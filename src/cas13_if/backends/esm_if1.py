"""Offline ESM-IF1 scoring and causal constrained sampling.

The implementation loads a user-fetched local checkpoint and never invokes the
hub loader. RNA chains are intentionally excluded: ESM-IF1 consumes protein
N/CA/C backbones only. Fixed tokens are inserted at their autoregressive decode
step, matching the upstream ``partial_seq`` semantics used by SynTnpB.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import math
import sys
import types
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from cas13_if.backends.base import InverseFoldingBackend
from cas13_if.schemas import (
    STANDARD_AA,
    BackendCapabilities,
    Candidate,
    EvidenceLevel,
    PositionTrace,
    SampleRequest,
    ScoreRequest,
    ScoreResult,
)
from cas13_if.structures.parser import (
    group_residues,
    parse_structure,
    protein_chain_sequence,
)

TRACE_ALPHABET = tuple(sorted(STANDARD_AA))


class EsmIf1Backend(InverseFoldingBackend):
    """Genuine, unconstrained ESM-IF1 backend using an offline checkpoint."""

    backend_name = "esm_if1"
    supports_constraints = False

    def __init__(self, checkpoint: Path, *, device: str = "auto") -> None:
        self.checkpoint = checkpoint.resolve()
        self.requested_device = device
        self._device = "not_loaded"
        self._model: Any = None
        self._alphabet: Any = None
        self._torch: Any = None
        self._batch_converter_class: Any = None
        self._checkpoint_sha256: str | None = None

    def capabilities(self) -> BackendCapabilities:
        return BackendCapabilities(
            scoring=True,
            sampling=True,
            protein_multichain=True,
            rna_atomic_context=False,
            hard_fixed=self.supports_constraints,
            allowed_residue_filter=self.supports_constraints,
            per_residue_probabilities=True,
        )

    def load(self) -> None:
        if not self.checkpoint.is_file():
            raise FileNotFoundError(
                f"local ESM-IF1 checkpoint is missing: {self.checkpoint}; "
                "run scripts/fetch_models.sh explicitly"
            )
        torch = importlib.import_module("torch")
        esm = importlib.import_module("esm")
        esm_data = importlib.import_module("esm.data")
        _install_inverse_folding_runtime_compatibility(esm_data, torch)
        if self.requested_device == "auto":
            device = "cuda" if bool(torch.cuda.is_available()) else "cpu"
        else:
            device = self.requested_device
        if device.startswith("cuda") and not bool(torch.cuda.is_available()):
            raise RuntimeError(
                "CUDA was requested but torch.cuda.is_available() is false"
            )
        # This is deliberately the local loader. The named pretrained helper is
        # not used because it is allowed to fetch weights implicitly.
        model, alphabet = esm.pretrained.load_model_and_alphabet_local(
            str(self.checkpoint)
        )
        model = model.eval().to(torch.device(device))
        self._torch = torch
        self._model = model
        self._alphabet = alphabet
        self._batch_converter_class = esm_data.BatchConverter
        self._device = device
        self._checkpoint_sha256 = _sha256(self.checkpoint)

    def _require_loaded(self) -> None:
        if self._model is None:
            raise RuntimeError(f"{self.__class__.__name__}.load() must be called first")

    def score(self, request: ScoreRequest) -> ScoreResult:
        self._require_loaded()
        coords, target_chain, conditioning_chains = self._load_conditioning_coords(
            request.structure_path, request.protein_chains
        )
        target_length = self._target_length(coords, target_chain)
        if len(request.sequence) != target_length:
            raise ValueError(
                f"sequence length {len(request.sequence)} does not match target "
                f"chain {target_chain!r} coordinate length {target_length}"
            )
        all_coords = self._concatenate(coords, target_chain)
        with self._torch.inference_mode():
            loss, padding_mask = self._get_sequence_loss(all_coords, request.sequence)
        valid_losses = [
            float(value)
            for value, is_padding in zip(
                loss.tolist(), padding_mask.tolist(), strict=True
            )
            if not bool(is_padding)
        ]
        if len(valid_losses) != len(request.sequence):
            raise RuntimeError(
                "ESM-IF1 returned a per-residue loss length inconsistent with "
                f"the target sequence: {len(valid_losses)} != {len(request.sequence)}"
            )
        per_residue = [-value for value in valid_losses]
        mean_log_likelihood = sum(per_residue) / len(per_residue)
        return ScoreResult(
            scaffold_id=request.scaffold_id,
            backend=self.backend_name,
            sequence=request.sequence,
            conditional_log_likelihood=sum(per_residue),
            perplexity=math.exp(-mean_log_likelihood),
            per_residue_log_probabilities=per_residue,
            is_mock=False,
            evidence_level=EvidenceLevel.INVERSE_FOLDING_COMPATIBILITY,
            metadata={
                "mean_conditional_log_likelihood": mean_log_likelihood,
                "target_chain": target_chain,
                "conditioning_protein_chains": conditioning_chains,
                "rna_atomic_context": False,
                "checkpoint_sha256": self._checkpoint_sha256,
                "device": self._device,
            },
        )

    def sample(self, request: SampleRequest) -> list[Candidate]:
        self._require_loaded()
        if (
            request.fixed_positions or request.allowed_residues
        ) and not self.supports_constraints:
            raise ValueError(
                "unconstrained ESM-IF1 backend does not accept residue constraints; "
                "use EsmIf1ConstrainedBackend"
            )
        coords, target_chain, conditioning_chains = self._load_conditioning_coords(
            request.structure_path, request.protein_chains
        )
        target_length = self._target_length(coords, target_chain)
        if len(request.parent_sequence) != target_length:
            raise ValueError(
                f"parent length {len(request.parent_sequence)} does not match "
                f"target chain coordinate length {target_length}"
            )
        all_coords = self._concatenate(coords, target_chain)
        candidates: list[Candidate] = []
        request_digest = _sampling_request_digest(self.backend_name, request)
        for sample_index in range(request.count):
            sample_seed = request.seed + sample_index
            sequence, traces = self._decode(
                all_coords=all_coords,
                target_length=target_length,
                fixed_positions=request.fixed_positions,
                allowed_residues=request.allowed_residues,
                temperature=request.temperature,
                seed=sample_seed,
            )
            candidates.append(
                Candidate(
                    candidate_id=(
                        f"{self.backend_name}-{request.scaffold_id}-"
                        f"{request_digest}-{sample_index:04d}"
                    ),
                    scaffold_id=request.scaffold_id,
                    backend=self.backend_name,
                    sequence=sequence,
                    parent_sequence=request.parent_sequence,
                    seed=sample_seed,
                    temperature=request.temperature,
                    is_mock=False,
                    evidence_level=EvidenceLevel.INVERSE_FOLDING_COMPATIBILITY,
                    fixed_positions={
                        index: token.upper()
                        for index, token in request.fixed_positions.items()
                    },
                    traces=traces,
                    metadata={
                        "semantics": "left_to_right_causal",
                        "future_fixed_tokens_visible": False,
                        "target_chain": target_chain,
                        "conditioning_protein_chains": conditioning_chains,
                        "rna_atomic_context": False,
                        "checkpoint_sha256": self._checkpoint_sha256,
                        "trace_alphabet": list(TRACE_ALPHABET),
                        "device": self._device,
                    },
                )
            )
        return candidates

    def _load_conditioning_coords(
        self, structure_path: str, chains: list[str]
    ) -> tuple[Any, str, list[str]]:
        path = Path(structure_path)
        if not path.is_file():
            raise FileNotFoundError(f"structure is missing: {path}")
        if not chains:
            raise ValueError(
                "protein_chains must explicitly name the target protein chain first; "
                "RNA chains must not be supplied to ESM-IF1"
            )
        parsed = parse_structure(path)
        chain_coords = {
            chain: self._extract_chain_coords(parsed, chain) for chain in chains
        }
        if len(chains) == 1:
            return chain_coords[chains[0]], chains[0], list(chains)
        return chain_coords, chains[0], list(chains)

    def _extract_chain_coords(self, atoms: Any, chain: str) -> NDArray[np.float32]:
        sequence, residue_keys = protein_chain_sequence(atoms, chain)
        if not sequence:
            raise ValueError(f"protein chain {chain!r} has no residues")
        residues = group_residues(atoms)
        coords = np.full((len(residue_keys), 3, 3), np.nan, dtype=np.float32)
        for residue_index, key in enumerate(residue_keys):
            by_name = {
                atom.name: atom
                for atom in residues[key]
                if atom.name in {"N", "CA", "C"}
            }
            for atom_index, atom_name in enumerate(("N", "CA", "C")):
                atom = by_name.get(atom_name)
                if atom is not None:
                    coords[residue_index, atom_index] = atom.coordinate
        return coords

    def _target_length(self, coords: Any, target_chain: str) -> int:
        if isinstance(coords, dict):
            return int(coords[target_chain].shape[0])
        return int(coords.shape[0])

    def _concatenate(self, coords: Any, target_chain: str) -> Any:
        if isinstance(coords, dict):
            padding = np.full((10, 3, 3), np.nan, dtype=np.float32)
            arrays = [coords[target_chain]]
            for chain, chain_coords in coords.items():
                if chain != target_chain:
                    arrays.extend((padding, chain_coords))
            return np.concatenate(arrays, axis=0)
        return coords

    def _convert_batch(
        self,
        all_coords: Any,
        sequence: str | None,
    ) -> tuple[Any, Any, Any, Any]:
        torch = self._torch
        dictionary = self._model.decoder.dictionary
        dictionary.cls_idx = dictionary.get_idx("<cath>")
        converter = self._batch_converter_class(dictionary)
        sequence_value = sequence if sequence is not None else "X" * len(all_coords)
        confidence_values = [1.0] * len(all_coords)
        labels, _, tokens = converter(
            [((all_coords, confidence_values), sequence_value)]
        )
        coordinates, confidence = labels[0]
        coordinate_tensor = torch.nn.functional.pad(
            torch.tensor(coordinates),
            (0, 0, 0, 0, 1, 1),
            value=np.inf,
        ).unsqueeze(0)
        confidence_tensor = torch.nn.functional.pad(
            torch.tensor(confidence),
            (1, 1),
            value=-1.0,
        ).unsqueeze(0)
        device = torch.device(self._device)
        coordinate_tensor = coordinate_tensor.to(device)
        confidence_tensor = confidence_tensor.to(device)
        tokens = tokens.to(device)
        padding_mask = torch.isnan(coordinate_tensor[:, :, 0, 0])
        coordinate_mask = torch.isfinite(coordinate_tensor.sum(-2).sum(-1))
        confidence_tensor = confidence_tensor * coordinate_mask + (-1.0) * padding_mask
        return coordinate_tensor, confidence_tensor, tokens, padding_mask

    def _get_sequence_loss(
        self,
        all_coords: Any,
        sequence: str,
    ) -> tuple[Any, Any]:
        torch = self._torch
        batch_coords, confidence, tokens, padding_mask = self._convert_batch(
            all_coords, sequence
        )
        previous_tokens = tokens[:, :-1]
        target = tokens[:, 1:]
        target_padding_mask = target == self._alphabet.padding_idx
        logits, _ = self._model.forward(
            batch_coords,
            padding_mask,
            confidence,
            previous_tokens,
        )
        loss = torch.nn.functional.cross_entropy(logits, target, reduction="none")
        return (
            loss[0].cpu().detach().numpy(),
            target_padding_mask[0].cpu().numpy(),
        )

    def _decode(
        self,
        *,
        all_coords: Any,
        target_length: int,
        fixed_positions: dict[int, str],
        allowed_residues: dict[int, set[str]],
        temperature: float,
        seed: int,
    ) -> tuple[str, list[PositionTrace]]:
        torch = self._torch
        device = torch.device(self._device)
        dictionary = self._model.decoder.dictionary
        batch_coords, confidence, _, padding_mask = self._convert_batch(
            all_coords, None
        )
        mask_index = dictionary.get_idx("<mask>")
        sampled_tokens = torch.full(
            (1, 1 + len(all_coords)),
            mask_index,
            dtype=torch.long,
            device=device,
        )
        sampled_tokens[0, 0] = dictionary.get_idx("<cath>")
        # Context chains are padding tokens. Only the first target_length
        # positions are decoded as amino acids.
        if len(all_coords) > target_length:
            sampled_tokens[0, 1 + target_length :] = dictionary.padding_idx
        for index, token in fixed_positions.items():
            sampled_tokens[0, index + 1] = dictionary.get_idx(token.upper())

        torch.manual_seed(seed)
        if self._device.startswith("cuda"):
            torch.cuda.manual_seed_all(seed)
        trace_indices = [dictionary.get_idx(token) for token in TRACE_ALPHABET]
        traces: list[PositionTrace] = []
        incremental_state: dict[str, Any] = {}
        with torch.inference_mode():
            encoder_out = self._model.encoder(batch_coords, padding_mask, confidence)
            for decoder_position in range(1, target_length + 1):
                output, _ = self._model.decoder(
                    sampled_tokens[:, :decoder_position],
                    encoder_out,
                    incremental_state=incremental_state,
                )
                raw_logits = output[0, :, -1]
                filtered_logits = raw_logits / temperature
                index_0 = decoder_position - 1
                allowed = allowed_residues.get(index_0)
                if allowed is not None:
                    allowed_indices = {
                        dictionary.get_idx(token.upper()) for token in allowed
                    }
                    disallowed = torch.ones_like(filtered_logits, dtype=torch.bool)
                    disallowed[list(allowed_indices)] = False
                    filtered_logits = filtered_logits.masked_fill(
                        disallowed, float("-inf")
                    )
                probabilities = torch.softmax(filtered_logits, dim=-1)
                is_fixed = index_0 in fixed_positions
                if not is_fixed:
                    sampled_tokens[0, decoder_position] = torch.multinomial(
                        probabilities, 1
                    ).squeeze(-1)
                selected_index = int(
                    sampled_tokens[0, decoder_position].detach().cpu().item()
                )
                selected_token = str(dictionary.get_tok(selected_index))
                if selected_token not in STANDARD_AA:
                    raise RuntimeError(
                        "ESM-IF1 sampled a non-standard token "
                        f"{selected_token!r} at target position {index_0}"
                    )
                traces.append(
                    PositionTrace(
                        index=index_0,
                        logits=[
                            float(raw_logits[index].detach().cpu().item())
                            for index in trace_indices
                        ],
                        probabilities=[
                            float(probabilities[index].detach().cpu().item())
                            for index in trace_indices
                        ],
                        selected_token=selected_token,
                        fixed=is_fixed,
                        temperature=temperature,
                        seed=seed,
                    )
                )
        sequence = "".join(trace.selected_token for trace in traces)
        return sequence, traces

    def metadata(self) -> dict[str, Any]:
        return {
            "backend": self.backend_name,
            "checkpoint": str(self.checkpoint),
            "checkpoint_sha256": self._checkpoint_sha256,
            "device": self._device,
            "loaded": self._model is not None,
            "is_mock": False,
            "rna_atomic_context": False,
            "offline_local_checkpoint": True,
        }


class EsmIf1ConstrainedBackend(EsmIf1Backend):
    """ESM-IF1 with decode-time hard-fixed/allowed-residue constraints."""

    backend_name = "esm_if1_constrained"
    supports_constraints = True


def _sampling_request_digest(backend_name: str, request: SampleRequest) -> str:
    payload = {
        "backend": backend_name,
        "scaffold_id": request.scaffold_id,
        "parent_sequence_sha256": hashlib.sha256(
            request.parent_sequence.encode("ascii")
        ).hexdigest(),
        "temperature": request.temperature,
        "seed": request.seed,
        "fixed_positions": sorted(request.fixed_positions.items()),
        "allowed_residues": [
            [index, sorted(residues)]
            for index, residues in sorted(request.allowed_residues.items())
        ],
        "protein_chains": request.protein_chains,
    }
    canonical = json.dumps(payload, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:12]


def _sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _install_inverse_folding_runtime_compatibility(
    esm_data: Any,
    torch: Any,
) -> None:
    """Install pinned ESM math utilities without its optional Biotite parser.

    ESM-IF1 imports sequence-independent tensor utilities from a module that
    also imports Biotite eagerly. This project supplies the same pinned tensor
    operations and coordinate batch converter, while structure I/O remains in
    ``cas13_if.structures.parser``. No model operation is approximated.
    """
    module_name = "esm.inverse_folding.util"
    if module_name in sys.modules:
        return
    functional = torch.nn.functional
    utility = types.ModuleType(module_name)

    def nan_to_num(ts: Any, val: float = 0.0) -> Any:
        replacement = torch.tensor(val, dtype=ts.dtype, device=ts.device)
        return torch.where(~torch.isfinite(ts), replacement, ts)

    def norm(
        tensor: Any,
        dim: int,
        eps: float = 1e-8,
        keepdim: bool = False,
    ) -> Any:
        return torch.sqrt(
            torch.sum(
                torch.square(tensor),
                dim=dim,
                keepdim=keepdim,
            )
            + eps
        )

    def normalize(tensor: Any, dim: int = -1) -> Any:
        return nan_to_num(
            torch.div(
                tensor,
                norm(tensor, dim, keepdim=True),
            )
        )

    def rbf(
        values: Any,
        v_min: float,
        v_max: float,
        n_bins: int = 16,
    ) -> Any:
        centers = torch.linspace(v_min, v_max, n_bins, device=values.device)
        centers = centers.view([1] * len(values.shape) + [-1])
        standard_deviation = (v_max - v_min) / n_bins
        z_score = (values.unsqueeze(-1) - centers) / standard_deviation
        return torch.exp(-(z_score**2))

    def rotate(vector: Any, rotation: Any) -> Any:
        rotation = rotation.unsqueeze(-3)
        vector = vector.unsqueeze(-1)
        return torch.sum(vector * rotation, dim=-2)

    def get_rotation_frames(coords: Any) -> Any:
        first = coords[:, :, 2] - coords[:, :, 1]
        second = coords[:, :, 0] - coords[:, :, 1]
        axis_one = normalize(first, dim=-1)
        residual = second - axis_one * torch.sum(
            axis_one * second, dim=-1, keepdim=True
        )
        axis_two = normalize(residual, dim=-1)
        axis_three = torch.cross(axis_one, axis_two, dim=-1)
        return torch.stack([axis_one, axis_two, axis_three], dim=-2)

    class CoordBatchConverter(esm_data.BatchConverter):  # type: ignore[misc]
        def __call__(self, raw_batch: Any, device: Any = None) -> Any:
            self.alphabet.cls_idx = self.alphabet.get_idx("<cath>")
            batch = []
            for coords, confidence, sequence in raw_batch:
                if confidence is None:
                    confidence = 1.0
                if isinstance(confidence, (float, int)):
                    confidence = [float(confidence)] * len(coords)
                if sequence is None:
                    sequence = "X" * len(coords)
                batch.append(((coords, confidence), sequence))
            labels, strings, tokens = super().__call__(batch)
            coordinates = [
                functional.pad(
                    torch.tensor(item),
                    (0, 0, 0, 0, 1, 1),
                    value=np.inf,
                )
                for item, _ in labels
            ]
            confidences = [
                functional.pad(
                    torch.tensor(item),
                    (1, 1),
                    value=-1.0,
                )
                for _, item in labels
            ]
            coordinate_tensor = self.collate_dense_tensors(coordinates, pad_v=np.nan)
            confidence_tensor = self.collate_dense_tensors(confidences, pad_v=-1.0)
            if device is not None:
                coordinate_tensor = coordinate_tensor.to(device)
                confidence_tensor = confidence_tensor.to(device)
                tokens = tokens.to(device)
            padding_mask = torch.isnan(coordinate_tensor[:, :, 0, 0])
            coordinate_mask = torch.isfinite(coordinate_tensor.sum(-2).sum(-1))
            confidence_tensor = (
                confidence_tensor * coordinate_mask + (-1.0) * padding_mask
            )
            return (
                coordinate_tensor,
                confidence_tensor,
                strings,
                tokens,
                padding_mask,
            )

        @staticmethod
        def collate_dense_tensors(samples: list[Any], pad_v: float) -> Any:
            if not samples:
                return torch.Tensor()
            if len({sample.dim() for sample in samples}) != 1:
                raise RuntimeError("samples have varying dimensions")
            maximum_shape = [
                max(dimensions)
                for dimensions in zip(
                    *(sample.shape for sample in samples), strict=True
                )
            ]
            result = torch.empty(
                len(samples),
                *maximum_shape,
                dtype=samples[0].dtype,
                device=samples[0].device,
            )
            result.fill_(pad_v)
            for index, sample in enumerate(samples):
                result[index][tuple(slice(0, length) for length in sample.shape)] = (
                    sample
                )
            return result

        def from_lists(
            self,
            coords_list: Any,
            confidence_list: Any = None,
            seq_list: Any = None,
            device: Any = None,
        ) -> Any:
            batch_size = len(coords_list)
            if confidence_list is None:
                confidence_list = [None] * batch_size
            if seq_list is None:
                seq_list = [None] * batch_size
            return self.__call__(
                zip(coords_list, confidence_list, seq_list, strict=True),
                device,
            )

    utility_any: Any = utility
    utility_any.nan_to_num = nan_to_num
    utility_any.norm = norm
    utility_any.normalize = normalize
    utility_any.rbf = rbf
    utility_any.rotate = rotate
    utility_any.get_rotation_frames = get_rotation_frames
    utility_any.CoordBatchConverter = CoordBatchConverter
    sys.modules[module_name] = utility
    sys.modules.setdefault(
        "esm.inverse_folding.multichain_util",
        types.ModuleType("esm.inverse_folding.multichain_util"),
    )
