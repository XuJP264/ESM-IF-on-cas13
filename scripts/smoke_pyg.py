#!/usr/bin/env python
"""Run a genuine, small PyTorch Geometric CPU operation and record it."""

from __future__ import annotations

import json
import time
from importlib.metadata import version
from pathlib import Path

import torch
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    torch.manual_seed(20260731)
    features = torch.tensor(
        [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        dtype=torch.float32,
    )
    edge_index = torch.tensor(
        [[0, 1, 1, 2], [1, 0, 2, 1]],
        dtype=torch.long,
    )
    graph = Data(x=features, edge_index=edge_index)
    layer = GCNConv(2, 3)
    started = time.perf_counter()
    with torch.inference_mode():
        output = layer(graph.x, graph.edge_index)
    elapsed = time.perf_counter() - started
    if output.shape != (3, 3) or not bool(torch.isfinite(output).all()):
        raise RuntimeError("PyTorch Geometric smoke produced an invalid tensor")
    result = {
        "schema_version": "1.0",
        "is_mock": False,
        "evidence_level": 0,
        "seed": 20260731,
        "device": str(output.device),
        "torch_version": torch.__version__,
        "torch_cuda_runtime": torch.version.cuda,
        "torch_cuda_available": torch.cuda.is_available(),
        "torch_geometric_version": version("torch-geometric"),
        "operation": "GCNConv(2,3) on a three-node bidirectional path",
        "output_shape": list(output.shape),
        "output_finite": True,
        "elapsed_seconds": elapsed,
    }
    output_path = repo / "artifacts/system/pytorch_geometric_real_smoke.json"
    output_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
