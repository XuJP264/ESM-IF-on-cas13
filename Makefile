SHELL := /usr/bin/env bash
ANALYSIS_ENV ?= .tools/envs/analysis
RUN := conda run -p $(ANALYSIS_ENV)

.PHONY: bootstrap fetch-references fetch-third-party fetch-models fetch-atlas
.PHONY: fetch-structures lint typecheck test smoke-cpu smoke-esm-if1
.PHONY: smoke-proteinmpnn smoke-ligandmpnn process-atlas cluster msa
.PHONY: conservation coevolution-smoke benchmark-experimental generate-pilot
.PHONY: report export-gpu-bundle verify-reproducibility

bootstrap:
	bash scripts/bootstrap_local.sh

fetch-references:
	bash scripts/fetch_references.sh

fetch-third-party:
	bash scripts/fetch_third_party.sh

fetch-models:
	bash scripts/fetch_models.sh

fetch-atlas:
	bash scripts/fetch_atlas.sh

fetch-structures:
	bash scripts/fetch_experimental_structures.sh

lint:
	$(RUN) ruff check src tests
	$(RUN) ruff format --check src tests

typecheck:
	$(RUN) mypy src

test:
	$(RUN) pytest -m "not real_model and not network and not slow" --cov=cas13_if --cov-report=term-missing

smoke-cpu:
	$(RUN) cas13-if preflight --config configs/benchmark_experimental.yaml --fixture

smoke-esm-if1:
	bash scripts/run_real_smoke.sh esm-if1

smoke-proteinmpnn:
	bash scripts/run_real_smoke.sh proteinmpnn

smoke-ligandmpnn:
	bash scripts/run_real_smoke.sh ligandmpnn

process-atlas:
	$(RUN) cas13-if build-dataset --config configs/atlas_processing.yaml

cluster:
	$(RUN) cas13-if cluster --config configs/atlas_processing.yaml

msa:
	$(RUN) cas13-if build-msa --config configs/atlas_processing.yaml

conservation:
	$(RUN) cas13-if conservation --config configs/atlas_processing.yaml

coevolution-smoke:
	$(RUN) cas13-if coevolution --config configs/atlas_processing.yaml --fixture

benchmark-experimental:
	bash scripts/run_experimental_benchmark.sh

generate-pilot:
	$(RUN) cas13-if sample --config configs/esm_if1_sampling.yaml

report:
	$(RUN) cas13-if report --config configs/benchmark_experimental.yaml

export-gpu-bundle:
	bash scripts/export_gpu_bundle.sh

verify-reproducibility:
	bash scripts/run_static_validation.sh
	bash scripts/verify_gpu_bundle.sh
