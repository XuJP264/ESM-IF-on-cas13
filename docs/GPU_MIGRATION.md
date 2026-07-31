# GPU/HPC 迁移与运行指南

## 适用范围与证据边界

迁移采用“小型代码/配置 bundle + 独立大资产同步”的方式。模型权重、Atlas
原始数据和预测结构不会被打进单个巨大 tar。GPU 节点输出仍遵守 Level 0–4
证据分级；结构预测或多模型一致性最高支持 Level 3，不能称为 wet-lab
验证有效的 Cas13。

本地 RTX 4060 Laptop 实测只有 8188 MiB 显存，而且当前 WSL 缺少
`/dev/dxg`。大规模采样、回折、正式 DCA 和领域适配应迁移；Atlas 下载、
流式解析、MMseqs2/MAFFT 及报告可以留在 CPU 节点。

## 1. 源节点准备

先确保待迁移代码已经 commit，且模型/数据 manifest 与实际文件一致：

```bash
git status --short
make lint
make typecheck
make test
make smoke-esm-if1
make smoke-proteinmpnn
make smoke-ligandmpnn
```

如已有候选回折任务，把包含 `candidates.fasta`、`jobs.jsonl`、
`expected_outputs.json` 和 `shards/` 的目录作为可选参数传入：

```bash
bash scripts/export_gpu_bundle.sh results/runs/<run_id>/refold_jobs
```

没有真实候选时也可导出基础 bundle：

```bash
make export-gpu-bundle
```

输出目录名包含当前 git short SHA 和配置 hash。脚本拒绝覆盖已有 bundle。
`bundle-manifest.json` 会记录：

- 完整 git commit 和导出时工作树是否 dirty；
- 环境规格与 exact/pip/conda locks；
- container、运行脚本、配置和 refold output schema；
- 每个已存在的大资产的相对路径、大小和 SHA256；
- 缺失资产清单；
- 可选输入 shard 清单。

内部文件由 `SHA256SUMS` 校验，大资产由 `ASSET_SHA256SUMS` 校验。正式迁移
应要求 `git_worktree_dirty_at_export=false`；dirty bundle 只用于诊断。

## 2. 传输 bundle 和大资产

先把小 bundle 复制到 GPU 节点，例如：

```bash
rsync -a artifacts/bundles/gpu-bundle-<sha>-<hash>/ \
  gpu-host:/work/ESM-IF/artifacts/bundles/gpu-bundle-<sha>-<hash>/
```

在 GPU 节点 clone bundle 指定的精确 commit：

```bash
git clone https://github.com/XuJP264/ESM-IF-on-cas13.git /work/ESM-IF
cd /work/ESM-IF
git checkout <full-commit-from-bundle-manifest>
```

大资产从源资产根目录单独同步。`SOURCE_ROOT` 可以是 rsync 的远端路径，
`TARGET_ROOT` 是 GPU 节点 repo 根目录，第三个参数是已经复制的 bundle：

```bash
bash scripts/sync_assets.sh \
  user@source-host:/home/junpeng/ESM-IF \
  /work/ESM-IF \
  /work/ESM-IF/artifacts/bundles/gpu-bundle-<sha>-<hash>
```

同步使用 `rsync --partial`，中断后可重复执行。完成后脚本会依据
`ASSET_SHA256SUMS` 对目标文件逐一校验；任何缺失或 hash 不一致都会非零
退出。不要把模型权重、Atlas JSON 或预测结构提交 Git。

## 3. GPU 节点预检和环境恢复

节点最低要求：

- `git`、`conda`、`tmux`、`rsync`；
- 可工作的 NVIDIA driver 和 `nvidia-smi`；
- 足够磁盘；
- 已同步 pinned `third_party` checkout 和 checkpoint。

先验证 bundle 内部，再验证同步后的大资产：

```bash
cd /work/ESM-IF
bash scripts/verify_gpu_bundle.sh \
  artifacts/bundles/gpu-bundle-<sha>-<hash> \
  /work/ESM-IF
```

建立四个项目内隔离环境：

```bash
bash scripts/bootstrap_gpu_node.sh all
```

也可分别建立：

```bash
bash scripts/bootstrap_gpu_node.sh analysis
bash scripts/bootstrap_gpu_node.sh esm-if1
bash scripts/bootstrap_gpu_node.sh ligandmpnn
bash scripts/bootstrap_gpu_node.sh bioinformatics
```

脚本强制 `PYTHONNOUSERSITE=1`，ESM 从 pinned 本地 checkout 安装，并重新
导出 conda explicit、conda list 和 pip freeze。完成后必须确认：

```bash
nvidia-smi
.tools/envs/esm_if1/bin/python -c \
  'import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())'
.tools/envs/ligandmpnn/bin/python -c \
  'import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())'
```

若 `torch.cuda.is_available()` 为 false，不得继续把任务记录为 GPU 运行；
先修复 driver/container/调度器映射。

## 4. GPU smoke 和长任务

先执行单结构 smoke：

```bash
make smoke-esm-if1
make smoke-proteinmpnn
make smoke-ligandmpnn
```

所有 GPU 长任务必须通过 tmux 启动：

```bash
bash scripts/launch_gpu_tmux.sh \
  configs/benchmark_experimental.yaml \
  benchmark-experimental
```

脚本生成唯一 session 和 `results/runs/<timestamp>-<task>-gpu/`，记录：

- `TASK`、`CONFIG`、`GIT_COMMIT`；
- `nvidia-smi.txt`、`STARTED_AT`、`FINISHED_AT`；
- `stdout.log`、`stderr.log`、`EXIT_CODE`；
- `SUCCESS` 或 `FAILED`。

查看和恢复：

```bash
cat results/latest_gpu_session.txt
cat results/latest_gpu_run.txt
tmux attach -t "$(cat results/latest_gpu_session.txt)"
tail -f "$(cat results/latest_gpu_run.txt)/stdout.log"
```

脚本不自动无限重试。失败时保留 run 目录，修复原因后创建新 run；不要覆盖
旧 run。

## 5. Refold shard 运行与结果回收

每个 shard 必须保持原始 `candidate_id`、provider、seed 和 `is_mock`。
provider 输出应满足 `src/cas13_if/refold/output_schema.json`，每个候选至少
提供：

- `result.json`；
- 结构文件；
- PAE 文件；
- mean pLDDT；
- provider、seed 和 mock 状态。

只同步结果目录和日志回源节点，例如：

```bash
rsync -a --partial \
  gpu-host:/work/ESM-IF/results/runs/<gpu_run_id>/ \
  results/runs/<gpu_run_id>/
```

回收后执行 ingest/QC。缺失结构、PAE 或 result JSON 必须进入
`failed_job_retry.jsonl`，不得静默丢失。US-align/TM-align 不可用时
TM-score 状态必须是 `not_run`，不能用近似分数替代。

## 6. 当前未完成的 GPU 验收

截至 2026-07-31，bundle/refold 接口和本指南已实现，但以下项目仍需在真实
GPU 节点取证后才能标记完成：

- GPU 节点四环境 bootstrap；
- ESM-IF1、ProteinMPNN、LigandMPNN GPU smoke；
- tmux 正式 benchmark 的完整退出码和日志审计；
- 真实 candidate shard 的 provider 执行与回收；
- 大规模 AF/ColabFold/AF3/Protenix/Boltz 回折；
- 正式 DCA、Foldseek 全量结构聚类和 ESM-IF1 domain adaptation。

早期 clean bundle
`artifacts/bundles/gpu-bundle-7dc0491d8441-6ad46d8577/` 和
`artifacts/bundles/gpu-bundle-ab6c9c5d011d-6ad46d8577/` 已完成内部
checksum 及当时已有模型/结构资产校验；它们是在 Atlas 完整下载前导出，
所以不应作为当前生产数据 bundle。官方 Atlas 现在已经完整下载并写入
manifest，大小为 5,267,508,328 bytes，SHA256 为
`5b4ba2fb99638d279e0c126100e19a4b77aba487b37b7df118e4bf4acd494720`。
本阶段最终验收会从最新 clean commit 再导出一次 bundle，并要求 Atlas
出现在 `ASSET_SHA256SUMS` 而不是 `missing_assets`。

当前状态是“操作指南、源节点 bundle 机制和本地校验完成，目标 GPU 节点
真实验收未完成”，而不是“GPU 实验已完成”。此外，Atlas 缺失
direct-repeat orientation 是源数据阻断；把 DCA 任务搬到 GPU 不能解决
0 个高置信配对的问题。
