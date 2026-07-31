# 阶段性工作、测试、缺口与 GPU 迁移审计

更新时间：2026-07-31（Asia/Shanghai）

## 结论

网络中断后的续跑已经完成了官方 CRISPR-Cas Atlas v1.0 下载、真实流式
解析、exact dedup、六阈值 MMseqs2 聚类、严格 40% cluster split、真实
subtype MSA、保守性统计，以及 72 条 ESM-IF1 pilot 候选的全 Atlas
新颖度审计。

当前可确认：

- Milestone 0 已完成；
- Milestone 1 的本地 CPU、环境、依赖和权重部分已完成，GPU 验收未完成；
- Milestone 2 的 ESM-IF1 pilot 和三模型真实 CPU smoke 已完成，匹配设计
  位置/新颖度的完整多模型矩阵未完成；
- Milestone 3 的官方数据、真实漏斗、dedup、六阈值聚类和严格 split 已
  完成，subtype-held-out/scaffold-held-out 辅助 split 尚未实现；
- Milestone 4 已取得真实 subtype MSA 和 coverage-gated conservation，
  但 scaffold 映射尚未完成；paired-repeat、真实 MI/APC 和 DCA 被 Atlas
  repeat orientation 缺失阻断；
- Milestone 5 已完成真正的 causal constrained ESM-IF1 pilot 和 Level 1
  候选新颖度审计，完整 baseline/ablation matrix 未完成；
- Milestone 6 的 provider-neutral 接口、mock E2E、迁移脚本和操作指南已
  完成源节点验证，目标 GPU 节点的实际传输、bootstrap 和 GPU run 未完成。

现有结果最高为 Level 2 inverse-folding compatibility。只有 14 条候选
通过当前预注册统计过滤并可标为 Level 1；这不等于功能有效。没有 Level 3
真实回折证据，没有 Level 4 wet-lab 证据，任何候选都不得称为“已验证有效
Cas13”。

## 1. 已完成的真实工作

### 1.1 仓库、规范与可复现基础

- `AGENTS.md`、`.agent/PLANS.md` 和首个 ExecPlan 已建立并持续更新；
- 项目树、MIT 自有代码许可、CITATION、CI、Snakemake、Makefile、配置、
  四套隔离环境、三类 lock、第三方/模型/数据 manifest 已建立；
- 正式 run 使用不可覆盖目录，记录 resolved config、seed、环境、硬件、
  git、输入/输出 hash、失败和 SUCCESS/FAILED；
- mock、fixture、real 和证据 Level 均显式记录；
- 当前本地 `main` 比 `origin/main` 超前；最终质量门和 bundle 验证通过后
  才 push。

### 1.2 本机实测资源

| 项目 | 实测 |
|---|---|
| GPU | NVIDIA GeForce RTX 4060 Laptop |
| 显存 | 8188 MiB，不是 32 GB |
| NVIDIA driver | 560.94 |
| `nvidia-smi` 驱动兼容 CUDA 标签 | 12.6 |
| CPU | 32 logical CPUs |
| 内存 | 15 GiB |
| 当前磁盘余量 | 约 710 GiB |
| 当前 GPU 状态 | OS 阻断；`/dev/dxg` 不存在 |

当前 `nvidia-smi` 返回 `GPU access blocked by the operating system`，
PyTorch `torch.cuda.is_available()` 为 false。因此本轮真实模型运行均
明确使用 CPU。

### 1.3 环境、第三方和权重

四个环境位于 `.tools/envs/`，不污染系统 Python：

- `analysis`；
- `esm_if1`；
- `ligandmpnn`；
- `bioinformatics`。

每个环境均有 conda explicit、conda list 和 pip freeze lock。真实可用工具
包括 MMseqs2 18.8cc5c、MAFFT 7.526、HMMER 3.4、Infernal 1.1.5、
seqkit 2.13.0、Foldseek 10.941cd33、TM-align 20240303 和 Git LFS 3.7.1。
PyTorch Geometric GCNConv CPU smoke 已通过。

固定上游：

- ESM：`2b369911bb5b4b0dda914521b9475cad1656b2ac`；
- SynTnpBs：`f3ea8e69c6f71baa56c4bb388e9df0489720f968`；
- ProteinMPNN：`8907e6671bfbfc92303b5f79c4b5e6ce47cdef57`；
- LigandMPNN：`26ec57ac976ade5379920dbd43c7f97a91cf82de`。

模型均从本地 checkpoint 加载且有 hash：

| 模型 | 状态 |
|---|---|
| ESM-IF1 142M | toy、6E9F、5XWP 真实 CPU score/sample 通过 |
| ProteinMPNN v48_020 | 6E9F 真实 CPU smoke 通过 |
| LigandMPNN ligand_mpnn | 6E9F RNA B/C atomic-context 真实 CPU smoke 通过 |
| LigandMPNN protein_mpnn | 权重已下载，独立 smoke 未运行 |
| LigandMPNN soluble_mpnn | 权重已下载，独立 smoke 未运行 |

### 1.4 实验结构和真实 inverse-folding

RCSB 真实结构 6E9F、5XWP、6E9E、5XWY 已下载、hash 和 QC；RNA 只用于
接触注释或 LigandMPNN atomic context，不会作为 ESM-IF1 蛋白链输入。

ESM-IF1 真实 CPU score：

- 6E9F chain A，864 aa：conditional log-likelihood
  `-1813.9333906933316`，perplexity `8.161760905453558`；
- 5XWP chain A，1125 aa：conditional log-likelihood
  `-2553.345790145657`，perplexity `9.675923652828677`。

真正的 decode-time constrained decoder 已验证：

- fixed token 在自回归到达该位置时被强制输入；
- 不是生成后覆盖；
- future fixed token 对早期位置不可见；
- all-fixed、partial-fixed、free、非法位置、长度不匹配和 seed 重现均有测试；
- 6E9F 真实 constrained sample 固定位置违反为 0。

权威 ESM pilot：

`results/runs/20260731-benchmark-experimental-a998ff40aa-ab6c9c5/`

- 2 scaffold × 3 constraint condition × 3 temperature × 4 sample；
- 72/72 candidate ID 唯一；
- fixed-position violations 为 0；
- `is_mock=false`，CPU，Level 2；
- 条件间 fixed-position 比例不同，因此 raw recovery 不能解释为方法胜负。

### 1.5 真实 Atlas 数据

官方源：

`data/raw/atlas/v1.0/crispr-cas-atlas-v1.0.json`

- size：5,267,508,328 bytes；
- SHA256：
  `5b4ba2fb99638d279e0c126100e19a4b77aba487b37b7df118e4bf4acd494720`；
- license：CC BY-NC 4.0；
- manifest status：`downloaded_verified`。

clean-provenance 真实解析 run：

`results/runs/20260731-atlas-processing-e8356ef7b5-eebc1a5-r001/`

| 漏斗项 | 数量 |
|---|---:|
| Atlas operons | 1,246,088 |
| Type VI operons | 11,707 |
| Cas effector annotations | 6,174,375 |
| Cas13 records | 12,353 |
| Cas13 exact unique | 4,070 |
| evolution-eligible exact unique | 3,500 |
| high-confidence Cas13–repeat pairs | 0 |
| ambiguous pairs | 11,727 |
| processing failures | 0 |

Cas13 record subtype：

- VI-B：5,163；
- VI-D：6,857；
- VI-F：166；
- VI-I：167；
- subtype conflict：40。

Atlas v1.0 没有可恢复的 direct-repeat orientation。项目没有猜测 strand，
所以 11,727 条进入 ambiguous 表，0 条进入高置信配对表。这是源数据阻断，
不是 GPU 或内存不足。

### 1.6 聚类、split、MSA 和 conservation

MMseqs2 使用 identity 1.0/0.9/0.7/0.5/0.4/0.3、coverage 0.8、
coverage mode 0、cluster mode 2、16 threads。

| identity | cluster 数 |
|---:|---:|
| 100% | 3,877 |
| 90% | 1,797 |
| 70% | 1,323 |
| 50% | 1,003 |
| 40% | 783 |
| 30% | 516 |

100% cluster 少于 4,070 exact unique 是因为 80% coverage 下可把完全一致的
短片段与较长序列聚在一起；exact hash 与 MMseqs cluster 是不同语义。

严格 40% cluster split：

- train：3,335 sequences；
- validation：160；
- test：575；
- leakage gate：passed。

inclusive MSA 暴露 `truncated=00` 仍包含 48–80 aa HMM fragment，且没有
90% coverage 列。该结果保留在 audit 目录，不能用于约束。随后在任何候选
test metric 前预注册 700–1600 aa 宽松全长门，并在每个 70% cluster 选择
最长合格成员。

正式 subtype MSA：

| subtype | sequences | columns | median length |
|---|---:|---:|---:|
| VI-B | 489 | 7,380 | 1,139 |
| VI-D | 182 | 3,524 | 954 |
| VI-F | 45 | 1,786 | 1,154 |
| VI-I | 50 | 2,164 | 1,084.5 |

coverage ≥ 0.8 的列：

- VI-B：696；
- VI-D：724；
- VI-F：1,095；
- VI-I：763。

这些只是“进入 scaffold mapping 的候选列”。在 scaffold-to-MSA mapping
和 mapping confidence 通过前，不会自动进入 hard/soft constraint。

### 1.7 候选新颖度

真实 run：

`results/runs/20260731-candidate-filtering-b14455d461-6d258de/`

对 72 条 pilot 候选使用全 4,070 exact-unique Atlas FASTA、MMseqs2
sensitivity 7.5 和 query coverage ≥ 0.8：

- 19 条有满足覆盖门的 Atlas hit；
- 观察到的最大 Atlas identity：0.457；
- parent identity 范围：0.1653–0.4213；
- 14 条通过全部预注册门，最高可标 Level 1；
- 41 条因 low-complexity window 失败；
- 18 条因 homopolymer > 7 失败；
- 53 条没有满足覆盖门的 Atlas hit，按 fail-closed 处理，不能把“未检出”
  自动解释为高度新颖。

5XWP 属 VI-A，而本次 Atlas Cas13 HMM 记录没有解析出 VI-A；其候选缺少
同 subtype Atlas 覆盖是数据库覆盖限制，不是候选优越性的证据。

## 2. 当前测试和报告状态

最近完整 CPU/fixture 质量门：

- Ruff lint：通过；
- Ruff format：通过；
- strict mypy：通过，43 个源文件；
- pytest：47/47 通过；
- branch-aware coverage：70.89%，门槛 70%；
- fixture MI/APC：通过，明确 `is_mock=true`，不是 DCA；
- mock refold E2E：通过，明确 `is_mock=true`。

当前项目报告：

`results/runs/20260731-benchmark-experimental-bcfd0be469-3ebd1c9/report/`

报告明确列出：

- available real：Atlas、cluster、MSA、conservation、candidate novelty、
  三模型 smoke、实验结构和 ESM benchmark；
- not run：matched multimodel benchmark、real refold；
- data-blocked：paired MSA、real MI/APC、formal DCA；
- maximum evidence：Level 2；
- 无 Level 4 功能声明。

## 3. 本节点仍可完成但尚未完成

以下工作不应因有 GPU 节点而推迟：

1. 把 VI-D scaffold 映射到正式 VI-D MSA，并生成 mapping confidence；
2. 实现 subtype-held-out 和 scaffold-held-out 辅助 split；
3. 对 MSA 最短长度 500/600/800 aa 做预注册 sensitivity；
4. 把 ProteinMPNN、LigandMPNN、MSA profile、matched random mutation
   扩展为相同设计位置/相同新颖度的 CPU 小规模矩阵；
5. 为 5XWP 运行真实 ProteinMPNN/LigandMPNN baseline；
6. 运行 LigandMPNN `protein_mpnn` 和 `soluble_mpnn` 独立 checkpoint smoke；
7. 实现当前仍是显式 `not_run` 的生产 CLI wrapper，包括通用
   `score/sample/sequence-qc`、真实 refold export/ingest 命令入口；
8. 完成候选 funnel、匹配统计、bootstrap/effect size 和 failure analysis；
9. 更新 manuscript 的实际方法数字和图表。

这些是尚未实现/运行，不应描述为已完成。

## 4. 本节点失败过的工作及原因

| 事项 | 原因 | 当前处理 |
|---|---|---|
| Atlas 首次下载 | DNS 无法解析 Google Storage | 网络恢复后续传完成，size/hash 验证通过 |
| fetch finalizer 首次幂等测试 | 系统 Python 3.10 没有 `datetime.UTC` | 改用 `timezone.utc`，系统 Python 与 analysis env 均通过 |
| clean Atlas 重跑首次启动 | 安全门拒绝覆盖已有 canonical 输出 | 保留旧目录后 clean run 成功，输出逐字节一致 |
| ESM 环境首次建成 | user-site `fair-esm` 泄漏 | pinned 本地源码安装并强制 `PYTHONNOUSERSITE=1` |
| ProteinMPNN 首次 validator | 把 893 numbering slots 错当 864 resolved residues | 显式保留 29 个 masked `X` slot 后通过 |
| LigandMPNN 前两次 | 缺 `dm-tree`；NumPy 1.26 移除 `np.int` | 按上游 requirements 加 `dm-tree`、pin NumPy 1.23.5 |
| inclusive MSA | `truncated=00` 仍含 48–80 aa fragment | 保留 audit；预注册全长门后重建正式 MSA |
| 本机 GPU smoke | WSL `/dev/dxg` 缺失，NVML 被 OS 阻断 | CPU 结果保留；真实 GPU smoke 迁移 |
| Cas13–repeat 共进化 | Atlas 缺 repeat orientation，0 高置信 pair | data-blocked；不猜 strand，不用 MI 冒充 DCA |
| 一篇 Cas13d 作者 PDF | 公开端点返回内容不是 PDF | 无效文件不入库，只保留元数据/合法链接 |

当前没有“测试失败但仍标记完成”的事项。

## 5. 应迁移到 GPU/HPC 的工作

优先或必须迁移：

- ESM-IF1、ProteinMPNN、LigandMPNN 的正式 GPU smoke；
- 大规模 ESM-IF1/LigandMPNN 采样和完整消融矩阵；
- 多状态大批量候选生成；
- AF2/ColabFold/AF3/Protenix/Boltz 多 seed 回折；
- 全 Atlas 结构预测和 Foldseek 结构聚类；
- ESM-IF1 Cas13 domain adaptation；
- 有合格 paired data 后的长序列 plmDCA/GREMLIN/CCMpred。

注意：当前 direct-repeat DCA 首先是源数据 orientation 阻断。迁到 GPU
不会自动解决 0 高置信 pair；必须先获得有授权、方向可信的 paired data。

## 6. GPU 迁移指南是否完成

源节点指南和接口已完成：

- `docs/GPU_MIGRATION.md`；
- `scripts/export_gpu_bundle.sh`；
- `scripts/verify_gpu_bundle.sh`；
- `scripts/bootstrap_gpu_node.sh`；
- `scripts/sync_assets.sh`；
- `scripts/launch_gpu_tmux.sh`；
- provider-neutral refold FASTA/JSONL、deterministic shard、retry、
  pLDDT/PAE/structure ingest、US-align comparison 和 missing-output audit；
- mock refold E2E，`is_mock=true`。

已验证的行为：

- bundle 不嵌入 GB 级资产；
- `SHA256SUMS` 校验内部文件；
- `ASSET_SHA256SUMS` 校验模型、结构和 Atlas；
- rsync 支持 `--partial`；
- target bootstrap 使用四个隔离环境；
- tmux launcher 记录唯一 session、run dir、完整日志、退出码和现场；
- 不无限静默重试。

本轮 Atlas-complete clean bundle：

`artifacts/bundles/gpu-bundle-3e53026923aa-7540febfb2/`

- manifest commit：
  `3e53026923aacdc9a87de1b7005dfa844d837934`；
- export dirty：false；
- bundle size：约 372 KiB；
- internal SHA256：passed；
- 5 个 checkpoint、Atlas JSON、8 个 PDB/mmCIF asset hash：passed；
- `missing_assets=[]`；
- 大资产没有嵌入 bundle。

尚未完成的是目标 GPU 节点验收：

- 第二台机器真实资产传输；
- 目标节点四环境 bootstrap；
- `torch.cuda.is_available()=true`；
- 三模型 GPU smoke；
- 实际 tmux 长任务；
- 真实 refold shard 执行和结果回收。

因此准确表述是：“迁移操作指南、脚本、源节点 bundle 机制和 mock E2E
已完成；目标 GPU 节点实测尚未完成。”

目标节点主流程：

```bash
git clone https://github.com/XuJP264/ESM-IF-on-cas13.git /work/ESM-IF
cd /work/ESM-IF
git checkout <bundle-manifest中的完整commit>
bash scripts/verify_gpu_bundle.sh <bundle-dir> /work/ESM-IF
bash scripts/bootstrap_gpu_node.sh all
make smoke-esm-if1
make smoke-proteinmpnn
make smoke-ligandmpnn
bash scripts/launch_gpu_tmux.sh \
  configs/benchmark_experimental.yaml \
  benchmark-experimental
```

## 7. 下一阶段顺序

1. 完成本轮 docs/STATUS/DECISIONS/ExecPlan 更新；
2. 重跑最终 lint、typecheck、test、CPU smoke 和 reproducibility gate；
3. push 并验证 GitHub Actions；
4. 本机继续 VI-D scaffold-to-MSA mapping 和小规模 matched baseline；
5. GPU 节点执行三模型 GPU smoke，再进入大规模采样/回折。
