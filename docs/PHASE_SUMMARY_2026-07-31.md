# 阶段性工作、测试与迁移审计（2026-07-31）

## 结论

当前项目已完成 Milestone 0，并完成了 Milestone 1 的本地 CPU/依赖部分及
Milestone 2 的 ESM-IF1 实验结构 pilot benchmark。真实 ESM-IF1 toy、
6E9F、5XWP CPU 推理，6E9F 受约束 CPU 生成，ProteinMPNN 6E9F CPU
smoke 和 LigandMPNN RNA-context 6E9F CPU smoke 已成功。修正候选 ID
后的正式 ESM-IF1 benchmark 产生 72/72 个唯一候选，固定位置违反为 0。
这些结果最高支持 Level 2 的模型兼容性，不构成 Level 4 功能验证。
Milestone 3 的真实 Atlas 下载已经从 DNS 中断点重新启动，当前为可续传
`.part`，所以真实数据漏斗、聚类、MSA 和共进化结果仍为 `not_run`。

本文件是网络中断后的审计快照，不把“代码已实现”写成“实验已完成”，也
不把 fixture/mock 结果写成真实科学结果。

## 1. 已完成并验证的工作

### 1.1 仓库、规范与研究边界

- 已建立 `AGENTS.md`、`.agent/PLANS.md` 和
  `docs/execplans/0001_bootstrap_and_real_baselines.md`。
- 已建立完整项目树、MIT 代码许可证、CITATION、README、ROADMAP、
  Snakemake、Makefile、CI、配置、环境、容器、数据/模型/第三方 manifest、
  manuscript 和报告目录。
- 全局采用 Level 0–4 证据分级。任何 Level 1–3 计算候选都不得描述成
  “已验证有效 Cas13”。
- 已提交 M0：
  `7599be0 chore: initialize publication-grade research repository`。
- 当前分支为 `main`，远端为
  `https://github.com/XuJP264/ESM-IF-on-cas13.git`。新增工作尚未提交或
  push。

### 1.2 机器审计

真实审计保存在：

- `artifacts/system/hardware.json`
- `artifacts/system/software_initial.txt`

初始实测为：

- NVIDIA GeForce RTX 4060 Laptop GPU；
- 显存 8188 MiB，不是 32 GB；
- NVIDIA driver 560.94，`nvidia-smi` 显示的驱动兼容 CUDA 为 12.6；
- 32 个逻辑 CPU；
- 15 GiB 系统内存；
- 当前磁盘约 723 GiB 可用；
- Docker 28.5.1、Apptainer 1.4.5、Conda 26.3.2 可用；
- `nvcc`、`gh`、可用 Git LFS 及多数生物信息工具在初始审计时缺失。

网络中断后的复查中，`nvidia-smi` 返回
`GPU access blocked by the operating system`。这与初始成功审计不同，
必须在真实 GPU 推理前重新诊断，不能宣称本机 CUDA 已通过。

### 1.3 当前 CPU/fixture 质量门

网络中断后重新实际执行：

```text
make lint
make typecheck
make test
```

结果：

- Ruff lint：通过；
- Ruff format check：通过，65 个 Python 文件已格式化；
- strict mypy：通过，42 个源文件无错误；
- pytest：42/42 通过；
- branch-aware coverage：70.08%，达到项目设定的 70% 门槛。

这些测试覆盖：

- Atlas 流式 parser、方向处理、保守配对、模糊配对和 exact dedup；
- MMseqs cluster mapping、cluster-level split 和泄漏阻断；
- FASTA、真实 alignment 合法性、gap、sequence weighting；
- conservation、entropy、MI/APC；
- PDB/mmCIF、蛋白/RNA 链、insertion code、缺失主链和 RNA contact；
- Shrake–Rupley 相对 SASA；
- hard/soft/free mask、位置索引；
- constrained decoder、固定位置保持、seed 重现；
- backend/candidate/score schema；
- novelty、refold export/ingest、provenance；
- fixture Atlas → candidate → mock refold → report；
- fixture subtype MSA → conservation。

Mock E2E 明确为 `is_mock=true`、Level 0，不支持科学性能结论。

### 1.4 第三方代码、论文与模型资产

以下仓库已获取到固定 commit，源码保持未修改：

- ESM：`2b369911bb5b4b0dda914521b9475cad1656b2ac`
- SynTnpBs：`f3ea8e69c6f71baa56c4bb388e9df0489720f968`
- ProteinMPNN：`8907e6671bfbfc92303b5f79c4b5e6ce47cdef57`
- LigandMPNN：`26ec57ac976ade5379920dbd43c7f97a91cf82de`

真实 checkpoint 已下载并记录：

| 模型 | 大小 | SHA256 | 当前状态 |
|---|---:|---|---|
| ESM-IF1 142M | 1,700,450,121 B | `be4ba36edec22a9bfaa4946ff6b2815f1f19d8a3d7e0eada8b796d5a0eae9fd4` | 真实 CPU toy、6E9F、5XWP smoke 通过；GPU 被 OS 阻断 |
| ProteinMPNN v48_020 | 6,681,301 B | `c9cb4a671d79604111231f8dbfc7c590e06f1197453b7a6854ac6661a642f5bd` | 真实 6E9F CPU smoke 通过；GPU 被 OS 阻断 |
| LigandMPNN protein_mpnn | 6,681,301 B | 同上 | 已下载，未运行真实 smoke |
| LigandMPNN ligand_mpnn | 10,541,943 B | `161cd264061fda9680cbb940255522ae42f2966c552d045d87913d9452a80970` | 真实 6E9F RNA-context CPU smoke 通过；GPU 被 OS 阻断 |
| LigandMPNN soluble_mpnn | 6,681,301 B | `7af52d090172c230c7f0e9d21e02203f6b3a38b16db58d3c7a3960e0a9a6e31a` | 已下载，未运行真实 smoke |

ESM-IF1 和 LigandMPNN 的开放论文 PDF 已下载并校验。Cas13d 某作者
PDF 端点返回的内容不是 PDF；无效临时文件未进入仓库，论文元数据和公开
链接仍保留。这是参考资料获取失败，不是模型或算法失败。

### 1.5 实验结构与结构 QC

已从 RCSB 官方接口真实下载并校验：

- 6E9F：VI-D ternary；
- 5XWP：VI-A ternary；
- 6E9E：VI-D binary matched state；
- 5XWY：VI-A binary matched state。

权威产物：

- `data/manifests/experimental_structures.yaml`
- `data/manifests/experimental_structure_funnel.json`
- `reports/experimental_structure_data_card.md`

当前真实 QC：

- 4/4 下载结构通过现有 Level 0 坐标/QC 门；
- 6E9F chain A：864 个坐标残基，对应 SEQRES 954，90 个未建模/未映射；
- 5XWP chain A：1125 个坐标残基，对应 SEQRES 1160，35 个未建模/未映射；
- 显式记录 chain break、MSE、RNA chain、crRNA/target chain、hash 和
  SEQRES 映射；
- RNA 只用于接触注释，不会传给标准 ESM-IF1 蛋白 N/CA/C 输入。

进一步的人工文献核查发现：

- 6E9F 为防止切割，把 R295/H300/R849/H854 全部突变为 A；
- 5XWP 的 R1048/H1053 在沉积结构中为 A，而 R472/H477 保留。

因此已建立 `data/manifests/cas13_functional_residues.yaml`，分别记录
沉积构建体氨基酸和文献支持的生物学 R/H 残基。后续设计必须保护后者，
而 construct recovery 与 biological recovery 必须分开报告。

### 1.6 模型后端与当前真实结果

- 统一 inverse-folding schema 和接口；
- 离线本地 checkpoint 的 `EsmIf1Backend`；
- 真正 decode-time hard-fixed 的 `EsmIf1ConstrainedBackend`；
- 固定 token 在自回归到达该位置时输入，不是生成后覆盖；
- 记录每位置 logits、probability、selected token、fixed/free、
  temperature 和 seed；
- 明确保持左到右因果语义，未来固定 token 对早期位置不可见；
- `MsaProfileBackend` 和 `MatchedRandomMutationBackend`；
- 真实 ESM-IF1、ProteinMPNN、LigandMPNN smoke 脚本；
- 6E9F/5XWP ESM-IF1 benchmark 脚本，含 construct/biological sequence、
  HEPN、RNA interface、second shell、buried/surface 和 temperature sweep；
- 报告聚合器，缺失结果会标记 `not_run`，不会填造性能数字。

ESM-IF1、ProteinMPNN 和 LigandMPNN 已取得本文件第 9 节所列真实结果。
完整 temperature sweep 与正式方法矩阵仍未完成，因此不能把已通过的 smoke
扩大解释为完整 benchmark。

## 2. 网络中断和当前失败项

### 2.1 CRISPR-Cas Atlas

下载前已真实核验：

- 官方 URL：
  `https://storage.googleapis.com/crispr-cas-atlas-xy7q13lmk9/crispr-cas-atlas-v1.0.json`
- Content-Length：5,267,508,328 B；
- 磁盘空间充足；
- 顶层为 JSON array；
- 2 MiB schema probe 显示真实字段可由当前 parser 处理；
- Atlas v1.0 的已检查记录没有 direct-repeat orientation 字段。

本次恢复下载时失败：

```text
curl: (6) Could not resolve host: storage.googleapis.com
```

失败发生在 DNS/连接建立阶段，正式目的目录只有 README，没有 `.part`
数据可统计。恢复命令仍为：

```bash
bash scripts/fetch_atlas.sh
```

脚本支持 `.part`、`curl --continue-at -`、大小检查、SHA256 和原子 rename。
网络恢复后已在本节点重新执行该命令；正式 `.part` 正在增长，初始 DNS
失败记录仍保留。只有完整 5,267,508,328 字节下载、SHA256 和原子 rename
完成后才会把 Atlas 状态标为 downloaded。该任务不需要迁移 GPU。

### 2.2 ESM-IF1 隔离环境

`.tools/envs/esm_if1` 的第一次 Conda transaction 暴露了 user-site 泄漏；
该问题已修复。现在从固定 commit 的本地 `third_party/esm` 构建并安装，
所有 bootstrap/smoke 强制 `PYTHONNOUSERSITE=1`，三类 lock 已重新导出，
真实 CPU load/score/sample 已通过。

同一复查中：

```text
torch=2.4.1
torch CUDA runtime=12.1
torch.cuda.is_available()=false
```

同时 `nvidia-smi` 当前被操作系统阻断。真实 CPU ESM smoke 已运行；GPU
smoke 仍必须在 `/dev/dxg` 恢复后单独通过或在 GPU 节点复跑，不能由 CPU
结果替代。

### 2.3 隔离环境与本地工具

四个环境现在均已建立在 `.tools/envs/`，未修改系统 Python：

- `analysis`：Ruff、mypy、pytest、Snakemake 及分析依赖；
- `esm_if1`：PyTorch 2.4.1、CUDA runtime 12.1、PyG 2.6.1 和 pinned
  ESM；
- `ligandmpnn`：pinned NumPy 1.23.5、`dm-tree` 和 MPNN 依赖；
- `bioinformatics`：MMseqs2 18.8cc5c、MAFFT 7.526、HMMER 3.4、
  Infernal 1.1.5、seqkit 2.13.0、Foldseek 10.941cd33、TM-align
  20240303、Git LFS 3.7.1。

每个环境均有 `*-linux-64.explicit.txt`、`*-conda-list.txt` 和
`*-pip-freeze.txt`。真实 PyTorch Geometric GCNConv CPU smoke 已通过，
结果在 `artifacts/system/pytorch_geometric_real_smoke.json`。

bioinformatics 第一次长安装在环境和 locks 已写完后以 shell syntax error
退出。原因是运行中的旧 shell 与此时刚更新的脚本内容不一致；当前脚本
`bash -n` 已通过，并在同一环境上重新完整执行成功。因此这是已恢复的
bootstrap 收尾失败，不是工具安装失败。

## 3. Milestone 状态

| Milestone | 当前状态 | 尚缺验收 |
|---|---|---|
| M0 仓库和长期规范 | complete | 无 |
| M1 环境和依赖 | local CPU complete; GPU pending | 四环境、三类 locks、PyG/ESM/ProteinMPNN/LigandMPNN CPU smoke 已通过；只缺真实 GPU smoke |
| M2 实验结构 benchmark | ESM pilot complete; full matrix pending | 72-candidate ESM temperature/constraint sweep 和正式报告通过；ProteinMPNN/LigandMPNN 尚未进入匹配方法矩阵 |
| M3 Atlas | official download in progress | 完整下载/hash；真实 stream parse；funnel；dedup；MMseqs 六阈值聚类；split audit |
| M4 进化约束 | implementation/fixture only | 真实 subtype MSA、conservation、paired data；MI/APC；结构验证；formal DCA |
| M5 约束生成 | real pilot, matrix pending | 6E9F 真实 constrained ESM 固定位置零违反；尚缺 baseline matrix、novelty/candidate report |
| M6 GPU/refold 迁移 | interface/fixture partial | 完整 bundle；实际 verify；GPU node bootstrap 验证；真实 candidate shards；回收测试 |

## 4. 哪些工作应继续在本节点完成

以下工作不应仅因为存在 GPU 节点就推迟：

1. 诊断 WSL/操作系统的 GPU 阻断；CPU model load 已完成；
2. 完成已恢复的 Atlas 下载并解析；
3. MMseqs2 exact/90/70/50/40/30% 聚类与泄漏审计；
4. subtype-specific MAFFT、conservation/entropy；
5. 高置信 paired-repeat 数据漏斗；如果 orientation 无法恢复，应如实得到
   很小或为零的高置信集合，而不是静默翻转 repeat；
6. 正式报告、bundle export 和 reproducibility verification；
7. 小规模 MPNN 候选生成、新颖度分析和匹配方法矩阵。

## 5. 哪些工作应迁移到 GPU/HPC 资源节点

必须或优先迁移的工作：

- 全 Atlas Cas13 的大规模结构预测；
- 大规模候选回折、多 seed AF2/ColabFold/AF3/Protenix/Boltz；
- 大规模 ESM-IF1/LigandMPNN 采样和完整消融矩阵；
- ESM-IF1 Cas13 domain adaptation；
- 长 Cas13 protein + direct-repeat 的正式 plmDCA/GREMLIN/CCMpred；
- Foldseek 全结构聚类；
- 本机 8 GiB 显存无法容纳的 batch 或多状态联合工作。

无需迁移但可在节点复跑以确认环境一致性的工作：

- 单结构 ESM-IF1/ProteinMPNN/LigandMPNN smoke；
- bundle verification；
- 少量候选 refold ingest；
- US-align/TM-align 结构比较。

## 6. GPU 迁移指南完成度

已存在：

- `docs/GPU_MIGRATION.md`
- `scripts/export_gpu_bundle.sh`
- `scripts/verify_gpu_bundle.sh`
- `scripts/bootstrap_gpu_node.sh`
- `scripts/sync_assets.sh`
- `scripts/launch_gpu_tmux.sh`
- provider-neutral refold FASTA/JSONL export、deterministic sharding、
  expected-output schema、retry manifest、pLDDT/PAE/structure ingest、
  missing-output audit；
- mock prediction E2E，`is_mock=true`。

中断恢复后，迁移路径已进一步补强：

- bundle 现在包含 git commit、dirty 状态、配置 hash、环境 locks、
  containers、运行脚本、clone instructions、输入 shard 清单和 expected
  output schema；
- 大资产不会嵌入 bundle；`ASSET_SHA256SUMS` 记录逐资产大小和 SHA256；
- `sync_assets.sh` 在 rsync 续传后对目标资产逐一验 hash；
- `verify_gpu_bundle.sh` 同时验证 bundle 内部和可选目标资产根目录；
- GPU bootstrap 覆盖 analysis、ESM-IF1、LigandMPNN、bioinformatics 四环境；
- tmux launcher 记录 git、GPU、配置、起止时间、日志、退出码以及
  SUCCESS/FAILED。

本地已经达到的验收：

- 已从 clean commit `7dc0491d84419771b4b5e14d2f8daf39f36e68d1`
  导出
  `artifacts/bundles/gpu-bundle-7dc0491d8441-6ad46d8577/`；
- `git_worktree_dirty_at_export=false`；
- bundle 内所有文件的 `SHA256SUMS` 已通过；
- 5 个模型 checkpoint 和 8 个实验结构 PDB/mmCIF 已按
  `ASSET_SHA256SUMS` 在本地资产根目录逐一通过；
- `make verify-reproducibility` 已通过 shell parsing、manifest、
  Ruff、format、strict mypy、42 个 pytest 和 bundle verification。

仍未达到的目标节点验收：

- `bootstrap_gpu_node.sh` 尚未在目标 GPU 节点验证完全隔离安装和 GPU
  imports；
- `sync_assets.sh` 的逐资产 hash 路径尚未在第二台机器完成真实传输验收；
- 尚无真实 candidate shards；
- 尚未验证 tmux launcher 对实际 benchmark 的退出码、日志和失败现场；
- 正式 DCA job export/ingest 未完成；
- 大规模结构预测 provider 只有通用交换接口，没有在目标软件上真实跑通。

当前基础 bundle 的 `missing_assets` 只有尚在下载的 Atlas 正式 JSON；
脚本正确地将其列为缺失而没有嵌入未完成 `.part`。因此答案是：迁移操作
指南、clean bundle 导出、本地内部校验和现有资产校验已完成；第二台机器的
真实传输、环境恢复和 GPU 运行仍未完成，不能称为目标 GPU 迁移验收完成。

## 7. 从中断点恢复的执行顺序

1. 完成已恢复的 Atlas 下载；
2. 真实 Atlas process → cluster → MSA → conservation；
3. 导出并验证 GPU bundle；
4. 更新 ExecPlan、STATUS、DECISIONS，按已验证功能分小 commit；
5. 全部质量门通过后 push；
6. 在 GPU 节点执行三模型 GPU smoke、正式大规模采样/回折和 DCA。

## 8. 证据边界

截至当前更新：

- Level 0：仓库、I/O、fixture workflow、真实结构下载/QC 已有证据；
- Level 1：尚无针对完整 Atlas 的正式候选新颖性结果；
- Level 2：已有真实 ESM-IF1 score/constrained sample、72-candidate ESM
  pilot benchmark、ProteinMPNN 和 RNA-context LigandMPNN smoke；完整
  匹配新颖度的多方法 benchmark 尚未完成；
- Level 3：尚无多模型/回折支持结果；
- Level 4：不在本项目当前范围内，且没有 wet-lab 结果。

## 9. 中断后继续推进的更新

本审计写入后，已从中断点继续执行并取得以下新结果：

- 修复了 ESM 环境的 user-site 泄漏；
- `fair-esm` 现在从固定 ESM commit 的本地源码构建，安装路径位于
  `.tools/envs/esm_if1`；
- bootstrap 和 smoke 脚本强制 `PYTHONNOUSERSITE=1`；
- 采用项目 backend 中与固定上游实现等价的 ESM tensor/coordinate batch
  适配，结构 I/O 使用本项目已测试的 PDB/mmCIF parser，未修改
  `third_party/esm`；
- 真实 CPU ESM-IF1 toy score/sample 通过；
- 全固定 toy constrained sampling 精确恢复，固定位置违反数为 0；
- 6E9F chain A（864 aa）真实 score：
  conditional log-likelihood `-1813.9333906933316`，
  perplexity `8.161760905453558`；
- 5XWP chain A（1125 aa）真实 score：
  conditional log-likelihood `-2553.345790145657`，
  perplexity `9.675923652828677`；
- 结果写入 `artifacts/system/esm_if1_real_smoke.json`，`is_mock=false`，
  最大证据为 Level 2 inverse-folding compatibility；
- 推理设备为 CPU。GPU smoke 仍因当前 WSL 缺少 `/dev/dxg` 而
  `not_run`。

这些结果更新了本文件前部“尚无 Level 2”的时间点状态：现在已有真实
ESM-IF1 Level 2 score 证据，但仍没有 Level 1 候选新颖性、Level 3
多模型/回折证据或 Level 4 wet-lab 证据。

ProteinMPNN 后续真实 CPU smoke 也已通过：

- checkpoint：
  `models/proteinmpnn/v_48_020.pt`，SHA256
  `c9cb4a671d79604111231f8dbfc7c590e06f1197453b7e0eada8b796d5a0eae9fd4`；
- 6E9F A 链严格坐标序列为 864 aa；
- 上游 ProteinMPNN residue-number tensor 为 893 位，其中 29 个内部
  缺坐标槽位为 `X` 并被 mask；
- 去除 `X` 后与严格坐标序列一致；
- 逐位置概率已保存；
- CPU 耗时约 9.90 秒，`is_mock=false`；
- 结果：`artifacts/system/proteinmpnn_real_smoke.json`。

第一次 ProteinMPNN smoke 失败是验证器错误：它要求上游 893 位编号张量
必须等于 864 位严格坐标序列。修正后不删除缺失槽位，而是同时审计两种
长度并验证 mask，真实模型运行随即通过。

LigandMPNN 真实 CPU smoke 随后通过：

- checkpoint：
  `models/ligandmpnn/ligandmpnn_v_32_010_25.pt`，SHA256
  `161cd264061fda9680cbb940255522ae42f2966c552d045d87913d9452a80970`；
- 设计链为 6E9F A（864 aa），上下文链为 RNA B/C；
- 上游 parser 保留 1680 个 B/C 非蛋白原子，残基/离子类型包括
  A、C、G、U、Mg；
- 固定残基 A58 保持；
- per-residue statistics 和包含 RNA context 的 backbone 输出均存在；
- CPU 耗时约 9.80 秒，`is_mock=false`；
- 结果：`artifacts/system/ligandmpnn_real_smoke.json`。

前两次 LigandMPNN smoke 均在推理前失败，原因依次为缺失 `dm-tree` 和
NumPy 1.26 已移除上游代码使用的 `np.int`。已依据 pinned 上游
requirements 将隔离环境补入 `dm-tree` 并固定 NumPy 1.23.5，重新导出
locks；没有修改 `third_party/LigandMPNN` 源码。

实验结构 ESM-IF1 pilot benchmark 已完成：

- 权威成功 run：
  `results/runs/20260731-benchmark-experimental-a998ff40aa-7599be0-r003`；
- 2 个 scaffold（6E9F、5XWP）× 3 个方法条件 × 3 个 temperature ×
  4 个 seed/sample，共 72 个真实 CPU 候选；
- 72/72 个 `candidate_id` 唯一，最大 ID 重数为 1；
- 所有方法合计固定位置违反数为 0；
- wall time 约 1409.0 秒，设备为 CPU，`is_mock=false`、Level 2；
- 输出包括每候选 trace、每位置 native score、区域注释、Markdown 和
  HTML benchmark；
- 项目汇总报告：
  `results/runs/20260731-benchmark-experimental-bcfd0be469-7599be0/report/`。

历史失败/被替代 run 均保留：

- `r001` 的模型计算本身成功，但旧 candidate ID 没有编码 temperature 和
  constraint 条件，72 行只有 16 个唯一 ID，因此只保留作审计，不用于
  后续 refold manifest；
- `r002` 在模型加载前因 wrapper 使用相对 config 路径失败；
- `r003` 修正 config 规范化、不可覆盖 run ID 和 sampling-condition
  digest 后通过。

当前 raw recovery 是不同固定位置比例下的 pilot 描述值，不能直接解释为
某方法优于另一方法；正式论文比较仍必须匹配设计位置和新颖度。

其他恢复后通过的本地检查：

- 真实 PyTorch Geometric GCNConv CPU smoke；
- fixture paired-MSA MI/APC、20 次 bootstrap、20 次 permutation smoke；
- MI/APC 输出明确标注 `is_mock=true`，正式 DCA 为 `not_run`；
- mock preflight 固定位置违反为 0；
- Ruff、format、strict mypy、42 个 pytest 和 70.08% branch coverage。
