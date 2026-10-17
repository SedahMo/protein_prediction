# 抗体工程化蛋白稳定性分析

将 **序列/结构预检 → FoldX RepairPDB → 配对点突变 BuildModel → 重复采样 → 可追溯报告** 连成可复现流程。输入校验、解析、可视化及测试仅依赖 Python 标准库；真实能量计算接入用户持有许可的 FoldX。

**当前交付的是分析流程，尚未计算两条实际工程化蛋白。** 内置数值与小型 PDB 均为软件测试用合成数据。

## 先选择比较类型

| 比较问题 | 处理方式 | 能支持的结论 |
| --- | --- | --- |
| 同一构型上的点替换/组合替换 | 同一参考结构；明确链、PDB 残基号和 WT；BuildModel 生成匹配 WT | 同一研究内 ΔΔG = G(mutant) − G(WT)，正值表示模型预测更不稳定 |
| 不同 linker，尤其长度不同 | 独立结构集合与连接区/接口分析模板 | 结合模型质量和实验判断；不直接按总能量排序 |
| VH–VL ↔ VL–VH | 独立建模各连接方式，检查接口及 linker 上下文 | 比较构型、表达、单体比例与热转变等证据 |
| 不同蛋白或独立预测结构 | 保留来源、残基覆盖和模型置信度 | 不生成跨蛋白 FoldX 总能量排行榜 |

同长度 linker 的替换虽然可以写成点突变，也须先判断骨架是否仍可共用。完整依据见 [方法学](docs/methodology.md)。

## 不安装 FoldX，先跑通

在仓库根目录，使用 Python ≥3.10：

```bash
python3 -m protein_stability doctor
python3 -m unittest discover -s tests -v
python3 -m protein_stability validate-fasta examples/toy_reference.fasta
python3 -m protein_stability plan examples/toy_plan.json --out runs/toy-plan
python3 -m protein_stability demo --out runs/demo
```

打开 `runs/demo/report.html`。报告含人工数据标识、均值±样本 SD 图、逐次记录、四类比较边界和数据来源。重新解析示例：

```bash
python3 -m protein_stability report runs/demo/synthetic_run --out runs/demo-reparsed
```

错误输入或不完整结果会产生非零退出状态。规划/执行目录必须是新目录。预生成文件：[演示报告](examples/demo_report.html)、[待计算模板](examples/report_template.html)。

## 接入真实输入与 FoldX

1. 将成熟链 FASTA 和完整参考 PDB 放入 `private_inputs/`。FASTA 第一段标题为记录 ID；序列接受 20 种标准氨基酸，小写会归一化。
2. 复制 [点突变配置](examples/point_mutation.json)，填写实际文件、链映射、替换清单和软件版本。路径相对于配置文件解析。配置中的 A1V 是格式占位符，需替换。
3. 查看预检和规划结果；完成结构质量、二硫键、侧链、质子化及预测置信度检查。
4. 按 [FoldX 环境说明](docs/foldx_setup.md) 提供许可覆盖当前机器的软件与资源，再执行。

```bash
python3 -m protein_stability inspect-pdb private_inputs/reference.pdb --output runs/pdb-inspection.json
python3 -m protein_stability plan examples/point_mutation.json --out runs/real-plan
export FOLDX_BINARY=/absolute/private/path/foldx
python3 -m protein_stability run examples/point_mutation.json --out runs/real-run --license-confirmed
python3 -m protein_stability report runs/real-run --out runs/real-report
```

`--license-confirmed` 是使用者对授权范围的声明，不会下载或认证许可。`doctor` 仅探测路径，不执行 FoldX。可独立准备 PDB，链映射支持多链，如 `{"H":"heavy","L":"light"}`：

```bash
python3 -m protein_stability prepare-pdb private_inputs/reference.pdb --fasta private_inputs/reference.fasta --chain-map '{"A":"reference"}' --output runs/prepared.pdb --audit runs/prepared-audit.json
```

预处理保留坐标、链 ID 和 PDB 残基号，要求选择链与 FASTA 完全一致。缺失主链原子、多模型、插入码或链片段歧义会阻止执行。选择 blank/A alternate location；其他情况须明确处理。水和其他 HETATM 会移除并记录；配体/修饰影响结果时应另建适用的准备方案。预检不补环、补原子或预测结构，侧链完整性和几何质量另行核实。

## 结果和追溯

- RepairPDB 执行一次；每次 BuildModel 使用同一修复参考结构，独立目录、`numberOfRuns=1`，默认重复 5 次。
- 只解析 `Dif_*.fxout` / `Diff_*.fxout`；不将 `Raw_`、`Average_` 绝对值当作 ΔΔG。
- 核对变体行号/模型名、全部重复、有限数值及文件哈希；不完整研究不输出部分排名。
- `plan.json`、`run_manifest.json` 保存比较类型、映射、输入/软件/产物哈希、参数、命令和状态；保留原始 stdout/stderr。
- 报告输出 `results.json`、`summary.csv`、自包含 `report.html`；图按配置顺序展示。

重新生成真实报告时需保留完整运行目录（计划、输入快照、修复结构、逐次原始产物及日志）。可以整体移动或归档该目录；解析不依赖原始输入和软件所在路径。只保留 manifest 和差值表不足以通过真实运行的完整性检查。

配置中可选 `foldx.resources_dir` 指向授权资源目录；规划会记录资源哈希，执行时复制至各工作目录。原始资源不改动，软件自动生成的资源也随原始运行文件保留。

SD 描述算法重复的取样变异，不能当作实验误差、准确度或结构集合不确定度；单次运行 SD 为空。流程不将 ΔΔG 转成 Tm、聚集倾向或亲和力。版本、哈希和原始产物支持审计与重新解析；随机搜索和软件差异可能使重新计算的数值不同。

## 其它构型的报告模板

```bash
python3 -m protein_stability comparison-template --type linker --study-id linker-study --out runs/linker-template
python3 -m protein_stability comparison-template --type domain_order --out runs/orientation-template
python3 -m protein_stability comparison-template --type cross_protein --out runs/cross-protein-template
```

[结构比较记录模板](examples/structural_comparison.json) 收集模型来源、域边界、置信度、接口及实验结果，不进入 BuildModel。上游结构预测使用用户选定工具，模型、版本/种子及置信度经核查后进入流程；仓库不自动调用付费服务。

## 环境与开发

无需安装包即可在根目录运行；也可用 `python3 -m pip install --no-build-isolation --no-deps -e .` 安装命令 `protein-stability`，需要本机已有 setuptools 构建工具。

```bash
docker build -t antibody-stability .
docker run --rm antibody-stability
```

镜像不包含 FoldX。默认基础镜像为 `python:3.12-slim`；长期固定环境时以 `--build-arg PYTHON_IMAGE=python@sha256:<已核验摘要>` 固定镜像。GitHub Actions 配置覆盖 Python 3.10 / 3.12 / 3.13 的模块测试与演示。

`protein_stability/` 核心模块；`tests/` unittest；`examples/` 合成夹具和模板；`docs/` 方法与环境说明；`runs/` 本地结果；`private_inputs/`、`foldx_private/` 私有输入与授权资源。后三者已被 Git 忽略。

当前 Codex 环境实测见 [environment_check.json](docs/environment_check.json) 和 [validation.md](docs/validation.md)。
