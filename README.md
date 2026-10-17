# protein_prediction

可复现的抗体工程化蛋白稳定性分析流程。运行依赖仅 Python 标准库（Python ≥3.10），FoldX 为可选、用户自行授权的外部程序。

## 快速验证（不需要 FoldX）

```bash
export PYTHONPATH="$PWD/src"
python -m unittest discover -s tests -v
python -m protein_prediction doctor
python -m protein_prediction validate-fasta examples/toy.fasta
python -m protein_prediction inspect-pdb examples/toy.pdb
python -m protein_prediction plan examples/config.json --out runs/demo
python -m protein_prediction report examples/Dif_synthetic.fxout --manifest runs/demo/manifest.json --out runs/demo/report
```

打开 `runs/demo/report/index.html`，同时得到 CSV、JSON。仓库内 `examples/report/` 提供同一报告模板预览。示例为两残基、非物理几何的测试夹具，能量完全合成；不能用于 FoldX 科学计算。所有输出目录首次创建，已有运行目录拒绝覆盖。

也可在联网环境 `python -m pip install -e .` 后使用 `protein-stability` 命令。

## 真实计算

1. 准备完整参考 FASTA、经过人工审查的单模型 PDB、点替换列表及结构来源记录。见 [结构准备](docs/structure-preparation.md)。
2. 复制 `examples/config.json`，替换 pdb/fasta 路径、链映射、variants；去掉 synthetic 标签，改用 `data_kind: "computed"`。路径相对于配置文件。
3. 安装合法授权的 FoldX，设置 `FOLDX_BIN` 为可执行文件路径。若版本需要 rotabase，配置 `rotabase` 为本地文件路径。不要提交二进制、授权文件或 rotabase。
4. 先运行 `plan` 检查输入；再用新目录执行：

```bash
python -m protein_prediction run config/project.json --out runs/project-001 --foldx /absolute/path/to/foldx
```

执行顺序：检查输入 → 保留原始编号的结构复制 → RepairPDB → 检查修复后的残基映射 → BuildModel（numberOfRuns）→ 配对 Dif 结果解析 → HTML/CSV/JSON。失败时记录 failed 状态与错误，保留日志。默认每个命令超时 3600 秒，可配置 `timeout_seconds`。

`manifest.json` 保存输入与二进制 SHA256、完整命令、环境、配置和顶层输出校验和；`foldx-help.txt` 保存安装版本的帮助输出（另请在项目记录中写明准确版本/发行日期）。相同输入与版本使过程可追踪；FoldX 随机优化不承诺逐字节复现。不要把算法重复 SD 当成预测置信区间。

## 比较边界

| 设计变化 | 当前实现 | 如何判读 |
|---|---|---|
| 同一结构中的点替换/多点联合替换 | 完整点突变执行接口 | 使用 BuildModel 配对 WT 的 ΔΔG |
| 不同 linker 长度/序列/构型 | 阻止进入点突变路线 | 独立建模、多构象及连接区几何评估 |
| VH–VL 改 VL–VH | 阻止进入点突变路线 | 独立建模，比较域完整性、界面与聚集倾向 |
| 不同蛋白/骨架/长度 | 阻止混合排序 | 不将绝对 FoldX 能量或能量/残基作为可靠稳定性名次 |

详见 [比较设计](docs/comparisons.md)、[Codex 云环境](docs/codex.md)。当前不包含结构预测引擎、分子动力学、结合自由能或实验 Tm 预测。此前讨论的两条真实序列未从聊天片段重建，用户应提供完整、确认过的 FASTA 与结构。

## FoldX 输出支持

只接受制表符分隔的 `Dif_*.fxout`，按 `Pdb` 和 `total energy` 表头解析。支持 `<stem>_<variant-index>_<run-index>.pdb`（末尾 `.pdb` 可省略）的命名；显式检查变体和重复覆盖，不按行顺序猜测身份。其他版本命名失败即停止；需先根据该版本 PdbList/Raw 文件确认配对规则再扩展解析器。Average/Raw 文件不能替代逐运行 Dif 数据。

ΔΔG = G(mutant) − G(WT)，kcal/mol；正值预测去稳定化。输出均值、样本 SD、最小值、最大值及运行次数，不设通用的实验成败阈值。

## 官方方法来源

- [FoldX BuildModel](https://foldxsuite.crg.eu/command/BuildModel)：局部配对 WT、重复计算、Dif/Average/Raw 输出。
- [FoldX RepairPDB](https://foldxsuite.crg.eu/command/RepairPDB)：结构修复的用途和输出。
- [individual_list 格式](https://foldxsuite.crg.eu/parameter/mutant-file)：WT 氨基酸 + 链 + PDB 残基编号 + mutant；多点以逗号分隔并以分号结束。

## 开发验证

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

GitHub Actions 对 Python 3.10 与 3.12 执行无商业软件测试，并生成合成示例报告。模拟 FoldX 测试只验证程序调用契约；真实 FoldX 端到端计算尚待有授权的环境验证。
