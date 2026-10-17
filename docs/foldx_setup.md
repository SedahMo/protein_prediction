# FoldX 环境与执行接口

本仓库的输入检查、解析、测试和 HTML 报告可以在没有 FoldX 的环境运行。真实 RepairPDB / BuildModel 需要用户自行提供合法可用的 FoldX 可执行文件、合适的结构和明确的替换清单；仓库不分发软件、许可或私有蛋白序列。

官方资料核对日期：2026-10-17 UTC。下面给出的是经官网核对的 FoldX 命令接口，不代表当前环境已经安装或执行过 FoldX。仓库封装命令与配置以 [README](../README.md) 和命令的 `--help` 为准。

## 许可与云端运行

FoldX 官方提供学术及非营利研究机构研究用途许可；公司用途应取得适用的商业许可。官方将商业许可区分为限定机器/地点的 On-Site 与可覆盖授权用户异地数据中心运行的 Global / Multi-Site。将二进制上传到云端前，确认实际许可覆盖该执行环境。两周评估许可仅供评估，不能据此将输出用于研发。[1]

许可协议、软件二进制与 rotabase 资源不应提交到 Git，也不放入公共镜像。通过环境变量或私有目录提供路径；输入和运行结果默认放入已忽略的本地目录。检查日志、报告与 git diff 后再决定是否共享真实项目资料。

## 环境检查

在仓库目录运行：

```bash
python -m protein_stability doctor
```

环境检查用于说明 Python、操作系统、架构和 FoldX 路径等是否可用。它不能认证许可，也不能证明软件在真实结构上已经完成计算。对于用户提供的二进制，另外核对系统架构、执行权限及动态库依赖，并保存实际软件版本与文件哈希。

例如，设置私有二进制路径后人工检查本机支持的参数：

```bash
export FOLDX_BINARY=/absolute/private/path/foldx
"$FOLDX_BINARY" --help
"$FOLDX_BINARY" --help --command=BuildModel
```

官网页面及本地二进制可能对应不同版本。本地 `--help` 和该版本发行说明决定可用参数；不要将旧版本名称或资源要求硬编码为所有版本都适用。

官网手册说明：4.x 需要工作目录中的 `rotabase.txt` 或软链接；5.x 可自动生成。[2] 因此“未提供 rotabase”不一定意味着 5.x 不能运行。使用与实际版本配套的资源，记录资源来源/哈希，并保存自动生成的资源以便追溯。

## 最小真实运行：独立工作目录

以下命令演示 FoldX 本体语法。先在一个新的、可写的运行目录准备 `reference.pdb` 与 `individual_list.txt`；`FOLDX_BINARY` 为绝对路径。旧结果不复用，避免输出覆盖或结果混淆。

`individual_list.txt` 使用 PDB 的链 ID 与残基号，例如：

```text
FA39L;
LA41F,KA42I;
```

第一行是 A 链 39 位 Phe→Leu，第二行是同一候选上的两个替换。示例仅用于说明格式，不对应本仓库虚构演示结构。官方要求文件名以 `individual_list` 开头、突变之间以逗号连接、每行分号结束且无空格。[3] 在执行前必须确认该位置确有预期的 WT 氨基酸；本仓库不将插入/缺失或结构域交换转换成此格式。

```bash
"$FOLDX_BINARY" --command=RepairPDB --pdb=reference.pdb --pH=7 --temperature=298
"$FOLDX_BINARY" --command=BuildModel --pdb=reference_Repair.pdb --mutant-file=individual_list.txt --numberOfRuns=5 --pH=7 --temperature=298 --out-pdb=true
```

温度单位为 K。官网列出的默认值为 pH 7、298 K；这里显式传入，以便记录并保持两个阶段一致。[4,5] 对具体配方还需考虑实际离子强度、辅因子、二硫键及其他设置。改变条件时应建立新运行，完整保存参数。

## 预期输出与解析

| 阶段 | 默认输出模式 | 用途 |
| --- | --- | --- |
| RepairPDB | `<stem>_Repair.pdb` | 修复后的参考坐标 |
| RepairPDB | `<stem>_Repair.fxout` | 修复残基的能量诊断 |
| BuildModel | `Dif_<stem>_BM.fxout` | 突变体与对应 WT 的能差；主解析入口 |
| BuildModel | `Average_<stem>_BM.fxout` | 多次运行的汇总 |
| BuildModel | `Raw_<stem>_BM.fxout` | WT 和突变体的逐项原始能量 |
| BuildModel | `PdbList_<stem>_BM.fxout` | 突变体及相应 WT 的结构名称映射 |

其中 BuildModel 的 `<stem>` 是其输入 PDB 的文件干名，可能已包含 `_Repair`。官网的 `output-file` 参数还会影响标签；解析应依据文件类型和表头，而不是只认某一个完整文件名。官方差值前缀为 `Dif_`，不是 `Diff_`。[6,7]

保留原始标题/列名；`Total Energy` 列在 `Raw_` 中是绝对评分，在 `Dif_` 中才是匹配能差。错误文件、空文件、非数值和不完整重复应报错，不能默认为 0。只含 `Average_` 的输出不足以重建逐次运行的变异。

如为每个重复单独调用 BuildModel，可将 `numberOfRuns=1`，每次放入独立目录。若一次调用 `numberOfRuns=5`，应检查真实产物中逐次结果和 WT 对应关系；不能把一个平均行重复记为 5 次。增加重复用于观察搜索收敛，不给出真实实验误差范围。

## 真实项目需要补充的输入

| 输入 | 需要说明的内容 |
| --- | --- |
| `FASTA` | 研究实体、结构域和 linker 边界；成熟链与标签定义 |
| 可信 `PDB` | 实验来源或预测工具/版本/种子、链与残基映射、缺失片段、寡聚状态 |
| 置信度/结构质控 | 预测模型的 pLDDT/PAE 及来源；几何、二硫键和 VH/VL 接口检查 |
| 替换清单 | WT/链/PDB 残基号/目标残基；组合突变按候选分组 |
| 软件与资源 | 可用二进制、版本、许可覆盖环境、必要的配套资源 |
| 执行参数 | pH、K 温度、重复数、比较类型及输入文件哈希 |

没有上述结构与软件时，先运行演示、输入预检和报告模块；结果应标记为演示或待计算。结构预测属于上游接口，当前流程不声称仅凭 FASTA 已经求得真实折叠 ΔΔG。不同 linker、VH/VL 顺序或独立蛋白比较的处理见 [方法学](methodology.md)。

## 官方资料

1. [Licensing and Services](https://foldxsuite.crg.eu/licensing-and-services)。
2. [FoldX Manual](https://foldxsuite.crg.eu/documentation)：命令帮助与 rotabase 版本差异。
3. [mutant-file](https://foldxsuite.crg.eu/parameter/mutant-file)：individual list 格式。
4. [pH](https://foldxsuite.crg.eu/parameter/pH)。
5. [temperature](https://foldxsuite.crg.eu/parameter/temperature)。
6. [RepairPDB](https://foldxsuite.crg.eu/command/RepairPDB)。
7. [BuildModel](https://foldxsuite.crg.eu/command/BuildModel)。
