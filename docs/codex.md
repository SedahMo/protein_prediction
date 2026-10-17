# Codex 云端执行说明

检查日期：2026-10-17 UTC。当前会话环境：Linux、Python 3.12.14、git 2.51.1；PATH 中未发现 foldx。无第三方库的验证、计划、解析和报告模块已运行。并未获取 FoldX、验证许可、预测三维结构或计算真实 ΔΔG。

可用于后续 Codex 任务的指令：

> 先阅读 README 和 docs，运行 doctor 与 unittest。使用我提供的完整 FASTA、单模型 PDB、structure-review.md 和点突变配置。验证链序列与残基编号，先运行 plan。只有在已提供合法授权的 FoldX 可执行文件及其所需 ancillary 文件后执行 run。不得下载未经授权的 FoldX，不要上传授权文件。保留 manifest、help 输出、全部日志和原始 fxout。若结构含缺失残基、插入码、HETATM 或 alternate conformer，停止自动处理并输出具体修复需求。不得以绝对能量比较 linker、VH/VL 顺序或不同骨架。把逐运行 ΔΔG、均值/SD/范围和结构来源整理到报告，并明确真实计算和合成夹具。

建议生产任务固定 Python/FoldX 发行版本和计算环境镜像 digest，并备份输入哈希与模型原始文件。当前标准库路线没有 pip 运行依赖，因此 PYTHONPATH 用法无需网络安装。FoldX 许可或额外结构预测引擎的可用性要在目标云环境逐项验证。
