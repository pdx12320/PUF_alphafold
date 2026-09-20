# 骨架分析输入缓存

这些文件用于 `architecture_validation/extract.py`、`prepare.py`、`combined12/train.py` 及第二批 repeat 映射。数据从清理前快照逐字节迁入，来源和 SHA256 见 provenance.json。保留的是序列、特征、矩阵和残基映射；旧阶段的模型、报告与图表已删除。

历史假设、验证结果、失败原因与预处理质控统一见 [Dry Lab DBTL](../../docs/DRY_LAB_DBTL.md)。当前统计/模型重现直接使用 architecture_validation 已合并缓存，无需先运行旧分析。
