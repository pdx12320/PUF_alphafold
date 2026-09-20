# 联合动态阈值复现

2026-09-20：原固定C388分析目录已删除。所需构建特征、实验率和分组逐字节迁入本目录inputs/；迁移校验见inputs/provenance.json。此处不依赖旧权重。

在仓库根目录运行：

```bash
python -m pip install -r c388_threshold/requirements.txt
python c388_threshold/thresholds.py
python c388_threshold/joint.py
python c388_threshold/validate_joint.py
```

thresholds.py保留固定50%实验标签下训练内选择模型分数阈值的动态流程。joint.py扫描实验标签与模型分数的联合阈值；validate_joint.py对30/35/40/50/55%候选进行训练内分数阈值验证。命令覆盖相应结果。

28种不等价标签划分、4个固定模型、两种验证的既有结果不变。选择标签和模型仍有全数据探索偏差，位置组合留出可能共享组成位置。40%候选尚无部署权重或独立新批次验证。

旧冻结模型迁移、固定分数扫描及单独固定50%模型的报告和权重已退出，历史结果见 [DBTL Cycle 5](../docs/DRY_LAB_DBTL.md)。源文件清单在inputs/raw_input_manifest.json；从大型AF3原始ZIP重新提取旧阶段特征，需要DBTL链接的清理前版本。当前缓存足以执行上述动态流程。
