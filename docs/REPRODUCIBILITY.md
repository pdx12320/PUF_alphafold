# 复现说明

所有命令从仓库根目录运行，建议Python 3.12及requirements.txt所列依赖。2026-09-20整理仅迁移必要输入与通用函数，未重训或更改现有主结果。退出当前目录的五轮基线、探索和质控见 [DBTL](DRY_LAB_DBTL.md)。

## 当前骨架结果

```bash
python tools/validate_snapshot.py
python tools/reproduce_current.py
```

后者按stats.py → models.py → plot_static.py → finish.py → write_report.py运行，读取architecture_validation的main_features/all_features/all_arrays/all_residue_mapping和combined12/manifest.json。无需先运行已删除阶段。运行会覆盖对应结果，建议独立checkout。

`python combined12/train.py`可从 `data/architecture_inputs/batch1_features.csv` 与 `new_batch/final_summary.csv` 重新训练24种PUF12模型。`architecture_validation/prepare.py`所需第一批矩阵和映射均在 `data/architecture_inputs/`；迁移前后SHA256见该目录provenance.json。

## 当前分类结果

- C388动态阈值：[独立输入与命令](../c388_threshold/REPRODUCE.md)。
- C871/C295：[22个固定入组突变体的报告](../trm22_offtarget/REPORT.md)与[复现命令](../trm22_offtarget/README.md)。
- C388 local/nonlocal：[快照核验与原包恢复](../c388_local_nonlocal_optimization/README.md)。仓库仅含摘录快照，完整重跑需要原始94文件ZIP；缺少的PP分析仍未完成。

## 原始AF3输入和保留范围

当前骨架原始提取入口为 `architecture_validation/extract.py`，共享CP函数为 `tools/cp_features.py`。需要upload目录的两份2026-09-07 ZIP及91–96 ZIP；第一批使用已验证的输入缓存。大型ZIP不随本仓库提供。

`new_batch/map_repeats.py`使用迁移后的repeat模板。第二批历史冻结评分依赖已删除权重，其analyze.py、finish.py和旧RF副本同步退出；第二批已保存结果仍供合并分析读取，`report_plot.py`可由保留的评分和输入缓存生成对应第二批图表。需要完整重跑已退出的历史流程时，使用DBTL链接的清理前Git版本。

原始文件清单仍见raw_input_inventory.json，其中曾记录的截断ZIP不可用于提取。先核验原始ZIP完整性，再运行耗时流程。此前N端范围、重复WT与RNA/seed独立性问题的诊断集中于DBTL Cycle 4。

`trm22_offtarget/` 使用自包含输入缓存，`SOURCE_MANIFEST.json` 记录输入在固定Git提交中的路径与SHA256，`sources/` 保留ContactSeek来源代码和许可证。旧的30突变体C871/C295目录不再作为当前执行入口。

## 快照核验

`tools/validate_snapshot.py`校验当前文件清单SHA256、24个独立PUF12和4个成功标签、预测覆盖及等价分组。重跑浮点和序列化结果可受环境影响，快照SHA256用于已发布文件核验。

预测分数未校准。原始model/seed为技术预测重复，泛化性能来自构建级折外预测；全数据重拟合训练评分不能替代验证。模型文件仅从可信来源加载。
