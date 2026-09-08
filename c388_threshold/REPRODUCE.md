# 复现顺序与过程记录

日期：2026-09-08。所有命令在仓库根目录运行，建议独立 Python 3.12 环境。

```bash
python -m pip install -r c388_threshold/requirements.txt
# 1. 从已保存的构建级特征、实验标签和冻结旧模型开始
python c388_analysis/analyze.py
python c388_analysis/summarize.py
python c388_analysis/finish.py
# 2. 固定 C388 50%，扫描模型分数阈值并做训练内阈值选择
python c388_threshold/thresholds.py
# 3. 重新定义标签并重训，联合扫描 C388 与分数阈值
python c388_threshold/joint.py
# 4. 冻结旧模型的阈值敏感性；不作为新模型推荐依据
python c388_threshold/frozen_joint.py
# 5. 固定候选 C388 分界，验证训练折内选择分类阈值
python c388_threshold/validate_joint.py
```

命令会重写结果文件。联合扫描不额外调 RF 参数。28 个不等价标签划分 × 4 个固定模型 × 2 种留出验证；训练集中出现单类的组合跳过，并未伪造预测。

先完成 C388 50% 的冻结模型迁移、固定模型重训和嵌套模型选择，再按用户要求进行分类阈值扫描，最后扩展为 C388 标签阈值和模型分数阈值的联合探索。候选 30%、35%、40%、50%、55% 的补充验证仍受事后挑选标签阈值的影响，不是完全独立的新批次验证。

代码归档只调整 finish.py 的历史参考路径，使其指向本仓库 architecture_validation/ 和 combined12/；没有修改模型参数、标签、数据或既有计算结果。模型 joblib 是 C388 50% 的最终重拟合模型，不能直接冒充 C388 40% 重训模型。40% 的结果为交叉验证输出，尚未发布其部署模型。

统计/机器学习复现所需特征缓存、模型、实验聚合数据、分组、全部结果表及图已包含。原始大型 AF3 ZIP/XLSX 和 units/ 不重复入库；原始输入文件名及 SHA256 见 c388_analysis/raw_input_manifest.json。construct_arrays.npz 为可由原始 AF3 重算的矩阵中间缓存，当前统计/模型/阈值脚本不读取它，因此未重复上传。

如需从原始 ZIP 重提特征，按 c388_analysis/README.txt 放置 upload/ 输入后运行 inventory.py 与 extract_features.py；严格的 mmCIF 解析仅针对本次上传格式。

固定 C388 50% 报告中“未搜索更有利阈值”是该阶段的历史描述；后续探索请以 c388_threshold/README.md 为准。mutation-group-out 是 TRM 突变位置组合留出，不能称为 architecture-out。
