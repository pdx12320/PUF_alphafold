# 相对WT动态分类重分析

## 本轮结果

| 位点 | WT编辑率 | Δ候选分界(pp) | 留一BA | 位置组合分组BA |
|---|---:|---:|---:|---:|
| C388 | 72.561% | −45 | 76.0% | 78.0% |
| C871 | 67.5015% | −25 | 78.5% | 79.4% |
| C295 | 32.958% | −20 | 82.5% | 52.5% |

[完整过程](PROCESS.md) · [结果汇总](results/wt_relative_summary.csv) · [逐突变体预测](results/wt_relative_per_construct_predictions.csv) · [HTML报告](WT_relative_classification_report.html) · [核验记录](results/validation_checks.json)

从仓库根目录运行：

```bash
cd wt_relative_classification
python -m pip install -r requirements.txt
python classify.py > run.log 2>&1
python report.py
```

运行 report.py 会从已保存结果重新生成HTML报告和核验记录；无需再次训练。若完整重跑，先运行 classify.py。


在此目录运行 `python classify.py`，完成后运行 `python report.py`。
依赖：numpy pandas scipy scikit-learn joblib threadpoolctl；来源版本见 versions.json。

Δ编辑效率（百分点）=突变体平均编辑率−该位点WT平均编辑率。
类别1：Δ≥动态分界；类别0：Δ<动态分界。阈值扫描为−100至100 pp的5 pp步长加全部相邻观测值中点，等价标签去重，每类至少5个构建。
模型：Total CP / Interface CP / S4-S12-S24 / Repository特征逻辑回归，ContactSeek Top3 / CCR浅层随机森林。
inputs保留已提取CP、实验记录和结构特征；sources保留本轮调用的ContactSeek函数与许可证；来源记录见source_provenance.json。源结构压缩包未包含；可从已提取输入复现本轮全部分类。
audit_reference保留上一轮折外分数，用于逐分组核验平移不变性。

编辑分界与模型经全体数据探索筛选，最终成绩未消除该层选择偏差。模型分数阈值和CCR在训练折内确定。组合突变的位置分组并非严格的所有组成位置隔离。
主要文件：results/wt_relative_summary.csv，results/wt_relative_per_construct_predictions.csv，results/validation_checks.json，WT_relative_classification_report.html。
