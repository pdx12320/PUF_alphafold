import os,json,collections
import joblib,numpy as np,pandas as pd
from scipy.special import expit
O=os.path.dirname(__file__);m=joblib.load(O+'/final_summary_XGB.joblib');cols=json.load(open(O+'/summary_model_manifest.json'))['columns'];rules={};constant=0.
for txt in m[-1].get_booster().get_dump(dump_format='json'):
 t=json.loads(txt)
 if 'leaf' in t:constant+=t['leaf'];continue
 j=int(t['split'][1:]);key=(j,t['split_condition']);leaf={v['nodeid']:v['leaf'] for v in t['children']}
 if key not in rules:rules[key]={'feature':cols[j],'feature_index':j,'z_threshold':t['split_condition'],'raw_threshold_approx':float(t['split_condition']*m[1].scale_[j]+m[1].mean_[j]),'leaf_if_less':0.,'leaf_if_greater_equal':0.,'tree_count':0}
 rule=rules[key];rule['leaf_if_less']+=leaf[t['yes']];rule['leaf_if_greater_equal']+=leaf[t['no']];rule['tree_count']+=1
out={'base_margin':constant,'rules':list(rules.values()),'input_means':m[1].mean_.tolist(),'input_scales':m[1].scale_.tolist(),'columns':cols,'prediction':'sigmoid(base_margin + sum(rule leaf)); standardized inputs cast to float32 before comparison','note':'Raw thresholds are for interpretation; use the saved pipeline at numerical boundaries.'};json.dump(out,open(O+'/XGBoost_scoring_formula.json','w'),indent=2);pd.DataFrame(out['rules']).to_csv(O+'/XGBoost_scoring_rules.csv',index=False)
a=pd.read_csv(O+'/all_construct_features.csv');xx=m[1].transform(a[cols].values).astype(np.float32);margin=np.full(len(a),constant)
for rule in out['rules']:margin+=np.where(xx[:,rule['feature_index']]<np.float32(rule['z_threshold']),rule['leaf_if_less'],rule['leaf_if_greater_equal'])
error=float(np.max(abs(expit(margin)-m.predict_proba(a[cols].values)[:,1])));assert error<1e-6
json.dump({'boost_formula_max_absolute_error':error,'raw_formula_scorer_parity_absolute_error':float(abs(pd.read_csv(O+'/scorer_raw_check.csv').exploratory_score[0]-pd.read_csv(O+'/final_formula_training_scores.csv').query('construct=="puf12_r123_r567_r567_r678"').exploratory_score.iloc[0]))},open(O+'/validation_checks.json','w'),indent=2)
print(pd.DataFrame(out['rules']).to_string(index=False))
