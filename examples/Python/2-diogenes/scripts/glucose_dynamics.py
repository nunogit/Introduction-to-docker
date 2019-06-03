
# coding: utf-8

import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_validate, cross_val_predict, cross_val_score, LeaveOneGroupOut, LeavePGroupsOut, KFold, LeaveOneOut
from sklearn.linear_model import Ridge, RidgeCV, ElasticNet, ElasticNetCV
from support_functions import znorm, load_data, encode_residuals, train_predictor

data_og, data_sim, var, abs_target, signed_target, res_0, res_30, res_60, res_90, res_120 = load_data('568__')

pattern = encode_residuals(data_og, res_0, res_30, res_60, res_90, res_120)
pattern_ = pd.DataFrame(pattern)

RESULTS = np.zeros((5,4))
mRESIDUALS = np.zeros((410,5)) #410 5 for pattern, 410 4 for residuals
pRESIDUALS = np.zeros((410,5)) #410 5 for pattern, 410 4 for residuals
features = []

for idx, i in enumerate(pattern_):
    print(i)
    X, y = data_og[var], pattern_.iloc[:,i] # different targets possible
    train_R2, train_R2_std, test_R2, test_MSE, CV_features, ypreds, ytests = train_predictor(X, y, var, 'GBR')
    
    
    
    RESULTS[0,idx] = train_R2
    RESULTS[1,idx] = train_R2_std
    RESULTS[2,idx] = test_R2
    RESULTS[3,idx] = test_MSE
    
    mRESIDUALS[:,idx] = ytests
    pRESIDUALS[:,idx] = ypreds
    
    features.append(CV_features)
features = pd.concat(features,axis=1)

#np.savetxt('hybrid/ML_568_independent_RESULTS.csv', RESULTS, delimiter=",", header='train_R2, train_R2_std, test_R2, test_MSE')
#np.savetxt('hybrid/ML_568_independent_RESIDUALSmes.csv', mRESIDUALS, delimiter=",", header='30, 60, 90, 120')
#np.savetxt('hybrid/ML_568_independent_RESIDUALSpred.csv', pRESIDUALS, delimiter=",", header='30, 60, 90, 120')
#features.to_csv('hybrid/features_per_timepoint.csv')
print('Done!')

