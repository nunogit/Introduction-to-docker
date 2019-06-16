def znorm(x):
    """
    Standardize variable to mean = 0, std = 1
    """
    import numpy as np
    
    x_np = np.asarray(x)
    return (x_np - x_np.mean()) / x_np.std()

def load_data(params):
    """
   
    Inputs: 
    params: name of the parameters used in model. e.g. 1568
    
    Outputs:
    data_og: dataframe with original Diogenes variables 
    data_sim: dataframe with interpolated curves at 1min resolution
    var: list with variable names
    abs_target: absolute residuals
    signed_target: residuals with original sign
    """
    import numpy as np
    import pandas as pd
    from support_functions import znorm
    
    data_og = pd.read_csv('diogenes80_noNaN.csv',index_col='Unnamed: 0')
    residuals = pd.read_csv('residuals_EDES_'+str(params)+'.csv', header=None) # 1257, 568, 1568
    data_sim = pd.read_csv('metric_data.csv', index_col='Unnamed: 0')

    residuals.index = data_og.index

    var = ['X3.OH.butyrate',
           'oxo.isovaleric.acid', 'proline', 'tyrosine', 'glucose',
           'triglycerides', 'leucine', 'creatine', 'valine', 'creatinine',
           'isoleucine', 'acetate', 'hydroxyisobutyrate', 'lactate', 'alanine',
           'glycine', 'n.acetyl', 'acetoacetate', 'glucose_CID1',
           'isoleucine_CID1', 'leucine_CID1', 'oxo.isovaleric.acid_CID1',
           'hydroxyisobutyrate_CID1', 'triglycerides_CID1', 'valine_CID1',
           'lactate_CID1', 'creatinine_CID1', 'tyrosine_CID1', 'proline_CID1',
           'alanine_CID1', 'acetate_CID1', 'n.acetyl_CID1', 'acetoacetate_CID1',
           'glycine_CID1', 'creatine_CID1', 'X3.OH.butyrate_CID1', 'ffa1',
           'dmatsu1', 'dfp1', 'crp1', 'trig1', 'fibri1', 'dldl1',
           'adipo1', 'chol1', 'fvii1', 'hdl1','aci10384', 'aci10381',
           'aci10353', 'aci10360', 'sc0250', 'gender',
           'BMI', 'WHR', 'dhomaB', 'dhomaIR', 'dhomaIS']

    var1 = ['ffa1','dmatsu1', 'dfp1', 'crp1', 'trig1', 'fibri1', 'dldl1',
           'adipo1', 'chol1', 'hdl1', 'aci10384', 'aci10381',
           'aci10353', 'sc0250','BMI', 'WHR', 'dhomaIS', 'gender']
    
    target = znorm(abs(residuals.iloc[:,0:5].sum(axis=1))) # residuals, the target variable -- skewed
    logtarget = znorm(np.log(abs(residuals.iloc[:,0:5].sum(axis=1)))) # log transformed residuals
    abs_target = abs(residuals).iloc[:,0:5].sum(axis=1)
    signed_target = residuals.iloc[:,0:5].sum(axis=1)
    
    res_0 = residuals.iloc[:,0]
    res_30 = residuals.iloc[:,1]
    res_60 = residuals.iloc[:,2]
    res_90 = residuals.iloc[:,3]
    res_120 = residuals.iloc[:,4]
    
    return data_og, data_sim, var, abs_target, signed_target, res_0, res_30, res_60, res_90, res_120

def encode_residuals(ogtt, res_0, res_30, res_60, res_90, res_120):
    
    import numpy as np
    
    pattern = np.zeros((410,5))
    pattern[:,0] = ogtt.iloc[:,1].values
    pattern[:,1] = (res_0-res_30/30)
    pattern[:,2] = (res_30-res_60/30)
    pattern[:,3] = (res_60-res_90/30)
    pattern[:,4] = (res_90-res_120/30)
    
    return pattern

def train_predictor(X, y, var, model):
    
    import numpy as np
    import pandas as pd    
    from sklearn.model_selection import cross_validate, cross_val_predict, cross_val_score, LeaveOneGroupOut, LeavePGroupsOut, KFold, LeaveOneOut
    from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
    from sklearn.linear_model import ElasticNet, Ridge
    from sklearn import metrics
    from scipy.stats import pearsonr
    
    ytests = []
    ypreds = []
    feats = []
    train_score = []
    cv = LeaveOneOut()
    for train_idx, test_idx in cv.split(X):
        X_train, X_test = X.values[train_idx], X.values[test_idx]
        y_train, y_test = y.values[train_idx], y.values[test_idx]
        
        if model == 'Elasticnet':
            M = ElasticNet(alpha=0.2, copy_X=True, fit_intercept=True, 
                           l1_ratio=0.5, max_iter=10000, normalize=False, positive=False, 
                           precompute=False, random_state=0, selection='cyclic', tol=0.0001, warm_start=False)
        elif model == 'GBR': # subsample=0.8, learning_rate=0.003, min_samples_leaf=8,
            M = GradientBoostingRegressor(n_estimators=2000, learning_rate=0.002, min_samples_leaf=4, random_state=1, loss='ls')
        elif model == 'Randomforest':
            M = RandomForestRegressor(bootstrap=True, criterion='mse', max_depth=2, max_features='auto', max_leaf_nodes=None, min_impurity_decrease=0.0, 
                                      min_impurity_split=None, min_samples_leaf=4, min_samples_split=2, min_weight_fraction_leaf=0.0, 
                                      n_estimators=1000, n_jobs=None, oob_score=False, random_state=0, verbose=0, warm_start=False)
            
        est = M 
        
        est.fit(X_train, y_train)
            
            
        y_pred = est.predict(X_test)
        
        if model == 'Elasticnet':
            feats.append(est.coef_)
            mean_imp = np.asarray(feats).mean(axis=0)
            std_imp = np.asarray(feats).std(axis=0)
            CV_features = pd.DataFrame(mean_imp,columns=['avg_importance'], index=var)
            CV_features['std_importance'] = std_imp
        elif model == 'GBR':
            feats.append(est.feature_importances_)
            mean_imp = np.asarray(feats).mean(axis=0)
            std_imp = np.asarray(feats).std(axis=0)
            CV_features = pd.DataFrame(mean_imp,columns=['avg_importance'], index=var)
            CV_features['std_importance'] = std_imp
        elif model == 'Randomforest':
            feats.append(est.feature_importances_)
            mean_imp = np.asarray(feats).mean(axis=0)
            std_imp = np.asarray(feats).std(axis=0)
            CV_features = pd.DataFrame(mean_imp,columns=['avg_importance'], index=var)
            CV_features['std_importance'] = std_imp

        ytests += list(y_test)
        ypreds += list(y_pred)

        train_score.append(est.score(X_train, y_train))

    train_R2 = np.asarray(train_score).mean(axis=0)
    train_R2_std = np.asarray(train_score).std(axis=0)

    test_R2 = metrics.r2_score(ytests, ypreds)
    test_MSE= metrics.mean_squared_error(ytests, ypreds)


    

        
    return train_R2, train_R2_std, test_R2, test_MSE, CV_features, ypreds, ytests