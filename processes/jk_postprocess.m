%% First, make db from jk_make_db, then run jk_master_file (suite2p)

%% Then, curate ROIs

%% then, process z-stack and C2
% process_zstack % This can be ignored (as in L4 mice)
% draw_C2(fn)

% draw_C2('027_9998_109')
% draw_C2('036_9998_101')
% draw_C2('041_9998_1014')

% draw_C2('052_6002_005')
draw_C2('052_6002_004')


% draw_C2('074_999_002')
% draw_C2('075_999_003')
% draw_C2('076_999_001')

%%

process_zstack 
% it includes assign_depth % this gives you _proc_final file
% 
% append_isC2(fn)

% results in _proc_final.mat files

%% then, sort active cells and make dF/F0

% sorting, it's included in datadFcorrection

datdFcorrection

% then, spike deconvolution
mlspike_on_suite2p

%% aligning FOVs




















%% touch model GLM

% run GLM (touch model)


% save the function result from touch model
glm_results_saving_touch
% which uses
% glm_results_cell_function_touch
% results in 
% 'JK0xxSyyglm_cell_function_lasso_NC.mat' files in each suite2p sub-directory
% and
% 'cellFunctionLasso_NC.mat' in suite2p directory

% save each DE just in case (and for comparison with whisker model later)
glm_results_dev_exp_saving_touch
% which uses
% glm_results_dev_exp_touch
% results in 
% 'glmResults_devExp_touch_NC.mat' in suite2p directory

%
% % ?? DE diff saving - obsolete?
% d190716_glm_DE_diff_saving
% % results in
% % 'glm_DE_cellFunction_NC_JK%30dS%02d.mat' files in each suit2p sub-directory

%% Angle tuning

% run angle tuning analysis
angle_tuning_predecision
% results in 
% 'JK041S03angle_tuning_lasso_predecision_NC.mat' in each suite2p sub-directory

% save the summary data
angle_tuning_summary_saving
% results in 
% 'angle_tuning_summary_predecision.mat' in suite2p directory

%% whisker model GLM

% run GLM (whisker model)
% it requires touch model results, since it's run only from touch
% responsive cells

% save the results
glm_results_dev_exp_saving_wkv_touchCell
% which uses
% glm_results_dev_exp_wkv_touchCell
% results in
% 'glmResults_devExp_WKV_touchCell_NC.mat' in suite2p directory

% save results of exclusion method
glm_results_WKV_exclusion_saving
% results in 
% 'glmResult_WKV_NC_JK%03dS%02d' in each suite2p sub-directory
% and in
% 'glmResults_WKV_touchCell_exclusion' in suite2p directory

%% angle tuning model from wkv GLM

wkv_angle_tuning_each
% requires
% 'cellFunctionLasso_NC.mat'
% results in
% 'angle_tuning_model_lasso_NC_JK%03dS%02d' in each suite2p sub-directory



