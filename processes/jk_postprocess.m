%% First, make db from jk_make_db, then run jk_master_file (suite2p)

%% Then, curate ROIs

%% then, process z-stack and C2
% process_zstack % This can be ignored (as in L4 mice)
% draw_C2(fn)

draw_C2('074_999_002')
draw_C2('075_999_003')
draw_C2('076_999_001')

% assign_depth % this gives you _proc_final file
% append_isC2(fn)

% results in _proc_final.mat files

%% then, sort active cells and make dF/F0

% sorting, it's included in datadFcorrection

% datdFcorrection

%% then, spike deconvolution
% mlspike_on_suite2p