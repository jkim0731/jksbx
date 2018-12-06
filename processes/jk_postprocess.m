%% First, make db from jk_make_db, then run jk_master_file (suite2p)

%% Then, curate ROIs

%% then, process z-stack and C2
% process_zstack
% draw_C2(fn)


% assign_depth % this gives you _proc_final file
% append_isC2(fn)
fn = '070_702_000';



%% then, sort active cells and make dF/F0

% sorting, it's included in datadFcorrection

% datdFcorrection

%% then, spike deconvolution
% mlspike_on_suite2p