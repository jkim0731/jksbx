%% SET ALL DEFAULT OPTIONS HERE

% UPDATE fall 2017: non-rigid and rigid registration scripts merged; red
% channel mean image can be computed while registering green channel; red
% channel binary can be computed during green channel registration
% (ly x lx x time like green channel)

% UPDATE end-of-summer 2017: default neuropil extraction is now "surround"
% and it's very fast. Cell extraction is on the raw data (no pixel-scaling or smoothing). 

% UPDATE summer 2017: default spike deconvolution changed to a customized version of
% OASIS (due to our results in this paper http://www.biorxiv.org/content/early/2017/06/27/156786). Please
% Please download the OASIS code from https://github.com/zhoupc/OASIS_matlab, and
% add the folder, with its subfolders, to your Matlab path. 

% UPDATE Christmas 2016: number of clusters determined automatically, but
% do specify the "diameter" of an average cell for best results. You can do this with either
% db(iexp).diameter, or ops0.diameter

% check out the README file for detailed instructions
% **** and for more options available ****

addpath('C:\Users\shires\Documents\GitHub\jksbx\Suite2P-master') % add the path to your make_db file

% overwrite any of these default options in your make_db file for individual experiments
jk_make_db; % RUN YOUR OWN MAKE_DB SCRIPT TO RUN HERE


ops0.optotune_ringing_time = 8; % in ms. To crop top portion of each frame.


ops0.toolbox_path = 'C:\Users\shires\Documents\GitHub\jksbx\Suite2P-master';
if exist(ops0.toolbox_path, 'dir')
	addpath(genpath(ops0.toolbox_path)) % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end

% mex -largeArrayDims SpikeDetection/deconvL0.c (or .cpp) % MAKE SURE YOU COMPILE THIS FIRST FOR DECONVOLUTION

ops0.useGPU                 = 1; % if you can use an Nvidia GPU in matlab this accelerates registration approx 3 times. You only need the Nvidia drivers installed (not CUDA).
ops0.fig                    = 1; % turn off figure generation with 0
% ops0.diameter               = 12; % most important parameter. Set here, or individually per experiment in make_db file
%                                   % being calculated in jk_run_pipeline according to the imaging magnification

% ---- root paths for files and temporary storage (ideally an SSD drive. my SSD is C:/)
ops0.RootStorage            = 'D:\TPM\JK\'; % Suite2P assumes a folder structure, check out README file
% ops0.temp_tiff              = 'D:\TPM\JK\temp.tif'; % copies each remote tiff locally first, into this file
ops0.RegFileRoot            = 'D:\TPM\JK\';  % location for binary file
ops0.DeleteBin              = 1; % set to 1 for batch processing on a limited hard drive
ops0.ResultsSavePath        = 'D:\TPM\JK\'; % a folder structure is created inside
ops0.RegFileTiffLocation    = ''; % leave empty to NOT save registered tiffs (slow)
% if you want to save red channel tiffs, also set ops0.REDbinary = 1

% ---- registration options ------------------------------------- %
ops0.doRegistration         = 1; % skip (0) if data is already registered
ops0.showTargetRegistration = 1; % shows the image targets for all planes to be registered
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
ops0.NimgFirstRegistration  = 500; % number of images to include in the first registration pass 
ops0.nimgbegend             = 0; % frames to average at beginning and end of blocks
ops0.dobidi                 = 1; % infer and apply bidirectional phase offset
ops0.nonrigid               = 0;

% ---- cell detection options ------------------------------------------%
ops0.ShowCellMap            = 0; % during optimization, show a figure of the clusters
ops0.sig                    = 1;  % spatial smoothing length in pixels; encourages localized clusters
ops0.nSVDforROI             = 1000; % how many SVD components for cell clustering
ops0.NavgFramesSVD          = 1000; % how many (binned) timepoints to do the SVD based on
ops0.signalExtraction       = 'surround'; % how to extract ROI and neuropil signals: 
% ops0.writeSVDroi            = 1;
%  'raw' (no cell overlaps), 'regression' (allows cell overlaps), 
%  'surround' (no cell overlaps, surround neuropil model)

% ----- neuropil options (if 'surround' option) ------------------- %
% all are in measurements of pixels
ops0.innerNeuropil  = 3; % padding around cell to exclude from neuropil
ops0.outerNeuropil  = Inf; % radius of neuropil surround
% if infinity, then neuropil surround radius is a function of cell size
if isinf(ops0.outerNeuropil)
    ops0.minNeuropilPixels = 200; % minimum number of pixels in neuropil surround
    ops0.ratioNeuropil     = 2; % ratio btw neuropil radius and cell radius
    % radius of surround neuropil = ops0.ratioNeuropil * (radius of cell)
end
ops0.saveNeuropil          = 1;

% ----- spike deconvolution and neuropil subtraction options ----- %
% ops0.imageRate              = 30; % imaging rate (cumulative over
% planes!). Approximate, for initialization of deconvolution kernel.
% Calculated later by info from sbx.
ops0.doDeconvolution        = 0;
ops0.sensorTau              = 1; % decay half-life (or timescale). Approximate, for initialization of deconvolution kernel.
ops0.maxNeurop              = 0.8; % for the neuropil contamination to be less than this (sometimes good, i.e. for interneurons)

% ----- if you have a RED channel ---------------------- ------------%
ops0.AlignToRedChannel      = 0; % compute registration offsets using red channel
ops0.REDbinary              = 0; % make a binary file of registered red frames
% if db.expred, then compute mean red image for green experiments with red
% channel available while doing registration
ops0.redMeanImg             = 1; 
% for red cell detection (identify_redcells_sourcery.m)
% redratio = red pixels inside / red pixels outside
% redcell = redratio > mean(redratio) + redthres*std(redratio)
% notred = redratio < mean(redratio) + redmax*std(redratio)
ops0.redthres               = 1.5; % the higher the thres the less red cells
ops0.redmax                 = 1; % the higher the max the more NON-red cells

db0 = db;
%% RUN THE PIPELINE HERE

running_times = zeros(length(db0),1);
running_times_deconv = zeros(length(db0),1);
running_times_red = zeros(length(db0),1);

for iexp = 1:length(db0)
    %%
    tic_start = tic;
    db = db0(iexp);
    
    if ~isfield(db, 'RootStorage') || isempty(db.RootStorage)
        db.RootStorage = ops0.RootStorage;
    end
    cd([db.RootStorage, filesep, db.mouse_name]) 
    % add info from scanbox .mat data and calculate ops0.imageRate
    matfnlist = dir([db.RootStorage, filesep, db.mouse_name, filesep, sprintf('%s_%03d_',db.mouse_name,db.session), '*.sbx']);
    for ifile = 1 : length(matfnlist)
        matfnlist(ifile).name = [matfnlist(ifile).name(1:end-4), '.mat'];
        if db.session < 6000 && db.session > 4000 % spontaneous imaging
            load(matfnlist(ifile).name)
            info.event_id = [];
            info.frame = [];
            info.line = [];
            save(matfnlist(ifile).name, 'info')
        end
    end
    try    
        load(fullfile(matfnlist(1).folder, matfnlist(1).name)) % loading 'info' variable from scanbox .mat file. using the first one of the session    
    catch
        try
            load([db.RootStorage, filesep, db.mouse_name, filesep, matfnlist(1).name])
        catch
            continue            
        end
    end
    
    db.info = info;
    ops0.diameter = getOr(db, 'diameter', 10*str2double(info.config.magnification_list(info.config.magnification,:)));
    ops0.imageRate = db.info.resfreq*(2-db.info.scanmode)/db.info.sz(1); % calculating imaging rate based on the scanning mode.     
        
    if isfield(db.info, 'deadband')
        ops0.useX = max((1-db.info.scanmode)*100, round(db.info.deadband(1)/2) + 10) + 1 : ...
            db.info.sz(2)-round(db.info.deadband(2)/2)-10; % crop first 100 columns in bidirectional scanning. and also considering deadband (10 as buffer for bidirectional alignment)
    else
        ops0.useX = (1-db.info.scanmode)*100: db.info.sz(2)-10;
    end
    if db.info.volscan % chop off first optotune_ringing_time ms of each plane because of optotune ringing
        ops0.useY = round(ops0.optotune_ringing_time/ (1000/db.info.resfreq) *(2-db.info.scanmode)) : db.info.sz(1);
    else
        db.useY = 1: db.info.sz(1);
    end
    
    curr_dir = pwd;
    db.sbxfnlist = cell(length(matfnlist),1);    
      
    for ifile = 1 : length(matfnlist)
        [~, db.sbxfnlist{ifile}, ~] = fileparts(matfnlist(ifile).name);
        if db.session > 8000
            onFrames = laser_on_frames(db.sbxfnlist{ifile});
            jksbxsplittrial_nobitcode(db.sbxfnlist{ifile}, onFrames); % temporary solution for not sending bitcode during piezo stimulation before 2018/02
        elseif db.session > 4000
            onFrames = laser_on_frames(db.sbxfnlist{ifile});
            jksbxsplittrial(db.sbxfnlist{ifile}, onFrames);
        end
        load([db.sbxfnlist{ifile},'.trials'],'-mat')        
        if ifile == 1
            db.frame_to_use = frame_to_use; % frame_to_use is in cell format. 
            db.trials = trials; % trials is a structure, containing trialnum, frames, and lines
            db.max_idx = sbx_maxidx(db.sbxfnlist{ifile}); % max_idx is a 1d array containing maximum index up to each file including previous index. Aim to help file reading with multiple files in one session
        else               
            for iplane = 1 : length(frame_to_use)                    
                db.frame_to_use{iplane} = [db.frame_to_use{iplane}, db.max_idx(end) + 1 + frame_to_use{iplane}];
            end
            db.trials = [db.trials, trials];
            db.max_idx = [db.max_idx, db.max_idx(end)+sbx_maxidx(db.sbxfnlist{ifile})+1];
        end
    end
    
    db.blockimaging = blockimaging; db.num_layer = num_layer;   db.num_plane = num_plane;
    db.nplanes = db.num_layer * db.num_plane; % nplanes meaning total # of planes that are imaged, while num_plane #of planes during a single imaging block (or within a layer)
    db.nsessions = length(db.session);
    
    % for red channel detection and treatment
    if info.channels == 2 % pmt0 only
        db.expred = [];
%         ops0.REDbinary = 0;
        ops0.redMeanImg = 0;
        ops0.AlignToRedChannel = 0; % just to make sure
    elseif info.channels == 1 % pmt0 & pmt1 
        db.expred = db.session;
%         ops0.REDbinary = 1;
        ops0.redMeanImg = 1;
    else
        error('only red channel')
    end
    
    cd(curr_dir)
    %%
    jk_run_pipeline(db, ops0);
    
    running_times(iexp) = toc(tic_start)/60; % in min
    
    % deconvolved data into st, and neuropil subtraction coef in stat
    if isfield(db,'doDeconvolution') && db.doDeconvolution
        jk_add_deconvolution(ops0, db);
    end
    
    running_times_deconv(iexp) = toc(tic_start)/60; % in min
    
    % add red channel information (if it exists)
    if isfield(db,'expred') && ~isempty(db.expred)
        % creates mean red channel image aligned to green channel
        % use this if you didn't get red channel during registration
        % OR you have a separate experiment with red and green just for this
        red_expts = ismember(db.session, getOr(db, 'expred', []));
        if ~ops0.redMeanImg || sum(red_expts)==0
            jkrun_REDaddon_sourcery(db, ops0);
        end
        
        % identify red cells in mean red channel image
        % fills dat.stat.redcell, dat.stat.notred, dat.stat.redprob
        jk_identify_redcells_sourcery(db, ops0);         
    end
    
    running_times_red(iexp) = toc(tic_start)/60; % in min
    
end
%% STRUCTURE OF RESULTS FILE
% cell traces are in dat.Fcell
% neuropil traces are in dat.FcellNeu
% manual, GUI overwritten "iscell" labels are in dat.cl.iscell
%  
% stat(icell) contains all other information:
% iscell: automated label, based on anatomy
% neuropilCoefficient: neuropil subtraction coefficient, based on maximizing the skewness of the corrected trace (ICA)
% st: are the deconvolved spike times (in frames)
% c:  are the deconvolved amplitudes
% kernel: is the estimated kernel
