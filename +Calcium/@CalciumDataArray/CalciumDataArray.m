classdef CalciumDataArray < handle
    % calcium data array from two-photon microscopy
    % sbx file from Neurolabware scanbox, processed by Suite2P (m
    % Pachitariu). 
    % Currently optimized for 2pad experiments, 2 blocks, 4 planes in each
    % block. 
    
    % assume dF/F is calculated already and saved in F_*_proc_final.mat
    % 'dat'
    
    % Inputs: 
%         suite2p files (Fmouse_session_plane#_proc.mat)
%         mat file from scanbox
%         .trials file made during suite2p (custom code)

    properties
        mouseName = '';
        sessionName = '';
        mimg = {}; 
        cellmap = {}; % cell roi of each plane. Each ROI has it's corresponding cell #, except for the first digit reserved for plane #.
        trials = {};
        dF = {}; % dF/F_0 with 5th percentile
        spk = {}; % inferred spikes from MLspike
        cellNums = []; % cell number as index        
        cellDepths = [];
        isC2 = []; % 0 or 1
        celly = []; % mid point in y axis. In whole FOV index
        cellx = []; % mid point in x axis. In whole FOV index
        c2ypoints = [];
        c2xpoints = [];
        frameRate = [];
        imagingDepth = [];
        pixResolution = [];
        rollingWindowForBaseF = 100; % in s
        baseFprctile = 5;
        noise = [];
        npcoeff = [];
        fovsize = [];
        fovyrange = {};
        fovxrange = {};
        fovdepth = {};
    end

    properties (Dependent = true)

    end

    methods (Access = public)
        function obj = CalciumDataArray(mouseName, sessionName)
            % Assume this function is called at the relevant directory
            if isnumeric(mouseName)
                mouseName = sprintf('%03d',mouseName);
            end
            if isnumeric(sessionName)
                sessionName = sprintf('%03d',sessionName);
            end
            obj.mouseName = mouseName;
            obj.sessionName = sessionName;
            
            cellNums{1} = []; cellDepths{1} = [];
            cellNums{2} = []; cellDepths{2} = [];

%             fnlist = dir(['F_', mouseName, '_', sessionName, '_plane*_proc_final_spikes_noiseCorrected.mat']);
            fnlist = dir(['F_', mouseName, '_', sessionName, '_plane*_proc_final.mat']);
            for i = 1 : length(fnlist)
                load(fnlist(i).name);
                %
                %
                %
                %

                if strcmp(mouseName, '027') && (strcmp(sessionName, '009') || strcmp(sessionName, '010'))
                    pi = i + 4; % only for JK027 S09 and S10
                else
                    pi = i; % all otherwise (the ones with both layers imaged)
                end
                %
                %
                %
                %
                %
                %
                
                obj.mimg{pi} = dat.ops.mimg1;
                obj.cellmap{pi} = zeros(size(dat.ops.mimg1),'single');
                tempCellmap = zeros([length(dat.ops.yrange), length(dat.ops.xrange)],'single');
                
                inds = find([dat.stat.iscell]);
                for j = 1 : length(inds)
                    fprintf('processing cells %d/%d of plane #%d\n', j, length(inds),pi)
                    tempCellmap(dat.stat(inds(j)).ipix) = j + pi*1000;
                    obj.cellNums(end+1) = j + pi*1000;
                    if isfield(dat.stat,'depth')
                        obj.cellDepths(end+1) = round(dat.stat(inds(j)).depth);
                    end
                    obj.isC2(end+1) = dat.isC2(j);
                    obj.celly(end+1) = dat.ops.yrange(round(dat.stat(inds(j)).med(1)));
                    obj.cellx(end+1) = dat.ops.xrange(round(dat.stat(inds(j)).med(2)));
                    obj.dF{end+1} = dat.dF(j,:);
                    obj.spk{end+1} = spk.n(j,:);
                    obj.noise(end+1) = dat.noise(j);
                    obj.npcoeff(end+1) = dat.npcoeffs(j);
                    if pi <= dat.ops.num_plane
                        cellNums{1} = [cellNums{1}, j + pi*1000];
                        if isfield(dat.stat,'depth')
                            cellDepths{1} = [cellDepths{1}, dat.stat(inds(j)).depth];
                        end
                    else
                        cellNums{2} = [cellNums{2}, j + pi*1000];
                        if isfield(dat.stat,'depth')
                            cellDepths{2} = [cellDepths{2}, dat.stat(inds(j)).depth];
                        end
                    end
                end
                obj.cellmap{pi}(dat.ops.yrange,dat.ops.xrange) = tempCellmap;
                obj.fovyrange{pi} = dat.ops.useY(dat.ops.yrange);
                obj.fovxrange{pi} = dat.ops.useX(dat.ops.xrange);
                if isfield(dat.stat,'depth')
                    obj.fovdepth{pi} = dat.depth.FOV;
                end
            end
            
            load(fnlist(1).name)
            obj.frameRate = dat.ops.imageRate / dat.ops.num_plane;
            numTrials = length(dat.ops.trials);
            obj.trials = cell(numTrials,1);
            for i = 1 : length(obj.trials)
                obj.trials{i}.trialNum = dat.ops.trials(i).trialnum;
                obj.trials{i}.frames = dat.ops.trials(i).frames(1):dat.ops.trials(i).frames(2);
                obj.trials{i}.startLine = dat.ops.trials(i).lines(1);
                if ~isempty(find(ismember(dat.ops.frame_to_use{4},obj.trials{i}.frames),1,'first'))
                    obj.trials{i}.planes = 1:4;
                    basePlane = 4;
                    obj.trials{i}.cellNums = cellNums{1};
                    obj.trials{i}.cellDepths = cellDepths{1};
                else
                    obj.trials{i}.planes = 5:8;
                    basePlane = 8;
                    obj.trials{i}.cellNums = cellNums{2};
                    obj.trials{i}.cellDepths = cellDepths{2};
                end
                obj.trials{i}.inds = find(ismember(dat.ops.frame_to_use{basePlane},obj.trials{i}.frames)); % conservative treatment for frame_to_use starts before TTL signal, just in case.
                obj.trials{i}.frameNum = length(obj.trials{i}.inds);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Very important. Calculating time points %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                frames = dat.ops.frame_to_use{basePlane}(obj.trials{i}.inds) - obj.trials{i}.frames(1); % frames are # of starting from the trial onset. Matching with TTL1 signal, 
                %   which is from BControl. It's because frame_to_use is calculated per volume, to match # of frames in each planes from the same layer.
                % This way of calculating "frames" is specific (and correct) for calculation next line. (for the part where (j-1) instead of (j) is used.)
                % It's bacause "frames" CANNOT be 0, because frame_to_use is calculated to have at least 1 frame after the TTL1 signal has arrived.
                for j = 1 : 4
                    obj.trials{i}.time{5-j} = (frames + (j-1)) / dat.ops.imageRate - (obj.trials{i}.startLine - 1) / dat.ops.imageRate / dat.ops.info.sz(1); % adjusting for start line. 
                    % offset can be as large as 30 ms. 
                    % This is the BEGINNING time point of each frame. For any event to be assigned to the corresponding frame,
                    %   that frame's time point should be should be LARGER THAN OR EQUAL TO (>=) the event timing.
                end
                obj.trials{i}.dF = zeros(length(obj.trials{i}.cellNums),obj.trials{i}.frameNum);
                obj.trials{i}.spk = zeros(length(obj.trials{i}.cellNums),obj.trials{i}.frameNum);
                for j = 1 : length(obj.trials{i}.cellNums) 
                    obj.trials{i}.dF(j,:) = obj.dF{obj.cellNums == obj.trials{i}.cellNums(j)}(obj.trials{i}.inds);
                    obj.trials{i}.spk(j,:) = obj.spk{obj.cellNums == obj.trials{i}.cellNums(j)}(obj.trials{i}.inds);
                end
            end
            if isfield(dat.stat,'depth')
                obj.imagingDepth = dat.depth.imagingDepth;
                obj.pixResolution = dat.depth.xyUmperpix;
            end
            obj.c2ypoints = dat.c2ypoints - dat.ops.useY(1)+1;
            obj.c2xpoints = dat.c2xpoints - dat.ops.useX(1)+1;
            obj.fovsize = dat.ops.info.sz;
            
%             obj.cellNums = [cellNums{1}, cellNums{2}];
%             obj.cellDepths = [cellDepths{1}, cellDepths{2}];
        end
    end

    methods % Dependent property methods; cannot have attributes.

    end
end