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
        cellInd = []; % cell number as index
        cellDepth = [];
        frameRate = [];
        imagingDepth = [];
        pixResolution = [];
        rollingWindowForBaseF = 100; % in s
        baseFprctile = 5;
        active = []; % records if the cell is active (1) or not (0). 
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
            
            cellNums{1} = []; cellDepth{1} = [];
            cellNums{2} = []; cellDepth{2} = [];

            fnlist = dir(['F_', mouseName, '_', sessionName, '_plane*_proc_final.mat']);
            for i = 1 : length(fnlist)
                load(fnlist(i).name);
                obj.mimg{i} = dat.ops.mimg1;
                obj.cellmap{i} = zeros(size(dat.ops.mimg1),'single');
                tempCellmap = zeros([length(dat.ops.yrange), length(dat.ops.xrange)],'single');
                
                numCells = 0;
                for j = 1 : length(dat.stat)
                    fprintf('processing cells %d/%d of plane #%d\n', j, length(dat.stat),i)
                    if dat.stat(j).iscell
                        numCells = numCells + 1;
                        tempCellmap(dat.stat(j).ipix) = numCells + i*1000;
                        obj.cellInd(end+1) = numCells + i*1000;
                        obj.cellDepth(end+1) = round(dat.stat(j).depth);                        
                        obj.dF{end+1} = dat.dF(j,:);
                        if i <= dat.ops.num_plane
                            cellNums{1} = [cellNums{1}, numCells + i*1000];
                            cellDepth{1} = [cellDepth{1}, dat.stat(j).depth];
                        else
                            cellNums{2} = [cellNums{2}, numCells + i*1000];
                            cellDepth{2} = [cellDepth{2}, dat.stat(j).depth];
                        end                        
                    end
                end
                obj.cellmap{i}(dat.ops.yrange,dat.ops.xrange) = tempCellmap;
            end
            
            load(fnlist(1).name)            
            obj.frameRate = dat.ops.imageRate / dat.ops.num_plane;
            numTrials = length(dat.ops.trials);
            obj.trials = cell(numTrials,1);
            for i = 1 : length(obj.trials)
                obj.trials{i}.trialNum = dat.ops.trials(i).trialnum;
                obj.trials{i}.frames = dat.ops.trials(i).frames(1):dat.ops.trials(i).frames(2);                
                if ~isempty(find(ismember(dat.ops.frame_to_use{4},obj.trials{i}.frames),1,'first'))
                    obj.trials{i}.planes = 1:4;
                    basePlane = 4;
                    obj.trials{i}.cellNums = cellNums{1};
                    obj.trials{i}.cellDepth = cellDepth{1};
                else                    
                    obj.trials{i}.planes = 5:8;
                    basePlane = 8;
                    obj.trials{i}.cellNums = cellNums{2};
                    obj.trials{i}.cellDepth = cellDepth{2};
                end
                obj.trials{i}.inds = find(ismember(dat.ops.frame_to_use{basePlane},obj.trials{i}.frames));
                obj.trials{i}.frameNum = length(obj.trials{i}.inds);                
                frames = dat.ops.frame_to_use{basePlane}(obj.trials{i}.inds) - obj.trials{i}.frames(1);                
                for j = 1 : 4
                    obj.trials{i}.time{5-j} = (frames + (j-1)) / dat.ops.imageRate;
                end
                obj.trials{i}.dF = zeros(length(obj.trials{i}.cellNums),obj.trials{i}.frameNum);
               for j = 1 : length(obj.trials{i}.cellNums) 
                    obj.trials{i}.dF(j,:) = obj.dF{obj.cellInd == obj.trials{i}.cellNums(j)}(obj.trials{i}.inds);
                end
            end
            
            obj.imagingDepth = dat.depth.imagingDepth;
            obj.pixResolution = dat.depth.xyUmperpix;
            
        end
    end

    methods % Dependent property methods; cannot have attributes.

    end
end     