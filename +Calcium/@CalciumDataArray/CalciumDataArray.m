classdef CalciumDataArray < handle
    % calcium data array from two-photon microscopy
    % sbx file from Neurolabware scanbox, processed by Suite2P (m
    % Pachitariu). 
    % Currently optimized for 2pad experiments, 2 blocks, 4 planes in each
    % block. 
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
        F = {}; % cell calcium trace
        cellInd = []; % cell number as index
        frameRate = [];
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
            
            cellNums{1} = [];
            cellNums{2} = [];

            fnlist = dir(['F_', mouseName, '_', sessionName, '_plane*_proc.mat']);            
            for i = 1 : length(fnlist)
                load(fnlist(i).name);
                obj.mimg{i} = dat.ops.mimg(dat.ops.yrange, dat.ops.xrange);
                obj.cellmap{i} = zeros([length(dat.ops.yrange), length(dat.ops.xrange)],'single');
                
                numCells = 0;
                for j = 1 : length(dat.stat)
                    if dat.stat(j).iscell
                        numCells = numCells + 1;
                        obj.cellmap{i}(dat.stat(j).ipix) = numCells;
                        obj.cellInd(end+1) = numCells + i*1000;
                        obj.F{end+1} = dat.Fcell{1}(j,:) - dat.FcellNeu{1}(j,:) * dat.stat(j).neuropilCoefficient;
                        if i <= 4
                            cellNums{1} = [cellNums{1}, numCells + i*1000];
                        else
                            cellNums{2} = [cellNums{2}, numCells + i*1000];
                        end
                    end
                end
            end
            
            load(fnlist(1).name)            
            obj.frameRate = dat.ops.imageRate / 4;
            numTrials = length(dat.ops.trials);
            obj.trials = cell(numTrials,1);
            for i = 1 : length(obj.trials)
                obj.trials{i}.trialNum = dat.ops.trials(i).trialnum;
                obj.trials{i}.frames = dat.ops.trials(i).frames(1):dat.ops.trials(i).frames(2);                
                if ~isempty(find(ismember(dat.ops.frame_to_use{4},obj.trials{i}.frames),1,'first'))
                    obj.trials{i}.planes = 1:4;
                    basePlane = 4;
                    obj.trials{i}.cellNums = cellNums{1};                    
                else                    
                    obj.trials{i}.planes = 5:8;
                    basePlane = 8;
                    obj.trials{i}.cellNums = cellNums{2};
                end
                obj.trials{i}.inds = find(ismember(dat.ops.frame_to_use{basePlane},obj.trials{i}.frames));
                obj.trials{i}.frameNum = length(obj.trials{i}.inds);                
                frames = dat.ops.frame_to_use{basePlane}(obj.trials{i}.inds) - obj.trials{i}.frames(1);                
                for j = 1 : 4
                    obj.trials{i}.time{5-j} = (frames + (j-1)) / dat.ops.imageRate;
                end
                obj.trials{i}.F = zeros(length(obj.trials{i}.cellNums),obj.trials{i}.frameNum);
                for j = 1 : length(obj.trials{i}.cellNums)
                    obj.trials{i}.F(j,:) = obj.F{obj.cellInd == obj.trials{i}.cellNums(j)}(obj.trials{i}.inds);
                end
            end            
        end
    end

    methods % Dependent property methods; cannot have attributes.

    end
end     