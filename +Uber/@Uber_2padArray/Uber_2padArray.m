classdef Uber_2padArray < handle
    properties
        mouseName = '';
        sessionName = '';
        trialNums = [];
        trials = {};
        cellNums = [];
        cellDepths = [];
        frameRate = [];
        mimg = {};
        cellmap = {};
        noise = [];
        npcoeff = [];
        isC2 = [];
        celly = [];
        cellx = [];
        c2ypoints = [];
        c2xpoints = [];
        fovsize = [];
        fovyrange = {};
        fovxrange = {};
        fovdepth = {}; 
        pixResolution = [];
        planeTrialInds = {}; % 1: indices of trials (of u.trials) for plane #1 (1-4)
                             % 2: indices of trials (of u.trials) for plane #2 (5-8)
    end
    
    properties (Dependent = true)
    end
    
    methods (Access = public)
%         function obj = Uber_2padArray(bArray,w3Array,caArray)
        function obj = Uber_2padArray(bArray,wfArray,caArray)
            obj.mouseName = bArray.mouseName;
            obj.sessionName = bArray.sessionName;
            obj.cellNums = caArray.cellNums;
            obj.cellDepths = caArray.cellDepths;
            
            btrialNums = cellfun(@(x) x.trialNum, bArray.trials);            
            wftrialNums = cellfun(@(x) x.trialNum, wfArray.trials);
            catrialNums = cellfun(@(x) x.trialNum, caArray.trials);
            obj.trialNums = intersect(btrialNums, intersect(wftrialNums,catrialNums));
            
            obj.trials = cell(length(obj.trialNums),1);
            for i = 1 : length(obj.trials)
                disp(['Building trial ', num2str(i), '/', num2str(length(obj.trials))]) 
                b = bArray.trials{cellfun(@(x) x.trialNum == obj.trialNums(i), bArray.trials)};
                wf = wfArray.trials{cellfun(@(x) x.trialNum == obj.trialNums(i), wfArray.trials)};
                ca = caArray.trials{cellfun(@(x) x.trialNum == obj.trialNums(i), caArray.trials)};
                obj.trials{i} = Uber.Uber_2pad(b,wf,ca,obj.trialNums(i));
            end
            
            obj.frameRate = caArray.frameRate;
            obj.mimg = caArray.mimg;
            obj.cellmap = caArray.cellmap;
            obj.noise = caArray.noise;
            obj.npcoeff = caArray.npcoeff;
            obj.isC2 = caArray.isC2;
            obj.celly = caArray.celly;
            obj.cellx = caArray.cellx;
            obj.c2ypoints = caArray.c2ypoints;
            obj.c2xpoints = caArray.c2xpoints;
            obj.fovsize = caArray.fovsize;
            obj.fovyrange = caArray.fovyrange;
            obj.fovxrange = caArray.fovxrange;
            obj.fovdepth = caArray.fovdepth;
            obj.pixResolution = caArray.pixResolution;
        end
    end
    
    methods % Dependent property methods; cannot have attributes.
        
    end
    
end