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
    end
    
    properties (Dependent = true)
    end
    
    methods (Access = public)
%         function obj = Uber_2padArray(bArray,w3Array,caArray)
        function obj = Uber_2padArray(bArray,wlArray,w3Array,caArray)
            obj.mouseName = bArray.mouseName;
            obj.sessionName = bArray.sessionName;
            obj.cellNums = caArray.cellInd;
            obj.cellDepths = caArray.cellDepth;
            
            btrialNums = cellfun(@(x) x.trialNum, bArray.trials);            
            wltrialNums = cellfun(@(x) x.trialNum, wlArray.trials);
            w3trialNums = cellfun(@(x) x.trialNum, w3Array.trials);
            catrialNums = cellfun(@(x) x.trialNum, caArray.trials);
            obj.trialNums = intersect(btrialNums, intersect(wltrialNums, intersect(w3trialNums,catrialNums)));
            
            obj.trials = cell(length(obj.trialNums),1);
            for i = 1 : length(obj.trials)
                disp(['Building trial ', num2str(i), '/', num2str(length(obj.trials))]) 
                b = bArray.trials{cellfun(@(x) x.trialNum == obj.trialNums(i), bArray.trials)};
                wl = wlArray.trials{cellfun(@(x) x.trialNum == obj.trialNums(i), wlArray.trials)};
                w3 = w3Array.trials{cellfun(@(x) x.trialNum == obj.trialNums(i), w3Array.trials)};
                ca = caArray.trials{cellfun(@(x) x.trialNum == obj.trialNums(i), caArray.trials)};
%                 obj.trials{i} = Uber.Uber_2pad(b,w3,ca,obj.trialNums(i));
                obj.trials{i} = Uber.Uber_2pad(b,wl,w3,ca,obj.trialNums(i));
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