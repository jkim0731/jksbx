classdef Uber_2padArray < handle
    properties
        mouseName = '';
        sessionName = '';
        trialNums = [];
        trials = {};
        cellNums = [];
        cellDepths = [];
        cellActive = [];
        frameRate = [];
        mimg = {};
        cellmap = {};
    end
    
    properties (Dependent = true)
    end
    
    methods (Access = public)
        function obj = Uber_2padArray(bArray,wlArray,caArray)            
            obj.mouseName = bArray.mouseName;
            obj.sessionName = bArray.sessionName;
            obj.cellNums = caArray.cellInd;            
            obj.cellDepths = caArray.cellDepth;
            obj.cellActive = caArray.active;
            
            btrialNums = cellfun(@(x) x.trialNum, bArray.trials);
            wltrialNums = cellfun(@(x) x.trialNum, wlArray.trials);
            catrialNums = cellfun(@(x) x.trialNum, caArray.trials);
            obj.trialNums = intersect(btrialNums, intersect(wltrialNums, catrialNums));
            
            obj.trials = cell(length(obj.trialNums),1);
            for i = 1 : length(obj.trials)
                disp(['Building trial ', num2str(i), '/', num2str(length(obj.trials))]) 
                b = bArray.trials{cellfun(@(x) x.trialNum == obj.trialNums(i), bArray.trials)};
                wl = wlArray.trials{cellfun(@(x) x.trialNum == obj.trialNums(i), wlArray.trials)};
                ca = caArray.trials{cellfun(@(x) x.trialNum == obj.trialNums(i), caArray.trials)};
                obj.trials{i} = Uber.Uber_2pad(b,wl,ca,obj.trialNums(i));
            end
            
            obj.frameRate = caArray.frameRate;
            obj.mimg = caArray.mimg;
            obj.cellmap = caArray.cellmap;            
        end                
    end
    
    methods % Dependent property methods; cannot have attributes.
        
    end
    
end