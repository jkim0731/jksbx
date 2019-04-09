classdef Uber_2pad < handle
    properties
        trialNum = [];
        
        % from the Task (Solo)
        leftLickTime = []; % in sec
        rightLickTime = []; % in sec
        answerLickTime = []; % in sec
        drinkingOnsetTime = []; % in sec. Actual drinking time, i.e., intersect(b.drinkingTime, b.left(or right)LickTime)
        poleUpOnsetTime = []; % in sec
        poleDownOnsetTime = []; % in sec
        response = []; % correct = 1, wrong = 0, miss = -1
        angle = []; % in degrees
        distance = []; % lateral distance
        position = []; % anterior-posterior position
        taskTarget = ''; % 'Angle' or 'Distance'
        distractor = ''; % 'On' or 'Off'
        trialType = ''; % rn, ln, rc, rf, lc, lf, ra, la, oo
        
        % from WF
        poleMovingTime = []; % in sec
        poleUpTime = []; % in sec
        whiskerTime = []; % in sec, all frames
        theta = [];
        nof = 0;
        frameDuration = 1/311;
        
        phi = [];
        kappaH = [];
        kappaV = [];
        
        protractionTouchChunks = {};
        retractionTouchChunks = {};
        protractionTouchDuration = [];
        protractionTouchSlideDistance = [];
        
%         protractionTouchChunksByWhisking = {};
%         retractionTouchChunksByWhisking = {};
        
        % from two-photon data
        planes = []; % 1:4 or 5:8 (for now, 2018/04/01)
        neuindSession = []; % index (number) of neurons observed in the current session, 
% which means that there is a single value assigned to each neuron every session (not same across sessions, yet. Going to be dealt in _animal array)
% This is matched with F matrix.
% 4 digits, assuming # of neurons can be larger than 100 in one plane
% (maybe in the future, low mag imaging), first digit represents the plane.
       
%         neuindAnimal = []; % index (number) of neurons observed across sessions 
% This should be in another array (Uber_2pad_animal array)

        cellDepth = [];
        dF = []; % F should be (length(neuind),timepoints), calculated by Fneu - neuropilCoefficient * Fneuropil. dF/F_0.
        tpmTime = []; % each frame has it's own time. plane from top to bottom
        spk = [];
        
    end
    
    properties (Dependent = true)
        protractionTouchDTheta
        protractionTouchDPhi
        protractionTouchDKappaV
        protractionTouchDKappaH
    end
    
    methods (Access = public)
        function obj = Uber_2pad(b,wf,ca,trialNum)
            obj.trialNum = trialNum;
            obj.leftLickTime = b.beamBreakTimesLeft;
            obj.rightLickTime = b.beamBreakTimesRight;
            obj.answerLickTime = b.answerLickTime;
            if ~isempty(b.drinkingTime)
                if b.trialType(1) == 'r'
                    drinkingInd = find(b.beamBreakTimesRight > b.drinkingTime(1) & b.beamBreakTimesRight < b.drinkingTime(end));
                    if ~isempty(drinkingInd)
                        obj.drinkingOnsetTime = b.beamBreakTimesRight(drinkingInd(1));
                    end
                elseif b.trialType(1) == 'l'
                    drinkingInd = find(b.beamBreakTimesLeft > b.drinkingTime(1) & b.beamBreakTimesLeft < b.drinkingTime(end));
                    if ~isempty(drinkingInd)
                        obj.drinkingOnsetTime = b.beamBreakTimesLeft(drinkingInd(1));
                    end
                end
            end
            obj.poleUpOnsetTime = b.poleUpOnsetTime;
            obj.poleDownOnsetTime = b.poleDownOnsetTime;
            obj.response = b.trialCorrect; % 1 hit, 0 wrong, -1 miss
            obj.angle = b.servoAngle;
            obj.distance = b.motorDistance;
            obj.position = b.motorApPosition;
            obj.taskTarget = b.taskTarget;
            obj.distractor = b.distractor;
            obj.trialType = b.trialType;
            
            obj.protractionTouchChunks = wf.protractionTFchunks;
            obj.retractionTouchChunks = wf.retractionTFchunks;
            obj.protractionTouchDuration = wf.protractionTouchDuration;
            obj.protractionTouchSlideDistance = cellfun(@(x) max(x) - x(1), wf.protractionSlide);
%             obj.protractionTouchChunksByWhisking = wf.protractionTFchunksByWhisking;
%             obj.retractionTouchChunksByWhisking = wf.retractionTFchunksByWhisking;

            
            obj.poleMovingTime = (wf.poleMovingFrames-1)*wf.framePeriodInSec;
            obj.poleUpTime = (wf.poleUpFrames-1)*wf.framePeriodInSec;
            obj.whiskerTime = wf.time;
            obj.theta = wf.theta;
            obj.phi = wf.phi;
            obj.kappaH = wf.kappaH;
            obj.kappaV = wf.kappaV;
            obj.nof = wf.nof;
            obj.frameDuration = wf.framePeriodInSec;
            
            obj.planes = ca.planes;
            obj.neuindSession = ca.cellNums;
            obj.cellDepth = ca.cellDepth;
            obj.dF = ca.dF;
            obj.spk = ca.spk;
            obj.tpmTime = cell(length(ca.time),1);
            for i = 1 : length(obj.tpmTime)
                obj.tpmTime{i} = ca.time{i};
            end
%             obj.tpmTime = ca.time{1};
            
        end

    end
    
    methods % Dependent property methods; cannot have attributes.
        function value = get.protractionTouchDTheta(obj)
            value = cellfun(@(x) max(obj.theta(x)) - obj.theta(x(1)), obj.protractionTouchChunks);
        end
        function value = get.protractionTouchDPhi(obj)
            if ~isempty(obj.protractionTouchChunks)
                value = zeros(1, length(obj.protractionTouchChunks));                
                for i = 1 : length(value)
                    [~, maxInd] = max(abs(obj.phi(obj.protractionTouchChunks{i}) - obj.phi(obj.protractionTouchChunks{i}(1))));
                    value(i) = obj.phi(obj.protractionTouchChunks{i}(maxInd)) - obj.phi(obj.protractionTouchChunks{i}(1));
                end
            else
                value = [];
            end
        end
        function value = get.protractionTouchDKappaV(obj)
            if ~isempty(obj.protractionTouchChunks)
                value = zeros(1, length(obj.protractionTouchChunks));                
                for i = 1 : length(value)
                    [~, maxInd] = max(abs(obj.kappaV(obj.protractionTouchChunks{i}) - obj.kappaV(obj.protractionTouchChunks{i}(1))));
                    value(i) = obj.kappaV(obj.protractionTouchChunks{i}(maxInd)) - obj.kappaV(obj.protractionTouchChunks{i}(1));
                end
            else
                value = [];
            end

        end
        function value = get.protractionTouchDKappaH(obj)
            if ~isempty(obj.protractionTouchChunks)
                value = zeros(1, length(obj.protractionTouchChunks));                
                for i = 1 : length(value)
                    [~, maxInd] = max(abs(obj.kappaH(obj.protractionTouchChunks{i}) - obj.kappaH(obj.protractionTouchChunks{i}(1))));
                    value(i) = obj.kappaH(obj.protractionTouchChunks{i}(maxInd)) - obj.kappaH(obj.protractionTouchChunks{i}(1));
                end
            else
                value = [];
            end

        end
    end
    
end