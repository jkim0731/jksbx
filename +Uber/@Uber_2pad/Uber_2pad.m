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
        
        % from the Whisker Behavior
        poleMovingTime = []; % in sec
        poleUpTime = []; % in sec
        whiskerTime = []; % in sec, all frames
        theta = [];
        nof = 0;
        frameDuration = 1/311;
        
        % from W3_2pad
        phi = [];
        kappaH = [];
        kappaV = [];
        
        % from WL_2pad
%         touchChunks = {}; % touch chunks, each converted into time (in sec)
        protractionTouchChunks = {};
        retractionTouchChunks = {};
        
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
    end
    
    methods (Access = public)
        function obj = Uber_2pad(b,w3,ca,trialNum)
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
            
            obj.poleMovingTime = (w3.poleMovingFrames-1)*w3.framePeriodInSec;
            obj.poleUpTime = (w3.poleUpFrames-1)*w3.framePeriodInSec;
            obj.whiskerTime = w3.time;
            obj.theta = w3.theta;
            obj.phi = w3.phi;
            obj.kappaH = w3.kappaH;
            obj.kappaV = w3.kappaV;
            obj.nof = w3.nof;
            obj.frameDuration = w3.framePeriodInSec;
            
            if ~isempty(w3.protractionTouchChunks)
                obj.protractionTouchChunks = cellfun(@(x) (x-1) * w3.framePeriodInSec, w3.protractionTouchChunks, 'uniformoutput', false);
            end
            if ~isempty(w3.retractionTouchChunks)
                obj.retractionTouchChunks = cellfun(@(x) (x-1) * w3.framePeriodInSec, w3.retractionTouchChunks, 'uniformoutput', false);
            end
            
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
        
    end
    
end