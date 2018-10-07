classdef Uber_2pad < handle
    properties
        trialNum = [];
        
        % from the Task (Solo)
        leftLickTime = []; % in sec
        rightLickTime = []; % in sec
        answerLickTime = []; % in sec
        drinkingTime = []; % in sec
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
        %
        % from thTouchFrames, for now 2018/05/01 JK. Needs to implement
        % thresholds from kappa and theta
        touchChunks = {}; % touch chunks, each converted into time (in sec)
        %
        % whisker variables, not for now. Need 3D reconstruction first
%         kappa = cell(1,2); % filling up all frames. NaN when there is no trace at that frame. if one whisker is lost, the other is treated as lost as well. (both have same entries of NaN)
%         theta = cell(1,2); % same as kappa
        
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
        tpmTime = {}; % each frame has it's own time. plane from top to bottom
        
        active = []; % 0 (inactive) or 1 (active) for each cell
        
    end
    
    properties (Dependent = true)
    end
    
    methods (Access = public)
        function obj = Uber_2pad(b,wl,ca,trialNum)        
                    %
                    %
                    %
                    %
                    %
                    %
                    %
                    
                    if wl.framePeriodInSec > 1
                        wl.framePeriodInSec = wl.framePeriodInSec/wl.framePeriodInSec/wl.framePeriodInSec;
                    end
                    %
                    %
                    %
                    %
                    %
                    %
                    % Temporary solution for mistake in framePeriodInSec
                    % with JK05X 2018/10/04
                    %
                    %
            
            obj.trialNum = trialNum;
            obj.leftLickTime = b.beamBreakTimesLeft;
            obj.rightLickTime = b.beamBreakTimesRight;
            obj.answerLickTime = b.answerLickTime;
            obj.drinkingTime = b.drinkingTime;
            obj.poleUpOnsetTime = b.poleUpOnsetTime;
            obj.poleDownOnsetTime = b.poleDownOnsetTime;
            obj.response = b.trialCorrect; % 1 hit, 0 wrong, -1 miss
            obj.angle = b.servoAngle;
            obj.distance = b.motorDistance;
            obj.position = b.motorApPosition;
            obj.taskTarget = b.taskTarget;
            obj.distractor = b.distractor;
            obj.trialType = b.trialType;
            
            obj.poleMovingTime = (wl.poleMovingFrames-1)*wl.framePeriodInSec;
            obj.poleUpTime = (wl.poleUpFrames-1)*wl.framePeriodInSec;
            obj.whiskerTime = (0:wl.nof-1)*wl.framePeriodInSec;
            if ~isempty(wl.protractionTouchFrames) % only for protraction, for now 2018/10/02 JK
                tf = wl.protractionTouchFrames;
                
                % single frame correction. Just for now. It should correct
                % 101->111 first, and then 010->000
                binaryFrames = zeros(max(tf),1);
                binaryFrames(tf) = 1;
                flipInds = find( diff(diff(binaryFrames)) > 1) + 1;
                binaryFrames(flipInds) = 1 - binaryFrames(flipInds);                
                tf = find(binaryFrames);
                % 010->000
                binaryFrames = zeros(max(tf),1);
                binaryFrames(tf) = 1;
                flipInds = find( diff(diff(binaryFrames)) < -1) + 1;
                binaryFrames(flipInds) = 1 - binaryFrames(flipInds);                
                tf = find(binaryFrames);
                
                % making chunks. It worked in wl function, but just in
                % case...
                if size(tf,1) > size(tf,2)
                    tf = tf';
                end
                chunkPoints = [1, find(diff(tf)>1) + 1, length(tf)+1]; % +1 because of diff. first 1 for the first chunk. So this is actually start points of each chunk.
                chunks = cell(1,length(chunkPoints)-1); % 
                for i = 1 : length(chunks)
                    chunks{i} = [tf(chunkPoints(i) : chunkPoints(i+1)-1)]';
                end
                
                obj.touchChunks = cell(length(chunks),1);
                for i = 1 : length(chunks)
                    obj.touchChunks{i} = (chunks{i}-1) * wl.framePeriodInSec;
                end
            else
                obj.touchChunks = {};
            end
%             topind = round(wl.time{1}*wl.framePeriodInSec) + 1;
%             frontind = round(wl.time{2}*wl.framePeriodInSec) + 1;
%             nonanind = intersect(topind, frontind);
%             indfin{1} = find(ismember(topind, nonanind));
%             indfin{2} = find(ismember(frontind, nonanind));
            % whisker variables, not for now. Need 3D reconstruction first
%             for i = 1 : 2 % 2pad, 2 views only (top and front)
%                 obj.kappa{i} = nan(wl.nof,1);
%                 obj.theta{i} = nan(wl.nof,1);
%                 obj.kappa{i}(nonanind) = wl.kappa{i}(indfin{i});
%                 obj.theta{i}(nonanind) = wl.theta{i}(indfin{i});
%             end
            
            obj.planes = ca.planes;
            obj.neuindSession = ca.cellNums;
            obj.cellDepth = ca.cellDepth;
            obj.dF = ca.dF;
            obj.tpmTime = cell(length(ca.time),1);
            for i = 1 : length(obj.tpmTime)
                obj.tpmTime{i} = ca.time{1};
            end            
        end
        
    end
    
    methods % Dependent property methods; cannot have attributes.
        
    end
    
end