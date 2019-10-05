function validOut = tune_prop_test(spkVal)

thresholdCategory = 0.05;
angles = 45:15:135;
anovactype = 'hsd';

multimodal = 0;
unimodalSingle = 0;
unimodalBroad = 0;
categorical = 0;
leaveOneOut = 0;

spkAnovaVal = cell2mat(spkVal);
if ~isempty(find(isnan(spkAnovaVal)))
    error('nan values')
end
groupAnova = zeros(size(spkAnovaVal));
angleLengths = [0;cumsum(cellfun(@length, spkVal))];
for ai = 1 : length(angles)
    groupAnova(angleLengths(ai)+1:angleLengths(ai+1)) = deal(ai);
end


[spkAnovaP, ~, spkAnovaStat] = anova1(spkAnovaVal, groupAnova, 'off');
spkPairComp = multcompare(spkAnovaStat, 'Ctype', anovactype, 'Display', 'off');
spkMeans = spkAnovaStat.means;
spkStd = cellfun(@std, spkVal);
%%
[~,testInd] = min(spkAnovaStat.means)
tempH = cellfun(@(x) ttest2(x,spkVal{testInd}), spkVal);
%%

tempH = cellfun(@(x) ttest(x), spkVal);
%%
tempH(isnan(tempH)) = deal(0);
sigInd = find(tempH); % significant indices

[~, maxind] = max(abs(spkMeans(sigInd)));
tunedAngleInd = sigInd(maxind);
tunedAngle = angles(tunedAngleInd);



% Categorization
ind__1 = find(spkPairComp(:,1) == tunedAngleInd);
ind__2 = find(spkPairComp(:,2) == tunedAngleInd);
testInd = union(ind__1, ind__2);
insigDiffInd = find(spkPairComp(testInd,6) >= thresholdCategory);
sigDiffInd = find(spkPairComp(testInd,6) < thresholdCategory);
temp = spkPairComp(testInd(insigDiffInd),1:2);
insigDiffIndGroup = unique(temp(:)); % sorted. Include tunedAngleInd, except when there's nothing

if isempty(insigDiffIndGroup)
    unimodalSingle = 1;
else
    broadInd = intersect(sigInd,insigDiffIndGroup);
    if length(broadInd) < 2
        unimodalSingle = 1;
    else                            
        broadNum = 1;
        for tunei = tunedAngleInd-1:-1:1
            if ismember(tunei, broadInd)
                broadNum = broadNum + 1;
            else 
                break
            end
        end
        for tunei = tunedAngleInd+1:length(angles)
            if ismember(tunei, broadInd)
                broadNum = broadNum + 1;
            else
                break
            end
        end
        if broadNum == length(broadInd)
            unimodalBroad = 1;
            % if broad, then it can be a categorical
            center = (length(angles)+1) / 2;
            compInd = union(find(spkPairComp(:,1) == tunedAngleInd), find(spkPairComp(:,2) == tunedAngleInd));
            indMat = spkPairComp(compInd,1:2);
            if tunedAngleInd < center
                withinInd = unique(mod( setdiff( find(indMat < center), find(indMat == tunedAngleInd) ) , size(indMat,1)));
                withinInd(withinInd==0) = size(indMat,1);
                betweenInd = unique(mod( find(indMat > center) , size(indMat,1) ));
                betweenInd(betweenInd==0) = size(indMat,1);
            else
                withinInd = unique(mod( setdiff( find(indMat > center), find(indMat == tunedAngleInd) ) , size(indMat,1)));
                withinInd(withinInd==0) = size(indMat,1);
                betweenInd = unique(mod( find(indMat < center) , size(indMat,1) ));
                betweenInd(betweenInd==0) = size(indMat,1);
            end
            if isempty(find(spkPairComp(compInd(withinInd),6) < thresholdCategory, 1)) && ... % nothing within the same half is different from the max ind
                    isempty(find(spkPairComp(compInd(betweenInd),6) >= thresholdCategory, 1)) % nothing between different half is same with the max ind
                categorical = 1; % categorical (>= 90 or <= 90)
            end
        else
            multimodal = 1;
        end
    end
    temp = spkPairComp(testInd(sigDiffInd),1:2);
    sigIndGroup = setdiff(temp(:), tunedAngleInd); % exclude tunedAngleInd. Any index that is significantly different from the tuned angle index.
    if ~isempty(find(diff(insigDiffIndGroup)>1,1))
        if sum(tempH(insigDiffIndGroup))>1 % to exclude tuned angle
            multimodal = 1; % multimodal. Including bipolar.
        end
        if length(sigIndGroup) == 1 && ... % only one bin is significantly different from the tuned bin. (can't be larger in response because of the way tuned bin is defined)
                all(tempH(insigDiffIndGroup)) % and all insignicant indices are different from 0
            leaveOneOut = 1 ; % leave-one-out. Part of multimodal in definition.
        end
    end
end

validOut.multimodal = multimodal;
validOut.unimodalSingle = unimodalSingle;
validOut.unimodalBroad = unimodalBroad;
validOut.categorical = categorical;
validOut.leaveOneOut = leaveOneOut;
validOut.tunedAngle = tunedAngle;