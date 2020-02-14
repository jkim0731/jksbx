function parfun_glmnet_whisker_perCell(info, tempSpike, tempPredictor, stratificationGroups, tindCell)

%% 
localDir = info.localDir;
mouse = info.mouse;
session = info.session;
cellInd = info.ci;
repetition = info.repetition;
posShift = info.posShift;
numFrames = info.numFrames;
trialNums = info.trialNums;
numCell = info.numCell;
glmnetOpt = info.glmnetOpt;
lambdaCV = info.lambdaCV;

%%
savedFnList = dir([sprintf('%sJK%03dS%02dci%04d_save_R',localDir,mouse,session,cellInd), '*']);
if isempty(savedFnList)
    Rinds = 1 : repetition;
else
    savedFnRinds = zeros(length(savedFnList),1);
    for i = 1 : length(savedFnList)
        [~, a] =strtok(savedFnList(i).name(end-10:end-4), 'R');
        savedFnRinds(i) = str2double(a(2:end));
    end
    if length(savedFnList) >= 10
        return
    end
    totalNum = 10 + length(savedFnList);
    Rinds = setdiff(1:totalNum, savedFnRinds);
end

parfor i = 1 : repetition
    saveFn = sprintf('%sJK%03dS%02dci%04d_save_R%02d',localDir,mouse,session,cellInd, Rinds(i));    
    fprintf('Mouse JK%03d session S%02d Running cell %d/%d (Loop %d) \n', mouse, session, cellInd, numCell, i);

    pfNumFrames = numFrames;
    pftindCell = tindCell;
    pfSpike = tempSpike;
    spkMedian = median(cellfun(@(x) sum(x), pfSpike));
    spkNumGroup = cell(2,1);
    tempInd = find(cellfun(@(x) sum(x) <= spkMedian, pfSpike));
    spkNumGroup{1} = trialNums(pftindCell(tempInd));
    spkNumGroup{2} = trialNums(pftindCell( setdiff(1:length(pftindCell),tempInd) ));

    ptouchGroup = stratificationGroups{1};
    choiceGroup = stratificationGroups{2};
    angleGroup = stratificationGroups{3};
    distanceGroup = stratificationGroups{4};
    tempTestTn = [];
    for pti = 1 : length(ptouchGroup)
        for ci = 1 : length(choiceGroup)
            for ai = 1 : length(angleGroup)
                for di = 1 : length(distanceGroup)
                    for spki = 1 : length(spkNumGroup)                                    
                        tempTn = intersect(ptouchGroup{pti}, intersect(choiceGroup{ci}, intersect(angleGroup{ai}, intersect(distanceGroup{di}, spkNumGroup{spki}))));
                        if ~isempty(tempTn)
                            tempTn = tempTn(randperm(length(tempTn)));
                            if length(tempTn) > 5
                                tempTestTn = [tempTestTn; tempTn(1:round(length(tempTn)*0.3))];
                            elseif length(tempTn) > 1
                                tempTestTn = [tempTestTn; tempTn(1:round(length(tempTn)*0.5))];                                        
                            end
                        end
                    end
                end
            end
        end
    end
    %
    [~,testInd] = ismember(tempTestTn, trialNums);

    tempTrainingTn = setdiff(trialNums, tempTestTn);
    [~,trainingInd] = ismember(tempTrainingTn, trialNums);


    iTrain = intersect(pftindCell, trainingInd);
    iTest = intersect(pftindCell, testInd);

    testTn = trialNums(testInd);
    trainingTn = trialNums(trainingInd);

    ratioi = length(iTest)/length(iTrain);

    if size(pfNumFrames,2) == 1
        pfNumFrames = pfNumFrames';
    end
    if size(pftindCell,2) == 1
        pftindCell = pftindCell';
    end
    trainingPredictorBinary = cell2mat(arrayfun(@(x,y) ones(1,x+posShift*2) * ismember(y, iTrain), pfNumFrames, pftindCell, 'un', 0));
    testPredictorBinary = cell2mat(arrayfun(@(x,y) ones(1,x+posShift*2) * ismember(y, iTest), pfNumFrames, pftindCell, 'un', 0));

    if sum(trainingPredictorBinary .* testPredictorBinary)
        error('Intersection between trainingPredictorInd and testPredictorInd')
    elseif sum(trainingPredictorBinary + testPredictorBinary) ~= size(tempPredictor,1)
        error('Number of total frames mismatch')
    end

    ratioInd = length(find(testPredictorBinary)) / length(find(trainingPredictorBinary));

    trainingInput = tempPredictor(find(trainingPredictorBinary),:);
    testInput = tempPredictor(find(testPredictorBinary),:);

    if size(pfSpike,2) == 1
        pfSpike = pfSpike';
    end
    spkTrain = cell2mat(cellfun(@(x) [nan(1,posShift), x, nan(1,posShift)], pfSpike(ismember(iTrain, pftindCell)),'uniformoutput',false));
    finiteIndTrain = intersect(find(isfinite(spkTrain)), find(isfinite(sum(trainingInput,2))));
    input = trainingInput(finiteIndTrain,:);
    spkTrain = spkTrain(finiteIndTrain)';
    tic
    cv = cvglmnet(input, spkTrain, 'poisson', glmnetOpt, [], lambdaCV);
    cellTime = toc;
    %% survived coefficients
    fitLambda = cv.lambda_1se;
    iLambda = find(cv.lambda == cv.lambda_1se);
    fitCoeffs = [cv.glmnet_fit.a0(iLambda);cv.glmnet_fit.beta(:,iLambda)];

    %% test
    spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x, nan(1,posShift)], pfSpike(ismember(iTest, pftindCell)),'uniformoutput',false));
    spkTest = spkTest';
    finiteIndTest = intersect(find(isfinite(spkTest)), find(isfinite(sum(testInput,2))));
    spkTest = spkTest(finiteIndTest)';
    %% (1) if the full model is significant
    model = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,:)]*[cv.glmnet_fit.a0(iLambda); cv.glmnet_fit.beta(:,iLambda)]);
    mu = mean(spkTest); % null poisson parameter
    nullLogLikelihood = sum(log(poisspdf(spkTest,mu)));
    fullLogLikelihood = sum(log(poisspdf(spkTest',model)));
    saturatedLogLikelihood = sum(log(poisspdf(spkTest,spkTest)));
    devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);
    fitDeviance = devianceFullNull;
    [fitCorrelation, fitCorrPval] = corr(spkTest', model);            
    fitDevExplained = 1 - (saturatedLogLikelihood - fullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
    fitCvDev = cv.glmnet_fit.dev(iLambda);
 
    %%
    parsave(saveFn, fitLambda, fitCoeffs, fitDeviance, fitCorrelation, fitCorrPval, fitDevExplained, fitCvDev, testTn, trainingTn, ratioi, ratioInd, ...
        ptouchGroup, choiceGroup, angleGroup, distanceGroup, spkNumGroup, cellTime);
end % end of parfor

end % end of function parfun_


function parsave(fn, fitLambda, fitCoeffs, fitDeviance, fitCorrelation, fitCorrPval, fitDevExplained, fitCvDev, testTn, trainingTn, ratioi, ratioInd, ...
    ptouchGroup, choiceGroup, angleGroup, distanceGroup, spkNumGroup, cellTime)
    save(fn, 'fit*', 'testTn', 'trainingTn', 'ratio*', '*Group', 'cellTime')
end
