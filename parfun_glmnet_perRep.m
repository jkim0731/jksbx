function parfun_glmnet_perRep(info, spikeAll, allPredictors, stratificationGroups, tindcellAll, planeIndAll, cindAll, ri)

%% 
localDir = info.localDir;
mouse = info.mouse;
session = info.session;
posShift = info.posShift;
trialNums = info.trialNums;
numCell = info.numCell;
glmnetOpt = info.glmnetOpt;
lambdaCV = info.lambdaCV;

%%
parfor ci = 1 : numCell
    saveFn = sprintf('%sJK%03dS%02dci%04d_save_R%02d',localDir,mouse,session,ci, ri);  
    if ~strcmp(dir(saveFn), saveFn) % Run only when it's not done yet
        fprintf('Mouse JK%03d session S%02d Loop %d: Running cell %d/%d \n', mouse, session, ri, ci, numCell);

        tempPredictor = allPredictors{planeIndAll(ci)};
        tindCell = tindcellAll{ci};
        cind = cindAll(ci);
        tempSpike = cellfun(@(x) x(cind,:), spikeAll(tindCell), 'un', 0);
    %     tempPlane = mod(planeIndAll(ci),4);
    %     if tempPlane ==0
    %         tempPlane = 4;
    %     end
    %     numFrames = cellfun(@(x) length(x.tpmTime{tempPlane}), u.trials(tindCell));
        numFrames = cellfun(@length, tempSpike);
        totalNumSpikeFrames = sum(cellfun(@(x) length(find(x)), tempSpike));
        if totalNumSpikeFrames < 10 % from the results of d190415_error_and_nonfit_cells_threshold_finding
            fprintf('Abort because of too little spike frames in cell index %d\n', cellnum)
        else
            spkMedian = median(cellfun(@(x) sum(x), tempSpike));
            spkNumGroup = cell(2,1);
            tempInd = find(cellfun(@(x) sum(x) <= spkMedian, tempSpike));
            spkNumGroup{1} = trialNums(tindCell(tempInd));
            spkNumGroup{2} = trialNums(tindCell( setdiff(1:length(tindCell),tempInd) ));

            ptouchGroup = stratificationGroups{1};
            choiceGroup = stratificationGroups{2};
            angleGroup = stratificationGroups{3};
            distanceGroup = stratificationGroups{4};
            tempTestTn = [];
            for pti = 1 : length(ptouchGroup)
                for choicei = 1 : length(choiceGroup)
                    for ai = 1 : length(angleGroup)
                        for di = 1 : length(distanceGroup)
                            for spki = 1 : length(spkNumGroup)                                    
                                tempTn = intersect(ptouchGroup{pti}, intersect(choiceGroup{choicei}, intersect(angleGroup{ai}, intersect(distanceGroup{di}, spkNumGroup{spki}))));
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


            iTrain = intersect(tindCell, trainingInd);
            iTest = intersect(tindCell, testInd);

            testTn = trialNums(testInd);
            trainingTn = trialNums(trainingInd);

            ratioi = length(iTest)/length(iTrain);

            if size(numFrames,2) == 1
                numFrames = numFrames';
            end
            if size(tindCell,2) == 1
                tindCell = tindCell';
            end
            trainingPredictorBinary = cell2mat(arrayfun(@(x,y) ones(1,x+posShift*2) * ismember(y, iTrain), numFrames, tindCell, 'un', 0));
            testPredictorBinary = cell2mat(arrayfun(@(x,y) ones(1,x+posShift*2) * ismember(y, iTest), numFrames, tindCell, 'un', 0));

            if sum(trainingPredictorBinary .* testPredictorBinary)
                error('Intersection between trainingPredictorInd and testPredictorInd')
            elseif sum(trainingPredictorBinary + testPredictorBinary) ~= size(tempPredictor,1)
                error('Number of total frames mismatch')
            end

            ratioInd = length(find(testPredictorBinary)) / length(find(trainingPredictorBinary));

            trainingInput = tempPredictor(find(trainingPredictorBinary),:);
            testInput = tempPredictor(find(testPredictorBinary),:);

            if size(tempSpike,2) == 1
                tempSpike = tempSpike';
            end
            spkTrain = cell2mat(cellfun(@(x) [nan(1,posShift), x, nan(1,posShift)], tempSpike(find(ismember(tindCell,iTrain))),'uniformoutput',false));
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
            spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x, nan(1,posShift)], tempSpike(find(ismember(tindCell,iTest))),'uniformoutput',false));
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
        end % end of if totalNumSpikeFrames < 10
    end % end of if ~strcmp
end % end of parfor

end % end of function parfun_


function parsave(fn, fitLambda, fitCoeffs, fitDeviance, fitCorrelation, fitCorrPval, fitDevExplained, fitCvDev, testTn, trainingTn, ratioi, ratioInd, ...
    ptouchGroup, choiceGroup, angleGroup, distanceGroup, spkNumGroup, cellTime)
    save(fn, 'fit*', 'testTn', 'trainingTn', 'ratio*', '*Group', 'cellTime')
end
