% save(savefnResultRe, 'fit*', 'allPredictors', '*InputMat', 'indPartial', '*Group', '*Tn', 'lambdaCV', '*Opt', 'done', '*Re', 'remainingCell', 'pThreshold*', '*Shift');

%%

poolobj = gcp('nocreate');
if poolobj.SpmdEnabled == 0
    error('SpmdEnabled changed to false')
end

if ~exist('pThresholdPartial', 'var')
    pThresholdPartial = 0.05;
end
if ~exist('pThresholdNull', 'var')
    pThresholdNull = 0.05;
end

if ~exist('posShift', 'var')
    posShift = 4;
end
if ~exist('negShift', 'var')
    negShift = 2;
end

mouse = 52;
session = 26;
repeat = 5;
restartingNum = 1;
glmPar = true;
savefnResult = sprintf('glmResponseType_JK%03dS%02d_m45_R%02d',mouse, session, repeat);

savefnResultRe = [savefnResult, '_02'];
% errorCellSession = [26, 82, 83, 87, 89, 104, 119, 120, 130, 246, 298, 319, 366, 391, 409, 411, 593, 615, 834, 1822]; % JK039 S23
% errorCellSession = [139,152, 163,259,453,1726]; % JK039 S24
errorCellSession = [160]; % JK052 S26
previousDone = done(find(done));

numCell = length(cIDAll);
remainingCell = setdiff(1 : numCell, previousDone);


startedRe = zeros(length(remainingCell),1);
fitLambdaRe = zeros(length(remainingCell),1);
fitCoeffsRe = cell(length(remainingCell),1);
fitDevExplainedRe = zeros(length(remainingCell),1);
ratioIndRe = zeros(length(remainingCell),1);
ratioiRe = zeros(length(remainingCell),1);
testTnRe = cell(length(remainingCell),1);
trainingTnRe = cell(length(remainingCell),1);

fitCvDevRe = zeros(length(remainingCell),1);
fitResultsRe = zeros(length(remainingCell),6);
doneRe = zeros(length(remainingCell),1);
cellTimeRe = zeros(length(remainingCell),1);

tindcellAll = cell(numCell,1);
cindAll = zeros(numCell,1);
planeIndAll = zeros(numCell,1);
for i = 1 : numCell
    tindcellAll{i} = find(cellfun(@(x) ismember(cIDAll(i), x.neuindSession), u.trials));
    cindAll(i) = find(u.trials{tindcellAll{i}(1)}.neuindSession == cIDAll(i));
    planeIndAll(i) = floor(cIDAll(i)/1000);
end
spikeAll = cellfun(@(x) x.spk, u.trials, 'uniformoutput', false);

poolobj = gcp('nocreate');
if poolobj.SpmdEnabled == 0
    error('SpmdEnabled turned to false at #2');
end

% clear u


poolobj = gcp('nocreate');
if poolobj.SpmdEnabled == 0
    error('SpmdEnabled changed to false')
end



if glmPar
    for cellnumInd = restartingNum : length(remainingCell)
        cellnum = remainingCell(cellnumInd);
        if ismember(cellnum, errorCellSession)
            fprintf('Abort because of previous error in cell index %d\n', cellnum)
        else

            startedRe(cellnumInd) = cellnum;
            cellTimeStart = tic;
            fprintf('Mouse JK%03d session S%02d Loop %d: Running cell %d/%d \n', mouse, session, repeat, cellnum, numCell);                

            cind = cindAll(cellnum);
            tindCell = tindcellAll{cellnum};

            totalNumSpikeFrames = sum(cellfun(@(x) length(find(x(cind,:))), spikeAll(tindCell)'));
            if totalNumSpikeFrames < 10 % from the results of d190415_error_and_nonfit_cells_threshold_finding
                fprintf('Abort because of too little spike frames in cell index %d\n', cellnum)
            else
                spkMedian = median(cellfun(@(x) sum(x(cind,:)), spikeAll(tindCell)'));                
                spkNumGroup = cell(2,1);
                spkNumGroup{1} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) sum(x.spk(cind,:)) <= spkMedian, u.trials(tindCell)) )));
                spkNumGroup{2} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) sum(x.spk(cind,:)) >  spkMedian, u.trials(tindCell)) )));

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
                totalTn = u.trialNums;
                [~,testInd] = ismember(tempTestTn, totalTn);

                tempTrainingTn = setdiff(totalTn, tempTestTn);
                [~,trainingInd] = ismember(tempTrainingTn, totalTn);

                iTrain = intersect(tindCell, trainingInd);
                iTest = intersect(tindCell, testInd);

                testTnRe{cellnumInd} = u.trialNums(testInd);
                trainingTnRe{cellnumInd} = u.trialNums(trainingInd);

                ratioiRe(cellnumInd) = length(iTest)/length(iTrain);

                planeInd = planeIndAll(cellnum);
                plane = mod(planeInd,4);
                if plane==0
                    plane = 4;
                end
                trainingPredictorInd = cell2mat(cellfun(@(x) (ones(1,length(x.tpmTime{plane})+posShift*2)) * ismember(x.trialNum, tempTrainingTn), u.trials(tindCell)','uniformoutput',false));
                testPredictorInd = cell2mat(cellfun(@(x) (ones(1,length(x.tpmTime{plane})+posShift*2)) * ismember(x.trialNum, tempTestTn), u.trials(tindCell)','uniformoutput',false));

                if (trainingPredictorInd .* testPredictorInd)
                    error('Intersection between trainingPredictorInd and testPredictorInd')
                elseif sum(trainingPredictorInd + testPredictorInd) ~= size(allPredictors{planeInd},1)
                    error('Number of total frames mismatch')
                end

                ratioIndRe(cellnumInd) = length(find(testPredictorInd)) / length(find(trainingPredictorInd));


                planeInd = planeIndAll(cellnum);
                trainingPredictorInd = cell2mat(cellfun(@(x) (ones(1,length(x.tpmTime{plane})+posShift*2)) * ismember(x.trialNum, tempTrainingTn), u.trials(tindCell)','uniformoutput',false));
                testPredictorInd = cell2mat(cellfun(@(x) (ones(1,length(x.tpmTime{plane})+posShift*2)) * ismember(x.trialNum, tempTestTn), u.trials(tindCell)','uniformoutput',false));

                trainingInput = allPredictors{planeInd}(find(trainingPredictorInd),:);
                testInput = allPredictors{planeInd}(find(testPredictorInd),:);

                spkTrain = cell2mat(cellfun(@(x) [nan(1,posShift), x(cind,:), nan(1,posShift)], spikeAll(iTrain)','uniformoutput',false));
                finiteIndTrain = intersect(find(isfinite(spkTrain)), find(isfinite(sum(trainingInput,2))));
                input = trainingInput(finiteIndTrain,:);
                spkTrain = spkTrain(finiteIndTrain)';

                cv = cvglmnet(input, spkTrain, 'poisson', glmnetOpt, [], lambdaCV);
                %% survived coefficients
                fitLambdaRe(cellnumInd) = cv.lambda_1se;
                iLambda = find(cv.lambda == cv.lambda_1se);
                fitCoeffsRe{cellnumInd} = [cv.glmnet_fit.a0(iLambda);cv.glmnet_fit.beta(:,iLambda)];
                coeffInds = find(cv.glmnet_fit.beta(:,iLambda));                

                %% test            
                spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x(cind,:), nan(1,posShift)], spikeAll(iTest)','uniformoutput',false));
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
                fitDevExplainedRe(cellnumInd) = 1 - (saturatedLogLikelihood - fullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
                fitCvDevRe(cellnumInd) = cv.glmnet_fit.dev(iLambda);
                cellTimeRe(cellnumInd) = toc(cellTimeStart);
                
                doneRe(cellnumInd) = cellnum;
            end
        end
    end % end of for cellnum
else
    parfor cellnumInd = 1 : length(remainingCell)
        cellnum = remainingCell(cellnumInd);
        if ismember(cellnum, errorCellSession)
            fprintf('Abort because of previous error in cell index %d\n', cellnum)
        else
        
            fitCoeffInd = zeros(1,6);
            startedRe(cellnumInd) = cellnum;
            cellTimeStart = tic;
            fprintf('Mouse JK%03d session S%02d Loop %d: Running cell %d/%d \n', mouse, session, repeat, cellnum, numCell);                

            cind = cindAll(cellnum);
            tindCell = tindcellAll{cellnum};

            totalNumSpikeFrames = sum(cellfun(@(x) length(find(x(cind,:))), spikeAll(tindCell)'));
            if totalNumSpikeFrames < 10 % from the results of d190415_error_and_nonfit_cells_threshold_finding
                fprintf('Abort because of too little spike frames in cell index %d\n', cellnum)
            else
                spkMedian = median(cellfun(@(x) sum(x(cind,:)), spikeAll(tindCell)'));                
                spkNumGroup = cell(2,1);
                spkNumGroup{1} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) sum(x.spk(cind,:)) <= spkMedian, u.trials(tindCell)) )));
                spkNumGroup{2} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) sum(x.spk(cind,:)) >  spkMedian, u.trials(tindCell)) )));

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
                totalTn = u.trialNums;
                [~,testInd] = ismember(tempTestTn, totalTn);

                tempTrainingTn = setdiff(totalTn, tempTestTn);
                [~,trainingInd] = ismember(tempTrainingTn, totalTn);

                testTnRe{cellnumInd} = u.trialNums(testInd);
                trainingTnRe{cellnumInd} = u.trialNums(trainingInd);

                iTrain = intersect(tindCell, trainingInd);
                iTest = intersect(tindCell, testInd);

                ratioiRe(cellnumInd) = length(iTest)/length(iTrain);

                planeInd = planeIndAll(cellnum);
                plane = mod(planeInd,4);
                if plane==0
                    plane = 4;
                end
                trainingPredictorInd = cell2mat(cellfun(@(x) (ones(1,length(x.tpmTime{plane})+posShift*2)) * ismember(x.trialNum, tempTrainingTn), u.trials(tindCell)','uniformoutput',false));
                testPredictorInd = cell2mat(cellfun(@(x) (ones(1,length(x.tpmTime{plane})+posShift*2)) * ismember(x.trialNum, tempTestTn), u.trials(tindCell)','uniformoutput',false));

                if (trainingPredictorInd .* testPredictorInd)
                    error('Intersection between trainingPredictorInd and testPredictorInd')
                elseif sum(trainingPredictorInd + testPredictorInd) ~= size(allPredictors{planeInd},1)
                    error('Number of total frames mismatch')
                end

                ratioIndRe(cellnumInd) = length(find(testPredictorInd)) / length(find(trainingPredictorInd));

                trainingInput = allPredictors{planeInd}(find(trainingPredictorInd),:);            
                testInput = allPredictors{planeInd}(find(testPredictorInd),:);

                spkTrain = cell2mat(cellfun(@(x) [nan(1,posShift), x(cind,:), nan(1,posShift)], spikeAll(iTrain)','uniformoutput',false));
                finiteIndTrain = intersect(find(isfinite(spkTrain)), find(isfinite(sum(trainingInput,2))));
                input = trainingInput(finiteIndTrain,:);
                spkTrain = spkTrain(finiteIndTrain)';

                cv = cvglmnet(input, spkTrain, 'poisson', glmnetOpt, [], lambdaCV);

                %% survived coefficients
                fitLambdaRe(cellnumInd) = cv.lambda_1se;
                iLambda = find(cv.lambda == cv.lambda_1se);
                fitCoeffsRe{cellnumInd} = [cv.glmnet_fit.a0(iLambda);cv.glmnet_fit.beta(:,iLambda)];

                %% test
                spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x(cind,:), nan(1,posShift)], spikeAll(iTest)','uniformoutput',false));
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
                fitDevExplainedRe(cellnumInd) = 1 - (saturatedLogLikelihood - fullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
                fitCvDevRe(cellnumInd) = cv.glmnet_fit.dev(iLambda);
            end

            doneRe(cellnumInd) = cellnum;
            cellTimeRe(cellnumInd) = toc(cellTimeStart);
        end

    end % end of parfor cellnum
end
%%
started(remainingCell) = startedRe;
fitLambda(remainingCell) = fitLambdaRe;
fitCoeffs(remainingCell) = fitCoeffsRe;
fitDevExplained(remainingCell) = fitDevExplainedRe;
fitCvDev(remainingCell) = fitCvDevRe;
fitResults(remainingCell,:) = fitResultsRe;
done(remainingCell) = doneRe;
cellTime(remainingCell) = cellTimeRe;
ratioInd(remainingCell) = ratioIndRe;
ratioi(remainingCell) = ratioiRe;
testTn(remainingCell) = testTnRe;
trainingTn(remainingCell) = trainingTnRe;
%             rtest(ri).fitInd = fitInd; % parameters surviving lasso in training set
%             rtest(ri).fitCoeffs = fitCoeffs;
%             rtest(ri).fitCoeffInds = fitCoeffInds; % first column is dummy
%             
%             rtest(ri).fitResults = fitResults;        
%             % fitResult(:,1) if full fitting is significant (compared to null model), 0 if not
%             % fitResult(:,2) for touchInd, compared to full fitting. if excluding touch
%             % is significantly less fit, then 1, 0 otherwise
%             % fitResult(:,3) for sound, (:,4) for reward, (:,5) for whisking, and (:,6) for licking
%     
%             rtest(ri).devExplained = devExplained;
%             rtest(ri).cvDev = cvDev;

% save(savefnResultRe, 'fit*', 'allPredictors', '*InputMat', 'indPartial', '*Group', '*Tn', 'lambdaCV', '*Opt', 'done', '*Re', 'remainingCell', 'pThreshold*', '*Shift', 'testInd', 'trainingInd', 'cIDAll', ...
%     'cellTime', 'ratio*');
save(savefnResultRe, 'fit*', 'allPredictors', 'indPartial', '*Group', 'testTn', 'trainingTn', 'lambdaCV', '*Opt', 'done', 'pThreshold*', '*Shift', 'cellTime', 'cIDAll', 'ratio*', 'restartingNum');
%%
%         end % of ri. random group selection index
% push_myphone(sprintf('GLM done for JK%03d S%02d', mouse, session))


%%
% myCluster = parcluster('local');
% delete(myCluster.Jobs)
% clear myCluster
% parpool(34, 'SpmdEnabled', true);