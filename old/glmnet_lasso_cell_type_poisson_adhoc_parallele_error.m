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
if ~exist('totalTn', 'var')
    totalTn = union(testTn,trainingTn);
end
if ~exist('cIDAll', 'var')
    cIDAll = u.cellNums;    
end
[~,testInd] = ismember(testTn, totalTn);

trainingTn = setdiff(totalTn, testTn);
[~,trainingInd] = ismember(trainingTn, totalTn);

mouse = 41;
session = 3;
restartingNum = 2;
savefnResult = sprintf('glmResponseType_JK%03dS%02d_glmnet_m14',mouse, session);

savefnResultRe = [savefnResult, '_05'];

previousDone = done(find(done));

numCell = length(cIDAll);
remainingCell = setdiff(1 : numCell, previousDone);


startedRe = zeros(length(remainingCell),1);
fitLambdaRe = zeros(length(remainingCell),1);
fitCoeffsRe = cell(length(remainingCell),1);
fitIndRe = cell(length(remainingCell),1);
fitDFRe = zeros(length(remainingCell),1);
fitDevExplainedRe = zeros(length(remainingCell),1);
fitCvDevRe = zeros(length(remainingCell),1);
fitResultsRe = zeros(length(remainingCell),6);
fitCoeffIndsRe = zeros(length(remainingCell),6);
doneRe = zeros(length(remainingCell),1);
cellTimeRe = zeros(length(remainingCell),1);


tindcellAll = cell(numCell,1);
cindAll = zeros(numCell,1);
planeIndAll = zeros(numCell,1);
iTrainAll = cell(numCell,1);
iTestAll = cell(numCell,1);
for i = 1 : numCell
    tindcellAll{i} = find(cellfun(@(x) ismember(cIDAll(i), x.neuindSession), u.trials));
    cindAll(i) = find(u.trials{tindcellAll{i}(1)}.neuindSession == cIDAll(i));
    planeIndAll(i) = floor(cIDAll(i)/1000);
    iTrainAll{i} = intersect(tindcellAll{i}, trainingInd);
    iTestAll{i} = intersect(tindcellAll{i}, testInd);
end
spikeAll = cellfun(@(x) x.spk, u.trials, 'uniformoutput', false);

% clear u


poolobj = gcp('nocreate');
if poolobj.SpmdEnabled == 0
    error('SpmdEnabled changed to false')
end

for cellnumInd = restartingNum : length(remainingCell)
    celltic = tic;
%     poolobj = gcp('nocreate')
    
%     fprintf('Is NumWorkers a prop = %d \n', isprop(poolobj,'NumWorkers'))
%     if poolobj.NumWorkers < parSetNum
%         push_myphone('Lasso GLM error')
%         error('Error occurred')
%     end    
    cellnum = remainingCell(cellnumInd);
%     if ismember(cellnum, doneRe)
%         error('Lasso GLM error')
%     end
    fitCoeffInd = nan(1,6);

%                 fprintf('Mouse JK%03d session S%02d Loop %d: Running cell %d/%d \n', mouse, session, ri,cellnum, length(u.cellNums));
    fprintf('Mouse JK%03d session S%02d: Running cell %d/%d \n', mouse, session,cellnum, numCell);
    startedRe(cellnumInd) = cellnum;
    iTrain = iTrainAll{cellnum};
    cind = cindAll(cellnum);
    planeInd = planeIndAll(cellnum);
    
    spkTrain = cell2mat(cellfun(@(x) [nan(1,posShift), x(cind,:), nan(1,posShift)], spikeAll(iTrain)','uniformoutput',false));                
    finiteIndTrain = intersect(find(isfinite(spkTrain)), find(isfinite(sum(trainingInputMat{planeInd},2))));
    input = trainingInputMat{planeInd}(finiteIndTrain,:);
    spkTrain = spkTrain(finiteIndTrain)';

    cv = cvglmnet(input, spkTrain, 'poisson', glmnetOpt, [], lambdaCV, [], true);
    
    %% survived coefficients
    fitLambdaRe(cellnumInd) = cv.lambda_1se;
    iLambda = find(cv.lambda == cv.lambda_1se);
    fitCoeffsRe{cellnumInd} = [cv.glmnet_fit.a0(iLambda);cv.glmnet_fit.beta(:,iLambda)];
    coeffInds = find(cv.glmnet_fit.beta(:,iLambda));                
%             rtest(ri).fitInd{cellnum} = coeffInds;
    fitIndRe{cellnumInd} = coeffInds;
    for i = 1 : length(indPartial)
        if sum(ismember(indPartial{i},coeffInds)>0)
            fitCoeffInd(i + 1) = 1;
        else
            fitCoeffInd(i + 1) = 0;
        end
    end

    %% test
    iTest = iTestAll{cellnum};                
    spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x(cind,:), nan(1,posShift)], spikeAll(iTest)','uniformoutput',false));
    spkTest = spkTest';
    finiteIndTest = intersect(find(isfinite(spkTest)), find(isfinite(sum(testInputMat{planeInd},2))));
    spkTest = spkTest(finiteIndTest)';
    %% (1) if the full model is significant
    fitResult = zeros(1,6);

    model = exp([ones(length(finiteIndTest),1),testInputMat{planeInd}(finiteIndTest,:)]*[cv.glmnet_fit.a0(iLambda); cv.glmnet_fit.beta(:,iLambda)]);
    mu = mean(spkTest); % null poisson parameter
    nullLogLikelihood = sum(log(poisspdf(spkTest,mu)));
    fullLogLikelihood = sum(log(poisspdf(spkTest',model)));
    saturatedLogLikelihood = sum(log(poisspdf(spkTest,spkTest)));
    devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);
    dfFullNull = length(coeffInds);                
    fitDFRe(cellnumInd) = dfFullNull;
    fitDevExplainedRe(cellnumInd) = 1 - (saturatedLogLikelihood - fullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
    fitCvDevRe(cellnumInd) = cv.glmnet_fit.dev(iLambda);

    if devianceFullNull > chi2inv(1-pThresholdNull, dfFullNull)
        fitResult(1) = 1;
        %% (2) test without each parameter (as a group)                
        for pi = 1 : 5
            if find(ismember(coeffInds, indPartial{pi}))
                if all(ismember(coeffInds, indPartial{pi}))
                    fitResult(pi+1) = 1;
                    break
                else
                    tempTrainInput = trainingInputMat{planeInd}(:,setdiff(coeffInds,indPartial{pi}));
                    tempTestInput = testInputMat{planeInd}(finiteIndTest,setdiff(coeffInds,indPartial{pi}));
                    cvPartial = cvglmnet(tempTrainInput(finiteIndTrain,:), spkTrain, 'poisson', partialGlmOpt, [], lambdaCV, [], true);
                    
                    iLambda = find(cvPartial.lambda == cvPartial.lambda_1se);
                    partialLogLikelihood = sum(log(poisspdf(spkTest', exp([ones(length(finiteIndTest),1), tempTestInput] * [cvPartial.glmnet_fit.a0(iLambda); cvPartial.glmnet_fit.beta(:,iLambda)]))));
                    devianceFullPartial = 2*(fullLogLikelihood - partialLogLikelihood);
                    dfFullPartial = dfFullNull - cvPartial.glmnet_fit.df(iLambda);
                    if devianceFullPartial > chi2inv(1-pThresholdPartial, dfFullPartial)
                        fitResult(pi+1) = 1;
                    end
                end
            end
        end
    end

    fitResultsRe(cellnumInd,:) = fitResult;
    fitCoeffIndsRe(cellnumInd,:) = fitCoeffInd;
    doneRe(cellnumInd) = cellnum;
    cellTimeRe(cellnumInd) = toc(celltic);
end % end of parfor cellnum

%%
started(remainingCell) = startedRe;
fitLambda(remainingCell) = fitLambdaRe;
fitCoeffs(remainingCell) = fitCoeffsRe;
fitInd(remainingCell) = fitIndRe;
fitDF(remainingCell) = fitDFRe;
fitDevExplained(remainingCell) = fitDevExplainedRe;
fitCvDev(remainingCell) = fitCvDevRe;
fitResults(remainingCell,:) = fitResultsRe;
fitCoeffInds(remainingCell,:) = fitCoeffIndsRe;
done(remainingCell) = doneRe;
cellTime(remainingCell) = cellTimeRe;
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

save(savefnResultRe, 'fit*', 'allPredictors', '*InputMat', 'indPartial', '*Group', '*Tn', 'lambdaCV', '*Opt', 'done', '*Re', 'remainingCell', 'pThreshold*', '*Shift', 'testInd', 'trainingInd', 'cIDAll', 'cellTime');
%%
%         end % of ri. random group selection index
push_myphone(sprintf('GLM done for JK%03d S%02d', mouse, session))


%%
% myCluster = parcluster('local');
% delete(myCluster.Jobs)
% clear myCluster
% parpool(34, 'SpmdEnabled', true);