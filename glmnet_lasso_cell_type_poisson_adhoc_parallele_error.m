 save(savefnResult, 'fit*', 'allPredictors', '*InputMat', 'indPartial', '*Group', '*Tn', 'lambdaCV', '*Opt', 'done');
 
 %%
mouse = 36;
session = 1;
savefnResult = sprintf('glmResponseType_JK%03dS%02d_glmnet_m14',mouse, session);

savefnResultRe = [savefnResult, '_02'];

previousDone = done(find(done));

remainingCell = setdiff(1 : length(u.cellNums), previousDone);


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
for cellnumInd = 2 : length(remainingCell)
    cellnum = remainingCell(cellnumInd);
%             for cellnum = 102, 127, (212 convergence error), 221, 658
%         ci = 0;
%         for cellnum = 1:division:length(u.cellNums)
%             ci = ci + 1;
%         for cellnum = 1
%     cellnum = 1;
    fitCoeffInd = zeros(1,6);

%                 fprintf('Mouse JK%03d session S%02d Loop %d: Running cell %d/%d \n', mouse, session, ri,cellnum, length(u.cellNums));
    fprintf('Mouse JK%03d session S%02d: Running cell %d/%d \n', mouse, session,cellnum, length(u.cellNums));
    startedRe(cellnumInd) = cellnum;

    cID = u.cellNums(cellnum);

    % find out trial indices for this specific cell
    tindcell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));

    iTrain = intersect(tindcell, trainingInd);
    % find out row number of this cell
    cind = find(u.trials{iTrain(1)}.neuindSession == cID);
    planeInd = floor(cID/1000);

    spkTrain = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cind,:), nan(1,posShift)], u.trials(iTrain)','uniformoutput',false));                
    %%
    finiteIndTrain = intersect(find(isfinite(spkTrain)), find(isfinite(sum(trainingInputMat{planeInd},2))));
    input = trainingInputMat{planeInd}(finiteIndTrain,:);
    spk = spkTrain(finiteIndTrain)';

    cv = cvglmnet(input, spk, 'poisson', glmnetOpt, [], lambdaCV);
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
    cID = u.cellNums(cellnum);

    % find out trial indices for this specific cell
    tindcell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));

    iTest = intersect(tindcell, testInd);
    % find out row number of this cell
    cind = find(u.trials{iTest(1)}.neuindSession == cID);
    planeInd = floor(cID/1000);

    spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cind,:), nan(1,posShift)], u.trials(iTest)','uniformoutput',false));
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
                    cvPartial = cvglmnet(tempTrainInput(finiteIndTrain,:), spk, 'poisson', partialGlmOpt, [], lambdaCV);
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

save(savefnResultRe, 'fit*', 'allPredictors', '*InputMat', 'indPartial', '*Group', '*Tn', 'lambdaCV', '*Opt', 'done', '*Re', 'remainingCell');
%%
%         end % of ri. random group selection index
push_myphone(sprintf('GLM done for JK%03d S%02d', mouse, session))


%%
% myCluster = parcluster('local');
% delete(myCluster.Jobs)
% clear myCluster
% parpool(34, 'SpmdEnabled', true);