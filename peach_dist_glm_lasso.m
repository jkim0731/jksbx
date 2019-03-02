function [fitInd, fitCoeffs, fitCoeffInd, fitResult, fitDevExplained, fitCvDev, fitLambda, dfFullNull] = ...
    peach_dist_glm_lasso(iTrainAll, cindAll, planeIndAll, posShift, spikeAll, trainingInputMat, glmnetOpt, lambdaCV, indPartial, ...
                iTestAll, testInputMat, pThresholdNull, partialGlmOpt, pThresholdPartial, cellnum)

    fitCoeffInd = nan(1,6);

    iTrain = iTrainAll{cellnum};
    cind = cindAll(cellnum);
    planeInd = planeIndAll(cellnum);

    spkTrain = cell2mat(cellfun(@(x) [nan(1,posShift), x(cind,:), nan(1,posShift)], spikeAll(iTrain)','uniformoutput',false));                
    finiteIndTrain = intersect(find(isfinite(spkTrain)), find(isfinite(sum(trainingInputMat{planeInd},2))));
    input = trainingInputMat{planeInd}(finiteIndTrain,:);
    spkTrain = spkTrain(finiteIndTrain)';

    cv = cvglmnet(input, spkTrain, 'poisson', glmnetOpt, [], lambdaCV);
    %% survived coefficients
    fitLambda = cv.lambda_1se;
    iLambda = find(cv.lambda == cv.lambda_1se);
    fitCoeffs = [cv.glmnet_fit.a0(iLambda);cv.glmnet_fit.beta(:,iLambda)];
    coeffInds = find(cv.glmnet_fit.beta(:,iLambda));
    fitInd = coeffInds;
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
    fitDevExplained = 1 - (saturatedLogLikelihood - fullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
    fitCvDev = cv.glmnet_fit.dev(iLambda);

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
                    cvPartial = cvglmnet(tempTrainInput(finiteIndTrain,:), spkTrain, 'poisson', partialGlmOpt, [], lambdaCV);
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

