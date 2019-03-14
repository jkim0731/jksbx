        startRepetition = 3;
        for ri = startRepetition : repetition % repetition index
            parfor cellnum = 1 : numCell
%             for cellnum = 102, 127, (212 convergence error), 221, 658
    %         ci = 0;
    %         for cellnum = 1:division:length(u.cellNums)
    %             ci = ci + 1;
    %         for cellnum = 1
            %     cellnum = 1;
            if ~ismember(numCell, errorCellSession)
                cellTimeStart = tic;
                fitCoeffInd = zeros(1,6);
                
                fprintf('Mouse JK%03d session S%02d Loop %02d Repeat %02d: Running cell %d/%d \n', mouse, session, loop, ri, cellnum, numCell);
%                 fprintf('Mouse JK%03d session S%02d: Running cell %d/%d \n', mouse, session,cellnum, numCell);
                started(cellnum) = cellnum;
                
                iTrain = iTrainAll{cellnum};
                cind = cindAll(cellnum);
                planeInd = planeIndAll(cellnum);
    
                spkTrain = cell2mat(cellfun(@(x) [nan(1,posShift), x(cind,:), nan(1,posShift)], spikeAll(iTrain)','uniformoutput',false));                
                finiteIndTrain = intersect(find(isfinite(spkTrain)), find(isfinite(sum(trainingInputMat{planeInd},2))));
                input = trainingInputMat{planeInd}(finiteIndTrain,:);
                spkTrain = spkTrain(finiteIndTrain)';
    
                cv = cvglmnet(input, spkTrain, 'poisson', glmnetOpt, [], lambdaCV);
                %% survived coefficients
                fitLambda(cellnum) = cv.lambda_1se;
                iLambda = find(cv.lambda == cv.lambda_1se);
                fitCoeffs{cellnum} = [cv.glmnet_fit.a0(iLambda);cv.glmnet_fit.beta(:,iLambda)];
                coeffInds = find(cv.glmnet_fit.beta(:,iLambda));                
    %             rtest(ri).fitInd{cellnum} = coeffInds;
                fitInd{cellnum} = coeffInds;
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
                fitDeviance(cellnum) = devianceFullNull;
                [fitCorrelation(cellnum), fitCorrPval(cellnum)] = corr(spkTest', model);            
                dfFullNull = length(coeffInds);                
                fitDF(cellnum) = dfFullNull;
                fitDevExplained(cellnum) = 1 - (saturatedLogLikelihood - fullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
                fitCvDev(cellnum) = cv.glmnet_fit.dev(iLambda);

                if devianceFullNull > chi2inv(1-pThresholdNull, dfFullNull)
                    fitResult(1) = 1;
                    
                    %% (2) test without each parameter (as a group)                
%                     for pi = 1 : 5
%                         if find(ismember(coeffInds, indPartial{pi}))
%                             if all(ismember(coeffInds, indPartial{pi}))
%                                 fitResult(pi+1) = 1;
%                                 break
%                             else
%                                 tempTrainInput = trainingInputMat{planeInd}(:,setdiff(coeffInds,indPartial{pi}));
%                                 tempTestInput = testInputMat{planeInd}(finiteIndTest,setdiff(coeffInds,indPartial{pi}));
%                                 cvPartial = cvglmnet(tempTrainInput(finiteIndTrain,:), spkTrain, 'poisson', partialGlmOpt, [], lambdaCV);
%                                 iLambda = find(cvPartial.lambda == cvPartial.lambda_1se);
%                                 partialLogLikelihood = sum(log(poisspdf(spkTest', exp([ones(length(finiteIndTest),1), tempTestInput] * [cvPartial.glmnet_fit.a0(iLambda); cvPartial.glmnet_fit.beta(:,iLambda)]))));
%                                 devianceFullPartial = 2*(fullLogLikelihood - partialLogLikelihood);
%                                 dfFullPartial = dfFullNull - length(setdiff(coeffInds, indPartial{pi}));
%                                 if devianceFullPartial > chi2inv(1-pThresholdPartial, dfFullPartial)
%                                     fitResult(pi+1) = 1;
%                                 end
%                             end
%                         end
%                     end
                end
                
                fitResults(cellnum,:) = fitResult;
                fitCoeffInds(cellnum,:) = fitCoeffInd;
                done(cellnum) = cellnum;
                cellTime(cellnum) = toc(cellTimeStart);
            end
            end % end of parfor cellnum
%%
            save(sprintf('%s_R%02d',savefnResult, ri), 'fit*', 'allPredictors', '*InputMat', 'indPartial', '*Group', '*Tn', 'lambdaCV', '*Opt', 'done', 'pThreshold*', '*Shift', 'cellTime', 'testInd', 'trainingInd', 'cIDAll');
            push_myphone(sprintf('Lasso GLM done for JK%03d S%02d Loop %02d repeat %02d', mouse, session, loop, ri))
        end % of ri. random group selection index
        push_myphone(sprintf('Lasso GLM done for JK%03d S%02d, Big loop %02d', mouse, session, loop))