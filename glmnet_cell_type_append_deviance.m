baseDir = 'C:\JK\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 

repetition = 10;

for mi = 1 : length(mice)
    mouse = mice(mi);
    cd(sprintf('%s%03d', baseDir, mouse))
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        fntemplate = sprintf('glmResponseType_JK%03dS%02d_glmnet_m18_R',mouse,session);
        for j = 1 : repetition
            if ~exist(sprintf('%s%02d.mat',fntemplate,j),'file')
                break
            end
        end
        if j < repetition
            break
        end
        
        load(sprintf('UberJK%03dS%02d',mouse,session))

        for ri = 1 : repetition
            fn = sprintf('%s%02d.mat',fntemplate,ri);
            dat = load(sprintf('%s%02d.mat',fntemplate,ri));
            if ~isfield(dat, 'fitDeviance')
                numCell = length(dat.cIDAll);
                
                fitDeviance = zeros(length(numCell),1);
                
                tindcellAll = cell(numCell,1);
                cindAll = zeros(numCell,1);
                planeIndAll = zeros(numCell,1);
                iTrainAll = cell(numCell,1);
                iTestAll = cell(numCell,1);
                for i = 1 : numCell
                    tindcellAll{i} = find(cellfun(@(x) ismember(dat.cIDAll(i), x.neuindSession), u.trials));
                    cindAll(i) = find(u.trials{tindcellAll{i}(1)}.neuindSession == dat.cIDAll(i));
                    planeIndAll(i) = floor(dat.cIDAll(i)/1000);
                    iTrainAll{i} = intersect(tindcellAll{i}, dat.trainingInd);
                    iTestAll{i} = intersect(tindcellAll{i}, dat.testInd);
                end
                spikeAll = cellfun(@(x) x.spk, u.trials, 'uniformoutput', false);
                testInputMat = dat.testInputMat;
                posShift = dat.posShift;
                fitCoeffs = dat.fitCoeffs;
                parfor cellnum = 1 : numCell
                    iTrain = iTrainAll{cellnum};
                    cind = cindAll(cellnum);
                    planeInd = planeIndAll(cellnum);

                    iTest = iTestAll{cellnum};                
                    spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x(cind,:), nan(1,posShift)], spikeAll(iTest)','uniformoutput',false));
                    spkTest = spkTest';
                    finiteIndTest = intersect(find(isfinite(spkTest)), find(isfinite(sum(testInputMat{planeInd},2))));
                    spkTest = spkTest(finiteIndTest)';
                    
                    model = exp([ones(length(finiteIndTest),1),testInputMat{planeInd}(finiteIndTest,:)]*[fitCoeffs{cellnum}]);
                    mu = mean(spkTest); % null poisson parameter
                    nullLogLikelihood = sum(log(poisspdf(spkTest,mu)));
                    fullLogLikelihood = sum(log(poisspdf(spkTest',model)));
                    saturatedLogLikelihood = sum(log(poisspdf(spkTest,spkTest)));
                    fitDeviance(cellnum) = 2*(fullLogLikelihood - nullLogLikelihood);                    
                end
                save(fn, 'fitDeviance', '-append')
            end
            dat = load(sprintf('%s%02d.mat',fntemplate,ri));
            if ~isfield(dat, 'fitCorrelation') || ~isfield(dat, 'fitCorrPval')
                numCell = length(dat.cIDAll);
                
                fitCorrelation = zeros(length(numCell),1);
                fitCorrPval = zeros(length(numCell),1);
                
                tindcellAll = cell(numCell,1);
                cindAll = zeros(numCell,1);
                planeIndAll = zeros(numCell,1);
                iTrainAll = cell(numCell,1);
                iTestAll = cell(numCell,1);
                for i = 1 : numCell
                    tindcellAll{i} = find(cellfun(@(x) ismember(dat.cIDAll(i), x.neuindSession), u.trials));
                    cindAll(i) = find(u.trials{tindcellAll{i}(1)}.neuindSession == dat.cIDAll(i));
                    planeIndAll(i) = floor(dat.cIDAll(i)/1000);
                    iTrainAll{i} = intersect(tindcellAll{i}, dat.trainingInd);
                    iTestAll{i} = intersect(tindcellAll{i}, dat.testInd);
                end
                spikeAll = cellfun(@(x) x.spk, u.trials, 'uniformoutput', false);
                testInputMat = dat.testInputMat;
                posShift = dat.posShift;
                fitCoeffs = dat.fitCoeffs;
                parfor cellnum = 1 : numCell
                    iTrain = iTrainAll{cellnum};
                    cind = cindAll(cellnum);
                    planeInd = planeIndAll(cellnum);

                    iTest = iTestAll{cellnum};                
                    spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x(cind,:), nan(1,posShift)], spikeAll(iTest)','uniformoutput',false));
                    spkTest = spkTest';
                    finiteIndTest = intersect(find(isfinite(spkTest)), find(isfinite(sum(testInputMat{planeInd},2))));
                    spkTest = spkTest(finiteIndTest)';                    
                    model = exp([ones(length(finiteIndTest),1),testInputMat{planeInd}(finiteIndTest,:)]*[fitCoeffs{cellnum}]);
                    [fitCorrelation(cellnum), fitCorrPval(cellnum)] = corr(spkTest', model);
                end
                save(fn, 'fitCorrelation', 'fitCorrPval', '-append')
            end
            fprintf('JK%03d S%02d rep%02d done.\n', mouse, session, ri)
        end        
    end
end
        