cellnum = 654
coeff = fitCoeffs{cellnum};
tn = testTn{cellnum};



cellID = u.cellNums(cellnum);
planeInd = floor(cellID/1000);

totalTind = find(cellfun(@(x) ismember(cellID, x.neuindSession), u.trials));

[~, testTind2] = ismember(tn, u.trialNums);
testTind3 = intersect(totalTind, testTind2);
obsInds = find(cell2mat(cellfun(@(x) ones(1,length(x.tpmTime{1})+2*posShift) * ismember(x.trialNum, tn), u.trials(totalTind)', 'uniformoutput', false)));
testInput2 = allPredictors{planeInd}(obsInds,:);

cind = find(u.trials{totalTind(1)}.neuindSession == cellID);


spkTest2 = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cind,:), nan(1,posShift)], u.trials(testTind3)','uniformoutput',false));
size(spkTest2)
size(testInput2)
%
spkTest2 = spkTest2';
finiteIndTest2 = intersect(find(isfinite(spkTest2)), find(isfinite(sum(testInput2,2))));
spkTest2 = spkTest2(finiteIndTest2)';

fitResult = zeros(1,6);

model = exp([ones(length(finiteIndTest2),1),testInput2(finiteIndTest2,:)]*coeff);
mu = mean(spkTest2); % null poisson parameter
nullLogLikelihood = sum(log(poisspdf(spkTest2,mu)));
fullLogLikelihood = sum(log(poisspdf(spkTest2',model)));
saturatedLogLikelihood = sum(log(poisspdf(spkTest2,spkTest2)));
devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);
devExp = (fullLogLikelihood- nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);

devExp 
fitDevExplained(cellnum)
