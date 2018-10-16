clear
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  
rollingWindowForBaseF = [100, 180, 300]; % in s % Tried varying timepoints and found out that baseline fluorescence value does not change
baseFprctile = 1;
baseDir = 'D:\TPM\JK\suite2p\';
for mi = 1 : length(mice)
% for mi = 1
    cd([baseDir, sprintf('%03d',mice(mi))])
    fList = dir('F_*_proc_final.mat');
    result = struct;

    for i = 1 : length(fList)
%     for i = 2
        fprintf('mouse %03d session %d\n', mice(mi), i)
        load(fList(i).name)
        
        for wi = 1 : length(rollingWindowForBaseF)
            inds = find([dat.stat.iscell]);
            len = size(dat.Fcell{1},2);
            window = round(rollingWindowForBaseF(wi)*(dat.ops.imageRate/dat.ops.num_plane));

            tempF = dat.Fcell{1}(inds,:) - dat.FcellNeu{1}(inds,:) .* [dat.stat(inds).neuropilCoefficient]';
            tempFF = [tempF, tempF(:,end-window+1:end)];
            ttbase = zeros(length(inds),size(tempF,2));
            parfor k = 1 : len
                ttbase(:,k) = prctile(tempFF(:,k:k+window),baseFprctile,2);
            end
            result(i).baseF{wi} = ttbase;
            result(i).dF{wi}  = (tempF - ttbase)./ttbase;            
            result(i).F{wi} = dat.Fcell{1}(inds,:) - (0.5 + wi * 0.1) * dat.FcellNeu{1}(inds,:);
            dF{wi} = dFtemp(find([dat.stat.iscell]),:);
            baseF{wi} = baseFtemp(find([dat.stat.iscell]),:);
            
            result(i).diffdF{wi} = diff(result(i).dF{wi},1,2);
        end
    end
end

%%

clear
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  
rollingWindowForBaseF = [100, 180, 300]; % in s % Tried varying timepoints and found out that baseline fluorescence value does not change
baseFprctile = 1;
baseDir = 'D:\TPM\JK\suite2p\';
% for mi = 1 : length(mice)
for mi = 1
    cd([baseDir, sprintf('%03d',mice(mi))])
    fList = dir('F_*_proc_final.mat');
    result = struct;

    for i = 1 : length(fList)
%     for i = 2
        fprintf('mouse %03d session %d\n', mice(mi), i)
        load(fList(i).name)
        
        for wi = 1 : length(rollingWindowForBaseF)
            inds = find([dat.stat.iscell]);
            len = size(dat.Fcell{1},2);
            window = round(rollingWindowForBaseF(wi)*(dat.ops.imageRate/dat.ops.num_plane));

            tempF = dat.Fcell{1}(inds,:) - dat.FcellNeu{1}(inds,:) .* [dat.stat(inds).neuropilCoefficient]';
            tempFF = [tempF, tempF(:,end-window+1:end)];
            ttbase = zeros(length(inds),size(tempF,2));
            parfor k = 1 : len
                ttbase(:,k) = prctile(tempFF(:,k:k+window),baseFprctile,2);
            end
            result(i).baseF{wi} = ttbase;
            result(i).dF{wi}  = (tempF - ttbase)./ttbase;            
            result(i).F{wi} = dat.Fcell{1}(inds,:) - (0.5 + wi * 0.1) * dat.FcellNeu{1}(inds,:);
            dF{wi} = dFtemp(find([dat.stat.iscell]),:);
            baseF{wi} = baseFtemp(find([dat.stat.iscell]),:);
            
            result(i).diffdF{wi} = diff(result(i).dF{wi},1,2);
        end
    end
end

%% How to deal with negative F (after neuropil correction)

% first, watch the distribution
% if it's just few cells having all the negative values, just ignore these
% cells.
% with different coefficient

ci = 3;
pi = 4;

[ii, jj] = ind2sub(size(result(pi).F{ci}),find(result(pi).F{ci}<0));
[aa, bb, cc] = (unique(ii));
negdist = zeros(length(aa),1);
for i = 1 : length(aa)
    negdist(i) = length(find(ii == aa(i)));
end

%%

clear
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  
rollingWindowForBaseF = [100, 180, 300]; % in s % Tried varying timepoints and found out that baseline fluorescence value does not change
baseFprctile = 1;
baseDir = 'D:\TPM\JK\suite2p\';
for mi = 1 : length(mice)
% for mi = 1
    cd([baseDir, sprintf('%03d',mice(mi))])
    fList = dir('F_*_proc_final.mat');
    result = struct;

    for i = 1 : length(fList)
        fprintf('mouse %03d session %d\n', mice(mi), i)
        load(fList(i).name)
        savefn = [fList(i).name(1:end-9), 'npcoeff.mat'];
        inds = find([dat.stat.iscell]);
        npcoeffs = min(dat.Fcell{1}(inds,:) ./ dat.FcellNeu{1}(inds,:), [], 2);
        save(savefn, 'npcoeffs')        
    end
end

%%
clear
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  
rollingWindowForBaseF = [100, 180, 300]; % in s % Tried varying timepoints and found out that baseline fluorescence value does not change
baseFprctile = 1;
baseDir = 'D:\TPM\JK\suite2p\';
for mi = 1 : length(mice)
% for mi = 1
    cd([baseDir, sprintf('%03d',mice(mi))])
    fList = dir('F_*_proc_npcoeff.mat');
    result = struct;

    for i = 1 : length(fList)
        npc(mi,i) = load(fList(i).name);
    end
end
%%
temp = [];
for i = 1 : length(mice)
    for j = 1 : size(npc,2)
        temp = [temp; npc(i,j).npcoeffs];
    end
end
figure, histogram(temp)
%%
basediff = baseF{3} - baseF{1};
dfdiff = dF{3} - dF{1};
%%
figure, histogram(dF{1}(4,:), 'normalization', 'probability')
%%
prctile(diffdF{1}(4,:),90)
prctile(diffdF{1}(4,:),10)
%%
val = zeros(100000,1);
for i = 1 : length(val)
    temp = diffdF{1}(4,:);
    inds = randperm(length(temp),7);
    val(i) = mean(cumsum(temp(inds)));
end
%%
figure, histogram(val,'normalization','probability')
prctile(val,95)
prctile(val,5)
%%
pct = zeros(size(diffdF{1},1),2,3);
for i = 1 : 3
    pct(:,1,i) = prctile(diffdF{i},5,2);
    pct(:,2,i) = prctile(diffdF{i},95,2);
end
%%
figure, 
for ii = 1 : 2
    subplot(1,2,ii), hold on
    for i = 1 : 3
        histogram(pct(:,ii,i), [-2:0.05:2],'normalization', 'probability')
    end
end
% 
% %%
% figure, hold on
% for i = 1 : 3
%     plot(dF