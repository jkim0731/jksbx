% clear

% par = tps_mlspikes('par');
% par.a = 0.3;
% par.tau = 0.6;
% par.pnonlin = [0.73 -0.05];
% par.drift.parameter = 0.01;

% calcium = dat.dF(2,:);
% par.finetune.sigma = dat.noise(2);
% par.dt = 0.1295;
% [n, fit, ~, ~, ~, drift] = tps_mlspikes(calcium,par);
% figure, plot(calcium), hold on, plot(fit,'k-'), plot(drift, 'k--'), yyaxis right, plot(n), ylim([0 20])

stdwindow = 5;
lowerprct = 5; % 5th percentile

mice = [74:76];
baseDir = 'D:\2p\JK\';
for mi = 1 : length(mice)
    cd(sprintf('%s%03d',baseDir, mice(mi)))
    flist = dir('F_*_proc_final.mat');
    for fi = 1 : length(flist)
%     for fi = 1
        load(flist(fi).name)
%         indCell = find([dat.stat.iscell]);
        dt = 1/dat.ops.imageRate*dat.ops.num_plane;
        savefn = flist(fi).name;
        %% test for smoothing dF 2018/11/25 JK
%         savefn = [flist(fi).name(1:end-4), '_smoothed.mat'];
%         savefn = flist(fi).name;
%         errorDiffs = cell(size(dat.dF,1),1);
%         errorResponse = zeros(size(dat.dF,1),2);
%         for j = 1 : size(dat.dF,1)
%             tempstd = zeros(size(dat.dF,2) - (stdwindow-1), 1);
%             temp = smooth(dat.dF(j,:));
%             temp = dat.dF(j,:);
%             tempsmooth = smooth(temp,5);
%             parfor k = 1 : length(tempstd)
%                 tempstd(k) = std(tempsmooth(k+2:k+(stdwindow-1)/2+2)); % chopping off flanking half-windows
%             end
%             tempstdInds = find(tempstd<=prctile(tempstd,lowerprct));
%             baseInds = zeros(length(tempstdInds),stdwindow);
%             diffs = zeros(length(tempstdInds),stdwindow);
%             parfor k = 1 : length(tempstdInds)
%                 baseInds(k,:) = tempstdInds(k) : tempstdInds(k)+(stdwindow-1);
%                 diffs(k,:) = temp(baseInds(k,:)) - mean(temp(baseInds(k,:)));
%             end            
%             errorDiffs{j} = diffs(:);
%         end
%         dat.noise = cellfun(@(x) sqrt(sum(x.^2)/length(x)), errorDiffs);
        %%
        
        n = zeros(size(dat.dF), 'single');
        fit = zeros(size(dat.dF), 'single');
        drift = zeros(size(dat.dF), 'single');
        
        parfor ci = 1 : size(dat.dF,1)
            par = tps_mlspikes('par');
            par.a = 0.3;
            %%
%             par.a = 0.16; % for smoothed data
            %%
            par.tau = 0.6;
            par.pnonlin = [0.73 -0.05];
            par.drift.parameter = 0.01;
            par.dt = dt;
            par.finetune.sigma = dat.noise(ci);
            par.algo.nspikemax = ceil(50 * dt); % 95th percentile ~ 50 Hz (ref. O'Connor 2010b Neuron)
            [n(ci,:), fit(ci,:), ~, ~, ~, drift(ci,:)] = tps_mlspikes(dat.dF(ci,:), par);
        end

        spk.n = n; spk.fit = fit; spk.drift = drift;
        
%         for ci = 1 : size(dat.dF,1)
%             par(ci) = tps_mlspikes('par');
%             par(ci).a = 0.3;
%             par(ci).tau = 0.6;
%             par(ci).pnonlin = [0.73 -0.05];
%             par(ci).drift.parameter = 0.01;
%             par(ci).dt = dt;
% %             par(ci).algo.dogpu = 1;
%             par(ci).finetune.sigma = dat.noise(ci);
% %             if dat.stat(indCell(ci)).depth > 350 % L4 cells
%                 par(ci).algo.nspikemax = ceil(50 * dt); % 95th percentile ~ 50 Hz (ref. O'Connor 2010b Neuron)
% %             else % L2/3 cells
% %                 par(ci).algo.nspikemax = ceil(20 * dt); % 95th percentile ~ 20 Hz (ref. O'Connor 2010b Neuron)
% %             end
%         end
%         [spk.n, spk.fit, ~, ~, ~, spk.drift] = tps_mlspikes(num2cell(dat.dF, 2), par);

        save(savefn, 'dat', 'spk')
    end
end