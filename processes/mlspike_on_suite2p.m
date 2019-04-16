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

mice = [52];
sessions = {[26]};  
baseDir = 'Y:\Whiskernas\JK\suite2p\';
for mi = 1 : length(mice)
    mouse = mice(mi);
    cd(sprintf('%s%03d',baseDir, mice(mi)))
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        flist = dir(sprintf('F_%03d_%03d_plane*_proc_final.mat', mouse, session));
%         for fi = 1 : length(flist)
        for fi = 1
            fprintf('mouse %03d session %d plane %d\n', mouse, session, fi)
            load(flist(fi).name)
            dt = 1/dat.ops.imageRate*dat.ops.num_plane;
            savefn = flist(fi).name;
            %%
            n = zeros(size(dat.dF), 'single');
            fit = zeros(size(dat.dF), 'single');
            drift = zeros(size(dat.dF), 'single');

            parfor ci = 1 : size(dat.dF,1)
                par = tps_mlspikes('par');
                par.a = 0.3;
                par.tau = 0.6;
                par.pnonlin = [0.73 -0.05];
                par.drift.parameter = 0.01;
                par.dt = dt;
                par.finetune.sigma = dat.noise(ci);
                par.algo.nspikemax = ceil(50 * dt); % 95th percentile ~ 50 Hz (ref. O'Connor 2010b Neuron)
                [n(ci,:), fit(ci,:), ~, ~, ~, drift(ci,:)] = tps_mlspikes(dat.dF(ci,:), par);
            end

            spk.n = n; spk.fit = fit; spk.drift = drift;
            save(savefn, 'dat', 'spk')
        end
    end
end