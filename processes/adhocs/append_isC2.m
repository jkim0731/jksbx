%%
clear
% mice = [25,27,30,36,37,39,52,53,54,56];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  

% mice = [70, 74,75,76];
% sessions = {[4], [4],[4],[4]};  
mice = [52];
sessions = {[26]};  

baseD = 'Y:\Whiskernas\JK\suite2p\';
for mi = 1 : length(mice)
    mouse = mice(mi);
    cd(sprintf('%s%03d',baseD,mouse));
    load(sprintf('JK%03dC2.mat',mouse))
%     flist = dir('F_0*_final.mat');
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        flist = dir(sprintf('F_%03d_%03d_plane*_proc.mat', mouse, session));    
        for fi = 1 : length(flist)
            fprintf('JK%03d S%02d plane%d processing\n', mouse, session, fi)
            load(flist(fi).name)        
            inds = find([dat.stat.iscell]);
            dat.isC2 = zeros(length(inds),1);
            for i = 1 : length(inds)
                if inpolygon(dat.ops.useY(dat.ops.yrange(round(dat.stat(inds(i)).med(1)))), dat.ops.useX(dat.ops.xrange(round(dat.stat(inds(i)).med(2)))), ypoints, xpoints)
                    dat.isC2(i) = 1;
                end
            end
            dat.c2ypoints = ypoints;
            dat.c2xpoints = xpoints;
    %         save(flist(fi).name,'dat')
            savefn = [flist(fi).name(1:end-4), '_final.mat'];
            save(savefn,'dat')
        end
    end
end
