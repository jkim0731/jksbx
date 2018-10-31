%%
clear
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  

baseD = 'C:\Data\suite2p\';
for mi = 1 : length(mice)
    cd(sprintf('%s%03d',baseD,mice(mi)));
    load(sprintf('JK%03dC2.mat',mice(mi)))
    flist = dir('F_0*_final.mat');
    for fi = 1 : length(flist)        
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
        save(flist(fi).name,'dat')
    end
end
