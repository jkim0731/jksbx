baseDir = 'D:\';
mouse = 37;
session = 5554;
sbxList = ls(sprintf('%s%03d\\%03d_%d_*.sbx',baseDir,mouse,mouse,session));
for sbxi = 1:length(sbxList)
    sbxFn = sbxList(sbxi,1:end-4);
    figure, plot(jksbx_mean_signal(sbxFn))
    drawnow;
end