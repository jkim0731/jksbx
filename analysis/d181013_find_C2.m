i = 1;
inds = find([dat.stat.iscell]);
ind = inds(i);

signal = (dat.Fcell{1}(ind,:) - dat.FcellNeu{1}(ind,:) * dat.stat(ind).neuropilCoefficient);
%%
sg = sgolayfilt(double(signal),3,9);
%
figure, plot(signal), hold on,
plot(sg)
%%
plot(signal)
%%
fn = '025_995_002';
show_C2_old(fn,1:4,5)

%%
a = sbxread(fn,732,1);
figure, imshow(mat2gray(squeeze(a)))

%%
a = sbxread(fn,0,jkget_maxidx(fn));
%%
figure, plot(1:4:size(a,4),squeeze(mean(mean(a(1,:,:,1:4:end)))))