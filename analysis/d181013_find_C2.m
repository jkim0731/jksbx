%% F, drawing and saving in suite2p folder
clear
mouse = '036';
sbxdir = 'D:\TPM\JK\';
suite2pdir = 'D:\TPM\JK\suite2p\';
cd([sbxdir, mouse])
%%
fn = '056_997_000';
diffim = show_C2_old(fn);
%
figure, imagesc(diffim{1}, [0 0.4]), axis image, axis off, hold on
%%
cd([sbxdir, mouse])
totalDiffim = cell(8,1);
for i = 1 : 8
    fn = sprintf('036_9998_10%d',i);
    totalDiffim{i} = show_C2_old(fn);
end

%%
temp = zeros(size(totalDiffim{1}{1}));
for i = 1 : 8
    for j = 3
        temp = temp + totalDiffim{i}{j}/8;
    end
end
figure, imagesc(temp, [0 0.4]), axis image, axis off, hold on


%%

xpoints = [];
ypoints = [];
while true
    [x, y] = ginput(1);
    if isempty(x) % enter
        break
    end
    xpoints = [xpoints, x];
    ypoints = [ypoints, y];
    plot(xpoints, ypoints, 'r-', 'linewidth', 5)
end
plot([xpoints, xpoints(1)], [ypoints, ypoints(1)], 'r-', 'linewidth', 5)

cd([suite2pdir,mouse])
savefn = sprintf('JK%sC2.mat',mouse);
save(savefn, 'xpoints', 'ypoints')

%%
mouse = '054';
fn = '054_997_000';
loadfn = 'zstack_054_995';
sbxdir = 'D:\TPM\JK\';
suite2pdir = 'D:\TPM\JK\suite2p\';

%%
cd([sbxdir, mouse])
load(fn)
numlayers = length(info.otwave);
temp = jksbxreadframes(fn,0:numlayers:jkget_maxidx(fn));
figure, plot(squeeze(mean(mean(temp))))



%%
cd([sbxdir, mouse])
totalDiffim = cell(8,1);
for i = 1 : 8
    fn = sprintf('054_6001_00%d',i);
    totalDiffim{i} = show_C2(fn);
end

%%
temp = zeros(size(totalDiffim{1}{1}));
for i = 1 : 8
    for j = 1       
        temp = temp + totalDiffim{i}{j}/8;        
    end
end

figure, imagesc(temp, [0 0.4]), axis image, axis off


%%
cd([sbxdir, mouse])
% fn = '053_999_008';
diffim = show_C2(fn);
%
figure, imagesc(diffim{1}, [0 0.4]), axis image, axis off

%%

% loadfn = 'zstack_041_998';
cd([suite2pdir,mouse])
load(loadfn)
% StackSlider(zstack)
%%
zstack = make_zstack('041_998',1);
%%
figure, imshow(mat2gray(mean(zstack(:,:,end-80:end-30),3)))


%%
cd(['Y:\Whiskernas\JK\ISI\JK', mouse])
%%
a = read_qcamraw('vas.qcamraw',1);
imtool(mat2gray(a'))

%%
isi_image('c2')

%%
isi_register('c2',0.045)


