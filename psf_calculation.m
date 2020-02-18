% calculating psf

%% load z-stack file
% it should have knobby info and processed z-stack and max projection
load('zstack_psf_191124_3_.mat')

%% pick points of the bead

figure, imagesc(zstack_maxproj)
%%
[x,y] = ginput(1);

regionSize = 20;
yrange = round(y-regionSize:y+regionSize);
xrange = round(x-regionSize:x+regionSize);

temp3D = zstack(yrange, xrange, :);

implay(mat2gray(temp3D))

%%
implay(zstack)

%%
zstackSel = zstack(:,:,30:45);
figure, imagesc(max(zstackSel,[],3))
%%
implay(zstackSel)
%%
figure, imagesc(max(zstackSel,[],3))
[x,y] = ginput(1);

regionSize = 20;
yrange = round(y-regionSize:y+regionSize);
xrange = round(x-regionSize:x+regionSize);

temp3D = zstack(yrange, xrange, :);

implay(mat2gray(temp3D))


%%
axialPsf = squeeze(max(temp3D));
figure, imagesc(squeeze(max(temp3D)))

%%
axialIntensity = max(axialPsf(20:22,:));
figure, plot(axialIntensity)


%%
minVal = mean(axialIntensity([1:20,60:107]));
normIntensity = axialIntensity-minVal;

%%
pixFWHM = length(find(normIntensity > max(normIntensity)/2));

%%
zstep = mean(diff([knobbyInfo.z]));

FWHM = zstep * pixFWHM