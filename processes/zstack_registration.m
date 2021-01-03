% Making registered z-stack 
% from z-stack files
% Translational registration based on intensity.
% Align each frame (2:end) to adjacent lower frame (1:end-1) first,
% and then apply consecutive imwarp through all layers.
% (This works better than first registering each and then calculating
% transformation. Maybe becuase of continuous x-y translation across depth
% and imMarin not catching it well. And also this is better in dealing 
% with errors.) 2020/05/12 JK

% Considered methods:
% (1) sbx fft align. Very fast, but only considers center part, only in 
% rectangular form. When I need to remove some margin, this is restricting
% too much. Not useful.
% (2) MATLAB imregtform. Blurring when using imwarp iteratively.
% (2)-1 modifying tform iteratively (how?), and then applying imwarp.

% Run each z-stack to check the quality of registration, and apply any
% modification necessary.

% Input:
%     z-stack file (containing 'zstack' and 'knobbyInfo')
%     
% Output:
%     Store 'zstackDepths', 'zstackReg', 'tform', and 'tformSerial' 
%     in a separate file (zstackReg_ file). 
%     zstackReg in gray scale. CLAHE is applied to zstackReg.

%% Basic setting

imMargin = 50; % variable setting

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,39,52];
numMice = length(mice);
zstackFilenum = [1000, 1000, 2000, 997, 998, 995];


%% load data for each mouse
mi = 3;
mouse = mice(mi);

% load z-stack
fn = sprintf('%s%03d\\zstack_%03d_%d',baseDir,mouse,mouse,zstackFilenum(mi));
load(fn, 'zstack', 'knobbyInfo')

savefn = sprintf('%s%03d\\zstackReg_%03d_%d',baseDir,mouse,mouse,zstackFilenum(mi));

% zstack depths from knobby info
zstackDepths = zeros(length(knobbyInfo), 1);
for di = 1 : length(knobbyInfo)
    zstackDepths(di) = - knobbyInfo(di).z / cosd(knobbyInfo(di).a); % conver to positive for depth.
end

%% Running frame-by-frame registration
% Only getting tform first
if mi == 6 % onl for JK052
    zstack = zstack(:,101:725,:);
end
zstack = mat2gray(zstack);
tformSerial = cell(size(zstack,3),1);
tformSerial{end} = affine2d(eye(3));
[opt, met] = imregconfig('monomodal');
tic
fprintf('Registerting JK%03d\n', mouse)
parfor i = 1 : size(zstack,3)-1
    fixed = adapthisteq(squeeze(zstack(:,:,i+1)));
    moving = adapthisteq(squeeze(zstack(:,:,i)));
    t = imregtform(moving(imMargin:end-imMargin,imMargin:end-imMargin), fixed(imMargin:end-imMargin,imMargin:end-imMargin), 'translation', opt, met);
    tformSerial{i} = t;
end
toc


%% tform QC

diffT = zeros(size(tformSerial,1),1);
for i = 1 : size(tformSerial,1)
    diffT(i) = sqrt(sum(sum((tformSerial{i}.T - eye(3)).^2)));
end

figure,
plot(diffT)

%%
ti = 82;
figure
subplot(121), imshowpair(adapthisteq(zstack(:,:,ti+1)), adapthisteq(zstack(:,:,ti)))
title('raw zstack')
ztemp = imwarp(zstack(:,:,ti), tformSerial{ti}, 'OutputView', imref2d([size(zstack,1),size(zstack,2)]));
subplot(122), imshowpair(adapthisteq(zstack(:,:,ti+1)), adapthisteq(ztemp))
title('registered')

%% When there's no error, run registration
zstackReg = zeros(size(zstack));
zstackReg(:,:,1) = adapthisteq(zstack(:,:,1));
tform = cell(size(zstack,3),1);
tform{end} = affine2d(eye(3));
tic
fprintf('Making registered zstack of JK%03d\n', mouse)
parfor i = 1 : size(zstack,3)-1
    ztemp = zstack(:,:,i);
    tempT = eye(3);
    for j = i:size(zstack,3)-1
        tempT = tempT * tformSerial{j}.T;
    end
    tform{i} = affine2d(tempT);
    ztemp = imwarp(ztemp, tform{i}, 'OutputView', imref2d(size(ztemp)));
    zstackReg(:,:,i) = adapthisteq(ztemp);
end
toc

%% registration QC
ti = 80;
tdepth = 50; % upward
figure
subplot(121), imshow(adapthisteq(mean(zstack(:,:,ti:ti+tdepth),3)))
subplot(122), imshow(mean(zstackReg(:,:,ti:ti+tdepth),3))

%%
ti = 149;
pdiff = 20; % 20 planes difference, upward
figure
subplot(121), imshowpair(adapthisteq(zstack(:,:,ti)), adapthisteq(zstack(:,:,ti+pdiff)))
title('raw zstack')
% ztemp = zstack(:,:,ti+pdiff);
% for i = ti+pdiff-1:-1:ti
%     ztemp = imwarp(ztemp, tform{i}, 'OutputView', imref2d([size(zstack,1),size(zstack,2)]));
% end
% subplot(122), imshowpair(adapthisteq(zstack(:,:,ti)), adapthisteq(ztemp))
subplot(122), imshowpair((zstackReg(:,:,ti)), (zstackReg(:,:,ti+pdiff)))
title('registered')

%% Save results
tic
fprintf('Saving JK%03d\n', mouse)
save(savefn,'zstackDepths','zstackReg','tform','tformSerial');
toc