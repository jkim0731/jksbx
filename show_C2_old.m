function diffImg = show_C2_old(fn,varargin)

% For JK025,027, and 030: baseline = 5
% For JK036~041 : baseline = 4

global info
fn = strtok(fn,'.');
if ~exist([fn,'.sbx'],'file')
    try
        cd(['D:\2p\JK\', fn(1:3)])
    catch
        error('Could not find the directory')
    end
    if ~exist([fn,'.sbx'],'file')
        error('Could not find the file')
    end
end
% No need to align because it's under anesthesia 2017/10/19 JK
% if ~exist([fn,'.align'],'file')
%     jksbxaligndir_optotune_top_piezo({fn})
% end
% 
% load([fn,'.align'],'-mat')
a = sbxread(fn,0,1);
if info.volscan == 0
    num_plane = 1;
else
    num_plane = length(info.otwave);
end
if nargin < 2
    plane = 1:num_plane;
    baseline_time = 3;
elseif nargin == 2
    plane = varargin{1};
    baseline_time = 3;
else
    plane = varargin{1};
    baseline_time = varargin{2};
end

if length(plane) < 2
    subplotnum1 = 1;
    subplotnum2 = 1;
elseif length(plane) < 3
    subplotnum1 = 1;
    subplotnum2 = 2;
elseif length(plane) < 10
    subplotnum1 = 2;
    subplotnum2 = ceil(length(plane)/2);
else
    subplotnum1 = 3;
    subplotnum2 = ceil(length(plane)/3);
end

max_idx = info.max_idx;

freq = info.resfreq/info.sz(1)*(2-info.scanmode);

baseline_frame_num_each = floor(baseline_time * freq);
stim_start_frames = info.frame(info.event_id == 3);
baseline_start_frames = stim_start_frames - baseline_frame_num_each;
stim_end_frames = info.frame(info.event_id == 2);
% stim_end_frames = stim_start_frames + freq;
if length(stim_end_frames) < length(stim_start_frames)
    stim_start_frames = stim_start_frames(1:end-1);
end
removeFrames = [];
for i = 1 : length(stim_end_frames)
    if stim_end_frames(i) < stim_start_frames(1)
        removeFrames = [removeFrames, i];
    else
        break
    end
end
stim_end_frames(removeFrames) = [];
% stim_end_frames = stim_start_frames + 60; % for passive pole.... temporary
% if stim_end_frames(1) < stim_start_frames(1)
%     stim_end_frames = stim_end_frames(2:end);
% end
baseline_end_frames = stim_start_frames - 1;

diffImg = cell(length(plane),1);
% avgImg = zeros(info.sz);
% figure,
parfor i_plane = 1 : length(plane)
    plane_frames = plane(i_plane)-1:num_plane:max_idx;
    diffImg{i_plane} = zeros(info.sz);
    for i = 1 : min(length(stim_start_frames), length(stim_end_frames))
        baseline_frames = intersect(plane_frames,baseline_start_frames(i):baseline_end_frames(i));
        stim_frames = intersect(plane_frames,stim_start_frames(i):stim_end_frames(i));
        if ~isempty(baseline_frames) && ~isempty(stim_frames)
            baseImg = mean(squeeze(jksbxreadframes(fn,baseline_frames)),3);            
            stimImg = mean(squeeze(jksbxreadframes(fn,stim_frames)),3);
            
            baseImg = imgaussfilt(baseImg);
            stimImg = imgaussfilt(stimImg);
            
            diffTemp = (stimImg - baseImg)./ baseImg;
            diffImg{i_plane} = diffImg{i_plane} + diffTemp / length(stim_start_frames);
        end
    end

%     temp = (diffImg{i_plane}-min(min(diffImg{i_plane})))/(max(max(diffImg{i_plane})) - min(min(diffImg{i_plane}))); % normalization
%     subplot(subplotnum1,subplotnum2,i_plane), imagesc(diffImg{i_plane}(101:end,101:end-70)), axis image, axis off
    
    
    

%     avgImg = avgImg + temp;
%     baseline_im = zeros(info.sz);
%     stim_im = zeros(info.sz);
%     for i = 1 : length(baseline_frames)
%         % No need to use alignment information because it's under anesthesia 2017/10/19 JK
% %         baseline_im = baseline_im + circshift(double(squeeze(sbxread(fn,baseline_frames(i),1))/length(baseline_frames)), T(floor(baseline_frames(i)/num_plane),:));
%         temp_im = sbxread(fn,baseline_frames(i),1);
%         baseline_im = baseline_im + double(squeeze(temp_im(1,:,:,:))/length(baseline_frames));
%     end
%     for i = 1 : length(stim_frames)
%         temp_im = sbxread(fn,stim_frames(i),1);
%         stim_im = stim_im + double(squeeze(temp_im(1,:,:,:))/length(stim_frames));
%     end
% 
%     diff_im = (stim_im - baseline_im)./ baseline_im * 100;
end

% figure, imshow(mat2gray(avgImg))
