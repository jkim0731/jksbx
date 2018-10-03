function show_C2(fn,varargin)

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

max_idx = jkget_maxidx(fn);

freq = info.resfreq/info.sz(1)*(2-info.scanmode);

baseline_frame_num_each = floor(baseline_time * freq);
udpDelayBuffer = 0.3; % in s
bufferFrames = ceil(freq * udpDelayBuffer);
baseline_start_frames = info.frame(info.event_id == 3)+bufferFrames;
stim_start_frames = baseline_start_frames + baseline_frame_num_each-bufferFrames;
stim_end_frames = info.frame(info.event_id == 2) - floor(freq*2);
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

figure,
for i_plane = 1 : length(plane)
    plane_frames = plane(i_plane)-1:num_plane:max_idx;
    diffImg = zeros(info.sz);
    for i = 1 : min(length(stim_start_frames), length(stim_end_frames))
        baseline_frames = intersect(plane_frames,baseline_start_frames(i):baseline_end_frames(i));
        stim_frames = intersect(plane_frames,stim_start_frames(i):stim_end_frames(i));
        if ~isempty(baseline_frames) && ~isempty(stim_frames)
            tempImg = double(sbxread(fn,baseline_frames(1),1));
            baseImg = squeeze(tempImg(1,:,:,:))/length(baseline_frames);
            for j = 2 : length(baseline_frames)
                tempImg = double(sbxread(fn,baseline_frames(j),1));
                baseImg = baseImg + squeeze(tempImg(1,:,:,:)/length(baseline_frames));            
            end

            tempImg = double(sbxread(fn,stim_frames(1),1));
            stimImg = squeeze(tempImg(1,:,:,:))/length(stim_frames);
            for j = 2 : length(stim_frames)
                tempImg = double(sbxread(fn,stim_frames(j),1));
                stimImg = stimImg + squeeze(tempImg(1,:,:,:)/length(stim_frames));            
            end
            diffTemp = (stimImg - baseImg)./ baseImg * 100;
            diffImg = diffImg + diffTemp / length(stim_start_frames);
            
        end
    end

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

    subplot(subplotnum1,subplotnum2,i_plane), imagesc(diffImg(101:end,101:end-70),[-50 100]), axis image, axis off
end
