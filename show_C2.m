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
num_plane = length(info.otwave);
if nargin < 2
    plane = 1:num_plane;
    baseline_time = 5;
elseif nargin == 2
    plane = varargin{1};
    baseline_time = 5;
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

if info.scanmode
    freq = 15;
else
    freq = 30;
end
baseline_frame_num_each = floor(baseline_time * freq);

stim_start_frames = info.frame(info.event_id == 3);
baseline_start_frames = stim_start_frames - baseline_frame_num_each;
stim_end_frames = info.frame(info.event_id == 2);
if stim_end_frames(1) < stim_start_frames(1)
    stim_end_frames = stim_end_frames(2:end);
end
baseline_end_frames = stim_start_frames - 1;
num_trial = length(stim_start_frames);

figure,
for i_plane = 1 : length(plane)
    plane_frames = plane(i_plane)-1:num_plane:max_idx;

    baseline_frames = [];
    stim_frames = [];
    for i = 1 : length(stim_start_frames)
        baseline_frames = [baseline_frames, intersect(plane_frames,baseline_start_frames(i):baseline_end_frames(i))];
        stim_frames = [stim_frames, intersect(plane_frames,stim_start_frames(i):stim_end_frames(i))];
    end

    baseline_im = zeros(info.sz);
    stim_im = zeros(info.sz);
    for i = 1 : length(baseline_frames)
        % No need to use alignment information because it's under anesthesia 2017/10/19 JK
%         baseline_im = baseline_im + circshift(double(squeeze(sbxread(fn,baseline_frames(i),1))/length(baseline_frames)), T(floor(baseline_frames(i)/num_plane),:));
        baseline_im = baseline_im + double(squeeze(sbxread(fn,baseline_frames(i),1))/length(baseline_frames));
    end
    for i = 1 : length(stim_frames)
        stim_im = stim_im + double(squeeze(sbxread(fn,stim_frames(i),1))/length(stim_frames));
    end

    diff_im = (stim_im - baseline_im)./ baseline_im * 100;

    subplot(subplotnum1,subplotnum2,i_plane), imagesc(diff_im), axis image
end
