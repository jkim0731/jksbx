function diffImg = find_C2(fn,varargin)

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
% if ~exist([fn,'.align'],'file')
%     jksbxaligndir_optotune_top_piezo({fn})
% end
% 
% load([fn,'.align'],'-mat')

a = sbxread(fn,0,1);
if info.volscan
    num_plane = length(info.otwave);
else
    num_plane = 1;
end

if nargin < 2
    plane = num_plane;
    baseline_time = 2;
elseif nargin == 2
    plane = varargin{1};
    baseline_time = 2;
else
    plane = varargin{1};
    baseline_time = varargin{2};
end

max_idx = jkget_maxidx(fn);
plane_frames = plane-1:num_plane:max_idx;

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
diffImg = zeros([info.sz,length(plane)]);

for i_plane = 1 : length(plane)
    plane_frames = plane(i_plane)-1:num_plane:max_idx;    
    for i = 1 : length(stim_start_frames)
        baseline_frames = intersect(plane_frames,baseline_start_frames(i):baseline_end_frames(i));
        stim_frames = intersect(plane_frames,stim_start_frames(i):stim_end_frames(i));
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
        
        diffImg(:,:,i_plane) = diffImg(:,:,i_plane) + (stimImg - baseImg)./ baseImg * 100 / length(stim_start_frames);
    end
end
diffImg = squeeze(diffImg);