function jksbxaligndir(varargin)
clear info
global info

if(nargin>=1) % cell with filenames to be aligned
    for(i=1:length(varargin{1}))
        d(i).name = varargin{1}{i};
    end
else
    d = dir('*.sbx');
end

% Align all *.sbx files in the list

for(i=1:length(d))
    fn = strtok(d(i).name,'.');
    if exist([fn,'.align'],'file')
        sprintf('File %s is already aligned',fn);
    elseif exist([fn,'.trials'],'file') % ignore frames outside of each trial, because those are blank. If not treated, these lead to weird interference pattern during bidirectional scanning.
        load([fn,'.trials'],'-mat')
        trial_frames = [];
        for tn = 1 : length(trials)
            trial_frames = [trial_frames, trials(tn).frames(1)+5:trials(tn).frames(2)-5]; % ignore first 5 and last 5 frames of each trial
        end

        sbxread(fn,0,1);            % read one frame to read the header of the image sequence
             % this contains the information about the structure of the image
        tic
        if info.volscan
            num_plane = length(info.otwave_um);
            for j = 1 : num_plane
                temp_start = num_plane + j -1; % to discard first frames of each plane. First frame is 0th. 
                frame_dif = mod(info.max_idx+1,num_plane) - mod(j,num_plane);
                if frame_dif < 0
                    temp_end = info.max_idx - (num_plane + frame_dif);
                else
                    temp_end = info.max_idx - frame_dif;
                end
                [frame_to_align, fta_ind] = intersect(temp_start:num_plane:temp_end,trial_frames);
                [m{j}, T_temp] = jksbxalignx(fn,frame_to_align);
                T{j} = zeros(length(temp_start:num_plane:temp_end),2);
                for Tind = 1 : length(fta_ind)
                    T{j}(fta_ind(Tind),:) = T_temp(Tind,:);
                end
                T{j} = [0 0;T{j}]; % To compensate for ignoring the first frame
            end
        else
            [frame_to_align, fta_ind] = intersect(0:info.max_idx,trial_frames);
            [m,T_temp] = jksbxalignx(fn,frame_to_align);   %
            T = zeros(info.max_idx+1,2);
            for Tind = 1 : length(fta_ind)
                T(fta_ind(Tind),:) = T_temp(Tind,:);
            end
        end
        save([fn '.align'],'m','T','frame_to_align','fta_ind');
        clear m T
        display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx,round(toc/60)));
    else
        sbxread(fn,1,1);            % read one frame to read the header of the image sequence
           % this contains the information about the structure of the image
        tic
        if info.volscan
            if isempty(info.otwave_um) % temporary remedy 2017/0714 JK
                num_plane = 30;
            else
                num_plane = length(info.otwave_um);
            end
%                 info.max_idx = 

            for j = 1 : num_plane
                temp_start = num_plane + j -1; % to discard first frames of each plane. First frame is 0th. 
                frame_dif = mod(info.max_idx+1,num_plane) - mod(j,num_plane);
                if frame_dif < 0
                    temp_end = info.max_idx - (num_plane + frame_dif);
                else
                    temp_end = info.max_idx - frame_dif;
                end
                [m{j}, T{j}] = jksbxalignx(fn,temp_start:num_plane:temp_end);
                T{j} = [0 0;T{j}]; % To compensate for ignoring the first frame
            end
        else
            [m,T] = jksbxalignx(fn,0:info.max_idx);   %
        end
        save([fn '.align'],'m','T');
        clear m T
        display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx,round(toc/60)));
    end
end
