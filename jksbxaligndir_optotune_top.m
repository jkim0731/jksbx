function jksbxaligndir_optotune_top(varargin)
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
    if exist([fn,'.align*'],'file')
        sprintf('File %s is already aligned',fn);
    elseif exist([fn,'.trials'],'file') % ignore frames outside of each trial, because those are blank. If not treated, these lead to weird interference pattern during bidirectional scanning.
        load([fn,'.trials'],'-mat')
        trial_frames = [];
        sbxread(fn,0,1);            % read one frame to read the header of the image sequence
             % this contains the information about the structure of the image
        tic
        if info.volscan
            num_plane = length(info.otwave_um);
           
            blockimaging = 0;            
            num_layer = 0;
            num_block = [];
            block_message_i = [];
            for message_i = 1 : length(info.messages)
                if ~isempty(strfind(info.messages{message_i},'objective'))
                    blockimaging = 1;
                    num_block = [num_block, str2double(info.messages{message_i}(end))]; % number of layer on each imaging block. (1, 2, 3, ..., num_layer, 1, 2, 3, ..., num_layer, 1, 2, 3, ...)
                    num_layer = max([num_layer, str2double(info.messages{message_i}(end))]); % total num of layer
                    block_message_i = [block_message_i, message_i]; % index of info.messages written with imaging layer. (like a split marker)
                end
            end
            if blockimaging % for now, assume that number of planes is same across different blocks                
                if isempty(strfind(info.messages{1},'objective'))
                    temp_num_layer = mod(num_block(1)-1,num_layer);
                    if temp_num_layer == 0
                        temp_num_layer = num_layer;
                    end
                    num_block = [temp_num_layer, num_block]; % if the info.messages did not start with a layer indication, add the first layer at the beginning
                    block_message_i = [0, block_message_i]; % the beginning index is 0 now.
                end
                if isempty(strfind(info.messages{end},'objective'))
                    temp_num_layer = mod(num_block(end)+1,num_layer);
                    if temp_num_layer == 0
                        temp_num_layer = num_layer;
                    end
                    num_block = [num_block, temp_num_layer]; % if the info.messages did not end with a layer indication, add the next layer at the end
                    block_message_i = [block_message_i, length(info.messages) + 1]; % adding index at the end (1 larger index than the actual info.messages length)
                end
                
                layer_trials = cell(1,num_layer); % allocate trials to each layer
                for ii = 1 : num_layer
                    layer_trials{ii} = [];
                end
                for ii = 1 : length(block_message_i)-1
                    layer_trials{num_block(ii)} = [layer_trials{num_block(ii)}, str2double(info.messages{block_message_i(ii)+1}) : str2double(info.messages{block_message_i(ii+1)-1})];
                end

                trial_frames = cell(1,num_layer); % allocate frames to each layer
                for ii = 1 : length(trial_frames)
                    trial_frames{ii} = []; 
                    for jj = 1 : length(layer_trials{ii})
                        for kk = 1 : length(trials)
                            if trials(kk).trialnum == layer_trials{ii}(jj)
                                trial_frames{ii} = [trial_frames{ii}, ...
                                    trials(kk).frames(1)+num_plane : trials(kk).frames(2)-num_plane];% ignore first and last frame of each plane, each trial 
                                break
                            end
                        end
                    end
                end                
                
                m = cell(1,num_layer*num_plane);
                T = cell(1,num_layer);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Very important variable
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                frame_to_align = cell(1,num_layer*num_plane); % this is going to be used for the rest of the analysis.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for ind_layer = 1 : num_layer
                    frame_to_align{ind_layer*num_plane} = intersect(num_plane-1:num_plane:info.max_idx, trial_frames{ind_layer});
                    [m{ind_layer*num_plane}, T{ind_layer}] = jksbxalignx(fn,frame_to_align{ind_layer*num_plane});
                    for j = 1 : num_plane-1
                        temp_ind = (ind_layer-1)*num_plane + j;
                        frame_to_align{temp_ind} = frame_to_align{ind_layer*num_plane} - (num_plane - j);
                        m{temp_ind} = zeros(info.sz);
                        for k = 1 : length(T{ind_layer})
                            m{temp_ind} = m{temp_ind} + ...
                                double(circshift(squeeze(sbxread(fn, frame_to_align{temp_ind}(k), 1)), T{ind_layer}(k,:)))/length(T{ind_layer});
                        end
                    end
                end
            else % not block imaging    
                m = cell(1,num_plane);                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Very important variable
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                frame_to_align = cell(1,num_plane); % this is going to be used for the rest of the analysis.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                for tn = 1 : length(trials)
                    trial_frames = [trial_frames, trials(tn).frames(1)+num_plane:trials(tn).frames(2)-num_plane]; % ignore first and last frame of each plane, each trial 
                end
                frame_to_align{num_plane} = intersect(num_plane-1:num_plane:info.max_idx,trial_frames);
                [m{num_plane}, T] = jksbxalignx(fn,frame_to_align);
                for j = 1 : num_plane-1
                    frame_to_align{j} = frame_to_align{num_plane} - (num_plane - j);
                    m{j} = zeros(info.sz);
                    for k = 1 : length(T)
                        m{j} = m{j} + circshift(squeeze(sbxread(fn, frame_to_align{j}(k), 1)), T(k,:))/length(T);
                    end
                end
            end
        else
            for tn = 1 : length(trials)
                trial_frames = [trial_frames, trials(tn).frames(1)+4:trials(tn).frames(2)-4]; % ignore first and last 4 frames of each trial 
            end
            frame_to_align = intersect(0:info.max_idx,trial_frames);
            [m,T] = jksbxalignx(fn,frame_to_align);   %            
        end
        save([fn '.align_ot'],'m','T','frame_to_align');
        clear m T
        display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx,round(toc/60)));
    else % no trial file, which means it was not a behavior session
        sbxread(fn,1,1);            % read one frame to read the header of the image sequence
           % this contains the information about the structure of the image
        tic
        if info.volscan
            if isempty(info.otwave_um) % temporary remedy 2017/0714 JK
                num_plane = 30;
            else
                num_plane = length(info.otwave_um);
            end

            for j = num_plane
                temp_start = num_plane + j -1; % to discard first frames of each plane. First frame is 0th. 
                frame_dif = mod(info.max_idx+1,num_plane) - mod(j,num_plane);
                if frame_dif < 0
                    temp_end = info.max_idx - (num_plane + frame_dif);
                else
                    temp_end = info.max_idx - frame_dif;
                end
                [m{j}, T] = jksbxalignx(fn,temp_start:num_plane:temp_end);
                T = [0 0;T]; % To compensate for ignoring the first frame
            end
            for j = 1 : num_plane-1                
                m{j} = zeros(size(m{num_plane}));
                for k = 1 : length(T)
                    if k == 1
                        m{j} = circshift(squeeze(sbxread(fn, j + (k-1)*num_plane - 1, 1)), T(k,:))/length(T);
                    else
                        m{j} = m{j} + circshift(squeeze(sbxread(fn, j + (k-1)*num_plane - 1, 1)), T(k,:))/length(T);
                    end
                end
            end

        else
            [m,T] = jksbxalignx(fn,0:info.max_idx);   %
        end
        save([fn '.align'],'m','T');
        clear m T
        display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx,round(toc/60)));
    end
end
