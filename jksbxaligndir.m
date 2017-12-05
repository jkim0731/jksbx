function jksbxaligndir(varargin)
clear info
global info

if(nargin>=1) % cell with filenames to be aligned
    if isempty(varargin{1})
        d = dir('*.sbx');
    else
        for(i=1:length(varargin{1}))
            d(i).name = varargin{1}{i};
        end
    end
else
    d = dir('*.sbx');
end

if nargin >=2
    ref = varargin{2}; % 'red' or 'green'
else
    ref = '';
end

% Align all *.sbx files in the list

for i = 1:length(d)
    fn = strtok(d(i).name,'.');
    temp = sbxread(fn,0,1); 
    if exist([fn,'.align'],'file')
        sprintf('File %s is already aligned',fn);
        continue
    elseif exist([fn,'.trials'],'file') % ignore frames outside of each trial, because those are blank. If not treated, these lead to weird interference pattern during bidirectional scanning.
        load([fn,'.trials'],'-mat')
        trial_frames = [];
        for tn = 1 : length(trials)
            trial_frames = [trial_frames, trials(tn).frames(1)+1:trials(tn).frames(2)-1]; % ignore the first and the last frames of each trial (because of possible clipping by laser blanking and re-opening)
        end
    else % no trial file
        trial_frames = 0:info.max_idx;
    end

    tic
    if info.volscan
        num_plane = length(info.otwave_um);
    else
        num_plane = 1;
    end        
    m1 = cell(1,num_plane); % for green. Empty if only pmt1 was used.
    m2 = cell(1,num_plane); % for red. Empty if only pmt0 was used.        
    T1 = cell(1,num_plane); % Empty if only pmt1 was used.
    T2 = cell(1,num_plane); % Empty if only pmt0 was used.
    
    frame_to_align = cell(1,num_plane); % frames subject for alignment, based on the actual total frames (no pre-allocating for each plane)    
    max_idx = info.max_idx;
    channels = info.channels;
    
    blockimaging = 0;            
    num_layer = 0;
    num_block = [];
    block_message_i = [];
    for message_i = 1 : length(info.messages)
        if ~isempty(strfind(info.messages{message_i},'objective')) % for now, assume that number of planes is same across different blocks 
            blockimaging = 1;                
            num_block = [num_block, str2double(info.messages{message_i}(end))]; % number of layer on each imaging block. (1, 2, 3, ..., num_layer, 1, 2, 3, ..., num_layer, 1, 2, 3, ...)
            num_layer = max([num_layer, str2double(info.messages{message_i}(end))]); % total num of layer
            block_message_i = [block_message_i, message_i]; % index of info.messages written with imaging layer. (like a split marker)
        end
    end
    
    if blockimaging
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
                            trials(kk).frames(1)+num_plane : trials(kk).frames(2)-num_plane];% ignore first and last frame of each plane, each trial (because of possible clipping by laser blanking and re-opening)
                        break
                    end
                end
            end
        end                

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Very important variable
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        frame_to_align = cell(1,num_layer*num_plane); % this is going to be used for the rest of the analysis.        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ind_layer = 1 : num_layer                        
            for ind_plane = 1 : num_plane
                frame_to_align{(ind_layer-1)*num_plane + ind_plane} = intersect(ind_plane-1:num_plane:max_idx, trial_frames{ind_layer});                  
            end
        end
% %                 [m{ind_layer*num_plane}, T{ind_layer*num_plane}] = jksbxalignx(fn,frame_to_align{ind_layer*num_plane});

% % %     %                             movements = sqrt(diff(T{temp_ind}(:,1)).^2 + diff(T{temp_ind}(:,2)).^2);
% % %     %                             movements_std = std(movements);                        
% % %     %                             if ~isempty(find(movements > 2*movements_std,1))
% % %     %                                 fix_ind = find(movements > 2 * movements_std) + 1; % because of diff
% % %     %                                 for jj = 1 : length(fix_ind)
% % %     %                                     if fix_ind(jj) < length(movements)
% % %     %                                         T{temp_ind}(fix_ind(jj),:) = round((T{temp_ind}(fix_ind(jj)-1,:) + T{temp_ind}(fix_ind(jj)+1,:)) / 2);
% % %     %                                     else
% % %     %                                         T{temp_ind}(fix_ind(jj),:) = T{temp_ind}(fix_ind(jj)-1,:);
% % %     %                                     end
% % %     %                                 end
% % %     %                                 m{temp_ind} = zeros(size(m{temp_ind}));
% % %     %                                 for jj = 1 : length(T{temp_ind})                                
% % %     %                                     m{temp_ind} = m{temp_ind} + circshift(double(squeeze(sbxread(fn,jj-1,1)))/length(T{temp_ind}),T{temp_ind}(jj,:));
% % %     %                                 end
% % %     %                             end

    else % if not block_imaging
        for j = 1 : num_plane
            temp_start = num_plane -1; 
            frame_dif = mod(max_idx+1,num_plane) - mod(j,num_plane);
            if frame_dif < 0
                temp_end = max_idx - (num_plane + frame_dif);
            else
                temp_end = max_idx - frame_dif;
            end
            frame_to_align{j} = intersect(temp_start:num_plane:temp_end,trial_frames);             
        end        
    end
    
    for j = 1 : length(frame_to_align)
        if channels > 1 % 2 for pmt0 (green), 3 for pmt1 (red)
            if channels == 2
                [m1{j}, T1{j}] = jksbxalignx(fn,frame_to_align{j});
            else % info.channels == 3
                [m2{j}, T2{j}] = jksbxalignx(fn,frame_to_align{j});
            end
        else % 2 channel imaging
            if strcmp(ref,'green')
                [m1{j}, T1{j}] = jksbxalignx(fn,frame_to_align{j});
                T2{j} = T1{j};  
                m2{j} = align2pre(fn,T2{j},frame_to_align{j},2);                     
            elseif strcmp(ref,'red')
                [m2{j}, T2{j}] = jksbxalignx(fn,frame_to_align{j}, 2);
                T1{j} = T2{j};  
                m1{j} = align2pre(fn,T1{j},frame_to_align{j},1); 
            elseif strcmp(ref,'no') % no alignment at all. Use blank T arrays
                T1{j} = zeros(length(frame_to_align{j}),2);    T2{j} = T1{j};
                m1{j} = align2pre(fn,T1{j},frame_to_align{j},1);
                m2{j} = align2pre(fn,T2{j},frame_to_align{j},2);
            else % empty ref, meaning each channel is aligned by its own - alignment can be differ between channels, but could be used to compare between those 2 channels
                [m1{j}, T1{j}] = jksbxalignx(fn,frame_to_align{j});
                [m2{j}, T2{j}] = jksbxalignx(fn,frame_to_align{j}, 2);                    
            end
        end
    end
    save([fn '.align'],'m1','m2','T1','T2','frame_to_align');        
    display(sprintf('Done %s: Aligned %d images in %d min',fn,max_idx,round(toc/60)));        
end
end

function m = align2pre(fn, T, frame_to_align, ch)
    a = sbxread(fn,0,1);
    m = zeros(size(a,2),size(a,3));
    for Tind = 1 : size(T,1)
        temp = sbxread(fn, frame_to_align(Tind), 1);
        temp = double(squeeze(temp(ch,:,:)));
        m = m + circshift(temp,T(Tind,:))/size(T,1);
    end
    m = uint16(m);
end