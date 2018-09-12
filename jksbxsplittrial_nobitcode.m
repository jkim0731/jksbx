function jksbxsplittrial_nobitcode(fn, varargin)

% same as jksbxsplittrial, but when there was no bitcode sent to scanbox
% (happened in experiments before 2018/02 during piezo stimulation)
% Just use info.messages

% splits an image (usually 30~60 min long) into each trials (use aligned-to-reference)
% number each trials with bitcode

%% Updates needed
% Laser on time point (frame & line) for each trial should be calculated
% (based on the fluorescence, maybe differecen from average image -
% randomly selected ~ 500 frames within trials)

%% Updates

% 2017/05/18: TTL1 ("durid" in the state matrix) is up during the whole trial (starts right after state
% 35, which has laser OFF and ON with 1+ sec blanking in-between and also
% whisker video trigger), and TTL0 ("slid" in the state matrix) is for
% trial number bitcode. TTL0 starts with stateBC (first line of the state
% matrix), so the start of a trial is with TTL Type 3, and the end of a
% trial is with TTL Type 2 (TTL1 going down, which means "durid" is off).

% % 1 sec before the pole rise, 1 sec after the pole down (32 frames)
% 2016/08/25: TTL2 is up for the entire trial except for state 35, matching
% with whisker video. -> Starting with video ON, end with next video ON
% (laser is on all the time, for now. Considering blanking between trials,
% only if the continuous laser on seems to affect viability)

% First gap time depends on TimeupSt of sBC in the state matrix. Bit time
% and gap time depends on those in "make_and_upload_state_matrix.m" in Behavior 
% Protocol. For 2port_ang_dist now, it's 10, 2, and 5 ms. 

% 2018/07/09
% including onFrames. For spontaneous and piezo imaging ~ 2018/06.
% splittrial includes laser off frames, so onFrames can be included as an
% input, and trials only within this onFrames are counted.
    %% check if already splitted or not
%     if exist([fn,'.trials'],'file')
%         fprintf('%s has already been split\n',fn)
%         return
%     else    

    %% input arguments
    if nargin > 1
        if isnumeric(varargin{1}) && min(varargin{1}) > 0 
            onFrames = varargin{1};
            laserOffIncluded = 1;
        else
            error('Input argument is either not numeric or includes negative values')
        end
    else
        laserOffIncluded = 0;
    end
    %% index sorting    
        load([fn,'.mat']);
        clear info
        a = squeeze(jksbxread(fn,0,1));
        global info    
        
        if isfield(info,'event_id') && ~isempty(info.event_id)
            % info.frame has limit at 2^16. Correct this
            % Don't save this for now. 2017/06/20 JK
            if info.max_idx > 2^16-1            
                if isfield(info,'frame')
                    overlimit_idx = find(diff(info.frame)<0);
                    for i = 1 : length(overlimit_idx)
                        if i == length(overlimit_idx)
                            info.frame(overlimit_idx(i)+1:end) = info.frame(overlimit_idx(i)+1:end) + i * 2^16;
                        else
                            info.frame(overlimit_idx(i)+1:overlimit_idx(i+1)) = info.frame(overlimit_idx(i)+1:overlimit_idx(i+1)) + i * 2^16;
                        end                    
                    end
                end
            end
            
            if isfield(info, 'blankstart') % blankstart is set manually. Sometimes during file transfer using windows, the files get breached and turns into white blank frames. 2018/03/03 JK
                info.max_idx = info.blankstart-1;
            end

            % event_id 3 (trial started - state 40) should follow directly after event_id 2
            % (end of a trial - state 35), but there are some error and event_id 3 is split
            % into two events of 1 and 2. Compensate.
            for i = 1 : size(info.event_id,1)-2
                if info.event_id(i) == 2
                    if info.event_id(i+1) == 1
                        if info.event_id(i+2) == 2
                            dframe = info.frame(i+2) - info.frame(i+1);
                            dline = info.line(i+2) - info.line(i+1) + info.sz(1) * dframe;
                            if dline < 3 % ~ 1 ms tolerance in 31Hz 512 line imaging
                                info.event_id(i+1) = 0;
                                info.event_id(i+2) = 3;
                            end
                        end
                    end
                end
            end
            % ? Need to consider 2 comes before 1 ?

            start_event = find(info.event_id == 3); % pole up is linked to both ttl0 & ttl1 up, making the event as 3. Refer to make_and_upload_state_matrix.m in the behavior protocol.
            end_event = find(info.event_id == 2); % pole down is linked to ttl1 down only. 

            while start_event(1) > end_event(1) % exception error for when the event started with pole_down
                end_event = end_event(2:end);
            end 
            while start_event(end) > end_event(end) % exception error for when the event ended with pole_up
                start_event = start_event(1:end-1);
            end

            num_event = min(length(start_event),length(info.messages)); % now that all the events were matched with pole_up/pole_down events, # of start_events are same as # of events                                      

            %% reading and saving
            trials = struct('trialnum',[],'frames',[], 'lines', []);
            for i = 1:num_event
                trials(i).trialnum = str2double(info.messages{i});
                trials(i).frames = [info.frame(start_event(i)),info.frame(end_event(i))];
                trials(i).lines = [info.line(start_event(i)),info.line(end_event(i))];
            end         
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [~,plane_sorted] = sort(info.otwave,'descend'); % sorting from the top. 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Assume objective is already sorted descending. (objective 1 higher, i.e., shallower, than objective 2)
            % Overall goal is to have all planes (including layers) sorted in descending order
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            if info.volscan
                num_plane = length(info.otwave_um);
            else
                num_plane = 1;
                plane_sorted = 1;
            end        
            max_idx = info.max_idx;
            blockimaging = 0;            
            num_layer = 1;
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
                frames_beginning = 0:num_plane:max_idx;
                frames_ending = num_plane-1:num_plane:max_idx;
                for ii = 1 : num_layer
                    trial_frames{ii} = []; 
                    for jj = 1 : length(layer_trials{ii})
                        for kk = 1 : length(trials)
                            if trials(kk).trialnum == layer_trials{ii}(jj)                                
                                begin_frame = frames_beginning(find(frames_beginning > trials(kk).frames(1), 1, 'first'));
                                end_frame = frames_ending(find(frames_ending < trials(kk).frames(2), 1, 'last'));
                                if laserOffIncluded
                                    currTrialFrames = intersect(onFrames, begin_frame:end_frame);
                                else
                                    currTrialFrames = begin_frame:end_frame;
                                end
                                trial_frames{ii} = [trial_frames{ii}, currTrialFrames];                                
%                                 trial_frames{ii} = [trial_frames{ii}, trials(kk).frames(1) : trials(kk).frames(2)];
%                                 Changed to include only the frames with full-FOV recording, and matching number of frames in each plane at the same layer in each trials. 2018/03/07 JK.
                                break
                            end
                        end
                    end
                end                

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Very important variable
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                frame_to_use = cell(1,num_layer*num_plane); % this is going to be used for the rest of the analysis.        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%

                for ind_layer = 1 : num_layer                        
                    for ind_plane = 1 : num_plane
                        frame_to_use{(ind_layer-1)*num_plane + plane_sorted(ind_plane)} = intersect(ind_plane-1:num_plane:max_idx, trial_frames{ind_layer});                  
                    end
                end

            else % if not block_imaging
                trial_frames = [];
                frame_to_use = cell(num_plane,1);
                frames_beginning = 0:num_plane:max_idx;
                frames_ending = num_plane-1:num_plane:max_idx;
                for trial_i = 1 : length(trials)
                    begin_frame = frames_beginning(find(frames_beginning > trials(trial_i).frames(1), 1, 'first'));
                    end_frame = frames_ending(find(frames_ending > trials(trial_i).frames(2), 1, 'last'));                    
                    if laserOffIncluded
                        currTrialFrames = intersect(onFrames, begin_frame:end_frame);
                    else
                        currTrialFrames = begin_frame:end_frame;
                    end
                    trial_frames = [trial_frames, currTrialFrames];
                end
                for ind_plane = 1 : num_plane
                    frame_to_use{plane_sorted(ind_plane)} = intersect(plane_sorted(ind_plane):num_plane:max_idx,trial_frames);             
                end        
            end
        else
            layer_trials = [];
            trials = []; 
            blockimaging = 0; num_layer = 1;
            if isfield(info, 'blankstart') % blankstart is set manually. Sometimes during file transfer using windows, the files get breached and turns into white blank frames. 2018/03/03 JK
                info.max_idx = info.blankstart-1;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Assume objective is already sorted descending. (objective 1 higher, i.e., shallower, than objective 2)
            % Overall goal is to have all planes (including layers) sorted in descending order
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            if info.volscan                
                num_plane = length(info.otwave_um);
                [~,plane_sorted] = sort(info.otwave,'descend'); % sorting from the top. 
                for ind_plane = 1 : num_plane
                    if laserOffIncluded
                        trial_frames = onFrames;
                        frame_to_use{plane_sorted(ind_plane)} = intersect(onFrames,(plane_sorted(ind_plane):num_plane:info.max_idx));
                    else
                        trial_frames = 0:info.max_idx;
                        frame_to_use{plane_sorted(ind_plane)} = (plane_sorted(ind_plane):num_plane:info.max_idx);
                    end 
                end 
            else
                num_plane = 1;
                frame_to_use = cell(num_plane,1);                
                if laserOffIncluded
                    trial_frames = onFrames;
                else
                    trial_frames = 0:info.max_idx;
                end
                frame_to_use{1} = trial_frames;
            end
        end
        clear global
%     end
    save([fn,'.trials'],'trials', 'frame_to_use', 'blockimaging', 'num_layer', 'num_plane', 'trial_frames');
end
