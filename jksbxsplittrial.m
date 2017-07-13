function jksbxsplittrial(fn)
% splits an image (usually 30~60 min long) into each trials (use aligned-to-reference)
% number each trials with bitcode

% % 1 sec before the pole rise, 1 sec after the pole down (32 frames)
% 2016/08/25: TTL2 is up for the entire trial except for state 35, matching
% with whisker video. -> Starting with video ON, end with next video ON
% (laser is on all the time, for now. Considering blanking between trials,
% only if the continuous laser on seems to affect viability)

% 2017/05/18: TTL1 ("durid" in the state matrix) is up during the whole trial (starts right after state
% 35, which has laser OFF and ON with 1+ sec blanking in-between and also
% whisker video trigger), and TTL0 ("slid" in the state matrix) is for
% trial number bitcode. TTL0 starts with stateBC (first line of the state
% matrix), so the start of a trial is with TTL Type 3, and the end of a
% trial is with TTL Type 2 (TTL1 going down, which means "durid" is off).

% First gap time depends on TimeupSt of sBC in the state matrix. Bit time
% and gap time depends on those in "make_and_upload_state_matrix.m" in Behavior 
% Protocol. For 2port_ang_dist now, it's 10, 2, and 5 ms. 
    %% check if already splitted or not
    if exist([fn,'.trials'],'file')
        display(sprintf('%s has already been split',fn))
        return
    else    
    %% index sorting    
        load([fn,'.mat']);
        a = squeeze(jksbxread(fn,0,1));
        global info    
        
        if isfield(info,'event_id')
        if ~isempty(info.event_id)
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
            
            % event_id 3 (pole_up) should follow directly after event_id 2
            % (pole_down), but there are some error and event_id 3 is split
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
            
            if start_event(1) > end_event(1) % exception error for when the event started with pole_down
                end_event = end_event(2:end);
            end 
            if start_event(end) > end_event(end) % exception error for when the event ended with pole_up
                start_event = start_event(1:end-1);
            end
            
            num_event = length(start_event); % now that all the events were matched with pole_up/pole_down events, # of start_events are same as # of events
            bc_chunk_idx = cell(num_event,1); % ttl1 is for the bitcode of trialnum. 
            for i = 1 : num_event
                if end_event(i)-start_event(i)  <= 1
                    disp(['No bitcode arrival at event #' num2str(i) ' in ' fn '.sbx'])
                else
                    bc_chunk_idx{i} = start_event(i):end_event(i)-1; % must have at least 2 indices
                end                
            end                  

            %% reading and saving
            trials = struct('trialnum',[],'frames',[], 'lines', []);
            for i = 1:num_event
                trials(i).trialnum = read_bitcode(bc_chunk_idx{i}, 10, 2, 5);                
                if isempty(find(cellfun(@(x) strcmp(x,num2str(trials(i).trialnum)),info.messages),1))
                    disp(['Bitcode mismatch with message received by scanbox in event #' num2str(i) ' in ' fn '.sbx'])
                end
                trials(i).frames = [info.frame(start_event(i)),info.frame(end_event(i))];
                trials(i).lines = [info.line(start_event(i)),info.line(end_event(i))];
            end            
            save([fn,'.trials'],'trials');
        end
        end
    end
end
%% reading the bitcode
% First_gaptm, bittm (bit time), gaptm (gap time) from
% "make_and_upload_start_matrix.m". 
function trialnum = read_bitcode(bit_idx,first_gaptm,bittm,gaptm)
    global info
    line_rate = 1000/info.resfreq/(2-info.scanmode); % duration of imaging a line in ms
    n = length(bit_idx);
    dline = zeros(n-1,1);
    for i = 1 : n-1
        dframe = info.frame(bit_idx(i)+1) - info.frame(bit_idx(i));
        dline(i) = info.line(bit_idx(i)+1) - info.line(bit_idx(i)) + info.sz(1) * dframe; 
    end    
    dbit = dline*line_rate; 
    dbit(1) = dbit(1) - first_gaptm;
    dbit = round(dbit/(bittm + gaptm));
    invbitcode = [];
    for i = 1 : n-1    
        invbitcode = [invbitcode; zeros(dbit(i)-1,1);1];
    end

    trialnum = invbit2num(invbitcode);
end
%% function for converting inverse bitcode to a number
function num = invbit2num(invbitcode)
    num = 0;
    for ii = 1 : length(invbitcode)
        num = num + invbitcode(ii) * pow2(ii-1);
    end
end