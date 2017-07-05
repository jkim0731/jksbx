function jksbxsplittrial(fn)
% splits an image (usually 30~60 min long) into each trials (use aligned-to-reference)
% number each trials with bitcode
% 1 sec before the pole rise, 1 sec after the pole down (32 frames)

    %% check if already splitted or not
    if exist([fn,'.trials'])
        display(sprintf('%s has already been split',fn))
        return
    else    
    %% index sorting
        a = squeeze(jksbxread(fn,0,1));
        global info

        load([fn,'.mat']);

        ttl2 = find(info.event_id == 2); % ttl2 = idx of every ttl2 event
        ttl2diff = ttl2(2:end) - ttl2(1:end-1); % ttl2diff = difference in ttl2 between neighboring idx
        ttl2_up_idx = find(ttl2diff == 1); % ttl2_up_idx = idx of every pole up movement
        ttl2_down_idx = ttl2_up_idx + 1; % ttl2_down_idx = idx of every pole down movement EXCEPT for a probable first down when it comes before first up.

        num_bitcode = length(find(ttl2diff(1:end-1)>1));
        ttl1chunk_idx = cell(num_bitcode,1); % ttl1 is for the bitcode of trialnum. 
        temp_numbitcode = 1;
        for i = 1 : length(ttl2diff)-1
            if ttl2diff(i) > 1
                for j = 1:ttl2diff(i)-1
                    ttl1chunk_idx{temp_numbitcode} = [ttl1chunk_idx{temp_numbitcode}, ttl2(i)+j]; % contains the information of i-th bitcode
                end
                temp_numbitcode = temp_numbitcode + 1;
            end
        end

            % disregard trials before the first down always. Because we cannot be sure about the
            % integrity of the first bitcode, and throwing away one trial will not be significant.
            if ttl2diff(1) == 1
                ttl2_up_idx = ttl2_up_idx(2:end);
                ttl2_down_idx = ttl2_down_idx(2:end);
            end 

        pole_up_frame = info.frame(ttl2(ttl2_up_idx));
        pole_down_frame = info.frame(ttl2(ttl2_down_idx));

        if pole_up_frame - 32*1 - 2 < 1
            pole_up_frame = pole_up_frame(2:end);
            pole_down_frame = pole_down_frame(2:end);
            ttl1chunk_idx = ttl1chunk_idx{2:end};
        end


        d = dir([fn,'.sbx']);
        info.max_frame =  d.bytes/info.recordsPerBuffer/info.sz(2)/4; % calculating total frame number. This might be modified (copied and modified from sbxread). This one starts with 1, so to read, it should be sbxread(fn,info.max_frame -1, 1)

        if pole_down_frame + 32*1 > info.max_frame
            pole_up_frame = pole_up_frame(1:end-1);
            pole_down_frame = pole_down_frame(1:end-1);
            ttl1chunk_idx = ttl1chunk_idx{1:end-1};
        end

        first_frames = pole_up_frame - 32*1 -2; % 30.98 Hz, 32 for safety + 2 more frames for Kalman;
        last_frames = pole_down_frame + 32 * 1;

        if (length(pole_up_frame) ~= length(pole_down_frame)) || length(pole_up_frame) ~= length(ttl1chunk_idx)
            error('index mismatch!')
            return
        end

        %% reading and saving
        trials = struct('trialnum',[],'frames',[]);
        for i = 1:length(first_frames)
            tempfirstframe = first_frames(i);
            templastframe = last_frames(i);
            trialnum = read_bitcode(ttl1chunk_idx{i});
            trials(i).trialnum = trialnum;
            trials(i).frames = [first_frames(i):last_frames(i)];
        end
        if ~isempty(trials.trialnum) && ~isempty(trials.frames)
            save([fn,'.trials'],'trials');
        else
            disp('no available trials')
        end
    end
end
%% reading the bitcode
% assume every image is 31 Hz, so one bit interval is approximately 114
% lines
% Least value bit arrives first
% maximum bit-length is 10 bit (1024), and the maximum bitcode-length is
% therefore 11 because...
% first one bit is for calculating timing
function trialnum = read_bitcode(bit_idx)
    global info
    n = length(bit_idx);
    dframe = zeros(n-1,1);
    dline = zeros(n-1,1);
    for i = 1 : n-1
        dframe(i) = info.frame(bit_idx(i)+1) - info.frame(bit_idx(i));
        dline(i) = info.line(bit_idx(i)+1) - info.line(bit_idx(i)) + info.sz(1) * dframe(i); 
    end
    dbit = round(dline/114); % this number 114 is empirical. there are variations.
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

