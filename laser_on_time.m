% function laser_on_time(sbxfn, varargin)

% frames and lines of laser ON of each trial

%% temporary for test
% sbxfn = '025_004_000';
padding = 20;

% if naragin > 1 && isnumeric(varargin{1}) && varargin{1} == 2 && size(a,1) == 2
%     ch = varagin{1};
% else
    ch = 1; % default channel
% end

%% main

for si = 1 : 25
    sbxfn = sprintf('025_%03d_000',si);
    %% set up 
    sbxfn = strtok(sbxfn, '.');
    trialfn = [sbxfn, '.trials'];
    if ~exist(trialfn, 'file')    
        jksbxsplittrial(sbxfn)
    end
    load(trialfn, '-mat') % loading trials, trial_frames, frame_to_use, num_plane

    global info
    a = sbxread(sbxfn,0,1);
    if isfield(info, 'blankstart') % blankstart is set manually. Sometimes during file transfer using windows, the files get breached and turns into white blank frames. 2018/03/03 JK
        info.max_idx = info.blankstart-1;
    end
    
    ondelay{si} = zeros(length(trials),1);
    offdelay{si} = zeros(length(trials),1);

    for ti = 1 : length(trials) % trial index
    % for ti = 11
        testFrames = trials(ti).frames(1) : trials(ti).frames(2) + padding;
        tempF = sbxread(sbxfn,testFrames(1),length(testFrames));
        tempF = squeeze(mean(mean(squeeze(tempF(ch,:,:,:)))));    

        threshold = (mean(tempF) - std(tempF) + min(tempF))/2;
    %     figure, plot(testFrames,tempF, 'k.'), hold on, plot(testFrames, repmat(threshold,length(tempF),1), 'b-')
    %     set(gca,'xticklabel',num2str(get(gca,'xtick')','%d'))

        laserOFFind = find(tempF < threshold);
        if ismember(1, laserOFFind)
            i = 1;
            while true
                i = i + 1;
                if ~ismember(i, laserOFFind)
                    i = i - 1;
                    break
                end
            end
            first = i;
        else
            first = 0;
        end

        if ismember(length(tempF), laserOFFind)
            i = length(tempF);
            while true
                i = i - 1;
                if ~ismember(i, laserOFFind)
                    i = i + 1;
                    break
                end
            end
            last = i;
        else
            last = length(tempF) + 1;
        end
        last = last - (trials(ti).frames(2) - trials(ti).frames(1) + 1);
        ondelay{si}(ti) = first;
        offdelay{si}(ti) = last;
    end
end
%%
% figure, subplot(1,2,1), histogram(ondelay), subplot(122), histogram(offdelay)