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
    frame_to_align = cell(1,num_plane);
    fta_ind = cell(1,num_plane);
    max_idx = info.max_idx;
    channels = info.channels;
    parfor j = 1 : num_plane
        temp_start = num_plane -1; 
        frame_dif = mod(max_idx+1,num_plane) - mod(j,num_plane);
        if frame_dif < 0
            temp_end = max_idx - (num_plane + frame_dif);
        else
            temp_end = max_idx - frame_dif;
        end
        [frame_to_align{j}, fta_ind{j}] = intersect(temp_start:num_plane:temp_end,trial_frames);            
        if channels > 1 % 2 for pmt0 (green), 3 for pmt1 (red)
            if channels == 2
                [m1{j}, T1{j}] = temp_align(fn,frame_to_align{j},fta_ind{j}, length(temp_start:num_plane:temp_end));
            else % info.channels == 3
                [m2{j}, T2{j}] = temp_align(fn,frame_to_align{j},fta_ind{j}, length(temp_start:num_plane:temp_end));
            end
        else % 2 channel imaging
            if strcmp(ref,'green')
                [m1{j}, T1{j}] = temp_align(fn,frame_to_align{j},fta_ind{j}, length(temp_start:num_plane:temp_end));
                T2{j} = T1{j};  
                m2{j} = align2pre(fn,T2{j},j,num_plane,2);                     
            elseif strcmp(ref,'red')
                [m2{j}, T2{j}] = temp_align(fn,frame_to_align{j},fta_ind{j}, length(temp_start:num_plane:temp_end), 2);
                T1{j} = T2{j};  
                m1{j} = align2pre(fn,T1{j},j,num_plane,1); 
            elseif strcmp(ref,'no') % no alignment at all. Use blank T arrays
                T1{j} = zeros(length(frame_to_align{j}),2);    T2{j} = T1{j};
                m1{j} = align2pre(fn,T1{j},j,num_plane,1);
                m2{j} = align2pre(fn,T2{j},j,num_plane,2);
            else % empty ref, meaning each channel is aligned by its own - alignment can be differ between channels, but could be used to compare between those 2 channels
                [m1{j}, T1{j}] = temp_align(fn,frame_to_align{j},fta_ind{j}, length(temp_start:num_plane:temp_end));
                [m2{j}, T2{j}] = temp_align(fn,frame_to_align{j},fta_ind{j}, length(temp_start:num_plane:temp_end), 2);                    
            end
        end
    end
    save([fn '.align'],'m1','m2','T1','T2','frame_to_align','fta_ind');        
    display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx,round(toc/60)));        
end
end

function [m, T] = temp_align(fn, frame_to_align, fta_ind, T_length, varargin)
    if nargin > 4
        ch = varargin{1};
    else
        ch = 1;
    end
    [m, T_temp] = jksbxalignx(fn,frame_to_align,ch);
    T = zeros(T_length,2);
    for Tind = 1 : length(fta_ind)
        T(fta_ind(Tind),:) = T_temp(Tind,:);
    end    
end

function m = align2pre(fn, T, j, num_plane, ch)
    a = sbxread(fn,0,1);
    m = zeros(size(a,2),size(a,3));
    for Tind = 1 : size(T,1)
        temp = sbxread(fn, j-1 + (Tind-1)*num_plane, 1);
        temp = double(squeeze(temp(ch,:,:)));
        m = m + circshift(temp,T(Tind,:))/size(T,1);
    end
    m = uint16(m);
end