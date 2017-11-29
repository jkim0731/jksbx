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

for(i=1:length(d))
    fn = strtok(d(i).name,'.');
    temp = sbxread(fn,0,1); 
    if isempty(ref)
        if exist([fn,'.align'],'file')
            sprintf('File %s is already aligned',fn);
        elseif exist([fn,'.trials'],'file') % ignore frames outside of each trial, because those are blank. If not treated, these lead to weird interference pattern during bidirectional scanning.
            load([fn,'.trials'],'-mat')
            trial_frames = [];
            for tn = 1 : length(trials)
                trial_frames = [trial_frames, trials(tn).frames(1)+1:trials(tn).frames(2)-1]; % ignore the first and the last frames of each trial
            end

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
            tic
            if info.volscan
                if isempty(info.otwave_um) % temporary remedy 2017/0714 JK
                    num_plane = 30;
                else
                    num_plane = length(info.otwave_um);
                end

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
    elseif strcmp(ref,'red')
        if size(temp,1) < 2
            error('The file should have 2 channels')
        else
            if exist([fn,'.align'],'file')
                sprintf('File %s is already aligned',fn);
            elseif exist([fn,'.trials'],'file') % ignore frames outside of each trial, because those are blank. If not treated, these lead to weird interference pattern during bidirectional scanning.
                load([fn,'.trials'],'-mat')
                trial_frames = [];
                for tn = 1 : length(trials)
                    trial_frames = [trial_frames, trials(tn).frames(1)+1:trials(tn).frames(2)-1]; % ignore the first and the last frames of each trial
                end

                tic
                if info.volscan
                    num_plane = length(info.otwave_um);
                    m1 = cell(1,num_plane);
                    m2 = cell(1,num_plane);
                    for j = 1 : num_plane
                        temp_start = num_plane + j -1; % to discard first frames of each plane. First frame is 0th. 
                        frame_dif = mod(info.max_idx+1,num_plane) - mod(j,num_plane);
                        if frame_dif < 0
                            temp_end = info.max_idx - (num_plane + frame_dif);
                        else
                            temp_end = info.max_idx - frame_dif;
                        end
                        [frame_to_align, fta_ind] = intersect(temp_start:num_plane:temp_end,trial_frames);
                        [m2{j}, T_temp] = jksbxalignx_red(fn,frame_to_align);
                        T2{j} = zeros(length(temp_start:num_plane:temp_end),2);
                        for Tind = 1 : length(fta_ind)
                            T2{j}(fta_ind(Tind),:) = T_temp(Tind,:);
                        end
                        T2{j} = [0 0;T2{j}]; % To compensate for ignoring the first frame
                        m1{j} = zeros(size(m2{j},1),size(m2{j},2));
                        for Tind = 1 : size(T2{j},1)
                            temp_green = sbxread(fn, j-1 + (Tind-1)*num_plane, 1);
                            temp_green = double(squeeze(temp_green(1,:,:)));
                            m1{j} = m1{j} + circshift(temp_green,T2{j}(Tind,:))/size(T2{j},1);
                        end
                    end
                else
                    [frame_to_align, fta_ind] = intersect(0:info.max_idx,trial_frames);
                    [m2,T_temp] = jksbxalignx_red(fn,frame_to_align);   %
                    T2 = zeros(info.max_idx+1,2);
                    for Tind = 1 : length(fta_ind)
                        T2(fta_ind(Tind),:) = T_temp(Tind,:);
                    end
                    
                    m1 = zeros(size(m2,1),size(m2,2));
                    for Tind = 1 : size(T2,1)
                        temp_green = sbxread(fn,Tind-1,1);
                        temp_green = double(squeeze(temp_green(1,:,:)));
                        m1 = m1 + circshift(temp_green,T2(Tind,:))/size(T2,1);
                    end
                end
                save([fn '.align'],'m1','m2','T2','frame_to_align','fta_ind');
                clear m1 m2 T2
                display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx,round(toc/60)));
            else        
                tic
                if info.volscan
                    num_plane = length(info.otwave_um);
                    m1 = cell(1,num_plane);
                    m2 = cell(1,num_plane);
                    for j = 1 : num_plane
                        temp_start = num_plane + j -1; % to discard first frames of each plane. First frame is 0th. 
                        frame_dif = mod(info.max_idx+1,num_plane) - mod(j,num_plane);
                        if frame_dif < 0
                            temp_end = info.max_idx - (num_plane + frame_dif);
                        else
                            temp_end = info.max_idx - frame_dif;
                        end
                        [m1{j}, T2{j}] = jksbxalignx_red(fn,temp_start:num_plane:temp_end);
                        T2{j} = [0 0;T2{j}]; % To compensate for ignoring the first frame
                        
                        m2{j} = zeros(size(m1{j},1),size(m1{j},2));
                        for Tind = 1 : size(T2{j},1)
                            temp_green = sbxread(fn, j-1 + (Tind-1)*num_plane, 1);
                            temp_green = double(squeeze(temp_green(1,:,:)));
                            m2{j} = m2{j} + circshift(temp_green,T2(Tind,:))/size(T2{j},1);
                        end                        
                    end
                else
                    [m1,T2] = jksbxalignx(fn,0:info.max_idx);
                    m2 = zeros(size(m1,1),size(m1,2));
                    for Tind = 1 : size(T2,1)
                        temp_green = sbxread(fn,Tind-1,1);
                        temp_green = double(squeeze(temp_green(1,:,:)));
                        m2 = m2 + circshift(temp_green,T2(Tind,:))/size(T2,1);
                    end
                end
                save([fn '.align'],'m1','m2','T2');
                clear m1 m2 T2
                display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx,round(toc/60)));
            end
        end
    elseif strcmp(ref,'green')
        if size(temp,1) < 2
            error('The file should have 2 channels')
        else
            if exist([fn,'.align'],'file')
                sprintf('File %s is already aligned',fn);
            elseif exist([fn,'.trials'],'file') % ignore frames outside of each trial, because those are blank. If not treated, these lead to weird interference pattern during bidirectional scanning.
                load([fn,'.trials'],'-mat')
                trial_frames = [];
                for tn = 1 : length(trials)
                    trial_frames = [trial_frames, trials(tn).frames(1)+1:trials(tn).frames(2)-1]; % ignore the first and the last frames of each trial
                end

                tic
                if info.volscan
                    num_plane = length(info.otwave_um);
                    m1 = cell(1,num_plane);
                    m2 = cell(1,num_plane);
                    for j = 1 : num_plane
                        temp_start = num_plane + j -1; % to discard first frames of each plane. First frame is 0th. 
                        frame_dif = mod(info.max_idx+1,num_plane) - mod(j,num_plane);
                        if frame_dif < 0
                            temp_end = info.max_idx - (num_plane + frame_dif);
                        else
                            temp_end = info.max_idx - frame_dif;
                        end
                        [frame_to_align, fta_ind] = intersect(temp_start:num_plane:temp_end,trial_frames);
                        [m1{j}, T_temp] = jksbxalignx(fn,frame_to_align);
                        T1{j} = zeros(length(temp_start:num_plane:temp_end),2);
                        for Tind = 1 : length(fta_ind)
                            T1{j}(fta_ind(Tind),:) = T_temp(Tind,:);
                        end
                        T1{j} = [0 0;T1{j}]; % To compensate for ignoring the first frame
                        
                        m2{j} = zeros(size(m1{j},1),size(m1{j},2));
                        for Tind = 1 : size(T1{j},1)
                            temp_red = sbxread(fn, j-1 + (Tind-1)*num_plane, 1);
                            temp_red = double(squeeze(temp_red(2,:,:)));
                            m2{j} = m2{j} + circshift(temp_red,T1{j}(Tind,:))/size(T1{j},1);
                        end
                    end
                else
                    [frame_to_align, fta_ind] = intersect(0:info.max_idx,trial_frames);
                    [m1,T_temp] = jksbxalignx(fn,frame_to_align);   %
                    T1 = zeros(info.max_idx+1,2);
                    for Tind = 1 : length(fta_ind)
                        T1(fta_ind(Tind),:) = T_temp(Tind,:);
                    end
                    
                    m2 = zeros(size(m1,1),size(m1,2));
                    for Tind = 1 : size(T1,1)
                        temp_red = sbxread(fn,Tind-1,1);
%                         temp_red = double(squeeze(temp_red(2,:,:)));
                        temp_red = squeeze(temp_red(2,:,:));
                        m2 = m2 + circshift(temp_red,T1(Tind,:))/size(T1,1);
                    end
                end
                save([fn '.align'],'m1','m2','T1','frame_to_align','fta_ind');
                clear m1 m2 T1
                display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx,round(toc/60)));
            else        
                tic
                if info.volscan
                    num_plane = length(info.otwave_um);
                    m1 = cell(1,num_plane);
                    m2 = cell(1,num_plane);
                    for j = 1 : num_plane
                        temp_start = num_plane + j -1; % to discard first frames of each plane. First frame is 0th. 
                        frame_dif = mod(info.max_idx+1,num_plane) - mod(j,num_plane);
                        if frame_dif < 0
                            temp_end = info.max_idx - (num_plane + frame_dif);
                        else
                            temp_end = info.max_idx - frame_dif;
                        end
                        [m1{j}, T1{j}] = jksbxalignx(fn,temp_start:num_plane:temp_end);
                        T1{j} = [0 0;T1{j}]; % To compensate for ignoring the first frame
                        
                        m2{j} = zeros(size(m1{j},1),size(m1{j},2));
                        for Tind = 1 : size(T1{j},1)
                            temp_red = sbxread(fn, j-1 + (Tind-1)*num_plane, 1);
                            temp_red = double(squeeze(temp_red(2,:,:)));
                            m2{j} = m2{j} + circshift(temp_red,T1{j}(Tind,:))/size(T1{j},1);
                        end                        
                    end
                else
                    [m1,T1] = jksbxalignx(fn,0:info.max_idx);
                    m2 = zeros(size(m1,1),size(m1,2));
                    for Tind = 1 : size(T1,1)
                        temp_red = sbxread(fn,Tind-1,1);
                        temp_red = double(squeeze(temp_red(2,:,:)));
                        m2 = m2 + circshift(temp_red,T1(Tind,:))/size(T1,1);
                    end
                end
                save([fn '.align'],'m1','m2','T1');
                clear m1 m2 T1
                display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx,round(toc/60)));
            end
        end
    else
        error('ref channel should be either red or green')
    end
        
end
