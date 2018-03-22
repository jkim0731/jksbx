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
        load([fn,'.trials'],'-mat', 'frame_to_use') 
    else % no trial file
        frame_to_use{1} = 0:info.max_idx;
    end

    tic
    
    if info.volscan
        num_plane = length(info.otwave_um);
    else
        num_plane = 1;
        plane_sorted = 1;
    end        
    
    m1 = cell(1,num_plane); % for green. Empty if only pmt1 was used.
    m2 = cell(1,num_plane); % for red. Empty if only pmt0 was used.        
    T1 = cell(1,num_plane); % Empty if only pmt1 was used.
    T2 = cell(1,num_plane); % Empty if only pmt0 was used.
    
    for j = 1 : length(frame_to_use)
        if info.channels > 1 % 2 for pmt0 (green), 3 for pmt1 (red)
            if info.channels == 2
                [m1{j}, T1{j}] = sbxalignx(fn,frame_to_use{j});
            else % info.channels == 3
                [m2{j}, T2{j}] = sbxalignx(fn,frame_to_use{j});
            end
        else % 2 channel imaging
            if strcmp(ref,'green')
                [m1{j}, T1{j}] = sbxalignx(fn,frame_to_use{j});
                T2{j} = T1{j};  
                m2{j} = align2pre(fn,T2{j},frame_to_use{j},2);                     
            elseif strcmp(ref,'red')
                [m2{j}, T2{j}] = jksbxalignx(fn,frame_to_use{j}, 2);
                T1{j} = T2{j};  
                m1{j} = align2pre(fn,T1{j},frame_to_use{j},1); 
            elseif strcmp(ref,'no') % no alignment at all. Use blank T arrays
                T1{j} = zeros(length(frame_to_use{j}),2);    T2{j} = T1{j};
                m1{j} = align2pre(fn,T1{j},frame_to_use{j},1);
                m2{j} = align2pre(fn,T2{j},frame_to_use{j},2);
            else % empty ref, meaning each channel is aligned by its own - alignment can be differ between channels, but could be used to compare between those 2 channels
                [m1{j}, T1{j}] = jksbxalignx(fn,frame_to_use{j});
                [m2{j}, T2{j}] = jksbxalignx(fn,frame_to_use{j}, 2);                    
            end
        end
    end
    save([fn '.align'],'m1','m2','T1','T2','frame_to_use');        
    fprintf('Done %s: Aligned %d images in %d min\n',fn,info.max_idx,round(toc/60));        
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