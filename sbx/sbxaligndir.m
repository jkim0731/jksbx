function sbxaligndir(varargin)

if(nargin>=1) % cell with filenames to be aligned
    for(i=1:length(varargin{1}))
        d(i).name = varargin{1}{i};
    end
else
    d = dir('*.sbx');
end

% Align all *.sbx files in the list

for(i=1:length(d))
    try
        fn = strtok(d(i).name,'.');
        if(exist([fn '.align'])==0)
            sbxread(fn,1,1);            % read one frame to read the header of the image sequence
            global info;                % this contains the information about the structure of the image
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
                    [m{j}, T{j}] = sbxalignx(fn,temp_start:num_plane:temp_end);
                    T{j} = [0 0;T{j}]; % To compensate for ignoring the first frame
                end
            else
                [m,T] = sbxalignx(fn,0:info.max_idx);   %
            end
            save([fn '.align'],'m','T');
            display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx,round(toc/60)));
        else
            sprintf('File %s is already aligned',fn)
        end
    catch
        sprintf('Could not align %s',fn)
    end
end
