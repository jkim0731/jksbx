function jkshowm(fn,varargin)

[a, b] = strtok(fn,'.');
if isempty(b)
    fn = [fn, '.align'];
elseif ~strcmp(b, '.align')
    fn = [a,'.align'];
end

try
    cd(['D:\2p\JK\', strtok(fn,'_')])
    load(fn,'-mat')
catch
    try 
        cd(['D:\TPM\', strtok(fn,'_')])
        load(fn,'-mat')
    catch
        try 
            cd(['Y:\JK_temp\2p\', strtok(fn,'_')])
            load(fn,'-mat')
        catch
            error('Wrong directory')
        end
    end
end
ind = 1;
if nargin > 1
    ind = varargin{1};    
end
if nargin > 2 && strcmp(varargin{2},'montage')
    montage_flag = 1;
else
    montage_flag = 0;
end

figure, hold on, 
% [mouse, remain] = strtok(fn,'_');
% pod = strtok(remain,'_');
% title(['JK', mouse, ', POD ', pod],'FontSize',15)
if iscell(m)
    if length(ind) == 1
        imshow(m{ind})
    else
        if montage_flag
            I = zeros([size(m{1}),1,length(ind)]);
            for i = 1 : length(ind)
                I(:,:,1,i) = m{ind(i)};
            end
            montage(mat2gray(I))
        else
            nrow = floor(sqrt(length(ind)));
            ncol = ceil(length(ind)/nrow);
            for i = 1 : length(ind)
                subplot(nrow, ncol, i), imshow(m{ind(i)}), axis image, axis off
            end
        end
    end
else
    imshow(m), axis image, axis off
end
