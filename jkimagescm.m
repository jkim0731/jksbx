function jkimagescm(fn,varargin)

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

figure, hold on,
% [mouse, remain] = strtok(fn,'_');
% pod = strtok(remain,'_');
% title(['JK', mouse, ', POD ', pod],'FontSize',15')
if iscell(m)
    if length(ind) == 1
        imagesc(m{ind}), axis image, axis off
    else
        nrow = floor(sqrt(length(ind)));
        ncol = ceil(length(ind)/nrow);
        for i = 1 : length(ind)
            subplot(nrow, ncol, i), imagesc(m{ind(i)}), axis image, axis off
        end
    end
else
    imagesc(m), axis image, axis off
end
