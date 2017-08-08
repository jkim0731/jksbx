function jkimagescm(fn,varargin)

cd(['D:\2p\JK\', strtok(fn,'_')])
ind = 1;
if nargin > 1
    ind = varargin{1};
end
[a, b] = strtok(fn,'.');
if isempty(b)
    fn = [fn, '.align'];
elseif ~strcmp(b, '.align')
    fn = [a,'.align'];
end
load(fn,'-mat')
figure,
if iscell(m)
    if length(ind) == 1
        imagesc(m{ind})
    else
        nrow = floor(sqrt(length(ind)));
        ncol = ceil(length(ind)/nrow);
        for i = 1 : length(ind)
            subplot(nrow, ncol, i), imagesc(m{ind(i)})
        end
    end
else
    imagesc(m)
end