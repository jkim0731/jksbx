function jkshowm_sbxcom(fn,varargin)

cd(['D:\TPM\', strtok(fn,'_')])
ind = 1;
if nargin > 1
    ind = varargin{1};    
end
if nargin > 2 && strcmp(varargin{2},'montage')
    montage_flag = 1;
else
    montage_flag = 0;
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
        imshow(m{ind})
    else
        if montage_flag
            I = zeros([size(m{1}),1,length(ind)]);
            for i = 1 : length(ind)
                I(:,:,1,i) = m{ind(i)};
            end
            montage(I)
        else
            nrow = floor(sqrt(length(ind)));
            ncol = ceil(length(ind)/nrow);
            for i = 1 : length(ind)
                subplot(nrow, ncol, i), imshow(m{ind(i)})
            end
        end
    end
else
    imshow(m)
end