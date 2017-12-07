function jksbx2tif_forsuite2p(fn, varargin)

tic
if nargin > 1
    ch = varargin{1};
else
    ch = 1;
end
curr_dir = pwd;
global info
sbxread(fn,0,1);
if ~exist(fn,'dir')
    mkdir(fn)
end
cd(fn)
for i = 1 : length(info.aligned.m)
    if ~exist(num2str(i),'dir')
        mkdir(num2str(i))
    end
end
cd(curr_dir)

for i = 1 : length(info.aligned.m)
    for j = 1 : length(info.aligned.frame_to_align{i})
        a = sbxread(fn,info.aligned.frame_to_align{i}(j),1);
        a = squeeze(a(ch,101:end,151:end));
        temp_fn = sprintf('%s_%d_%05d.tif',fn,i,j);
        imwrite(a,[pwd,filesep,fn,filesep,num2str(i),filesep,temp_fn],'tif');
    end            
end    
toc