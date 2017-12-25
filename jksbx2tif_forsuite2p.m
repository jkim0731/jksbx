function jksbx2tif_forsuite2p(fn, varargin)

tic

if ~exist([fn,'.trials'],'file')
    jksbxsplittrial(fn)
end
if ~exist([fn,'.align'],'file')
    jksbxaligndir({fn},'green')
end

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
if isfield(info.aligned, 'm')
    m = info.aligned.m;
elseif isfield(info.aligned, 'm1')
    m = info.aligned.m1;
else
    error('no aligned average image')
end
for i = 1 : length(m)    
    if ~exist(num2str(i),'dir')
        mkdir(num2str(i))
    end
end
cd(curr_dir)

fta = info.aligned.frame_to_align;

for i = 1 : length(m)
    parfor j = 1 : length(fta{i})
        a = sbxread(fn,fta{i}(j),1);
        a = squeeze(a(ch,101:end,101:end));
        temp_fn = sprintf('%s_%d_%05d.tif',fn,i,j);
        imwrite(a,[pwd,filesep,fn,filesep,num2str(i),filesep,temp_fn],'tif');
    end            
end    
toc