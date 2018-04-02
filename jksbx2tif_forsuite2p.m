function jksbx2tif_forsuite2p(fn, varargin)

tic

if ~exist([fn,'.trials'],'file')
    jksbxsplittrial(fn)
end

load([fn,'.trials'],'-mat') % loading 'frame_to_use'
ftu = frame_to_use;

if nargin > 1
    ch = varargin{1};
else
    ch = 1;
end

if nargin > 2
    fnummax = varargin{2};
else
    fnummax = 4000;
end

curr_dir = pwd;
global info
sbxread(fn,0,1);
if ~exist(fn,'dir')
    mkdir(fn)
end
cd(fn)
for i = 1 : length(frame_to_use)    
    if ~exist(num2str(i),'dir')
        mkdir(num2str(i))
    end
end

cd(curr_dir)

for i = 1 : length(ftu)
    for j = 1 : length(ftu{i})        
        a = sbxread(fn,ftu{i}(j),1);
        a = squeeze(a(ch,101:end,101:end));
        if ch == 1
            temp_fn = sprintf('%s_%d_green_%02d.tif',fn,i,ceil(j/fnummax));
        elseif ch == 2
            temp_fn = sprintf('%s_%d_red_%02d.tif',fn,i,ceil(j/fnummax));
        end        
        imwrite(a,[pwd,filesep,fn,filesep,num2str(i),filesep,temp_fn],'WriteMode', 'append');
    end            
end    

toc
