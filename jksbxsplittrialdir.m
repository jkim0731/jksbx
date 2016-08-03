function jksbxsplittrialdir(varargin)

if(nargin>=1) % cell with filenames to be split
    for(i=1:length(varargin{1}))
        d(i).name = varargin{1}{i};
    end
else
    d = dir('*.sbx');
end

for(i=1:length(d))
    try
        fn = strtok(d(i).name,'.');
        jksbxsplittrial(fn)
    catch
        display(sprintf('Could not split %s',fn))
    end
end