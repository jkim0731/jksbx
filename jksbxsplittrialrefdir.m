function jksbxsplittrialdir(varargin)

d = dir('*.sbx');

if(nargin > 1) % cell with filenames to be aligned
    reffn = varargin{1};
end

for(i=1:length(d))
    try
        fn = strtok(d(i).name,'.');
        jksbxsplittrial(fn)
    catch
        display(sprintf('Could not split %s',fn))
    end
end