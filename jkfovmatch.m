
% pwd

ref_fn = '650_018_000';
load([ref_fn,'.align'], '-mat')

if iscell(m)
    m = m{4};
end

mmfile = memmapfile('scanbox.mmap','Writable',true, ...
    'Format', { 'int16' [1 16] 'header' } , 'Repeat', 1);
fh = figure;
WinOnTop(fh)
 
while(true)
     
    while(mmfile.Data.header(1)<0) % wait for a new frame...
        if(mmfile.Data.header(1) == -2) % exit if Scanbox stopped
            return;
        end
    end
    mmfile.Format = {'int16' [1 16] 'header' ; ...
        'uint16' double([mmfile.Data.header(2) mmfile.Data.header(3)]) 'chA'};
    mchA = double(intmax('uint16')-mmfile.Data.chA);
        
    imshowpair(mchA,m)
    
    mmfile.Data.header(1) = -1; % signal Scanbox that frame has been consumed!
    drawnow limitrate;
     
end