num_iter = 10;
mmfile = memmapfile('scanbox.mmap','Writable',true, ...
    'Format', { 'int16' [1 16] 'header' } , 'Repeat', 1);
fh = figure;
WinOnTop(fh);

while(true)
    while(mmfile.Data.header(1)<0) % wait for a new frame...
        if(mmfile.Data.header(1) == -2) % exit if Scanbox stopped
            return;
        end
    end
       
    mchA = [];
    for i = 1 : num_iter
        mmfile.Format = {'int16' [1 16] 'header' ; ...
            'uint16' double([mmfile.Data.header(2) mmfile.Data.header(3)]) 'chA'};
%         aaa = mmfile.Data;
        if i == 1 
            mchA = double(intmax('uint16')-mmfile.Data.chA)/num_iter;
        else
            mchA = mchA + double(intmax('uint16')-mmfile.Data.chA)/num_iter;    
        end
        mmfile.Data.header(1) = -1; % signal Scanbox that frame has been consumed!
    end    
    imagesc(mchA)
    drawnow limitrate;
     
end