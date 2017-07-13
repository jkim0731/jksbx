
% pwd

num_iter = 8;
mouse_num = '018';

% ref_fn_list = {'648_000_006';'650_000_001';'651_000_004';'652_000_004';'653_000_004'};
ref_fn_list = {'017_033_000';'018_032_000';'020_032_000'};
for i = 1 : length(ref_fn_list)
    if strfind(ref_fn_list{i},mouse_num) == 1
        ref_fn = ref_fn_list{i};
        break
    end
end

cd(['d:\2p\JK\',ref_fn(1:3)])
load([ref_fn,'.align'], '-mat')
if iscell(m)
    ref = mat2gray(m{end});
else
    ref = mat2gray(m);
end

mmfile = memmapfile('scanbox.mmap','Writable',true, ...
    'Format', { 'int32' [1 16] 'header' } , 'Repeat', 1);
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
        mmfile.Format = {'int32' [1 16] 'header' ; ...
            'uint16' double([mmfile.Data.header(2) mmfile.Data.header(3)]) 'chA'};
%         aaa = mmfile.Data;
        if i == 1 
            mchA = double(intmax('uint16')-mmfile.Data.chA)/num_iter;
        else
            mchA = mchA + double(intmax('uint16')-mmfile.Data.chA)/num_iter;    
        end
        mmfile.Data.header(1) = -1; % signal Scanbox that frame has been consumed!
    end    
    imshowpair(mchA,ref)
    drawnow limitrate;
     
end