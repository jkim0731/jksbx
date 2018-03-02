
% pwd
close all

num_iter = 100;
mouse_num = '049';
layer = 1;

ref_fn_list1 = {'049_998_001'};
ref_fn_list2 = {'036_5554_001','037_5554_001','038_5554_001','039_5554_001','041_5554_001'};
if layer == 1
    for i = 1 : length(ref_fn_list1)
        if strfind(ref_fn_list1{i},mouse_num) == 1            
            ref_fn = ref_fn_list1{i};
            break            
        end
    end
elseif layer == 2    
    for i = 1 : length(ref_fn_list2)
        if strfind(ref_fn_list2{i},mouse_num) == 1
            ref_fn = ref_fn_list2{i};
            break
        end
    end
else
    error('No filename list')
end
cd(['d:\2p\JK\',ref_fn(1:3)])
load([ref_fn,'.align'], '-mat')

ref1 = mat2gray(m1{end}); % all m's are cell from 2017/11/29. Every .align file now has m1 for green and m2 for red. empty if not taken. 

% figure(1), imagesc(ref1(101:end-10,101:end-10)), axis image, 

mmfile = memmapfile('scanbox.mmap','Writable',true, ...
    'Format', { 'int16' [1 16] 'header' } , 'Repeat', 1);
fh1 = figure; h1 = axes('Parent',fh1);
WinOnTop(fh1);

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
        if i == 1 
            mchA = double(intmax('uint16')-mmfile.Data.chA);
            mfinal = mchA;
        else
            prev_mchA = mfinal;
            mchA = double(intmax('uint16')-mmfile.Data.chA);
            
            [u1, v1] = jkfftalign(mchA(50:end-50,50:end-50),prev_mchA(50:end-50,50:end-50));
            mchA = circshift(mchA,[u1,v1]);
%             Ar = circshift(mchA,[u1,v1]);
%             mchA = (Ar+prev_mchA)/2;            
            mfinal = mfinal + mchA/num_iter;
        end
        
        mmfile.Data.header(1) = -1; % signal Scanbox that frame has been consumed!
    end  
    
%     set(0,'CurrentFigure',fh1), imshowpair(mchA(101:end-10,101:end-10),ref1(101:end-10,101:end-10)), axis image,
%     set(0,'CurrentFigure',fh1), imshowpair(mchB(101:end-10,101:end-10),ref2(101:end-10,101:end-10)), axis image,
    imshowpair(mchA(101:end-10,101:end-10),ref1(101:end-10,101:end-10),'Parent',h1), axis image,
%     imshowpair(mchB(101:end-10,101:end-10),ref2(101:end-10,101:end-10),'Parent',h2), axis image,
    drawnow limitrate;     
end