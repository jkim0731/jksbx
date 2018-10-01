
% pwd
close all

num_iter = 15;
mouse_num = '053';
layer = 1;

ref_fn_list = {'052_S01_000','053_S01_000','054_S01_000','056_S01_000'};
for i = 1 : length(ref_fn_list1)
    if strfind(ref_fn_list{i},mouse_num) == 1            
        ref_fn = ref_fn_list{i};
        break
    end
end
cd(['d:\2p\JK\',ref_fn(1:3)])
if ~exist([ref_fn,'.align'], 'file')
    jksbxaligndir(ref_fn)
end
load([ref_fn,'.align'], '-mat')

if layer == 1
    ref1 = mat2gray(m1{1}); 
    ref2 = mat2gray(m2{1});
else
    ref1 = mat2gray(m1{5});
    ref2 = mat2gray(m2{5});
end
% figure(1), imagesc(ref1(101:end-10,101:end-10)), axis image, 
figure, imshow(mat2gray(ref1)), axis image

mmfile = memmapfile('scanbox.mmap','Writable',true, ...
    'Format', { 'int16' [1 16] 'header' } , 'Repeat', 1);
fh1 = figure; 

while(true)
    while(mmfile.Data.header(1)<0) % wait for a new frame...
        if(mmfile.Data.header(1) == -2) % exit if Scanbox stopped
            return;
        end
    end
       
    mchA = []; 
    mchB = [];
    for i = 1 : num_iter 
        mmfile.Format = {'int16' [1 16] 'header' ; ...
            'uint16' double([mmfile.Data.header(2) mmfile.Data.header(3)]) 'chA';...
            'uint16' double([mmfile.Data.header(2) mmfile.Data.header(3)]) 'chB'};
        if i == 1 
            mchA = double(intmax('uint16')-mmfile.Data.chA);
            mchB = double(intmax('uint16')-mmfile.Data.chB);
            mfinalA = mchA;
            mfinalB = mchB;
        else
            prev_mchA = mfinalA;
            mchA = double(intmax('uint16')-mmfile.Data.chA);
            
            [u1, v1] = jkfftalign(mchA(50:end-50,50:end-50),prev_mchA(50:end-50,50:end-50));
            mchA = circshift(mchA,[u1,v1]);
            mfinalA = mfinalA + mchA/num_iter;
            
            prev_mchB = mfinalB;
            mchB = double(intmax('uint16')-mmfile.Data.chB);
            
            [u1, v1] = jkfftalign(mchB(50:end-50,50:end-50),prev_mchB(50:end-50,50:end-50));
            mchB = circshift(mchB,[u1,v1]);
            mfinalB = mfinalB + mchB/num_iter;

        end
        mmfile.Data.header(1) = -1; % signal Scanbox that frame has been consumed!
    end  
    
    subplot(121),imshowpair(mfinalA(101:end-10,101:end-10),ref1(101:end-10,101:end-10)), axis image,
    subplot(122),imshowpair(mfinalB(101:end-10,101:end-10),ref2(101:end-10,101:end-10)), axis image,
    drawnow limitrate;     
end