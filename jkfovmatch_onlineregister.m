
% pwd
close all

num_iter = 15;
mouse_num = '070';
layer = 1;

ref_fn_list1 = {'070_002_000'};
ref_fn_list2 = {''};
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
if ~exist([ref_fn,'.align'], 'file')
    jksbxaligndir({ref_fn})
end
load([ref_fn,'.align'], '-mat')

ref1 = circshift(mat2gray(m1{1}),[0,0]); % all m's are cell from 2017/11/29. Every .align file now has m1 for green and m2 for red. empty if not taken. 
% ref2 = mat2gray(m2{end});
% figure(1), imagesc(ref1(101:end-10,101:end-10)), axis image, 
figure, imshow(mat2gray(ref1)), axis image

mmfile = memmapfile('scanbox.mmap','Writable',true, ...
    'Format', { 'int16' [1 16] 'header' } , 'Repeat', 1);
figure,
% h1 = subplot(121);
% h2 = subplot(122);

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
        
%         mmfile.Format = {'int16' [1 16] 'header' ; ...
%             'uint16' double([mmfile.Data.header(2) mmfile.Data.header(3)]) 'chB'};
%         if i == 1 
%             mchB = double(intmax('uint16')-mmfile.Data.chB);
%             mfinalB = mchB;
%         else
%             prev_mchB = mfinal;
%             mchB = double(intmax('uint16')-mmfile.Data.chB);
%             
%             [u1, v1] = jkfftalign(mchB(50:end-50,50:end-50),prev_mchB(50:end-50,50:end-50));
%             mchB = circshift(mchB,[u1,v1]);
% %             Ar = circshift(mchA,[u1,v1]);
% %             mchA = (Ar+prev_mchA)/2;            
%             mfinalB = mfinalB + mchB/num_iter;
%         end
% 
%         
        mmfile.Data.header(1) = -1; % signal Scanbox that frame has been consumed!
    end  
    
%     imshowpair(mfinalA(101:end-10,101:end-10),ref1(101:end-10,101:end-10),'Parent', h1), axis image,
%     imshowpair(mfinalB(101:end-10,101:end-10),ref2(101:end-10,101:end-10), 'Parent', h2), axis image,
    imshowpair(mfinalA(101:end-10,101:end-10),ref1(101:end-10,101:end-10)), axis image,

    drawnow limitrate;     
end