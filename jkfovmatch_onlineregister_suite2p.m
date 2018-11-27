
% pwd
close all

num_iter = 30;
mouse_num = '070';
layer = 1;
redchannel = false;

ref_fn_list1 = {'070_002'};
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
if ~exist(['regops_', ref_fn, '.mat'], 'file')
    error('Run suite2p first')
end
load(['regops_', ref_fn, '.mat'])

ref1 = mat2gray(ops1{1}.mimg1); % 1 is the top FOV
if redchannel
    ref2 = mat2gray(ops1{1}.REDmimg1);
end
figure, imshow(mat2gray(ref1)), axis image

ref1fov = zeros(ops1{1}.info.sz);
ref1fov(ops1{1}.useY, ops1{1}.useX) = ref1;
if redchannel
    ref2fov = zeros(ops1{1}.info.sz);
    ref2fov(ops1{1}.useY, ops1{1}.useX) = ref2;
end

mmfile = memmapfile('scanbox.mmap','Writable',true, ...
    'Format', { 'int16' [1 16] 'header' } , 'Repeat', 1);
figure,
if redchannel
    h1 = subplot(121);
    h2 = subplot(122);
end

while(true)
    while(mmfile.Data.header(1)<0) % wait for a new frame...
        if(mmfile.Data.header(1) == -2) % exit if Scanbox stopped
            return;
        end
    end
       
    mchA = []; 
    mchB = [];
    for i = 1 : num_iter
        if redchannel
            mmfile.Format = {'int16' [1 16] 'header' ; ...
                'uint16' double([mmfile.Data.header(2) mmfile.Data.header(3)]) 'chA';...
                'uint16' double([mmfile.Data.header(2) mmfile.Data.header(3)]) 'chB'};
        else
            mmfile.Format = {'int16' [1 16] 'header' ; ...
                'uint16' double([mmfile.Data.header(2) mmfile.Data.header(3)]) 'chA'};
        end
        if i == 1 
            mchA = double(intmax('uint16')-mmfile.Data.chA);
            mfinalA = mchA;
            if redchannel
                mchB = double(intmax('uint16')-mmfile.Data.chB);
                mfinalB = mchB;
            end
        else
            prev_mchA = mfinalA;
            mchA = double(intmax('uint16')-mmfile.Data.chA);
            
            [u1, v1] = jkfftalign(mchA(50:end-50,50:end-50),prev_mchA(50:end-50,50:end-50));
            mchA = circshift(mchA,[u1,v1]);
            mfinalA = mfinalA + mchA/num_iter;
            
            if redchannel
                prev_mchB = mfinalB;
                mchB = double(intmax('uint16')-mmfile.Data.chB);

                [u1, v1] = jkfftalign(mchB(50:end-50,50:end-50),prev_mchB(50:end-50,50:end-50));
                mchB = circshift(mchB,[u1,v1]);
                mfinalB = mfinalB + mchB/num_iter;
            end
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
    
    if redchannel
        imshowpair(mfinalA(101:end-10,101:end-10),ref1fov(101:end-10,101:end-10),'Parent', h1), axis image,    
        imshowpair(mfinalB(101:end-10,101:end-10),ref2(101:end-10,101:end-10), 'Parent', h2), axis image,
    else
        imshowpair(mfinalA(101:end-10,101:end-10),ref1fov(101:end-10,101:end-10)), axis image,
    end

    drawnow limitrate;     
end