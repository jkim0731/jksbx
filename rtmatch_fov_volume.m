% method 01: matching imaging plane with pre-imaged volume
% Real-time correlation analysis with the volume centered at target plane.
% Volume should be imaged at the very beginning of experiment, 0.5 um resolution with 50 um thickness
% Assume I start 'Focus' first.

mouse = 'JK049';
layer = 1; % or 2. 1 upper, 2 lower.

% Load volumetric imaging data.
data = load([mouse, '_target_volume']); % need to contain target_img(:,:,:,1:2), target_img_maxproj(:,:,1:2)

% Set # of frames you want to use
num_frames = 100;

% Real-time imaging registration between frames

mmfile = memmapfile('scanbox.mmap','Writable',true, ...
    'Format', { 'int16' [1 16] 'header' } , 'Repeat', 1);
fh1 = figure; h1 = axes('Parent',fh1);

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
        mfinal = mfinal + mchA/num_iter;
    end

    mmfile.Data.header(1) = -1; % signal Scanbox that frame has been consumed!
end  

figure, imshow(mat2gray(mfinal))

% Register the imaged plane with maximum projection of the volumetric data
% (for rotation)



% Calculate correlation between the imaged plane with every slice of the
% volumetric data. 



% First compare with a threshold to see if the plane is within the volume.
% And then according to the correlation values, assign which slice is the best match.
% Move the objective to the desired plane.




% Repeat imaging until the maximum correlation value is converged (need to
% have a threshold to testing convergence).
