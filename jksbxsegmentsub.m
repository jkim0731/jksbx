
function plist = jksbxsegmentsub(img,mag,roinum)

% segment a subimage

switch mag
    case 1
%         W = 50;
        W = 10; % 03/03/2016 JK
    case 2
        W = 100;
    case 4
        W = 150;
end
W = W + 1;

if ~isa(img,'double'), img = double(img); end
mimg = mean(img,3);

mimg = mat2gray(mimg);
subplot(2,3,1), imagesc(mimg), axis image, axis off, drawnow

radius = 5;

[col_im, row_im] = meshgrid(1:size(img,1), 1:size(img,2));
ref_rg = (col_im - W).^2 + (row_im - W).^2 <= radius.^2;
ref_rg = ref_rg';
% subplot(2,3,2), imagesc(ref_rg), axis image, axis off, drawnow

ref_rg_im = bsxfun(@times, img, ref_rg);
ref_ts = squeeze(mean(mean(ref_rg_im,1),2));
corr_im = zeros(size(img,1),size(img,2));
p_im = zeros(size(img,1),size(img,2));

[max_f, max_idx] = sort(ref_ts, 'descend');
max_f_std = max_f > 3*std(ref_ts(:));
max_idx = max_idx(1:min(10000,sum(max_f_std(:))));

for i = 1 : size(img,1)
    for j = 1 : size(img,2)
        temp_ts = squeeze(img(i,j,max_idx));
        [corr_im(i,j), p_im(i,j)] = corr(ref_ts(max_idx), temp_ts);
    end
end

subplot(2,3,2), imagesc(corr_im), axis image, axis off, drawnow, colorbar
subplot(2,3,3), imagesc(p_im, [0 10^(-14)]), axis image, axis off, drawnow, colorbar

sigcorr_im = corr_im .* ((10^(-14) - (double(p_im < 10^-14)).*p_im)/10^(-14));

corr_th = 0.1;
sigcorr_im = corr_im > corr_th;
subplot(2,3,4), imagesc(double(sigcorr_im)), axis image, axis off, drawnow



L = zeros(size(img,1),size(img,2));
bw = sigcorr_im;
bw = imfill(bw,4,'holes');
cc = regionprops(bw,'Area','Solidity','PixelIdxList','Eccentricity','PixelList');
plist = [];
for i = 1 : length(cc)
if(cc(i).Area<1000*mag^2 && cc(i).Area > 25*mag^2 && cc(i).Solidity>0.7 && cc(i).Eccentricity<0.9)
    L(cc(i).PixelIdxList)=1;
    plist = cc(i).PixelList;
end
end
subplot(2,3,5), imagesc(double(L)), axis image, axis off, drawnow
title(['cell ', num2str(roinum), ': ', num2str(length(max_idx)), 'time points'])


subplot(2,3,6), imhist(corr_im)

