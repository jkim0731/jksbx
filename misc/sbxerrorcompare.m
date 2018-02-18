frame = 139083;

fn = '030_011_000';

dir1 = 'E:\JKbackup_20180113\TPM\030';
dir2 = 'D:\TPM\JK\030';

cd(dir1)
global info
fig1 = squeeze(sbxread(fn,frame,1));
cd(dir2)
fig2 = squeeze(sbxread(fn,frame,1));

figure, 
subplot(121), imagesc(fig1), axis image, axis off
