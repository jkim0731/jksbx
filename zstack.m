%%
base_d = 'D:\2p\JK\test\tdTomato_L4\';
for i = 36:3:39
    cd(sprintf('%s%03d',base_d,i));
    delete *.align *.tif
    jksbxaligndir({},'green')

%     flist = dir('*.align');
%     load(flist(1).name,'-mat');    
%     if exist('m','var')
%         im = zeros(size(m,1),size(m,2),length(flist));        
%         
%         for j = 1 : length(flist)
%             load(flist(j).name,'-mat')
%             im(:,:,j) = m;
%         end
%         [~,T] = jksbxalignx_im(im,1:length(flist));
%         tifname = sprintf('%03d_000.tif',i);
%         for j = 1 : length(flist)
%             im(:,:,j) = circshift(im(:,:,j),T(j,:));
%             imwrite(im(:,:,j),tifname,'WriteMode','append')
%         end
%     elseif exist('m1','var') && exist('m2','var')
%         im1 = zeros(size(m1,1),size(m1,2),length(flist));
%         im2 = zeros(size(m2,1),size(m2,2),length(flist));
%         for j = 1 : length(flist)
%             load(flist(j).name,'-mat')
%             im1(:,:,j) = m1;
%             im2(:,:,j) = m2;
%         end                
%         [~,T1] = jksbxalignx_im(im1,1:length(flist));
%         [~,T2] = jksbxalignx_im(im2,1:length(flist));
%         gtifname = sprintf('%03d_000g.tif',i);
%         rtifname = sprintf('%03d_000r.tif',i);
%         for j = 1 : length(flist)
%             im1(:,:,j) = circshift(im1(:,:,j),T1(j,:));
%             im2(:,:,j) = circshift(im2(:,:,j),T2(j,:));
%         end
%         
%         for j = 1 : length(flist)
%             im1(:,:,j) = im1(:,:,j)/max(max(max(im1)));
%             im2(:,:,j) = im2(:,:,j)/max(max(max(im2)));
%             imwrite(im1(:,:,j),gtifname,'WriteMode','append')
%             imwrite(im2(:,:,j),rtifname,'WriteMode','append')
%         end
%     end
end
