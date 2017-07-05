function sbxmovie_reg(fname,idx)
global info;

clf;
colormap gray;
writerObj = VideoWriter([fname '.m4v'],'MPEG-4');
writerObj.Quality = 96;
writerObj.FrameRate = 30;
open(writerObj);
load([fname,'.regnxc'],'-mat')

for(i=1:length(idx))
    z = squeeze(sbxread(fname,idx(i),1));
    z = circshift(z,T(idx(i)+1,:));
    
    imshow(z(:,70:end-70));
    text(20,20,sprintf('%3.1f sec frame #%05d',idx(i)/(info.resfreq/512*(2-info.scanmode)),idx(i)),'color',[1 1 1],'fontsize',16);
    drawnow;
    frame = getframe;
    writeVideo(writerObj,frame);
end

close(writerObj);