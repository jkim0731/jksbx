function movie_jk(fname,idx)

global info;

clf;
colormap gray;
writerObj = VideoWriter([fname '.m4v'],'MPEG-4');
writerObj.Quality = 96;
writerObj.FrameRate = 30;
open(writerObj);

fnalign = strcat(fname,'.align'); 
if(exist(fnalign))
    sprintf('%f already aligned. Using that information.', fname);
    load(fnalign,'-mat','T');
end

for i = idx
    z = sbxread(fname,i,1);
    z = squeeze(z);
    if(exist(fnalign))
        if i ~= 0
            z = circshift(z,T(i,:));
        end
    end
    
    imshow(z(:,70:end-70));
    text(20,20,sprintf('%3.1f sec frame #%05d',idx(i)/(info.resfreq/512*(2-info.scanmode)),idx(i)),'color',[1 1 1],'fontsize',16);
    drawnow;
    frame = getframe;
    writeVideo(writerObj,frame);
    
end

close(writerObj);