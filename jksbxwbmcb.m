function jksbxwbmcb(src,callbackdata)

global segmenttool_h data nhood frames options th_corr bgimg mode_h
global zimg hline vline nm

try
    if(mode_h.Value == 1)
        p = gca;
        z = round(p.CurrentPoint);
        z = z(1,1:2);
        if(z(1)>0 && z(2)>0 && z(1)< size(bgimg.CData,2) && z(2)<size(bgimg.CData,1))
            cm = squeeze(sum(bsxfun(@times,data(z(2),z(1),:),data),3))/size(data,3);
            imgth = gather(cm>th_corr);
            bgimg.CData(:,:,2) = uint8(255*imgth);
            D = bwdistgeodesic(imgth,z(1,1),z(1,2));
            bw = imdilate(isfinite(D),strel('disk',1));
            bgimg.CData(:,:,3) = uint8(255*bw);
            
            % move axis...
            zimg.CData(:,:,1) = uint8(255*(nm + edge(bw)/3));
            set(zimg.Parent,'xlim',[z(1)-32 z(1)+32],'ylim',[z(2)-32 z(2)+32]);
            set(zimg.Parent,'CLimMode','manual');
            ztemp = 255*nm(sub2ind(size(nm),[z(1)-32 z(1)+32],[z(2)-32 z(2)+32]));
            set(zimg.Parent,'CLim',[min(ztemp), max(ztemp)]);
            set(hline,'xdata',[z(1)-32 z(1)+32],'ydata', [z(2) z(2)]);
            set(vline,'ydata',[z(2)-32 z(2)+32],'xdata', [z(1) z(1)]);            
            
            drawnow;
        end
    end
catch
    return
end