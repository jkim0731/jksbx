
function redraw_meanimg(h)

if h.dat.procmap
%     I = h.dat.mimg_proc(:,:,h.dat.map);
    I = h.dat.mimg_proc(h.dat.ylim(1):h.dat.ylim(2), h.dat.xlim(1):h.dat.xlim(2), h.dat.map);
else    
%     I = h.dat.mimg(:,:,h.dat.map);
%     I = h.dat.mimg(h.dat.ylim(1):h.dat.ylim(2), h.dat.xlim(1):h.dat.xlim(2), h.dat.map);
    I = h.dat.mimg(max(h.dat.ylim(1),1):h.dat.ylim(2), max(h.dat.xlim(1),1):h.dat.xlim(2), h.dat.map);
end

mu = median(I(:));
sd1 = mean(abs(I(I<mu+1e-7) - mu));
sd2 = 1e-7 + mean(abs(I(I>mu-1e-7) - mu));

axes(h.axes2); imagesc(I, mu + 5*[-sd1 sd2]);
% xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off
colormap('gray')
axes(h.axes3); imagesc(I, mu + 5*[-sd1 sd2]);
% xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off
colormap('gray')

drawnow