function ops = jk_write_reg_to_tiff(fid, ops, iplane, isRED)
%%
Ly = ops.Ly;
Lx = ops.Lx;
bitspersamp = 16;

frewind(fid);

ix = 0;    
nframesleft = ops.Nframes;

while nframesleft>0
     ix = ix + 1;
    nfrtoread = min(nframesleft, 2000);
    data = fread(fid,  Ly*Lx*nfrtoread, 'int16');                
    nframesleft = nframesleft - nfrtoread;
    data = reshape(data, Ly, Lx, []);        

     % this is now done in a separate function
%         if ix==1
%             navg = min(size(data,3), ops.nimgbegend);
%             ops.mimg_beg(:,:,k) = mean(data(:,:,1:navg), 3);
%         end
%         if nframesleft<ops.nimgbegend
%             nconc = min(ops.nimgbegend - nframesleft, nfrtoread);
%             datend = cat(3, datend, data(:,:,(nfrtoread-nconc+1):nfrtoread));
%         end
%         if nframesleft<=0
%              ops.mimg_end(:,:,k) = mean(datend, 3);
%         end

    foldr = fullfile(ops.RegFileTiffLocation, ops.mouse_name, sprintf('%03d',ops.session), sprintf('Plane%d', iplane));
    if ~exist(foldr, 'dir')
        mkdir(foldr)
    end
    if isRED
        partname = sprintf('%s_%03d_2P_plane%d_%d_RED.tif', ops.mouse_name, ops.session, iplane, ix);

    else
        partname = sprintf('%s_%03d_2P_plane%d_%d.tif', ops.mouse_name, ops.session, iplane, ix);
    end
    fname = fullfile(foldr, partname);

    TiffWriter(int16(data),fname,bitspersamp);
end
