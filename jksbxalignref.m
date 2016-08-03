function jksbxalignref(fn, varargin)
% align everyframe that is aligned (.align exists) again to a reference
% image (in reffn, i.e., ref301.bmp)

if nargin > 1
    reffn = varargin{1};
else
    try
        rx_reffn = '?*ref\d*.(tif|tiff|bmp|jpg|jpeg|png)';
        d = dir; % assume that reffn has a form of refxxx.yyy (xxx being the number, yyy being the image|file format - tif, tiff, bmp, jpg, jpge, png.)
        for i = 1 : length(d)
            if ~isempty(regexp(d(i).name,rx_reffn,'match'))
                reffn = d(i).name;
                break
            end
        end
    catch
        error('error in reference image name')
        return
    end
end
      
if exist([fn, '.alignref'])
    display(sprintf('%s has already been aligned to a reference file %s',fn, reffn))
    return
else
    if ~exist([fn, '.align'])
        disp(sprintf('Aligning %s...', fn))
        a = sbxread(fn,1,1);            % read one frame to read the header of the image sequence
        global info                % this contains the information about the structure of the image
        tic
        [m,T] = sbxalignx(fn,0:info.max_idx-1);   %
        save([fn '.align'],'m','T');
        display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx,round(toc/60)));
    else
        disp(sprintf('%s is already aligned', fn))
        load([fn,'.align'], '-mat');
    end

    im_ref = imread(reffn);
    a = squeeze(jksbxread(fn,0,1)); 
    global info
    rep = squeeze(jksbxread(fn,info.max_idx-1,1)); % representative image is the last one that did not move. sbxalignx makes every frame aligned compared to the last frame
    [u, v] = fftalign(rep,im_ref);
    uv = repmat([u, v], size(T,1),1);
    T = T + uv;
    m = circshift(m,[u, v]);
    save([fn '.alignref'],'m','T');
    return
end
end




%% test code
% load('301_004_000.align','-mat','m')
% m_align = m;
% load('301_004_000.alignref','-mat','m')
% m_ref = m;
% ref_im = imread('ref301.bmp');
% figure, imshowpair(ref_im,m_align)
% figure, imshowpair(ref_im,m_ref)