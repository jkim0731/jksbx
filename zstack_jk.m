function [] = zstack_jk(fntag,fstart,fend,fnum)
% Converting z-stack to a single tif file
% Assume that each .sbx file has just one channel
% starting from fstart, until fend
% Assume the name is xxx_yyy_zzz (zzz includes fstart to fend (fnum)
% Reading at the frame # (fnum)

tifname = strcat(fntag,'_zstack.tif');
if exist(tifname)
    error('tif file exists')
    return
end
for i = fstart : fend
    if i < 10
        fn = strcat(fntag,'_00',num2str(i));
    elseif i < 100
        fn = strcat(fntag,'_0',num2str(i));
    elseif i < 1000
        fn = strcat(fntag,'_',num2str(i));
    else 
        error('too many files'); % can a filenum be over 999?
        return
    end
    temp = squeeze(sbxread(fn,fnum-1,1));
    imwrite(temp,tifname,'WriteMode', 'append');
end
return