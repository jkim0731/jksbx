function maxidx = jkget_maxidx(fn)
load([fn,'.mat']);

if(info.scanmode==0)
    info.recordsPerBuffer = info.recordsPerBuffer*2;
end
switch info.channels
    case 1
        factor = 1;
    case 2
        factor = 2;
    case 3
        factor = 2;
end
d = dir([fn, '.sbx']);    
maxidx =  d.bytes/info.recordsPerBuffer/info.sz(2)*factor/4 - 1;