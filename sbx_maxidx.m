function maxidx = sbx_maxidx(fn)

a = sbxread(fn,0,1);
global info
maxidx = info.max_idx;
clear info
end

