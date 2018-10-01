fn = '054_008_000';
global info
sbxread(fn,0,1);
laser = zeros(info.max_idx,1);
framePerChunk = 10000;
numChunk = ceil(info.max_idx/framePerChunk);
%%
for i = 1 : numChunk - 1    
    img = sbxread(fn,(i-1)*framePerChunk,framePerChunk);
    img = squeeze(img(1,:,:,:));
    laser((i-1)*framePerChunk+1:i*framePerChunk) = squeeze(mean(mean(img)));
end
%%
img = sbxread(fn,(numChunk-1)*framePerChunk,info.max_idx-(numChunk-1)*framePerChunk);
img = squeeze(img(1,:,:,:));
laser((numChunk-1)*framePerChunk+1:end) = squeeze(mean(mean(img)));
%%
ttlOn = info.frame(info.event_id == 3);

ttlOff = info.frame(info.event_id == 2);
ttlOn(249:end) = ttlOn(249:end) + 2^16;
ttlOff(250:end) = ttlOff(250:end) + 2^16;
%%
figure, plot(laser), hold on, plot(ttlOn, laser(ttlOn),'r.'), hold on, 
plot(ttlOff(2:end-1), laser(ttlOff(2:end-1)), 'k.')

