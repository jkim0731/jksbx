function rf = sbxsparsenoise(fname)

% analyze sparse noise experiment

% example experimnet: fname = 'd:\gm8\gm8_000_002_nonrigid'


% edit Luis 1/10/17. Allow BOTH rigid and nonrigid .signals files as input..
    % -----
    if contains(fname,'rigid') % search for rigid in filename
        si = strfind(fname,'_'); 
        fnamelog = fname( 1:si(end)-1); % remove it
    else 
        fnamelog = fname;
    end
    % -----
    
log = sbxreadsparsenoiselog(fnamelog); % read log
load([fname, '.signals'],'-mat');
sig = medfilt1(sig,11);             % median filter
sig = zscore(sig);
dsig = diff(sig);    
p = prctile(dsig,65);
dsig = bsxfun(@minus,dsig,p);
dsig = dsig .* (dsig>0);
dsig = zscore(dsig);

ncell = size(dsig,2);
nstim = size(log,1);

r = zeros(max(log.xpos),max(log.ypos),13,2,ncell);

h = waitbar(0,'Processing...');
for(i=1:nstim)
    if(log.mean(i)==0)  % rf map for dark spots
        r(log.xpos(i),log.ypos(i),:,1,:) =  squeeze(r(log.xpos(i),log.ypos(i),:,1,:)) + ...
            dsig(log.sbxborn(i)-2:log.sbxborn(i)+10,:);
    else                % rf map for brigths spots
        r(log.xpos(i),log.ypos(i),:,2,:) =  squeeze(r(log.xpos(i),log.ypos(i),:,2,:)) + ...
            dsig(log.sbxborn(i)-2:log.sbxborn(i)+10,:);
    end
    waitbar(i/nstim,h);
end
delete(h);

rf = cell(1,ncell);

h = fspecial('gauss',250,40);
k = 0;

hh = waitbar(0,'Filtering...');

for(n = 1:ncell)
    dark   = imresize(filter2(h,squeeze(mean(r(:,:,6:8,1,n),3)),'same'),0.25);
    bright = imresize(filter2(h,squeeze(mean(r(:,:,6:8,2,n),3)),'same'),0.25);
    rf{n} = {dark' bright'};
    waitbar(n/ncell,hh);
end
delete(hh);

save([fname '.sparsenoise'],'rf','-v7.3');


