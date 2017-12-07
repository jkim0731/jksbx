% compute bidirectional phase offsets from poor line scanning timing
function BiDiPhase = jk_BiDiPhaseOffsets(data)

[Ly, Lx, ~, NT] = size(data);

% lines scanned one direction
yr1 = 101:2:floor(Ly/2)*2; % Because of optotune ringing
% lines scanned in other direction
yr2 = 101:2:floor(Ly/2)*2; % Because of optotune ringing

% compute phase correlation between lines in x-direction
eps0 = single(1e-6);
Nmax = min(50, NT);
d1 = fft(data(yr1,101:end,:,1:Nmax),[],2); % 101:end because of bidirecitonal biproduct
d2 = conj(fft(data(yr2,101:end,:,1:Nmax),[],2)); % 101:end because of bidirecitonal biproduct
d1 = d1./(abs(d1) + eps0);
d2 = d2./(abs(d2) + eps0);

cc = ifft(d1 .* d2,[],2);
cc = fftshift(cc, 2);
cc = mean(mean(mean(cc,1),3),4);

[~, ix] = max(cc);
ix       = ix - (floor(Lx/2) + 1);

BiDiPhase = -1 * ix;



    