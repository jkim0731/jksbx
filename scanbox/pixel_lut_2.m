
function S = pixel_lut_2

global sbconfig;

ncol = 796;   
% nsamp =  round(sbconfig.lasfreq / sbconfig.resfreq);
nsamp =  floor(sbconfig.lasfreq / sbconfig.resfreq); % JK
M = ncol+2;

n = acos(linspace(1,-1,M))*nsamp/(2*pi);
n = n(2:end-1);

S = floor(n)-1;
S = [S; S+1; S+2 ; S+3];
S = reshape(S,1,[]);


