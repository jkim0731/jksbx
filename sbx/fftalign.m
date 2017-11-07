
function [u,v] = fftalign(A,B)

N = min(size(A))-80;    % leave out margin

yidx = round(size(A,1)/2-N/2) + 1 : round(size(A,1)/2+ N/2);
xidx = round(size(A,2)/2-N/2) + 1 : round(size(A,2)/2+ N/2);

A = A(yidx,xidx);
B = B(yidx,xidx);

C = fftshift(real(ifft2(fft2(A).*fft2(rot90(B,2)))));
[~,i] = max(C(:));
[ii jj] = ind2sub(size(C),i);

u = round(N/2-ii);
v = round(N/2-jj);