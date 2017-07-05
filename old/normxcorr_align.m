function [u,v] = normxcorr_align(A,B)

N = min(size(A))-20;    % leave margin

yidx = round(size(A,1)/2)-N/2 + 1 : round(size(A,1)/2)+ N/2;
xidx = round(size(A,2)/2)-N/2 + 1 : round(size(A,2)/2)+ N/2;

A = A(yidx,xidx);
B = B(yidx,xidx);


C = normxcorr2(A,B);
[~,i] = max(C(:));
% [~,i] = max(abs(C(:))); % why absolute?
[ii jj] = ind2sub(size(C),i(1));

u = N/2-ii;
v = N/2-jj;
