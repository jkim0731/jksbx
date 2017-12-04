function [m,T] = jksbxalignx(fname,idx,varargin)

% Aligns images in fname for all indices in idx
% 
% m - mean image after the alignment
% T - optimal translation for each frame

if nargin < 3
    ch = 1; % green (pmt0) when 2 channel imaging. 
else
    ch = varargin{1}; % 1 for green (pmt0), 2 for red (pmt1) when 2 channel imaging
end

if(length(idx)==1)
    
    A = sbxread(fname,idx(1),1);
    A = squeeze(A(ch,:,:));
    m = A;
    T = [0 0];
    
elseif (length(idx)==2)
    
    A = sbxread(fname,idx(1),1);
    B = sbxread(fname,idx(2),1);
    A = squeeze(A(ch,:,:));
    B = squeeze(B(ch,:,:));
    Ap = A(151:end-50,151:end-50);
    Bp = B(151:end-50,151:end-50);
    
    [u, v] = jkfftalign(Ap,Bp);
    
    Ar = circshift(A,[u,v]);
    m = (Ar+B)/2;
    T = [[u v] ; [0 0]];
    
else
    
    idx0 = idx(1:floor(end/2));
    idx1 = idx(floor(end/2)+1 : end);
    try
        [A,T0] = jksbxalignx(fname,idx0,ch);
        [B,T1] = jksbxalignx(fname,idx1,ch);
    catch
        disp(['fname = ', fname])
        disp(['idx0 = ', num2str(idx0)])
        [A,T0] = jksbxalignx(fname,idx0);
        [B,T1] = jksbxalignx(fname,idx1);
    end
   
    [u, v] = jkfftalign(A,B);
     
    Ar = circshift(A,[u, v]);
    m = (Ar+B)/2;
    T = [(ones(size(T0,1),1)*[u v] + T0) ; T1];
    
end
