function [m2,T2] = sbxalignx2iter(fname,idx)

% Aligns images in fname for all indices in idx
% 
% m - mean image after the alignment
% T - optimal translation for each frame


load([fname, '.align'], '-mat') % contains m and T

if(length(idx)==1)
    
    A = sbxread(fname,idx(1),1);
    A = squeeze(A(1,:,:));
        A = circshift(A,T(idx(1)+1,:));
    
    m2 = A;
    T2 = [0 0];
    
elseif (length(idx)==2)
    
    A = sbxread(fname,idx(1),1);
    B = sbxread(fname,idx(2),1);
    A = squeeze(A(1,:,:));
    B = squeeze(B(1,:,:));
        A = circshift(A,T(idx(1)+1,:));
        B = circshift(B,T(idx(2)+1,:));    

    [u v] = fftalign(A,B);
    
    Ar = circshift(A,[u,v]);
    m2 = (Ar+B)/2;
    T2 = [[u v] ; [0 0]];
    
else
    
    idx0 = idx(1:floor(end/2));
    idx1 = idx(floor(end/2)+1 : end);
        
    [A,T0] = sbxalignx(fname,idx0);
    [B,T1] = sbxalignx(fname,idx1);
    
    [u v] = fftalign(A,B);
     
    Ar = circshift(A,[u, v]);
    m2 = (Ar+B)/2;
    T2 = [(ones(size(T0,1),1)*[u v] + T0) ; T1];
    
end
