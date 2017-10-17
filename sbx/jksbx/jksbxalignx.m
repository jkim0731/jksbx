function [m,T] = jksbxalignx(fname,idx)

% Aligns images in fname for all indices in idx
% 
% m - mean image after the alignment
% T - optimal translation for each frame

if(length(idx)==1)
    
    A = sbxread(fname,idx(1),1);
    A = squeeze(A(1,:,:));
    m = A;
    T = [0 0];
    
elseif (length(idx)==2)
    
    A = sbxread(fname,idx(1),1);
    B = sbxread(fname,idx(2),1);
    A = squeeze(A(1,:,:));
    B = squeeze(B(1,:,:));
    Ap = A(151:end-50,151:end-50);
    Bp = B(151:end-50,151:end-50);
    
    [u v] = jkfftalign(Ap,Bp);
    
    Ar = circshift(A,[u,v]);
    m = (Ar+B)/2;
    T = [[u v] ; [0 0]];
    
else
    
    idx0 = idx(1:floor(end/2));
    idx1 = idx(floor(end/2)+1 : end);
    try
        [A,T0] = jksbxalignx(fname,idx0);
        [B,T1] = jksbxalignx(fname,idx1);
    catch
        disp(['fname = ', fname])
        disp(['idx0 = ', num2str(idx0)])
        [A,T0] = jksbxalignx(fname,idx0);
        [B,T1] = jksbxalignx(fname,idx1);
    end
   
    [u v] = jkfftalign(A,B);
     
    Ar = circshift(A,[u, v]);
    m = (Ar+B)/2;
    T = [(ones(size(T0,1),1)*[u v] + T0) ; T1];
    
end
