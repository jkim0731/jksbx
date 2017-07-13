function [m,T] = jkaligncell(mg)

% mg is a cell containing series of images (m group) to be aligned
% Aligns images in fname for all indices in idx
% 
% m - mean image after the alignment
% T - optimal translation for each frame

if(length(mg)==1)
    
    m = mg{1};
    T = [0 0];
    
elseif (length(mg)==2)
    
    A = mg{1};
    B = mg{2};
    Ap = A(:,150:end);
    Bp = B(:,150:end);
    
    [u v] = fftalign(Ap,Bp);
    
    Ar = circshift(A,[u,v]);
    m = (Ar+B)/2;
    T = [[u v] ; [0 0]];
    
else
    
    idx0 = mg(1:floor(end/2));
    idx1 = mg(floor(end/2)+1 : end);
    
    [A,T0] = jkaligncell(idx0);
    [B,T1] = jkaligncell(idx1);
   
    [u v] = fftalign(A,B);
     
    Ar = circshift(A,[u, v]);
    m = (Ar+B)/2;
    T = [(ones(size(T0,1),1)*[u v] + T0) ; T1];
    
end
