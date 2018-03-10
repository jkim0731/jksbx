function [m,T] = align_zstack(zstack)

% Aligns images in all frames in a zstack
% 
% m - mean image after the alignment (similar to maximum projection)
% T - optimal translation for each frame


if size(zstack,3) == 1
    
    m = squeeze(zstack(:,:,1));    
    T = [0 0];
    
elseif size(zstack,3) == 2
    
    A = squeeze(zstack(:,:,1));
    B = squeeze(zstack(:,:,2));
    Ap = adapthisteq(A(50:end-50,151:end-50));
    Bp = adapthisteq(B(50:end-50,151:end-50));
    
    [u, v] = jkfftalign(Ap,Bp);
    
    Ar = circshift(A,[u,v]);
    m = (Ar+B)/2;
    T = [[u v] ; [0 0]];
    
else
    
    zstack1 = zstack(:,:,1:floor(end/2));
    zstack2 = zstack(:,:,floor(end/2)+1 : end);
    
    [A,T0] = align_zstack(zstack1);
    [B,T1] = align_zstack(zstack2); 

    Ap = adapthisteq(A(50:end-50,151:end-50));
    Bp = adapthisteq(B(50:end-50,151:end-50));
    
    [u, v] = jkfftalign(Ap,Bp);    
     
    Ar = circshift(A,[u, v]);
    m = (Ar+B)/2;
    T = [(ones(size(T0,1),1)*[u v] + T0) ; T1];
    
end
