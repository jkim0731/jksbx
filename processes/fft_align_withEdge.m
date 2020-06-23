function [u,v] = fft_align_withEdge(A,B,N)
    % Copied from Scanbox (Dario)
    N = 50;
    
    
    C = fftshift(real(ifft2(fft2(A).*fft2(rot90(B,2)))));
    [~,i] = max(C(:));
    [ii, jj] = ind2sub(size(C),i);

    u = -ii;
    v = -jj;

end
