function C_hat = compute_C_hat(L1,L2,N,Phis)
%COMPUTE_C_HAT computes the fourier transform of the C subspace and vectorizes.
%   We have basis vectors for C, which are matrices with. So we compute the
%   2D fourier transform of these and vectorize them into columns of a new 
%   matrix C_hat.
%   INPUTS:
%       L1      -       First dimesion of "in"
%       L2      -       Second dimension of "in"
%       N       -       Number of basis vectors in the subspace
%
%   OUTPUTS:
%       C_hat   -       (L1xL2)-by-K matrix. Each column is the vectorized
%                       2D Fourier transform of one L1-by-L2 basis vector
    
    C_hat = zeros(L1*L2,N);
    
    %pick basis vectors as gaussian matrices
    
    for i = 1:N
        C_hat(:,i) = reshape((1/sqrt(L1*L2))*fft2(Phis(:,:,i)),L1*L2,1);
    end

end

