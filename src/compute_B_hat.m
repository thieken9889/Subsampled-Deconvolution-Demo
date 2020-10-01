function B_hat = compute_B_hat(L1,L2,K)
%COMPUTE_B_HAT computes the fourier transform of the B subspace and vectorizes.
%   We know that the basis vectors for B are just the matrices with one
%   nonzero entry. So we compute the fourier transform of these and
%   vectorize them into columns of a new matrix B_hat.
%   INPUTS:
%       L1      -       First dimesion of "in"
%       L2      -       Second dimension of "in"
%       K       -       Number of basis vectors in the subspace
%
%   OUTPUTS:
%       B_hat   -       (L1xL2)-by-K matrix. Each column is the vectorized
%                       2D Fourier transform of one L1-by-L2 basis vector
    
    %assume the 2D filter is square and causal
    %i.e. only the top left corner is filled
    
    if(floor(sqrt(K))~=sqrt(K))
        disp('K does not correspond to square filter');
        B_hat = [];
        return
    end
    
    %square filter is sqrt(K) by sqrt(K)
    filter_size = sqrt(K);
    
    B_hat = zeros(L1*L2,K);
    
    for i = 1:K
        %single point in FIR filter
        temp = zeros(L1,L2);
        temp2 = zeros(K,1);
        temp2(i) = 1;
        temp(1:filter_size,1:filter_size) = reshape(temp2,filter_size,filter_size);
        
        B_hat(:,i) = reshape((1/sqrt(L1*L2))*fft2(temp), L1*L2, 1);
    end
    
end

