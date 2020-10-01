function [out] = A_op_2D(L1,L2,K,N,B_hat,C_hat,in,mode)
%A_OP Function to perform the actions of the linear operator B.
%
%   A is a linear operator that goes from R^(K-by-N) to C^(L1-by-L2). This
%   function impliments that linear operater in the format that TFOCS needs
%   for its optimization algorithms. This operator corresponds to our
%   obeservation of the convolution of 2 signals
%
%   INPUTS:
%       L1      -       First dimesion of "in"
%       L2      -       Second dimension of "in"
%       K       -       Number of basis vectors in the B subspace
%       N       -       Number of basis vectors in the C subspace
%       in      -       K-by-N input (forward operation), or L1-by-L2 input
%                       (adjoint operation)
%       mode    -       Indicates whether to perform the forward or adjoint
%                       operation (or just return input and output size)
%
%   OUTPUTS:
%       out     -       if mode is 0 then returns the input and output
%                       dimensions, if mode is 1 then performs the forward 
%                       opertion, if mode is 2 performs the adjoint operation

    switch mode
        case 0
            %return the input size and output size for forward operation
            out = {[K,N], [L1,L2]};
        case 1
            %forward operation - need B and C which are the subspaces that
            %the convolved vectors come from.
            %This generates the observations from the outer product that we
            %are trying to find
            out = zeros(L1*L2,1);
            
            out = sum((sqrt(L1*L2)*B_hat*in).*C_hat,2);
            
            out = reshape(out,L1,L2);
            
        case 2
            % adjoint operator
            out = zeros(K,N);

            
            out = B_hat'*diag(reshape(in,L1*L2,1))*conj(C_hat)*sqrt(L1*L2);
    end
end

