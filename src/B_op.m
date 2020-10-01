function out = B_op(L1,L2,K,in,mode)
%B_OP Function to perform the actions of the linear operator B.
%
%   B is a linear operator that goes from R^K to C^(L1 by L2). This
%   function impliments that linear operater in the format that TFOCS needs
%   for its optimization algorithms. This subspace is the subspace of 2-D
%   FIR filters, so we already have the basis vectors.
%
%   INPUTS:
%       L1      -       First dimesion of "in"
%       L2      -       Second dimension of "in"
%       K       -       Number of basis vectors in the subspace
%       in      -       L1 by L2 input
%       mode    -       Indicates whether to perform the forward or adjoint
%                       operation (or just return input and output size)
%
%   OUTPUTS:
%       out     -       if mode is 0 then returns the input and output
%                       dimensions, if mode is 1 then performs the forward 
%                       opertion, if mode is 2 performs the adjoint operation
    
    if(floor(sqrt(K))~=sqrt(K))
        disp('K does not correspond to square filter');
        out = [];
        return
    end

    %assume the 2D filter is square and causal
    %i.e. only the top left corner is filled
    
    %square filter is sqrt(K) by sqrt(K)
    filter_size = sqrt(K);
    
    switch mode
        case 0
            %input and output size for forward operation
            out = {[K,1],[L1,L2]};
        case 1
            %case for forward operation (synthesis)
            out = zeros(L1,L2);
            
            %faster way
            out(1:filter_size,1:filter_size) = reshape(in,filter_size,filter_size);
            
        case 2
            %adjoint operation (analysis)
            
            out = zeros(K,1);
            
            out = reshape(in(1:filter_size,1:filter_size),K,1);
end

