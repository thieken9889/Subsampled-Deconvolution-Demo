function out = C_op(L1,L2,N,Phis,in,mode)
%C_OP Function to perform the actions of the linear operator C.
%
%   C is a linear operation that goes from R^N to C^(L1 by L2). This
%   function impliments that linear operater in the format that TFOCS needs
%   for its optimization algorithms. This one is typically the subspace of
%   an image set, not the filter.
%
%   INPUTS:
%       L1      -       First dimesion of "in"
%       L2      -       Second dimension of "in"
%       N       -       Number of basis vectors in the subspace
%       Phis    -       The basis vectors of the subspace
%       in      -       L1 by L2 input
%       mode    -       Indicates whether to perform the forward or adjoint
%                       operation (or just return input and output size)
%
%   OUTPUTS:
%       out     -       if mode is 0 then returns the input and output
%                       dimensions, if mode is 1 then performs the forward 
%                       opertion, if mode is 2 performs the adjoint operation

    
    %subspace has matrices as its basis vectors

    %will be used to multiply along third direction of basis vectors
    mtx_vec_multplr = dsp.ArrayVectorMultiplier('Dimension',3);

    switch mode
        case 0
            %input and output size for forward operation
            out = {[N,1],[L1,L2]};
            
        case 1
            %case for forward operation
            %scale each basis vector by input then sum
            out = sum(mtx_vec_multplr(Phis,in),3);
            
        case 2
            %adjoint operation
            out = zeros(N,1);
            
            for i = 1:N
                out(i) = trace(Phis(:,:,i)'*in);
            end

    end
    
end

