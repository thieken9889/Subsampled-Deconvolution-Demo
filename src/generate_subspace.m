function [scores, coeffs] = generate_subspace(L1, L2, n_samp, n_comps, showplots)
    
    %check that L1 = L2, which is required for the phantom() function
    if L1 ~= L2
        error('L1 must equal L2 to generate this subspace')
    end
    
    %get parameters for generic shepp-logan then we can alter them slightly
    [~,E_shepp] = phantom('Modified Shepp-Logan',L1);

    %generate set of phantoms
    data = zeros(L1,L2,n_samp);
    for i = 1:n_samp
        %parameters for phantom
        E = E_shepp;
        E(3:end,4:5) = E(3:end,4:5)+(0.01)*randn(8,2); %shift some of the elipses around randomly
        data(:,:,i) = phantom(E,L1) +(0.01)*randn(L1,L1); %noise is 1% relative to total intensity

    end
    
    % PCA - principal component analysis 
    data_vec = reshape(data,L1*L2,n_samp);
    data_vec = data_vec.';

    [coeffs,scores,~] = pca(data_vec,'Centered',false);

    data_vec_recon = scores(:,1:n_comps+1)*coeffs(:,1:n_comps+1).'; %only use n_comps # of components to reconstruct
    data_recon = reshape(data_vec_recon.',L1,L2,n_samp);
    
    if showplots
        
        figure(1)
        for i = 1:10
            subplot(2,5,i)
            imshow(kron(data(:,:,i),ones(10,10)),[])
            title(num2str(i))
        end
        pause(0.5)

        figure(2)
        for i = 1:10
            subplot(2,5,i)
            imshow(kron(data_recon(:,:,i),ones(10)),[]);
            title(num2str(i))
        end
    end
    
end

