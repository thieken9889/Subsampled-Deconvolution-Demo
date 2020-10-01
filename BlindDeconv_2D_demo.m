clear,clc;
close all

% add src folder with all of our functions to the path
addpath('src');

%generate a the signals w=Bh and x=Cm
L1 = 50;
L2 = L1; % must equal L1
K = 9;  % must be a perfect square
N = 10;

%import basis for linear operator C, this was generated in other script
[scores, coeffs] = generate_subspace(L1, L2, 1000, N, 0);
Phis = reshape(coeffs(:,1:N),L1,L2,N); %use loaded principle components as subspace

% create ground truth
h = randn(K,1); %random FIR filter coefficients
m = scores(1,1:N).'; %use first sample from generated data and first N principal components
X0 = h*m';

%precomupte the bases for subspaces B and C
B_hat = compute_B_hat(L1,L2,K);
C_hat = compute_C_hat(L1,L2,N,Phis);

% the two vectors we are convolving
w = B_op(L1,L2,K,h,1);
x = C_op(L1,L2,N,Phis,m,1);


%fully sampled in fourier domain
y_hat = (1/sqrt(L1*L2))*(fft2(w).*fft2(x));
y = sqrt(L1*L2)*ifft2(y_hat);

figure(1);
imshow(kron(y,ones(10)),[]);
title('Convolution')
pause(1);

%parameters for subsampling
M = round(0.75*L1*L2);
Omega = randperm(L1*L2,M)';

%operators for subsampling
ifftOp = linop_handles({[L1,L2],[L1,L2]}, @(x)ifft2(x)*sqrt(L1*L2), @(x)fft2(x)/sqrt(L1*L2) ,'C2C');
sampOp   = linop_compose( linop_subsample({[L1*L2,1],[M,1]},Omega), linop_reshape([L1,L2],[L1*L2,1]), ifftOp );

%inal lin_op
lin_op = linop_compose(sampOp, @(x,mode)A_op_2D(L1,L2,K,N,B_hat,C_hat,x,mode));

%subsampled data
z = sampOp(y_hat,1);

%adjust options for TFOCS
opts = tfocs_SCD;
opts.maxIts = 10000;
opts.tol = 1e-5;
opts.errFcn = @(f,z,x)error_measure(X0,f,z,x); % measure error between interations


%use TFOCS SCD solver to preform nuclear norm minimization
X1 = tfocs_SCD(prox_nuclear,{lin_op, -z},prox_l2(1e-3),0.01,[],[],opts);

%check the residual between the ground truth and reconstruction (MSE)
err_MSE = norm(X1-X0).^2 /numel(X0)

%check cosine of the angle between ground truth and reconstruction
corr_coeff = trace(X1'*X0)/(norm(X0)*norm(X1))

%% visualize our solution for h and m
% note that we will always be off by some unknown scaling factor

%svd to separate (almost) rank-1 matrix, remember X0 = h*m.'
[U, S, V] = svd(X1);
hnew = U(:,1); %left singular vector with largest singular value is h
mnew = V(:,1); %right singular vector with largest singular value is m

xnew = C_op(L1,L2,N,Phis,real(mnew),1); %reconstruct image with our solution
wnew = B_op(L1,L2,K,real(hnew),1); %reconstruct impulse response with our solution


%things are off by a scaling factor, so showing the normalized coeffcients
note = 'Coeffcients are Normalized';

% compare reconstructed filter coefficients to the ground truth
figure(1);
subplot(1,2,2)
stem((-1*real(hnew))/max(abs(hnew)));
title('Reconstructed Filter Coeffs (Normalized)')
subplot(1,2,1)
stem(h/max(abs(h)));
title('Orignial Filter Coeffs (Normalized)')
pause(1);

% calculate cosine of angle between filter coefficients
h_cor_coeff = (h'*hnew)/(norm(h)*norm(hnew))

% Compare the reconstructed PCA weights to the ground truth
figure(2);
subplot(1,2,2)
stem(-1*real(mnew)/max(abs(mnew)));
ylim([-0.15 0.15])
title('Reconstructed PCA weights (Normalized)')
subplot(1,2,1)
stem(m/max(abs(m)));
ylim([-0.15 0.15])
title('Orignial PCA weights (Normalized)')
pause(1);

% calculate the cosine of the angle between the PCA weights
m_cor_coeff = (m'*mnew)/(norm(m)*norm(mnew))

% Show the images created using the PCA weights
figure(3);
subplot(1,2,2)
imshow(kron(-1*xnew,ones(10)),[]);
title('Reconstructed Image')
subplot(1,2,1)
imshow(kron(x,ones(10)),[]);
title('Orignial Image')
pause(1);
