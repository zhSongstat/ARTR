%% Exp_2_Gaussian_noise
clc,clear
rng(42)
%%
noise_ratio = 0.2;
%%
L0 = double(imread('1.jpg'));
O = L0;maxP = max(abs(O(:)));
[n1,n2,n3] = size(L0);totalNum=n1*n2*n3;
Allindex = 1:totalNum;
%% Generate the noise term by Gaussian distribution N(5,1)
meanmu = 0; sigma2 = 30;
E0 = zeros(n1,n2,n3);
index_E0 = randperm(length(Allindex),floor(noise_ratio*length(Allindex)));
noise_Omega = Allindex(index_E0);
E0(noise_Omega) = normrnd(meanmu,sigma2,1,length(noise_Omega));
%% Generate data following P_Omega(X) = P_Omega(L+S+E)
X = L0+E0;
corruption_rate = 0.5; 
sample_rate = 1 - corruption_rate;
Omegac =  randperm(totalNum,floor(totalNum*corruption_rate));
X(Omegac(:)) = 0;
Omega = setdiff(Allindex,Omegac);
%% 
itmax = 500;
W = ones(size(L0));
W(Omegac(:)) = 0;
eps = 1e-6;
etamax = log(50);
eta = log(0.01);
S = zeros(size(X));
Z = X;
N = sum(W(:));
%%
opt = 'skinny';
transform.L = 'fft';
[U,Sig,V] = H_tsvd(X,transform,opt);
K = ceil(0.35*min(n1,n2));
U = H_tprod(U(:,1:K,:),Sig(1:K,1:K,:));
V = H_tran(V);V = V(1:K,:,:);
%%
lambda = 1/sqrt(sample_rate*max(n1,n2)*n3);
mu1 = 0.05*sqrt(lambda);
%%
Lh = N*exp(etamax)./(2*sqrt(pi))+4*N*exp(3*etamax);
Lfk = N*exp(3*etamax)*sqrt(2./pi)+mu1;
alpha = 0.99 * 2./(max(Lh,Lfk));
%%
beta = 0.8;
[Lk1,Sk1,E,Z] = ARTR(X,U,V,S,Z,etamax,itmax,beta,lambda,mu1,W,eps,alpha);
X_hat = max(Lk1,0);
X_hat = min(X_hat,maxP);
psnr_ours = PSNR(O,X_hat,maxP)
SSIM_FRTR_L = ssim(X_hat,O)