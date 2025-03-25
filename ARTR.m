
function [Lk1,Sk1,E,Z]=ARTR(X,U,V,S,Z,etamax,itmax,beta,lambda,mu,W,eps,alpha)
L = H_tprod(U,V);N = sum(W(:));E = zeros(size(L));
sizeX = size(X);alpha1 = alpha;alpha2 = alpha;
flag = 0;eta = log(0.01);
Z0 = L.*double(~W)+S.*double(~W)+E.*double(~W);
Z = Z0+ X.*W;
while flag == 0
    for k = 1:itmax
       Lk = L;Sk = S;

       D = Z-S-E;
       U = Update_U(U,V,D,beta,mu);
       V = Update_V(U,V,D,beta,mu);
       L = H_tprod(U,V);



       if lambda == 0
           S = zeros(sizeX);
       else
           S = my_prox_l1(Z-L-E,lambda./mu);
       end

       for iE = 1:10
       D_E = exp(3*eta).*sqrt(2/pi).*W.*exp(-0.5*exp(2*eta).*(E.^2)).*E + mu.*(L+S+E-Z);
       E = E - alpha2.*D_E;
       end

        E0 = W.*exp(-0.5.*exp(2*eta).*(E.^2));
        E1 = W.*exp(-0.5.*exp(2*eta).*(E.^2)).*E.^2;
        D_eta = 0.5*N*exp(eta)./sqrt(pi) - exp(eta)*sqrt(2/pi)*sum(E0(:)) + exp(3*eta)*sqrt(2/pi)*sum(E1(:));
        eta0 = eta - alpha1.*D_eta;
        if eta0 <etamax
            eta = eta0;
        else
            eta = etamax;
        end
        
       Z0 = L.*double(~W)+S.*double(~W)+E.*double(~W);
       Z = Z0+ X.*W;
       
        Lk1 = L;Sk1 = S;
        if max(Lk1(:)-Lk(:)) < eps && max(Sk1(:)-Sk(:)) < eps
            flag = 1;
        elseif k == itmax
            flag = 1;
        end
    end 
end
end
