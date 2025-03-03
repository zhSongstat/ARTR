function U = Update_U(U,V,W,beta,mu)
%%
Ndim = length(size(U));
Nway= size(U);
Iu = eye(Nway(2));
L = ones(1,Ndim);
for x = 3:Ndim
    U = fft(U,[],x);
    V = fft(V,[],x);
    W = fft(W,[],x);
    L(x) = L(x-1) * Nway(x);
end
%%
U(:,:,1) = mu.*W(:,:,1)*(V(:,:,1)')*pinv(beta.*Iu+mu.*V(:,:,1)*V(:,:,1)');
%V(:,:,1) = pinv(mu.*U(:,:,1).'*U(:,:,1)+beta.*I)*U(:,:,1).'*(mu.*W(:,:,1));
for j = 3 : Ndim
    for i = L(j-1)+1 : L(j)
        I = unfoldi(i,j,L);
        halfnj = floor(Nway(j)/2)+1;
        if I(j) <= halfnj && I(j) >= 2
            U(:,:,i) = mu.*W(:,:,i)*(V(:,:,i)')*pinv(beta.*Iu+mu.*V(:,:,i)*V(:,:,i)');
        elseif I(j) > halfnj
            n_ = nc(I,j,Nway);
            i_ = foldi(n_,j,L);
            U(:,:,i) = conj(U(:,:,i_));
        end
        
    end
end
for jj = Ndim:-1:3
    U = ifft(U,[],jj);
end
U = real(U);
end