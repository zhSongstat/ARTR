function V = Update_V(U,V,W,beta,mu)
%%
Ndim = length(size(V));
Nway= size(V);
Iv = eye(Nway(1));
L = ones(1,Ndim);
for x = 3:Ndim
    U = fft(U,[],x);
    V = fft(V,[],x);
    W = fft(W,[],x);
    L(x) = L(x-1) * Nway(x);
end
%%
%U(:,:,1) = mu.*W(:,:,1)*(V(:,:,1).')*pinv(beta.*Iu+mu.*V(:,:,1)*V(:,:,1).');
V(:,:,1) = pinv(mu.*U(:,:,1).'*U(:,:,1)+beta.*Iv)*U(:,:,1).'*(mu.*W(:,:,1));
for j = 3 : Ndim
    for i = L(j-1)+1 : L(j)
        I = unfoldi(i,j,L);
        halfnj = floor(Nway(j)/2)+1;
        if I(j) <= halfnj && I(j) >= 2
            V(:,:,i) = pinv(mu.*U(:,:,i).'*U(:,:,i)+beta.*Iv)*U(:,:,i).'*(mu.*W(:,:,i));
        elseif I(j) > halfnj
            n_ = nc(I,j,Nway);
            i_ = foldi(n_,j,L);
            V(:,:,i) = conj(V(:,:,i_));
        end
        
    end
end
for jj = Ndim:-1:3
    V = ifft(V,[],jj);
end
V = real(V);

end