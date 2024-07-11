function B=fcgltran4d(A)
% 4D Fast Chebyshev Transform from nodal to spectral
%
% Performs the fast transform of data sampled at the
% Chebyshev-Gauss-Lobatto Nodes x=cos(k*pi/N);
%
% A - original data in columns
% B - transformed data in columns
%
% Written by Yulia Peet 07/05/22  
% Contact: ypeet@asu.edu
[N,M,P,Q]=size(A);
    F=ifft([A(1:N,:,:,:);A(N-1:-1:2,:,:,:)]);
    B=real([F(1,:,:,:); 2*F(2:N,:,:,:)]);
    G=permute(B,[2 1 3 4]);
    F=ifft([G(1:M,:,:,:);G(M-1:-1:2,:,:,:)]);
    B=permute(real([F(1,:,:,:); 2*F(2:M,:,:,:)]),[2 1 3 4]);
    G=permute(B,[3 2 1 4]);
    F=ifft([G(1:P,:,:,:);G(P-1:-1:2,:,:,:)]);
    B=permute(real([F(1,:,:,:); 2*F(2:P,:,:,:)]),[3 2 1 4]);
    G=permute(B,[4 2 3 1]);
    F=ifft([G(1:Q,:,:,:);G(Q-1:-1:2,:,:,:)]);
    B=permute(real([F(1,:,:,:); 2*F(2:Q,:,:,:)]),[4 2 3 1]);
end
