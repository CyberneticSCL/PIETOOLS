%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fcgltran2dperm(.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding unknown

function B=fcgltran2dperm(A,direction)
% 2D Fast Chebyshev Transform with matrix permutation 
%
% Performs the fast transform of data sampled at the
% Chebyshev-Gauss-Lobatto Nodes x=cos(k*pi/N);
%
% A - original data in columns
% B - transformed data in columns
% direction - set equal to 1 for nodal to spectral
%             anything else for spectral to nodal
%
% Written by Greg von Winckel 03/08/04  
% Contact: gregvw@chtm.unm.edu
[N,M]=size(A);
if direction==1 % Nodal-to-spectral
    F=ifft([A(1:N,:);A(N-1:-1:2,:)]);
    B=real([F(1,:); 2*F(2:N,:)]);
    G=permute(B,[2 1]);
    F=ifft([G(1:M,:);G(M-1:-1:2,:)]);
    B=permute(real([F(1,:); 2*F(2:M,:)]),[2 1]);
else            % Spectral-to-nodal
    F=fft([A(1,:); [A(2:N,:);A(N-1:-1:2,:)]/2]);
    B=real(F(1:N,:));
    G=B.';
    F=fft([G(1,:); [G(2:M,:);G(M-1:-1:2,:)]/2]);
    B=real(F(1:N,:)).';
end
