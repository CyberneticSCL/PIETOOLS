function [Zth,H] = poly_decomp_1D(Q)
% poly_decomp_1D creates a monomial-matrix decomposition of the matrix-valued polynomial
% Q in a single spatial variable such that Q(s)=H*Zth(s) where Z is the
% vector of n-dimensional monomial bases 
%
% Note: This representation is unique for the monomial basis Zth. The use 
% of alternative bases for the polynomial matrices is not currently
% supported.
%
% Note: To get a left decomposition, use [ZthT,HT] = poly_decomp_1D(Q') and 
% then Q=ZthT'*HT'
%
% INPUTS 
%   Q: A pvar object of dimension n times m with a single independent variable (there is currently no error
%   checking to ensure the absence of other independent variables.
%
% OUTPUT 
%   Zth - a pvar object of dimension q times m which is a matrix of vector-valued monomial bases
%   H   - a non-sparse matrix of dimension n times q

% NOTES:
% Distributed with DelayTOOLS
% Compatable with MULTIPOLY and SOSTOOLS as of June 2013
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial Coding:
% 6/26/20 - MMP  -  titled poly_decomp_1D
if ~isa(Q,'polynomial')
    [nd md]=size(Q);
    Zth=polynomial(eye(nd));
    H = Q;
    return
end
ss = Q.Varname(1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first step is to construct the separable representation of the
% Polynomial Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Q, we need to first convert to the forms
% Q=H Z(s) or for us this is Q=HZ(theta)
% however, we need to be consistent about our Z

Qt=Q;
[nd md]=size(Qt);
dQ=Qt.degmat;
dm=max(max(dQ));
smallZth=polynomial(eye(dm+1),[0:dm]',ss,[dm+1 1]); % common vector of monomials
bigZth=[];
for i=1:md
    bigZth=blkdiag(bigZth,smallZth);
end

% Now use bigZth to find H!
nZ=length(smallZth);
for i=1:nd
    for j=1:md
        [CQij,Ztemp,etemp] = poly2basis(Qt(i,j),smallZth); % examine each element of NQ
        bigC(i,(nZ*(j-1)+1):(nZ*j))=CQij.';
    end
end
H=full(bigC);
Zth=bigZth;
%H=bigC;
%Zth=bigZth;

% Q'=H*Zth so Q=Zth'*H'
% Zth*H-Q %uncomment to verify the representation (should be 0)
