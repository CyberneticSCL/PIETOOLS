function [Zth,H] = poly_decomp_1D_econ(Q)
% poly_decomp_1D_econ creates a monomial-matrix decomposition of the matrix-valued polynomial
% Q in a single spatial variable such that Q(s)=H*Zth(s) where Z is the
% smallest vector of n-dimensional monomial bases needed
%
% Note: This representation is unique for the monomial basis Zth. The use 
% of alternative bases for the polynomial matrices is not currently
% supported.
%
% Note: To get a left decomposition, use [ZthT,HT] = poly_decomp_1D_econ(Q') and 
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
% 6/26/20 - MMP  -  titled poly_decomp_1D_econ

if ~isa(Q,'polynomial')
    [nd md]=size(Q);
    Zth=polynomial(eye(md));
    H = Q;
    return
end
if length(Q.Varname)>1
    error('Q has more than 1 independent variable')
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
bZth=[];
for i=1:md
    sZth{i}=monomials(Qt(:,i));
    bZth=blkdiag(bZth,sZth{i});
    nZ(i)=length(sZth{i});
end


% Now use bigZth to find H!
for i=1:nd
    for j=1:md
        [CQij,Ztemp,etemp] = poly2basis(Qt(i,j),sZth{j}); % examine each element of NQ
        bigC(i,(sum(nZ(1:(j-1)))+1):sum(nZ(1:j)))=CQij.';
    end
end
H=full(bigC);
Zth=bZth;

% Q'=H*Zth so Q=Zth'*H'
% Zth*H-Q %uncomment to verify the representation (should be 0)
