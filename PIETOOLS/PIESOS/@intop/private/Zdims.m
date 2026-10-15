function [nmons,strides,nZ] = Zdims(Zin)
% [NMONS,STRIDES,NZ] = ZDIMS(ZIN) returns, for the Kronecker product basis
% of monomials defined by the cell ZIN, the number of monomials NMONS(i)
% along each direction, the stride STRIDES(i) of each direction in the
% Kronecker product (first direction slowest), and the total number of
% monomials NZ = prod(NMONS).
%
% INPUTS
% - Zin:     1 x N cell, with element i a nmons(i) x deps(i) array, each
%            row specifying the degrees of one monomial along direction i;
%
% OUTPUTS
% - nmons:   1 x N array of integers, with nmons(i) the number of monomials
%            along direction i;
% - strides: 1 x N array of integers, with strides(i) the number of
%            monomials along the directions i+1,...,N, so that the index of
%            the monomial in the Kronecker product basis increases by
%            strides(i) if the monomial along direction i is incremented;
% - nZ:      scalar integer specifying the total number of monomials,
%            nZ = prod(nmons);
%
% NOTES
% - If N=0, nZ=1.
%

% DJ, 10/15/2026: Initial coding;

N = numel(Zin);
nmons = zeros(1,N);
for i=1:N
    nmons(i) = size(Zin{i},1);
end
strides = ones(1,N);
for i=N-1:-1:1
    strides(i) = strides(i+1)*nmons(i+1);
end
nZ = prod(nmons);

end
