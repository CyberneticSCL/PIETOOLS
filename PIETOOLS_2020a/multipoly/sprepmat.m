function b = sprepmat(a,m,n)
% function B = sprepmat(A,M,N)
%
% DESCRIPTION 
%   B is a sparse matrix consisting of an MxN tiling of the 
%   sparse matrix A.
%   [repmat fails on scalar, sparse matrices in Matlab 6.5:
%     >>repmat(sparse(5),2,3) 
%    This yields an error.]
%   
% INPUTS 
%   A: sparse matrix
%   M,N: Copies of A in the row and column directions
%
% OUTPUTS  
%   B: sparse matrix
%  
% SYNTAX 
%   B = sprepmat(A,M,N);
%

% 1/30/2003: PJS  Initial Coding  

[nra,nca] = size(a);
ridx = (1:nra)';
ridx = ridx(:,ones(1,m));
cidx = (1:nca)';
cidx = cidx(:,ones(1,n));
b = a(ridx,cidx);
