function [Z,C] = monomial_outer(Z1,Z2)
% Given two monomial vectors Z1, Z2 in multiple variables (but same order), this returns
% Z1*Z2' = (I_n1\otimes Z')*C, where ni =length(Zi), C size n1*n x n2
% n = length(Z).

nvars = numel(Z1);
Z = cell(1,nvars);

C = 1;
for i=1:nvars
    [Zi,Ci] = monomial_outer_single(Z1{i},Z2{i});
    Z{i} = Zi;
    C = kron(C,Ci);
end
% this is now ((I\otimes Z(s1)') \otimes(I\otimes Z(s2)')\otimes  ... \otimes (I\otimes Z(sn)'))*C
% we need (I\otimes Z(s1)'\otimes Z(s2)' \otimes.... Z(sn)')*Cnew. so
% rearrange rows of C.
n1 = cellfun(@numel,Z1);
nz = cellfun(@numel,Z);

% Current row layout of C:
%   [i1,z1,i2,z2,...,id,zd]
%
% Desired row layout:
%   [i1,i2,...,id,z1,z2,...,zd]

rowShape = reshape([n1; nz],1,[]);          % [n1(1),nz(1),...,n1(d),nz(d)]
newOrder = [1:2:2*nvars, 2:2:2*nvars];     % all i's first, then all z's

C = reshape(C,[rowShape,size(C,2)]);
C = permute(C,[newOrder,2*nvars+1]);
C = reshape(C,prod(n1)*prod(nz),[]);
end

function [Z,C] = monomial_outer_single(Z1,Z2)
% Given two monomial vectors Z1, Z2 in same variable, this returns
% Z1*Z2' = (I_n1\otimes Z')*C, where ni =length(Zi), C size n1*n x n2
% n = length(Z).

n1 = length(Z1);
n2 = length(Z2);
[Z,~,map] = unique(Z1+Z2');
n = length(Z);

% Build sparse coefficient matrix C
[I,J] = ndgrid(1:n1,1:n2);

rows = (I(:)-1)*n + map(:);
cols = J(:);
vals = ones(numel(rows),1);

C = sparse(rows,cols,vals,n1*n,n2);
end