function C = plus(A,B)
% This performs addition assuming inputs are compatible
%
% Inputs: A, B to be added
%
% Outputs: C = A+B. 

% initialize added output class
C = B;

% Error handling: Checks to ensure A and B are compatible
if any(A.dims~=B.dims)
    error('Summands A and B have different dimensions')
end


% dimensions of params in A and B, vars_in and vars_out

for i=1:numel(A.params)  % linear indexing of multi-dimensional cell array
    C.params{i} = A.params{i}+C.params{i};  % adding quadpoly objects.
end
end
