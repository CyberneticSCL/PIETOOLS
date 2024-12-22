%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elem_order.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 8_1_2024

function elem_order = vec_order2elem_order(vecs_order,vec_sizes)
% ELEM_ORDER = VEC_ORDER2ELEM_ORDER(VECS_ORDER,VEC_SIZES) takes a list of
% indices specifying an order of vectors, and determines the associated
% order of the elements in each vector
%
% INPUT
% - vecs_order:     nx1 array of corresponding to a permutation of the
%                   integers 1:n, specifying some order of vector-valued
%                   objects labeled 1 to n;
% - vec_sizes:      nx1 array specifying for each of the vector-valued
%                   objects 1 through n how many elements this vector
%                   contains.
%
% OUTPUT
% - elem_order:     mx1 array where m=sum(vec_sizes), specifying the order
%                   of the individual elements appearing in the
%                   vector-valued objects as per the specified vector
%                   order;
%
% EXAMPLE
% Suppose we have vectors a of size 1, b of size 2, and c of size three,
% building a total vector v = [a;b;c] = [a1;b1;b2;c1;c2;c3];
% Suppose we re-order them as v_new = [b;c;a] = [b1;b2;c1;c2;c3;a1];
% Then vecs_order = [2;3;1],    vec_sizes = [1;2;3], and
%   elem_order = [2;3;4;5;6;1];
% 

% Extract the inputs;
if ~isnumeric(vecs_order) || ~isnumeric(vec_sizes)
    error("The order and size of the vector-valued objects should be integer-valued.")
end
if numel(vecs_order)~=numel(vec_sizes)
    error("The number of vector sizes should match the number of vectors.")
end
n_comps = numel(vecs_order);

vec_sizes = vec_sizes(:);

% Determine the number of state variables in each component;
vec_size_cum = [0;cumsum(vec_sizes)];
% Set first and last index of variables in each component;
idcs = [vec_size_cum(1:end-1),vec_size_cum(2:end)];
% Determine all indices of variables that appear in each component;
idcs_c = mat2cell(idcs,ones(n_comps,1),2);
idcs_c = cellfun(@(x) (x(1)+1:x(2))',idcs_c,'UniformOutput',false);
% Reorder the variable indices to match the new order of the components;
idcs_reordered_c = idcs_c(vecs_order);
% Concatenate variable indices to get order of variables.
elem_order = cell2mat(idcs_reordered_c);

end