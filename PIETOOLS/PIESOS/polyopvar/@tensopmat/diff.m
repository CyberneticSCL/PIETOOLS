function dCop = diff(Cop,vars,d)
% DCOP = DIFF(COP,VARS,D) takes the composition of the block matrix COP of
% tensor-PI operators with a Dth order differential operator in variables
% VARS
%   DCOP = d^D/(d vars^D) o COP
% VARS and D can be row vectors to simultaneously differentiate with
% respect to multiple variables. D can also be a column vector or matrix to
% take multiple derivatives, concatenated vertically.
%
% INPUTS
% - Cop:    m x n 'tensopmat' object representing a block matrix of
%           tensor-PI operators;
% - vars:   1 x r 'cellstr' or 'polynomial' (pvar) array specifying the
%           names of the spatial variables with respect to which to take
%           the derivative;
% - d:      (optional) p x r array of nonnegative integers specifying the
%           order of the derivatives to be taken with respect to each of
%           the spatial variables. Defaults to 1;
%
% OUTPUTS
% - dCop:   m*p x n 'tensopmat' representing the composition of the
%           requested differential operator with the operator Cop. The
%           derivatives specified by the rows of d are stacked vertically,
%           so that block (i,j) of the kth derivative is block
%           (i+(k-1)*P,j) of dCop, for P the number of row-blocks of Cop;
%
% NOTES
% - Each block is differentiated by its own 'diff' function: the blocks of
%   type 'tensopvar' by the product and Leibniz rules, those of type
%   'intop' (functionals) are mapped to zero, and blocks M of type
%   'double', 'polynomial' or 'dpvar', representing the operator x --> M*x
%   for x in \R, are differentiated elementwise. Empty blocks remain empty;
% - A block M(s) is differentiated as the function M(s)*x of s, so the
%   output variables of a row consisting only of such blocks are those
%   remaining in the derivative;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - diff
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 10/20/2026: Initial coding;

% Check that the operator is properly specified
if ~isa(Cop,'tensopmat')
    error("Block operator must be specified as object of type 'tensopmat'.")
elseif ~Cop.chkval
    Cop = initialize(Cop);
end

% Check that the variables and orders of differentiation are properly
% specified
if nargin<=2
    [vars,d] = parse_diff_inputs(vars);
else
    [vars,d] = parse_diff_inputs(vars,d);
end
p = size(d,1);

% Differentiation with no variables returns the original operator
[P,Q] = size(Cop.ops);
if isempty(vars)
    dCop = repmat(Cop,p,1);
    return
end

% Differentiate each block, for each of the rows of d
dops = cell(p*P,Q);
for k=1:p
    for i=1:P
        for j=1:Q
            op_ij = Cop.ops{i,j};
            if isempty(op_ij)
                continue
            elseif is_native_op(op_ij)
                % The block M is a 'double', 'polynomial' or 'dpvar' acting
                % on R: differentiate elementwise, one order at a time. A
                % 'double' is constant, so its derivative is zero. For a
                % single variable, diff of a 'dpvar' is the elementwise
                % derivative.
                M = op_ij;
                for v=find(d(k,:)>0)
                    if isa(M,'double')
                        M = zeros(size(M));
                    else
                        for l=1:d(k,v)
                            M = diff(M,vars{v});
                        end
                    end
                end
                dops{(k-1)*P+i,j} = M;
            else
                dops{(k-1)*P+i,j} = diff(op_ij,vars,d(k,:));
            end
        end
    end
end
dCop = tensopmat(dops,Cop.vars,Cop.dom);

end
