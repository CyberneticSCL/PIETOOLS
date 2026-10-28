function [prog,Pop,Qcell,alpha_list] = possopvar(prog,dim,vars,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP,QCELL,ALPHA_LIST] = POSSOPVAR(PROG,DIM,VARS,DOM,DEG,OPTIONS)
% declares a positive semidefinite, self-adjoint 'sdopvar' decision
% operator on L_2^m[S3], S3 the spatial variables VARS, as the quadratic
% form
%
%       Pop = sum_{i,j} (Z_i)* (g Q_{ij}) (Z_j),     Q >= 0,
%
% in the basis operators
%
%       (Z_i x)(theta) = int I_alpha(theta-s) (I_m kron Z^i(theta,s))
%                            x(s) ds,
%
%       I_alpha(r) = prod_k I_{alpha_k}(r_k),
%       I_0(r) = delta(r),  I_1(r) = 1 for r >= 0,  I_-1(r) = 1 for r <= 0,
%
% the integral over the box of S3, I_1 and I_-1 zero elsewhere, alpha =
% alpha_i the multi-index of basis operator i, and Pop = Pop* >= 0. A
% multi-index stores 1, 2, 3 for the entries 0 (multiplier), 1 (lower
% integral) and -1 (upper integral), and 4 for the full-domain integral
% I_1 + I_-1 = 1. Z^i(theta,s) is the vector of the monomials in theta and
% s within the degree caps of basis operator i (degree 0 in s_k where
% alpha_k is a multiplier), and g(theta) the weight of OPTIONS.psatz, a
% multiplier on L_2[theta].
%
% INPUT
% - prog:       'struct', the LPI program to which the decision variables
%               of Q and the constraint Q >= 0 are added;
% - dim:        'double', m or [m,m], a positive integer;
% - vars:       1 x n3 'cellstr' of the distinct names of the variables S3;
%               or a 'char' name; or a 'polynomial' object (its variables);
%               or a 'double' scalar n3, for the names 's1',...,'sn3'.
%               Names ending in '_int' or '_dum' are an error;
% - dom:        'double', n3 x 2 array, dom(k,:) = [a_k,b_k], a_k < b_k,
%               the domain of vars{k}; or a 1 x 2 row for every variable;
% - deg:        the degree caps of the monomials Z^i(theta,s), over the
%               variables [theta_1..theta_n3, s_1..s_n3] in the order of
%               VARS: a 'double' scalar d, the cap of every theta_k and
%               every s_k; or a 'struct' with the optional fields
%                 int     scalar or 1 x n3, the cap in each theta_k
%                         (default 1);
%                 mult    scalar or 1 x n3, the cap in each s_k (default
%                         int);
%                 joint   scalar, the cap of the total degree (default
%                         none);
%                 subset  2^(2*n3) array, the cap of the total degree in
%                         each subset of [theta_1..theta_n3, s_1..s_n3],
%                         the subset with members b at index
%                         1+sum(2.^(b-1)) (default none);
%               or a 'cell' with one such scalar or 'struct' per row of
%               ALPHA_LIST;
% - options:    (optional) 'struct' with the optional fields
%                 psatz    'double', the weight g: 0 (default) g = 1; 1 g =
%                          prod_k (theta_k-a_k)(b_k-theta_k); 2j+1 or
%                          2j+2, j = 1..n3, for the j-th variable v of
%                          sort(VARS), g = (theta_v-a_v)/(b_v-a_v) or
%                          (b_v-theta_v)/(b_v-a_v). Any other value is an
%                          error;
%                 sep      'logical' scalar or 1 x n3 (default false):
%                          where sep(k), alpha_k takes the values 1 and 4,
%                          elsewhere 1, 2 and 3;
%                 include  the basis operators, among all multi-indices
%                          (enumerated with direction 1 fastest): an N x n3
%                          array of multi-indices, a 'logical' mask, or
%                          linear indices into that enumeration, at least
%                          one and all distinct (default all);
%                 type     'char', 'pos' or absent; any other value is an
%                          error.
%
% OUTPUT
% - prog:       'struct', PROG with the decision variables of Q and the
%               constraint Q >= 0;
% - Pop:        'sdopvar' object, m x m, on L_2^m[S3]; Pop = Pop* >= 0;
% - Qcell:      N x N 'cell', N the number of basis operators: Qcell{i,j}
%               the 'cellstr' array of the names of the entries of Q_{ij},
%               rows (component, monomial of Z^i) and columns (component,
%               monomial of Z^j), monomial fastest;
% - alpha_list: N x n3 'double', row i the multi-index alpha_i, in the
%               order of VARS and of Qcell.
%
% See also SOPQUADVAR, POSCOPVAR, EQ_OPTS_SOPVAR, DEGBALANCE,
% SETTINGS2POSSOPVAR, LPI_EQ_SDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - possopvar
%
% Copyright (C)2026 PIETOOLS Team
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
% MP, 08/22/2026:  Initial coding
% MMP, 09/21/2026: Body moved to 'sopquadvar' and this reduced to a wrapper.
%                  The body never imposed positivity itself -- it assembles
%                  a quadratic form, and the whole of the positivity
%                  question was the literal 'pos' passed to 'sosquadvar' --
%                  so it belonged in a constructor that Sec. 8.4's
%                  'poscopvar' can share, rather than being duplicated. The
%                  earlier header entries for the sep option, deg.subset,
%                  the sorted-S3 cell order and the block accumulation now
%                  live with the code in 'sopquadvar'. Behaviour here is
%                  unchanged, checked against the Gram dimensions measured
%                  beforehand and the 'test_possopvar' suite.
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  posmopvar -> poscopvar.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), the licence block in the toolbox form
%                  (title line added), step comments in the body; code
%                  unchanged. The deleted header text was the statement
%                  that this wraps sopquadvar, with the reference to it for
%                  the inputs and outputs, and the reason options.type is
%                  refused.

% OPTIONS: a 'struct', struct() when omitted; type 'pos' only
if nargin<6 || isempty(options)
    options = struct();
end
if ~isa(options,'struct')
    error("Options should be specified as a 'struct' object.")
end
if isfield(options,'type') && ~isempty(options.type) ...
        && ~strcmp(char(options.type),'pos')
    error("'possopvar' declares a positive operator; for type '%s' call "...
          +"'sopquadvar' directly.",char(options.type))
end
options.type = 'pos';

% Five inputs required; the quadratic form with Q >= 0
if nargin<5
    error("Not enough input arguments.")
end
[prog,Pop,Qcell,alpha_list] = sopquadvar(prog,dim,vars,dom,deg,options);

end
