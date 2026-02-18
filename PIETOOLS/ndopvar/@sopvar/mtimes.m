function C = mtimes(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% C = mtimes(A,B) composes two sopvar operators C=AB where the output
% domain of B must be the input domain of A
% 
% Date: 1/16/26
% Version: 1.0
% 
% INPUT
% A: sopvar class object L_2^q[S4,S5] to L_2^r[S6,S5]
% B: sopvar class object L_2^p[S1,S3] to L_2^q[S2,S3]
%
% OUTPUT
% C = AB:  L_2^p[S7,S8] to L_2^r[S9,S8]
% 
% Where S7= S1, S3 cap S4                   %S1=S1, S3b= S3 cap S4 
%       S8= S3 cap S5                       %S3a= S3 cap S5 
%       S9= S6, S2 cap S5                   %S4=S6 ,S2a= S2 cap S5 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - mtimes
%
% Copyright (C)2026  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 1_16_2026
%
% Error handling: Checks to ensure A and B are compatible
if A.dims(1)~=B.dims(2)
    error('number of output components of B is different than number of input components of A');
end
if any(~strcmp(A.vars_in,B.vars_out))
    error('names of B output variables differ from names of A input variables');
end
if any(any(A.dom_out~=B.dom_in))
    error('output domain of B differs from input domain of A');
end

if numel(A.params)~=numel(B.params)   % why does this have to be true? -Sachin
    error('number of terms in summands is not equal -- one of them is probably malformed');
end


C = sopvar(B.vars_in,A.vars_out,[B.dims(1),A.dims(2)],B.dom_in,A.dom_out);


% Outline
% 1. Extract A_\alpha parameters from A and B_\beta parameters from B 
% 2. Extract Monomials and Coefficients from each set of parameters: C^A_\alpha and C^B_\beta
% 3a. Compute the D_{\alpha_a,\beta_b} 
% 3b. Compute the G_\alpha,\beta terms and express as Z C_{\alpha,\beta}^G Z 
% 4. Redistribute C^A_\alpha G_\alpha,\beta C^B_\beta as Z C^C_{\alpha,\beta} Z
% 5. Map \sum_{\alpha,\beta} C^C_{\alpha,\beta} to \sum_\gamma C^C_\gamma
% End: construct quadpolys from Z1,Z2, C^C_\gamma

% Now, just need to compute the parameters!
% I would suggest looping over the gammas and creating a set of terms to
% add for each gamma, since for some alpha,betas they get added to two of
% the gammas
alphaIdx = cell(1,ndims(A.params));
betaIdx = cell(1,ndims(B.params));
for i=1:numel(A.params)
    for j=1:numel(B.params)
        [alphaIdx{:}] = ind2sub(size(A.params),i);
        [betaIdx{:}] = ind2sub(size(B.params),j);
        [gammaIdx,Cparamstemp] = ...
            termCompose(A.params{i},B.params{j},alphaIdx,betaIdx,B.vars_in,A.vars_out,A.vars_in);
        C.params{gammaIdx} = C.params{gammaIdx}+Cparamstemp;
    end
end
end
