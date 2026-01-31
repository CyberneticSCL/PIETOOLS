function At = ctranspose(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pt = ctranspose(P) transposes a sopvar operator P: L_2^p[S1,S3] to L_2^q[S2,S3]
% Date: 1/16/26
% Version: 1.0
% 
% INPUT
% A: sopvar class object
%   of the form 
%         sum_alpha int_S1 int_S3dum  I_alpha(S_3-S_3dum) Z_d(S_2,S_3) C{\alpha} Z_d(S_3dum,S_1)
%                         alpha \in \{0,1,-1\}^{n_3}
%
%   where S_1 is the variables in S_i (output) not in S_j (input)
%         S_2 is the variables in S_j (input) not in S_i (output)
%         S_3 is the variables common to S_j (input) and S_i (output)
%         S_3dum are dummy versions of the variables in S_3% 
% OUTPUT
% At: transpose of the input opvar with the same matlab structure as
%   At: L_2^q[S_2,S_3] to L_2^p[S1,...Sn]
%   of the form 
%         sum_alpha int_S1 int_S3dum  I_alpha(S_3-S_3dum) Z_d(S_2,S_3) C{\alpha} Z_d(S_3dum,S_1)
%                         alpha \in \{0,1,-1\}^{n_3}
%
%   where S_1 is the variables in S_i (output) not in S_j (input)
%         S_2 is the variables in S_j (input) not in S_i (output)
%         S_3 is the variables common to S_j (input) and S_i (output)
%         S_3dum are dummy versions of the variables in S_3% 

%

% If A has the form sum_alpha int_S1 int_S3dum  I_alpha(S_3-S_3dum)
% Z_d(S_2,S_3) C Z_d(S_3dum,S_1) where alpha \in {0,1,-1}^n3
% Then At has the form sum_alpha int_S1 int_S3dum  I_alpha(S_3-S_3dum) Z_d(S_2,S_3) C Z_d(S_3dum,S_1)

% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - ctranspose
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

% initialize the transpose object
At = A;  
At.vars_S2 = A.vars_S1;  % input S1 vars become output vars S2
At.vars_S1 = A.vars_S2;  % output S2 vars become input vars S1
At.dom_2 = A.dom_1;  % input S1 doms become output S2 doms
At.dom_1 = A.dom_2;  % output S2 doms become input S1 doms

At.dims = A.dims';       % matrix dimension is transposed

% collect all vars
%varsMain = union(A.vars_in,A.vars_out); 

% create dummy vars corresponding to vars
varsDummy = strrep(varsMain,'s','t');

% NOTE, need verification that A is in correct format

% total vars in old S1,S2,S3 
nvars1 = numel(A.vars_S1);
nvars2 = numel(A.vars_S2);
nvars3 = numel(A.vars_S3);

% parameter cell dimension should not change

cellsize = ones(1,nvars3);
tmpsize = size(A.params);
%cellsize(1:length(tmpsize)) = tmpsize;  


% For the adjoint, the parameters in P.params(alpha) are swapped to
% positions  P.params(-alpha) using indexing alpha \in {0,1,-1}^n3. However,
% the actual indexing is \alpha \in {1,2,3}^n3
% so, position [1 2 3 3] is swapped with [1 3 2 2]

idx = repmat({[1 3 2]}, 1, nvars3);
C = C(idx{:});

% for each parameter in cell structure repeat
for i=1:numel(A.params)
    % find multiindex from linear index
    Aidx = cell(1,numel(cellsize));
    [Aidx{:}] = ind2sub(cellsize,i);
    
    % create transposed index to place in correct location
    % 2 and 3 are swapped; 1 is multiplier, 2 is lower int, 3 is upper int
    Atidx = cellfun(@(x) x+(x==2)-(x==3), Aidx, UniformOutput=false);
    % note, if x==1, then x + (x==2) - (x==3) = 1
    %          x==2, then x + (x==2) - (x==3) = 3
    %          x==3, then x + (x==2) - (x==3) = 2
    
    % transpose Parameter
    tmp = A.params{i}';
    
    % for each spatial variable si in S3, check if we need to swap
    % si and si_dum 
    for j=1:nvars3
        % dimension along si is 1, 
        % must be multiplier or full integral between different spaces, swap (si,si_dum)
        if cellsize(j)==1   
            tmp = var_swap(tmp,varsMain{j},varsDummy{j});
        elseif any(Atidx{j}==[2,3]) % must be semisep term, swap (si,si_dum)
            tmp = var_swap(tmp,varsMain{j},varsDummy{j});
        else 
            % multiplier term along si, appearing in both domain and range
            % do not swap
        end
    end
    % place the parameter in the transposed location
    At.params{Atidx{:}} = tmp;
end
end
