function At = ctranspose(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pt = ctranspose(P) transposes a sopvar operator P: L_2^p[S1,S3] to L_2^q[S2,S3]
% Date: 1/16/26
% Version: 1.0
% 
% INPUT
% A: sopvar class object
%   of the form 
%         sum_alpha int_S1 int_S3dum  I_alpha(S_3-S_3dum) Z_1(S_2,S_3) C{\alpha} Z_2(S1,S_3dum)
%                         alpha \in \{0,1,-1\}^{n_3}
%
%   where S_2 is the variables in S_i (output) not in S_j (input)
%         S_1 is the variables in S_j (input) not in S_i (output)
%         S_3 is the variables common to S_j (input) and S_i (output)
%         S_3dum are dummy versions of the variables in S_3
% OUTPUT
% At: transpose of the input opvar with the same matlab structure as
%   At: L_2^q[S_2,S_3] to L_2^p[S1,...Sn]
%   of the form 
%         sum_alpha int_S2 int_S3dum  I_alpha(-S_3+S_3dum) Z_2(S_1,S_3) C{\alpha} Z_1(S_2, S_3dum)
%                         alpha \in \{0,1,-1\}^{n_3}
%
%   where S_1 is the variables in S_i (output) not in S_j (input)
%         S_2 is the variables in S_j (input) not in S_i (output)
%         S_3 is the variables common to S_j (input) and S_i (output)
%         S_3dum are dummy versions of the variables in S_3
%

% If A has the form 
% sum_alpha int_S1 int_S3dum  I_alpha(S_3-S_3dum) Z_d(S_2,S_3) C Z_d(S_1, S_3dum) 
% where alpha \in {0,1,-1}^n3
% Then At has the form 
% sum_alpha int_S2dum int_S3dum  I_alpha(-S_3+S_3dum) Z_d(S_1,S_3) C Z_d(S_2, S_3dum)
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

% NOTE, need verification that A is in correct format

% total vars in old S1,S2,S3 
nvars3 = numel(A.vars_S3);
S3Var = A.vars_S3;
S3Vardummy = strrep(S3Var,'s','t');

% For the adjoint, the parameters in P.params(alpha) are swapped to
% positions  P.params(-alpha) using indexing alpha \in {0,1,-1}^n3. However,
% the actual indexing is \alpha \in {1,2,3}^n3
% so, position [1 2 3] is swapped with [1 3 2]

idx = repmat({[1 3 2]}, 1, nvars3);
C = A.params;  % assuming params is indexed by alpha in {1,2,3]^n3
C = C(idx{:});  % swap the parameters along 2 and 3. 

% for each parameter in cell structure repeat
for i=1:numel(C)
    % C{i} has form Z_1(S_2,S_3) C Z_2(S_1, S_3dum) 
    matC = C{i}.C;
    leftZ = C{i}.Zs;
    rightZ = C{i}.Zt;
    leftVar = C{i}.ns;
    rightVar = C{i}.nt;
    dim = C{i}.dim;
    
    newrightVar = strrep(leftVar,'s','t');
    newleftVar = strrep(rightVar,'t','s');

    % % we need form Z_2(S_1,S_3) C^T Z_1(S_2,S_3dum)
    % [tf,loc]=ismember(newrightVar,S3Var); 
    % newrightVar(tf)=S3Vardummy(loc(tf)); % changes (S_2,S_3) to (S_2,S_3dum)
    % [tf,loc]=ismember(newleftVar,S3Vardummy); 
    % newleftVar(tf) = S3Var(loc(tf)); % changes (S_1,S_3dum) to (S_1,S_3)

    C{i} = quadPoly(matC',rightZ,leftZ,dim',newleftVar,newrightVar);
end
At.params = C;
end
