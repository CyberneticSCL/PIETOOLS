function At = ctranspose(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pt = ctranspose(P) transposes a sopvar operator P: L_2^p[Si,...Sn] to L_2^q[Sj,...,Sm]
% Date: 1/16/26
% Version: 1.0
% 
% INPUT
% A: sopvar class object
% 
% OUTPUT
% At: transpose of the input opvar with the same matlab structure as
%   At: L_2^q[Sj,...,Sm] to L_2^p[Si,...Sn]
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
At.vars_out = A.vars_in;  % input vars become output vars
At.vars_in = A.vars_out;  % output vars become input vars

At.dims = A.dims';       % matrix dimension is transposed

% collect all vars
varsMain = union(A.vars_in,A.vars_out); 
% create dummy vars corresponding to vars
varsDummy = strrep(varsMain,'s','t');

% total vars
nvars = numel(varsMain);
% parameter cell dimension
cellsize = ones(1,nvars);
tmpsize = size(A.params);
cellsize(1:length(tmpsize)) = tmpsize;  
% matlab size(A.params) function drops trailing 1s, i.e, 1x3x1 is truncated to 1x3. 
% In the above cellsize, we explicitly append them back for transpose logic.

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
    
    % for each spatial variable si, check if we need to swap
    % si and si_dum 
    for j=1:nvars
        % dimension along si is 1, 
        % must be multiplier or full integral, swap (si,si_dum)
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
