function obj = sdopvar2ndopvar(objsdopvar)
% OBJ = SDOPVAR2NDOPVAR(OBJSDOPVAR) converts an 'sdopvar' object to an
% 'ndopvar' object representing the same decision operator.
%
% INPUTS
% - objsdopvar: 'sdopvar' object whose input and output variables and
%               domains coincide (it maps L2^n[S] to L2^m[S]) and whose
%               monomial bases are the complete sets 0:deg in every
%               variable, with the same degree on both sides;
%
% OUTPUTS
% - obj:        'ndopvar' object representing the same operator, on the
%               decision variables objsdopvar.Zd;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - sdopvar2ndopvar
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
% AT, 08/05/2026: Initial coding 
% MMP, 08/29/2026: The monomial bases must be the complete sets 0:deg on
%                  both sides (an ndopvar holds one degree per variable),
%                  and the coefficients are mapped entry by entry into the
%                  ndopvar row layout, decision variable inner.
% MMP, 10/04/2026: Variable lists compared with 'isequal' on the cellstr
%                  rather than 'strcmp' on 'char' of it, as in 'plus'
%                  (09/07/2026): the char form depends on padding. Not a
%                  live fault. Ported from commit ca42b7f0 (09/21/2026).
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): help block
%                  corrected (input name, conditions), a log entry added for
%                  the 08/29/2026 lines, stale and commented-out lines
%                  removed, step comments; code unchanged.


if ~isa(objsdopvar, 'sdopvar') 
    error('The input should be sdopvar')
end



% An ndopvar has one variable list and one domain for both sides
P = objsdopvar;
dims = [P.dims(1), P.dims(2)];
dvarname = P.Zd;
if ~isequal(P.dom.in, P.dom.out)
    error('Input/Output domains dismatch')
end
% 'isequal' on the cellstr rather than 'strcmp' on 'char' of it, which      % MMP, 10/04/2026
% pads a ragged name list into a matrix (the test then depends on padding). % MMP, 10/04/2026
% if any(~strcmp(char(P.vars.in),char(P.vars.out)))                         % MMP, 10/04/2026 (was)
if ~isequal(P.vars.in(:),P.vars.out(:))                                     % MMP, 10/04/2026
    error('Input/Output vars dismatch')
end
% An ndopvar carries each variable with its dummy: vars(:,1) the primary
% variables, vars(:,2) the dummies
vars_old = string(char(P.vars.in));

for idx = 1:length(vars_old)
    [var1, var1_dummy] = pvar(vars_old{idx}, [vars_old{idx}, '_dum']);
    if idx == 1
        vars_new = [var1, var1_dummy];
    else
        vars_new = [vars_new; [var1, var1_dummy]];
    end
end
% The monomial bases: one degree per variable, equal on both sides
ZL = P.ZL;
ZR = P.ZR;
degL = cellfun(@(x) max(x), ZL);
degR = cellfun(@(x) max(x), ZR);                                            % MMP, 08/29/2026
if ~isequal(degR, degL)
    error('left and right monomials have different degrees')
end
deg = degL;
% An ndopvar carries a single degree per variable and always uses the
% complete monomial basis 0:deg, so a gapped basis, or a left basis
% differing from the right one, cannot be represented faithfully.
for kk = 1:numel(ZL)                                                        % MMP, 08/29/2026
    if ~isequal(reshape(ZL{kk},1,[]),0:degL(kk)) ...                        % MMP, 08/29/2026
            || ~isequal(reshape(ZR{kk},1,[]),0:degR(kk))                    % MMP, 08/29/2026
        error("Monomial bases must be the complete set 0:deg in every "...  % MMP, 08/29/2026
              +"variable for conversion to an ndopvar; use change_degree "... % MMP, 08/29/2026
              +"or pad the basis first.")                                   % MMP, 08/29/2026
    end                                                                     % MMP, 08/29/2026
end                                                                         % MMP, 08/29/2026
% now transform coefficient matrix 
% Cell j of the ndopvar is C_j with dim(1)*numel(Z_L)*(nd+1) rows, one
% block of nd+1 consecutive rows (constant first) per (matrix row,
% monomial), and one column per (matrix column, monomial in the dummy
% variables of its integral directions)

N = size(P.dom.in,1);
left_size_block = dims(1)*prod(cellfun(@(x) length(x), ZL));
left_size_full = left_size_block*(length(dvarname) + 1);
right_size = dims(2)*prod(cellfun(@(x) length(x), ZR));

A = P.params.A;
B= P.params.B;
sz_C = size(A);
n_dvarnames = length(dvarname);
C_new = cell(size(A));
for ii=1:numel(A) 
    % Multi-index of cell ii: entry 1 is a multiplier, 2 or 3 an integral
    idcs = cell(1,N);
    [idcs{:}] = ind2sub(sz_C,ii);
    idcs = cell2mat(idcs); 
    % An integral direction carries the full dummy basis 0:deg, a
    % multiplier direction only the constant: mark the columns of the
    % sdopvar coefficient that the ndopvar cell keeps
    is_int = logical(idcs-1);
    column_idx = 1; 
    for dim_idx = 1:length(deg)
        if is_int(dim_idx)
            column_idx = kron(column_idx, ones(1, deg(dim_idx) + 1));
        else
            var_temp= zeros(1, deg(dim_idx) + 1);
            var_temp(1) = 1;
            column_idx = kron(column_idx, var_temp); 
        end
    end
    column_idx = kron(ones(1, dims(2)), column_idx);
    cols_mon = find(column_idx(:));                                         % MMP, 08/29/2026
    col_map = zeros(right_size,1);   col_map(cols_mon) = 1:numel(cols_mon); % MMP, 08/29/2026

    % Map each stored coefficient to its place in the ndopvar cell.         % MMP, 08/29/2026
    % The rows of C_j are grouped as (I_k kron [1;d]): one block of q+1     % MMP, 08/29/2026
    % consecutive rows per (matrix row, monomial) slot, constant first, so  % MMP, 08/29/2026
    % the decision variable is the INNER index and not a trailing block.    % MMP, 08/29/2026
    Ak = A{ii};     Bk = B{ii};                                             % MMP, 08/29/2026
    if isempty(Ak),  Ak = sparse(left_size_block*right_size,1);  end        % MMP, 08/29/2026
    if isempty(Bk),  Bk = sparse(n_dvarnames,left_size_block*right_size); end % MMP, 08/29/2026
    if numel(Ak)~=left_size_block*right_size ...                            % MMP, 08/29/2026
            || size(Bk,1)~=n_dvarnames || size(Bk,2)~=numel(Ak)             % MMP, 08/29/2026
        error("Parameter "+num2str(ii)+" has inconsistent dimensions.")     % MMP, 08/29/2026
    end                                                                     % MMP, 08/29/2026
    [va,~,sa] = find(Ak(:));                                                % MMP, 08/29/2026
    [lb,vb,sb] = find(Bk);                                                  % MMP, 08/29/2026
    vpos = [va(:); vb(:)];                                                  % MMP, 08/29/2026
    lvl = [zeros(numel(va),1); lb(:)];      % 0 = constant term            % MMP, 08/29/2026
    val = [sa(:); sb(:)];                                                   % MMP, 08/29/2026
    ccol = floor((vpos-1)/left_size_block) + 1;                             % MMP, 08/29/2026
    rslot = vpos - (ccol-1)*left_size_block;                                % MMP, 08/29/2026
    % A multiplier direction carries no dummy monomials, so the columns     % MMP, 08/29/2026
    % outside cols_mon must be empty for the operator to be representable.  % MMP, 08/29/2026
    newcol = col_map(ccol);                                                 % MMP, 08/29/2026
    if any(newcol==0)                                                       % MMP, 08/29/2026
        error("Parameter "+num2str(ii)+" depends on a dummy variable in a "... % MMP, 08/29/2026
              +"multiplier direction, which an 'ndopvar' cannot represent; "... % MMP, 08/29/2026
              +"contract that direction first.")                            % MMP, 08/29/2026
    end                                                                     % MMP, 08/29/2026
    C_new{ii} = sparse((rslot-1)*(n_dvarnames+1)+1+lvl, newcol, val, ...    % MMP, 08/29/2026
                       left_size_full, numel(cols_mon));                    % MMP, 08/29/2026
end


% Assemble the ndopvar
obj = ndopvar(); % empty ndopvar
obj.dom =  P.dom.in;
obj.deg =  deg;
obj.vars = vars_new;
obj.dvarname = dvarname;
obj.C = C_new;

end
