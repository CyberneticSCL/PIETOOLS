function obj = sopvar2nopvar(objsopvar)
% OBJ = SOPVAR2OPVAR(OBJSOPVAR) takes a sopvar object representing a 4-PI
% operator component and returns an nopvar object representing the same
% operator.
%
% INPUTS
% - objSopvar:  'sopvar' object representing a 1D PI operator. It can not
%               map between different function space, 
%               (i.e. it maps L2^n to L2^n);
%
% OUTPUTS
% - obj:        'nopvar' object representing the same operator as the input;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - sopvar2opvar
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
% MMP, 09/29/2026: degR is now read from ZR. It was read from ZL, so the
%   degree check compared degL with itself and never fired: ZL = 0:1 with
%   ZR = 0:2 returned a degree-1 nopvar without the s'^2 column. Also reject
%   any basis other than the complete set 0:deg, since the columns are
%   selected by position in that basis (ZL = [0 2], ZR = 0:2 returned the s^2
%   row as s^1). And reject content on dummy monomials of a multiplier
%   direction (params assigned past the canonical form), which was dropped
%   silently. All three mirror 'sdopvar2ndopvar', fixed 08/29/2026.


if ~isa(objsopvar, 'sopvar') 
    error('The input should be sopvar')
end



P = objsopvar;
dims = [P.dims(1), P.dims(2)];
if ~isequal(P.dom.in, P.dom.out)
    error('Input/Output domains dismatch')
end
if any(~strcmp(char(P.vars.in),char(P.vars.out)))
    error('Input/Output vars dismatch')
end
% ndopvar include dummy variables
% P.vars(:, 2) -- only dummy variables
vars_old = string(char(P.vars.in));

for idx = 1:length(vars_old)
    [var1, var1_dummy] = pvar(vars_old{idx}, [vars_old{idx}, '_dum']);
    if idx == 1
        vars_new = [var1, var1_dummy];
    else
        vars_new = [vars_new; [var1, var1_dummy]];
    end
end
% now construct monomial basis
% ZL = ZR;
ZL = P.ZL;
ZR = P.ZR;
degL = cellfun(@(x) max(x), ZL);
% degR = cellfun(@(x) max(x), ZL);                                          % MMP, 09/29/2026 (was)
% The degree of the right basis; ZL here made the check below vacuous.      % MMP, 09/29/2026
degR = cellfun(@(x) max(x), ZR);                                            % MMP, 09/29/2026
if ~isequal(degR, degL)
    error('left and right monomials have different degrees')
end
deg = degL;
% A nopvar has one degree per variable and the complete basis 0:deg on      % MMP, 09/29/2026
% both sides, and the columns below are picked by position in it, so any    % MMP, 09/29/2026
% other basis drops or misplaces coefficients. As in sdopvar2ndopvar.       % MMP, 09/29/2026
for kk = 1:numel(ZL)                                                        % MMP, 09/29/2026
    if ~isequal(reshape(ZL{kk},1,[]),0:degL(kk)) ...
            || ~isequal(reshape(ZR{kk},1,[]),0:degR(kk))                    % MMP, 09/29/2026
        error("Monomial bases must be the complete set 0:deg in every "...
              +"variable for conversion to an nopvar; use change_degree "...
              +"or pad the basis first.")                                   % MMP, 09/29/2026
    end                                                                     % MMP, 09/29/2026
end                                                                         % MMP, 09/29/2026
% now transform coefficient matrix 
% in ndopvar it is (Ik o [1;d])^T Cj in R^{k times }
% Cj has the size dim(1)*len(Zl)*(len(dvarnames)+1) times (len(Zr)*dim(2))
% Cj(1:dim(1)*len(Zl), :) are our A

N = size(P.dom.in,1);
left_size  = dims(1)*prod(cellfun(@(x) length(x), ZL));
right_size = dims(2)*prod(cellfun(@(x) length(x), ZR));

C_new = cell(size(P.params));
sz_C = size(C_new);
for ii=1:numel(P.params) 
    % new C is dim(1)*len(Zl)*(len(dvarnames)+1) times (len(Zr)*dim(2))
    % need to determinated which columns to choose
    % Determine the index of element ii along each dimension of the cell C
    idcs = cell(1,N);
    [idcs{:}] = ind2sub(sz_C,ii);
    idcs = cell2mat(idcs); % array of indeces 
    % If element ii corresponds to an integral, we need to account for the
    % monomial basis in the associated dummy variable
    is_int = logical(idcs-1);
    column_idx = 1;% column_idx indicates 
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
    [~, cols_mon, ~] = find(column_idx);

    C_sopvar = P.params{ii};
    % A multiplier direction carries no dummy monomials, so the columns     % MMP, 09/29/2026
    % outside cols_mon must be empty (as in sdopvar2ndopvar). Only a params % MMP, 09/29/2026
    % assignment that bypassed the canonical form can put content there,    % MMP, 09/29/2026
    % which the selection below would drop silently.                        % MMP, 09/29/2026
    other = true(1,size(C_sopvar,2));   other(cols_mon) = false;            % MMP, 09/29/2026
    if nnz(C_sopvar(:,other))                                               % MMP, 09/29/2026
        error("Parameter "+num2str(ii)+" depends on a dummy variable in a "...
              +"multiplier direction, which an 'nopvar' cannot represent; "...
              +"contract that direction first.")                            % MMP, 09/29/2026
    end                                                                     % MMP, 09/29/2026
    C_new{ii} = C_sopvar(:, cols_mon);
    % cdim = n*prod(deg(is_int)+1);
    % % Set sparse coefficients of dimension rdim x cdim
    % rho = (q+10)/(rdim*cdim);
    % Pop.C{ii} = sprand(rdim,cdim,rho);
end


obj = nopvar(); % empty nopvar
obj.dom =  P.dom.in;
obj.deg =  deg;
obj.vars = vars_new;
obj.C = C_new;

end