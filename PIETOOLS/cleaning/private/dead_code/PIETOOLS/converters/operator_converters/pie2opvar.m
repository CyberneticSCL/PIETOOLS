function PIE_old = pie2opvar(PIE)
% PIE_OLD = PIE2OPVAR(PIE) takes a 'pie_struct' object, representing a
% linear PIE, and returns a struct with the fields of the PIE structure
% used before the introduction of the current 'pie_struct', with the
% operators specified as 'opvar' (1D) or 'opvar2d' (2D) objects instead of
% 'copvar' objects.
%
% INPUT
% - PIE:    'pie_struct' object representing a PIE with at most 2 spatial
%           variables. If the PIE has not been initialized, it is
%           initialized first;
%
% OUTPUT
% - PIE_OLD:    struct with fields
%       dom, vars:      spatial domain, and spatial variables with
%                       associated dummy variables (N x 2 'polynomial'),
%                       as PIE.dom, PIE.vars;
%       T, Tw, Tu, A, B1, B2, C1, D11, D12, C2, D21, D22:
%                       'opvar' or 'opvar2d' objects representing the
%                       operators of the PIE, so that B1=PIE.Bw, B2=PIE.Bu,
%                       C1=PIE.Cz, D11=PIE.Dzw, D12=PIE.Dzu, C2=PIE.Cy,
%                       D21=PIE.Dyw, and D22=PIE.Dyu;
%       F:              [];
%       K:              the operator K defining an algebraic constraint, as
%                       'opvar' or 'opvar2d' object. Only a field if the
%                       PIE has such a constraint;
%       misc:           PIE.misc;
%       x_tab, w_tab, u_tab, y_tab, z_tab:
%                       tables with a row for each component, as PIE.x_tab,
%                       etc.;
%       dim:            number of spatial variables, as PIE.dim;
%
% NOTES
% - An 'opvar' or 'opvar2d' object maps between R^n x L2^m[s], or
%   R^n x L2^m1[s1] x L2^m2[s2] x L2^m3[s1,s2], so spaces of the same type
%   are combined, and the rows and columns of each operator are ordered as
%   [R; L2[s1]; L2[s2]; L2[s1,s2]]. The rows of each table are reordered
%   in the same manner, with the order within each type unchanged.
% - Operators for an algebraic constraint with input terms (Kw or Ku
%   nonzero) cannot be represented in the old format, and an error is
%   thrown. If the PIE has no constraint (K without rows), there is no
%   field K, as for the old PIE structure, so that e.g.
%   isfield(PIE_old,'K') is false.
% - An 'opvar' or 'opvar2d' object involves a dummy variable for each of
%   the spatial variables, named by adding '_dum' to the spatial variable.
% - A PIE that involves no spatial variables is represented using the
%   default variable and domain of an 'opvar' object.
%
% See also PIE2COPVAR, OPVAR2COPVAR, COPVAR2OPVAR, COPVAR2OPVAR2D.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pie2opvar
%
% Copyright (C)2026  PIETOOLS Team
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
% DJ, 10/29/2026: Initial coding;

if ~isa(PIE,'pie_struct')
    error("The PIE should be a linear PIE, of type 'pie_struct'.")
elseif ~PIE.is_initialized
    PIE = initialize(PIE);
end
vnames = PIE.space.vars;
if numel(vnames)>2
    error("PIEs in more than 2 spatial variables are not supported in the old format.")
end

% % % The spatial variables and domain
if isempty(vnames)
    % No spatial variables --> use the default of an opvar object.
    tmp = opvar();
    vars = [tmp.var1,tmp.var2];
    dom = tmp.I;
else
    vars = PIE.vars;
    dom = PIE.dom;
end

% % % The operators, using the order of the fields of the old pie_struct.
PIE_old = struct('dom',dom,'vars',vars);
opnames = {'T','T';  'Tw','Tw'; 'Tu','Tu'};
for ii=1:size(opnames,1)
    PIE_old.(opnames{ii,1}) = assemble_op(PIE.(opnames{ii,2}),vnames,vars,dom);
end
PIE_old.F = [];
if size(PIE.K,1)>0
    if ~(PIE.Kw==0) || ~(PIE.Ku==0)
        error("PIEs with an algebraic constraint involving inputs cannot be represented in the old format.")
    end
    PIE_old.K = assemble_op(PIE.K,vnames,vars,dom);
end
opnames = {'A','A'; 'B1','Bw'; 'B2','Bu'; 'C1','Cz'; 'D11','Dzw'; 'D12','Dzu';
           'C2','Cy'; 'D21','Dyw'; 'D22','Dyu'};
for ii=1:size(opnames,1)
    PIE_old.(opnames{ii,1}) = assemble_op(PIE.(opnames{ii,2}),vnames,vars,dom);
end

% % % The remaining fields.
PIE_old.misc = PIE.misc;
for type = {'x','u','w','y','z'}
    % Order the components by the type of space (R, L2[s1], L2[s2],
    % L2[s1,s2]), as for the rows of the operators.
    tab = PIE.([type{1},'_tab']);
    kind = ones(size(PIE.space.(type{1}).depmat,1),1);
    for kk=1:numel(vnames)
        kind = kind + kk*PIE.space.(type{1}).depmat(:,kk);
    end
    [~,perm] = sort(kind);
    PIE_old.([type{1},'_tab']) = tab(perm,:);
end
PIE_old.dim = size(vars,1);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function Pop = assemble_op(Pc,vnames,vars,dom)
% POP = ASSEMBLE_OP(PC,VNAMES,VARS,DOM) converts the 'copvar' object Pc to
% an 'opvar' (for at most 1 variable) or 'opvar2d' object, combining the
% spaces of the same type (R, L2[s1], L2[s2], L2[s1,s2]).
%
% INPUT
% - Pc:     'copvar' object;
% - vnames: 1 x N cellstr with the names of the spatial variables (N<=2);
% - vars:   N x 2 'polynomial' with the spatial and dummy variables;
% - dom:    N x 2 array with the domain of each variable;
%
% OUTPUT
% - Pop:    'opvar' or 'opvar2d' object on the specified variables and
%           domain, representing the same operator as Pc.

nvars = numel(vnames);
if nvars<=1
    Pop = opvar();
    nm = {'P','Q1';'Q2','R'};
else
    Pop = opvar2d();
    nm = {'R00','R0x','R0y','R02';
          'Rx0','Rxx','Rxy','Rx2';
          'Ry0','Ryx','Ryy','Ry2';
          'R20','R2x','R2y','R22'};
end
Pop.I = dom;
Pop.var1 = vars(:,1);       Pop.var2 = vars(:,2);

% Establish the type of each of the spaces: 1 for R, 2 for L2[s1], 3 for
% L2[s2], and 4 for L2[s1,s2].
[M,N] = size(Pc);
kind_out = ones(M,1);       kind_in = ones(N,1);
for kk=1:nvars
    is_var = strcmp(Pc.vars,vnames{kk});
    kind_out = kind_out + kk*any(Pc.space_out(:,is_var),2);
    kind_in = kind_in + kk*any(Pc.space_in(:,is_var),2);
end
% Combine the spaces of the same type, with the offset in the rows and
% columns of each type specifying where each space starts.
d = zeros(size(nm,1),2);
off_out = zeros(M,1);       off_in = zeros(N,1);
for ii=1:M
    off_out(ii) = d(kind_out(ii),1);
    d(kind_out(ii),1) = d(kind_out(ii),1) + Pc.dim_out(ii);
end
for jj=1:N
    off_in(jj) = d(kind_in(jj),2);
    d(kind_in(jj),2) = d(kind_in(jj),2) + Pc.dim_in(jj);
end
Pop.dim = d;

% Convert each block, and place it in the appropriate slot of the operator.
for ii=1:M
    for jj=1:N
        if isempty(Pc.C{ii,jj})
            continue
        end
        Pblk = copvar(Pc.C(ii,jj));
        if nvars<=1
            Bop = copvar2opvar(Pblk);
        else
            Bop = copvar2opvar2d(Pblk,vnames);
        end
        slot = nm{kind_out(ii),kind_in(jj)};
        rows = off_out(ii) + (1:Pc.dim_out(ii));
        cols = off_in(jj) + (1:Pc.dim_in(jj));
        Pop.(slot) = place(Pop.(slot),Bop.(slot),rows,cols);
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function dest = place(dest,val,rows,cols)
% DEST = PLACE(DEST,VAL,ROWS,COLS) sets the elements ROWS,COLS of the
% matrix DEST equal to VAL, or does so for each of the elements of DEST and
% VAL, if these are cells or structs with the same structure.
if isstruct(val)
    for fname = fieldnames(val)'
        dest.(fname{1}) = place(dest.(fname{1}),val.(fname{1}),rows,cols);
    end
elseif iscell(val)
    for kk=1:numel(val)
        dest{kk} = place(dest{kk},val{kk},rows,cols);
    end
else
    if isa(val,'polynomial') && ~isa(dest,'polynomial')
        dest = polynomial(dest);
    end
    dest(rows,cols) = val;
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
