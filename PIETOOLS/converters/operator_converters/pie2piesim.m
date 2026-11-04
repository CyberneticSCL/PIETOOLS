function PIE_out = pie2piesim(PIE)
% PIE_OUT = PIE2PIESIM(PIE) returns the linear PIE PIE as the struct that
% PIESIM and the executive witness (PIE_WITNESS_SOP) discretize: each
% operator a struct of its polynomial kernels, with the components of each
% kind of space (R, L2[s1], L2[s2], L2[s1,s2]) stacked in that order.
%
% INPUT
% - PIE:        'pie_struct' object, a linear PIE in at most 2 spatial
%               variables; initialized here if it is not;
%
% OUTPUT
% - PIE_out:    'struct' with fields
%   dom:        N x 2 'double', the interval of each of the N spatial
%               variables; [0,1] if N = 0;
%   vars:       N x 2 'polynomial', each spatial variable and its dummy
%               variable; [s1, s1_dum] if N = 0;
%   T, Tw, Tu, A, B1, B2, C1, D11, D12, C2, D21, D22:
%               'struct', the operators T, Tw, Tu, A, Bw, Bu, Cz, Dzw, Dzu,
%               Cy, Dyw and Dyu of the PIE. For N <= 1 the map
%               R^n0 x L2^n1[a,b] -> R^m0 x L2^m1[a,b],
%                 [x0; x1] -> [P x0 + int_a^b Q1(t) x1(t) dt;
%                              Q2(s) x0 + R.R0(s) x1(s)
%                                + int_a^s R.R1(s,t) x1(t) dt
%                                + int_s^b R.R2(s,t) x1(t) dt],
%               with fields P, Q1, Q2 and R (fields R0, R1, R2), each
%               'double' or 'polynomial' in s = var1 and t = var2;
%               dim = [m0 n0; m1 n1]; I = [a b]. For N = 2 the map between
%               R x L2[x] x L2[y] x L2[x,y], x and y the two variables in the
%               order of PIE.vars, with a field Rab per output space a and
%               input space b (a, b in {0, x, y, 2}); along each variable on
%               which both spaces depend the field is a cell, index 1 the
%               multiplier, 2 the integral from the lower bound to the
%               variable, 3 from the variable to the upper bound (Rxx, Rx2,
%               R2x 3 x 1; Ryy, Ry2, R2y 1 x 3; R22 3 x 3); dim 4 x 2, rows
%               R, L2[x], L2[y], L2[x,y], columns output, input; I 2 x 2;
%               var1, var2 2 x 1;
%   F:          [];
%   K:          'struct' as T, the operator K of the algebraic constraint;
%               a field only if the PIE has a constraint;
%   misc:       PIE.misc;
%   x_tab, w_tab, u_tab, y_tab, z_tab:
%               the component tables of the PIE, rows sorted by kind of
%               space as the rows and columns of the operators;
%   dim:        'double', the number N of spatial variables.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pie2piesim
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
% Initial coding MMP, 11/04/2026: the assembly of PIE2OPVAR (DJ, 10/29/2026)
%   and the kernel formulas of SOPVAR2OPVAR, SOPVAR2OPVAR2D and the
%   orientation rule of COPVAR2OPVAR2D, writing plain structs, so that PIESIM
%   needs no operator class.
% MMP, 11/04/2026: Moved from PIESIM/PIESIM_routines/PIESIM_Setup to
%   converters/operator_converters, since PIE_WITNESS_SOP calls it too.

% Operators as 'copvar' objects; at most 2 spatial variables
if ~isa(PIE,'pie_struct')
    error("The PIE should be a linear PIE, of type 'pie_struct'.")
elseif ~PIE.is_initialized
    PIE = initialize(PIE);
end
vnames = PIE.space.vars;
if numel(vnames)>2
    error("PIESIM simulates PIEs in at most 2 spatial variables.")
end

% Spatial variables and domain; without spatial variables, s1 on [0,1]
if isempty(vnames)
    vars = [polynomial({'s1'}),polynomial({'s1_dum'})];
    dom = [0,1];
else
    vars = PIE.vars;
    dom = PIE.dom;
end

% Each operator as a struct of kernels; B1 = Bw, B2 = Bu, C1 = Cz,
% D11 = Dzw, D12 = Dzu, C2 = Cy, D21 = Dyw, D22 = Dyu
PIE_out = struct('dom',dom,'vars',vars);
opnames = {'T','T';  'Tw','Tw'; 'Tu','Tu'};
for ii=1:size(opnames,1)
    PIE_out.(opnames{ii,1}) = assemble_op(PIE.(opnames{ii,2}),vnames,vars,dom);
end
PIE_out.F = [];
% Constraint K x = 0: a field only if present; input terms Kw, Ku have no
% field in this layout
if size(PIE.K,1)>0
    if ~(PIE.Kw==0) || ~(PIE.Ku==0)
        error("PIEs with an algebraic constraint involving inputs cannot be simulated.")
    end
    PIE_out.K = assemble_op(PIE.K,vnames,vars,dom);
end
opnames = {'A','A'; 'B1','Bw'; 'B2','Bu'; 'C1','Cz'; 'D11','Dzw'; 'D12','Dzu';
           'C2','Cy'; 'D21','Dyw'; 'D22','Dyu'};
for ii=1:size(opnames,1)
    PIE_out.(opnames{ii,1}) = assemble_op(PIE.(opnames{ii,2}),vnames,vars,dom);
end

% Component tables, rows sorted by kind of space (R, L2[s1], L2[s2],
% L2[s1,s2]) as the operators, the order within a kind kept
PIE_out.misc = PIE.misc;
for type = {'x','u','w','y','z'}
    tab = PIE.([type{1},'_tab']);
    kind = ones(size(PIE.space.(type{1}).depmat,1),1);
    for kk=1:numel(vnames)
        kind = kind + kk*PIE.space.(type{1}).depmat(:,kk);
    end
    [~,perm] = sort(kind);
    PIE_out.([type{1},'_tab']) = tab(perm,:);
end
PIE_out.dim = size(vars,1);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Pop = assemble_op(Pc,vnames,vars,dom)
% POP = ASSEMBLE_OP(PC,VNAMES,VARS,DOM) returns the 'copvar' operator PC as
% a struct of kernels, the spaces of each kind stacked in the order R,
% L2[s1], L2[s2], L2[s1,s2].
%
% INPUT
% - Pc:     'copvar' object;
% - vnames: 1 x N cellstr, the names of the spatial variables (N <= 2);
% - vars:   N x 2 'polynomial', the spatial and dummy variables;
% - dom:    N x 2 'double', the interval of each variable;
%
% OUTPUT
% - Pop:    'struct', the operator in the layout of PIE2PIESIM.

nvars = numel(vnames);
if nvars<=1
    nm = {'P','Q1';'Q2','R'};
else
    nm = {'R00','R0x','R0y','R02';
          'Rx0','Rxx','Rxy','Rx2';
          'Ry0','Ryx','Ryy','Ry2';
          'R20','R2x','R2y','R22'};
end

% Kind of each output and input space: 1 for R, 2 for L2[s1], 3 for
% L2[s2], 4 for L2[s1,s2]
[M,N] = size(Pc);
kind_out = ones(M,1);       kind_in = ones(N,1);
for kk=1:nvars
    is_var = strcmp(Pc.vars,vnames{kk});
    kind_out = kind_out + kk*any(Pc.space_out(:,is_var),2);
    kind_in = kind_in + kk*any(Pc.space_in(:,is_var),2);
end
% Spaces of one kind stacked: dimension of each kind, and the offset at
% which each space starts within its kind
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

% Zero operator of these dimensions; each block's kernels placed in the
% slot of its two kinds, at its offsets
Pop = zero_op(d,vars,dom);
for ii=1:M
    for jj=1:N
        if isempty(Pc.C{ii,jj})
            continue
        end
        slot = nm{kind_out(ii),kind_in(jj)};
        if nvars<=1
            val = kernels_1d(Pc.C{ii,jj});
        else
            val = kernels_2d(Pc.C{ii,jj},vnames,slot);
        end
        rows = off_out(ii) + (1:Pc.dim_out(ii));
        cols = off_in(jj) + (1:Pc.dim_in(jj));
        Pop.(slot) = place(Pop.(slot),val,rows,cols);
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Pop = zero_op(d,vars,dom)
% POP = ZERO_OP(D,VARS,DOM) returns the zero operator with dimensions D, in
% the layout of PIE2PIESIM.
%
% INPUT
% - d:      2 x 2 'double' (at most one spatial variable) or 4 x 2 (two),
%           the dimension of each kind of space, columns output, input;
% - vars:   N x 2 'polynomial', the spatial and dummy variables;
% - dom:    N x 2 'double', the interval of each variable;
%
% OUTPUT
% - Pop:    'struct', every kernel a zero of its size: 'double' for one
%           variable; 'polynomial' for two, and empty where both its
%           dimensions are zero (R00 then []).

if size(d,1)==2
    % One variable: P, Q1, Q2 and R.R0-R2 zero matrices
    R = struct('R0',zeros(d(2,1),d(2,2)),'R1',zeros(d(2,1),d(2,2)),...
               'R2',zeros(d(2,1),d(2,2)));
    Pop = struct('P',zeros(d(1,1),d(1,2)),'Q1',zeros(d(1,1),d(2,2)),...
                 'Q2',zeros(d(2,1),d(1,2)),'R',R);
else
    % Two variables: slot (a,b) a cell along x if spaces a and b both
    % depend on x (kinds 2 and 4), along y if both depend on y (3 and 4)
    nm = {'R00','R0x','R0y','R02';
          'Rx0','Rxx','Rxy','Rx2';
          'Ry0','Ryx','Ryy','Ry2';
          'R20','R2x','R2y','R22'};
    hasx = [0 1 0 1];       hasy = [0 0 1 1];
    Pop = struct();
    for a=1:4
        for b=1:4
            dab = [d(a,1),d(b,2)];
            if all(dab==0) && a==1 && b==1
                z = [];
            elseif all(dab==0)
                z = polynomial([]);
            else
                z = polynomial(zeros(dab));
            end
            shp = [1+2*(hasx(a)&&hasx(b)), 1+2*(hasy(a)&&hasy(b))];
            if isequal(shp,[1,1])
                Pop.(nm{a,b}) = z;
            else
                Pop.(nm{a,b}) = repmat({z},shp);
            end
        end
    end
end
Pop.I = dom;
Pop.var1 = vars(:,1);       Pop.var2 = vars(:,2);
Pop.dim = d;

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function val = kernels_1d(B)
% VAL = KERNELS_1D(B) returns the kernels of the 'sopvar' block B on at
% most one spatial variable s, as polynomials in s and its dummy s_dum.
%
% INPUT
% - B:      'sopvar' object, a map between R^n and L2^n[s];
%
% OUTPUT
% - val:    the matrix P ('double') of a map R -> R; Q1(s) of L2 -> R; Q2(s)
%           of R -> L2; for L2 -> L2, the 'struct' with fields R0(s), the
%           multiplier, R1(s,s_dum) on [a,s] and R2(s,s_dum) on [s,b].

% Variable of the block (s1 if none) and its dummy
sv = unique([B.vars.in,B.vars.out]);
if isempty(sv),     vname = 's1';   else,   vname = char(sv{1});    end
dname = [vname,'_dum'];

% Kernel (I_q x Z_L)' C (I_p x Z_R) of each parameter C, Z_L and Z_R the
% output and input monomials
if isempty(B.vars.in) && isempty(B.vars.out)
    % R -> R: the matrix itself
    val = B.params{1};
elseif isempty(B.vars.out)
    % L2 -> R: input monomials only
    Q1 = combine(quadPoly(B.params{1},{},B.ZR,B.dims,{},{vname},0));
    val = quadPoly.quadPoly2polynomial(Q1);
elseif isempty(B.vars.in)
    % R -> L2: output monomials only
    Q2 = combine(quadPoly(B.params{1},B.ZL,{},B.dims,{vname},{},0));
    val = quadPoly.quadPoly2polynomial(Q2);
else
    % L2 -> L2: the multiplier has its input variable identified with s;
    % the two integral kernels act on the dummy s_dum
    R0 = combine(quadPoly(B.params{1},B.ZL,B.ZR,B.dims,{vname},{vname},0));
    R1 = combine(quadPoly(B.params{2},B.ZL,B.ZR,B.dims,{vname},{dname},0));
    R2 = combine(quadPoly(B.params{3},B.ZL,B.ZR,B.dims,{vname},{dname},0));
    val = struct('R0',combine(quadPoly.quadPoly2polynomial(R0)),...
                 'R1',combine(quadPoly.quadPoly2polynomial(R1)),...
                 'R2',combine(quadPoly.quadPoly2polynomial(R2)));
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function val = kernels_2d(B,vnames,slot)
% VAL = KERNELS_2D(B,VNAMES,SLOT) returns the kernels of the 'sopvar' block
% B on at most two spatial variables, with x = VNAMES{1} and y = VNAMES{2},
% for the field SLOT of the 2-D layout of PIE2PIESIM.
%
% INPUT
% - B:      'sopvar' object, a map between any two of R, L2[x], L2[y] and
%           L2[x,y];
% - vnames: 1 x 2 cellstr, the names of x and y;
% - slot:   char, the field Rab of B's output and input spaces;
%
% OUTPUT
% - val:    'polynomial' in the variables and the dummies <name>_dum, or a
%           cell of them with axis d along direction d (shape of SLOT).

vars = B.vars;
% Directions of B by its own rule: its variables in sorted order, a lone
% variable named 's2' or 'y' along y
allv = unique([vars.in(:).',vars.out(:).']);
if isempty(allv)
    xy = {'',''};
elseif isscalar(allv) && (strcmp(allv{1},'s2') || strcmp(allv{1},'y'))
    xy = {'',allv{1}};
elseif isscalar(allv)
    xy = {allv{1},''};
else
    xy = {allv{1},allv{2}};
end
[mapsx, x_idx_in ] = name_at(vars.in ,xy{1});
[maps2x,x_idx_out] = name_at(vars.out,xy{1});
[mapsy, y_idx_in ] = name_at(vars.in ,xy{2});
[maps2y,y_idx_out] = name_at(vars.out,xy{2});

% Field of B in its own directions (kind 1 R, 2 L2[x], 3 L2[y], 4 L2[x,y])
nm = {'R00','R0x','R0y','R02';
      'Rx0','Rxx','Rxy','Rx2';
      'Ry0','Ryx','Ryy','Ry2';
      'R20','R2x','R2y','R22'};
cidx = [mapsx,mapsy]*[1;2]+1;
ridx = [maps2x,maps2y]*[1;2]+1;
hasx = [0 1 0 1];       hasy = [0 0 1 1];
shp = [1+2*(hasx(ridx)&&hasx(cidx)), 1+2*(hasy(ridx)&&hasy(cidx))];

% Variable names in the stored order of vars.out and vars.in (ZL{k} belongs
% to vars.out(k)): an input variable that is also an output variable is
% the dummy copy <name>_dum, any other keeps its name (integrated out)
var1 = vars.out(:).';
var2 = cell(1,numel(vars.in));
for k = 1:numel(vars.in)
    if ismember(vars.in{k},vars.out)
        var2{k} = [char(vars.in{k}),'_dum'];
    else
        var2{k} = char(vars.in{k});
    end
end

% Kernel of each parameter; index 1 along a direction is the multiplier,
% whose input variable is identified with the output one (params index
% (i1,i2) along the directions of B)
Rcell = cell(size(B.params));
for i=1:numel(B.params)
    [i1,i2] = ind2sub(size(B.params),i);
    var2_i = var2;
    if i1==1
        if mapsx && maps2x
            var2_i(x_idx_in) = var1(x_idx_out);
        elseif mapsy && maps2y
            var2_i(y_idx_in) = var1(y_idx_out);
        end
    end
    if i2==1 && (mapsx && maps2x)
        % With both directions shared, i2 indexes y
        if mapsy && maps2y
            var2_i(y_idx_in) = var1(y_idx_out);
        end
    end
    Ri = quadPoly(B.params{i},B.ZL,B.ZR,B.dims,var1,var2_i,0);
    Rcell{i} = combine(quadPoly.quadPoly2polynomial(Ri));
end
if isequal(shp,[1,1])
    val = Rcell{1};
else
    val = reshape(Rcell,shp);
end

% B's directions against (x,y): where B takes them in reverse order, its
% field is the mirror of SLOT (x and y exchanged) and its cell is transposed
sv = {'s1','s2'};
if mapsx,   sv{1} = char(vars.in{x_idx_in});    end
if maps2x,  sv{1} = char(vars.out{x_idx_out});  end
if mapsy,   sv{2} = char(vars.in{y_idx_in});    end
if maps2y,  sv{2} = char(vars.out{y_idx_out});  end
v = [vars.in, vars.out];
sw = ~isempty(v) && find(strcmp(sv,v{1}),1) ~= find(strcmp(vnames,v{1}),1);
mir = [1 3 2 4];
[ro,co] = find(strcmp(nm,slot));
if sw
    own = nm{mir(ro),mir(co)};
    if iscell(val),     val = val.';    end
else
    own = slot;
end
% A block whose own field is not SLOT would be dropped: refuse it
if ~strcmp(own,nm{ridx,cidx})
    error("Block of field %s does not map between the spaces of field %s.",nm{ridx,cidx},slot)
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [tf,idx] = name_at(list,nm)
% [TF,IDX] = NAME_AT(LIST,NM) returns whether the char NM is in the cellstr
% LIST and its position (false and 0 if NM is empty or absent).
if isempty(nm)
    tf = false;     idx = 0;
else
    [tf,idx] = ismember(nm,list);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function dest = place(dest,val,rows,cols)
% DEST = PLACE(DEST,VAL,ROWS,COLS) returns DEST with its elements ROWS,COLS
% set to VAL, element by element if DEST and VAL are cells or structs of
% the same structure.
%
% INPUT
% - dest:   'double', 'polynomial', cell or struct of these;
% - val:    the same structure as DEST, with blocks of size
%           numel(ROWS) x numel(COLS);
% - rows, cols: 'double' index vectors;
%
% OUTPUT
% - dest:   DEST with the blocks of VAL in place, 'polynomial' wherever VAL
%           is.
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
