function Pop = pi4_copvar_sop(P,Q1,Q2,R,dims,var1,var2,dom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POP = PI4_COPVAR_SOP(P,Q1,Q2,R,DIMS,VAR1,VAR2,DOM) returns the container
% of the 4-PI operator
%
%   Pop: R^p0 x L2^p1[s] -> R^q0 x L2^q1[s],   s in DOM = [a,b],
%   Pop*[v; x] = [ P v + int_a^b Q1(s) x(s) ds ;
%                  Q2(s) v + R0(s) x(s) + int_a^s R1(s,t) x(t) dt
%                                       + int_s^b R2(s,t) x(t) dt ],
%
% with one 'sopvar' block per nonzero component.
%
% INPUT
% - P:      q0 x p0 'double';
% - Q1:     q0 x p1 'double' or 'polynomial' in VAR1, or [] for zero;
% - Q2:     q1 x p0 'double' or 'polynomial' in VAR1, or [] for zero;
% - R:      1 x 3 cell {R0,R1,R2} of q1 x p1 'double' or 'polynomial' in
%           VAR1 (s) and VAR2 (t), [] for zero;
% - dims:   2 x 2 'double', [q0 p0; q1 p1];
% - var1:   'polynomial', the variable s;
% - var2:   'polynomial', the dummy variable t, named as VAR1 with '_dum'
%           appended (a 'sopvar' records no dummy name);
% - dom:    1 x 2 'double', the interval [a,b].
%
% OUTPUT
% - Pop:    'copvar' with a row per output space and a column per input
%           space of positive dimension, R^q0 and R^p0 before L2. A zero
%           component is a [] block, except one zero block in a row or
%           column that would have none (the container reads the space of
%           a row or column off a block). With no space on one side, the
%           container has no space there (an empty channel of a PIE).
%
% NOTES
% The blocks, monomial bases and zero blocks are those 'opvar2copvar'
% builds from the opvar with fields P, Q1, Q2, R.R0-R2, dim = DIMS, var1,
% var2 and I = DOM. Cost: a fixed operator, O(blocks x monomials).
%
% EXAMPLE - the Volterra operator x -> int_0^s x(t) dt on L2[0,1]:
%   pvar s1 s1_dum
%   V = pi4_copvar_sop([],[],[],{0,1,0},[0 0; 1 1],s1,s1_dum,[0 1]);
%
% See also MAT2COPVAR_SOP, INT2COPVAR, EYE_COPVAR_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pi4_copvar_sop
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
% Initial coding MMP, 11/04/2026. The operator builder of
%                convert_PIETOOLS_DDF and convert_PIETOOLS_DDE2PIE_legacy,
%                which built these operators as opvar objects.

v1 = var1.varname;      v2 = var2.varname;
if ~strcmp(v2{1},[v1{1},'_dum'])
    error('pi4_copvar_sop:dummy',['The dummy variable ''%s'' should be named ',...
          '''%s_dum'': a sopvar records no dummy name.'],v2{1},v1{1});
end
comp = {P, Q1; Q2, R};
rows = [dims(1,1)>0, dims(2,1)>0];      cols = [dims(1,2)>0, dims(2,2)>0];
% Blocks: the nonzero components, plus one zero block in each row, then
% each column, they leave empty
keep = false(2,2);
for i = find(rows),     for j = find(cols),     keep(i,j) = ~pi4_iszero(comp{i,j});    end,    end
ir = find(rows);        jc = find(cols);
for i = ir,     if ~isempty(jc) && ~any(keep(i,:)),     keep(i,jc(1)) = true;   end,    end
for j = jc,     if ~isempty(ir) && ~any(keep(:,j)),     keep(ir(1),j) = true;   end,    end
C = cell(nnz(rows),nnz(cols));
for a = 1:numel(ir)
    for b = 1:numel(jc)
        i = ir(a);      j = jc(b);
        if ~keep(i,j),  continue,   end                 % structurally zero block
        q = dims(i,1);  p = dims(j,2);
        c = comp{i,j};
        vars = struct('in',{{}},'out',{{}});
        sdom = struct('in',zeros(0,2),'out',zeros(0,2));
        ZL = cell(1,0);     ZR = cell(1,0);
        % An empty component is a zero matrix of the block size
        if i==1 && j==1                                 % R -> R: P v
            if isempty(c),  c = zeros(q,p);     end
            params = {double(polynomial(c))};
        elseif i==1                                     % L2 -> R: Q1
            if isempty(c),  c = zeros(q,p);     end
            Qf = quadPoly.polynomial2quadPoly(c,{},v1);
            params = {full(Qf.C)};      ZR = Qf.Zt;
            vars.in = v1;               sdom.in = dom;
        elseif j==1                                     % R -> L2: Q2
            if isempty(c),  c = zeros(q,p);     end
            Qf = quadPoly.polynomial2quadPoly(c,v1,{});
            params = {full(Qf.C)};      ZL = Qf.Zs;
            vars.out = v1;              sdom.out = dom;
        else                                            % L2 -> L2: R0, R1, R2
            Rf = cell(3,1);
            for k = 1:3
                ck = c{k};      if isempty(ck),     ck = zeros(q,p);    end
                Rf{k} = quadPoly.polynomial2quadPoly(ck,v1,v2);
            end
            [params,ZL,ZR] = pi4_params(Rf);
            vars.out = v1;      vars.in = v1;
            sdom.out = dom;     sdom.in = dom;
        end
        C{a,b} = sopvar(params,vars,ZL,ZR,sdom,[q,p]);
    end
end
% The spaces, stated from DIMS: a row or column without a block has none
% to read them off
isL2 = [false, true];
if any(rows & isL2) || any(cols & isL2)
    reg = {char(v1{1})};        rdom = double(dom);
    so = isL2(rows)';           si = isL2(cols)';
else
    reg = cell(1,0);            rdom = zeros(0,2);
    so = false(nnz(rows),0);    si = false(nnz(cols),0);
end
meta = struct('vars',{reg},'dom',rdom,'space_out',so,'space_in',si,...
              'dim_out',dims(rows,1),'dim_in',dims(cols,2));
Pop = copvar(C,meta);
end

function tf = pi4_iszero(c)
% TF = PI4_ISZERO(C) is true when the component C, a 'double', a
% 'polynomial' or a cell of them (R), is empty or all zero.
if iscell(c),   tf = all(cellfun(@pi4_iszero,c));   return,     end
if isempty(c),  tf = true;  return,     end
c = polynomial(c);
tf = isempty(c.coefficient) || ~any(c.coefficient(:));
end

function [C,Zs,Zt] = pi4_params(Q)
% [C,ZS,ZT] = PI4_PARAMS(Q) writes the kernels of one block over common
% monomial bases.
%
% INPUT
% - Q:      n x 1 cell of 'quadPoly' of equal size, in one primary and one
%           dummy variable;
%
% OUTPUT
% - C:      n x 1 cell of 'double', the coefficients of Q{k} on ZS, ZT;
% - Zs, Zt: 1 x 1 cells, the monomials in the primary and the dummy variable.
if isscalar(Q)
    Q = Q{1};
    C = Q.C;        C = {full(C)};
    Zs = Q.Zs;      Zt = Q.Zt;
    return
end
Ccell = cell(size(Q));
Zscell = cell(size(Q));     Ztcell = cell(size(Q));
for i=1:numel(Q)
    Ccell{i} = Q{i}.C;
    Zscell{i} = Q{i}.Zs;   Ztcell{i} = Q{i}.Zt;
end
% Common bases: the union of the monomials, Ps*Zs and Pt*Zt the originals
Zsfull = zeros(0,1);        Ztfull = zeros(0,1);
nZs = zeros(numel(Q),1);    nZt = zeros(numel(Q),1);
for i=1:numel(Q)
    nZs(i) = size(Zscell{i}{1},1);
    nZt(i) = size(Ztcell{i}{1},1);
    Zsfull = [Zsfull; Zscell{i}{1}];                                        %#ok<AGROW>
    Ztfull = [Ztfull; Ztcell{i}{1}];                                        %#ok<AGROW>
end
[Ps,Zs] = uniquerows_integerTable(Zsfull);
[Pt,Zt] = uniquerows_integerTable(Ztfull);
C = Ccell;
[m,n] = size(Q{1});
nnZs = [0;cumsum(nZs)];
nnZt = [0;cumsum(nZt)];
for i=1:numel(C)
    Ps_i = spIkron(m,Ps(nnZs(i)+1:nnZs(i+1),:));
    Pt_i = spIkron(n,Pt(nnZt(i)+1:nnZt(i+1),:));
    C{i} = full(Ps_i'*Ccell{i}*Pt_i);
end
Zs = {Zs};
Zt = {Zt};
end
