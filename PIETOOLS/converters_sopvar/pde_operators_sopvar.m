function ops = pde_operators_sopvar(PDE,Top)
% Compute the 'copvar' operators defining the PIE associated to a PDE,
%   (d/dt)*(Tx*xf + Tw*w + Tu*u) = A  * xf + B1  * w + B2  * u,
%                              z = C1 * xf + D11 * w + D12 * u,
%                              y = C2 * xf + D21 * w + D22 * u,
% for a PDE in any number of spatial variables.
%
% INPUTS:
% - PDE:    An initialized 'pde_struct' with ordered state components, as
%           passed to the conversion routine after "reorder_comps".
% - Top:    Struct with 'copvar' fields x, u, w, such that the PDE state x
%           is expressed in terms of the fundamental state xf and the inputs
%           as x = Top.x*xf + Top.u*u + Top.w*w, see "tmap_separable_sopvar".
%
% OUTPUTS:
% - ops:    Struct with fields A, B1, B2, C1, C2, D11, D12, D21, D22,
%           containing the 'copvar' operators defining the PIE. Operators
%           that have no input or no output space are [].
%
% NOTES:
% - Each term in the PDE (or output equation) takes the form
%       C * int [ d^D/ds^D x(s)|_{s=loc} ] ds
%   for a coefficient C, a derivative D, a location loc, and integral
%   domain I. We express the term as C*(Delta_loc * D_D * Tx) * xf, with the
%   composition of the derivative and evaluation operators with the 'copvar'
%   operator Tx a 'copvar' operator again, and apply the coefficient C using
%   the 'sopvar' parameters of the shared variables.
% - Derivatives of the inputs are not supported.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pde_operators_sopvar
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
% Initial coding, 'copvar' version of the PIE operators computed in
% convert_PIETOOLS_PDE, 10/28/2026

vars = PDE.vars;
dom = PDE.dom;
nvars = size(vars,1);
x_tab = PDE.x_tab;
names = pvar2varname(vars(:,1));    names = names(:)';

% Layout of the spaces of the state and inputs
[x_masks,x_dims,x_cols] = space_layout(x_tab(:,2),logical(x_tab(:,3:2+nvars)),nvars);
% Layout of the spaces of the fundamental state, i.e. the columns of Top.x,
% which may differ from that of the PDE state, e.g. for periodic BCs
[~,order_x] = sort(names);
xf_masks = zeros(1,numel(Top.x.dim_in));
for jj=1:numel(xf_masks)
    xf_masks(jj) = sum(2.^(order_x(Top.x.space_in(jj,:))-1));
end
xf_dims = Top.x.dim_in(:)';
lay_i.w = layout_of(PDE.w_tab,nvars);
lay_i.u = layout_of(PDE.u_tab,nvars);

ops = struct('A',[],'B1',[],'B2',[],'C1',[],'C2',[],'D11',[],'D12',[],'D21',[],'D22',[]);
objs = {'x','y','z'};
names_x = {'A','B1','B2'};  names_y = {'C2','D21','D22'};   names_z = {'C1','D11','D12'};
names_obj = {names_x,names_y,names_z};

% Keep track of the (differentiated, evaluated) rows of Top for each state
% component, to avoid recomputing them
cache = containers.Map('KeyType','char','ValueType','any');

for io=1:3
    obj = objs{io};
    if isempty(PDE.(obj)),  continue,   end
    tab = PDE.([obj,'_tab']);
    lay_o = layout_of(tab,nvars);
    o_dep = logical(tab(:,3:2+nvars));
    empty_e = struct('r',{},'c',{},'rows',{},'cols',{},'sub',{},'C',{});

    % Collect the direct input terms, and the state terms
    E_in.w = empty_e;       E_in.u = empty_e;
    sum_op = struct('x',[],'u',[],'w',[]);
    for eqnum=1:numel(PDE.(obj))
        in_eq = o_dep(eqnum,:);
        eq_mask = in_eq*(2.^(0:nvars-1)');
        rows = lay_o.cols{eqnum};
        for jj=1:numel(PDE.(obj){eqnum}.term)
            term = PDE.(obj){eqnum}.term{jj};
            Cval = polynomial(term.C);
            if isfield(term,'w') || isfield(term,'u')
                % % Term involving an input
                if isfield(term,'w'),   nm = 'w';   comp = term.w;     tab_i = PDE.w_tab;
                else,                   nm = 'u';   comp = term.u;     tab_i = PDE.u_tab;
                end
                in_i = logical(tab_i(comp,3:2+nvars));
                for v=1:nvars
                    if in_i(v) && ~in_eq(v)
                        Cval = subs(Cval,vars(v,1),vars(v,2));
                    else
                        Cval = subs(Cval,vars(v,2),vars(v,1));
                    end
                end
                E_in.(nm)(end+1) = struct('r',eq_mask,'c',in_i*(2.^(0:nvars-1)'), ...
                            'rows',rows,'cols',lay_i.(nm).cols{comp},'sub',ones(1,nvars),'C',Cval); %#ok<AGROW>
                continue
            end

            % % Term involving the state
            comp = term.x;
            in_x = logical(x_tab(comp,3:2+nvars));
            nxc = x_tab(comp,2);
            Dval_full = zeros(1,nvars);     Dval_full(in_x) = term.D;
            % Evaluation at the lower (-1), upper (1), or in the interior (0)
            loc_val = zeros(1,sum(in_x));
            Rloc = term.loc;
            Rdom = dom(in_x,:);
            if ~isa(Rloc,'double') && ~isdouble(Rloc)
                for kk=1:sum(in_x)
                    if isdouble(Rloc(kk)) && double(Rloc(kk))==Rdom(kk,1)
                        loc_val(kk) = -1;
                    elseif isdouble(Rloc(kk)) && double(Rloc(kk))==Rdom(kk,2)
                        loc_val(kk) = 1;
                    end
                end
            elseif ~isempty(Rloc)
                is_lower = double(Rloc)==Rdom(:,1)';
                loc_val(is_lower) = -1;
                loc_val(~is_lower) = 1;
            end
            loc_full = zeros(1,nvars);      loc_full(in_x) = loc_val;
            in_S = in_x & (loc_full==0);

            % Take the integrals and multipliers in the coefficient operator
            idx_xvar = cumsum(in_x);
            Idoms = term.I;
            tuples = ones(1,nvars);
            for v=1:nvars
                if in_S(v) && ~in_eq(v)
                    Cval = subs(Cval,vars(v,1),vars(v,2));
                elseif in_S(v) && in_eq(v)
                    ll = idx_xvar(v);
                    Iv = [];
                    if ~isempty(Idoms) && numel(Idoms)>=ll,     Iv = Idoms{ll};     end
                    if isempty(Iv)
                        continue
                    elseif isa(Iv,'double') || (isa(Iv,'polynomial') && isdouble(Iv))
                        t1 = tuples;    t2 = tuples;
                        t1(:,v) = 2;    t2(:,v) = 3;
                        tuples = [t1;t2];
                    elseif isa(Iv,'polynomial') && isdouble(Iv(1))
                        tuples(:,v) = 2;
                    elseif isa(Iv,'polynomial') && isdouble(Iv(2))
                        tuples(:,v) = 3;
                    end
                elseif ~in_S(v) && in_eq(v)
                    Cval = subs(Cval,vars(v,2),vars(v,1));
                end
            end
            cmask = in_S*(2.^(0:nvars-1)');
            ents = empty_e;
            for it=1:size(tuples,1)
                ents(end+1) = struct('r',eq_mask,'c',cmask,'rows',rows, ...
                            'cols',1:nxc,'sub',tuples(it,:),'C',Cval); %#ok<AGROW>
            end
            Coef = assemble_copvar_blocks(ents,lay_o.masks,lay_o.dims,cmask,nxc,vars,dom);

            % Compose the coefficient operator with the evaluated derivative
            % of the operators T
            key_base = sprintf('%d|%s|%s',comp,mat2str(Dval_full),mat2str(loc_full));
            for nm = {'x','u','w'}
                Tnm = Top.(nm{1});
                if isempty(Tnm) || is_zero_copvar(Tnm),     continue,   end
                key = [key_base,'|',nm{1}];
                if ~isKey(cache,key)
                    cache(key) = diff_eval_rows(Tnm,comp,x_masks,x_cols{comp},in_x, ...
                                        names,term.D,loc_val);
                end
                Dt = cache(key);
                if isempty(Dt),     continue,   end
                part = Coef*Dt;
                if isempty(sum_op.(nm{1})),     sum_op.(nm{1}) = part;
                else,   sum_op.(nm{1}) = sum_op.(nm{1}) + part;
                end
            end
        end
    end

    % Assemble the operators from the state to the object
    onames = names_obj{io};
    cols_x = {x_masks,x_dims};
    cols_n = {'x',xf_masks,xf_dims; 'w',lay_i.w.masks,lay_i.w.dims; 'u',lay_i.u.masks,lay_i.u.dims};
    for kk=1:3
        nm = cols_n{kk,1};
        cm = cols_n{kk,2};      cd = cols_n{kk,3};
        if isempty(cm),     continue,   end
        if strcmp(nm,'x'),  Pdir = assemble_copvar_blocks(empty_e,lay_o.masks,lay_o.dims,cm,cd,vars,dom);
        else,               Pdir = assemble_copvar_blocks(E_in.(nm),lay_o.masks,lay_o.dims,cm,cd,vars,dom);
        end
        if ~isempty(sum_op.(nm)),   Pdir = Pdir + sum_op.(nm);  end
        ops.(onames{kk}) = Pdir;
    end
end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function lay = layout_of(tab,nvars)
% Layout of the spaces in which the components in a table exist.
[lay.masks,lay.dims,lay.cols] = space_layout(tab(:,2),logical(tab(:,3:2+nvars)),nvars);
end

function Dt = diff_eval_rows(Top,comp,x_masks,rows,in_x,names,Dterm,loc_val)
% Compose the rows of the 'copvar' operator Top associated to the state
% component "comp" with a differential operator and evaluation operator,
%   Dt = Delta_loc * d^D/ds^D * Top(rows,:).
% The result has a single row space, in which the variables evaluated at a
% boundary are no longer present.
nvars = numel(in_x);
cmask = in_x*(2.^(0:nvars-1)');
i_c = find(x_masks==cmask);
ncols = numel(Top.dim_in);
blocks = cell(1,ncols);
xnames = names(in_x);
is_eval = loc_val~=0;
for j=1:ncols
    blk = Top.C{i_c,j};
    if isempty(blk),    continue,   end
    blk = sopvar_rows(blk,rows);
    if any(Dterm>0)
        blk = diff(blk,xnames,Dterm(:)');
    end
    for kk=find(is_eval)
        if is_zero_sopvar(blk),     break,      end
        if loc_val(kk)<0,   blk = trace_sopvar(blk,xnames{kk},'lower');
        else,               blk = trace_sopvar(blk,xnames{kk},'upper');
        end
    end
    % Blocks that vanish are structurally zero
    if is_zero_sopvar(blk),     continue,   end
    blocks{1,j} = blk;
end
% Space after evaluation
in_S = in_x;    in_S(in_x) = ~is_eval;
[~,order] = sort(names);
space_out = in_S(order);
meta = struct('vars',{Top.vars},'dom',Top.dom,'space_out',space_out,'space_in',Top.space_in, ...
                'dim_out',numel(rows),'dim_in',Top.dim_in);
Dt = copvar(blocks,meta);
end

function tf = is_zero_sopvar(P)
% True if the 'sopvar' operator is zero, or has no monomials in its kernel.
tf = any(cellfun(@isempty,P.ZL)) || any(cellfun(@isempty,P.ZR));
if tf,  return,     end
tf = true;
for m=1:numel(P.params)
    if any(any(P.params{m}~=0)),    tf = false;     return,     end
end
end

function tf = is_zero_copvar(P)
% True if all blocks of the container are structurally or numerically zero.
tf = true;
for k=1:numel(P.C)
    blk = P.C{k};
    if isempty(blk),    continue,   end
    for m=1:numel(blk.params)
        if any(any(blk.params{m}~=0)),  tf = false;     return,     end
    end
end
end
