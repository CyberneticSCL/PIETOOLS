function [Pop_b2BC, Pop_f2BC, Pop_u2BC, Pop_w2BC, layout] = bc_operators_sopvar(PDE,bc_state_tab)
% Define 'copvar' objects to represent the BCs of the PDE in terms of the
% fundamental state xf and a set of core boundary states xb (as defined by
% bc_state_tab) as
%   0 = Pop_b2BC * xb + Pop_f2BC * xf + Pop_u2BC * u + Pop_w2BC * w;
% for a PDE in any number of spatial variables.
%
% INPUTS:
% - PDE:            An initialized 'pde_struct' whose BCs have been expanded
%                   in terms of core boundary and fundamental state
%                   components using "expand_PIETOOLS_PDE_BCs".
% - bc_state_tab:   Table defining the core boundary state, as returned by
%                   "ftc_expansion_sopvar".
%
% OUTPUTS:
% - Pop_b2BC, Pop_f2BC, Pop_u2BC, Pop_w2BC:
%               'copvar' objects such that the BCs are represented as
%                   0 = Pop_b2BC*xb + Pop_f2BC*xf + Pop_u2BC*u + Pop_w2BC*w;
%               The rows are ordered by the space the BC lives in, the
%               columns by the space the state lives in, as for the
%               operators of "ftc_expansion_sopvar". Pop_u2BC and
%               Pop_w2BC are returned as [] if there are no inputs u or w.
% - layout:     Struct with the ordered spaces (binary indices "masks") and
%               their dimensions ("dims") for the fields x, b, BC, u, w.
%
% NOTES:
% - This is the 'copvar' version of "construct_PI_ops_BCs" in
%   convert_PIETOOLS_PDE, in which the parameters of each term are placed
%   in the 'sopvar' parameter associated to the shared variables of the BC
%   and the state, so that any number of variables is supported.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - bc_operators_sopvar
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
% Initial coding, adapting construct_PI_ops_BCs to 'copvar', 10/28/2026

dom = PDE.dom;
vars = PDE.vars;
nvars = size(vars,1);
x_tab = PDE.x_tab;
pow2 = 2.^(0:nvars-1)';

% Layouts of the spaces in which the state, core boundary state, BCs and
% inputs exist
[x_masks,x_dims,x_cols] = space_layout(x_tab(:,2),logical(x_tab(:,3:2+nvars)),nvars);
[bc_masks,bc_dims,bc_rows] = space_layout(PDE.BC_tab(:,2),logical(PDE.BC_tab(:,3:2+nvars)),nvars);
[u_masks,u_dims,u_cols] = space_layout(PDE.u_tab(:,2),logical(PDE.u_tab(:,3:2+nvars)),nvars);
[w_masks,w_dims,w_cols] = space_layout(PDE.w_tab(:,2),logical(PDE.w_tab(:,3:2+nvars)),nvars);
% The core boundary state is already ordered by space
b_dep = logical(bc_state_tab(:,3:2+nvars));
nb_tab = size(bc_state_tab,1);
b_rmask = zeros(nb_tab,1);
if nb_tab>0 && nvars>0,     b_rmask = double(b_dep)*pow2;   end
b_masks = unique(b_rmask(:))';
b_dims = zeros(size(b_masks));
b_cols = cell(nb_tab,1);
for ii=1:numel(b_masks)
    for kk=find(b_rmask==b_masks(ii))'
        b_cols{kk} = b_dims(ii) + (1:bc_state_tab(kk,2));
        b_dims(ii) = b_dims(ii) + bc_state_tab(kk,2);
    end
end

layout = struct('x',struct('masks',x_masks,'dims',x_dims), ...
                'b',struct('masks',b_masks,'dims',b_dims), ...
                'BC',struct('masks',bc_masks,'dims',bc_dims), ...
                'u',struct('masks',u_masks,'dims',u_dims), ...
                'w',struct('masks',w_masks,'dims',w_dims));

% % % Loop over all BCs and terms, collecting the kernels
empty_entry = struct('r',{},'c',{},'rows',{},'cols',{},'sub',{},'C',{});
E_f = empty_entry;      E_b = empty_entry;
E_u = empty_entry;      E_w = empty_entry;
for BCnum = 1:numel(PDE.BC)
    in_BC = logical(PDE.BC_tab(BCnum,3:2+nvars));
    rmask = in_BC*pow2;
    rows = bc_rows{BCnum};
    for jj=1:numel(PDE.BC{BCnum}.term)
        term = PDE.BC{BCnum}.term{jj};
        Cval = polynomial(term.C);
        if isfield(term,'w') || isfield(term,'u')
            % % Input term: a finite-dimensional input
            if isfield(term,'w')
                comp = term.w;      tab = PDE.w_tab;    cols_all = w_cols;
            else
                comp = term.u;      tab = PDE.u_tab;    cols_all = u_cols;
            end
            cmask = double(tab(comp,3:2+nvars)>0)*pow2;
            Cval = to_primary(Cval,vars);
            e = struct('r',rmask,'c',cmask,'rows',rows,'cols',cols_all{comp}, ...
                        'sub',ones(1,nvars),'C',Cval);
            if isfield(term,'w'),   E_w(end+1) = e; else, E_u(end+1) = e; end %#ok<AGROW>
            continue
        end

        % % State term: fundamental or core boundary state
        comp = term.x;
        in_x = logical(x_tab(comp,3:2+nvars));
        Dmax = x_tab(comp,3+nvars:2+2*nvars);
        Dmax = Dmax(in_x);
        is_f = all(term.D==Dmax);
        if is_f
            % Fundamental state
            in_S = in_x;
            cmask = in_x*pow2;
            cols = x_cols{comp};
        else
            % Core boundary state
            is_intrr_full = zeros(1,nvars);
            if ~(isa(term.loc,'double') || isdouble(term.loc))
                is_intrr = zeros(1,sum(in_x));
                for ll=1:sum(in_x),     is_intrr(ll) = ~isdouble(term.loc(ll));     end
                is_intrr_full(in_x) = is_intrr;
            end
            Dval_full = zeros(1,nvars);
            Dval_full(in_x) = term.D;
            xcomp_info = [comp,x_tab(comp,2),is_intrr_full,Dval_full];
            kk = find(all(bc_state_tab==xcomp_info,2),1);
            if isempty(kk)
                error('A core boundary state component in the BCs could not be identified.')
            end
            in_S = logical(is_intrr_full);
            cmask = in_S*pow2;
            cols = b_cols{kk};
        end

        % Establish how each variable enters the operator: integrated over
        % (only in the state), shared, or only in the BC
        idx_xvar = cumsum(in_x);        % position of each variable among those of the state
        Idoms = term.I;
        tuples = ones(1,nvars);
        for v=1:nvars
            if in_S(v) && ~in_BC(v)
                % Integrated over the full domain: kernel in dummy variable
                Cval = subs(Cval,vars(v,1),vars(v,2));
            elseif in_S(v) && in_BC(v)
                % Shared variable: multiplier or partial integral
                ll = idx_xvar(v);
                Iv = [];
                if ~isempty(Idoms) && numel(Idoms)>=ll,     Iv = Idoms{ll};     end
                if isempty(Iv)
                    continue
                elseif isa(Iv,'double') || (isa(Iv,'polynomial') && isdouble(Iv))
                    % Full integral, as the sum of two partial integrals
                    t1 = tuples;    t2 = tuples;
                    t1(:,v) = 2;    t2(:,v) = 3;
                    tuples = [t1;t2];
                elseif isa(Iv,'polynomial') && isdouble(Iv(1))
                    tuples(:,v) = 2;
                elseif isa(Iv,'polynomial') && isdouble(Iv(2))
                    tuples(:,v) = 3;
                end
            elseif ~in_S(v) && in_BC(v)
                % Variable only in the BC: kernel uses the primary variable
                Cval = subs(Cval,vars(v,2),vars(v,1));
            end
        end
        for it=1:size(tuples,1)
            e = struct('r',rmask,'c',cmask,'rows',rows,'cols',cols, ...
                        'sub',tuples(it,:),'C',Cval);
            if is_f
                E_f(end+1) = e; %#ok<AGROW>
            else
                E_b(end+1) = e; %#ok<AGROW>
            end
        end
    end
end

Pop_f2BC = assemble_copvar_blocks(E_f,bc_masks,bc_dims,x_masks,x_dims,vars,dom);
Pop_b2BC = assemble_copvar_blocks(E_b,bc_masks,bc_dims,b_masks,b_dims,vars,dom);
if isempty(u_masks),    Pop_u2BC = [];
else,   Pop_u2BC = assemble_copvar_blocks(E_u,bc_masks,bc_dims,u_masks,u_dims,vars,dom);
end
if isempty(w_masks),    Pop_w2BC = [];
else,   Pop_w2BC = assemble_copvar_blocks(E_w,bc_masks,bc_dims,w_masks,w_dims,vars,dom);
end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function Cval = to_primary(Cval,vars)
% Replace all dummy variables by the primary variables.
for v=1:size(vars,1)
    Cval = subs(Cval,vars(v,2),vars(v,1));
end
end
