function dV = Liediff(V,PIE)
% DV = LIEDIFF(V,PIE) computes the Lie derivative of a polynomial
% functional V along the PIE defined by PIE
%
% INPUTS
% - V:      1 x 1 'polyopvar' object representing the distributed
%           polynomial Lyapunov functional in terms of the PDE state, u
% - PIE:    'struct' with field 'T', a 'nopvar' object representing the map
%           from fundamental to PDE state, and field 'f', a 'polyopvar'
%           object representing the right-hand side of the PIE, so that
%               d/dt u = d/dt T*x = f(x)
%
% OUTPUTS
% - dV:     1 x 1 'polyopvar' object representing the derivative of the
%           Lyapunov functional along the PIE, now expressed in terms of
%           the PIE state.
%
% NOTES
% The Lyapunov functional must be specified in terms of the PDE state.
% If V.varname = {'x1'; ...; 'xn'}, this means PIE.T must be an n x n
% 'nopvar' object, representing the map [x1;...;xn]=T*[x1_f;...;xn_f] from
% fundamental state to PDE state. Note that dV.varname = V.varname, but now
% 'x1' in dV represents the fundamental state associated with 'x1' in V.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Liediff
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
% DJ, 01/26/2026: Initial coding

Top = PIE.T;
RHS = PIE.f;
Vdegs = V.degmat;
fdegs = RHS.degmat;
nvars = numel(V.varname);

if ~isequal(RHS.varname,V.varname)
    error("State variables in PIE and functional must match.")
end

% Decompose the T operator into individual operators acting on each state
% variable (this is really not ideal, but the only option for now...)
Top_cell = cell(nvars,nvars);
if isa(Top,'nopvar') || isa(Top,'ndopvar')
    Top_tmp = ndopvar2dopvar(Top);
else
    Top_tmp = Top;
end
for lindx=1:numel(Top_cell)
    [ridx,cidx] = ind2sub(size(Top_cell),lindx);
    Top_cell{lindx} = dopvar2ndopvar(Top_tmp(ridx,cidx));
end

% Replace the variable x by T*x in each of the monomials defining V. Since
% the monomials are expressed in terms of elements x_{i}, but T(i,:)*x may
% involve all variables again, this is somewhat cumbersome
old_state_arr = [];
new_state_arr = [];
old_monnum = [];
dtot = max(sum(Vdegs,2));
for ii=1:size(Vdegs,1)
    degs_ii = Vdegs(ii,:);
    % For each factor in the monomial, establish which PDE state variable
    % it corresponds to
    state_idcs_ii = (1:numel(degs_ii));
    state_idcs_ii = repelem(state_idcs_ii,degs_ii);
    % Apply the map ui = T(i,:)*v = sum_{j=1}^{n} T(i,j)*vj
    % splitting the monomial into j^dtot terms
    new_state_idcs_jj = state_idcs_ii;
    for jj=1:numel(state_idcs_ii)
        nr = size(new_state_idcs_jj,1);
        new_state_idcs_jj = repmat(new_state_idcs_jj,[nvars,1]);
        new_state_idcs_jj(:,jj) = kron((1:nvars)',ones(nr,1));
    end
    [nr,nc] = size(new_state_idcs_jj);
    new_state_arr = [new_state_arr; [new_state_idcs_jj, zeros(nr,dtot-nc)]];
    old_state_arr = [old_state_arr; [repmat(state_idcs_ii,[nr,1]), zeros(nr,dtot-nc)]];
    old_monnum = [old_monnum; ii*ones(nr,1)];
end

% Declare an empty polynomial 
tmp_poly = RHS;
tmp_poly.C.ops = {};
tmp_poly.degmat = zeros(0,nvars);
dV = tmp_poly;
for ii=1:size(new_state_arr,1)
    % Loop over each of the monomials in V, replacing x by T*x, and
    % evaluating the derivative d/dt T*x = f(x)
    degs_ii = Vdegs(old_monnum(ii),:);
    Kop_ii = V.C.ops{old_monnum(ii)};
    dtot_ii = sum(degs_ii);
    old_state_idcs = old_state_arr(ii,1:dtot_ii);
    state_idcs = new_state_arr(ii,1:dtot_ii);
    % Apply the product rule, d/dt(Tx*Tx) = Ax*Tx + Tx*Ax, looping over
    % each factor in the monomial and taking the temporal derivative only
    % along this factor
    dV_ii = tmp_poly;
    for jj=1:dtot_ii
        % Indicate the state variable for which to take the derivative by 0
        state_idcs_jj = state_idcs;
        state_idcs_jj(jj) = 0;
        % Build the product (T*x)*...*(T*x)*(A*x)*(T*x*...*T*x) for A in 
        % factor j
        Fx_jj = cell(1,dtot_ii);
        skip_jj = false;
        for state_num = 1:nvars
            Fx_tmp = tmp_poly;
            Fx_tmp.degmat = [zeros(1,state_num-1),1,zeros(1,nvars-state_num)];
            Fx_tmp.C.ops{1} = Top_cell{old_state_idcs(jj),state_num};
            is_state_num = state_idcs_jj==state_num;
            if any(is_state_num) && Top_cell{old_state_idcs(jj),state_num}==0
                skip_jj = true;
                break
            end
            Fx_jj(is_state_num) = {Fx_tmp};
        end
        if skip_jj
            continue
        end
        % Set factor jj as the right-hand side of the PIE
        dV_jj = tmp_poly;
        for ll=1:size(fdegs,1)
            % Loop over all terms in the PIE
            Fx_ll = RHS;
            op_ll = Fx_ll.C.ops{old_state_idcs(jj),ll};
            if isempty(op_ll) || (isa(op_ll,'nopvar') && op_ll==0)
                continue
            else
                Fx_ll.C.ops = Fx_ll.C.ops(old_state_idcs(jj),ll);
            end
            Fx_ll.degmat = Fx_ll.degmat(ll,:);
            Fx_jj{jj} = Fx_ll;
            % Take the composition of the functional operator with this
            % particular term
            dV_ll = mtimes_functional(Kop_ii,Fx_jj);
            dV_ll = combine_terms(dV_ll);
            dV_jj = dV_jj + dV_ll;
        end
        % Add the derivative for this particular factor to the full
        % derivative
        dV_ii = dV_ii + dV_jj;
    end
    dV = dV + dV_ii;
end

end