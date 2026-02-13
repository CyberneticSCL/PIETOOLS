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
nvars = numel(V.varname);

% Declare a fundamental state variable for each PDE state variable
if ~isequal(RHS.varname,V.varname)
    error("State variables in PIE and functional must match.")
end
varname_pde = V.varname;
varname_pie = varname_pde;
for i=1:nvars
    varname_pie{i} = [varname_pie{i},'_f'];
end
if any(ismember(varname_pie,V.varname))
    error("Proposed fundamental state variable names match PDE state variable names...")
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

% Split the left- and right-hand side of the PIE into a separate equation 
% for each state variable
if size(RHS.C.ops,1)~=nvars
    error("Number of equations should match number of state variables.")
end
RHS_arr = cell(nvars,1);
LHS_arr = cell(nvars,1);
for i=1:nvars
    % For each PDE state variable x_i, define a polynomial expressing the
    % dynamics of this variable in terms of fundamental state variables,
    %       d/dt x_i = RHS_{i}(xf)
    RHS_arr{i} = RHS;
    RHS_arr{i}.C.ops = RHS.C.ops(i,:);
    RHS_arr{i}.varname = [varname_pde; varname_pie];
    RHS_arr{i}.degmat = [zeros(size(RHS.degmat,1),nvars),RHS.degmat];
    RHS_arr{i}.varmat = [RHS.varmat; RHS.varmat];

    % For each PDE state variable x_i, define a polynomial expressing this
    % state variable in terms of fundamental state variables,
    %       x_i = sum_{j=1}^{n} Top(i,j)*xf_{j}
    LHS_arr{i} = RHS_arr{i};
    LHS_arr{i}.degmat = [zeros(nvars),eye(nvars)];
    LHS_arr{i}.C.ops = Top_cell(i,:);
end

% Declare an empty polynomial 
tmp_poly = RHS;
tmp_poly.C.ops = {};
tmp_poly.degmat = zeros(0,nvars);
tmp_poly.varname = varname_pie;
dV = tmp_poly;

% For each of the monomials in V, take the Lie derivative along the PIE
for i=1:size(Vdegs,1)
    % Extract the ith term from the function V, expressing it in terms of
    % both the PDE and PIE state variables
    degi = Vdegs(i,:);
    dVi_tmp = tmp_poly;
    dVi_tmp.C.ops = V.C.ops(i);
    dVi_tmp.degmat = [degi,zeros(1,nvars)];
    dVi_tmp.varname = [varname_pde; varname_pie];
    dVi_tmp.varmat = [V.varmat;V.varmat];
    % For each factor in the monomial, replace only that factor by the
    % right-hand side of the PIE, replacing all other factors by T*xf
    dVi = tmp_poly;
    nfctrs = cumsum(degi);
    di = nfctrs(end);
    for j=1:di
        % Replace factor j-1 by T*xf
        if j>1
            % Determine which PDE state variable appears in factor j-1
            state_idx = find(j-1<=nfctrs,1,'first');
            % Express this PDE state in terms of the PIE state
            dV0 = 0;
            for trm_num=1:size(dVi_tmp.degmat,1)
                dVk = dVi_tmp;
                dVk.degmat = dVk.degmat(trm_num,:);
                dVk.C.ops = dVk.C.ops(:,trm_num);
                dVk = subs(dVk,1,LHS_arr{state_idx});
                dV0 = dV0 + dVk;
            end
            dVi_tmp = dV0;
        end
        dVj = dVi_tmp;
        % Replace remaining factors by T*xf
        for k=j+1:di
            % Determine which PDE state variable appears in factor j-1
            state_idx = find(k<=nfctrs,1,'first');
            % Express this PDE state in terms of the PIE state
            dV0 = 0;
            for trm_num=1:size(dVj.degmat,1)
                dVk = dVj;
                dVk.degmat = dVk.degmat(trm_num,:);
                dVk.C.ops = dVk.C.ops(:,trm_num);
                dVk = subs(dVk,2,LHS_arr{state_idx});
                dV0 = dV0 + dVk;
            end
            dVj = dV0;
        end
        % Determine which PDE state variable appears in factor j
        state_idx = find(j<=nfctrs,1,'first');
        % Substitute for right-hand side of the corresponding PIE 
        dV0 = 0;
        for trm_num=1:size(dVj.degmat,1)
            dVk = dVj;
            dVk.degmat = dVk.degmat(trm_num,:);
            dVk.C.ops = dVk.C.ops(:,trm_num);
            dVk = subs(dVk,1,RHS_arr{state_idx});
            dV0 = dV0 + dVk;
        end
        dVj = dV0;
        % Remove the PDE state variable names, and add to the full
        % derivative
        dVj.degmat = dVj.degmat(:,nvars+1:end);
        dVj.varname = varname_pie;
        dVj.varmat = dVj.varmat(1:nvars,:);
        dVi = dVi + dVj;
    end
    % Add to the full Lie derivative
    dV = dV+dVi;
end

% Replace the PIE variable names with the PDE variable names for
% consistency with the input
dV.varname = V.varname;

end