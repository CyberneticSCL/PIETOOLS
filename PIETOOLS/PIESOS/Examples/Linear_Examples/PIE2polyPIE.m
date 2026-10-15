function PIE_out = PIE2polyPIE(PIE_in)
% PIE_OUT = PIE2POLYPIE(PIE_IN) takes a linear PIE in 'pie_struct' format
% and returns a structure representing the same PIE, with the right-hand
% side expressed as a (linear) distributed polynomial.
%
% INPUTS
% - PIE_in:     'pie_struct' object representing a linear PIE, with the
%               operators T and A in 'opvar' format (1D PIEs);
%
% OUTPUTS
% - PIE_out:    'struct' representing the same PIE as the input, with the
%               right-hand side of the PIE defined by a 'polyopvar'
%               object, PIE_out.f, so that
%                   d/dt T*x = f(x) = A*x;
%               It has fields T (a 'tensopmat'), f, x, varname_pde,
%               varname_pie, x_tab, vars, dom and dim, as the output of
%               `convert_PIETOOLS_PDE_nonlinear`, needed for e.g.
%               `piesos_LF`. The fundamental state is split into a state
%               collecting the finite-dimensional (ODE) components (if
%               any), followed by a scalar state for each of the
%               infinite-dimensional components, in the order of the PIE.
%               The state number i is named 'xi_pie' (PDE state 'xi'),
%               with i zero-padded to a common length;
%
% NOTES
% - Only the operators T and A of PIE_in are converted. The other
%   operators of PIE_in (Tw, Tu, B1, B2, C1, C2, D11-D22) are ignored.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIE2polyPIE
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
% DJ, 10/15/2026: Initial coding;

% Check that the PIE is a 1D PIE with operators in opvar format
Aop = PIE_in.A;     Top = PIE_in.T;
if ~isa(Aop,'opvar') || ~isa(Top,'opvar')
    error("PIE must have 1D operators T and A of type 'opvar'.")
end

% Declare the state components: one for the finite-dimensional states (if
% any), and a scalar state for each of the infinite-dimensional ones.
nx_ode = Aop.dim(1,1);      nx_pde = Aop.dim(2,1);
sizes = [nx_ode(nx_ode>0); ones(nx_pde,1)];
is_pde = [false(nx_ode>0,1); true(nx_pde,1)];
nstates = numel(sizes);
nnr = [0;cumsum(sizes)];

% Declare the PIE state variables
vars = PIE_in.vars;     dom = PIE_in.dom;
PIE_out = struct();
PIE_out.dom = dom;
PIE_out.vars = vars;
PIE_out.dim = size(vars,1);
PIE_out.x_tab = [(1:nstates)',sizes,double(is_pde)];
PIE_out.varname_pde = cell(nstates,1);
PIE_out.varname_pie = cell(nstates,1);
ndig = numel(num2str(nstates));       % zero-padded names, sorted in the order of the states
for i=1:nstates
    PIE_out.varname_pde{i} = ['x',num2str(i,['%0',num2str(ndig),'d'])];
    PIE_out.varname_pie{i} = ['x',num2str(i,['%0',num2str(ndig),'d']),'_pie'];
end
x = polyopvar(PIE_out.varname_pie,sizes,vars(:,1),dom,is_pde);
PIE_out.x = x;

% Declare the operators T and A as 'tensopmat', with a block for each pair
% of state components
Tcell = cell(nstates,nstates);      Acell = cell(nstates,nstates);
for i=1:nstates
    ridcs = nnr(i)+1:nnr(i+1);
    for j=1:nstates
        cidcs = nnr(j)+1:nnr(j+1);
        Tcell{i,j} = tensopmat(tensopvar(opvar2sopvar(clean_opvar(Top(ridcs,cidcs)))));
        Acell{i,j} = tensopmat(tensopvar(opvar2sopvar(clean_opvar(Aop(ridcs,cidcs)))));
    end
end
Tmat = vertcat_blocks(Tcell);
Amat = vertcat_blocks(Acell);

% Declare the PIE
PIE_out.T = Tmat;
PIE_out.f = Amat*x;

% nx = size(Aop,1);
% x_name = [repmat('x',[nx,1]),num2str((1:nx)')];
% x_name = mat2cell(x_name,ones(nx,1),size(x_name,2));
%
% % Declare the right-hand side of the PDE
% f = polyopvar(x_name,Aop.var1,Aop.dom);
% for i=1:nx
%     for j=1:nx
%         f.C.ops{i,j} = dopvar2ndopvar(Aop(i,j));
%     end
% end
%
% % Declare the PIE
% PIE_out = struct();
% PIE_out.T = dopvar2ndopvar(PIE_in.T);
% PIE_out.f = f;
% PIE_out.vars = PIE_in.vars;
% PIE_out.dom = PIE_in.dom;

end

function Cmat = vertcat_blocks(Ccell)
% Concatenate a cell of 'tensopmat' blocks into a single 'tensopmat'
rows = cell(size(Ccell,1),1);
for i=1:size(Ccell,1)
    row_i = [];
    for j=1:size(Ccell,2)
        row_i = [row_i,Ccell{i,j}]; %#ok<AGROW>
    end
    rows{i} = row_i;
end
Cmat = [];
for i=1:numel(rows)
    Cmat = [Cmat;rows{i}]; %#ok<AGROW>
end
end
