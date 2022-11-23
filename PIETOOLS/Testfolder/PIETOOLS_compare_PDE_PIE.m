function ismatch_arr = PIETOOLS_compare_PDE_PIE(PDE,PIE,tol)
% ismatch_arr = PIETOOLS_compare_PDE_PIE(PDE,PIE) checks that the PDE and
% PIE correspond to the same system.
% 
% INPUTS:
% - PDE:    A pde_struct object defining a PDE in the terms format
%           (see also the "@pde_struct/initialize" function).
% - PIE:    A pie_struct object, generated using the appropriate (1D or 2D)
%           PDE to PIE converter function on the PDE structure "PDE".
% - tol:    Allowed error in the PIE versus the PDE. That is, comparing
%           e.g. the x dynamics of the PDE to those of the PIE, the
%           polynomial x_pde - x_pie can have coefficients of which the
%           largest absolute value should not be greater than tol.
%
% OUTPUTS:
% - ismatch_arr:    A 5x1 logical array indicating whether:
%
%                   1 The PDE state xp := PIE.T*xf + PIE.Tw*w + PIE.Tu*u
%                     for a random fundamental state xf, and input w and u
%                     has an associated fundamental state xf_alt = D*xp
%                     that matches the randomly generated fundamental state
%                     xf. Here, the fundamental state should be the
%                     derivative of the PDE state in each spatial variable
%                     up to the maximally allowed degree, as per the
%                     continuity constraints in the PDE.
%                   
%                   2 The PDE state xp := PIE.T*xf + PIE.Tw*w + PIE.Tu*u
%                     satisfies the boundary conditions imposed in the PDE.
%
%                   3 The PIE dynamics xf_dyn = PIE.A*xf + PIE.B1*w + PIE.B2*u
%                     match the PDE dynamics of the PDE state
%                     xp := PIE.T*xf + PIE.Tw*w + PIE.Tu*u.
%                     as described by PDE.x.
%
%                   4 The regulated outputs z = PIE.C1*xf + PIE.D11*w + PIE.D12*u
%                     of the PIE match those of the PDE, as described by
%                     PDE.z.
%
%                   5 The observed outputs y = PIE.C2*xf + PIE.D21*w + PIE.D22*u
%                     of the PIE match those of the PDE, as described by
%                     PDE.y.
%
% NOTES:
% To check that the PDE and PIE match, a fundamental state xf, exogenous
% input w, and actuator input u are randomly generated. They are generated
% as polynomials, depending only on the spatial variables that the
% respective state component and input are allowed to depend on as
% indicated in the PDE and PIE structures. Coefficients in each polynomial
% are all integer.
% Using the randomly generated xf, w, and u, the PIE is applied to derive
% an associated PDE state xp, and it is checked that this PDE state returns
% the fundamental state xf, and satisfies the BCs. Then, the PIE is applied
% to (xf, w, u) to derive dynamics for xf, z, and y, and the PDE is
% applied to (xp, w, u) to derive dynamics for xp, z and y, which are then
% compared to verify the PDE and PIE match.
%
% The function should not be used to check equality of arbitrary PDE and
% PIE systems. Any mismatch in the dimensions, variables, domains, etc.
% will not be checked by the function, but will instead likely cause an
% internal error. As such, the PIE structure should be that generated from
% the converter files, applied to the input PDE structure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_compare_PDE_PIE(PDE,PIE)
%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 08/24/2022

% If no arguments are provided, try calling the function as a script.
if nargin==0
    if evalin('base',['exist(''PDE'',''var'');']) && evalin('base',['exist(''PIE'',''var'');'])
        evalin('base','ismatch_arr = PIETOOLS_compare_PDE_PIE(PDE,PIE);');
        return
    else
        error('No PDE and PIE have been specified!')
    end
elseif nargin<2
    error('No PDE and PIE have been specified!')
end
if nargin<3
    tol = 1e-12;    % Allowed error in the PDE and PIE results.
end

% Set some parameters.
maxdeg_comps = 2;   % Maximal degree of the monomials appearing in the randomly generated states/inputs.       



% Extract the spatial variables of the PIE.
% ! We're not checking that these match those of the PDE !
vars = PIE.vars;
nvars = size(vars,1);
nvars_pde = size(PDE.vars,1);   % The 2D converter might augment 1D systems to 2D ones, introducing variables that are not used.

% Extract the number of state and input components of the PIE.
% ! We're not checking that these match those of the PDE !
nx = size(PIE.x_tab,1);
nu = size(PIE.u_tab,1);
nw = size(PIE.w_tab,1);


% % % Generate a random fundamental state, exogenous input, and actuator
% % % input.
% Initialize an empty fundamental state.
xf = polynomial(zeros(sum(PIE.x_tab(:,2)),1));
kk = 0;
for ii=1:nx
    % For each state component, construct a random fundamental state
    % component, that depends on the allowed variables.
    nx_ii = PIE.x_tab(ii,2);
    has_vars_Lcomp = logical(PIE.x_tab(ii,3:2+nvars));
    if ~any(has_vars_Lcomp)
        xf_ii = polynomial(randi([1,10],[nx_ii,1]));
    else
        var1_Lcomp = vars(has_vars_Lcomp,1);
        xf_ii = PIETOOLS_random_poly_generator(nx_ii,1,var1_Lcomp,maxdeg_comps);
    end
    xf(kk+1:kk+length(xf_ii)) = xf_ii;
    kk = kk+length(xf_ii);
end
% Initialize an empty exogenous input.
w_pie = polynomial(zeros(sum(PIE.w_tab(:,2)),1));
kk = 0;
for ii=1:nw
    % Construct random exogenous inputs that depend on the allowed 
    % variables.
    nw_ii = PIE.w_tab(ii,2);
    has_vars_Lcomp = logical(PIE.w_tab(ii,3:2+nvars));
    if ~any(has_vars_Lcomp)
        w_ii = polynomial(randi([1,10],[nw_ii,1]));
    else
        var1_Lcomp = vars(has_vars_Lcomp,1);
        w_ii = PIETOOLS_random_poly_generator(nw_ii,1,var1_Lcomp,maxdeg_comps);
    end
    w_pie(kk+1:kk+length(w_ii)) = w_ii;
    kk = kk+length(w_ii);
end
% Initialize and empty actuator input.
u_pie = polynomial(zeros(sum(PIE.u_tab(:,2)),1));
kk = 0;
for ii=1:nu
    % Construct random actuator inputs that depend on the allowed 
    % variables.
    nu_ii = PIE.u_tab(ii,2);
    has_vars_Lcomp = logical(PIE.u_tab(ii,3:2+nvars));
    if ~any(has_vars_Lcomp)
        u_ii = polynomial(randi([1,10],[nu_ii,1]));
    else
        var1_Lcomp = vars(has_vars_Lcomp,1);
        u_ii = PIETOOLS_random_poly_generator(nu_ii,1,var1_Lcomp,maxdeg_comps);
    end
    u_pie(kk+1:kk+length(u_ii)) = u_ii;
    kk = kk+length(u_ii);
end

% Using the PIE states and inputs, generate an associated PDE state 
%   xp = T*x + Tw*w + Tu*u.
xp_pie = apply_opvar(PIE.T,xf) + apply_opvar(PIE.Tw,w_pie) + apply_opvar(PIE.Tu,u_pie);

% Keep track of which state, input and output variable in the PDE is
% associated to which state, input and output variable in the PIE.
% (The order is slightly adjusted in the PDE to PIE converter).
[~, pie2pde_indcs_x] = compute_PDE2PIE_conversion_indices(PDE,PIE,'x');     % xf(f2p_indcs_x) ~ xp;
[~, pie2pde_indcs_w] = compute_PDE2PIE_conversion_indices(PDE,PIE,'w');
[~, pie2pde_indcs_u] = compute_PDE2PIE_conversion_indices(PDE,PIE,'u');
[~, pie2pde_indcs_z] = compute_PDE2PIE_conversion_indices(PDE,PIE,'z');
[~, pie2pde_indcs_y] = compute_PDE2PIE_conversion_indices(PDE,PIE,'y');

% Convert the PIE states and inputs to associated pde versions.
xp_pde = xp_pie(pie2pde_indcs_x);
w_pde = w_pie(pie2pde_indcs_w);
u_pde = u_pie(pie2pde_indcs_u);



% ----------------------------------------------------------------------- %
%                       Start the checking process                        %
% ----------------------------------------------------------------------- %                       

% % % First, check that the composition of the differential operator and
% % % the PI operator T returns the fundamental state:
% % % xf = D * (T*x + Tw*w + Tu*u);

% Initialize an empty copy of the fundamental state.
xf_alt = polynomial(zeros(size(xf,1),1));
kk = 0;
for ii=1:nx
    % Differentiate the PDE state component up to the maximally allowed
    % order to retrieve the fundamental state component.
    nx_ii = PIE.x_tab(ii,2);
    x_indcs = kk+(1:nx_ii)';
    xdiff_ii = PIE.x_tab(ii,3+nvars:2+2*nvars);
    xf_alt_ii = xp_pie(x_indcs);
    for kk=1:nvars
        for dd=1:xdiff_ii(kk)
            xf_alt_ii = diff(xf_alt_ii,vars(kk,1));
        end
    end
    xf_alt(x_indcs) = xf_alt_ii;
    kk = x_indcs(end);
end
% Check that the obtained fundamental state matches the actual fundamental
% state.
err_Tops = cleanpoly(xf_alt-xf,tol);
if any(any(err_Tops.coeff))
    fprintf(2,'\n The observed fundamental state does not match the actual fundamental state for the constructed PIE. \n')
    max_err = max(max(abs(err_Tops.coeff)));
    fprintf(2,[' - A maximal discrepancy of ',num2str(max_err),' occurs in the coefficients. \n'])
    ismatch_xf_xp = false;
else
    ismatch_xf_xp = true;
end



% % % Next, check that the PDE state satisfies the BCs
BC_PDE = polynomial(zeros(sum(PDE.BC_tab(:,2)),1));
rr = 0;
for ii = 1:numel(PDE.BC)
    rr_indcs = rr+(1:PDE.BC_tab(ii,2))';
    has_vars_Lcomp = logical(PDE.BC_tab(ii,3:2+nvars_pde));
    for jj=1:numel(PDE.BC{ii}.term)
        term_jj = PDE.BC{ii}.term{jj};
        term_val_jj = apply_term(PDE,term_jj,has_vars_Lcomp,xp_pde,u_pde,w_pde);
        BC_PDE(rr_indcs) = BC_PDE(rr_indcs) + term_val_jj;
    end
    rr = rr_indcs(end);
end
err_BCs = cleanpoly(BC_PDE,tol);    % Boundary functions should be equal to 0.
if any(any(err_BCs.coeff))
    fprintf(2,'\n The PDE state associated to the PIE does not satisfy the desired boundary conditions. \n')
    max_err = max(max(abs(err_BCs.coeff)));
    fprintf(2,[' - A maximal discrepancy of ',num2str(max_err),' occurs in the coefficients. \n'])
    ismatch_BCs = true;
else
    ismatch_BCs = true;
end



% % % Next, check that the RHS of the PDE matches the RHS of the PIE.
x_dyn_pie = apply_opvar(PIE.A,xf) + apply_opvar(PIE.B1,w_pie) + apply_opvar(PIE.B2,u_pie);
x_dyn_pie = x_dyn_pie(pie2pde_indcs_x);
x_dyn_pde = polynomial(zeros(size(x_dyn_pie,1),1));
rr = 0;
for ii = 1:numel(PDE.x)
    rr_indcs = rr+(1:PDE.x_tab(ii,2))';
    has_vars_Lcomp = logical(PDE.x_tab(ii,3:2+nvars_pde));
    for jj=1:numel(PDE.x{ii}.term)
        term_jj = PDE.x{ii}.term{jj};
        term_val_jj = apply_term(PDE,term_jj,has_vars_Lcomp,xp_pde,u_pde,w_pde);
        x_dyn_pde(rr_indcs) = x_dyn_pde(rr_indcs) + term_val_jj;
    end
    rr = rr_indcs(end);
end
err_RHS_x = cleanpoly(x_dyn_pde-x_dyn_pie,tol);
if any(any(err_RHS_x.coeff))
    fprintf(2,'\n The observed PDE dynamics do not match the PIE dynamics. \n')
    max_err = max(max(abs(err_RHS_x.coeff)));
    fprintf(2,[' - A maximal discrepancy of ',num2str(max_err),' occurs in the coefficients. \n'])
    ismatch_x_dyn = false;
else
    ismatch_x_dyn = true;
end



% % % Next, check that regulated output of the PDE matches that of the PIE.
% Compute the value of the PIE output.
z_dyn_pie = apply_opvar(PIE.C1,xf) + apply_opvar(PIE.D11,w_pie) + apply_opvar(PIE.D12,u_pie);
% Order outputs in accordance with the PDE order.
z_dyn_pie = z_dyn_pie(pie2pde_indcs_z);
% Compute the value of the PDE output.
z_dyn_pde = polynomial(zeros(size(z_dyn_pie,1),1));
rr = 0;
for ii = 1:numel(PDE.z)
    rr_indcs = rr+(1:PDE.z_tab(ii,2))';
    has_vars_Lcomp = logical(PDE.z_tab(ii,3:2+nvars_pde));
    for jj=1:numel(PDE.z{ii}.term)
        term_jj = PDE.z{ii}.term{jj};
        term_val_jj = apply_term(PDE,term_jj,has_vars_Lcomp,xp_pde,u_pde,w_pde);
        z_dyn_pde(rr_indcs) = z_dyn_pde(rr_indcs) + term_val_jj;
    end
    rr = rr_indcs(end);
end
err_RHS_z = cleanpoly(z_dyn_pde-z_dyn_pie,tol);
if any(any(err_RHS_z.coeff))
    fprintf(2,'\n The regulated output of the PDE does not match that of the PIE. \n')
    max_err = max(max(abs(err_RHS_z.coeff)));
    fprintf(2,[' - A maximal discrepancy of ',num2str(max_err),' occurs in the coefficients. \n'])
    ismatch_z_dyn = false;
else
    ismatch_z_dyn = true;
end



% % % Finally, check that observed output of the PDE matches that of the PIE.
% Compute the value of the PIE output.
y_dyn_pie = apply_opvar(PIE.C2,xf) + apply_opvar(PIE.D21,w_pie) + apply_opvar(PIE.D22,u_pie);
% Order outputs in accordance with the PDE order.
y_dyn_pie = y_dyn_pie(pie2pde_indcs_y);
% Compute the value of the PDE output.
y_dyn_pde = polynomial(zeros(size(y_dyn_pie,1),1));
rr = 0;
for ii = 1:numel(PDE.y)
    rr_indcs = rr+(1:PDE.y_tab(ii,2))';
    has_vars_Lcomp = logical(PDE.y_tab(ii,3:2+nvars_pde));
    for jj=1:numel(PDE.y{ii}.term)
        term_jj = PDE.y{ii}.term{jj};
        term_val_jj = apply_term(PDE,term_jj,has_vars_Lcomp,xp_pde,u_pde,w_pde);
        y_dyn_pde(rr_indcs) = y_dyn_pde(rr_indcs) + term_val_jj;
    end
    rr = rr_indcs(end);
end
err_RHS_y = cleanpoly(y_dyn_pde-y_dyn_pie,tol);
if any(any(err_RHS_y.coeff))
    fprintf(2,'\n The observed output of the PDE does not match that of the PIE. \n')
    max_err = max(max(abs(err_RHS_y.coeff)));
    fprintf(2,[' - A maximal discrepancy of ',num2str(max_err),' occurs in the coefficients. \n'])
    ismatch_y_dyn = false;
else
    ismatch_y_dyn = true;
end



% Indicate in which aspects the PDE and PIE might differ.
ismatch_arr = [ismatch_xf_xp;
               ismatch_BCs;
               ismatch_x_dyn;
               ismatch_z_dyn;
               ismatch_y_dyn];

end



%%
function term_val = apply_term(PDE,term_jj,has_vars_Lcomp,xp,u,w)

vars = PDE.vars;
%dom = PDE.dom;
nvars = size(vars,1);

if isfield(term_jj,'w')
    % Extract the desired input component.
    Robj = term_jj.w;
    nnw_arr = cumsum([0;PDE.w_tab(:,2)]);
    Robj_indcs = (nnw_arr(Robj)+1:nnw_arr(Robj+1));
    term_val = polynomial(w(Robj_indcs));
    
    % Determine which variables the input component depends on.
    has_vars_Rcomp = logical(PDE.u_tab(Robj,3:2+nvars));
    var1_Rcomp = vars(has_vars_Rcomp,1);
    var2_Rcomp = vars(has_vars_Rcomp,2);
    
elseif isfield(term_jj,'u')
    % Extract the desired input component.
    Robj = term_jj.u;
    nnu_arr = cumsum([0;PDE.u_tab(:,2)]);
    Robj_indcs = (nnu_arr(Robj)+1:nnu_arr(Robj+1));
    term_val = polynomial(u(Robj_indcs));
    
    % Determine which variables the input component depends on.
    has_vars_Rcomp = logical(PDE.u_tab(Robj,3:2+nvars));
    var1_Rcomp = vars(has_vars_Rcomp,1);
    var2_Rcomp = vars(has_vars_Rcomp,2);
    
else
    % Extract the desired state component.
    Robj = term_jj.x;
    nnx_arr = cumsum([0;PDE.x_tab(:,2)]);
    Robj_indcs = (nnx_arr(Robj)+1:nnx_arr(Robj+1));
    term_val = polynomial(xp(Robj_indcs));
    
    % Determine which variables the state component depends on.
    has_vars_Rcomp = logical(PDE.x_tab(Robj,3:2+nvars));
    var1_Rcomp = vars(has_vars_Rcomp,1);
    var2_Rcomp = vars(has_vars_Rcomp,2);
    
    % Apply the desired derivative and delta operators
    Dval = term_jj.D;
    locval = term_jj.loc;
    for kk=1:size(var1_Rcomp,1)
        var1_kk = var1_Rcomp(kk);
        % Apply the derivative wrt var kk.
        for dd=1:Dval(kk)
            term_val = diff(term_val,var1_kk);
        end
        % Apply the delta operator.
        term_val = polynomial(subs(term_val,var1_kk,locval(kk)));        
    end
end
    
% % % Perform any specified integration.
% First, check for which variables no integration is performed.
use_int = true(1,length(term_jj.I));
for kk=1:size(term_jj.I,1)
    if isempty(term_jj.I{kk})
        use_int(kk) = false;
    end
end
use_int_full = false(1,nvars);
use_int_full(has_vars_Rcomp) = use_int;
% If both the LHS and the term depend on a variable, then any integral is
% (comprised of) a partial integral, and we have to use dummy variables to
% allow the use of a kernel.
Cval = term_jj.C;
use_dummy_vars = has_vars_Lcomp & use_int_full;
if ~nvars==0 && any(use_dummy_vars)
    term_val = subs(term_val,vars(use_dummy_vars,1),vars(use_dummy_vars,2));
    %Cval = subs(Cval,vars(~use_dummy_vars,2),vars(~use_dummy_vars,1));
end
term_val = Cval * term_val;
% Apply the desired integrals.
use_dummy_vars_Rcomp = use_dummy_vars(has_vars_Rcomp);
for kk=1:size(term_jj.I,1)
    if ~isempty(term_jj.I{kk})
        if use_dummy_vars_Rcomp(kk)
            var_kk = var2_Rcomp(kk);
        else
            var_kk = var1_Rcomp(kk);
        end
        term_val = int(term_val,var_kk,term_jj.I{kk}(1),term_jj.I{kk}(2));
    end    
end
% Finally, check that the term indeed does not depend on any illegal
% variables.
Lvars = vars(has_vars_Lcomp,1);
if isa(term_val,'polynomial') && any(~ismember(term_val.varname,Lvars.varname))
    error('The term appears to depend on a variable that the LHS does not depend on... Is something wrong in the initializer?')
end
    
end



%%
function [pde2pie_indcs, pie2pde_indcs] = compute_PDE2PIE_conversion_indices(PDE,PIE,obj)
% Determine indices to relate PIE state (or input or output) variables to
% PDE state (or input or output) variables.

ncomps = size(PDE.([obj,'_tab']),1);

nnpde_arr = cumsum([0;PDE.([obj,'_tab'])(:,2)]);
pde2pie_indcs = mat2cell([nnpde_arr(1:end-1),nnpde_arr(2:end)],ones(ncomps,1),2);
pde2pie_indcs = cellfun(@(x) (x(1)+1:x(2))',pde2pie_indcs,'UniformOutput',false);
pde2pie_indcs = cell2mat(pde2pie_indcs(PIE.([obj,'_tab'])(:,1)));    % xf ~ xp(p2f_indcs);

[~,pie_order] = sort(PIE.([obj,'_tab'])(:,1));
nnpie_arr = cumsum([0;PIE.([obj,'_tab'])(:,2)]);
pie2pde_indcs = mat2cell([nnpie_arr(1:end-1),nnpie_arr(2:end)],ones(ncomps,1),2);
pie2pde_indcs = cellfun(@(x) (x(1)+1:x(2))',pie2pde_indcs,'UniformOutput',false);
pie2pde_indcs = cell2mat(pie2pde_indcs(pie_order));    % xp ~ xf(f2p_indcs);

end