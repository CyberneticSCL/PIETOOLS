


% Generate a random PDE, and store in workspace in case of error.
%PDE = PIETOOLS_random_PDE_generator(2,3,[3,1,1,1,1]);
%assignin('base','PDE_test',PDE);
% Convert the PDE to a PIE, and also store in workspace in case of error.
%PIE = convert_PIETOOLS_PDE_2D(PDE);
%assignin('base','PIE_test',PIE);

maxdeg_comps = 2;
vars = PIE.vars;
dom = PIE.dom;
nvars = size(vars,1);

nx = size(PIE.x_tab,1);
nu = size(PIE.u_tab,1);
nw = size(PIE.w_tab,1);

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



% % % First, check that the composition of the differential operator and
% % % the PI operator T returns the fundamental state:
% % % xf = D * (T*x + Tw*w + Tu*u);
% Compute the (time-independent) PDE state.
xp_pie = apply_opvar(PIE.T,xf) + apply_opvar(PIE.Tw,w_pie) + apply_opvar(PIE.Tu,u_pie);
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

tol = 1e-12;
% Check that the obtained fundamental state matches the actual fundamental
% state.
err_Tops = cleanpoly(xf_alt-xf,tol);
if any(any(err_Tops.coeff))
    error('The observed fundamental state does not match the actual fundamental state for the constructed PIE')
end

% The fundamental state components in the PIE are ordered slightly
% differently from the PDE state components.
[p2f_indcs_x, f2p_indcs_x] = compute_PDE2PIE_conversion_indices(PDE,PIE,'x');
[p2f_indcs_w, f2p_indcs_w] = compute_PDE2PIE_conversion_indices(PDE,PIE,'w');
[p2f_indcs_u, f2p_indcs_u] = compute_PDE2PIE_conversion_indices(PDE,PIE,'u');
[p2f_indcs_z, f2p_indcs_z] = compute_PDE2PIE_conversion_indices(PDE,PIE,'z');
[p2f_indcs_y, f2p_indcs_y] = compute_PDE2PIE_conversion_indices(PDE,PIE,'y');

xp = xp_pie(f2p_indcs_x);
%w_pde = 0*w_pde;        u_pde = 0*u_pde;
w_pde = w_pie(f2p_indcs_w);
u_pde = u_pie(f2p_indcs_u);

% % % Next, check that the RHS of the PDE matches the RHS of the PIE.
RHS_PIE = apply_opvar(PIE.A,xf) + apply_opvar(PIE.B1,w_pie) + apply_opvar(PIE.B2,u_pie);
RHS_PIE_pde = RHS_PIE(f2p_indcs_x);
RHS_PDE = polynomial(zeros(size(RHS_PIE,1),1));
rr = 0;
for ii = 1:numel(PDE.x)
    rr_indcs = rr+(1:PDE.x_tab(ii,2))';
    has_vars_Lcomp = logical(PDE.x_tab(ii,3:2+nvars));
    for jj=1:numel(PDE.x{ii}.term)
        term_jj = PDE.x{ii}.term{jj};
        term_val_jj = apply_term(PDE,term_jj,has_vars_Lcomp,xp,u_pde,w_pde);
        RHS_PDE(rr_indcs) = RHS_PDE(rr_indcs) + term_val_jj;
    end
    rr = rr_indcs(end);
end
err_RHS_x = cleanpoly(RHS_PDE-RHS_PIE_pde,tol);
if any(any(err_RHS_x.coeff))
    error('The observed fundamental state does not match the actual fundamental state for the constructed PIE')
end


% % % Next, check that regulated output of the PDE matches that of the PIE.
% Compute the value of the PIE output.
z_PIE = apply_opvar(PIE.C1,xf) + apply_opvar(PIE.D11,w_pie) + apply_opvar(PIE.D12,u_pie);
% Order outputs in accordance with the PDE order.
z_PIE_pde = z_PIE(f2p_indcs_z);
% Compute the value of the PDE output.
z_PDE = polynomial(zeros(size(z_PIE,1),1));
rr = 0;
for ii = 1:numel(PDE.z)
    rr_indcs = rr+(1:PDE.z_tab(ii,2))';
    has_vars_Lcomp = logical(PDE.z_tab(ii,3:2+nvars));
    for jj=1:numel(PDE.z{ii}.term)
        term_jj = PDE.z{ii}.term{jj};
        term_val_jj = apply_term(PDE,term_jj,has_vars_Lcomp,xp,u_pde,w_pde);
        z_PDE(rr_indcs) = z_PDE(rr_indcs) + term_val_jj;
    end
    rr = rr_indcs(end);
end
err_RHS_z = cleanpoly(z_PDE-z_PIE_pde,tol);
if any(any(err_RHS_z.coeff))
    error('The observed fundamental state does not match the actual fundamental state for the constructed PIE')
end


% % % Next, check that observed output of the PDE matches that of the PIE.
% Compute the value of the PIE output.
y_PIE = apply_opvar(PIE.C2,xf) + apply_opvar(PIE.D21,w_pie) + apply_opvar(PIE.D22,u_pie);
% Order outputs in accordance with the PDE order.
y_PIE_pde = y_PIE(f2p_indcs_y);
% Compute the value of the PDE output.
y_PDE = polynomial(zeros(size(y_PIE,1),1));
rr = 0;
for ii = 1:numel(PDE.y)
    rr_indcs = rr+(1:PDE.y_tab(ii,2))';
    has_vars_Lcomp = logical(PDE.y_tab(ii,3:2+nvars));
    for jj=1:numel(PDE.y{ii}.term)
        term_jj = PDE.y{ii}.term{jj};
        term_val_jj = apply_term(PDE,term_jj,has_vars_Lcomp,xp,u_pde,w_pde);
        y_PDE(rr_indcs) = y_PDE(rr_indcs) + term_val_jj;
    end
    rr = rr_indcs(end);
end
err_RHS_y = cleanpoly(y_PDE-y_PIE_pde,tol);
if any(any(err_RHS_y.coeff))
    error('The observed fundamental state does not match the actual fundamental state for the constructed PIE')
end


% % % Finally, check that the PDE state also satisfies the BCs
BC_PDE = polynomial(zeros(sum(PDE.BC_tab(:,2)),1));
rr = 0;
for ii = 1:numel(PDE.BC)
    rr_indcs = rr+(1:PDE.BC_tab(ii,2))';
    has_vars_Lcomp = logical(PDE.BC_tab(ii,3:2+nvars));
    for jj=1:numel(PDE.BC{ii}.term)
        term_jj = PDE.BC{ii}.term{jj};
        term_val_jj = apply_term(PDE,term_jj,has_vars_Lcomp,xp,u_pde,w_pde);
        BC_PDE(rr_indcs) = BC_PDE(rr_indcs) + term_val_jj;
    end
    rr = rr_indcs(end);
end
err_BCs = cleanpoly(BC_PDE,tol);    % Boundary functions should be equal to 0.
if any(any(err_BCs.coeff))
    error('The observed fundamental state does not match the actual fundamental state for the constructed PIE')
end




function term_val = apply_term(PDE,term_jj,has_vars_Lcomp,xp,u,w)

vars = PDE.vars;
dom = PDE.dom;
nvars = size(vars,1);

if isfield(term_jj,'w')
    % Extract the desired input component.
    Robj = term_jj.w;
    nnw_arr = cumsum([0;PDE.w_tab(:,2)]);
    Robj_indcs = (nnw_arr(Robj)+1:nnw_arr(Robj+1));
    term_val = w(Robj_indcs);
    
    % Determine which variables the input component depends on.
    has_vars_Rcomp = logical(PDE.u_tab(Robj,3:2+nvars));
    var1_Rcomp = vars(has_vars_Rcomp,1);
    var2_Rcomp = vars(has_vars_Rcomp,2);
    
elseif isfield(term_jj,'u')
    % Extract the desired input component.
    Robj = term_jj.u;
    nnu_arr = cumsum([0;PDE.u_tab(:,2)]);
    Robj_indcs = (nnu_arr(Robj)+1:nnu_arr(Robj+1));
    term_val = u(Robj_indcs);
    
    % Determine which variables the input component depends on.
    has_vars_Rcomp = logical(PDE.u_tab(Robj,3:2+nvars));
    var1_Rcomp = vars(has_vars_Rcomp,1);
    var2_Rcomp = vars(has_vars_Rcomp,2);
    
else
    % Extract the desired state component.
    Robj = term_jj.x;
    nnx_arr = cumsum([0;PDE.x_tab(:,2)]);
    Robj_indcs = (nnx_arr(Robj)+1:nnx_arr(Robj+1));
    term_val = xp(Robj_indcs);
    
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
        term_val = subs(term_val,var1_kk,locval(kk));        
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
use_dummy_vars = has_vars_Lcomp & use_int_full;
if ~nvars==0 && any(use_dummy_vars)
    term_val = subs(term_val,vars(use_dummy_vars,1),vars(use_dummy_vars,2));
    Cval = subs(term_jj.C,vars(~use_dummy_vars,2),vars(~use_dummy_vars,1));
else
    Cval = term_jj.C;
end
term_val = Cval * term_val;
% Apply the desired integrals.
for kk=1:size(term_jj.I,1)
    if ~isempty(term_jj.I{kk})
        if use_dummy_vars(kk)
            var_kk = var2_Rcomp(kk);
        else
            var_kk = var1_Rcomp(kk);
        end
        term_val = int(term_val,var_kk,term_jj.I{kk}(1),term_jj.I{kk}(2));
    end    
end
% Finally, check that the term indeed does not depend on any illegal
% variables.
if isa(term_val,'polynomial') && any(~ismember(term_val.varname,vars(has_vars_Lcomp,1).varname))
    error('The term appears to depend on a variable that the LHS does not depend on... Is something wrong in the initializer?')
end
    
end


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