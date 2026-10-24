% Test_polyopvar_subs_ode_rows.m -- substitution of the PDE states u in the
% monomials Z(u) by the rows of Tx=T*x, for PIEs with ODE and PDE states
% (rows of different spatial dependence: the ODE rows depend on no variable
% and the PDE rows on s), also with vector-valued states:
%   - library example 9:  ODE state (size 1) and PDE state (size 1);
%   - library example 8:  ODE state (size 4) and PDE state (size 2);
% for monomials of degree 1 (all rows) and of degree 1 in a subset of the
% states, by name (@polyopvar/subs.m, subs_var with val_rows), and of the
% inner product <Tx,Tx> of the rows (sosprogramming/innerprod.m: ODE rows
% contribute u'*u, PDE rows the integral of u'*u).
%
% Verification is independent of subs: the PDE states u0=Tx(x0) and the
% monomials Z(u0) are evaluated by apply_op (eval_polyopvar), for random
% PIE states x0, and compared with the substituted polynomial Z(Tx)
% evaluated at x0.
%
% DJ (CC), 10/21/2026: Initial coding;

tol = 1e-8;
sv = polynomial({'s'});     sv = sv(1);
ncases = 0;
for ex = [9,8]
    rng(500+ex);
    PDE = examples_PDE_library_PIETOOLS_nouinput(ex);
    evalc('PIE = convert_PIETOOLS_PDE_nonlinear(PDE,''silent'');');
    sz = PIE.x_tab(:,2);
    nst = numel(sz);
    % Random PIE states: constants for ODE states, polynomials in s for PDE states
    X0 = cell(nst,1);
    for i=1:nst
        if PIE.x_tab(i,3)==0
            X0{i} = polynomial(randn(sz(i),1));
        else
            X0{i} = polynomial(randn(sz(i),1)) + randn(sz(i),1)*sv + randn(sz(i),1)*sv^2;
        end
    end
    Tx = PIE.T*PIE.x;
    u0 = eval_polyopvar(Tx,X0);
    row0 = cumsum([0;sz(:)]);
    U0 = cell(nst,1);
    for i=1:nst
        U0{i} = u0(row0(i)+1:row0(i+1));
    end
    xu = polyopvar(PIE.varname_pde,sz,{'s'},PIE.dom,logical(PIE.x_tab(:,3)));
    Zu = dmonomials(xu,1);
    % (a) all states
    Zref = eval_polyopvar(Zu,U0);
    Znew = subs(Zu,Zu.varname,Tx);
    Y = eval_polyopvar(Znew,state_values(Znew,PIE,X0,U0));
    check_equal(Y,Zref,tol,sprintf('example %d, all states',ex));
    ncases = ncases+1;
    % (b) substitution of a subset of the states: the others remain
    for i=1:nst
        vals = U0;      vals(i) = [];
        names = PIE.varname_pde(i);
        ridx = row0(i)+1:row0(i+1);
        % substitute only state i by its row(s) of Tx: the remaining PDE states stay
        Tx_i = Tx(ridx,:);
        Zi = subs(Zu,names,Tx_i);
        % evaluate: remaining PDE states are given by their values u0, the PIE states by x0
        Yi = eval_polyopvar(Zi,state_values(Zi,PIE,X0,U0));
        check_equal(Yi,Zref,tol,sprintf('example %d, only state %d substituted',ex,i));
        ncases = ncases+1;
    end
    % (c) inner product of the rows of Tx: ODE rows u'*u, PDE rows int u'*u ds
    V = innerprod(Tx,Tx);
    Vref = polynomial(0);
    for i=1:nst
        q = U0{i}.'*U0{i};
        if PIE.x_tab(i,3)~=0
            q = int(q,sv,0,1);
        end
        Vref = Vref + q;
    end
    vals = struct();
    for i=1:nst,     vals.(PIE.varname_pie{i}) = X0{i};     end
    Y = eval_V(V,vals);
    check_equal(Y,Vref,tol,sprintf('example %d, innerprod',ex));
    ncases = ncases+1;
end

fprintf('Test_polyopvar_subs_ode_rows passed (%d cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function xv = state_values(P,PIE,X0,U0)
% Values of the states of P, by name: PIE states x0, PDE states u0 (states of P that are not used in P
% are given a value of the right size)
xv = cell(numel(P.varname),1);
for k=1:numel(xv)
    [tf,jx] = ismember(P.varname{k},PIE.varname_pie);
    if tf
        xv{k} = X0{jx};
    else
        [~,jx] = ismember(P.varname{k},PIE.varname_pde);
        xv{k} = U0{jx};
    end
end
end

function Y = eval_V(V,vals)
% Evaluate a 1x1 polyopvar whose coefficients are integral operators acting on d-fold tensor products
% of the states (as in Test_piesos_LF): the dummy variables of s_k are numbered over the factors that
% depend on s_k, as s_k_j, if the monomial has more than 1 factor
vnames = V.vars.in;
dep_state = double(V.depmat.in)>0;
Y = 0;
Zf = full(V.Z);
for k=1:size(V.Z,1)
    op = V.C.ops{1,k};
    if isempty(op),     continue,   end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    slot = zeros(1,numel(vnames));
    v = 1;
    for j=1:d
        xj = vals.(V.varname{idcs(j)});
        for mm=1:numel(vnames)
            if dep_state(idcs(j),mm)
                slot(mm) = slot(mm)+1;
                if d>1
                    isv = strcmp(xj.varname,vnames{mm});
                    xj.varname(isv) = {[vnames{mm},'_',num2str(slot(mm))]};
                end
            end
        end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
end
function check_equal(Y,Yref,tol,label)
if ~isequal(size(Y),size(Yref))
    error("Test_polyopvar_subs_ode_rows failed (%s): size [%s] differs from reference [%s].",label,num2str(size(Y)),num2str(size(Yref)))
end
Yc = polynomial(Yref).coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("Test_polyopvar_subs_ode_rows (%s): reference is zero; test would be vacuous.",label)
end
err = polynomial(Y-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol*max(1,max(abs(full(Yc(:)))))
    error("Test_polyopvar_subs_ode_rows failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end
