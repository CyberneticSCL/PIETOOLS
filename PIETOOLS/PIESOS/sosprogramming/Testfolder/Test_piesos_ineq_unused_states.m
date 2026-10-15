% Test_piesos_ineq_unused_states.m -- tests piesos_ineq for a functional
% with state variables that do not appear in any monomial, as returned by
% piesos_LF (names of the PDE state and of the copies of the PIE state,
% besides the PIE state).
%
% For the transport equation x_t = -x_s, x(0)=0, the Lie derivative of the
% functional V = <u,P*u> + eps*<u,u> along the PIE is nonpositive for P=0
% (dV = -eps*u(1)^2), so the program "find P>=0 such that -dV is a sum of
% squares" is feasible. Infeasibility (pinf=1) is therefore a failure.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - A psd variable (sosposmatr) is added so that SeDuMi's preprocessing
%   does not fail on small programs;
% - Feasible means pinf=dinf=0 and feasratio>=0.5.

clear stateNameGenerator
pvar s t
dom = [0,1];
u = pde_var(s,dom);
PDE = [diff(u,t)==-diff(u,s);   subs(u,s,dom(1))==0];
evalc('PIE = convert_PIETOOLS_PDE_nonlinear(PDE,''silent'');');
x = PIE.x;
nfeas = 0;

opts_list = {struct('degs',2,'exclude',[0;0;0]), struct('degs',0,'exclude',[0;1;1]), struct('degs',struct('in',0,'out',2),'exclude',[0;1;1])};
for io = 1:numel(opts_list)
    prog = piesos_program(x);
    Popts = opts_list{io};      Popts.pos = true;
    evalc('[prog,Vx,dV] = piesos_LF(prog,PIE,''pde'',1,Popts);');
    % Strict positivity: eps*<u,u>, with derivative 2*eps*<u,f(x)>
    Tx = PIE.T*x;
    Vx = Vx + 1e-2*innerprod(Tx,Tx);            %#ok<NASGU>
    dV = dV + 2e-2*innerprod(Tx,PIE.f);
    evalc('[prog,W] = piesos_ineq(prog,-dV);');
    prog = sosposmatr(prog,20);
    evalc('prog = piesos_solve(prog);');
    info = prog.solinfo.info;
    if info.pinf || info.dinf || info.feasratio<0.5
        error("Test_piesos_ineq_unused_states failed (options %d): program infeasible (pinf=%d, feasratio=%g).",io,info.pinf,info.feasratio)
    end
    nfeas = nfeas+1;
end

% The variables of dV from piesos_LF are the PIE state only
prog = piesos_program(x);
evalc('[~,~,dV] = piesos_LF(prog,PIE,''pde'',1,struct(''degs'',0,''exclude'',[0;1;1],''pos'',true));');
if ~isequal(reshape(dV.varname,[],1),reshape(x.varname,[],1))
    error("Test_piesos_ineq_unused_states failed: dV has unused state variables.")
end

fprintf('Test_piesos_ineq_unused_states passed (%d programs feasible).\n',nfeas)
