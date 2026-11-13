function mx_piesim_eig2()
% MX_PIESIM_EIG2() is MX_PIESIM_EIG for S6 (g = 1), whose PIE has an
% auxiliary block in L2[s1] and a nonempty K; exact rightmost eigenvalue
% -pi^2 + g.
% INPUT
% none.
% OUTPUT
% none; one line per N.

% Initial coding MMP, 11/12/2026: diagnosis of D3 and D4.
syms sx sy real
S = mixed_domain_systems('S6_indomain_1d_2d',struct('g',1));
PDE = initialize(S.PDE);
for N=[8 16]
    opts = struct('N',N,'tf',0.1,'intScheme',1,'Norder',4,'dt',1e-5,'plot','no','ploteig','no');
    uinput = struct('ic',[sin(pi*sx), sin(pi*sx)*(1+cos(pi*sy)/2)]);
    try
        [opts2,uinput2] = PIESIM_options_check({PDE,opts,uinput});
        [PDE2,uinput2,opts2,psize] = PIESIM_input_check(PDE,uinput2,opts2);
        PIE = convert_PIETOOLS_PDE(PDE2,'pie');
        PIE = pie2piesim(PIE);
        PIE = rescalePIE_2D(PIE,[-1,1;-1,1]);
        [uinput2,psize] = PIESIM_initial_setup(uinput2,psize,PIE,opts2.type);
        [Dop,~,~] = PIESIM_discretize_all_2D(PIE,uinput2,psize);
        lam = eig(Dop.Atotal);  [~,o] = sort(real(lam),'descend');  lam = lam(o);
        fprintf('S6 aux/K        N %2d: size %4d, cond(Tcheb) %9.2e, max Re eig(Atotal) %+11.4e (next %+11.4e), exact %+9.4f\n', ...
            N,size(Dop.Atotal,1),cond(full(Dop.Tcheb_inv)),real(lam(1)),real(lam(2)),-pi^2+1);
    catch ME
        fprintf('S6 aux/K        N %2d: ERROR %s (%s:%d)\n',N,ME.message,ME.stack(1).name,ME.stack(1).line);
    end
end
end
