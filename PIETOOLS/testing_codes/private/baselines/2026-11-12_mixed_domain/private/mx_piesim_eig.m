function mx_piesim_eig()
% MX_PIESIM_EIG() repeats the steps of PIESIM.m up to
% PIESIM_discretize_all_2D for four 2-D PDEs (a control with ordinary BCs,
% S31 with an integral BC, PIESIM 2-D library examples 9 and 8) and prints,
% at N = 8 and 16, cond(Tcheb) (computed as cond(Tcheb_inv)) and the
% rightmost eigenvalues of Atotal = Tcheb_inv*Acheb, against the exact
% rightmost eigenvalue.
% INPUT
% none.
% OUTPUT
% none; one line per case and N.

% Initial coding MMP, 11/12/2026: diagnosis of D4 (MIXED_DOMAIN_REPORT.md).
syms st sx sy real
pvar s1 s2 t
cases = {};
% (1) control: Neumann at s1=0, Dirichlet at s1=1, Dirichlet in s2;
%     exact rightmost -(1/4 + 1) pi^2 = -12.337
clear stateNameGenerator
x = pde_var(1,[s1;s2],[0,1;0,1]);
cases(end+1,:) = {'control N-D', initialize([diff(x,t)==diff(x,s1,2)+diff(x,s2,2);
       subs(diff(x,s1),s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0]), ...
       struct('ic',cos(pi*sx/2)*sin(pi*sy)), -(1/4+1)*pi^2};
% (2) S31: int x ds1 = 0, x_s1(1,s2) = 0, Dirichlet in s2; exact -2 pi^2
clear stateNameGenerator
x = pde_var(1,[s1;s2],[0,1;0,1]);
cases(end+1,:) = {'S31 integral BC', initialize([diff(x,t)==diff(x,s1,2)+diff(x,s2,2);
       int(x,s1,[0,1])==0; subs(diff(x,s1),s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0]), ...
       struct('ic',cos(pi*sx)*sin(pi*sy)), -2*pi^2};
% (3) PIESIM 2-D library ex 9 and ex 8 (exact rightmost not computed here)
[PDE9,u9] = examples_pde_library_PIESIM_2D(9);
cases(end+1,:) = {'library ex 9', PDE9, u9, NaN};
[PDE8,u8] = examples_pde_library_PIESIM_2D(8);
cases(end+1,:) = {'library ex 8', PDE8, u8, NaN};
for i=1:size(cases,1)
    for N=[8 16]
        try
            opts = struct('N',N,'tf',0.1,'intScheme',1,'Norder',4,'dt',1e-5,'plot','no','ploteig','no');
            Dop = setup_dop(cases{i,2},opts,cases{i,3});
            lam = eig(Dop.Atotal);
            [~,o] = sort(real(lam),'descend');  lam = lam(o);
            fprintf('%-16s N %2d: size %4d, cond(Tcheb) %9.2e, max Re eig(Atotal) %+11.4e (next %+11.4e), exact %+9.4f\n', ...
                cases{i,1},N,size(Dop.Atotal,1),cond(full(Dop.Tcheb_inv)),real(lam(1)),real(lam(2)),cases{i,4});
        catch ME
            fprintf('%-16s N %2d: ERROR %s (%s:%d)\n',cases{i,1},N,ME.message,ME.stack(1).name,ME.stack(1).line);
        end
    end
end
end

function Dop = setup_dop(PDE,opts,uinput)
% The steps of PIESIM.m for a PDE input in 2-D, up to the discretization
[opts,uinput] = PIESIM_options_check({PDE,opts,uinput});
[PDE,uinput,opts,psize] = PIESIM_input_check(PDE,uinput,opts);
PIE = convert_PIETOOLS_PDE(PDE,'pie');
if isa(PIE,'pie_struct'), PIE = pie2piesim(PIE); end
PIE = rescalePIE_2D(PIE,[-1,1;-1,1]);
[uinput,psize] = PIESIM_initial_setup(uinput,psize,PIE,opts.type);
[Dop,~,~] = PIESIM_discretize_all_2D(PIE,uinput,psize);
end
