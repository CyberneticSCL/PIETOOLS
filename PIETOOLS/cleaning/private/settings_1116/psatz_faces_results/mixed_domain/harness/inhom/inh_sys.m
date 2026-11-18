function PIE = inh_sys(name,lam)
% PIE = INH_SYS(NAME,LAM) the variable-coefficient test systems on [0,1]^2
% (Dirichlet on all sides), NAME 'A', 'B', 'C', and the mixed cascade 'M'
% (LAM = [r1 r2]). SCRATCH, 11/16/2026.
clear stateNameGenerator
pvar s1 s2 t
switch name
    case 'A'    % x_t = (1+s1^2) x_s1s1 + (1+s2) x_s2s2 + lam x
        x = pde_var(1,[s1;s2],[0,1;0,1]);
        PDE = [diff(x,t)==(1+s1^2)*diff(x,s1,2)+(1+s2)*diff(x,s2,2)+lam*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0];
    case 'B'    % x_t = x_s1s1 + x_s2s2 + (lam - 3 s1 s2) x
        x = pde_var(1,[s1;s2],[0,1;0,1]);
        PDE = [diff(x,t)==diff(x,s1,2)+diff(x,s2,2)+(lam-3*s1*s2)*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0];
    case 'C'    % x_t = (1+s1) x_s1s1 + x_s2s2 + s2 x_s1 + lam x
        x = pde_var(1,[s1;s2],[0,1;0,1]);
        PDE = [diff(x,t)==(1+s1)*diff(x,s1,2)+diff(x,s2,2)+s2*diff(x,s1)+lam*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0];
    case 'M'    % S3 cascade, variable coefficients: x1 on the edge s2 = 0
        r1 = lam(1);    r2 = lam(2);
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        PDE = [diff(x1,t)==(1+s1)*diff(x1,s1,2)+r1*x1;
               diff(x2,t)==diff(x2,s1,2)+(1+s2^2)*diff(x2,s2,2)+r2*x2;
               subs(x1,s1,0)==0; subs(x1,s1,1)==0;
               subs(x2,s1,0)==0; subs(x2,s1,1)==0;
               subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
end
evalc('PIE = initialize(convert(PDE,{''pie'',''silent''}));');
end
