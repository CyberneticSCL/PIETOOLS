function S = mx_sys2(name,p)
% S = MX_SYS2(NAME,P) builds the control C1 (one 2-D heat state, Dirichlet
% on every edge, scalar w as a uniform source, z = int int x) and passes
% every other name to MX_SYS.
% INPUT
% - name:   char, the system name;
% - p:      (optional) 'struct' of parameters.
% OUTPUT
% - S:      'struct' as in MX_SYS; for C1 only the field PDE.

% Initial coding MMP, 11/12/2026.
if nargin<2, p = struct(); end
switch name
    case 'C1_heat2d_scalar_w'
        % Control for the H-infinity and H2 executives: one 2-D heat state,
        % Dirichlet on every edge, scalar w as a uniform source, z = int int x.
        % No lower-dimensional state, no auxiliary state (K empty).
        pvar s1 s2 t
        clear stateNameGenerator
        x = pde_var(1,[s1;s2],[0,1;0,1]);  w = pde_var('in',1);  z = pde_var('out',1);
        S = struct('name',name,'PDE',[],'par',p);
        S.PDE = [diff(x,t)==diff(x,s1,2)+diff(x,s2,2)+w;
                 z==int(x,[s1;s2],[0,1;0,1]);
                 subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0];
    otherwise
        S = mx_sys(name,p);
end
end
