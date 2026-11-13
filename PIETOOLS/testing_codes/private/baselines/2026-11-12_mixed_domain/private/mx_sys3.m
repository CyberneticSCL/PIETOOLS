function S = mx_sys3(name,p)
% S = MX_SYS3(NAME,P) builds the subsystems D1-D6 that localize the
% light-tier failures of S17 (transport edge state into a state of order 0
% in s1) and S18 (library PDE 43, the delay a 2-D transport state), and
% passes every other name to MX_SYS2.
% INPUT
% - name:   char, the system name;
% - p:      (optional) 'struct' of parameters.
% OUTPUT
% - S:      'struct' as in MX_SYS; for D1-D6 only the field PDE.

% Initial coding MMP, 11/12/2026.
if nargin<2, p = struct(); end
pvar s1 s2 t
clear stateNameGenerator
S = struct('name',name,'PDE',[],'par',p);
switch name
    case 'D1_transport_1d'
        % S17's x1 alone: x_t = -x_s, x(0) = 0; nilpotent, so stable
        x1 = pde_var(1,s1,[0,1]);
        S.PDE = [diff(x1,t)==-diff(x1,s1); subs(x1,s1,0)==0];
    case 'D2_s17_x2_alone'
        % S17's x2 alone: order 0 in s1, heat in s2, Dirichlet in s2;
        % spectrum r2 - n^2 pi^2 (each infinitely repeated)
        r2 = getp(p,'r2',0);
        x2 = pde_var(1,[s1;s2],[0,1;0,1],[0;2]);
        S.PDE = [diff(x2,t)==diff(x2,s2,2)+r2*x2; subs(x2,s2,0)==0; subs(x2,s2,1)==0];
    case 'D3_s17_decoupled'
        % S17 with zero edge data: D1 and D2 side by side
        r2 = getp(p,'r2',0);
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1],[0;2]);
        S.PDE = [diff(x1,t)==-diff(x1,s1); diff(x2,t)==diff(x2,s2,2)+r2*x2;
                 subs(x1,s1,0)==0; subs(x2,s2,0)==0; subs(x2,s2,1)==0];
    case 'D4_s17_heat_source'
        % S17 with a Dirichlet 1-D heat state as the edge data: S3 with a
        % target of order 0 in s1; margin r2* = pi^2
        r2 = getp(p,'r2',0);
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1],[0;2]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s2,2)+r2*x2;
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
    case 'D5_dde_scalar'
        % y' = a0 y + a1 y(t-tau), the sin(s1) mode of S18; with a0 = 0.9,
        % a1 = -1: tau* = acos(0.9)/sqrt(0.19) = 1.034726
        tau = getp(p,'tau',0.5);  a0 = getp(p,'a0',0.9);  a1 = getp(p,'a1',-1);
        y = pde_var(1);
        S.PDE = [diff(y,t)==a0*y+a1*subs(y,t,t-tau)];
    case 'D6_heat1d_reaction'
        % S18's x1 without the delay: x_t = x_ss + a0 x on [0,pi], Dirichlet;
        % rightmost eigenvalue a0 - 1
        a0 = getp(p,'a0',0.9);
        x1 = pde_var(1,s1,[0,pi]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2)+a0*x1; subs(x1,s1,0)==0; subs(x1,s1,pi)==0];
    otherwise
        S = mx_sys2(name,p);
end
end

function v = getp(p,f,d)
if isfield(p,f), v = p.(f); else, v = d; end
end
