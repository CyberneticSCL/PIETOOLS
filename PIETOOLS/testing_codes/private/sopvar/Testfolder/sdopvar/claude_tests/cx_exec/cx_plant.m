function PIE = cx_plant(kind,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE = CX_PLANT(KIND,...) returns the PIE of one baseline plant, built by
% the SAME PDE as the baseline builders in PIETOOLS_demos/cuadmm/private
% (bl_b_stab1, bl_b_io1, bl_b_syn1, bl_b_stab2, bl_b_io2), so container
% results line up with the Mosek reference table of the baseline suite.
%
%   cx_plant('rd',frac[,n])    u_t = u_ss + frac*pi^2*u, Dirichlet; lam* = pi^2
%   cx_plant('tr',frac[,n])    u_t = -u_s - frac*u, u(0) = 0
%   cx_plant('wave',frac[,n])  2n-state damped system, Dirichlet
%   cx_plant('io1'[,n])        rd at 0.5*pi^2 with distributed w, z = int u
%   cx_plant('syn1')           io1 plus control u and sensed y, z = [int u; u]
%   cx_plant('rd2',frac)       2-D heat, Dirichlet on [0,1]^2; lam* = 2*pi^2
%   cx_plant('io2')            2-D rd at r = 15 with w, z = int u; gain known
%
% The builders are not called: they run the stock executive themselves.
% Kept in step with them by hand; a change there must be copied here.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear stateNameGenerator
switch kind
    case {'rd','tr','wave'}
        frac = varargin{1};
        n = 1;  if numel(varargin)>=2 && ~isempty(varargin{2}),  n = varargin{2};   end
        pvar s t
        switch kind
            case 'rd'
                x = pde_var(n,s,[0,1]);
                PDE = [diff(x,t,1)==diff(x,s,2)+frac*pi^2*x; subs(x,s,0)==0; subs(x,s,1)==0];
            case 'tr'
                x = pde_var(n,s,[0,1]);
                PDE = [diff(x,t,1)==-diff(x,s,1)-frac*x; subs(x,s,0)==0];
            case 'wave'
                x = pde_var(2*n,s,[0,1]);
                PDE = [diff(x,t,1)==diff(x,s,2)-frac*x; subs(x,s,0)==0; subs(x,s,1)==0];
        end
    case 'io1'
        n = 1;  if ~isempty(varargin),  n = varargin{1};    end
        pvar s t
        x = pde_var(n,s,[0,1]);
        w = pde_var('in',n,s,[0,1]);
        z = pde_var('out',n);
        PDE = [diff(x,t,1)==diff(x,s,2)+0.5*pi^2*x+w;
               z==int(x,s,[0,1]);
               subs(x,s,0)==0;  subs(x,s,1)==0];
    case 'syn1'
        pvar s t
        x = pde_var(1,s,[0,1]);
        w = pde_var('in',1,s,[0,1]);
        u = pde_var('control',1);
        y = pde_var('sense',1);
        z = pde_var('out',2);
        PDE = [diff(x,t,1)==diff(x,s,2)+0.5*pi^2*x+w+u;
               z==[int(x,s,[0,1]); u];
               y==int(x,s,[0,1]);
               subs(x,s,0)==0;  subs(x,s,1)==0];
    case 'rd2'
        frac = varargin{1};
        pvar s1 s2 t
        x = pde_var(1,[s1;s2],[0,1;0,1]);
        PDE = [diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+frac*2*pi^2*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0;
               subs(x,s2,0)==0; subs(x,s2,1)==0];
    case 'io2'
        pvar s1 s2 t
        nu = 1;  rr = 15;
        x = pde_var(1,[s1;s2],[0,1;0,1]);
        w = pde_var('in',1,[s1;s2],[0,1;0,1]);
        z = pde_var('out',1);
        PDE = [diff(x,t,1)==rr*x+nu*(diff(x,s1,2)+diff(x,s2,2))+w;
               z==int(x,[s1;s2],[0,1;0,1]);
               subs(x,s1,0)==0; subs(x,s1,1)==0;
               subs(x,s2,0)==0; subs(x,s2,1)==0];
    otherwise
        error('cx_plant:kind','Unknown plant ''%s''.',kind)
end
evalc('PIE = initialize(convert(PDE));');
end
