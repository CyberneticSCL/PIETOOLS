function [sol,M] = bl_b_io2(exec,setname)
% 2-D plant with disturbance and regulated output, transcribed from the shipped
% Ex_2D_ReactionDiffusion_DDDD. THIS ROW HAS A GROUND TRUTH: that example file
% carries a closed-form L2 gain, recomputed below, so the returned gamma can be
% checked against an analytic value instead of only against another solver.
%
% Output count comes from nargout() rather than being hardcoded: the 2-D
% executives are not uniform (Hinf_gain_2D returns 4, H2_norm_2D_c returns 3),
% and asking for more outputs than a function declares aborts the case AFTER the
% solve has already been paid for.
%
% M is assembled field by field. struct() with a cell value expands into a
% struct ARRAY, which is how the first version of this file threw away four
% completed 2-D solves at m=14005 apiece.
clear stateNameGenerator
pvar s1 s2 t
nu = 1;  rr = 15;
x = pde_var(1,[s1;s2],[0,1;0,1]);
w = pde_var('in',1,[s1;s2],[0,1;0,1]);
z = pde_var('out',1);
PDE = [diff(x,t,1)==rr*x+nu*(diff(x,s1,2)+diff(x,s2,2))+w;
       z==int(x,[s1;s2],[0,1;0,1]);
       subs(x,s1,0)==0; subs(x,s1,1)==0;
       subs(x,s2,0)==0; subs(x,s2,1)==0];
PIE = initialize(convert(PDE));
st  = bl_settings(setname);
fn  = ['PIETOOLS_' exec];
nout = abs(nargout(fn));
out = cell(1,nout);
[out{1:nout}] = feval(fn,PIE,st);
sol = out{1};

M = struct();
M.Top = PIE.T;
M.Aop = PIE.A;
M.aux = {out(2:end)};
for k = 2:nout
    if isnumeric(out{k}) && isscalar(out{k}), M.gam = double(out{k}); break; end
end
% closed-form L2 gain, from the shipped example file's own tail
Mm = 10;  Nn = 10;
mu_mn   = nu*pi^2*((2*(1:Mm)'-1).^2 + (2*(1:Nn)-1).^2);
mn_fact = ((2*(1:Mm)'-1)).*((2*(1:Nn)-1));
fct     = sqrt(sum((1./((mu_mn-rr).*mn_fact)).^2,'all'));
M.gam_exact = (8/(pi^2))*fct;
end
