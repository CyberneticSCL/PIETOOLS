% PF_RUN_C3B  Phase 4 again: the io2 gain through the executive with the
% slack weight pinned at the tier's own level and one above (dw_P = dw_Q =
% 2, 3; the first run pinned 0 and 1, far below the tier default 2, and
% every configuration returned an uncertified gain). SCRATCH, MMP, 11/15/2026.
nu = 1;     rr = 15;    Mm = 10;    Nn = 10;
mu_mn   = nu*pi^2*((2*(1:Mm)'-1).^2 + (2*(1:Nn)-1).^2);
mn_fact = (2*(1:Mm)'-1).*(2*(1:Nn)-1);
gex = (8/pi^2)*sqrt(sum((1./((mu_mn-rr).*mn_fact)).^2,'all'));
fprintf('io2 closed-form gain %.7f\n',gex);
C = pf_cfg({'none','box','product','faces','facesS','corners'},2);
pf_hinf(@io2_pie,gex,C,[2 3],{'def','stock'},'light','c3b_io2');
disp('## C3B DONE');

function PIE = io2_pie()
pvar s1 s2 t
nu = 1;  rr = 15;
x = pde_var(1,[s1;s2],[0,1;0,1]);
w = pde_var('in',1,[s1;s2],[0,1;0,1]);
z = pde_var('out',1);
PDE = [diff(x,t,1)==rr*x+nu*(diff(x,s1,2)+diff(x,s2,2))+w;
       z==int(x,[s1;s2],[0,1;0,1]);
       subs(x,s1,0)==0; subs(x,s1,1)==0;
       subs(x,s2,0)==0; subs(x,s2,1)==0];
PIE = initialize(convert(PDE,{'pie','silent'}));
end
