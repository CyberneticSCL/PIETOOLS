% t5_opnorm.m -- validate opnorm_pi on opvars with known norms, then on a real
% PIETOOLS operator using structural invariants (self-adjointness, positivity)
% that the scalar tests cannot check.
cuadmm_path;
pvar s th t

mk = @(R0,R1,R2,P,Q1,Q2) mkop(R0,R1,R2,P,Q1,Q2,s,th);

fprintf('T5 --- identity on L2: norm 1, M=I\n');
M = cuadmm_private('opnorm_pi',mk(polynomial(1),[],[],[],[],[]),12);
fprintf('T5 ||M||=%.12f  ||M-I||=%.3e\n', norm(M), norm(M-eye(size(M))));

fprintf('T5 --- Volterra on L2: 2/pi=%.12f\n',2/pi);
M = cuadmm_private('opnorm_pi',mk(polynomial(0),polynomial(1),[],[],[],[]),24);
fprintf('T5 ||M||=%.12f err=%+.3e\n', norm(M), norm(M)-2/pi);

fprintf('T5 --- block diag(P=2, I on L2): norm 2\n');
M = cuadmm_private('opnorm_pi',mk(polynomial(1),[],[],2,polynomial(0),polynomial(0)),12);
fprintf('T5 ||M||=%.12f err=%+.3e  size=%d\n', norm(M), norm(M)-2, size(M,1));

fprintf('T5 --- REAL PIETOOLS operator: Top''*Top for the Dirichlet heat PIE\n');
x = pde_var('state',1,s,[0,1]);
PIE = convert([diff(x,t,1)==diff(x,s,2)+2*x; subs(x,s,0)==0; subs(x,s,1)==0]);
G = PIE.T'*PIE.T;
for D = [8 16 24]
    M = cuadmm_private('opnorm_pi',G,D);
    asym = norm(M-M','fro')/max(norm(M,'fro'),eps);
    ev   = eig((M+M')/2);
    fprintf('T5 D=%-3d ||M||=%.10f  rel_asym=%.3e  min_eig=%+.3e  (self-adjoint & PSD expected)\n', ...
            D, norm(M), asym, min(ev));
end
fprintf('T5DONE\n');

function Op = mkop(R0,R1,R2,P,Q1,Q2,s,th)
opvar Op;
Op.I = [0,1];  Op.var1 = s;  Op.var2 = th;
n0 = 0; if ~isempty(P), n0 = size(P,1); end
Op.dim = [n0 n0; 1 1];
if ~isempty(P),  Op.P  = P;  end
if ~isempty(Q1), Op.Q1 = Q1; end
if ~isempty(Q2), Op.Q2 = Q2; end
Op.R.R0 = R0;
if ~isempty(R1), Op.R.R1 = R1; end
if ~isempty(R2), Op.R.R2 = R2; end
end
