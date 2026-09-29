function out = stability_1d_sop(plant,level)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = STABILITY_1D_SOP(PLANT,LEVEL) solves the 1-D PIE stability LPI of
% executives/PIETOOLS_PIE2PDEstability.m (Q form, the default 'stability'
% executive) on the container path, extracts the Lyapunov operator with
% lpigetsol_sop and certifies it:
%
%   P = P1 (+P2) + eppos2*T'*T,  T'*Q - P = 0 (Q free),
%   D = A'*Q + Q'*A + epneg*P,   D + N = 0,  N >= 0,
%
% which certifies stability of the PIE T z_t = A z. The program is
% cx_PIE2PDEstability's (sopvar/Testfolder/sdopvar/claude_tests/cx_exec),
% which also returns its operators.
%
% CHECKS (each asserted; error id stability_1d_sop:check):
% (1) solve certified: MOSEK numerr 0, pinf 0; rel_b <= 1e-6 on the
%     solinfo.RRx rows; Gram blocks PSD (min eig >= -1e-8 max);
% (2) the extracted constraint operators reproduce the SDP's rows: the
%     rows of each relation, rebuilt by lpi_eq_sop on an empty copy of the
%     program and evaluated at solinfo.RRx, have norm equal to the
%     coefficient norm of the extracted T'Q - P (every coefficient is a
%     row), and between 1 and sqrt(2) times it for D + N ('symmetric');
%     the two row sets together are the program's rows (norms add up);
% (3) semantically, by action on random test functions (opcheck_sop):
%     T'Q = P in weak form <T y, Q x> = <y, P x>, and D = A'Q + Q'A (+
%     epneg P) in weak form <y, D x> = <A y, Q x> + <Q y, A x>, neither
%     computing an adjoint or a composition; D + N = 0; each to 1e-5
%     relative; P >= 0 and N >= 0 on the test functions;
% (4) control: solinfo.x in RRx's place breaks (3) by more than 1e3 times.
%
% INPUT
% - plant: cx_plant argument cell (default {'rd',0.5}, u_t = u_ss +
%          0.5 pi^2 u, Dirichlet: stable);
% - level: lpisettings name (default 'light').
% OUTPUT struct: rel_b, psd_relmin, the ratios of (2), the residuals and
%   margins of (3), the control, and the times (s).
%
% Initial coding MMP, 09/29/2026. Tier 1 end-to-end example E2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(plant),  plant = {'rd',0.5};     end
if nargin<2 || isempty(level),  level = 'light';        end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
chk = @(c,varargin) assert(c,'stability_1d_sop:check',varargin{:});
st = cx_settings(level,'mosek');
PIE = cx_plant(plant{:});

t0 = tic;
[prog,ops] = cx_PIE2PDEstability(PIE,st);
t_build = toc(t0);
if ~isfield(ops,'P'),   error('stability_1d_sop:plant','A 1-D PIE is required.'),   end
t0 = tic;
prog = lpisolve(prog,st.sos_opts);
t_solve = toc(t0);
info = prog.solinfo.info;
[rb,~,prel] = cx_resid(prog);
fprintf('E2 %s/%s: numerr %d, pinf %d, rel_b %.2e, psd_relmin %.2e; %d decision variables; build %.1f s, solve %.1f s\n',...
    strjoin(cellfun(@num2str,plant,'uni',0),','),level,info.numerr,info.pinf,rb,prel,...
    numel(prog.decvartable),t_build,t_solve);
chk(info.numerr==0 && info.pinf==0 && rb<=1e-6 && prel>=-1e-8,...
    '(1) solve not certified: numerr %d, pinf %d, rel_b %.3g, psd_relmin %.3g',info.numerr,info.pinf,rb,prel)

% % % Extraction.
t0 = tic;
Ps = lpigetsol_sop(prog,ops.P);     Qs = lpigetsol_sop(prog,ops.Q);
Ds = lpigetsol_sop(prog,ops.D);     Ns = lpigetsol_sop(prog,ops.N);
E1 = ops.T'*ops.Q - ops.P;          E2 = ops.D + ops.N;     % as imposed
E1s = lpigetsol_sop(prog,E1);       E2s = lpigetsol_sop(prog,E2);
t_extract = toc(t0);
chk(all(cellfun(@(X) isa(X,'copvar'),{Ps,Qs,Ds,Ns,E1s,E2s})),'(1) extracted operators not fixed')

% (2) extraction against the SDP rows, per relation.
r1 = rows_at(prog,E1,{});   r2 = rows_at(prog,E2,{'symmetric'});
rall = row_res(prog,1:prog.expr.num);
c1 = opcheck_sop('coef',E1s);   c2 = opcheck_sop('coef',E2s);
fprintf('   ||coef E(sol)|| / ||r_rows||: T''Q-P %.6f (rows %.2e), D+N %.4f (rows %.2e); rows %.3e = %.3e (all)\n',...
    c1/r1,r1,c2/r2,r2,hypot(r1,r2),rall);
chk(abs(hypot(r1,r2)/rall-1)<=1e-6,'(2) the rebuilt rows are not the program''s rows (%.6g vs %.6g)',hypot(r1,r2),rall)
chk(abs(c1/r1-1)<=1e-6,'(2) T''Q-P: coefficient norm %.6g vs row residual %.6g',c1,r1)
chk(c2>=r2*(1-1e-6) && c2<=sqrt(2)*r2*(1+1e-6),'(2) D+N: coefficient norm %.6g vs row residual %.6g',c2,r2)

% (3) semantic certificate.
o = struct('nt',6,'nqi',12,'nqo',14);
T = ops.T;  A = ops.A;
w1 = opcheck_sop('weak',{1,T,Qs; -1,[],Ps},Ps,o);                  % <Ty,Qx> - <y,Px>
wD = {1,A,Qs; 1,Qs,A; -1,[],Ds};                                    % <Ay,Qx> + <Qy,Ax> - <y,Dx>
if st.epneg~=0,     wD(end+1,:) = {st.epneg,[],Ps};     end
wd = opcheck_sop('weak',wD,Ds,o);
rN = opcheck_sop('res',Ds,-Ns,o);                                   % D = -N
pP = opcheck_sop('psd',Ps,o);   pN = opcheck_sop('psd',Ns,o);
fprintf('   by action: T''Q=P weak %.2e | D=A''Q+Q''A weak %.2e | D+N %.2e | min <x,Px> %.2e, <x,Nx> %.2e\n',...
    w1,wd,rN,pP,pN);
chk(w1<=1e-5 && rN<=1e-5,'(3) the extracted operators violate the equalities (%.3g, %.3g)',w1,rN)
chk(wd<=1e-10,'(3) getsol(D) differs from A''getsol(Q)+getsol(Q)''A in weak form (%.3g)',wd)
chk(pP>=-1e-8 && pN>=-1e-8,'(3) P or N not positive on test functions (%.3g, %.3g)',pP,pN)

% (4) control: the solver's cone vector in RRx's place.
bad = prog;     nd = numel(prog.decvartable);
xx = zeros(nd,1);   nx = min(nd,numel(prog.solinfo.x));     xx(1:nx) = prog.solinfo.x(1:nx);
bad.solinfo.RRx = xx;
w1x = opcheck_sop('weak',{1,T,lpigetsol_sop(bad,ops.Q); -1,[],lpigetsol_sop(bad,ops.P)},Ps,o);
rNx = opcheck_sop('res',lpigetsol_sop(bad,ops.D),-lpigetsol_sop(bad,ops.N),o);
fprintf('   control solinfo.x: T''Q=P weak %.2e, D+N %.2e (||x - RRx|| = %.2e)\n',w1x,rNx,...
    norm(xx-prog.solinfo.RRx(1:nd)));
chk(max(w1x,rNx)>1e3*max([w1,rN,1e-12]),'(4) reading solinfo.x is not detected')

out = struct('rel_b',rb,'psd_relmin',prel,'ndec',nd,'coef_rows',[c1/r1, c2/r2],'rows_res',[r1 r2],...
    'weak_TQP',w1,'weak_D',wd,'res_DN',rN,'psd',[pP pN],'control',[w1x rNx],...
    't',struct('build',t_build,'solve',t_solve,'extract',t_extract));
fprintf('E2 PASS\n');
end


% ------------------------------------------------------------------------
function r = rows_at(prog,E,o)
% Norm of the rows lpi_eq_sop(.,E,o{:}) imposes, at solinfo.RRx: E imposed
% on a copy of PROG with no expressions, then evaluated as row_res does.
pe = prog;
pe.expr.num = 0;
fn = fieldnames(pe.expr);
for k = 1:numel(fn),    if iscell(pe.expr.(fn{k})),    pe.expr.(fn{k}) = {};   end,    end
pe = lpi_eq_sop(pe,E,o{:});
r = row_res(pe,1:pe.expr.num);
end

function r = row_res(prog,rows)
% ||At_i' x - b_i|| over the expressions ROWS, x = solinfo.RRx
% (decvartable order, the order of At rows).
x = prog.solinfo.RRx(:);    r2 = 0;
for i = rows
    At = prog.expr.At{i};   b = prog.expr.b{i};
    r2 = r2 + sum((At'*x(1:size(At,1)) - b).^2);
end
r = sqrt(full(r2));
end
