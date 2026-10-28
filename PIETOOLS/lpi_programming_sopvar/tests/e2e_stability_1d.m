function out = e2e_stability_1d(plant,level)                                % MMP, 10/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = E2E_STABILITY_1D(PLANT,LEVEL) builds with cx_PIE2PDEstability, on
% the container path, the stability LPI of the 1-D PIE  T x_t = A x  of
% PLANT,
%
%   P = P1 (+P2) + eppos2*T'*T,  T'*Q - P = 0 (Q free),
%   D = A'*Q + Q'*A + epneg*P,   D + N = 0,  N >= 0,
%
% with P1, P2 and N the positive operators of lpisettings(LEVEL), solves
% it with MOSEK, extracts P, Q, D, N, T'*Q - P and D + N with
% lpigetsol_sop, and checks them, erroring at the first failed check:
%   (1) MOSEK numerr 0 and pinf 0; rel_b = ||At'x - b||/||b|| <= 1e-6 at
%       x = solinfo.RRx; the least eigenvalue of the Gram blocks at least
%       -1e-8 times the largest; the extracted operators are 'copvar'
%       objects;
%   (2) the rows of T'*Q - P = 0 and of D + N = 0 ('symmetric'), imposed
%       by lpi_eq_sop on a copy of the program without constraints, have
%       residual norms r1 and r2 at solinfo.RRx with sqrt(r1^2 + r2^2)
%       equal to the residual norm of all rows of the program to 1e-6
%       relative; the coefficient norm of the extracted T'*Q - P equals r1
%       to 1e-6 relative, that of D + N lies in [r2, sqrt(2) r2] (1e-6
%       slack);
%   (3) by the action on 6 random polynomial test functions x, y
%       (opcheck_sop, quadrature of the kernels): <T y, Q x> = <y, P x>
%       and D + N = 0 to 1e-5 relative; <A y, Q x> + <Q y, A x>
%       + epneg <y, P x> = <y, D x> to 1e-10 relative; <x, P x> and
%       <x, N x> at least -1e-8 (normalized);
%   (4) with the first entries of solinfo.x in place of solinfo.RRx, the
%       larger of the residuals of <T y, Q x> = <y, P x> and D + N = 0
%       exceeds 1e3 times the larger one at solinfo.RRx.
%
% INPUT
% - plant:  (optional) 'cell', the arguments of cx_plant naming a 1-D
%           plant; default {'rd',0.5}, the PIE of u_t = u_ss + 0.5 pi^2 u
%           on [0,1], u(0) = u(1) = 0;
% - level:  (optional) 'char', the name of the LPI settings (lpisettings);
%           default 'light'.
%
% OUTPUT
% - out:    'struct' with fields rel_b and psd_relmin ('double', of the
%           solve), ndec ('double', the number of decision variables),
%           coef_rows (1x2 'double', the coefficient norms over r1 and r2),
%           rows_res ([r1 r2]), weak_TQP, weak_D and res_DN ('double', the
%           residuals of (3)), psd (1x2 'double', the least <x,Px> and
%           <x,Nx>), control (1x2 'double', the residuals of (4)), t
%           ('struct' of seconds: build, solve, extract); printed lines
%           with the values.
%
% NOTES
% A 2-D plant is an error ('stability_1d_sop:plant'); a failed check
% errors with identifier 'stability_1d_sop:check'. Requires MOSEK.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - e2e_stability_1d
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 09/29/2026. Tier 1 end-to-end example E2.
% MMP, 10/21/2026: moved from examples/stability_1d_sop.m to
%                tests/e2e_stability_1d.m: an end-to-end check built on test
%                helpers (cx_exec, opcheck_sop), not an example of standard
%                usage; examples/stability_1d_sop.m is now that. Code
%                unchanged; the error identifiers keep the former name.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and headers of the local functions; code
%                unchanged. The deleted header text was the stock executive
%                the LPI reproduces and the location of the helper that
%                builds it. The tolerance of the weak form of D, 1e-10 in
%                the code, was given as 1e-5.
%                The (was) line of the former signature moved here from
%                under the function line, so that help starts with the
%                description:
% function out = stability_1d_sop(plant,level)                              % MMP, 10/21/2026 (was)

% Defaults; no warnings for non-canonical multipliers; a check asserts
% with one identifier
if nargin<1 || isempty(plant),  plant = {'rd',0.5};     end
if nargin<2 || isempty(level),  level = 'light';        end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
chk = @(c,varargin) assert(c,'stability_1d_sop:check',varargin{:});
% Settings with MOSEK; the PIE of the plant
st = cx_settings(level,'mosek');
PIE = cx_plant(plant{:});

% Build the LPI and its operators (1-D only), then solve
t0 = tic;
[prog,ops] = cx_PIE2PDEstability(PIE,st);
t_build = toc(t0);
if ~isfield(ops,'P'),   error('stability_1d_sop:plant','A 1-D PIE is required.'),   end
t0 = tic;
prog = lpisolve(prog,st.sos_opts);
t_solve = toc(t0);
% (1) solver flags, rel_b and the PSD margin
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

% Output
out = struct('rel_b',rb,'psd_relmin',prel,'ndec',nd,'coef_rows',[c1/r1, c2/r2],'rows_res',[r1 r2],...
    'weak_TQP',w1,'weak_D',wd,'res_DN',rN,'psd',[pP pN],'control',[w1x rNx],...
    't',struct('build',t_build,'solve',t_solve,'extract',t_extract));
fprintf('E2 PASS\n');
end


% ------------------------------------------------------------------------
function r = rows_at(prog,E,o)
% R = ROWS_AT(PROG,E,O) returns the residual norm, at prog.solinfo.RRx, of
% the rows that lpi_eq_sop(.,E,O{:}) imposes on a copy of PROG with no
% expressions.
% INPUT  prog: 'struct', a solved LPI program; E: 'cdopvar' object, the
%        constraint operator; o: 'cell', the options of lpi_eq_sop ({} or
%        {'symmetric'}).
% OUTPUT r: 'double'.
% Copy of PROG with no expressions
pe = prog;
pe.expr.num = 0;
fn = fieldnames(pe.expr);
for k = 1:numel(fn),    if iscell(pe.expr.(fn{k})),    pe.expr.(fn{k}) = {};   end,    end
% Impose E = 0; the residual of its rows
pe = lpi_eq_sop(pe,E,o{:});
r = row_res(pe,1:pe.expr.num);
end

function r = row_res(prog,rows)
% R = ROW_RES(PROG,ROWS) returns sqrt(sum_i ||At_i'*x - b_i||^2) over the
% expressions i in ROWS, x = prog.solinfo.RRx (decvartable order, the
% order of the rows of At_i).
% INPUT  prog: 'struct', a solved LPI program; rows: 'double' vector of
%        indices into prog.expr.At and prog.expr.b.
% OUTPUT r: 'double'.
x = prog.solinfo.RRx(:);    r2 = 0;
for i = rows
    At = prog.expr.At{i};   b = prog.expr.b{i};
    r2 = r2 + sum((At'*x(1:size(At,1)) - b).^2);
end
r = sqrt(full(r2));
end
