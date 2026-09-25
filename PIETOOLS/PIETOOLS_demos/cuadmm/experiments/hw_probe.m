% hw_probe.m -- localise the Hinf BUILD wall.
%
% bl_b_io1('Hinf_gain','heavy',2) reached 19.7 GB resident with no solve in
% sight and had to be killed. n=1 heavy is trivial (hinf_rd1_hv: m=267, build
% and solve in well under a second), so something between n=1 and n=2 is
% superlinear by a lot.
%
% Three candidate drivers, separated before any profiling:
%   n            the state/disturbance/output count
%   settings     light vs heavy (degree of the positive operators)
%   disturbance  DISTRIBUTED w (pde_var('in',n,s,dom)) as bl_b_io1 declares it,
%                versus a finite-dimensional w. A distributed disturbance makes
%                Bw a PI operator rather than a matrix, which changes the block
%                structure of the KYP form, and the earlier hinf_bisect harness
%                used a finite w -- so this has never been varied deliberately.
%
% BUILD ONLY: the executive is not called, because it solves internally and the
% question is whether the wall is before the solver. The build is mirrored from
% PIETOOLS_Hinf_gain.m with gamma pinned to a constant, which removes the
% objective and the lpi_ineq but leaves every operator construction identical.
%
% Guarded: memory is read after each stage and the case is abandoned if MATLAB
% passes CAP_GB, so a probe cannot repeat the thrash it is diagnosing.

cuadmm_path;
CAP_GB = 12;
fprintf('HW n|set|wdist|stage|t_s|mem_GB|m|note\n');
for wdist = [true false]
  for setname = {'light','heavy'}
    for n = [1 2 3]
      try
        run_one(n,setname{1},wdist,CAP_GB);
      catch ME
        fprintf('HW %d|%s|%d|ABORT|-|-|-|%s\n',n,setname{1},wdist, ...
                regexprep(ME.message,'[\t\r\n]+',' '));
      end
    end
  end
end
fprintf('HWDONE\n');

function run_one(n,setname,wdist,CAP_GB)
clear stateNameGenerator
pvar s t
x = pde_var(n,s,[0,1]);
if wdist, w = pde_var('in',n,s,[0,1]); else, w = pde_var('in',n); end
z = pde_var('out',n);
PDE = [diff(x,t,1)==diff(x,s,2)+0.5*pi^2*x+w;
       z==int(x,s,[0,1]);
       subs(x,s,0)==0;  subs(x,s,1)==0];
PIE = initialize(convert(PDE));
st = lpisettings(setname);
gam = 1;                                  % pinned: same operators, no objective

Top=PIE.T; Aop=PIE.A; Bwop=PIE.Bw; Czop=PIE.Cz; Dzwop=PIE.Dzw;
mark(n,setname,wdist,'convert',tic,NaN,'');

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
t0=tic; [prog,R1op] = poslpivar(prog,Top.dim,st.dd1,st.options1);
if st.override1~=1
    [prog,P2op] = poslpivar(prog,Top.dim,st.dd12,st.options12);
    Rop = R1op+P2op;
else, Rop = R1op; end
g=mark(n,setname,wdist,'poslpivar_R',t0,NaN,''); if g>CAP_GB, error('hw:cap','over cap at poslpivar_R'); end

t0=tic; Qdeg = get_lpivar_degs(Rop,Top);
[prog,Qop] = lpivar(prog,Top.dim,Qdeg);
g=mark(n,setname,wdist,'lpivar_Q',t0,NaN,''); if g>CAP_GB, error('hw:cap','over cap at lpivar_Q'); end

t0=tic; prog = lpi_eq(prog,Top'*Qop-Rop);
g=mark(n,setname,wdist,'lpi_eq_TQ',t0,NaN,''); if g>CAP_GB, error('hw:cap','over cap at lpi_eq_TQ'); end

t0=tic;
Iw = mat2opvar(eye(size(Bwop,2)),Bwop.dim(:,2),PIE.vars,PIE.dom);
Iz = mat2opvar(eye(size(Czop,1)),Czop.dim(:,1),PIE.vars,PIE.dom);
Dop = [-gam*Iw,    Dzwop',   Bwop'*Qop;
        Dzwop,     -gam*Iz,  Czop;
        Qop'*Bwop, Czop',    Aop'*Qop+Qop'*Aop];
g=mark(n,setname,wdist,'Dop_compose',t0,NaN,''); if g>CAP_GB, error('hw:cap','over cap at Dop_compose'); end

t0=tic; [prog,De1op] = poslpivar(prog,Dop.dim,st.dd2,st.options2);
if st.override2~=1
    [prog,De2op] = poslpivar(prog,Dop.dim,st.dd3,st.options3);
    Deop = De1op+De2op;
else, Deop = De1op; end
g=mark(n,setname,wdist,'poslpivar_De',t0,NaN,''); if g>CAP_GB, error('hw:cap','over cap at poslpivar_De'); end

t0=tic; prog = lpi_eq(prog,Deop+Dop,'symmetric');
S = cuadmm_private('sdpshape',prog);
mark(n,setname,wdist,'lpi_eq_neg',t0,S.m,sprintf('Ks=%s nnzAt=%d',mat2str(S.Ks),S.nnzAt));
end

function g = mark(n,setname,wdist,stage,t0,m,note)
M = memory; g = M.MemUsedMATLAB/2^30;
fprintf('HW %d|%s|%d|%s|%.1f|%.2f|%s|%s\n', n,setname,wdist,stage,toc(t0),g, ...
        mstr(m),note);
end

function s = mstr(m)
if isnan(m), s = ''; else, s = num2str(m); end
end
