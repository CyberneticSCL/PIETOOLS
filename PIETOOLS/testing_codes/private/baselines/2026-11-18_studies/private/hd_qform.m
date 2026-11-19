function hd_qform()
% HD_QFORM, 10/28/2026: build-only probe, PDE 45 (2-D H-infinity gain), of
% the PIE2PDE (Q) form sized to contain the PDE2PDE (P) form, against the P
% form, at lpisettings_sop stripped and light:
%   P form   PIETOOLS_Hinf_gain_sop 'P' (P + eppos I, KYP slack 'like' P);
%   Q form   Q free with the degree caps of P*T (P the storage operator of
%            the settings plus eppos I); R >= 0 sized from the support of
%            T'*Q (the rule of LPI_INEQ_SOP for a target read from itself,
%            so R can equal T'*Q for every Q of that basis); T'*Q == R in
%            every coefficient, as PIETOOLS_Hinf_gain_sop imposes it; KYP
%            slack 'like' R with weight degree max(dw_Q, dw_P).
% Every program goes to the shadow Sedumi2Mosek, which records it and stops
% before MOSEK. A positive operator with sum N^2 above LIMIT is reported,
% not built.
here = fileparts(mfilename('fullpath'));
lbdir = fullfile(fileparts(here),'library_baseline');
addpath(lbdir);     lb_setup();
addpath(fullfile(here,'shadow'),'-begin');
w = which('Sedumi2Mosek','-all');
assert(startsWith(w{1},fullfile(here,'shadow')),'hd_qform: shadow Sedumi2Mosek is not first on the path');
S = load(fullfile(lbdir,'catalog.mat'));    C = S.C;
PIE = lb_build_system(C(find(strcmp({C.key},'PDE:45'),1)), [], true);
% Run 1 (10/28): P form, then the Q form with the KYP slack 'like' R; that
% slack exceeded LIMIT (Gram 2705 at stripped, 4226 at light). Run 2: the
% Q form with the KYP slack 'like' P, the P form's own rule, which contains
% the P form's slack since K is the same operator when Q = P*T.
for tier = {'stripped','light'}
    St = lpisettings_sop(tier{1},0,'','mosek');
%   record(sprintf('%s P form',tier{1}), @() PIETOOLS_Hinf_gain_sop(PIE,St,'P'));
    record(sprintf('%s Q form, slack like P',tier{1}), @() qform(PIE,St,'P'));
end
end


function record(lbl,fn)
% RECORD(LBL,FN) runs FN to the shadow and prints the recorded SDP.
global HD_FP
HD_FP = [];     t0 = tic;
try
    fn();
    fprintf('##FP %-28s returned without reaching the solver\n', lbl);
catch ME
    t = toc(t0);
    if ~strcmp(ME.identifier,'hdshadow:stop')
        fprintf('##FP %-28s ERROR after %.1f s: %s (%s)\n', lbl, t, ME.message, ME.stack(1).name);
        return
    end
    f = HD_FP;
    fprintf('##FP %-28s build %6.1f s | m %6d | n %9d | scal %9d | Kf %d | Ks [%s] | nnzA %d\n', ...
            lbl, t, f.m, f.n, f.scal, f.Kf, num2str(f.Ks), f.nnzA);
end
end


function qform(PIE,settings,like)
% QFORM(PIE,SETTINGS,LIKE) builds the Q form of the header and calls
% LPISOLVE; LIKE 'R' sizes the KYP slack from R (dw max(dw_Q,dw_P)), 'P'
% from P with dw_P, as the P form does.
LIMIT = 1e7;
st = exec_sop_settings(settings);
tl = exec_tools_sop(PIE,st);
Tm = op2copvar_sop(PIE.T);  Am = op2copvar_sop(PIE.A);
Bw = op2copvar_sop(PIE.B1); Cz = op2copvar_sop(PIE.C1);  Dzw = op2copvar_sop(PIE.D11);
Iw = tl.ident(PIE.B1,2);    Iz = tl.ident(PIE.C1,1);
[sp,dm] = copvar_space_list(Tm,'out');
% Degree caps of P*T, P the storage of the settings with its margin,
% declared in a throwaway program
prog0 = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[~,P0] = tl.storage(prog0,Tm,'out');
Pm = P0 + tl.eye(1,tl.eppos);
dP = get_lpivar_degs_sop(Pm);   dQ = get_lpivar_degs_sop(Pm*Tm);
fprintf('##DEG P: %s | P*T: %s\n', dstr(dP), dstr(dQ));
% gam >= 0 as the objective; Q free; R >= 0 with T'*Q == R; the KYP slack
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,gamv] = lpidecvar(prog,'gam');
prog = lpi_ineq(prog,gamv);     prog = lpisetobj(prog,gamv);
[prog,Qm] = lpivar_cdopvar(prog,dm,sp,tl.dom,dQ);
fprintf('##Q  free coefficients %d\n', numel(prog.decvartable)-1);
TQ = Tm'*Qm;
[prog,Rm] = pos_sized(prog,TQ,struct(),LIMIT,'R',true);
prog = lpi_eq_sop(prog,TQ-Rm);
Km = [-(gamv*Iw),   Dzw',       Bw'*Qm;
       Dzw,         -(gamv*Iz), Cz;
       Qm'*Bw,      Cz',        Am'*Qm + Qm'*Am];
if strcmp(like,'P')
    o = struct('like',Pm,'dw',st.sop.dw_P,'margin',st.sop.margin);
else
    o = struct('like',Rm,'dw',max(st.sop.dw_Q,st.sop.dw_P),'margin',st.sop.margin);
end
pos_sized(prog,-Km,o,LIMIT,'KYP slack',false);        % size check only
prog = lpi_ineq_sop(prog,-Km,o);
lpisolve(prog,st.sos_opts);
end


function [prog,Pop] = pos_sized(prog,X,opts,LIMIT,lbl,build)
% [PROG,POP] = POS_SIZED(PROG,X,OPTS,LIMIT,LBL,BUILD): the positive
% operator LPI_INEQ_SOP declares for X >= 0 with OPTS (degrees and terms of
% GET_LIFT_DEGS, bases pruned to the support of the diagonal blocks of X),
% its Gram sizes printed; declared by POSCOPVAR_DIRECT, without an
% equality, when BUILD; an error when sum N^2 > LIMIT.
[sp,dm] = copvar_space_list(X,'out');
dom = struct('vars',{X.vars},'dom',X.dom);
[deg,terms] = get_lift_degs(X,opts);
include = cell(1,numel(sp));
for k = 1:numel(sp)
    if isempty(sp{k}) || isempty(X.C{k,k}),    continue,   end
    eo = eq_opts_sopvar(X.C{k,k});      inc = eo.include;
    for t = 1:size(inc,2)               % full-integral index 4 -> 2 and 3
        r = inc(:,t)==4;
        if any(r)
            lo = inc(r,:);  up = lo;    lo(:,t) = 2;    up(:,t) = 3;
            inc = unique([inc(~r,:); lo; up],'rows');
        end
    end
    include{k} = inc;
end
popts = struct('psatz',terms.codes,'psatz_offset',terms.offsets);
if any(~cellfun(@isempty,include)),     popts.include = include;    end
eopts = popts;  eopts.estimate = true;
[~,~,~,einfo] = poscopvar_direct(prog,dm,sp,dom,deg,eopts);
nd = sum(einfo.N.^2);
fprintf('##POS %-9s Gram blocks %s | sum N^2 %d | term codes %s\n', lbl, mat2str(einfo.N), nd, mat2str(terms.codes));
if nd > LIMIT
    error('hd_qform:limit','%s: sum N^2 = %d exceeds %d; not built.',lbl,nd,LIMIT);
end
Pop = [];
if build,   [prog,Pop] = poscopvar_direct(prog,dm,sp,dom,deg,popts);   end
end


function s = dstr(d)
% S = DSTR(D): a degree-cap struct of GET_LIFT_DEGS_SOP as text.
if isstruct(d)
    s = sprintf('mult %d, int [%d %d], out %d, in %d', d.mult, d.int(1), d.int(2), d.out, d.in);
else
    s = mat2str(d);
end
end
