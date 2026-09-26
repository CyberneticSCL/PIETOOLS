function prog = cx_stability_dual_2D(PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_STABILITY_DUAL_2D(PIE,ST) builds, with copvar/cdopvar only, the
% LPI of executives/2D/PIETOOLS_stability_dual_2D.m (2-D direct form, dual
% PIE T' xdot = A' x):
%
%   P = P0 (+ psatz terms) + diag(eppos) >= 0,
%   Q = ((T*P)*A')' + (T*P)*A' (+ 2*epneg*(T*P)*T') = -Qe,   Qe >= 0,
%
% returned UNSOLVED (structure only). Mirrors the stock settings block (same
% as PIETOOLS_stability_2D.m:56-117), 128 (program), 139-155 (P), 166-173
% (Q), 193-229 (pruned slack, psatz terms, equality). Not mirrored: line 176
% clean_opvar (no container version; the pruning test uses the same
% tolerance), line 188 use_sosineq=1 (no lpi_ineq_2d counterpart, errors),
% 196-215 dead code (toggle=0), and the solve_val block.
%
% P lives on T's INPUT space (T*P*A'). Scope and pruning as in
% cx_stability_2D (state purely in L2[s1,s2]; cx_stability_eq_opts_2D).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(PIE,'pie_struct')
    error('The PIE should be a ''pie_struct''.')
end
PIE = initialize(PIE);
% Settings (the primal file's lines 56-117, repeated in the stock dual).
if nargin==1
    st = settings_PIETOOLS_light_2D;
    st.sos_opts.simplify = 1;
    st.eppos = 1e-2*ones(4,1);  st.epneg = 0;
elseif ~isfield(st,'is2D') || ~st.is2D
    sos_opts = st.sos_opts;     st = st.settings_2d;    st.sos_opts = sos_opts;
end
if ~isfield(st,'eppos'),    eppos = [1e-4;1e-6;1e-6;1e-6];
else,                       eppos = st.eppos;   if numel(eppos)==1,  eppos = eppos*ones(4,1);  end
end
epneg = 0;  if isfield(st,'epneg'),    epneg = st.epneg;   end
LF_use_psatz = st.LF_use_psatz;
LF_deg_psatz = psatz_cell(st.LF_deg_psatz,LF_use_psatz);
LF_opts_psatz = psatz_cell(st.LF_opts_psatz,LF_use_psatz);
if st.use_sosineq
    error('cx_stability:gap',['use_sosineq=1 needs lpi_ineq_2d, which has no container '...
          'counterpart; only the equality path (the default) is transcribed.'])
end
eq_opts = st.eq_opts;   eq_deg = st.eq_deg;     eq_use_psatz = st.eq_use_psatz;
eq_deg_psatz = psatz_cell(st.eq_deg_psatz,eq_use_psatz);
eq_opts_psatz = psatz_cell(st.eq_opts_psatz,eq_use_psatz);

Tm = opvar2d2copvar(PIE.T);     Am = opvar2d2copvar(PIE.A);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);        % line 128

% P, lines 139-147: on T's input space.
[prog,Pm] = cx_stability_pos2d(prog,Tm,'in',PIE.dom,st.LF_deg,st.LF_opts);
for j = 1:length(LF_use_psatz)
    if LF_use_psatz(j)~=0
        [prog,P2m] = cx_stability_pos2d(prog,Tm,'in',PIE.dom,LF_deg_psatz{j},LF_opts_psatz{j});
        Pm = Pm + P2m;
    end
end
% Strict positivity, lines 150-155.
if ~all(eppos==0)
    np = PIE.T.dim(:,2);
    Ip = blkdiag(eppos(1)*eye(np(1)),eppos(2)*eye(np(2)),eppos(3)*eye(np(3)),eppos(4)*eye(np(4)));
    Pm = Pm + opvar2d2copvar(opvar2d(Ip,[np,np],PIE.dom,PIE.vars));
end

% Q, lines 166-173.
TPm = Tm*Pm;    TPAm = TPm*Am';
if epneg==0
    Qm = TPAm' + TPAm;
else
    Qm = TPAm' + TPAm + 2*epneg*(TPm*Tm');
end

% Qe >= 0 with Q = -Qe, lines 193-229.
ztol = 1e-12;                                       % line 175
eq_opts = cx_stability_eq_opts_2D(Qm,eq_opts,ztol);
[prog,Qem] = cx_stability_pos2d(prog,Qm,'out',PIE.dom,eq_deg,eq_opts);
for j = 1:length(eq_use_psatz)
    if eq_use_psatz(j)~=0
        eq_opts_psatz{j}.exclude = eq_opts_psatz{j}.exclude | eq_opts.exclude;
        eq_opts_psatz{j}.sep = eq_opts_psatz{j}.sep | eq_opts.sep;
        [prog,Qe2m] = cx_stability_pos2d(prog,Qm,'out',PIE.dom,eq_deg_psatz{j},eq_opts_psatz{j});
        Qem = Qem + Qe2m;
    end
end
prog = lpi_eq_cdopvar(prog,Qem+Qm,'symmetric');
end

function outcell = psatz_cell(incell,use_psatz)
% extract_psatz_deg / extract_psatz_opts of the stock file (identical bodies).
if all(use_psatz==0),                   outcell = {};   return, end
if isa(incell,'struct'),                outcell = repmat({incell},1,length(use_psatz));
elseif numel(incell)==1,                outcell = repmat(incell(1),1,length(use_psatz));
elseif numel(incell)==length(use_psatz),outcell = incell;
elseif numel(incell)>=max(use_psatz),   outcell = incell(use_psatz);
else,   error('For each element of ''use_psatz'', a field should be defined.')
end
end
