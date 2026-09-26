function prog = cx_PIE2PDEstability_dual(PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_PIE2PDESTABILITY_DUAL(PIE,ST) builds, with copvar/cdopvar only,
% the LPI of executives/PIETOOLS_PIE2PDEstability_dual.m (Q form, dual PIE):
%
%   P = P1 (+P2) + eppos2*T*T',     T*Q - P = 0   (Q free, shape of T'),
%   D = A*Q + Q'*A' + epneg*P = -N,  N = N1 (+N2) >= 0,
%
% returned UNSOLVED. Mirrors lines 63-70 (2-D hand-off: the stock dispatches
% to the DIRECT-form PIETOOLS_stability_dual_2D), 84-102 (settings), 108
% (program), 115-125 (P), 124-130 (Q and T*Q = P), 142 (D), 153-168
% (equality). Not mirrored: line 155 (sosineq_on=1: no container lpi_ineq,
% errors) and line 175 (no getsol).
%
% Qdeg = get_lpivar_degs(Pop,Top') (line 128) read off the cdopvar P by
% cx_stability_lpivar_degs (Top' is unused by the stock rule). Q maps T's
% output spaces to its input spaces (Top2.dim, Top2 = Top'); P lives on T's
% OUTPUT spaces (T*T', T*Q: out -> out).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(PIE,'pie_struct')
    error('The PIE should be a ''pie_struct''.')
end
PIE = initialize(PIE);
if PIE.dim==2                                       % stock lines 63-70
    if nargin<2,    prog = cx_stability_dual_2D(PIE);
    else,           prog = cx_stability_dual_2D(PIE,st);
    end
    return
end
if nargin<2                                         % stock lines 76-82
    st = settings_PIETOOLS_heavy;
    st.sos_opts.simplify = 1;
    st.eppos = 1e-4;    st.eppos2 = 1e-6;   st.epneg = 0;
end
dd1 = st.dd1;   dd12 = st.dd12;     options1 = st.options1;     options12 = st.options12;
override1 = st.override1;   epneg = st.epneg;   eppos2 = st.eppos2;
if st.sosineq_on                                    % stock lines 153-155
    error('cx_stability:gap',['sosineq_on=1 needs lpi_ineq, which has no container '...
          'counterpart; only the equality path (the default) is transcribed.'])
end
override2 = st.override2;   options2 = st.options2;     options3 = st.options3;
dd2 = st.dd2;   dd3 = st.dd3;

Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);        % line 108

% P = P1 (+P2) + eppos2*T*T', lines 115-125.
[prog,Pm] = cx_stability_pos(prog,Tm,'out',PIE.dom,dd1,options1);
if override1~=1
    [prog,P2m] = cx_stability_pos(prog,Tm,'out',PIE.dom,dd12,options12);
    Pm = Pm + P2m;
end
Pm = Pm + eppos2*Tm*Tm';

% Q: T.out -> T.in (the shape of Top2 = Top'), T*Q = P, lines 128-130.
Qdeg = cx_stability_lpivar_degs(Pm);
[spO,dmO] = cx_space_list(Tm,'in');     [spI,dmI] = cx_space_list(Tm,'out');
[prog,Qm] = lpivar_cdopvar(prog,struct('out',dmO,'in',dmI), ...
                           struct('out',{spO},'in',{spI}),PIE.dom,Qdeg);
prog = lpi_eq_cdopvar(prog,Tm*Qm - Pm);             % NOT symmetric (line 130)

% D, line 142.
Dm = Am*Qm + Qm'*Am';
if epneg~=0
    Dm = Dm + epneg*Pm;
end

% D = -N, lines 159-167.
[prog,Nm] = cx_stability_pos(prog,Dm,'out',PIE.dom,dd2,options2);
if override2~=1
    [prog,N2m] = cx_stability_pos(prog,Dm,'out',PIE.dom,dd3,options3);
    Nm = Nm + N2m;
end
prog = lpi_eq_cdopvar(prog,Dm+Nm,'symmetric');
end
