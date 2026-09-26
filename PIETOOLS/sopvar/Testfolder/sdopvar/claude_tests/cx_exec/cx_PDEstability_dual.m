function prog = cx_PDEstability_dual(PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_PDESTABILITY_DUAL(PIE,ST) builds, with copvar/cdopvar only, the
% LPI of executives/PIETOOLS_PDEstability_dual.m (direct form, dual PIE):
%
%   P = P1 (+P2) + diag(eppos I, eppos2 I) >= 0,
%   D = T*P*A' + A*P*T' + epneg*T*P*T' = -N,   N = N1 (+N2) >= 0,
%
% returned UNSOLVED. Mirrors lines 62-70 (2-D hand-off to the
% stability_dual_2D transcription), 84-102 (settings), 109 (program),
% 116-127 (P), 138 (D), 150-165 (equality). Not mirrored: line 152
% (sosineq_on=1: no container lpi_ineq, errors) and line 171 (no getsol).
%
% P lives on T's INPUT spaces (T*P*A'); for a PIE these equal the output
% spaces. Everything else as in cx_PDEstability.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(PIE,'pie_struct')
    error('The PIE should be a ''pie_struct''.')
end
PIE = initialize(PIE);
if PIE.dim==2                                       % stock lines 62-70
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
override1 = st.override1;   eppos = st.eppos;   epneg = st.epneg;   eppos2 = st.eppos2;
if st.sosineq_on                                    % stock lines 150-152
    error('cx_stability:gap',['sosineq_on=1 needs lpi_ineq, which has no container '...
          'counterpart; only the equality path (the default) is transcribed.'])
end
override2 = st.override2;   options2 = st.options2;     options3 = st.options3;
dd2 = st.dd2;   dd3 = st.dd3;

Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);        % line 109

% P = P1 (+P2), lines 116-123.
[prog,Pm] = cx_stability_pos(prog,Tm,'in',PIE.dom,dd1,options1);
if override1~=1
    [prog,P2m] = cx_stability_pos(prog,Tm,'in',PIE.dom,dd12,options12);
    Pm = Pm + P2m;
end
% Strict positivity, lines 126-127.
n = PIE.T.dim(:,2);
Imat = blkdiag(eppos*eye(n(1)),eppos2*eye(n(2)));
Pm = Pm + opvar2copvar(mat2opvar(Imat,PIE.T.dim(:,2),PIE.vars,PIE.dom));

% D, line 138, stock association ((T*P)*A').
Dm = Tm*Pm*Am' + Am*Pm*Tm';
if epneg~=0
    Dm = Dm + epneg*(Tm*Pm*Tm');
end

% D = -N, lines 156-164.
[prog,Nm] = cx_stability_pos(prog,Dm,'out',PIE.dom,dd2,options2);
if override2~=1
    [prog,N2m] = cx_stability_pos(prog,Dm,'out',PIE.dom,dd3,options3);
    Nm = Nm + N2m;
end
prog = lpi_eq_cdopvar(prog,Nm+Dm,'symmetric');                 % line 164 order
end
