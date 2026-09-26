function prog = cx_PDEstability(PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_PDESTABILITY(PIE,ST) builds, with copvar/cdopvar only, the LPI
% that executives/PIETOOLS_PDEstability.m builds (direct form, Cor. 7.3):
%
%   P = P1 (+P2) + diag(eppos I, eppos2 I) >= 0,
%   D = T'*P*A + A'*P*T + epneg*T'*P*T = -N,   N = N1 (+N2) >= 0,
%
% and returns it UNSOLVED, for cx_compare. Mirrors lines 65-73 (2-D hand-off
% to the stability_2D transcription), 86-104 (settings), 110 (program),
% 117-128 (P), 139 (D), 150-165 (negativity by equality). Not mirrored:
% line 152 (sosineq_on=1, lpi_ineq: no container counterpart, errors here)
% and line 172 (getsol_lpivar: no getsol for cdopvar; nothing is solved).
%
% Container specifics:
% - P lives on T's OUTPUT spaces (T'*P*A), read from the converted Tm, so a
%   pure-L2 PIE gives a 1x1 container (opvar2copvar.m:129-133).
% - poslpivar -> poscopvar through cx_stability_pos (degrees, psatz, sep,
%   exclude); the slack N is declared on D's own output spaces.
% - The strictness identity is opvar2copvar(mat2opvar(...)): no container
%   identity constructor exists.
% - epneg*T'*P*T is added only when epneg~=0: the stock adds 0*(...) at
%   epneg=0, the zero operator. Measured: identical SDP shape (ndv, K.s, m)
%   with the term skipped (cx_stability_check), and at epneg=0.1 with it.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(PIE,'pie_struct')
    error('The PIE should be a ''pie_struct''.')
end
PIE = initialize(PIE);
if PIE.dim==2                                       % stock lines 65-73
    if nargin<2,    prog = cx_stability_2D(PIE);
    else,           prog = cx_stability_2D(PIE,st);
    end
    return
end
if nargin<2                                         % stock lines 78-84
    st = settings_PIETOOLS_heavy;
    st.sos_opts.simplify = 1;
    st.eppos = 1e-4;    st.eppos2 = 1e-6;   st.epneg = 0;
end
dd1 = st.dd1;   dd12 = st.dd12;     options1 = st.options1;     options12 = st.options12;
override1 = st.override1;   eppos = st.eppos;   epneg = st.epneg;   eppos2 = st.eppos2;
if st.sosineq_on                                    % stock line 150-152
    error('cx_stability:gap',['sosineq_on=1 needs lpi_ineq, which has no container '...
          'counterpart; only the equality path (the default) is transcribed.'])
end
override2 = st.override2;   options2 = st.options2;     options3 = st.options3;
dd2 = st.dd2;   dd3 = st.dd3;

Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);        % line 110

% P = P1 (+P2), lines 117-124.
[prog,Pm] = cx_stability_pos(prog,Tm,'out',PIE.dom,dd1,options1);
if override1~=1
    [prog,P2m] = cx_stability_pos(prog,Tm,'out',PIE.dom,dd12,options12);
    Pm = Pm + P2m;
end
% Strict positivity, lines 127-128: n0 = R^n count, n1 = L2 count of T's output.
n = PIE.T.dim(:,1);
Imat = blkdiag(eppos*eye(n(1)),eppos2*eye(n(2)));
Pm = Pm + opvar2copvar(mat2opvar(Imat,PIE.T.dim(:,1),PIE.vars,PIE.dom));

% D, line 139, same association as the stock ((T'*P)*A).
Dm = Tm'*Pm*Am + Am'*Pm*Tm;
if epneg~=0
    Dm = Dm + epneg*(Tm'*Pm*Tm);
end

% D = -N with N = N1 (+N2) >= 0, lines 156-164.
[prog,Nm] = cx_stability_pos(prog,Dm,'out',PIE.dom,dd2,options2);
if override2~=1
    [prog,N2m] = cx_stability_pos(prog,Dm,'out',PIE.dom,dd3,options3);
    Nm = Nm + N2m;
end
prog = lpi_eq_cdopvar(prog,Dm+Nm,'symmetric');
end
