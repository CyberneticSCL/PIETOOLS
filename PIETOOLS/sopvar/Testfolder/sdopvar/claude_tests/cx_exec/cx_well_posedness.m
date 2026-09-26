function prog = cx_well_posedness(PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_WELL_POSEDNESS(PIE,ST) builds, with copvar/cdopvar only, the LPI
% of executives/PIETOOLS_well_posedness.m (dissipativity + surjectivity,
% omega = -epneg):
%
%   P = P1 (+P2) + I1 >= 0,  I1 = diag(eppos I, eppos2 I),
%   R = R1 (+R2) + I  >= 0,
%   D1 = T'*P*A + A'*P*T - 2*omega*T'*P*T = -N1,      N1 >= 0,
%   D2 = (T-A)*R*(T-A)' - I1             = +N2,      N2 >= 0,
%
% returned UNSOLVED. Mirrors lines 60-62 (2-D refused, as the stock does),
% 82-100 (settings), 106 (program), 113-131 (P, R), 142-143 (D1, D2),
% 158-175 (equality). Not mirrored: lines 154-157 (sosineq_on=1: no
% container lpi_ineq, errors) and 182-183 (no getsol).
%
% Container specifics:
% - eye(size(Rop)) (line 130) is the TOTAL dimension; @cdopvar/size returns
%   the block grid, so the count is taken from PIE.T.dim instead.
% - P on T's output spaces (T'*P*A), R on T's input spaces ((T-A)*R*(T-A)').
% - The omega term is added only when omega~=0 (the stock adds 0*(...), the
%   zero operator; measured identical SDP shape either way).
% - Declaration order P1, R1, (P2, R2), N1a, N2a, (N1b, N2b) as in the stock.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(PIE,'pie_struct')
    error('The PIE should be a ''pie_struct''.')
end
PIE = initialize(PIE);
if PIE.dim==2                                       % stock lines 60-62
    error("Well-posedness analysis of 2D PDEs is currently not supported.")
end
if nargin<2                                         % stock lines 74-80
    st = settings_PIETOOLS_heavy;
    st.sos_opts.simplify = 1;
    st.eppos = 1e-2;    st.eppos2 = 1e-2;   st.epneg = 0;
end
dd1 = st.dd1;   dd12 = st.dd12;     options1 = st.options1;     options12 = st.options12;
override1 = st.override1;   eppos = st.eppos;   omega = -st.epneg;  eppos2 = st.eppos2;
if st.sosineq_on                                    % stock lines 154-157
    error('cx_stability:gap',['sosineq_on=1 needs lpi_ineq, which has no container '...
          'counterpart; only the equality path (the default) is transcribed.'])
end
override2 = st.override2;   options2 = st.options2;     options3 = st.options3;
dd2 = st.dd2;   dd3 = st.dd3;

Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);        % line 106

% P, R, lines 113-124.
[prog,Pm] = cx_stability_pos(prog,Tm,'out',PIE.dom,dd1,options1);
[prog,Rm] = cx_stability_pos(prog,Tm,'in',PIE.dom,dd1,options1);
if override1~=1
    [prog,P2m] = cx_stability_pos(prog,Tm,'out',PIE.dom,dd12,options12);
    Pm = Pm + P2m;
    [prog,R2m] = cx_stability_pos(prog,Tm,'in',PIE.dom,dd12,options12);
    Rm = Rm + R2m;
end
% Strictness, lines 127-131: P >= I1, R >= I (identity on the TOTAL state).
n = PIE.T.dim(:,1);
Imat = blkdiag(eppos*eye(n(1)),eppos2*eye(n(2)));
I1 = opvar2copvar(mat2opvar(Imat,PIE.T.dim(:,1),PIE.vars,PIE.dom));
Pm = Pm + I1;
nR = sum(PIE.T.dim(:,2));
Rm = Rm + opvar2copvar(mat2opvar(eye(nR),PIE.T.dim(:,2),PIE.vars,PIE.dom));

% D1, D2, lines 142-143 (stock association).
Dm1 = Tm'*Pm*Am + Am'*Pm*Tm;
if omega~=0
    Dm1 = Dm1 - 2*omega*(Tm'*Pm*Tm);
end
TA = Tm - Am;
Dm2 = TA*Rm*TA' - I1;

% D1 = -N1, D2 = N2, lines 161-174.
[prog,N1m] = cx_stability_pos(prog,Dm1,'out',PIE.dom,dd2,options2);
[prog,N2m] = cx_stability_pos(prog,Dm2,'out',PIE.dom,dd2,options2);
if override2~=1
    [prog,N1b] = cx_stability_pos(prog,Dm1,'out',PIE.dom,dd3,options3);
    N1m = N1m + N1b;
    [prog,N2b] = cx_stability_pos(prog,Dm2,'out',PIE.dom,dd3,options3);
    N2m = N2m + N2b;
end
prog = lpi_eq_cdopvar(prog,Dm1+N1m,'symmetric');
prog = lpi_eq_cdopvar(prog,Dm2-N2m,'symmetric');
end
