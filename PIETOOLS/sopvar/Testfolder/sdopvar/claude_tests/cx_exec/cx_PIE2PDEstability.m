function prog = cx_PIE2PDEstability(PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_PIE2PDESTABILITY(PIE,ST) builds, with copvar/cdopvar only, the
% LPI of executives/PIETOOLS_PIE2PDEstability.m (Q form, Cor. 7.1, the
% default 'stability' executive):
%
%   P = P1 (+P2) + eppos2*T'*T,     T'*Q - P = 0   (Q free),
%   D = A'*Q + Q'*A + epneg*P = -N,  N = N1 (+N2) >= 0,
%
% returned UNSOLVED. Mirrors lines 60-68 (2-D hand-off: the stock dispatches
% to the DIRECT-form PIETOOLS_stability_2D, so this goes to cx_stability_2D),
% 81-99 (settings), 105 (program), 112-120 (P), 123-125 (Q and T'Q = P),
% 137 (D), 148-163 (equality). Not mirrored: line 150 (sosineq_on=1, no
% container lpi_ineq, errors) and line 170 (no getsol for cdopvar).
%
% Container specifics:
% - Qdeg = get_lpivar_degs(Pop,Top) (line 123) has no container overload;
%   cx_stability_lpivar_degs reads the same maxima off the cdopvar P.
% - Q from lpivar_cdopvar with lpivar's [d1 d2 d3] convention, which
%   reproduces lpivar exactly in 1-D (lpivar_cdopvar.m:47-50, 429-430).
% - P lives on T's INPUT spaces (T'*T, T'*Q: in -> in).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(PIE,'pie_struct')
    error('The PIE should be a ''pie_struct''.')
end
PIE = initialize(PIE);
if PIE.dim==2                                       % stock lines 60-68
    if nargin<2,    prog = cx_stability_2D(PIE);
    else,           prog = cx_stability_2D(PIE,st);
    end
    return
end
if nargin<2                                         % stock lines 73-79
    st = settings_PIETOOLS_heavy;
    st.sos_opts.simplify = 1;
    st.eppos = 1e-4;    st.eppos2 = 1e-6;   st.epneg = 0;
end
dd1 = st.dd1;   dd12 = st.dd12;     options1 = st.options1;     options12 = st.options12;
override1 = st.override1;   epneg = st.epneg;   eppos2 = st.eppos2;
if st.sosineq_on                                    % stock lines 148-150
    error('cx_stability:gap',['sosineq_on=1 needs lpi_ineq, which has no container '...
          'counterpart; only the equality path (the default) is transcribed.'])
end
override2 = st.override2;   options2 = st.options2;     options3 = st.options3;
dd2 = st.dd2;   dd3 = st.dd3;

Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);        % line 105

% P = P1 (+P2) + eppos2*T'*T, lines 112-120.
[prog,Pm] = cx_stability_pos(prog,Tm,'in',PIE.dom,dd1,options1);
if override1~=1
    [prog,P2m] = cx_stability_pos(prog,Tm,'in',PIE.dom,dd12,options12);
    Pm = Pm + P2m;
end
Pm = Pm + eppos2*Tm'*Tm;                            % stock association (eppos2*T')*T

% Q free with T'*Q = P, lines 123-125. Q has T's shape (Top.dim).
Qdeg = cx_stability_lpivar_degs(Pm);
[spO,dmO] = cx_space_list(Tm,'out');    [spI,dmI] = cx_space_list(Tm,'in');
[prog,Qm] = lpivar_cdopvar(prog,struct('out',dmO,'in',dmI), ...
                           struct('out',{spO},'in',{spI}),PIE.dom,Qdeg);
prog = lpi_eq_cdopvar(prog,Tm'*Qm - Pm);            % NOT symmetric (line 125)

% D, line 137.
Dm = Am'*Qm + Qm'*Am;
if epneg~=0
    Dm = Dm + epneg*Pm;
end

% D = -N, lines 154-162.
[prog,Nm] = cx_stability_pos(prog,Dm,'out',PIE.dom,dd2,options2);
if override2~=1
    [prog,N2m] = cx_stability_pos(prog,Dm,'out',PIE.dom,dd3,options3);
    Nm = Nm + N2m;
end
prog = lpi_eq_cdopvar(prog,Dm+Nm,'symmetric');
end
