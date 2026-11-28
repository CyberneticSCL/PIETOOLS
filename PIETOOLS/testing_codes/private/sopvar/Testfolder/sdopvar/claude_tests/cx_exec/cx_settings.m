function st = cx_settings(name,solver)
% ST = CX_SETTINGS(NAME,SOLVER) is lpisettings(NAME) with the solver set
% and simplification off, as the baseline suite's bl_settings does (which is
% in a private/ folder and so not callable from here). SOLVER defaults to
% 'mosek', the baseline's; use 'sedumi' in 2-D, where Mosek reports
% false-feasible stability verdicts at light settings.
%
% Initial coding MMP, 09/25/2026
% MMP, 11/08/2026: lpisettings_sop(NAME) (was the stock lpisettings(NAME),
%                whose folder leaves the path); help text updated.
% MMP, 11/12/2026: ST carries sop.stock, the stock slack of the tier
%                (LPISETTINGS_STOCK_SOP), which left the tiers on 11/12/2026
%                and which the cx transcriptions of the stock executives
%                read (POSLPIVAR_SETTINGS_SOP, role 'slack'); sop.slack stays
%                'like', so the container executives pose what they did.
% MMP, 11/23/2026: ST carries lpi.pos in the stock-field form (mult, int,
%                joint, sep) that the transcriptions read (role 'lf'):
%                entry 1 of each tier on 11/12/2026 (019b383e). The tiers
%                hold the storage row since 11/18/2026 (bb367fe1), on which
%                'lf' errors.
if nargin<2 || isempty(solver),     solver = 'mosek';   end
% st = lpisettings(name);                                                   % MMP, 11/08/2026 (was)
st = lpisettings(name);                                                     % MMP, 11/08/2026
[~,st.sop.stock] = lpisettings_stock(name);         % the stock slack       % MMP, 11/12/2026
% lpi.pos of role 'lf': mult, int, joint of each tier on 11/12/2026         % MMP, 11/23/2026
f = struct('stripped',[1 1 2],'light',[1 1 2],'custom',[1 1 2], ...
           'heavy',[2 1 2],'veryheavy',[3 3 6],'extreme',[1 -1 0]);         % MMP, 11/23/2026
f = f.(lower(name));                                                        % MMP, 11/23/2026
st.lpi.pos = struct('mult',f(1),'int',f(2),'joint',f(3),'sep',0);           % MMP, 11/23/2026
st.sos_opts.solver   = solver;
st.sos_opts.simplify = false;
end
