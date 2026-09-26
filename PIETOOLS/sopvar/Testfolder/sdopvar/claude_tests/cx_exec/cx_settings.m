function st = cx_settings(name,solver)
% ST = CX_SETTINGS(NAME,SOLVER) is lpisettings(NAME) with the solver set and
% simplification off, as the baseline suite's bl_settings does (which is in
% a private/ folder and so not callable from here). SOLVER defaults to
% 'mosek', the baseline's; use 'sedumi' in 2-D, where Mosek reports
% false-feasible stability verdicts at light settings.
%
% Initial coding MMP, 09/25/2026
if nargin<2 || isempty(solver),     solver = 'mosek';   end
st = lpisettings(name);
st.sos_opts.solver   = solver;
st.sos_opts.simplify = false;
end
