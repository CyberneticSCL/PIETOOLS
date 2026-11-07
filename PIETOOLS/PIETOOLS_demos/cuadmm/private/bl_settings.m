function st = bl_settings(name)
% CC, 11/07/2026: max_raise 0, so a case is one SDP, as with the stock
%   executives. The container executives re-solve at a raised slack degree
%   after 'unknown' or 'infeasible' (lpi_solve_loop_sop), so the dumped SDP
%   would depend on the solver's first verdict, and t_mosek (cpusec of the
%   returned program) would cover the last pass only.
st = lpisettings(name);
st.sos_opts.solver   = 'mosek';
st.sos_opts.simplify = false;
st.sop.max_raise = 0;   st.sop.max_raise_feas = 0;                          % CC, 11/07/2026
end


