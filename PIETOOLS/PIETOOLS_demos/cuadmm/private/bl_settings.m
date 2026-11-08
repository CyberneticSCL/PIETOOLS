function st = bl_settings(name)
% CC, 11/07/2026: max_raise 0, so a case is one SDP, as with the stock
%   executives. The container executives re-solve at a raised slack degree
%   after 'unknown' or 'infeasible' (lpi_solve_loop_sop), so the dumped SDP
%   would depend on the solver's first verdict, and t_mosek (cpusec of the
%   returned program) would cover the last pass only.
% CC, 11/08/2026: lpisettings_sop(NAME), the container tier, in place of the
%   stock lpisettings, whose folder leaves the path. Its sop options change
%   only the feasibility LPIs: stability and well-posedness take dw_feas and
%   eppos_stab 1e-2 (were dw_Q or dw_P, and eppos 1e-4). Measured over the
%   24 registry cases (MOSEK): objective cases identical; heavy stability
%   identical up to b x 100; light stability, well-posedness and 2-D stability
%   smaller (stab_rd1 m 126 -> 116, stab2_rd 2968 -> 2044) and certified, but
%   stab2_rd_psz90 (0.9 lam*) infeasible at this first pass. The tier's raise,
%   off here (max_raise 0), certifies it at the old program.
% st = lpisettings(name);                                                   % CC, 11/08/2026 (was)
st = lpisettings_sop(name);                                                 % CC, 11/08/2026
st.sos_opts.solver   = 'mosek';
st.sos_opts.simplify = false;
st.sop.max_raise = 0;   st.sop.max_raise_feas = 0;                          % CC, 11/07/2026
end


