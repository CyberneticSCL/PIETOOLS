function solcoeff = PIESIM_time_integrate(psize, opts, uinput, coeff, gridall, Dop)
% Phase-0 capture shadow of PIESIM_time_integrate. P0_F1 puts this folder
% first on the path, after storing a handle to the original in P0_CAPTURE.h;
% this function records the discretized operators Dop, the initial
% coefficients, the grids and the sizes that PIESIM passes to its time
% integrator, then calls the original. It records data only; the solution
% PIESIM returns is unchanged.
global P0_CAPTURE
P0_CAPTURE.n = P0_CAPTURE.n + 1;
P0_CAPTURE.psize = psize;
P0_CAPTURE.coeff = coeff;
P0_CAPTURE.gridall = gridall;
P0_CAPTURE.Dop = Dop;
solcoeff = P0_CAPTURE.h(psize, opts, uinput, coeff, gridall, Dop);
end
