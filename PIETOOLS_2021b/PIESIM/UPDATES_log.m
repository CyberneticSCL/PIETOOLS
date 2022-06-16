% PIESIM update of 06.16.2022

% examples_pde_library_PIESIM.m - updated to include terms format, 
% added more examples in terms format 

% solver_PIESIM.m - updated example count 

% PIESIM.m - changed arguments in transform to solution to allow for different 
% LHS of PDE (T) and map (T0) operators (between fundamental and primary states) -
% needed, for example, for Orr-Sommerfeld equation

% In folder "Discretization_of_Operators"

% PIESIM_discretize_icf.m -Renamed uinput.B21_nonpol to uinput.Bpw_nonpol, updated description of
% outputs

% PIESIM_4PI2Mat_cheb.m - added flag = 2 and corresponding discretizations - now supports
% computing of observed and regulated outputs

% PIESIM_discretize_all.m - updated comments to include C1, C2 operators

% PIESIM_discretize_ops.m - added discretization of C1, C2 operators to allow for
% computation of observed and regulated outputs. Separated T from T0 (LHS
% PDE from map operator).

% In folder PIESIM_setup

% PIESIM_input_check - Separated inputs check for batch and terms format
% of PDE. Enhanced the input check for terms format. Terms format has now
% its own structure ('PDT'), while batch format has 'PDE' structure. 
% Renamed uinput.B21_nonpol to uinput.Bpw_nonpol

% rescale_PIE.m - Separated PDE LHS operator (T) from the state map operator
% (T0)

% PIESIM_options_check.m - Added a clause -  "If piesize is not found, default to last argument"

% In folder "Postprocessing"

% PIESIM_plot_solution.m - added outputs for observed and regulated outputs 

% PIESIM_transform_to_solution.m - Added reconstruction of observed and regulated outputs.
% Added a support of case when LHS PDE operator is not the same as the
% state map operator. Note: input arguments changed.
