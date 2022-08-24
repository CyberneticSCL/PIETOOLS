%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_PDE_terms.m     PIETOOLS 2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIE_out = convert_PIETOOLS_PDE(PDE)
% convert_PIETOOLS_PDE is the new version of the PDE converter 
% file which performs the following two tasks.
% 1) It verifies the dimension compatibility of input parameters of ODE-PDE
% and sets any missing parameters to zero.
% 2) It converts the input ODE-PDE representation to a PIE
% representation.

% A Partial Integral Equation is defined by 12 PI operators as
%
% Tw \dot{w}(t) + Tu \dot{u}(t) + T \dot{x}(t) = A  x(t) + B1 w(t)  + B2 u(t);
%                                         z(t) = C1 x(t) + D11 w(t) + D12op u(t);
%                                         y(t) = C2 x(t) + D21 w(t) + D22op u(t);
%
% This script takes a user-defined PDE system in the format outlined in the
% header of "pde_struct" and converts it to a PIE by defining the 11
% PI operators {T,  Tw,  Tu;
%               A,  B1,  B2;
%               C1, D11, D12;
%               C2, D21, D22}; 
%
% Depending on the structure of the input "PDE", the "terms", "batch", or
% "terms_old" conversion function will be called. See the separate
% functions and the manual for more information on the input format.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, SS, DJ, 08/18/2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call the appropriate converter function.
if isa(PDE,'pde_struct') || isfield(PDE,'x')
    if PDE.dim<=1
        PIE_out = convert_PIETOOLS_PDE_terms(PDE);
    elseif PDE.dim==2
        PIE_out = convert_PIETOOLS_PDE_2D(PDE);
    else
        error('Conversion of PDEs in more than 2 spatial variables to PIEs is currently not supported.')
    end
elseif isfield(PDE,'n0') || isfield(PDE,'n1') || isfield(PDE,'n2')
    PIE_out = convert_PIETOOLS_PDE_batch(PDE);
elseif isfield(PDE,'n') || isfield(PDE,'ODE')  || isfield(PDE,'PDE')
    PIE_out = convert_PIETOOLS_PDE_terms_legacy(PDE);
else
    error('The input PDE is not appropriately specified. Please define your PDE as a "pde_struct" class object, and consult the manual and examples for illustration of the structure.')
end

end