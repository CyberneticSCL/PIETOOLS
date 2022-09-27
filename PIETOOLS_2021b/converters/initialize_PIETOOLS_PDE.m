%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize_PIETOOLS_PDEs.m     PIETOOLS 2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PDE_out = initialize_PIETOOLS_PDE(PDE,silent_initialize)
% Initializes and checks the PDE formulation.
% Depending on the structure of the input "PDE", the "terms", "batch", or
% "terms_old" initialization function will be called. See the separate
% functions and the manual for more information on the input format.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ, 08/18/2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call the appropriate initialization function.
if isa(PDE,'pde_struct') || isa(PDE,'struct') && isfield(PDE,'x')
    if isa(PDE,'struct')
        PDE = pde_struct(PDE);
    end
    if nargin==1
        PDE_out = initialize(PDE);
    else
        PDE_out = initialize(PDE,silent_initialize);
    end
elseif isfield(PDE,'n0') || isfield(PDE,'n1') || isfield(PDE,'n2')
    PDE_out = initialize_PIETOOLS_PDE_batch(PDE);
elseif isfield(PDE,'n') || isfield(PDE,'ODE')  || isfield(PDE,'PDE')
    if nargin==2 && silent_initialize
        evalin('base','silent_initialize_pde=true;');
    end
    PDE_out = initialize_PIETOOLS_PDE_terms_legacy(PDE);
else
    error('The input PDE is not appropriately specified. Please define your PDE as a "pde_struct" class object, and consult the manual and examples for illustration of the structure.')
end

end