function PIE = initialize_PIETOOLS_PIE(PIE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE = initialize_PIETOOLS_PIE(PIE) takes a PIE data structure and
% checks that all the necessary fields are appropriately specified, and
% assigns a default value to all the optional fields that have not been
% specified.
% 
% INPUT
%   PIE: "struct" or "pie_struct" class object.
%
% OUTPUT
%   PIE: "pie_struct" class object describing the same system as the input.
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding SS - 08/01/2022

if nargin==0
    PIE = pie_struct();
    return
end

if ~isa(PIE,'struct') && ~isa(PIE,'pie_struct')
    error('Input must be a ''struct'' or ''pie_struct'' class object.')
end

if ~isfield(PIE,'dom')
    PIE.dom = [0,1];
end
if ~isfield(PIE,'dim')
    PIE.dim = 1;
end
if ~isfield(PIE,'vars')
    PIE.vars = [pvar('s'),pvar('theta')];
end

% store dimension information in the structure
dimstruct = struct('nx',0,'nw',0,'nu',0,'nz',0,'ny',0);
% first figure out the dimensions from existing fields
if ~isfield(PIE,'T')
    error("PIE.T is a mandatory element of PIE structure and cannot be defaulted to zero");
else
    if any(PIE.T.dim(:,1)~=PIE.T.dim(:,2))
        error("PIE.T operator must be a map between symmetric spaces");
    end
    dimstruct.nx = PIE.T.dim(:,1);
end
if ~isfield(PIE,'Tw')
    dimstruct.nw = 0;
else
    if any(PIE.T.dim(:,1)~=PIE.Tw.dim(:,1))
        warning('PIE.Tw and PIE.T map to different spaces. Defaulting Tw to zero operator');
    else
    end
end
if ~isfield(PIE,'Tu')
    
end
if ~isfield(PIE,'A')
    
end
if ~isfield(PIE,'B1')
    
end
if ~isfield(PIE,'B2')
    
end
if ~isfield(PIE,'C1')
    
end
if ~isfield(PIE,'C2')
    
end
if ~isfield(PIE,'D11')
    
end
if ~isfield(PIE,'D12')
    
end
if ~isfield(PIE,'D21')
    
end
if ~isfield(PIE,'D22')
    
end
end