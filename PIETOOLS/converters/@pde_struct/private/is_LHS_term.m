function [is_LHS,obj,idx,tdiff] = is_LHS_term(PDE_term)
%IS_LHS_TERM Summary of this function goes here
%   Detailed explanation goes here

if ~isa(PDE_term,'struct')
    error("Input term should be of type 'struct'")
end

% Assume the input corresponds to just a standard term.
is_LHS = false;
obj = '';
idx = 0;
tdiff = 0;

% Check if the term corresponds to a temporal derivative of a state, or to
% an output.
if isfield(PDE_term,'x') && isfield(PDE_term,'tdiff') && PDE_term.tdiff>0
    is_LHS = true;
    obj = 'x';
    idx = PDE_term.x;
    tdiff = PDE_term.tdiff;
elseif isfield(PDE_term,'y')
    is_LHS = true;
    obj = 'y';
    idx = PDE_term.y;
elseif isfield(PDE_term,'z')
    is_LHS = true;
    obj = 'z';
    idx = PDE_term.z;
end

end