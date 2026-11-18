function [PDE_out] = terms2pde_struct(terms_in)
% PDE_OUT = TERMS2PDE_STRUCT(TERMS_IN) takes an object of type 'terms'
% and converts it to an object of type 'pde_struct', representing the same 
% terms.
%
% INPUT
% - terms_in:   'terms' class object, representing a set of terms to define
%               a PDE.
%
% OUTPUT
% - PDE_out:    'pde_struct' object representing the same terms as the
%               input "terms_in".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 06/30/2024

% % % Process the inputs
% Check that the input terms are indeed of type 'terms'
if ~isa(terms_in,'terms')
    error("Input must be an object of type 'terms'.")
end

% % 
% Convert the state vector in the terms to a PDE structure.
PDE_states = state2pde_struct(terms_in.statevec);
% Extract the operator defining the scaling/integrals of the terms.
Cop = terms_in.operator;
if ~isa(Cop,'opvar')
    error("Conversion of 'terms' objects to 'pde_struct' object is currently only supported for 1D systems.")
end

% % Determine which components in the vector belong on the left-hand side
% % of the PDE.
vec_size = zeros(numel(PDE_states.free),1);
is_LHS_obj = false(numel(PDE_states.free),1);
for ii=1:numel(PDE_states.free)
    vec_size(ii) = PDE_states.free{ii}.size;
    is_LHS_obj(ii) = is_LHS_term(PDE_states.free{ii}.term{1});
end

% % % If no LHS objects appear, we can just multiply.
if ~any(is_LHS_obj)
    PDE_out = Cop*PDE_states;
    return
end

% % % Otherwise, extract columns of Cop that correspond to LHS objects.
% % Extract elements of "PDE_states" corresponding to LHS and RHS objects.
PDE_states_LHS = PDE_states;    PDE_states_LHS.free = PDE_states_LHS.free(is_LHS_obj);
PDE_states_RHS = PDE_states;    PDE_states_RHS.free = PDE_states_RHS.free(~is_LHS_obj);

% % Create a list of indices of columns in Cop corresponding to RHS and LHS
% % objects.
% Create a matrix with starting and ending columns in Cop corresponding to
% each state component.
vec_size_cum = [0;cumsum(vec_size)];
full_idcs = [vec_size_cum(1:end-1)+1,vec_size_cum(2:end)];
% Extract start and end indices corresponding to LHS objects
LHS_idcs = full_idcs(is_LHS_obj,:)';
% Convert to a list of columns in Cop corresponding to LHS objects.
LHS_idcs = mat2cell(LHS_idcs,2,ones(1,size(LHS_idcs,2)));
LHS_idcs = cellfun(@(A) (A(1,1):A(2,1)),LHS_idcs,'UniformOutput',false);
LHS_idcs = cell2mat(LHS_idcs);
% % Set a list of columns in Cop corresponding to RHS objects.
RHS_idcs = 1:vec_size_cum(end);
RHS_idcs(LHS_idcs) = [];

% % Check that no integral is taken of LHS objects.
Cop_LHS = Cop(1:size(Cop,1),LHS_idcs);
if ~isempty(RHS_idcs)
    Cop_RHS = Cop(1:size(Cop,1),RHS_idcs);
else
    Cop_RHS = 0;
end
if any(Cop.dim(:,2))
    Cint_LHS = polynomial([Cop_LHS.Q1; Cop_LHS.R.R1 +Cop_LHS.R.R2]);
    if any(any(abs(Cint_LHS.C)>1e-10))
        error("Integration of outputs or temporal derivatives of state components is not supported.")
    end
end
% Conver Cop_LHS to just a matrix;
Cop_LHS = cleanpoly(polynomial([Cop_LHS.P, Cop_LHS.Q1; Cop_LHS.Q2, Cop_LHS.R.R0]),1e-10);
if ~isdouble(Cop_LHS)
    error("Multiplication of outputs or temporal derivatives of state components with polynomials is not supported.")
end
Cop_LHS = double(Cop_LHS);

% % % Set a separate equation for each LHS object.
% % Determine which rows in Cop_LHS correspond to which object.
vec_size_LHS = vec_size(is_LHS_obj);
vec_size_LHS_cum = [0;cumsum(vec_size_LHS)];
LHS_comp_num = zeros(size(Cop_LHS,1),1);
for ii=1:size(Cop_LHS,1)
    % Check which columns in Cop have nonzero elements.
    colnums_ii = find(Cop_LHS(ii,:)~=0);
    % Check which LHS component these columns correspond to.
    is_LHS_comp_ii = any((colnums_ii>=vec_size_LHS_cum(1:end-1)+1) & (colnums_ii<=vec_size_LHS_cum(2:end)),2);
    if sum(is_LHS_comp_ii,1)>1
        error("Some equations appear to involve more than one output or or temporal derivative of a state component; this is not supported.")
    elseif any(is_LHS_comp_ii)
        LHS_comp_num(ii) = find(is_LHS_comp_ii,1,'first');
    end
end

% % Finally, for each LHS object, declare an equation in 'pde_struct'
% % format.
PDE_out = [];
for kk=0:sum(is_LHS_obj)
    % Determine which rows in Cop_LHS correspond to LHS object kk.
    ridcs_kk = find(LHS_comp_num==kk);
    if isempty(ridcs_kk)
        continue
    end
    Cop_LHS_kk = Cop_LHS(ridcs_kk',1:end);
    if isa(Cop_RHS,'opvar')
        Cop_RHS_kk = Cop_RHS(ridcs_kk',1:size(Cop_RHS,2));
        PDE_out = [PDE_out; Cop_LHS_kk*PDE_states_LHS+Cop_RHS_kk*PDE_states_RHS];
    else
        PDE_out = [PDE_out; Cop_LHS_kk*PDE_states_LHS];
    end
end

end



function [is_LHS] = is_LHS_term(PDE_term)
% IS_LHS_TERM Checks if a 'struct' object "PDE_term" corresponds to a term
% in a 'pde_struct' object, belonging on the left-hand side of a PDE (i.e.
% involving an output signal or a temporal derivative of a state.
%
% INPUT
% - PDE_term:   'struct' object, representing a term in a 'pde_struct'
%               object.
% 
% OUTPUT
% - is_LHS:     Boolean variable specifying whether the term corresponds to
%               the left-hand side of the PDE, i.e. involves an output
%               signal or a temporal derivative of a state.

if ~isa(PDE_term,'struct')
    error("Input term should be of type 'struct'")
end

% Assume the input corresponds to just a standard term.
is_LHS = false;

% Check if the term corresponds to a temporal derivative of a state, or to
% an output.
if isfield(PDE_term,'x') && isfield(PDE_term,'tdiff') && PDE_term.tdiff>0
    is_LHS = true;
elseif isfield(PDE_term,'y')
    is_LHS = true;
elseif isfield(PDE_term,'z')
    is_LHS = true;
end

end