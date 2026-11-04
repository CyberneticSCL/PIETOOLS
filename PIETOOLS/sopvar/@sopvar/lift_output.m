function Pop_new = lift_output(Pop,vars,dom,dims)
% POP_NEW = LIFT_OUTPUT(POP,VARS,DOM,DIMS) returns the 'sopvar'
% object POP_NEW representing the PI operator POP with its output lifted
% to the function space in the variables VARS.out: the output is extended
% with the variables in VARS.out that it is not yet defined in, and is
% constant along these new variables,
%   (Pop_new*x)(s,t) = (Pop*x)(s),
% for s the output variables of Pop, and t the new output variables.
%
% INPUTS
% - Pop:    m x n 'sopvar' object representing an operator
%               Pop: L2^n[dom_in] --> L2^m[dom_out];
% - vars:   1 x M cellstr specifying the names of the output variables of
%           the lifted operator. Must contain all output variables of Pop.
%           Can also be specified as struct with field 'out' containing
%           this cellstr (other fields, e.g. 'in', are ignored);
% - dom:    M x 2 'double' array specifying the intervals on which these
%           variables are defined. Must match the domains of the output
%           variables of Pop. Can also be specified as struct with field
%           'out' (other fields are ignored);
% - dims:   (optional) number of rows of the lifted operator, specified as
%           an integer, as 1 x 2 array [rows,columns], or as struct with
%           field 'out'. Defaults to the number of rows of Pop. If Pop is
%           scalar-valued (m=1), it can be lifted to a vector-valued
%           operator, repeating the scalar output in each row;
%
% OUTPUTS
% - Pop_new: m x n 'sopvar' object (or dims(1) x n, if specified)
%           representing the lifted operator
%               Pop_new: L2^n[dom_in] --> L2^m[dom_out_new];
%
% NOTES
% - A new output variable that is not an input variable of Pop is
%   introduced as a multiplier variable, for which the output is the
%   constant monomial.
% - A new output variable t that is also an input variable of Pop, so
%   that Pop integrates x over t, becomes a variable along which Pop_new
%   acts as a PI operator. The (full) integral along t is then expressed
%   as the sum of the integral from a to t and from t to b, with the same
%   parameters, so that the output does not depend on t. The domain of t
%   must match that of the input variable.
% - Lifting changes only the output space. The input space is unchanged.
% - Since the parameters of a 'sopvar' object are constant coefficients,
%   the cost is independent of any decision variables.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lift_output
%
% Copyright (C) 2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 11/03/2026: Initial coding;

% Check that the inputs are properly specified
if ~isa(Pop,'sopvar')
    error("Operator must be specified as object of type 'sopvar'.")
end
if nargin<=3
    dims = [];
end
[m,n] = size(Pop);
[vars_out,dom_out,m_new] = parse_inputs(vars,dom,dims,Pop.vars,Pop.dom,[m,n]);

% Determine the new output variables
vars_old = reshape(Pop.vars.out,1,[]);
vars_in = reshape(Pop.vars.in,1,[]);
is_new = ~ismember(vars_out,vars_old);
vars_add = vars_out(is_new);
dom_add = dom_out(is_new,:);
[is_in,idcs_in] = ismember(vars_add,vars_in);
if any(is_in) && any(any(dom_add(is_in,:)~=Pop.dom.in(idcs_in(is_in),:)))
    error("A new output variable that is also an input variable must have the same domain.")
end

Pop_new = Pop;
if ~isempty(vars_add)
    % Declare the output as a constant function of the new variables: the
    % monomial in each new variable is s^0, so that the monomial vector
    % ZL, and with it the rows of the parameters, is unchanged
    vars_new = struct('in',{Pop.vars.in},'out',{[vars_old,vars_add]});
    dom_new = struct('in',Pop.dom.in,'out',[Pop.dom.out;dom_add]);
    ZL_new = [reshape(Pop.ZL,1,[]),num2cell(zeros(1,numel(vars_add)))];
    params_new = lift_params(Pop,vars_in,vars_old,vars_add(is_in),[m,n]);
    % The constructor puts the variables in the canonical order, and
    % permutes the monomial indices of the parameters to match
    Pop_new = sopvar(params_new,vars_new,ZL_new,Pop.ZR,dom_new,[m,n]);
end

% Lift a scalar output to a vector: since the parameters are stacked with
% the matrix row as outer index (see sopvar), repeating the scalar
% operator along the rows repeats the parameters
if m_new~=m
    Pop_new = repmat(Pop_new,m_new,1);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function params_new = lift_params(Pop,vars_in,vars_old,vars_S1to3,matdim)
% PARAMS_NEW = LIFT_PARAMS(POP,VARS_IN,VARS_OLD,VARS_S1TO3,MATDIM) returns
% the cell of parameters of the lifted operator, in which the input
% variables VARS_S1TO3, so far integrated out, become pass-through
% variables.
%
% INPUTS
% - Pop:    'sopvar' object that is lifted;
% - vars_in, vars_old: input and output variables of Pop;
% - vars_S1to3: 1 x K cellstr of input variables of Pop that are new output
%           variables;
% - matdim: 1 x 2 array of the matrix dimensions of Pop;
%
% OUTPUTS
% - params_new: 3 x ... x 3 cell of sparse matrices, for the pass-through
%           variables S3 = intersect(vars_in,[vars_old,vars_S1to3]), in the
%           order intersect returns (see sopvar);
%
params_old = Pop.params;
if isempty(vars_S1to3)
    % No new pass-through variables, so the parameters are unchanged
    params_new = params_old;
    return
end

% Set the parameters along the pass-through directions, in the order of
% intersect, as for sopvar
S3_old = intersect(vars_in,vars_old);               N3_old = numel(S3_old);
S3_new = intersect(vars_in,[vars_old,vars_S1to3]);  N3_new = numel(S3_new);
[~,pos_old] = ismember(S3_old,S3_new);
pos_new = ~ismember(S3_new,S3_old);
sz_old = [3*ones(1,N3_old),1];
sz_new = [3*ones(1,N3_new),1];

% Zero parameter, for the multiplier term along a new direction
nZL = prod(cellfun(@numel,Pop.ZL));
nZR = prod(cellfun(@numel,Pop.ZR));
Z0 = sparse(matdim(1)*nZL,matdim(2)*nZR);

% Along a new direction t, the full integral of the operator is the sum of
% the integrals over [a,t] (index 2) and [t,b] (index 3): the parameters of
% Pop are used for both, and the multiplier term (index 1) is zero. The
% result is independent of t.
params_new = cell(sz_new);
for lidx=1:numel(params_new)
    idcs = cell(1,N3_new+1);
    [idcs{:}] = ind2sub(sz_new,lidx);
    idcs = cell2mat(idcs);
    if any(idcs(pos_new)==1)
        params_new{lidx} = Z0;
    elseif N3_old==0
        params_new{lidx} = params_old{1};
    else
        idcs_old = num2cell([idcs(pos_old),1]);
        params_new{lidx} = params_old{sub2ind(sz_old,idcs_old{:})};
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [vars_out,dom_out,dim_out] = parse_inputs(vars,dom,dims,vars0,dom0,matdim0)
% [VARS_OUT,DOM_OUT,DIM_OUT] = PARSE_INPUTS(VARS,DOM,DIMS,VARS0,DOM0,
% MATDIM0) checks the inputs specifying the output space to which the
% operator is lifted, and returns them in standard form.
%
% INPUTS
% - vars, dom: the output variables (cellstr) and their domains (M x 2
%           array) of the lifted operator, or structs with field 'out'
%           containing these;
% - dims:   [], or the number of rows of the lifted operator, specified as
%           a scalar, as 1 x 2 array [rows,columns], or as struct with
%           field 'out';
% - vars0, dom0: structs with field 'out', the output variables and
%           domains of the operator that is lifted;
% - matdim0: 1 x 2 array of the matrix dimensions of that operator;
%
% OUTPUTS
% - vars_out: 1 x M cellstr of output variables of the lifted operator;
% - dom_out: M x 2 'double' array of the associated domains;
% - dim_out: number of rows of the lifted operator: matdim0(1), or dims if
%           specified, in which case matdim0(1) must be 1 or equal dims;
%
% Check that the output variables are properly specified, as cellstr or as
% struct with field 'out'
if isa(vars,'struct') && isfield(vars,'out')
    vars = vars.out;
end
if ~iscellstr(vars)
    error("Output variable names should be specified as 'cellstr' object, or struct with field 'out'.")
end
vars_out = reshape(vars,1,[]);          M = numel(vars_out);
if numel(unique(vars_out))~=M
    error("Output variable names must be unique.")
end

% Check that the domains are properly specified, as numeric array or as
% struct with field 'out'
if isa(dom,'struct') && isfield(dom,'out')
    dom = dom.out;
end
if ~isnumeric(dom) || ~all(size(dom)==[M,2])
    error("Number of output domains should match number of output variables, and be specified as M x 2 array or struct with field 'out'.")
end
dom_out = dom;

% The output variables of the operator must be kept, on the same domain,
% as lifting extends the output function space but does not restrict it
vars_old = reshape(vars0.out,1,[]);
[is_old,idcs_old] = ismember(vars_old,vars_out);
if ~all(is_old)
    error("The output variables of the operator must be among the output variables to lift to.")
end
if ~isempty(vars_old) && any(any(dom_out(idcs_old,:)~=dom0.out & ~(isnan(dom_out(idcs_old,:)) & isnan(dom0.out))))
    error("The domains of the output variables of the operator cannot be changed.")
end

% Check that the number of rows is properly specified
dim_out = matdim0(1);
if ~isempty(dims)
    if isa(dims,'struct') && isfield(dims,'out')
        dims = dims.out;
    end
    if ~isnumeric(dims) || ~(isscalar(dims) || numel(dims)==2) || any(dims~=round(dims)) || any(dims<0)
        error("Dimensions must be specified as nonnegative integer, 1x2 array of nonnegative integers, or struct with field 'out'.")
    elseif numel(dims)==2 && dims(2)~=matdim0(2)
        error("Number of columns cannot be changed when lifting the output.")
    end
    dim_out = dims(1);
    % Only a scalar output can be lifted to a vector: all components of
    % the lifted output are then equal to this scalar
    if dim_out~=matdim0(1) && matdim0(1)~=1
        error("Only a scalar-valued output can be lifted to a vector-valued output.")
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
