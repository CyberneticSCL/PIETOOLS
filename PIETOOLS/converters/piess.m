function PIE = piess(Top,Aop,Bop,Cop,Dop)
% PIE = PIESS(TOP,AOP,BOP,COP,DOP) takes a set of PI operators, and 
% constructs the associated PIE rerpresentation. Declaring Top through Dop
% as 'opvar' or 'opvar2d' objects, the output PIE is a 'pie_struct' object
% representing a "5-PIE"
%   d/dt Top*v(t) = Aop*v(t) + Bop*w(t)
%            z(t) = Cop*v(t) + Dop*w(t)
% where v(t) is the fundamental state, w(t) the exogenous input
% (disturbance), and z(t) the regulated output, all possibly vector-valued
% and involving coupled finite- and infinite-dimensional components,
% determined by the dimensions of the specified operators.
% Alternatively, if cells of operators are specified as
%       Top = {Tvop,Twop,Tuop},     Bop = {Bwop,Buop};
%       Cop = {Czop;Cuop},          Dop = {Dzwop,Dzuop; Dywop,Dyup};
% the output PIE is a 'pie_struct' object representing the "12-PIE"
%   d/dt (Top*v(t)+Twop*w(t)+Tuop*u(t)) = Aop*v(t) +Bwop*w(t) +Buop*u(t)
%                                  z(t) = Czop*v(t)+Dzwop*w(t)+Dzuop*u(t)
%                                  y(t) = Cyop*v(t)+Dywop*w(t)+Dyuop*u(t)
% where now u(t) is the controlled input, and y(t) the observed output,
% with dimensions and domain again determined by the input operators.
%
% INPUT
% - Top:    nv x nv 'opvar' or 'opvar2d' object representing the
%           PI operator defining the left-hand side of a 5-PIE.
%           Alternatively, a 1x3 cell {Tvop,Twop,Tuop} with all elements 
%           either 'opvar' or 'opvar2d' objects of dimensions 
%           {nv x nv, nv x nw, nv x nu}, representing the left-hand
%           side of a 12-PI. The operators Twop and Tuop can also be set
%           equal to [], in which case they will default to 0 operators
%           of suitable dimensions;
% - Aop:    nv x nv 'opvar' or 'opvar2d' object representing the state
%           operator in both 5-PIE and 12-PIE dynamics.
% - Bop:    nv x nw 'opvar' or 'opvar2d' object representing the
%           input PI operator defining the contribution of disturbances w 
%           to the PDE dynamics.
%           Alternatively, a 1x2 cell {Bwop,Buop} with both elements 
%           either 'opvar' or 'opvar2d' objects of dimensions 
%           {nv x nw, nv x nu}, representing the contribution of 
%           disturbances and controlled inputs to the PDE dynamics;
%           Defaults to [], returning 0 operator(s) of suitable dimensions;
% - Cop:    nz x nv 'opvar' or 'opvar2d' object representing the
%           output PI operator defining the contribution of the state v 
%           to the regulated output z.
%           Alternatively, a 2x1 cell {Czop;Cyop} with both elements 
%           either 'opvar' or 'opvar2d' objects of dimensions 
%           {nz x nv, ny x nv}, representing the contribution of 
%           the state to the regulated and controlled output.
%           Defaults to [], returning 0 operator(s) of suitable dimensions;
% - Dop:    nz x nw 'opvar' or 'opvar2d' object representing the
%           feedthrough PI operator defining the effect of the 
%           disturbance w to the regulated output z.
%           Alternatively, a 2x2 cell {Dzwop, Dzuop; Dywop, Dyuop} with 
%           all elements either 'opvar' or 'opvar2d' objects of dimensions 
%           {nz x nw, nz x nu; ny x nw, ny x nu}, representing the
%           feedthrough from each type of input to each type of output.
%           Defaults to [], returning 0 operator(s) of suitable dimensions;
%
% OUTPUT
% - PIE:    'pie_struct' object representing the PIE defined by the
%           specified operators. The specified operators will be stored in
%           the PIE structure under the appropriate fields, so that e.g.
%           PIE.T = Top, PIE.Tw = Twop, PIE.A = Aop, etc..
%
% NOTES
% - Passing a single argument "PIE = piess(Top)", the input Top will be
% interpreted as operator Dzwop from disturbance to regulated output, and a
% PIE will be constructed as just z(t) = Dzwop*w(t). However, PIETOOLS
% currently offers no functionality for actually working with such systems.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESS
%
% Copyright (C)2024  PIETOOLS Team
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
% DJ, 12/20/2024: Initial coding;


% % % Process the inputs
% % Check how may input arguments are passed.
if nargin==0
    % Return empty pie_var structure.
    PIE = pie_var();
elseif nargin==1
    % Just as Matlab's `ss`, interpret single argument as direct map from
    % input to output.
    Dop = Top;
    Top = [];   Aop = [];   Bop = [];   Cop = [];
elseif nargin==2
    % Interpret two arguments as autonomous system d/dt(Top*v) = Aop*v.
    Bop = [];   Cop = [];   Dop = [];
elseif nargin==3
    % Interpret three arguments as system with disturbance,
    %   d/dt(Top*v) = Aop*v + Bop*w
    Cop = [];   Dop = [];
elseif nargin==4
    % Interpret three arguments as system with disturbance and output,
    %   d/dt(Top*v) = Aop*v + Bop*w
    %             z = Cop*v
    Dop = [];
elseif nargin>5
    error('Too many input arguments.')
end

% % Extract operators defining the PIE.
% Operators Top s.t u = Top*v + Twop*w + Tuop*u
Twop = [];      Tuop = [];
if isa(Top,'cell')
    if isscalar(Top)
        Top = Top{1};
    elseif numel(Top)==2
        Twop = Top{2};
        Top = Top{1};
    elseif numel(Top)==3
        Tuop = Top{3};
        Twop = Top{2};
        Top = Top{1};
    else
        error("Operators Top defining the fundamental map should be defined as a single 'opvar' object or 1x3 cell {Top,Twop,Tuop} of 'opvar' objects.")
    end
end
if isequal(Top,0) || isempty(Top)
    Top = [];
elseif ~isa(Top,'opvar') && ~isa(Top,'opvar2d') && ~isa(Top,'polynomial') && ~isa(Top,'double')
    error("Operators Top defining the fundamental map should be defined as a single 'opvar' object or 1x3 cell {Top,Twop,Tuop} of 'opvar' objects.")
end
if (isa(Twop,'double') && isequal(Twop,0)) || isempty(Twop)
    Twop = [];
elseif ~isa(Twop,'opvar') && ~isa(Twop,'opvar2d') && ~isa(Twop,'polynomial') && ~isa(Twop,'double')
    error("Operators Top defining the fundamental map should be defined as a single 'opvar' object or 1x3 cell {Top,Twop,Tuop} of 'opvar' objects.")
end
if (isa(Tuop,'double') && isequal(Tuop,0)) || isempty(Tuop)
    Tuop = [];
elseif ~isa(Tuop,'opvar') && ~isa(Tuop,'opvar2d') && ~isa(Tuop,'polynomial') && ~isa(Tuop,'double')
    error("Operators Top defining the fundamental map should be defined as a single 'opvar' object or 1x3 cell {Top,Twop,Tuop} of 'opvar' objects.")
end
% Operator Aop
if isa(Aop,'cell')
    if ~isscalar(Aop)
        error("Operator Aop defining the PIE dynamics should be defined as a single 'opvar' object.")
    else
        Aop = Aop{1};
    end
end
if (isa(Aop,'double') && isequal(Aop,0)) || isempty(Aop)
    Aop = [];
elseif ~isa(Aop,'opvar') && ~isa(Aop,'opvar2d') && ~isa(Aop,'polynomial') && ~isa(Aop,'double')
    error("Operator Aop defining the PIE dynamics should be defined as a single 'opvar' object.")
end
% Operator Bop
Bwop = Bop;     Buop = [];      % assume no controlled input by default
if isa(Bop,'cell')
    if isscalar(Bop)
        Bwop = Bop{1};
    elseif numel(Bop)==2
        Bwop = Bop{1};
        Buop = Bop{2};
    else
        error("Input operators Bop should be defined as a single 'opvar' object Bwop or 1x2 cell {Bwop,Buop} of 'opvar' objects.")
    end
end
if (isa(Bwop,'double') && isequal(Bwop,0)) || isempty(Bwop)
    Bwop = [];
elseif ~isa(Bwop,'opvar') && ~isa(Bwop,'opvar2d') && ~isa(Bwop,'polynomial') && ~isa(Bwop,'double')
    error("Input operators Bop should be defined as a single 'opvar' object Bwop or 1x2 cell {Bwop,Buop} of ''opvar'' objects.")
end
if (isa(Buop,'double') && isequal(Buop,0)) || isempty(Buop)
    Buop = [];
elseif ~isa(Buop,'opvar') && ~isa(Buop,'opvar2d') && ~isa(Buop,'polynomial') && ~isa(Buop,'double')
    error("Input operators Bop should be defined as a single 'opvar' object Bwop or 1x2 cell {Bwop,Buop} of ''opvar'' objects.")
end
% Operator Cop
Czop = Cop;     Cyop = [];      % assume no controlled input by default
if isa(Cop,'cell')
    if isscalar(Cop)
        Czop = Cop{1};
    elseif numel(Cop)==2
        Czop = Cop{1};
        Cyop = Cop{2};
    else
        error("Output operators Cop should be defined as a single 'opvar' object Czop or 2x1 cell {Czop;Cyop} of 'opvar' objects.")
    end
end
if (isa(Czop,'double') && isequal(Czop,0)) || isempty(Czop)
    Czop = [];
elseif ~isa(Czop,'opvar') && ~isa(Czop,'opvar2d') && ~isa(Czop,'polynomial') && ~isa(Czop,'double')
    error("Output operators Cop should be defined as a single 'opvar' object Czop or 2x1 cell {Czop;Cyop} of 'opvar' objects.")
end
if (isa(Cyop,'double') && isequal(Cyop,0)) || isempty(Cyop)
    Cyop = [];
elseif ~isa(Cyop,'opvar') && ~isa(Cyop,'opvar2d') && ~isa(Cyop,'polynomial') && ~isa(Cyop,'double')
    error("Output operators Cop should be defined as a single 'opvar' object Czop or 2x1 cell {Czop;Cyop} of 'opvar' objects.")
end
% Operator Dop
Dzwop = Dop;    Dzuop = [];
Dywop = [];     Dyuop = [];
if isa(Dop,'cell')
    if isscalar(Dop)
        Dzwop = Dop{1};
    elseif all(size(Dop)==[1,2])
        Dzwop = Dop{1};     Dzuop = Dop{2};
    elseif all(size(Dop)==[2,1])
        Dzwop = Dop{1};     Dywop = Dop{2};
    elseif all(size(Dop)==[2,2])
        Dzwop = Dop{1,1};   Dzuop = Dop{1,2};
        Dywop = Dop{2,1};   Dyuop = Dop{2,2};
    else
        error("Feedthrough operators Dop should be declared as a single 'opvar' object Dzwop or 2x2 cell {Dzwop,Dzuop;Dywop,Dyuop} of 'opvar' objects.")
    end
elseif isa(Dop,'opvar') || isa(Dop,'opvar2d') || isa(Dop,'double') || isa(Dop,'polynomial')
    % If a single operator is specified, check if only a single input or
    % output operator is specified as well.
    if isempty(Bwop) && ~isempty(Buop)
        if isempty(Czop) && ~isempty(Cyop)
            Dyuop = Dop;    Dzwop = [];
        else
            if ~isempty(Cyop)
                warning("Assuming feedthrough operator maps to regulated output y.")
            end
            Dzuop = Dop;    Dzwop = [];
        end
    else 
        if ~isempty(Buop)
            warning("Assuming feedthrough operator maps disturbance w.")
        end
        if isempty(Czop) && ~isempty(Cyop)
            Dywop = Dop;    Dzwop = [];
        elseif ~isempty(Cyop)
            warning("Assuming feedthrough operator maps to regulated output y.")
        end
    end
else
    error("Feedthrough operators Dop should be declared as a single 'opvar' object Dzwop or 2x2 cell {Dzwop,Dzuop;Dywop,Dyuop} of 'opvar' objects.")
end
if isa(Dzwop,'double') && (isempty(Dzwop) || isequal(Dzwop,0))
    Dzwop = [];
elseif ~(isa(Dzwop,'opvar') || isa(Dzwop,'opvar2d') || isa(Dzwop,'double') || isa(Dzwop,'polynomial'))
    error("Feedthrough operators Dop should be declared as a single 'opvar' object Dzwop or 2x2 cell {Dzwop,Dzuop;Dywop,Dyuop} of 'opvar' objects.")
end
if isa(Dzuop,'double') && (isempty(Dzuop) || isequal(Dzuop,0))
    Dzuop = [];
elseif ~(isa(Dzuop,'opvar') || isa(Dzuop,'opvar2d') || isa(Dzuop,'double') || isa(Dzuop,'polynomial'))
    error("Feedthrough operators Dop should be declared as a single 'opvar' object Dzwop or 2x2 cell {Dzwop,Dzuop;Dywop,Dyuop} of 'opvar' objects.")
end
if isa(Dywop,'double') && (isempty(Dywop) || isequal(Dywop,0))
    Dywop = [];
elseif ~(isa(Dywop,'opvar') || isa(Dywop,'opvar2d') || isa(Dywop,'double') || isa(Dywop,'polynomial'))
    error("Feedthrough operators Dop should be declared as a single 'opvar' object Dzwop or 2x2 cell {Dzwop,Dzuop;Dywop,Dyuop} of 'opvar' objects.")
end
if isa(Dyuop,'double') && (isempty(Dyuop) ||  isequal(Dyuop,0))
    Dyuop = [];
elseif ~(isa(Dyuop,'opvar') || isa(Dyuop,'opvar2d') || isa(Dyuop,'double') || isa(Dyuop,'polynomial'))
    error("Feedthrough operators Dop should be declared as a single 'opvar' object Dzwop or 2x2 cell {Dzwop,Dzuop;Dywop,Dyuop} of 'opvar' objects.")
end


% % % Make sure dimensions, domain, and variables of operators match.
% % Extract the dimensions of the operators.
[nr_T,nc_T,vars_T,dom_T] = get_op_props(Top);
if any(nr_T~=nc_T)
    error("Operator Top defining the fundamental map should be symmetric.")
end
[nr_A,nc_A,vars_A,dom_A] = get_op_props(Aop);
if any(nr_A~=nc_A)
    error("Operator Top defining the fundamental map should be symmetric.")
end
[nr_Tw,nc_Tw,vars_Tw,dom_Tw] = get_op_props(Twop);
[nr_Tu,nc_Tu,vars_Tu,dom_Tu] = get_op_props(Tuop);
[nr_Bw,nc_Bw,vars_Bw,dom_Bw] = get_op_props(Bwop);      
[nr_Bu,nc_Bu,vars_Bu,dom_Bu] = get_op_props(Buop);
[nr_Cz,nc_Cz,vars_Cz,dom_Cz] = get_op_props(Czop);      
[nr_Cy,nc_Cy,vars_Cy,dom_Cy] = get_op_props(Cyop);
[nr_Dzw,nc_Dzw,vars_Dzw,dom_Dzw] = get_op_props(Dzwop);   
[nr_Dzu,nc_Dzu,vars_Dzu,dom_Dzu] = get_op_props(Dzuop);
[nr_Dyw,nc_Dyw,vars_Dyw,dom_Dyw] = get_op_props(Dywop);   
[nr_Dyu,nc_Dyu,vars_Dyu,dom_Dyu] = get_op_props(Dyuop);

% % Determine dimensions of the state, inputs, and outputs
nv_op = check_dims(nr_T,nr_Tw,nr_Tu,nr_A,nr_Bw,nr_Bu,nc_Cz,nc_Cy);
nw_op = check_dims(nc_Tw,nc_Bw,nc_Dzw,nc_Dyw);
nu_op = check_dims(nc_Tu,nc_Bu,nc_Dzu,nc_Dyu);
nz_op = check_dims(nr_Cz,nr_Dzw,nr_Dzu);
ny_op = check_dims(nr_Cy,nr_Dyw,nr_Dyu);
if any(isnan(nv_op))
    error("State dimensions of the specified operators are not consistent.")
end
if any(isnan(nw_op))
    error("Exogenous input dimensions of the specified operators are not consistent.")
end
if any(isnan(nu_op))
    error("Controlled input dimensions of the specified operators are not consistent.")
end
if any(isnan(nz_op))
    error("Regulated output dimensions of the specified operators are not consistent.")
end
if any(isnan(ny_op))
    error("Observed output dimensions of the specified operators are not consistent.")
end
% Determine the dimension of the spatial domain of functions
dim = log(max([size(nv_op,1),size(nw_op,1),size(nu_op,1),size(nz_op,1),size(ny_op,1)]))/log(2);

% % Determine the spatial variables and domain of the system
if dim>0
    vars = check_vars({vars_T,dom_T},{vars_Tw,dom_Tw},{vars_Tu,dom_Tu},...
                        {vars_A,dom_A},{vars_Bw,dom_Bw},{vars_Bu,dom_Bu},...
                         {vars_Cz,dom_Cz},{vars_Dzw,dom_Dzw},{vars_Dzu,dom_Dzu},...
                          {vars_Cy,dom_Cy},{vars_Dyw,dom_Dyw},{vars_Dyu,dom_Dyu});
    if ~isa(vars,'cell') && isnan(vars)
        error("Spatial variables and domains of operators must match.")
    end
    var1 = vars{1}(:,1);
    var2 = vars{1}(:,2);
    dom = vars{2};
    if any(isequal(var2,0))
        warning("Secondary spatial variables of the operators do not match; replacing with default variables (s1_dum,s2_dum,...)")
        var2 = polynomial(zeros(dim,1));
        for ii=1:dim
            var2(ii) = polynomial({[var1(ii).varname{1},'_dum']});
        end
    end
else
    % Aparantly all operators are finite-dimensional
    % --> augment to 1D system, with default domain and variables;
    tmp = opvar();
    var1 = tmp.var1;    var2 = tmp.var2;    dom = tmp.I;
end

% % Convert all operators to same class and on same domain,
% % and using same variables.
Top = rebuild_op(Top,nv_op,nv_op,[var1,var2],dom);
Twop = rebuild_op(Twop,nv_op,nw_op,[var1,var2],dom);
Tuop = rebuild_op(Tuop,nv_op,nu_op,[var1,var2],dom);
Aop = rebuild_op(Aop,nv_op,nv_op,[var1,var2],dom);
Bwop = rebuild_op(Bwop,nv_op,nw_op,[var1,var2],dom);
Buop = rebuild_op(Buop,nv_op,nu_op,[var1,var2],dom);
Czop = rebuild_op(Czop,nz_op,nv_op,[var1,var2],dom);
Dzwop = rebuild_op(Dzwop,nz_op,nw_op,[var1,var2],dom);
Dzuop = rebuild_op(Dzuop,nz_op,nu_op,[var1,var2],dom);
Cyop = rebuild_op(Cyop,ny_op,nv_op,[var1,var2],dom);
Dywop = rebuild_op(Dywop,ny_op,nw_op,[var1,var2],dom);
Dyuop = rebuild_op(Dyuop,ny_op,nu_op,[var1,var2],dom);


% % % Build the PIE representation.
pie_struct PIE;
PIE.vars = [var1,var2];     PIE.dom = dom;

PIE.T = Top;    PIE.Tw = Twop;      PIE.Tu = Tuop;
PIE.A = Aop;    PIE.Bw = Bwop;      PIE.Bu = Buop;
PIE.Cz = Czop;  PIE.Dzw = Dzwop;    PIE.Dzu = Dzuop;
PIE.Cy = Cyop;  PIE.Dyw = Dywop;    PIE.Dyu = Dyuop;

end




%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [nr_op,nc_op,vars,dom] = get_op_props(Pop)
% [NR_OP,NC_OP,VARS,DOM] = GET_OP_PROPS(POP) takes an object Pop 
% representing a PI operators, and returns the row dimension nr_op, column 
% dimension nc_op, spatial domain dom and spatial variables vars of the
% operator.
%
% INPUT
% - Pop:    Object of type 'opvar', 'opvar2d', 'polynomial', or 'double',
%           of which to return the properties.
%
% OUTPUT
% - nr_op:  Numeric array indicating the row dimension of Pop. If Pop
%           is of type 'double' or 'polynomial', then nr_op is a scalar,
%           given by nr_op=size(Pop,1). If Pop is of type 'opvar', or
%           'opvar2d', then nr_op is a 2x1 or 4x1 array, respectively,
%           given by Pop.dim(:,1).
% - nc_op:  Numeric array indicating the column dimension of Pop. If Pop
%           is of type 'double' or 'polynomial', then nc_op is a scalar,
%           given by nc_op=size(Pop,2). If Pop is of type 'opvar', or
%           'opvar2d', then nc_op is a 2x1 or 4x1 array, respectively,
%           given by Pop.dim(:,2).
% - vars:   nx2 array of type 'polynomial', specifying the primary spatial
%           variables (first column) and dummy variables (second column) in
%           which the operator is defined. If Pop is of type 'double' or
%           'polynomial', vars = zeros(0,2);
% - dom:    nx2 array indicating the spatial domain of functions on which
%           the operator acts. If Pop is of type 'double' or 'polynomial',
%           then dom = zeros(0,2). Otherwise, if Pop is of type 'opvar' or
%           'opvar2d', then dom is of dimension 1x2 or 2x2, respectively.
%
% NOTES
% If Pop is an empty array or equal to the scalar 0, it is assumed that the
% dimensions of Pop are not necessarily fixed, and nr_op=nc_op=zeros(0,1)
% is returned.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 12/19/2024: Initial coding;

if isa(Pop,'opvar') || isa(Pop,'opvar2d')
    nr_op = Pop.dim(:,1);
    nc_op = Pop.dim(:,2);
    dom = Pop.I;
    vars = [Pop.var1,Pop.var2];
elseif ~isempty(Pop)
    nr_op = size(Pop,1);
    nc_op = size(Pop,2);
    dom = zeros(0,1);
    vars = polynomial(zeros(0,2));
else
    nr_op = zeros(0,1);
    nc_op = zeros(0,1);
    dom = zeros(0,1);
    vars = polynomial(zeros(0,2));
end
end


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function n_op = check_dims(varargin)
% N_OP = CHECK_DIMS(n_op1,n_op2,...) takes a list of row or column
% dimensions of operators, and checks whether these dimensions match,
% returning a new operator dimension n_op equal to all specified dimensions
% or nan if dimensions don't match.
%
% INPUT
% - varargin:   A list of arrays (at least 1) of size 1x1, 2x1, or 4x1,
%               representing the row or column dimensions of matrices 
%               (1x1), 'opvar' objects (2x1), or 'opvar2d' objects (4x1).
%               Elements of each arrays must be real nonnegative integers, 
%               but different arrays may have different sizes. Arryas
%               may also be empty to indicate an unknown dimension, which
%               will be skipped. 
%
% OUTPUT
% - n_op:   CASE1: If all input arguments could feasibly represent the
%           row or column dimension of the same operator, then n_op will be
%           a numeric array of dimension nx1 indicating this row or column
%           dimension. Here, n will be the maximal number of elements of
%           all input arrays, representing the dimension of a matrix (n=1),
%           'opvar' object (n=2), or 'opvar2d' object (n=4), noting that
%           an 'opvar' object could be used to represent a matrix, and a
%           'opvar2d' object could be used to represent a 'opvar' object.
%           CASE2: If any of the input arguments defines a row or column
%           dimension that could not match that specified by the other
%           input arguments, n_op=nan is returned.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 12/19/2024: Initial coding;

if nargin==0
    error("Insufficient input arguments.")
elseif nargin==1
    n_op = varargin{1};
    return
else
    n_op1 = varargin{1};
    n_op2 = varargin{2};
end

% Check whether dimensions n_op2 match earlier value n_op1
if isempty(n_op1) && ~isempty(n_op2)
    n_op1 = n_op2;
elseif ~isempty(n_op1) && ~isempty(n_op2)
    if isscalar(n_op2) 
        if n_op2~=sum(n_op1)
            % Pass nan to indicate row dimensions do not match.
            n_op = nan;
            return
        end
    elseif numel(n_op2)==2
        % Dimensions correspond to 1D operator
        if isscalar(n_op1)
            if sum(n_op2)~=n_op1
                % Pass nan to indicate row dimensions do not match.
                n_op = nan;
                return
            else
                % Set current dimension as new default.
                n_op1 = n_op2;
            end
        elseif numel(n_op1)==2 && any(n_op2~=n_op1)
            % Pass nan to indicate row dimensions do not match.
            n_op = nan;
            return
        elseif numel(n_op1)==4
            if nnz(n_op1)~=nnz(n_op2) || any(n_op1(n_op1~=0)~=n_op2(n_op2~=0))
                % Pass nan to indicate row dimensions do not match.
                n_op = nan;  
                return
            end
        end
    elseif numel(n_op2)==4
        % Dimensions correspond to 2D operator
        if isscalar(n_op1)
            if sum(n_op2)~=n_op1
                % Pass nan to indicate row dimensions do not match.
                n_op = nan;
                return
            else
                % Set current dimension as new default.
                n_op1 = n_op2;
            end
        elseif numel(n_op1)==2
            if nnz(n_op1)~=nnz(n_op2) || any(n_op1(n_op1~=0)~=n_op2(n_op2~=0))
                % Pass nan to indicate row dimensions do not match.
                n_op = nan; 
                return
            else
                % Set current dimension as new default.
                n_op1 = n_op2;
            end
        elseif any(n_op2~=n_op1)
            n_op = nan;
            return
        end
    end
end

% If n_op==nan, the array elements already do not match, so there is no
% point in proceeding
if nargin>2
    n_op = check_dims(n_op1,varargin{3:end});
else
    n_op = n_op1;
end

end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function vars = check_vars(varargin)
% VARS = CHECK_DIMS(vars1,vars2,...) takes a list of arrays of variables
% used in defining PI operators, as well as the domains on which these
% variables exist, and checks whether these variables and domains match,
% returning a unique list of variables and domains on which all operators 
% are defined, or nan if variables don't match.
%
% INPUT
% - varargin:   A list of 1x2 cell arrays (at least 1), with the first
%               element of each cell being a 1x2 or 2x2 array of type
%               'polynomial', representing primary and dummy variables of 
%               'opvar' objects (1x1) or 'opvar2d' objects (2x1). The
%               second element of each cell should be a 1x2 or 2x2 array
%               specifying the lower limits (first column) and upper limits
%               (second column) of the domain in which each of the
%               variables exists. 
%
% OUTPUT
% - vars:   CASE1: If all input arguments could feasibly represent the
%           spatial or dummy variables of the same operator, then n_op will
%           be a 'polynomial' array of dimension nx1 indicating these
%           variables. Here, n will be the maximal number of elements of
%           all input arrays, representing the dimension of the space on
%           which the operator acts, being n=1 for 'opvar' objects' and
%           'n=2' for opvar2d objects.
%           CASE2: If any of the input arguments defines a variable
%           dimension that could not match that specified by the other
%           input arguments, vars=nan is returned.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 12/19/2024: Initial coding;

if nargin==0
    error("Insufficient input arguments.")
elseif nargin==1
    vars = varargin{1};
    return
else
    vars1 = varargin{1}{1};     dom1 = varargin{1}{2};
    vars2 = varargin{2}{1};     dom2 = varargin{2}{2};
end

% Check whether variables vars2 match earlier variables vars1
if isempty(vars1) && ~isempty(vars2)
    vars1 = vars2;
    dom1 = dom2;
elseif ~isempty(vars1) && ~isempty(vars2)
    if size(vars2,1)==1 
        if size(vars1,1)==1
            % Check that the domains match.
            if any(dom1~=dom2) || ~isequal(vars1(1),vars2(1))
                vars = nan;
                return
            end
            if ~isequal(vars1(2),vars2(2))
                % Use 0 to indicate that dummy variables don't match.
                vars1(1,2) = 0;
            end
        elseif size(vars1,1)==2
            idx = isequal(vars1(:,1),vars2(:,1));
            if any(dom1(idx,:)~=dom2)
                vars = nan;
                return
            end
            if ~isequal(vars1(idx,2),vars2(1,2))
                % Use 0 to indicate that dummy variables don't match.
                vars1(idx,2) = 0;
            end
        end
    elseif size(vars2,1)==2
        % Dimensions correspond to 2D operator
        if size(vars1,1)==1
            idx = isequal(vars1(:,1),vars2(:,1));
            if any(dom1~=dom2(idx,:))
                vars = nan;
                return
            end
            % Set current variables as new default.
            vars1 = vars2;
            dom1 = dom2;
            if ~isequal(vars2(idx,2),vars1(1,2))
                % Use 0 to indicate that dummy variables don't match.
                vars1(idx,2) = 0;
            end
        elseif size(vars1,1)==2
            % Make sure both variables are the same, and in the same order.
            if any(any(dom1~=dom2)) || ~all(isequal(vars1(:,1),vars2(:,1)))
                vars = nan;
                return
            end
            % Use 0 to indicate that dummy variables don't match.
            if any(~isequal(vars1(:,2),vars2(:,2)))
                vars1(~isequal(vars1(:,2),vars2(:,2)),2) = 0;
            end
        end
    end
end

% Repeat for remaining variables.
if nargin>2
    vars = check_vars({vars1,dom1},varargin{3:end});
else
    vars = {vars1,dom1};
end

end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Pop_out = rebuild_op(Pop_in,rdim,cdim,vars,dom)
% POP_OUT = REBUILD_OP(POP_IN,RDIM,CDIM,VARS,DOM) takes a matrix or 'opvar'
% object Pop_in, and converts it to an opvar object Pop_out of dimensions
% [rdim,cdim] and in variables vars on domain dom. If this is not possible,
% it throws an error.
%
% INPUT
% - Pop_in: object of type 'double', 'polynomial', 'opvar' or 'opvar2d',
%           representing some PI operator.
% - rdim:   nx1 array of type 'double' specifying the desired/expected row
%           dimension of the operator, where n=1 if Pop_in is a matrix, n=2
%           if 'opvar', and n=4 if 'opvar2d'. If set to [], the dimension
%           is assumed to be 0.
% - cdim:   nx1 array of type 'double' specifying the desired column
%           dimension of the operator, where n=1 if Pop_in is a matrix, n=2
%           if 'opvar', and n=4 if 'opvar2d'. If set to [], the dimension
%           is assumed to be 0.
% - vars:   Nx2 array of type 'polynomial' specifying the primary (first
%           column) and dummy (second column) variables desired/expected in
%           the PI operator. Here, N=0 for finite-dimensional systems,
%           N=1 for 1D systems, and N=2 for 2D systems.
% - dom:    Nx2 array of type 'double' specifying the lower (first column)
%           and upper (second column) limits of the intervals on which each
%           of the variables is defined.
%
% OUTPUT
% - Pop_out:    'opvar' (if N=1) or 'opvar2d' (if N=2) object of dimension
%           [rdim,cdim], representing the same operator as Pop_in, defined
%           in terms of variables vars on domain dom. If Pop_in is already
%           of type 'opvar' or 'opvar2d', the function merely checks that
%           the dimensions, variables, and domain defining Pop_in match the
%           expected values, throwing an error otherwise. If dummy
%           variables don't match, these are replaced with the expected
%           values (vars(:,2)). If Pop_in is of type 'double' or
%           'polynomial', it is augmented to an 'opvar' class object of
%           specified dimensions, if possible.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 12/19/2024: Initial coding;

% Check the dimensionality of the system.
N = size(vars,1);

% If the size of the object is not known, assume it is zero.
if isempty(rdim) || isequal(rdim,0)
    rdim = zeros(2^N,1);
end
if isempty(cdim) || isequal(cdim,0)
    cdim = zeros(2^N,1);
end

% Convert dimensions of matrix to those of an opvar
if isscalar(rdim)
    rdim_op = zeros(2^N,1);
    rdim_op(1) = rdim;
    rdim = rdim_op;
end
if isscalar(cdim)
    cdim_op = zeros(2^N,1);
    cdim_op(1) = cdim;
    cdim = cdim_op;
end

% Do not support conversion of opvars to opvar2d objects for now.
if (N==1 && isa(Pop_in,'opvar2d')) || (N==2 && isa(Pop_in,'opvar'))
    error("Inputs must all be of type 'opvar', or all be of type 'opvar2d'.")
end

if N==0
    % All operators are matrices; convert to 'opvar' objects.
    Pop_out = opvar();
    Pop_out.dim = [rdim,cdim;0,0];
    Pop_out.var1 = vars(:,1);   Pop_out.var2 = vars(:,2);
    Pop_out.I = dom;
    Pop_out.P = Pop_in;
elseif N==1
    % 1D PIE
    if isempty(Pop_in)
        % Build 0 opvar object.
        Pop_out = opvar();
        Pop_out.dim = [rdim,cdim];
        Pop_out.var1 = vars(:,1);   Pop_out.var2 = vars(:,2);
        Pop_out.I = dom;
    elseif isa(Pop_in,'opvar')
        % Make sure dimensions of the operator match.
        if any(any(Pop_in.dim~=[rdim,cdim]))
            error("Construction of the PIE failed: Dimensions of the operators don't match.")
        end
        % Make sure variables and domain match.
        if any(Pop_in.I~=dom)
            error("All PI operators must be defined on the same spatial domain.")
        elseif any(~isequal(Pop_in.var1,vars(:,1)))
            error("All PI operators must be defined in terms of the same spatial variables.")
        end
        Pop_out = Pop_in;
        % We don't care about dummy vars, replace them if necessary.
        if any(~isequal(Pop_in.var2,vars(:,2)))
            Pop_out = subs_vars_op(Pop_out,Pop_out.var2,vars(:,2));
        end
    else
        % Try to convert non-opvar object to opvar.
        try Pop_out = mat2opvar(Pop_in,[rdim,cdim],vars,dom);
        catch
            error("Construction of the PIE failed: Please specify all operators as objects of type 'opvar' or 'opvar2d'.")
        end
    end
elseif N==2
    % 2D PIE
    if isempty(Pop_in) || isequal(Pop_in,0)
        % Build 0 opvar2d object.
        Pop_out = opvar2d();
        Pop_out.dim = [rdim,cdim];
        Pop_out.var1 = vars(:,1);   Pop_out.var2 = vars(:,2);
        Pop_out.I = dom;
    elseif isa(Pop_in,'opvar2d')
        % Make sure dimensions of the operator match.
        if any(any(Pop_in.dim~=[rdim,cdim]))
            error("Construction of the PIE failed: Dimensions of the operators don't match.")
        end
        % Make sure variables and domain match.
        if any(any(Pop_in.I~=dom))
            error("All PI operators must be defined on the same spatial domain.")
        elseif any(~isequal(Pop_in.var1,vars(:,1)))
            error("All PI operators must be defined in terms of the same spatial variables.")
        end
        Pop_out = Pop_in;
        % We don't care about dummy vars, replace them if necessary.
        if any(~isequal(Pop_in.var2,vars(:,2)))
            Pop_out = subs_vars_op(Pop_out,Pop_out.var2,vars(:,2));
        end
    else
        % Try to convert non-opvar object to opvar.
        try Pop_out = mat2opvar(Pop_in,[rdim,cdim],vars,dom);
        catch
            error("Construction of the PIE failed: Please specify all operators as objects of type 'opvar' or 'opvar2d'.")
        end
    end
end

end