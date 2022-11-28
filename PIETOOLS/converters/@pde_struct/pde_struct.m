classdef (InferiorClasses={?opvar2d,?dpvar,?polynomial}) pde_struct
% pde_struct is a class of objects used to describe PDE systems in
% PIETOOLS 2022a. It has the following fields:
%
% - PDE.dim:    Integer "p" specifying the dimension of the
%               spatial domain of the system. Currently, only p<=2 is
%               supported.
% - PDE.vars:   A px2 pvar (polynomial class) object specifying
%               the p spatial variables that appear in the system in the
%               first column, and the associated dummy variables in the
%               second column.
% - PDE.dom:    A px2 "double" array specifying the interval on
%               which each of p spatial variables exists.
% - PDE.tau:    A qx2 "polynomial" array specifying temporal variables in
%               the first column, and the maximal value these delays can
%               assume in the second column.
%
% - PDE.x:      A cell with each element i defining a state 
%               component x_i and associated differential equation
%               \dot{x}_i = ....
% - PDE.y:      A cell with each element i defining an observed
%               output y_i and associated equation y_i = .... If not
%               specified, it is assumed that the system has no observed
%               outputs.
% - PDE.z:      A cell with each element i defining a regulated
%               output z_i and associated equation z_i = .... If not
%               specified, it is assumed that the system has no regulated
%               outputs.
% - PDE.u:      A cell with each element i defining an actuator
%               input u_i that appears in the system. If not specified, it
%               is assumed that the system has no actuator inputs.
% - PDE.w:      A cell with each element i defining an exogenous
%               input w_i that appears in the system. If not specified, it
%               is assumed that the system has no exogenous inputs.
% - PDE.BC:     A cell with each element i defining a boundary
%               condition 0 = ... of the system.
%
% For more information on how each of these objects is defined, see the
% "initialize_PIETOOLS_PDE" function.
%
% In addition, a pde_struct object also has hidden fields x_tab, y_tab,
% z_tab, u_tab, w_tab, and BC_tab. These are tables of which each row
% corresponds to a particular component specified in the PDE (e.g.
% x_tab(ii,:) corresponds to PDE.x{ii}). The first column provides an index
% associated to the considered component. The second column provides the
% size (e.g. number of state variables or inputs) associated to this
% component. The remaining columns 2+j for j=1:p are binary indices,
% setting e.g. x_tab(ii,2+j)=1 if state component x{ii} depends on the jth
% spatial variable in PDE.vars. For the x_tab, p additional columns 2+p+j
% for j=1:p provide the order up to which each state component is
% differentiable with each spatial variable it depends on, where always
% x_tab(ii,2+p+j)=0 if x_tab(ii,2+j)=0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 07/08/2022

    
properties
    % Primary and dummy spatial variables that appear in the system.
    vars = polynomial(zeros(0,2));
    % Interval on which each spatial variable exists.
    dom = zeros(0,2);
    % Delay variables (1st column) and associated delays (2nd column).
    tau = polynomial(zeros(0,2));
    
    x = cell(0);
    u = cell(0);
    w = cell(0);
    y = cell(0);
    z = cell(0);
    BC = cell(0);
    
end
properties (Hidden)
    % The tables with standardized information on each of the components
    % should not be specified by the user, but should be added to the
    % structure in "initialize_PIETOOLS_PDE".
    x_tab = zeros(0,2); % [state_num, state_size, isvar1, isvar2, vardiff1, vardiff2]
    u_tab = zeros(0,2);
    w_tab = zeros(0,2);
    y_tab = zeros(0,2);
    z_tab = zeros(0,2);
    BC_tab = zeros(0,2);
    
    is_initialized = false; % Has the system already been initialized.
    has_delay = false;      % Does the system involve any delayed states/inputs?
    has_hotd = false;       % Does the system involve any higher-order temporal derivatives?
end
properties (Dependent)
    % Dimension of the system
    dim = 0;
end

    
methods
    function obj = pde_struct(varargin)
        %PDE_STRUCT Construct an instance of this class
        %   Detailed explanation goes here
        
        if nargin==1
            if ischar(varargin{1})
                if nargout==0
                    assignin('caller', varargin{1}, pde_struct());
                end
            elseif isa(varargin{1},'pde_struct')
                obj = varargin{1};
            elseif isa(varargin{1},'struct')
                % If a struct is provided, convert to a "pde_struct".
                obj = pde_struct();
                % Copy information from any field in the struct to the
                % pde_struct.
                if isfield(varargin{1},'vars')
                    obj.vars = varargin{1}.vars;
                end
                if isfield(varargin{1},'dom')
                    obj.dom = varargin{1}.dom;
                end
                if isfield(varargin{1},'tau')
                    obj.dom = varargin{1}.tau;
                end
                
                if isfield(varargin{1},'x')
                    obj.x = varargin{1}.x(:);
                end
                if isfield(varargin{1},'u')
                    obj.u = varargin{1}.u(:);
                end
                if isfield(varargin{1},'w')
                    obj.w = varargin{1}.w(:);
                end
                if isfield(varargin{1},'y')
                    obj.y = varargin{1}.y(:);
                end
                if isfield(varargin{1},'z')
                    obj.z = varargin{1}.z(:);
                end
                if isfield(varargin{1},'BC')
                    obj.BC = varargin{1}.BC(:);
                end
                
                if isfield(varargin{1},'x_tab')
                    obj.x_tab = varargin{1}.x_tab;
                end
                if isfield(varargin{1},'u_tab')
                    obj.u_tab = varargin{1}.u_tab;
                end
                if isfield(varargin{1},'w_tab')
                    obj.w_tab = varargin{1}.w_tab;
                end
                if isfield(varargin{1},'y_tab')
                    obj.y_tab = varargin{1}.y_tab;
                end
                if isfield(varargin{1},'z_tab')
                    obj.z_tab = varargin{1}.z_tab;
                end
                if isfield(varargin{1},'BC_tab')
                    obj.BC_tab = varargin{1}.BC_tab;
                end
            else
                error("Input must be strings");
            end
        else
            for i=1:nargin
                if ischar(varargin{i})
                    % Allow new pde_structs to be defined using
                    % > pde_struct PDE1 PDE2 PDE3 ...
                    if nargout==0
                        assignin('caller', varargin{i}, pde_struct());
                    end
                else
                    error("Input must be strings");
                end
            end
        end
    end

    function dim_val = get.dim(obj)
        % Get spatial dimension of PDE.
        if size(obj.vars,1)==size(obj.dom,1)
            % Each variable should be assigned a domain.
            dim_val = size(obj.vars,1);
        else
            dim_val = nan;
        end
    end
end

end