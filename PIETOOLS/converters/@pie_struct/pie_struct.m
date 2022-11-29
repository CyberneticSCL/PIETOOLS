classdef (InferiorClasses={?opvar2d,?dpvar,?polynomial}) pie_struct
% pie_struct is a class of objects used to describe PIE systems in
% PIETOOLS 2022. It has the following fields:
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
%
% - PIE.T, PIE.Tw, PIE.Tu:
%               opvar or opvar2d objects describing the map from a
%               fundamental state xf in the PIE to an associated state x in
%               the equivalent PDE representation, as:
%                   x = T*xf + Tw*w + Tu*u;
%               where w are exogenous inputs, and u are actuator inputs.
%               The operators are uniquely defined by the boundary
%               conditions in the PDE structure, taking the form
%                   0 = F(x,w,u).
%
% - PIE.A,  PIE.B1,  PIE.B2,
%   PIE.C1, PIE.D11, PIE.D12,
%   PIE.C2, PIE.D21, PIE.D22:      
%               opvar or opvar2d objects defining a Partial Integral
%               Equation (PIE) of the form
%
%       Tu * u + Tw * w + T * x = A  * x + B1  * w + B2  * u;
%                             z = C1 * x + D11 * w + D12 * u;
%                             y = C2 * x + D21 * w + D22 * u;
%
% - PIE.x_tab, PIE.w_tab, PIE.u_tab, PIE.z_tab, PIE.y_tab:
%
%               ncomp x (2+nvars) arrays providing information on the state
%               components, exogenous inputs, actuator inputs, regulated
%               outputs, and observed outputs in the PIE structure. 
%               For example, if the fundamental state xf consists of 
%               ncomp_x components, the first column of PIE.x_tab provides 
%               for each of these state components an index, indicating
%               which state component in the original PDE is associated to
%               this fundamental state component. 
%               The second column of PIE.x_tab provides the size of each 
%               state component.
%               Columns 3:2+nvars are binary values, indicating for each of
%               the nvars variables in the PIE whether the state component
%               depends on this spatial variable.
%               NOTE: PIE.x_tab has nvars extra columns (3+nvars:2+2*nvars)
%               indicating for each spatial variable to what order the PDE
%               state compopnent is differentiated wrt this variable in
%               defining the associated fundamental state component in the
%               PIE.

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
% Initial coding SS, DJ - 08/12/2022
    
    properties
        dom = zeros(0,2);
        vars = polynomial(zeros(0,2));
        T;
        Tw;
        Tu;
        A;
        B1;
        B2;
        C1
        D11;
        D12;
        C2;
        D21;
        D22;
    end
    
    properties (Hidden)
        % The tables with standardized information on each of the components
        % should not be specified by the user, but should be added to the
        % structure by the conversion function".
        x_tab = zeros(0,2); % [state_num, state_size, isvar1, isvar2, vardiff1, vardiff2]
        u_tab = zeros(0,2);
        w_tab = zeros(0,2);
        y_tab = zeros(0,2);
        z_tab = zeros(0,2);
    end

    properties (Dependent)
        % Spatial dimension of the PIE.
        dim = 0;
    end
    
methods
    function obj = pie_struct(varargin)
        %PIE_STRUCT Construct an instance of this class
        %   Detailed explanation goes here
        
        if nargin==1
            if ischar(varargin{1})
                if nargout==0
                    assignin('caller', varargin{1}, pie_struct());
                end
            elseif isa(varargin{1},'pie_struct')
                obj = varargin{1};
            elseif isa(varargin{1},'struct')
                % If a struct is provided, convert to a "pie_struct".
                obj = pie_struct();
                % Copy information from any field in the struct to the
                % pie_struct.
                if isfield(varargin{1},'vars')
                    obj.vars = varargin{1}.vars;
                end
                if isfield(varargin{1},'dom')
                    obj.dom = varargin{1}.dom;
                end
                
                if isfield(varargin{1},'T')
                    obj.T = varargin{1}.T;
                end
                if isfield(varargin{1},'Tu')
                    obj.Tu = varargin{1}.Tu;
                end
                if isfield(varargin{1},'Tw')
                    obj.Tw = varargin{1}.Tw;
                end
                if isfield(varargin{1},'A')
                    obj.A = varargin{1}.A;
                end
                if isfield(varargin{1},'B1')
                    obj.B1 = varargin{1}.B1;
                end
                if isfield(varargin{1},'B2')
                    obj.B2 = varargin{1}.B2;
                end
                if isfield(varargin{1},'C1')
                    obj.C1 = varargin{1}.C1;
                end
                if isfield(varargin{1},'D11')
                    obj.D11 = varargin{1}.D11;
                end
                if isfield(varargin{1},'D12')
                    obj.D12 = varargin{1}.D12;
                end
                if isfield(varargin{1},'C2')
                    obj.C2 = varargin{1}.C2;
                end
                if isfield(varargin{1},'D21')
                    obj.D21 = varargin{1}.D21;
                end
                if isfield(varargin{1},'D22')
                    obj.D22 = varargin{1}.D22;
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
            else
                error("Input must be strings");
            end
        else
            for i=1:nargin
                if ischar(varargin{i})
                    % Allow new pde_structs to be defined using
                    % > pie_struct PIE1 PIE2 PIE3 ...
                    if nargout==0
                        assignin('caller', varargin{i}, pie_struct());
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