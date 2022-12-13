classdef (InferiorClasses={?state}) sys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys class is the container class to store equations and parameters for
% pde, pie, time-delay systems and does not take any input arguments for
% initialization. To initialize, use X = sys(). The type of sys can be
% modified using the commands X.type = 'pde' or 'pie'. Likewise, the domain
% can be changed using the command X.dom = [0,5]. The other properties of
% system object can not be accessed or modified directly to avoid
% unintended errors.
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications

    properties
        equation {validateEquation(equation)} = [];
        type char {mustBeMember(type,{'pde','dde','ddf','pie'})} = 'pde';
        params {mustBeA(params,{'pde_struct','pie_struct'})} = pde_struct();
        ControlledInputs;
        ObservedOutputs;
        dom = [0,1];
    end
    properties (Dependent)
        states;
    end
    methods
        function obj = sys(type)
            if nargin==0
                type = 'pde';
            end
            obj.equation = [];
            obj.params = pde_struct();
            obj.ControlledInputs = [];
            obj.ObservedOutputs = [];
            obj.type = type;
            fprintf('Initialized sys() object of type "%s"\n',obj.type);
        end
        function prop = get.states(obj)
            prop = getStatesFromEquations(obj);
        end
        function prop = get.params(obj)
            if isempty(obj.params.dom)
                obj = getParams(obj);
                prop = obj.params;
            else
                prop = obj.params;
            end
        end
        function obj = set.params(obj,params)
            obj.params = params;
        end
        function obj = set.dom(obj,dom)
            obj.dom = dom;
        end
        function out = get.ObservedOutputs(obj)
            if isempty(obj.ObservedOutputs)
                statelist = getStatesFromEquations(obj);
                if ~isempty(statelist)
                    out = zeros(length(statelist.veclength),1);
                else
                    out = [];
                end
            else
                out = obj.ObservedOutputs;
            end
        end
        function obj = set.ObservedOutputs(obj,OO)
            obj.ObservedOutputs = OO;
        end
        function out = get.ControlledInputs(obj)
            if isempty(obj.ControlledInputs)
                statelist = getStatesFromEquations(obj);
                if ~isempty(statelist)
                    out = zeros(length(statelist.veclength),1);
                else
                    out = [];
                end
            else
               out = obj.ControlledInputs;
            end
        end
        function obj = set.ControlledInputs(obj,CI)
            obj.ControlledInputs = CI;
        end
        
        prop = getStatesFromEquations(obj);
        obj = addequation(obj,eqn);
        obj = getParams(obj);    
        obj = convert(obj,convertto);
        obj = removeequation(obj,eqnNumber);
        obj = setControl(obj,input);
        obj = removeControl(obj,input);
        obj = setObserve(obj,output);
        obj = removeObserve(obj,output);
        disp(obj);
    end
end
function validateEquation(prop)
% if ~iscell(prop)
%     error('Equations must be stored in cell column array');
% end
eqntype = (~isa(prop,'terms'))&& (~isempty(prop));
% eqntype = cell2mat(eqntype);
if ~isempty(eqntype)&&any(eqntype(:))
    error("Equation entries should be terms type object.");
end
end