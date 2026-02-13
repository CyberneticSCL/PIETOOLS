classdef (InferiorClasses={?polynomial,?dpvar,?nopvar,?ndopvar,?tensopvar})polyopvar
% This function defines the class of 'polyopvar' objects, representing
% polynomial functions on distributed states.
%
% CLASS properties
% - P.varname:  p x 1 cellstr, specifying the names of the p distributed
%               state variables on which the polynomial is defined,
%                   (x1,x2,...,xp);
% - P.degmat:   nZ x p array of integers, representing the basis of
%               distributed monomials in terms of which the polynomial is
%               defined. If P.degmat(i,:) = [d1,...,dp], then the ith
%               distributed monomial is defined by
%                   [Z(x)]_{i} = x1^d1 o ... o xp^dp,
%               where o denotes the tensor product, and where
%                   xi^di = xi o xi o ... o xi } <-- di products
% - P.C:        m x n 'tensopvar' object, representing the coefficient
%               operators acting on the distributed monomial basis, so that
%                   P(x) = C*(Z(x) o I_{n})
%               where I_{n} is the n x n identity matrix;
% - P.pvarname: N x 1 cellstr, specifying the names of the N independent
%               variables (s1,...,sN) in which the distributed state 
%               variables are defined;
% - P.dom:      N x 2 array, with row i specifying the interval [a,b] on
%               which indepenednet variable si is defined;
% - P.varmat:   p x N boolean array, with the element in row i and column j
%               indicating whether state variable xi depends on independent
%               variable sj
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - polyopvar
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
% DJ, 01/15/2026: Initial coding

    properties
        C = tensopvar();
        degmat = zeros(0,0);
        varname = cell(0,1);
        pvarname = cell(0,1);
        varmat = zeros(0,0);
        dom = zeros(0,2);
    end

    properties(Dependent)
        matdim;
    end

    methods

        function [P] = polyopvar(varargin) %constructor
            if nargout==0
                % Declare state variables as
                %   polyopvar x1(s1,s2) x2(s1) x3(s3)
                for i=1:nargin
                    if ischar(varargin{i})
                        obj = varargin{i};
                        % Split "x1(s1,s2)" into "x1" and "(s1,s2)"
                        is_bracket = ismember(obj,'(');
                        if any(is_bracket)
                            if sum(is_bracket)>1 || ~strcmp(obj(end),')')
                                error("Distributed polynomial variables must be specified as e.g. 'x1(s1,s2)'")
                            end
                            % Determine the name of the state variable
                            xname = obj(1:find(is_bracket,1)-1);
                            % Determine the names of the indepdnent
                            % variables
                            pnames = obj(find(is_bracket,1)+1:end-1);
                            pnames = strsplit(pnames,{',',' '});
                            % Declare the independent variables as pvar
                            % objects in the workspace
                            for j=1:numel(pnames)
                                pname = pnames{j};
                                if ~isnan(str2double(pname))
                                    error("Independent variable names must include a letter.")
                                else
                                    assignin('caller',pname,pvar(pname));
                                end
                            end
                            % Declare the state variable as polyopvar
                            % object in the workspace
                            assignin('caller', xname, polyopvar(xname,pnames'));
                        else
                            assignin('caller', obj, polyopvar(obj));
                        end
                    else
                        error("Input must be strings");
                    end
                end
            elseif nargout==1
                % Declare state variable as
                %   x = polyopvar('x',s1,[0,1]);
                if nargin==0
                    return
                end
                varname = varargin{1};
                if ischar(varname) && size(varname,1)==1
                    varname = {varname};
                elseif ~iscellstr(varname)
                    error("Variable names must be specified as string of char objects.")
                end
                p = numel(varname);
                % Check which spatial variables x depends on
                if nargin<=1
                    pvarname = {};
                else
                    if iscellstr(varargin{2})
                        pvarname = varargin{2};
                    elseif isa(varargin{2},'polynomial')
                        if ~ispvar(varargin{2})
                            error("Spatial variables must be specified as array of 'pvar' objects.")
                        end
                        pvarname = cell(numel(varargin{2}),1);
                        for jj=1:numel(pvarname)
                            pvarname(jj) = varargin{2}(jj).varname;
                        end
                    else
                        error("Spatial variables must be specified as array of 'pvar' objects.")
                    end
                end
                N = numel(pvarname);
                % Check the domain on which the spatial variables exist
                if nargin<=2
                    dom = [zeros(N,1),ones(N,1)];
                else
                    if all(size(varargin{3})==[1,2])
                        dom = repmat(varargin{3},[N,1]);
                    elseif all(size(varargin{3})==[N,2])
                        dom = varargin{3};
                    else
                        error("The spatial domain should be specified as px2 array for 'p' variables.")
                    end
                end
                % % Declare an identity operator mapping x to x
                % C = nopvar();
                % C.dom = dom;
                % C.vars = polynomial(pvarname);
                % C.deg = zeros(N,1);
                % C.C = cell([3*ones(1,N),1]);
                % C.C{1} = 1;
                % Declare the distributed polynomial variable
                %P.C = tensopvar(C);
                P.degmat = eye(p);
                P.varname = varname;
                P.pvarname = pvarname;
                P.varmat = true(p,N);
                P.dom = dom;
            end
        end

        function matdim = get.matdim(obj)
            % % Determine the dimensions of the matrix-valued distributed
            % % polynomial
            matdim = obj.C.matdim;
        end

    end
end