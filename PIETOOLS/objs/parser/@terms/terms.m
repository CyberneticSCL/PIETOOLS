classdef (InferiorClasses={?state})terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an internal-use only class that stores equations in the form of
% P*x = 0, where P is PI operator (opvar class object) and x is a
% state-vector (state class object). There is no way to initialize or
% modify this class object directly.
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

properties (SetAccess=protected)
    operator opvar;
    statevec state;
end
methods
    function obj = terms(varargin)
        if nargin>2
            error('Terms object constructor takes only two inputs');
        end
        if nargin>1
            obj.statevec = varargin{2};
        end
        if nargin>0
            obj.operator = varargin{1};    
        end
    end
    
%     obj = delta(objA,var,val);
%     obj = diff(objA,var,order);
%     obj = horzcat(varargin);
    obj = int(objA,var,lim);
    obj = length(obj);
    obj = minus(objA,objB);
    obj = eq(objA,objB);
    obj = isequal(objA,objB);
    obj = mtimes(objA,objB); % note only one of two inputs can be terms object
    obj = plus(objA,objB);
%     obj = times(objA,objB); % note only one of two inputs can be terms object
    obj = uplus(objA);
    obj = uminus(objA);
    obj = vertcat(varargin);
end
end