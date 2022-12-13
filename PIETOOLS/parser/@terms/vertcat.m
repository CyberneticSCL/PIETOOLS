function obj = vertcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms-class method that performs [objA;objB] concatenation 
% Input: 
% varargin - terms class objects
% Output:
% obj - terms class object [varargin{1};...varargin{n}]
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

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
if nargin==1
    obj = varargin{1};
else
    objA = varargin{1};
    objB = varargin{2};
    
    if isa(objA,'state')
        objA = state2terms(objA);
    end
    if isa(objB,'state')
        objB = state2terms(objB);
    end
    
    if isempty(objA)
        obj = objB;
        return
    elseif isempty(objB)
        obj = objA;
        return
    end
    
    opvar zeroAB zeroBA;
    zeroAB.dim = [0 0; objA.operator.dim(2,1) objB.operator.dim(2,2)]; 
    zeroBA.dim = [0 0; objB.operator.dim(2,1) objA.operator.dim(2,2)];
    
    tempoperator = [objA.operator zeroAB; zeroBA objB.operator];
    tempstatevec = vertcat(objA.statevec, objB.statevec);
    
    obj = terms(tempoperator,tempstatevec);
    if nargin>2
        obj = vertcat(obj,varargin{3:end});
    end
end
end