function out = getStatesFromEquations(obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys-class method that identifies unique states in equations list
% Input: 
% obj - sys class object
% Output:
% out - vector of state class objects
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
%
% Initial coding SS - 06/28/2022
% MMP, 11/19/2026: Release comment pass (comment-only).

% Collect the states of each equation, each state name once
equations = obj.equation;
eqnNum = length(equations);
if eqnNum==0
    out = [];
else
out = [];
for i=1:eqnNum
    tempterms = equations(i);
    out = combine_statename(out,tempterms.statevec);
end
end
end
function out = combine_statename(varargin)
% OUT = COMBINE_STATENAME(VARARGIN) appends to the states VARARGIN{1}
% each state of VARARGIN{2:end} whose name (statename) is not yet in them.
%
% INPUT
% - varargin:   'state' objects (vectors), the first possibly empty;
%
% OUTPUT
% - out:        column of 'state' objects with distinct names.
out = varargin{1};
for i=2:nargin
    tmp = varargin{i};
    for j=1:length(tmp)
        s.type = '()'; s.subs = {j};
        temp = subsref(tmp,s); 
        if isempty(out)
            out = temp;
        else
            if ~ismember(temp.statename,out.statename)
                out = [out; temp];
            end
        end
    end
end
end