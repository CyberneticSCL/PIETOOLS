function obj = subs(obj,old,new)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that performs replaces obj(old) by obj(new)
% Input: 
% obj - state class object, obj(old)
% old, new - pvar variable
% Output:
% obj - state class object, obj(new)
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
if nargin~=3
    error('Substitution requires 3 inputs: state object, variable to be substituted, value to be substituted.');
end
if ~isa(new,'polynomial')&&((new~=0)&&(new~=1))
    error('Subs operation can only be performed with pvar variable or at the boundaries 0 and 1.');
end
for i=1:length(obj)
    if isequal(obj(i).var(1),old)
        if ~poly2double(obj(i).var(1)-new)
            error('Subs on time can only be performed to add delays: must be of the form t-tau, where tau is positive real');
        end
        if poly2double(obj(i).var(1)-new) && double(obj(i).var(1)-new)<0
            error('Positive delays are not allowed since the system is non-causal');
        end
    end
%     if ismember(obj(i).type,{'ode','in','out'})
%         error('Subs operation on ODE states, inputs, and outputs is not supported')
%     end
    idx = find(isequal(obj(i).var,old));
    obj(i).var(idx) = new;
end
end