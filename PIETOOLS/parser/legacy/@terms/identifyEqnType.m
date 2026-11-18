function out = identifyEqnType(equation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms-class method that identifies if equation has time derivatives,
% outputs, or just boundary values
% Input: 
% equation - terms class objects
% Output:
% out - string, 'ode', 'pde', 'out', 'bc'
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
statevec = equation.statevec;
classifiedeqns = zeros(length(equation),1);
for i=1:length(statevec)
    if strcmp(statevec(i).type,'ode')
    elseif strcmp(statevec(i).type,'pde')
    elseif strcmp(statevec(i).type,'out')
    end
end

if any(strcmp(statevec.state.type,'out')) %if equation has output
    out = 'out';
elseif any(cellfun(@(x) x(1)~=0, statevec.diff,'un',1))%if term has derivative of time
    if strcmp(statevec.state(find(cellfun(@(x) x(1)~=0, statevec.diff,'un',1))).type,'ode')% ode with derivative of time
        out = 'ode';
    else %pde
        out = 'pde';
    end
else % boundary condition
    out = 'bc';
end
end