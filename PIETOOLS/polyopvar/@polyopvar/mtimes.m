function p = mtimes_old(T,v)
% p = mtimes(T,v) returns the 'polyopvar' object p representing the 'nopvar'
% operator T acting on the 'polyopvar' object v, representing a single state 
% variable of degree 1.
%
% INPUTS
% - T:     nopvar operator with dimension (m,n); 
% - v:     polyopvar object with p=1 and v.degmat = 1;
%
% OUTPUS
% - p:     polyopvar object representing T*v.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - plus
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
% CR, 01/21/2026: Initial coding
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check v is a polyopvar with size(v.degmat) = [1 1].
[nz, p] = size(v.degmat);
if ~isequal([1, 1], [nz, p])
    error("size of v.degmat does not match expected size [1, 1].");
end

% check v.degmat=1.
if ~isequal(1, v.degmat)
    error("v.degmat != 1.");
end

% check v.varmat is all ones - no reason to include unused spatial terms.
Nv = numel(v.pvarname);
if ~isequal(ones(p,Nv), v.varmat)
    error("v includes unused spatial variables.");
end

% check if m=n.
[mv, nv] = size(v.C.ops);
if mv ~= nv
    error("v.C.dim does not have equal dimensions.")
end

% check n is same for T and v.
nt = T.dim(2);
if nt ~= nv
    error("dimensions of T and v do not match.")
end

% check number of spatial variables match.
Nt = size(T.vars,1);
if Nt ~= Nv
    error("number of spatial variables in T and v do not match.")
end

% check name of v and T's spatial variables all match.
if ~isequal(T.vars(:,1), v.pvarname)
    error("names of spatial variables do not match.");
end

% check domains of spatial variables match.
if ~isequal(T.dom, v.dom)
    error("domains of sqatial variables do not match.");
end


% specify polyopvar for T*v.
p = polyopvar();
p.varname = v.varname;
p.degmat = v.degmat;
p.C = tensopvar();
p.C.ops(1) = T;
p.pvarname = v.pvarname;
p.dom = v.dom;
p.varmat = v.varmat;

end


