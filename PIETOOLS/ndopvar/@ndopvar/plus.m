function Cop = plus(Aop,Bop)
% COP = PLUS(AOP,BOP) returns the 'ndopvar' object COP representing the sum
% of the PI operators defined by 'ndopvar' objects AOP and BOP
%
% INPUTS
% - Aop:    m x n 'ndopvar' object
% - Bop:    m x n 'ndopvar' object
%
% OUTPUS
% - Cop:    m x n 'ndopvar' object representing the sum of the operators
%           defined by Aop and Bop
%
% NOTES
% The operators defined by Aop and Bop must act on functions on the same
% spatial domain, i.e. Aop.dom = Bop.dom.

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
% DJ, 01/15/2026: Initial coding
% AT, 01/21/2026: allow operators with different degree

% Check that the operators can indeed be added
if any(Aop.dim~=Bop.dim)
    error("Dimensions of the operators must match.")
end
if size(Aop.dom,1)~=size(Bop.dom,1) || any(any(Aop.dom~=Bop.dom))
    error("Spatial domains on which operators are defined should match.")
end

% convert to the same degree
if any(Aop.deg~=Bop.deg)
    max_degree = max(Aop.deg, Bop.deg);
    Aop_new = change_degree(Aop, max_degree);
    Bop_new = change_degree(Bop, max_degree);
    Aop = Aop_new;
    Bop = Bop_new;
    % error("Addition of operators with different monomial degrees is currently not supported.")
end

if isa(Aop, 'ndopvar')
    Aop_dvarname = Aop.dvarname;
else
    Aop_dvarname = {};
end
if isa(Bop, 'ndopvar')
    Bop_dvarname = Bop.dvarname;
else
    Bop_dvarname = {};
end
if numel(Aop_dvarname) ~= numel(Bop_dvarname) || ~isequal(Aop_dvarname,Bop_dvarname)
    dvars1 = string(Aop_dvarname); % convert array to char array
    dvars2 = string(Bop_dvarname);
    % numberOfCharacters = max(size(dvars1{1}, 2), size(dvars2{1}, 2));
    % dvars1 = pad(dvars1, numberOfCharacters); % pad with ' ' if needed
    % dvars2 = pad(dvars2, numberOfCharacters); % pad with ' ' if needed 

    if isempty(dvars1) || isempty(dvars2)
        common_dvar = [];
        full_dvars = [dvars1; dvars2];
    else
        common_dvar = intersect(dvars1, dvars2, 'rows');
        new_dvars   = setdiff(dvars2, common_dvar, 'rows');
        full_dvars = char([dvars1; new_dvars]);
    end
    Aop = change_dec_var(Aop, full_dvars);
    Bop = change_dec_var(Bop, full_dvars);
end
    

% Assuming the same monomial degrees and decision variables, the sum of the 
% operators is just defined by the sum of the coefficients.
Cop = Aop;
for ii=1:numel(Cop.C)
    Cop.C{ii} = Aop.C{ii} + Bop.C{ii};
end

end