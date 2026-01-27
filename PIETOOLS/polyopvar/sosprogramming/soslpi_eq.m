function prog = soslpi_eq(prog,f)
% PROG = SOSLPI_EQ(PROG,F) takes an LPI program PROG with distributed
% polynomial variable F and adds the constraint F=0 to the program.
%
% INPUTS
% - prog:   lpiprogram structure;
% - f:      'polyopvar' object representing a distributed polynomial
%           variable;
%
% OUTPUTS
% - prog:   lpiprogram structure specifying the same optimization program
%           as the input, but with the added constraint f = 0;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - soslpi_eq
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
% DJ, 01/26/2026: Initial coding


if ~isa(f,'polyopvar')
    error("Distributed polynomial functional must be specified as object of type 'polyopvar'.")
end

for ii=1:size(f.degmat,1)
    % Consider only the term involving monomial ii
    Kfun_ii = f;
    Kfun_ii.degmat = f.degmat(ii,:);
    Kfun_ii.C.ops = f.C.ops(ii);
    % Make sure to enforce uniqueness of the monomials
    Kfun_ii = combine_terms(Kfun_ii);
    % Add the independent variables to the optimization program
    Kop_vars = Kfun_ii.C.ops{1}.vars;
    Kop_varname = Kop_vars.varname;
    is_missing = ~ismember(Kop_varname,prog.vartable.varname);
    prog.vartable = [prog.vartable; polynomial(Kop_varname(is_missing))];
    % Require the coefficient operator acting on the monomial to be zero
    Kop_ii = Kfun_ii.C.ops{1};
    prog = soseq(prog,Kop_ii.params);
end

end