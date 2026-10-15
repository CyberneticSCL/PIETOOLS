function prog = piesos_eq(prog,f,tol)
% PROG = PIESOS_EQ(PROG,F,TOL) adds the equality constraint F=0 to the
% PIESOS program PROG, for a distributed polynomial F.
%
% INPUTS
% - prog:   'struct' specifying a PIESOS program (see `piesos_program`);
% - f:      'polyopvar' object representing a distributed polynomial,
%           with coefficients of type 'dpvar' or 'intop';
% - tol:    (optional) scalar value specifying the tolerance below which
%           coefficients of F are taken to be zero. Defaults to 1e-14;
%
% OUTPUTS
% - prog:   'struct' specifying the same PIESOS program as the input, but
%           with the added constraint F=0;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_eq
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
% DJ, 10/15/2026: Initial coding;


if ~isa(f,'polyopvar')
    error("Distributed polynomial functional must be specified as object of type 'polyopvar'.")
end
if ~f.chkval && ~isvalid(f)
    error("Polynomial is not properly declared.")
end

% Combine coefficients acting on the same monomials, and remove coefficients
% below the tolerance (round-off residues)
if nargin<3 || isempty(tol)
    tol = 1e-14;
end
f = clean_poly(f,tol);

for ii=1:size(f.degmat,1)
    % Consider only the term involving monomial ii
    Kfun_ii = f;
    Kfun_ii.Z = f.Z(ii,:);
    Kfun_ii.C.ops = f.C.ops(ii);
    if isa(Kfun_ii.C.ops{1},'dpvar')
        % Require coefficients to be zero
        prog = soseq(prog,Kfun_ii.C.ops{1});
    elseif isa(Kfun_ii.C.ops{1},'intop')
        % The variables are the dummy variables of the kernels; {:} is applied to
        % a copy, as chained {:} on an empty cell errors in subsref, and
        % [cell(1,0),..] keeps a cellstr when there are no variables
        Kop_dum = Kfun_ii.C.ops{1}.varnames.in_dum;
        Kop_varname = reshape([cell(1,0),Kop_dum{:}],[],1);
        if isa(prog.vartable,'polynomial')
            is_missing = ~ismember(Kop_varname,prog.vartable.varname);
            prog.vartable = [prog.vartable; polynomial(Kop_varname(is_missing))];
        else
            is_missing = ~ismember(Kop_varname,prog.vartable);
            prog.vartable = [prog.vartable; Kop_varname(is_missing)];
        end
        % Require the coefficient operator acting on the monomial to be zero
        Kop_ii = Kfun_ii.C.ops{1};
        % The kernels are a cell of polynomial/dpvar objects; soseq takes one expression at a time
        Kop_kernels = kernels(Kop_ii);
        for k=1:numel(Kop_kernels)
            prog = soseq(prog,Kop_kernels{k});
        end
    else
        error("Coefficients acting on monomials should be specified as 'dpvar' or 'intop' objects.")
    end
end

end