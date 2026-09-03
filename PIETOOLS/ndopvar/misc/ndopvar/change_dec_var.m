function [Cop, dec_var_NEW] = change_dec_var(Aop,dec_var)
% COP = change_dec_var(Aop,dec_var) returns the 'nopvar/ndopvar' object Cop
% with new set of dec_var
%
% INPUTS
% - Aop:     'nopvar' object
% --dec_var:  q x 1 array of decvar names of new decision variables. 
% Aop may have the same names or different.
%
% OUTPUS
% - Cop:      'ndopvar' copy of Aop with new dec var
% - dec_var_NEW:     new names of dec_var -- 
% first elements of dec_var_NEW are dec_var
% Last elements of dec_var_NEW are Aop.dvarname (not in dec_var)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - change_dec_var
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
% AT, 01/21/2026: Initial coding



% Check that the operators can indeed be added
if ~isa(Aop, 'nopvar') && ~isa(Aop, 'ndopvar')
    error('The first argument must be "nopvar" or "ndopvar"')
end

if ~isa(dec_var, 'char')
    try
        % try to convert dec var to char
        % may work if dec_var is a cell array of names as in Aop
        dec_var = char(dec_var);        
    catch
        error('The second argument should be names of dec var "char"')
    end
end
if isa(Aop, 'nopvar')
    dec_var_OLD = [];
else
    % Aop is ndopvar
    dec_var_OLD = char(Aop.dvarname);
end
if isempty(dec_var_OLD)
    loc = [];
    dec_var_NEW = dec_var;
else
    [ind, loc] = ismember(dec_var_OLD, dec_var, 'rows');   
    dec_var_NEW = [dec_var; dec_var_OLD(ind == 0, :)];
end

% [~, loc] = ismember(dec_var_OLD, dec_var_NEW, 'rows');
loc_dec_var = [1; loc + 1]; % shift by 1, since we have [1; xi] dec var
Cop = ndopvar();
Cop.dvarname = num2cell(dec_var_NEW, 2);
Cop.dom = Aop.dom;
Cop.deg = Aop.deg;
Cop.vars = Aop.vars; 
N = size(Aop.dom,1);
Cop.C = cell([3*ones(1,N),1]);
binStr = dec2base(0:(3^N-1), 3);

size_of_monom = Aop.dim(1)*prod(Aop.deg+1); % size of monomimials 

Loc_L_C = ((0:(size_of_monom - 1))*(1+ length(dec_var_NEW)) + 1)'; % starting location of C rows
Loc_L_C = [kron(Loc_L_C, ones(length(loc_dec_var), 1)), kron(ones(length(Loc_L_C), 1), loc_dec_var - 1)];
Loc_L_C = sum(Loc_L_C, 2);

for iter_idx = 1:size(binStr, 1)

    substr = binStr(iter_idx, :);
    Qindx = str2num(substr')'+1; % array of .C indeces
    Qindx_cell = num2cell(Qindx); % cell array of .C indeces

    index_for_monom_in_theta = Qindx == 1; % include monomials for theta

    subdegree1 = Aop.deg;
    subdegree1(index_for_monom_in_theta) = 0;

    size_of_monom_r = Aop.dim(2)*prod(subdegree1 + 1);

    % change Cop values
    full_matrix = sparse(size_of_monom*(length(dec_var_NEW) + 1), size_of_monom_r);
    full_matrix(Loc_L_C, :) = Aop.C{Qindx_cell{:}};
    Cop.C{Qindx_cell{:}} = full_matrix;
end


end