function P_clean = clean_ndopvar(P1,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P1 = clean_opvar(P1,tol) removes decvar from ndopvar P1 with coefficients
% with a value less than a tolerance tol
% 
% INPUT
% P1: ndopvar class objects
% tol: acceptable tolerance value, defaults to 1e-14. 
%       -Any term (of any component) of P1 with a coefficient value less 
%        than or equal to tol gets removed.
% 
% OUTPUT
% P_clean: ndopvar class object with only nontrivial terms retained
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - clean_ndopvar
%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding AT - 01/28/2026

if nargin<2
    tol=1e-14;
end    

if ~isa(P1, 'ndopvar')
    error('Currently implemented only for ndopvar')
end


dec_var_OLD = char(P1.dvarname);
N = size(P1.dom,1);
binStr = dec2base(0:(3^N-1), 3);
size_of_monom = P1.dim(1)*prod(P1.deg+1);


indeces_of_decvar = zeros(length(dec_var_OLD)+1, 1); %exclude all decvar
indeces_of_decvar(1) = 1;% always include 1 
indeces_of_decvar = logical(indeces_of_decvar);
% choose new dec var 
for iter_idx = 1:size(binStr, 1)
    substr = binStr(iter_idx, :);
    Qindx = str2num(substr')'+1; % array of .C indeces
    Qindx_cell = num2cell(Qindx); % cell array of .C indeces
    
    full_matrix = max(abs(P1.C{Qindx_cell{:}}), [], 2); 
    full_matrix = reshape(full_matrix, length(dec_var_OLD)+1, []);
    full_matrix =  max(abs(full_matrix), [], 2);
    
    indeces_of_decvar = indeces_of_decvar + (full_matrix > tol   ); 
    indeces_of_decvar = logical(indeces_of_decvar);
    % check if there is a dec var
end

loc_dec_var = find(abs(indeces_of_decvar - 1)<1.e-8); 
dec_var_NEW = dec_var_OLD(loc_dec_var(2:end) - 1, :); % first is just for constant 
% no dec var
cell_newdvar = num2cell(dec_var_NEW, 2); % convert arr of dec var to cell
cell_newdvar = strtrim(cell_newdvar);    % remove whitespace 

P_clean = ndopvar();
P_clean.dvarname = cell_newdvar;
P_clean.dom = P1.dom;
P_clean.deg = P1.deg;
P_clean.vars = P1.vars; 
P_clean.C = cell([3*ones(1,N),1]);
 

indx_array = 1:(size_of_monom*(1+ length(dec_var_OLD)));
indx_array = mod(indx_array, 1 + length(dec_var_OLD));
chosen_idx = ismember(indx_array, mod(loc_dec_var, 1 + length(dec_var_OLD)));
% Loc_L_C = ((0:(size_of_monom - 1))*(1+ length(dec_var_NEW)) + 1)'; % starting location of C rows
% Loc_L_C = [kron(Loc_L_C, ones(length(loc_dec_var), 1)), kron(ones(length(Loc_L_C), 1), loc_dec_var - 1)];
% Loc_L_C = sum(Loc_L_C, 2);


for iter_idx = 1:size(binStr, 1)

    substr = binStr(iter_idx, :);
    Qindx = str2num(substr')'+1; % array of .C indeces
    Qindx_cell = num2cell(Qindx); % cell array of .C indeces

    index_for_monom_in_theta = Qindx == 1; % include monomials for theta

    % change P1 values
    full_matrix = P1.C{Qindx_cell{:}};
    P_clean.C{Qindx_cell{:}} = full_matrix(chosen_idx == 1, :);
end



end