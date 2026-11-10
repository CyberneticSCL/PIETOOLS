function disp_struct = get_display(Kop,xname)
% DISP_STRUCT = GET_DISPLAY(KOP,XNAME) constructs the terms of the
% expression
%   int ... int K(s) x(s) ds ... ds + ...
% of the functional KOP acting on a function x, for use in the 'display'
% method of 'intop', and in the 'get_display' methods of the classes of
% operators that have 'intop' objects as blocks.
%
% INPUTS
% - Kop:    'intop' object representing a functional, see 'intop';
% - xname:  (optional) char, the name of the function that the functional
%           acts on. Defaults to x;
%
% OUTPUTS
% - disp_struct:      struct with the fields
%       terms:  1 x n cell of char, the terms of the sum, without their
%               sign, see 'display_term';
%       signs:  1 x n array of 1 and -1, the sign of each term;
%       use_K:  logical, true if a kernel is displayed by name. The kernel
%               K_i can then be extracted by calling Kop.kernels{i};
%       ok:     logical, true if the functional is displayed explicitly (no
%               kernel by name), with at most 6 terms of at most 100
%               characters, so that it may be part of a larger display;
%       xname:  char, the name of the function that the functional acts on;
%       is_id:  logical, true if Kop is an identity operator, in which case
%               the only term is the name of the function, x;
%
% NOTES
% - The integrals involve only dummy variables, displayed as gamma, with
%   the number of the integrated copy as subscript if there are multiple,
%   see 'display_term'. The kernels Kop.kernels that can be extracted for a
%   kernel displayed by name use dummy variables named as the variables in
%   Kop.vars with _dum appended.
% - A kernel with at most 4 nonzero coefficients is displayed explicitly,
%   if Kop is scalar-valued. Kernels of a matrix-valued functional are
%   always displayed by name. The terms are numbered like the elements of
%   Kop.params, so a term with a zero kernel is not displayed, but does
%   take an index.
% - An identity operator is displayed as x, see 'is_id_op'.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - get_display
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
% DJ, 11/06/2026: Initial coding;

% Default settings
max_terms = 6;          % maximal number of terms of an operator displayed explicitly
max_chars = 100;        % maximal number of characters of each such term

N = numel(Kop.vars.in);
params_size = [size(Kop.params,1:N),1];       % trailing 1 for ind2sub, if N=1

if nargin<2 || isempty(xname)
    xname = 'x';
end
disp_struct = struct('terms',{cell(1,0)},'signs',zeros(1,0),'use_K',false,'ok',false,'xname',xname,'is_id',false);

if is_id_op(Kop)
    disp_struct.terms = {xname};      disp_struct.signs = 1;
    disp_struct.ok = true;            disp_struct.is_id = true;
    return
end

for k=1:numel(Kop.params)
    if nnz(Kop.params{k})==0
        continue
    end
    sub_idcs = cell(1,N);
    if N>0
        [sub_idcs{:}] = ind2sub(params_size,k);
    end
    [term_str,term_sign,use_K] = display_term(Kop,cell2mat(sub_idcs),k,xname);
    disp_struct.terms{end+1} = term_str;      disp_struct.signs(end+1) = term_sign;
    disp_struct.use_K = disp_struct.use_K || use_K;
end
disp_struct.ok = ~disp_struct.use_K && numel(disp_struct.terms)<=max_terms && all(cellfun(@numel,disp_struct.terms)<=max_chars);

end
