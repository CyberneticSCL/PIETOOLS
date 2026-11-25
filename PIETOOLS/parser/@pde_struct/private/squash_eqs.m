function PDE_new = squash_eqs(PDE,sz_eqs)
% PDE_NEW = SQUASH_EQS(PDE,SZ_EQS) takes a PDE structure describing free
% terms (no completed equations), and squashes all rows of terms (each row
% to be added to a separate equation) into a single vector-valued row of
% terms.
%
% INPUT
% - PDE:        'pde_struct' object representing a set of n rows of free
%               terms to be used to define equations, specified through the
%               fields PDE.free{ii}.term, for ii=1:n;
% - sz_eqs:     nx1 array specifying for each row of terms the vector size
%               of those terms, sz_eqs(ii) = PDE.free{ii}.size;
%
% OUTPUT
% - PDE_new:    'pde_struct' object representing the same equations, but
%               now with PDE_new.free only containing a single element,
%               with PDE_new.free{1}.term including all terms from all
%               elements PDE.free{ii}, premultiplied with a suitable
%               permutation matrix;
%
% NOTES
% A nonlinear term, the elementwise product of factors C_k*x_k, is
% premultiplied as a 'pde_struct' product, which splits it into one term
% per column of the permutation matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% DJ, 01/06/2025: Initial coding;
% DJ, 11/24/2026: Move from eq.m to be shared with plus; also squash
%                   nonlinear terms;

% Initialize a new object with just 1 vector-valued equation
PDE_new = PDE;
PDE_new.free = PDE.free(1);
PDE_new.free{1}.term = {};
var_list = cell(0,1);
for ii=1:numel(sz_eqs)
    if ~isfield(PDE.free{ii},'term') || isempty(PDE.free{ii}.term)
        % A row of zeros, as in [x(2);0], adds no terms.
        continue
    end
    % Augment all terms in row ii to new size by premultiplying coefficients.
    C_aug = [zeros(sum(sz_eqs(1:ii-1)),sz_eqs(ii));
             eye(sz_eqs(ii));
             zeros(sum(sz_eqs(ii+1:end)),sz_eqs(ii))];
    for jj=1:numel(PDE.free{ii}.term)
        trm = PDE.free{ii}.term{jj};
        if ~isscalar(trm)
            % C_aug times a nonlinear term is a sum of nonlinear terms.
            row_jj = PDE;   row_jj.free = PDE.free(ii);     row_jj.free{1}.term = {trm};
            row_jj = C_aug*row_jj;
            PDE_new.free{1}.term = [PDE_new.free{1}.term,row_jj.free{1}.term];
        elseif isfield(trm,'C')
            trm.C = C_aug*trm.C;
            PDE_new.free{1}.term{end+1} = trm;
        else
            trm.C = C_aug;
            PDE_new.free{1}.term{end+1} = trm;
        end
    end
    % Keep track of the variables that appear in all rows.
    for kk=1:size(PDE.free{ii}.vars,1)
        varname_kk = PDE.free{ii}.vars(kk).varname{1};
        if ~ismember(varname_kk,var_list)
            var_list = [var_list;{varname_kk}];                             %#ok<AGROW>
        end
    end
end
% Set the variables and size of the new squashed row of terms
PDE_new.free{1}.size = sum(sz_eqs);
PDE_new.free{1}.vars = polynomial(var_list);

end
