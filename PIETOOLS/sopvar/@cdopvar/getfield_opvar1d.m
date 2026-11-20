function val = getfield_opvar1d(P,fname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% val = getfield_opvar1d(P,fname) returns the value that the field fname
% would have if the PI operator
%   P: R^m x L2^n[a,b] -> R^p x L2^q[a,b]
% were a 'dopvar' object, acting as
%   [ P    Q1 ] [ u ]       (P u + int_a^b Q1(s) x(s) ds )
%   [ Q2   R  ] [ x ] (s) = (Q2(s) u + R0(s) x(s) + int_a^s R1(s,t) x(t) dt
%                                    + int_s^b R2(s,t) x(t) dt ).
%
% INPUT
% P:      cdopvar class object with the structure of a 'dopvar': at most
%         one spatial variable, and every block row (column) of R before
%         every block row (column) of L2
% fname:  char or string, the name of a 'dopvar' field: 'P', 'Q1', 'Q2',
%         'R', 'R.R0', 'R.R1', 'R.R2', 'I', 'var1', 'var2' or 'dim'
%
% OUTPUT
% val:    the value of the field:
%         'P':      p x m 'dpvar', constant in var1;
%         'Q1':     p x n 'dpvar' in var1;
%         'Q2':     q x m 'dpvar' in var1;
%         'R':      struct with fields R0, R1, R2, each a q x n 'dpvar',
%                   R0 in var1, R1 and R2 in var1 and var2;
%         'R.R0', 'R.R1', 'R.R2': that field of the struct 'R';
%         'I':      1 x 2 array [a,b], the domain P.dom of the variable;
%         'var1':   1 x 1 'polynomial', the spatial variable P.vars{1};
%         'var2':   1 x 1 'polynomial', the dummy variable, named as var1
%                   with '_dum' appended;
%         'dim':    2 x 2 array [p, m; q, n];
%         The 'dpvar' components are affine in the decision variables of
%         P.Zd that appear in them;
%
% NOTES:
% - The block rows of R of P are stacked in order into the rows of P and
%   Q1, those of L2 into the rows of Q2 and R; likewise for the block
%   columns. The structure of P ensures this keeps the order of the
%   components.
% - A component is assembled from the kernels of the blocks, so it is in
%   the output variable var1 and the dummy variable var2. For R, the
%   parameters 1, 2, 3 of a block are the multiplier and the integrals over
%   [a,s] and [s,b], giving R0, R1, R2. The kernel of Q1 is in the dummy variable
%   alone, which is renamed to var1, as a 'dopvar' defines Q1 in var1.
% - A block may be 'sdopvar', 'sopvar' (no decision variables, kernels
%   'polynomial') or []. Each component is returned as 'dpvar', as the
%   components of a 'dopvar' are, also if no block of it has decision
%   variables.
% - A structurally zero block ([]) gives zeros. An absent space gives a
%   component with zero rows or columns, as 'dopvar' does for a zero
%   dimension.
% - If P has no spatial variable, P: R^m -> R^p, then 'I' is 0 x 2, and
%   'var1', 'var2' are s and s_dum, the default variables of a new 'dopvar'
%   (P does not depend on them), so that s = P.var1 works on a new cdopvar.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - getfield_opvar1d (cdopvar)
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% DJ, 11/15/2026: Initial coding;

% Check the structure of P and the field name
if ~is_opvar1d(P)
    error('cdopvar:notOpvar1d',['The operator does not have the structure of a ''dopvar'', '...
        'R^m x L2^n[a,b] -> R^p x L2^q[a,b]; see ''is_opvar1d''.'])
end
if ~(ischar(fname) || (isstring(fname) && isscalar(fname)))
    error('cdopvar:badField','The field name must be a char or string.')
end
fname = char(fname);

% Block rows (columns) of L2; the others are R
isL2_out = any(P.space_out,2);          isL2_in = any(P.space_in,2);
% Names of the variable and its dummy, as 'sdopvar' names them in kernels
vname = reshape(P.vars,[],1);           dname = strcat(vname,'_dum');

% Extract the desired field. Use the assemble_kernels helper to extract
% the kernels associated with a particular field.
switch fname
    case 'I'
        val = P.dom;
    case 'var1'
        % Without a variable, the default s of a new 'dopvar'
        val = polynomial(vname);
        if isempty(vname),  val = polynomial({'s'});    end
    case 'var2'
        val = polynomial(dname);
        if isempty(dname),  val = polynomial({'s_dum'});    end
    case 'dim'
        val = [sum(P.dim_out(~isL2_out)), sum(P.dim_in(~isL2_in));
               sum(P.dim_out(isL2_out)),  sum(P.dim_in(isL2_in))];
    case 'P'
        K = assemble_kernels(P,~isL2_out,~isL2_in,1);
        val = K{1};
    case 'Q1'
        % The kernel is in the dummy variable only; rename it to var1
        K = assemble_kernels(P,~isL2_out,isL2_in,1);
        val = K{1};
        if ~isempty(val.varname)
            vnames = val.varname;
            vnames(strcmp(vnames,dname)) = vname;
            val = dpvar(val.C,val.degmat,vnames,val.dvarname,val.matdim);
        end
    case 'Q2'
        % The multiplier from R into L2, in var1
        K = assemble_kernels(P,isL2_out,~isL2_in,1);
        val = K{1};
    case {'R','R.R0','R.R1','R.R2'}
        % The kernels R0 (multiplier), R1 (int_a^s) and R2 (int_s^b) of the
        % L2 blocks; one of them for 'R.R0', 'R.R1' or 'R.R2'
        K = assemble_kernels(P,isL2_out,isL2_in,3);
        val = struct('R0',K{1},'R1',K{2},'R2',K{3});
        if numel(fname)>1
            val = val.(fname(3:end));
        end
    otherwise
        error('cdopvar:badField',['''%s'' is not a field of a ''dopvar''; use one of '...
            '''P'', ''Q1'', ''Q2'', ''R'', ''R.R0'', ''R.R1'', ''R.R2'', ''I'', '...
            '''var1'', ''var2'', ''dim''.'],fname)
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function K = assemble_kernels(P,rows,cols,nk)
% K = assemble_kernels(P,rows,cols,nk) stacks the kernels of the blocks
% P.C(rows,cols) into one 'dpvar' per parameter.
%
% INPUT
% P:    cdopvar class object, opvar-like
% rows: M x 1 logical, the block rows to stack
% cols: N x 1 logical, the block columns to stack
% nk:   number of parameters of each block: 1 for P, Q1, Q2, and 3 for R
%
% OUTPUT
% K:    1 x nk cell of sum(P.dim_out(rows)) x sum(P.dim_in(cols)) 'dpvar',
%       K{k} the kernel of parameter k of every block, with zeros for a
%       structurally zero block
%
% NOTES
% The kernels are computed using the '@sdopvar/kernels' method ('dpvar'),
% or '@sopvar/kernels' for a block without decision variables
% ('polynomial'), converted to 'dpvar' with 'poly2dpvar', so that every
% kernel is 'dpvar'. The 'dpvar' concatenation merges the decision
% variables of the blocks.

ir = find(rows);    jc = find(cols);
K = cell(1,nk);
% No block: zeros of the size of the absent (or empty) space
if isempty(ir) || isempty(jc)
    K(:) = {dpvar(zeros(sum(P.dim_out(rows)),sum(P.dim_in(cols))))};
    return
end
Kblk = cell(numel(ir),numel(jc),nk);
for a=1:numel(ir)
    for b=1:numel(jc)
        blk = P.C{ir(a),jc(b)};
        if isempty(blk)
            Kblk(a,b,:) = {dpvar(zeros(P.dim_out(ir(a)),P.dim_in(jc(b))))};
        else
            Kb = blk.kernels;
            if isa(blk,'sopvar')
                % 'polynomial' kernels; {}: no variable is a decision
                % variable (poly2dpvar(K) would take any 'coeff*' name as one)
                Kb = cellfun(@(Kk) poly2dpvar(Kk,{}),Kb,'UniformOutput',false);
            end
            Kblk(a,b,:) = Kb(1:nk);
        end
    end
end
for k=1:nk
    Krows = cell(numel(ir),1);
    for a=1:numel(ir)
        Krows{a} = horzcat(Kblk{a,:,k});
    end
    K{k} = vertcat(Krows{:});
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
