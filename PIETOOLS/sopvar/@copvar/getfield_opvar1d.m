function val = getfield_opvar1d(P,fname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% val = getfield_opvar1d(P,fname) returns the value that the field fname
% would have if the PI operator
%   P: R^m x L2^n[a,b] -> R^p x L2^q[a,b]
% were an 'opvar' object, acting as
%   [ P    Q1 ] [ u ]       (P u + int_a^b Q1(s) x(s) ds )
%   [ Q2   R  ] [ x ] (s) = (Q2(s) u + R0(s) x(s) + int_a^s R1(s,t) x(t) dt
%                                    + int_s^b R2(s,t) x(t) dt ).
%
% INPUT
% P:      copvar class object with the structure of an 'opvar', see
%         'is_opvar1d'
% fname:  char or string, the name of an 'opvar' field: 'P', 'Q1', 'Q2',
%         'R', 'R.R0', 'R.R1', 'R.R2', 'I', 'var1', 'var2' or 'dim'
%
% OUTPUT
% val:    the value of the field:
%         'P':      p x m 'polynomial', constant;
%         'Q1':     p x n 'polynomial' in var1;
%         'Q2':     q x m 'polynomial' in var1;
%         'R':      struct with fields R0, R1, R2, each a q x n 'polynomial',
%                   R0 in var1, R1 and R2 in var1 and var2;
%         'R.R0', 'R.R1', 'R.R2': that field of the struct 'R';
%         'I':      1 x 2 array [a,b], the domain P.dom of the variable;
%         'var1':   1 x 1 'polynomial', the spatial variable P.vars{1};
%         'var2':   1 x 1 'polynomial', the dummy variable, named as var1
%                   with '_dum' appended;
%         'dim':    2 x 2 array [p, m; q, n];
%
% NOTES:
% - The block rows of R of P are stacked in order into the rows of P and
%   Q1, those of L2 into the rows of Q2 and R; likewise for the block
%   columns. 'is_opvar1d' ensures this keeps the order of the components.
% - A component is assembled from the kernels of the blocks, 'kernels' of
%   'sopvar', so it is in the output variable var1 and the dummy variable
%   var2. For R, the parameters 1, 2, 3 of a block are the multiplier and
%   the integrals over [a,s] and [s,b] ('decode_params' of 'sopvar'), 
%   giving R0, R1, R2. The kernel of Q1 is in the dummy variable alone, 
%   which is renamed to var1, as an 'opvar' defines Q1 in var1.
% - A structurally zero block ([]) gives zeros. An absent space gives a
%   component with zero rows or columns, as 'opvar' does for a zero
%   dimension.
% - If P has no spatial variable, P: R^m -> R^p, then 'I' is 0 x 2, and
%   'var1', 'var2' are s and s_dum, the default variables of a new 'opvar'
%   (P does not depend on them), so that s = P.var1 works on a new copvar.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - getfield_opvar1d (copvar)
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
% DJ, 11/15/2026: Initial coding;

if ~is_opvar1d(P)
    error('copvar:notOpvar1d',['The operator does not have the structure of an ''opvar'', '...
        'R^m x L2^n[a,b] -> R^p x L2^q[a,b]; see ''is_opvar1d''.'])
end
if ~(ischar(fname) || (isstring(fname) && isscalar(fname)))
    error('copvar:badField','The field name must be a char or string.')
end
fname = char(fname);

% Block rows (columns) of L2; the others are R
isL2_out = any(P.space_out,2);          isL2_in = any(P.space_in,2);
% Names of the variable and its dummy, as 'sopvar' names them in kernels
vname = reshape(P.vars,[],1);           dname = strcat(vname,'_dum');

% Extract the desired field. Use the assemble_kernels helper to extract
% the kernels associated with a particular field.
switch fname
    case 'I'
        val = P.dom;
    case 'var1'
        % Without a variable, the default s of a new 'opvar'
        val = polynomial(vname);
        if isempty(vname),  val = polynomial({'s'});    end                 % DJ, 11/15/2026
    case 'var2'
        val = polynomial(dname);
        if isempty(dname),  val = polynomial({'s_dum'});    end             % DJ, 11/15/2026
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
            val = polynomial(val.coefficient,val.degmat,vnames,val.matdim);
        end
    case 'Q2'
        K = assemble_kernels(P,isL2_out,~isL2_in,1);
        val = K{1};
    case {'R','R.R0','R.R1','R.R2'}
        K = assemble_kernels(P,isL2_out,isL2_in,3);
        val = struct('R0',K{1},'R1',K{2},'R2',K{3});
        if numel(fname)>1
            val = val.(fname(3:end));
        end
    otherwise
        error('copvar:badField',['''%s'' is not a field of an ''opvar''; use one of '...
            '''P'', ''Q1'', ''Q2'', ''R'', ''R.R0'', ''R.R1'', ''R.R2'', ''I'', '...
            '''var1'', ''var2'', ''dim''.'],fname)
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function K = assemble_kernels(P,rows,cols,nk)
% K = assemble_kernels(P,rows,cols,nk) stacks the kernels of the blocks
% P.C(rows,cols) into one 'polynomial' per parameter.
%
% INPUT
% P:    copvar class object, opvar-like
% rows: M x 1 logical, the block rows to stack
% cols: N x 1 logical, the block columns to stack
% nk:   number of parameters of each block: 1 for P, Q1, Q2, and 3 for R
%
% OUTPUT
% K:    1 x nk cell of sum(P.dim_out(rows)) x sum(P.dim_in(cols))
%       'polynomial', K{k} the kernel of parameter k of every block, with
%       zeros for a structurally zero block
%
% NOTES
% The polynomial kernels are computed using the '@sopvar/kernels' method.

ir = find(rows);    jc = find(cols);
K = cell(1,nk);
% No block: zeros of the size of the absent (or empty) space
if isempty(ir) || isempty(jc)
    K(:) = {polynomial(zeros(sum(P.dim_out(rows)),sum(P.dim_in(cols))))};
    return
end
Kblk = cell(numel(ir),numel(jc),nk);
for a=1:numel(ir)
    for b=1:numel(jc)
        blk = P.C{ir(a),jc(b)};
        if isempty(blk)
            Kblk(a,b,:) = {polynomial(zeros(P.dim_out(ir(a)),P.dim_in(jc(b))))};
        else
            Kb = blk.kernels;
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
