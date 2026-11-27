function R = dpvar_op_opvar(op,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = DPVAR_OP_OPVAR(OP,...) evaluates the arithmetic OP with a 'dpvar'
% operand (scalar or matrix, constant in space) and an operator X of class
% 'opvar', 'dopvar', 'sopvar' or 'sdopvar':
%
%   'mtimes'  d*X, X*d: a scalar scales every block. A matrix D multiplies
%             the output (D*X) or input (X*D) components, and needs all
%             output (input) spaces of X to be ONE space, which is then the
%             output (input) space of the result;
%   'plus'    d+X, X+d: a scalar is d*I, for X square (the same spaces and
%             dimensions in and out, row i = column i); a matrix of the
%             total size of X is a multiplier split by the dimensions of X,
%             and one that would integrate is refused. A difference is
%             plus(A,-B);
%   'horzcat' a matrix entry maps R^{cols} into the spaces of the rows;
%   'vertcat' a matrix entry maps the spaces of the columns into R^{rows}:
%             for an L2 column that is the integral over the variables of
%             that space.
%
% INPUTS
% - op:     'mtimes' | 'plus' (two operands) | 'horzcat' | 'vertcat' (any);
% - ...:    the operands: at least one an 'opvar', 'dopvar', 'sopvar' or
%           'sdopvar' object, and a 'dpvar' among them (for 'plus', a
%           numeric scalar may take the place of the 'dpvar');
%
% OUTPUTS
% - R:      mtimes/plus on a block: a block ('sdopvar' when a decision
%           variable is involved). On a container, and for any
%           concatenation: a container. A 'dpvar' operand with decision
%           variables makes the result a decision operator.
%
% NOTES
% A 'dpvar' with no decision variables is its constant value, a 'double',
% so dpvar(M)*X, dpvar(c) + X, [dpvar(M), X] give FIXED results.
% Decision x decision is refused (dopvar:decisionTimesDecision): a 'dpvar'
% with a decision variable times a 'dopvar' or 'sdopvar' is quadratic,
% decided by the argument types. A dpvar that depends on a spatial
% variable is refused (mat2opvar_grid:spatial); polynomial operands are
% not handled.
%
% Cost: a scalar product is O(nnz) per block with a one-row B per decision
% variable of d, flat in the q of the container. A matrix product or a sum
% goes through the mtimes/plus of the class against the constant
% container, so it costs what they cost; a sum or concatenation with a
% container on a DIFFERENT decision list pays the merge of the lists
% (O(q log q)). 'plus' puts X first, so that the decision list of X heads
% the union and its blocks keep their rows.
%
% See also MAT2OPVAR_GRID, DPVAR2SDVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - dpvar_op_opvar
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
% Initial coding MMP, 09/29/2026

% Dispatch on the operation
switch op
    case 'mtimes',              R = op_mtimes(varargin{1},varargin{2});
    case 'plus',                R = op_plus(varargin{1},varargin{2});
    case {'horzcat','vertcat'}, R = op_cat(op,varargin);
    otherwise
        error('dpvar_op_opvar:badOp',"Unknown operation '%s'.",op)
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = op_mtimes(A,B)
% R = OP_MTIMES(A,B) is A*B with one operand a 'dpvar' and the other an
% operator: a constant scalar scales through the numeric product of the
% class, a dpvar scalar scales every block (fixed operator only), a matrix
% multiplies the output or input components through its constant
% container.
% INPUTS  A, B: the two factors, one a 'dpvar'.
% OUTPUT  R: block or container; 'sdopvar'/'dopvar' when decision
%         variables enter.
left = ~isoperator(A);              % the matrix is the left factor
if left,    D = A;  X = B;  else,   D = B;  X = A;      end
D = const_value(D);
if isnumeric(D) && isscalar(D)
    % A constant scalar: the numeric product of the class.
    R = mtimes(D,X);
    return
end
if all(size(D)==1)
    % dpvar scalar with decision variables: d*X = X*d, scale each block.
    % decision by type, not by content: a dopvar or sdopvar
    if isa(X,'dopvar') || isa(X,'sdopvar')
        error('dopvar:decisionTimesDecision',['A dpvar with decision '...
            'variables times a %s is quadratic in the decision variables and '...
            'is not representable. Multiply a decision scalar only into a '...
            'fixed operator (opvar/sopvar).'],class(X))
    end
    % d = a + b'*z: a the constant term, b one row per decision variable z;
    % a d that depends on a spatial variable is refused
    [Pc,ZLc] = dpvar2sdvar(D,struct('out',{reshape(D.varname,1,[])},'in',{cell(1,0)}));
    if any(cellfun(@(z) any(z(:)~=0),ZLc))
        error('mat2opvar_grid:spatial',['The dpvar depends on a spatial '...
            'variable; only matrices constant in space are supported here.'])
    end
    a = full(Pc.A);     b = Pc.B;       Zd = reshape(D.dvarname,[],1);
    % A block is scaled; a container has each nonempty block scaled
    if isa(X,'sopvar')
        R = scale_block(X,a,b,Zd);
        return
    end
    Cc = X.C;
    for ii = 1:numel(Cc)
        if ~isempty(Cc{ii}),    Cc{ii} = scale_block(Cc{ii},a,b,Zd);    end
    end
    meta = metadata(X);     meta.Zd = Zd;
    R = dopvar(Cc,meta);
    return
end
% A matrix: the constant container between one side of X and one space,
% then the composition of the class, which refuses decision x decision.
[meta,blk] = side_meta(X);
if left
    % D*X: the output spaces of X must be one space; D maps the output
    % components of X to size(D,1) components on that space
    so = meta.space_out(1,:);
    if any(any(xor(meta.space_out,so)))
        error('mtimes:ambiguousMatrix',['A matrix times an operator whose '...
            'output spaces differ is ambiguous; '...
            'build the multiplier with mat2opvar instead.'])
    end
    if size(D,2)~=sum(meta.dim_out)
        error('mtimes:dimMismatch',['The matrix has %d columns; the operator '...
            'has %d output components.'],size(D,2),sum(meta.dim_out))
    end
    mm = meta;  mm.space_out = so;  mm.dim_out = size(D,1);
    mm.space_in = meta.space_out;   mm.dim_in = meta.dim_out;
    Mop = mat2opvar_grid(D,mm);
    if blk,     R = Mop.C{1,1}*X;   else,   R = Mop*X;      end
else
    % X*D: likewise on the input side
    si = meta.space_in(1,:);
    if any(any(xor(meta.space_in,si)))
        error('mtimes:ambiguousMatrix',['An operator whose input spaces '...
            'differ times a matrix is ambiguous; '...
            'build the multiplier with mat2opvar instead.'])
    end
    if size(D,1)~=sum(meta.dim_in)
        error('mtimes:dimMismatch',['The matrix has %d rows; the operator '...
            'has %d input components.'],size(D,1),sum(meta.dim_in))
    end
    mm = meta;  mm.space_in = si;   mm.dim_in = size(D,2);
    mm.space_out = meta.space_in;   mm.dim_out = meta.dim_in;
    Mop = mat2opvar_grid(D,mm);
    if blk,     R = X*Mop.C{1,1};   else,   R = X*Mop;      end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = op_plus(A,B)
% R = OP_PLUS(A,B) is A + B with one operand a 'dpvar': a scalar is d*I on a
% square operator, a matrix of the total size a multiplier split by the
% dimensions, both built as a constant container and added with X first.
% INPUTS  A, B: the two summands, one a 'dpvar'.
% OUTPUT  R: block or container of the class the sum takes.
if isoperator(A),   X = A;  D = B;  else,   X = B;  D = A;      end
D = const_value(D);
[meta,blk] = side_meta(X);
m = sum(meta.dim_out);  n = sum(meta.dim_in);
square = isequal(meta.space_out,meta.space_in) && isequal(meta.dim_out(:),meta.dim_in(:));
if all(size(D)==1) && square
    % d*I on every diagonal block; off-diagonal slices are zero and stay [].
    Mat = D*speye(m);
elseif isequal(size(D),[m,n])
    Mat = D;
else
    error('plus:dimMismatch',['A %dx%d matrix cannot be added to an operator '...
        'with %d output and %d input components; a scalar needs a square '...
        'operator (same spaces in and out).'],size(D,1),size(D,2),m,n)
end
Mop = mat2opvar_grid(Mat,meta,true);
if blk,     Mop = Mop.C{1,1};   end
% X first: its decision list heads the union, so its blocks keep their rows.
R = X + Mop;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = op_cat(op,args)
% R = OP_CAT(OP,ARGS) concatenates ARGS with 'horzcat' or 'vertcat': each
% matrix entry becomes a container on the spaces of the first operator
% among the operands, its rows (horzcat) or columns (vertcat).
% INPUTS  op: 'horzcat' or 'vertcat'; args: cell of operands.
% OUTPUT  R: the concatenated container.
iref = find(cellfun(@isoperator,args),1);
meta = side_meta(args{iref});
nv = numel(meta.vars);
for k = 1:numel(args)
    a = args{k};
    % Only dpvar entries; anything else goes to the concatenation of the
    % class as it is, which drops [] and refuses a double (BadOperand).
    if ~isa(a,'dpvar'),     continue,   end
    a = const_value(a);
    mm = meta;
    if strcmp(op,'horzcat')
        % horzcat: the entry maps R^{cols} into the row spaces
        if size(a,1)~=sum(meta.dim_out)
            error('horzcat:dimMismatch',['Entry %d has %d rows; the row spaces '...
                'have %d components.'],k,size(a,1),sum(meta.dim_out))
        end
        mm.space_in = false(1,nv);      mm.dim_in = size(a,2);
    else
        % vertcat: the entry maps the column spaces into R^{rows}
        if size(a,2)~=sum(meta.dim_in)
            error('vertcat:dimMismatch',['Entry %d has %d columns; the column '...
                'spaces have %d components.'],k,size(a,2),sum(meta.dim_in))
        end
        mm.space_out = false(1,nv);     mm.dim_out = size(a,1);
    end
    % An empty entry is dropped; any other becomes a constant container
    if any(size(a)==0),     args{k} = [];   continue,   end
    args{k} = mat2opvar_grid(a,mm);
end
% Concatenate with the method of the operator class
if strcmp(op,'horzcat'),    R = horzcat(args{:});   else,   R = vertcat(args{:});   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = isoperator(x)
% TF = ISOPERATOR(X) is true for a container or a block of the four classes.
% INPUT   x: any value.
% OUTPUT  tf: logical, true for an 'opvar', 'dopvar', 'sopvar' or 'sdopvar'.
tf = isa(x,'opvar') || isa(x,'dopvar') || isa(x,'sopvar') || isa(x,'sdopvar');
end


function D = const_value(D)
% D = CONST_VALUE(D): a dpvar with no decision variables is the double it
% equals. The conversion errors if it still depends on a spatial variable,
% which is intended.
% INPUT   D: 'dpvar', or any other value.
% OUTPUT  D: 'double' for a dpvar without decision variables, else D.
if isa(D,'dpvar') && isempty(D.dvarname)
    D = double(D);
end
end


function [meta,blk] = side_meta(X)
% [META,BLK] = SIDE_META(X) returns the container metadata of X; for a
% single block, that of the 1 x 1 container built from the vars and dom of
% the block, O(nv), with no Zd check.
% INPUT   X: 'opvar', 'dopvar', 'sopvar' or 'sdopvar' object.
% OUTPUT  meta: 'struct' with fields vars, dom, space_out, space_in,
%         dim_out, dim_in; blk: logical, true for a block.
blk = isa(X,'sopvar') || isa(X,'sdopvar');
if ~blk
    meta = metadata(X);
    return
end
% Registry: the sorted variables of both sides, with their domains
v = unique([reshape(X.vars.out,1,[]),reshape(X.vars.in,1,[])]);
v = reshape(v,1,[]);        nv = numel(v);
dom = zeros(nv,2);
[tf,loc] = ismember(v,X.vars.out);      dom(tf,:) = X.dom.out(loc(tf),:);
[tf,loc] = ismember(v,X.vars.in);       dom(tf,:) = X.dom.in(loc(tf),:);
meta = struct('vars',{v},'dom',dom,'space_out',ismember(v,X.vars.out),...
    'space_in',ismember(v,X.vars.in),'dim_out',X.dims(1),'dim_in',X.dims(2));
end


function S = scale_block(P,a,b,Zd)
% S = SCALE_BLOCK(P,A,B,ZD) is (a + b'*d)*P for a FIXED block P:
% A = a*vec(C), B = b*vec(C)', per gamma cell. nnz(B) = nnz(b)*nnz(C); b
% has one entry per decision variable of d.
% INPUT   P: 'sopvar' object; a: scalar; b: nd x 1 coefficients of the
%         decision variables; Zd: nd x 1 cellstr of their names.
% OUTPUT  S: 'sdopvar' object on the bases of P.
S = sopvar2sdopvar(P,Zd);
nc = P.dims(1)*prod(cellfun(@numel,P.ZL))*P.dims(2)*prod(cellfun(@numel,P.ZR));
for k = 1:numel(S.params.A)
    v = S.params.A{k};
    if numel(v)~=nc
        % [] or scalar 0, the class zero-cell shorthand: expand, or B would
        % be nd x 1 and A 1 x 1, a block of inconsistent size.
        if nnz(v)
            error('dpvar_op_opvar:badCell','Gamma cell %d of the block has %d entries; expected %d.',k,numel(v),nc)
        end
        v = sparse(nc,1);
    end
    S.params.B{k} = b*v.';
    S.params.A{k} = a*v;
end
end
