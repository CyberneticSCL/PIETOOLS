function R = dpvar_op_copvar(op,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = DPVAR_OP_COPVAR(OP,...) evaluates the arithmetic OP with a 'dpvar'
% operand (scalar or matrix, constant in space) and an operator X of class
% 'copvar', 'cdopvar', 'sopvar' or 'sdopvar', with the semantics of the
% legacy @opvar/@dopvar methods (plus.m, mtimes.m, horzcat.m, vertcat.m):
%
%   'mtimes'  d*X, X*d: a scalar scales every block. A matrix D multiplies
%             the output (D*X) or input (X*D) components, and needs all
%             output (input) spaces of X to be ONE space, as legacy needs
%             nnz(dim(:,1)) <= 1: that space is then the result's;
%   'plus'    d+X, X+d: a scalar is d*I, for X square (the same spaces and
%             dimensions in and out, row i = column i); a matrix of X's
%             total size is a multiplier split by X's dimensions, and one
%             that would integrate is refused, as 'mat2opvar' refuses Q1.
%             A difference reaches here through minus = plus(A,-B);
%   'horzcat' a matrix entry maps R^{cols} into the rows' spaces;
%   'vertcat' a matrix entry maps the columns' spaces into R^{rows}: for
%             an L2 column that is the integral, legacy's Q1 (see
%             'mat2copvar_grid' for the one meaning a constant has).
%
% INPUTS
% - op:     'mtimes' | 'plus' (two operands) | 'horzcat' | 'vertcat' (any);
% - ...:    the operands, as the class method received them. At least one
%           is a copvar/cdopvar/sopvar/sdopvar; the caller routes here only
%           when a 'dpvar' is among them. Numeric operands keep the classes'
%           own behavior: a numeric scalar scales, anything else numeric is
%           refused as before (use dpvar(M) or mat2copvar_sop for a matrix);
%
% OUTPUTS
% - R:      mtimes/plus on a block: a block ('sdopvar' when a decision
%           variable is involved). On a container, and for any
%           concatenation: a container. Class follows legacy: a 'dpvar'
%           operand makes the result a decision operator.
%
% NOTES
% A 'dpvar' with no decision variables is its constant value, a 'double',
% so dpvar(M)*X, dpvar(c) + X, [dpvar(M), X] give FIXED results.
% Decision x decision stays refused (cdopvar:decisionTimesDecision): a
% 'dpvar' with a decision variable times a 'cdopvar' or 'sdopvar' is
% quadratic. Decided by the argument types, as @cdopvar/mtimes decides it.
% A dpvar that depends on a spatial variable is refused
% (mat2copvar_grid:spatial); polynomial operands are not handled.
%
% Cost: a scalar product is O(nnz) per block with a one-row B per decision
% variable of d, flat in the container's q. A matrix product or a sum goes
% through the class's own mtimes/plus against the constant container, so
% it costs what they cost; a sum or concatenation with a container on a
% DIFFERENT decision list pays 'merge_dvar_lists' (O(q log q)). 'plus' puts
% X first so that X's own list is the union's head and its blocks keep
% their rows.
%
% See also MAT2COPVAR_GRID, DPVAR2SDVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - dpvar_op_copvar
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
% Initial coding MMP, 09/29/2026. Tier 1b of the container parity map: gamma
%                as a decision variable (gam*I, gam - T'*T, [-gam, D'; ..])
%                for the H-inf/H2 executives and DEMO2/3/5/6/9, one branch
%                shared by the four classes' mtimes/plus/horzcat/vertcat.

switch op
    case 'mtimes',              R = op_mtimes(varargin{1},varargin{2});
    case 'plus',                R = op_plus(varargin{1},varargin{2});
    case {'horzcat','vertcat'}, R = op_cat(op,varargin);
    otherwise
        error('dpvar_op_copvar:badOp',"Unknown operation '%s'.",op)
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = op_mtimes(A,B)
left = ~isoperator(A);              % the matrix is the left factor
if left,    D = A;  X = B;  else,   D = B;  X = A;      end
D = const_value(D);
if isnumeric(D) && isscalar(D)
    % A constant scalar is the classes' own numeric path.
    R = mtimes(D,X);
    return
end
if all(size(D)==1)
    % dpvar scalar with decision variables: d*X = X*d, scale each block.
    if isdecision(X)
        error('cdopvar:decisionTimesDecision',['A dpvar with decision '...
            'variables times a %s is quadratic in the decision variables and '...
            'is not representable. Multiply a decision scalar only into a '...
            'fixed operator (copvar/sopvar).'],class(X))
    end
    [Pc,ZLc] = dpvar2sdvar(D,struct('out',{reshape(D.varname,1,[])},'in',{cell(1,0)}));
    if any(cellfun(@(z) any(z(:)~=0),ZLc))
        error('mat2copvar_grid:spatial',['The dpvar depends on a spatial '...
            'variable; only matrices constant in space are supported here.'])
    end
    a = full(Pc.A);     b = Pc.B;       Zd = reshape(D.dvarname,[],1);
    if isa(X,'sopvar')
        R = scale_block(X,a,b,Zd);
        return
    end
    Cc = X.C;
    for ii = 1:numel(Cc)
        if ~isempty(Cc{ii}),    Cc{ii} = scale_block(Cc{ii},a,b,Zd);    end
    end
    meta = metadata(X);     meta.Zd = Zd;
    R = cdopvar(Cc,meta);
    return
end
% A matrix: the constant container between X's side and one space, then
% the class's own composition, which refuses decision x decision.
[meta,blk] = side_meta(X);
if left
    so = meta.space_out(1,:);
    if any(any(xor(meta.space_out,so)))
        error('mtimes:ambiguousMatrix',['A matrix times an operator whose '...
            'output spaces differ is ambiguous (legacy: nnz(dim(:,1)) > 1); '...
            'build the multiplier with mat2copvar_sop instead.'])
    end
    if size(D,2)~=sum(meta.dim_out)
        error('mtimes:dimMismatch',['The matrix has %d columns; the operator '...
            'has %d output components.'],size(D,2),sum(meta.dim_out))
    end
    mm = meta;  mm.space_out = so;  mm.dim_out = size(D,1);
    mm.space_in = meta.space_out;   mm.dim_in = meta.dim_out;
    Mop = mat2copvar_grid(D,mm);
    if blk,     R = Mop.C{1,1}*X;   else,   R = Mop*X;      end
else
    si = meta.space_in(1,:);
    if any(any(xor(meta.space_in,si)))
        error('mtimes:ambiguousMatrix',['An operator whose input spaces '...
            'differ times a matrix is ambiguous (legacy: nnz(dim(:,2)) > 1); '...
            'build the multiplier with mat2copvar_sop instead.'])
    end
    if size(D,1)~=sum(meta.dim_in)
        error('mtimes:dimMismatch',['The matrix has %d rows; the operator '...
            'has %d input components.'],size(D,1),sum(meta.dim_in))
    end
    mm = meta;  mm.space_in = si;   mm.dim_in = size(D,2);
    mm.space_out = meta.space_in;   mm.dim_out = meta.dim_in;
    Mop = mat2copvar_grid(D,mm);
    if blk,     R = X*Mop.C{1,1};   else,   R = X*Mop;      end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = op_plus(A,B)
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
Mop = mat2copvar_grid(Mat,meta,true);
if blk,     Mop = Mop.C{1,1};   end
% X first: its decision list heads the union, so its blocks keep their rows.
R = X + Mop;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = op_cat(op,args)
% Each matrix entry becomes a container on the spaces of the first operator
% among the operands: its rows (horzcat) or columns (vertcat).
iref = find(cellfun(@isoperator,args),1);
meta = side_meta(args{iref});
nv = numel(meta.vars);
for k = 1:numel(args)
    a = args{k};
    % Only dpvar entries; anything else goes to the class's concatenation
    % as it is, which drops [] and refuses a double (BadOperand).
    if ~isa(a,'dpvar'),     continue,   end
    a = const_value(a);
    mm = meta;
    if strcmp(op,'horzcat')
        if size(a,1)~=sum(meta.dim_out)
            error('horzcat:dimMismatch',['Entry %d has %d rows; the row spaces '...
                'have %d components.'],k,size(a,1),sum(meta.dim_out))
        end
        mm.space_in = false(1,nv);      mm.dim_in = size(a,2);
    else
        if size(a,2)~=sum(meta.dim_in)
            error('vertcat:dimMismatch',['Entry %d has %d columns; the column '...
                'spaces have %d components.'],k,size(a,2),sum(meta.dim_in))
        end
        mm.space_out = false(1,nv);     mm.dim_out = size(a,1);
    end
    if any(size(a)==0),     args{k} = [];   continue,   end
    args{k} = mat2copvar_grid(a,mm);
end
if strcmp(op,'horzcat'),    R = horzcat(args{:});   else,   R = vertcat(args{:});   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = isoperator(x)
tf = isa(x,'copvar') || isa(x,'cdopvar') || isa(x,'sopvar') || isa(x,'sdopvar');
end


function tf = isdecision(X)
% By type, not by content: a cdopvar or sdopvar is decision, as in
% @cdopvar/mtimes. 'sdopvar' is a subclass of nothing here, so isa is exact.
tf = isa(X,'cdopvar') || isa(X,'sdopvar');
end


function D = const_value(D)
% A dpvar with no decision variables is the double it equals. 'double'
% errors if it still depends on a spatial variable, which is intended.
if isa(D,'dpvar') && isempty(D.dvarname)
    D = double(D);
end
end


function [meta,blk] = side_meta(X)
% Container metadata of X; for a single block, the 1 x 1 container's,
% built from the block's own vars and dom, O(nv), with no Zd check.
blk = isa(X,'sopvar') || isa(X,'sdopvar');
if ~blk
    meta = metadata(X);
    return
end
v = unique([reshape(X.vars.out,1,[]),reshape(X.vars.in,1,[])]);
v = reshape(v,1,[]);        nv = numel(v);
dom = zeros(nv,2);
[tf,loc] = ismember(v,X.vars.out);      dom(tf,:) = X.dom.out(loc(tf),:);
[tf,loc] = ismember(v,X.vars.in);       dom(tf,:) = X.dom.in(loc(tf),:);
meta = struct('vars',{v},'dom',dom,'space_out',ismember(v,X.vars.out),...
    'space_in',ismember(v,X.vars.in),'dim_out',X.dims(1),'dim_in',X.dims(2));
end


function S = scale_block(P,a,b,Zd)
% (a + b'*d)*P for a FIXED block P: A = a*vec(C), B = b*vec(C)', per gamma
% cell. nnz(B) = nnz(b)*nnz(C); b has one entry per decision variable of d.
S = sopvar2sdopvar(P,Zd);
nc = P.dims(1)*prod(cellfun(@numel,P.ZL))*P.dims(2)*prod(cellfun(@numel,P.ZR));
for k = 1:numel(S.params.A)
    v = S.params.A{k};
    if numel(v)~=nc
        % [] or scalar 0, the class zero-cell shorthand: expand, or B would
        % be nd x 1 and A 1 x 1 (getsol and mtimes then refuse the block).
        if nnz(v)
            error('dpvar_op_copvar:badCell','Gamma cell %d of the block has %d entries; expected %d.',k,numel(v),nc)
        end
        v = sparse(nc,1);
    end
    S.params.B{k} = b*v.';
    S.params.A{k} = a*v;
end
end
