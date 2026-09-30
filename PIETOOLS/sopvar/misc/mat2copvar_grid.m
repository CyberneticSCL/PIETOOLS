function P = mat2copvar_grid(Mat,meta,mult_only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = MAT2COPVAR_GRID(MAT,META,MULT_ONLY) returns the container of the
% spatially constant matrix MAT acting between the spaces META describes,
%
%   (P*x)_i = sum_j K_ij x_j,   K_ij the (i,j) block of MAT, cut by
%                               META.dim_out (rows) and META.dim_in (cols).
%
% The kernel form (sopvar document Sec. 4) gives a constant kernel exactly
% one meaning per pair of spaces: a multiplier in every variable the two
% spaces share (gamma = 1), constant in an output-only variable, and
% integrated over an input-only variable. In 1-D that is 'mat2opvar' (P, Q2,
% R0) plus the Q1 integral that 'dopvar/vertcat' makes of a matrix stacked
% on an L2 operator.
%
% INPUTS
% - Mat:       sum(dim_out) x sum(dim_in) 'double' or 'dpvar', constant in
%              the spatial variables;
% - meta:      container metadata, fields vars (sorted registry), dom,
%              space_out, space_in, dim_out, dim_in, as 'metadata' returns;
% - mult_only: (optional, default false) true refuses a nonzero block that
%              would integrate, i.e. whose input space has a variable the
%              output space lacks, as 'mat2opvar' refuses Q1;
%
% OUTPUTS
% - P:         'copvar' for a double, 'cdopvar' for a 'dpvar' (Zd = its
%              decision variables), by input CLASS, as 'mat2opvar' returns
%              opvar or dopvar.
%
% NOTES
% A zero block of MAT stays [], except one explicit zero block per row or
% column that would otherwise have none, placed as 'lpivar_cdopvar' places
% it: 'verify' and cdopvar(C) read a row's space off a block.
%
% Canonical by construction: every basis is degree 0, so the multiplier
% cell has content only at right degree 0 (CANONICAL MULTIPLIER FORM,
% sopvar.m), and the integral cells are zero.
%
% Cost: O(M*N) blocks. A block is the slice A(lin), B(:,lin) of the
% coefficients, O(nnz of the slice + number of columns), plus the column
% pointers of 3^n3 cells of q_i*p_j columns; nothing is O(q) but the one
% 'dpvar2sdvar' call, which is O(nnz(Mat.C)).
%
% See also DPVAR_OP_COPVAR, DPVAR2SDVAR, MAT2OPVAR, LPIVAR_CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mat2copvar_grid
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
% Initial coding MMP, 09/29/2026. The shared core of the Tier-1 dpvar
%                operator branches (dpvar_op_copvar) and the constructors
%                eye_copvar_sop, zeros_copvar_sop and mat2copvar_sop, which
%                replace the opvar2copvar(mat2opvar(...)) detour; that
%                detour has no N-D route.

if nargin<3 || isempty(mult_only),  mult_only = false;  end
isdec = isa(Mat,'dpvar');
if ~isdec && ~isnumeric(Mat)
    error('mat2copvar_grid:badClass',...
        "The matrix must be 'double' or 'dpvar', not '%s'.",class(Mat))
end
m = sum(meta.dim_out);      n = sum(meta.dim_in);
if ~isequal(size(Mat),[m,n])
    error('mat2copvar_grid:size',['The matrix is %dx%d, but the spaces have %d '...
        'output and %d input components.'],size(Mat,1),size(Mat,2),m,n)
end

% % % vec(Mat) = A + B'*d, column-major, d the decision variables.
if isdec
    % All of Mat's spatial variables go on the output side, so a degree
    % above 0 shows up in ZL: the matrix is not constant and is refused.
    vconv = struct('out',{reshape(Mat.varname,1,[])},'in',{cell(1,0)});
    [Pc,ZLc] = dpvar2sdvar(Mat,vconv);
    if any(cellfun(@(z) any(z(:)~=0),ZLc))
        error('mat2copvar_grid:spatial',['The dpvar depends on a spatial '...
            'variable; only matrices constant in space are supported here.'])
    end
    A = Pc.A;       B = Pc.B;       Zd = reshape(Mat.dvarname,[],1);
else
    A = sparse(double(Mat(:)));     B = sparse(0,m*n);      Zd = cell(0,1);
end

% % % The blocks, one slice of (A,B) each.
M = numel(meta.dim_out);    N = numel(meta.dim_in);
ro = cumsum([0;meta.dim_out(:)]);   co = cumsum([0;meta.dim_in(:)]);
C = cell(M,N);
for i = 1:M
    for j = 1:N
        r = ro(i)+(1:meta.dim_out(i));      c = co(j)+(1:meta.dim_in(j));
        lin = reshape(r(:) + m*(c(:).'-1),[],1);
        Aij = A(lin);       Bij = B(:,lin);
        if nnz(Aij)==0 && nnz(Bij)==0,      continue,   end
        C{i,j} = const_block(Aij,Bij,Zd,isdec,meta,i,j,mult_only);
    end
end
% One explicit zero block per empty row, then per empty column.
occ = ~cellfun(@isempty,C);
for i = find(~any(occ,2))'
    nc = meta.dim_out(i)*meta.dim_in(1);
    C{i,1} = const_block(sparse(nc,1),sparse(numel(Zd),nc),Zd,isdec,meta,i,1,false);
    occ(i,1) = true;
end
for j = find(~any(occ,1))
    nc = meta.dim_out(1)*meta.dim_in(j);
    C{1,j} = const_block(sparse(nc,1),sparse(numel(Zd),nc),Zd,isdec,meta,1,j,false);
end

% Metadata is the caller's, already validated, so the trusted constructor.
keep = {'vars','dom','space_out','space_in','dim_out','dim_in'};
mt = struct();
for k = 1:numel(keep),  mt.(keep{k}) = meta.(keep{k});  end
mt.dim_out = mt.dim_out(:);     mt.dim_in = mt.dim_in(:);
if isdec
    mt.Zd = Zd;
    P = cdopvar(C,mt);
else
    P = copvar(C,mt);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = const_block(Aij,Bij,Zd,isdec,meta,i,j,mult_only)
% The constant block from input space j to output space i, degree-0 bases,
% content in the all-multiplier cell only (gamma = 1 in every shared
% direction, linear index 1).

vars = meta.vars;
vo = vars(meta.space_out(i,:));     vi = vars(meta.space_in(j,:));   % sorted
S3 = vo(ismember(vo,vi));   S2 = vo(~ismember(vo,vi));  S1 = vi(~ismember(vi,vo));
if mult_only && ~isempty(S1) && (nnz(Aij)>0 || nnz(Bij)>0)
    error('mat2copvar_grid:integral',['Block (%d,%d) maps out of a space with '...
        'variable(s) {%s} the output space lacks, so a constant block there is '...
        'an integral, not a multiplier; refused here as mat2opvar refuses Q1.'],...
        i,j,strjoin(S1,','))
end
vout = [S2,S3];     vin = [S3,S1];      % the classes' canonical order
[~,io] = ismember(vout,vars);   [~,ii] = ismember(vin,vars);
bvars = struct('out',{vout},'in',{vin});
bdom = struct('out',meta.dom(io,:),'in',meta.dom(ii,:));
ZL = repmat({0},1,numel(vout));     ZR = repmat({0},1,numel(vin));
dims = [meta.dim_out(i),meta.dim_in(j)];    nc = dims(1)*dims(2);
shape = [3*ones(1,numel(S3)),1,1];
if isdec
    Ac = repmat({sparse(nc,1)},shape);      Bc = repmat({sparse(numel(Zd),nc)},shape);
    Ac{1} = Aij;    Bc{1} = Bij;
    b = sdopvar(struct('A',{Ac},'B',{Bc}),bvars,Zd,ZL,ZR,bdom,dims);
else
    Pc = repmat({sparse(dims(1),dims(2))},shape);
    Pc{1} = reshape(Aij,dims(1),dims(2));
    b = sopvar(Pc,bvars,ZL,ZR,bdom,dims);
end

end
