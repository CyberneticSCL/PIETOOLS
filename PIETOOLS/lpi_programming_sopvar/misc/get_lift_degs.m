function [deg,terms,info] = get_lift_degs(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DEG,TERMS,INFO] = GET_LIFT_DEGS(P,OPTS) returns the degrees and the
% Positivstellensatz terms of a positive container operator Pop for the
% constraint P = Pop, Pop >= 0, with P self-adjoint, by the following
% lift/weight rule of the separated form. In each registry
% direction d, with i and j the exponents of s_d and s'_d in a monomial of
% a kernel of P (the upper cells, adjoints of the lower ones, not read),
%
%   Dmin_d = the largest min(i,j) of the lower cells;
%   Dmax_d = the largest max(i,j) of the lower cells and the largest
%            exponent in d of a block between two L2 spaces of which one
%            lacks d;
%   Mdeg_d = the largest i of the multiplier cells (-1 if none);
%   J_d    = the largest exponent in d of a block between an R^q space and
%            an L2 space (-1 if none).
%
% A positive operator of lift degree D_d and weight degree w_d reaches the
% lower kernels with min(i,j) <= D_d and max(i,j) <= D_d + 2 w_d + 1, the
% multiplier degree 2 w_d and, between R^q and L2, the degree
% D_d + w_d + wR + 1; reaching the support is necessary for P = Pop, not
% sufficient. The rule sets D_d = Dmin_d + OPTS.dD and w_d to the smallest
% weight degree that reaches the support plus OPTS.dw,
%
%   w_d = max(ceil((Dmax_d-D_d-1)/2), ceil(Mdeg_d/2), J_d-D_d-wR-1, 0)
%         + OPTS.dw,
%
% raised to OPTS.wmin when the support is read from P itself (no
% OPTS.like).
%
% INPUT
% - P:     'dopvar' or 'opvar' object, square (the same output and input
%          spaces), or an 'sdopvar' or 'sopvar' object from a space to
%          itself (vars.out = vars.in);
% - opts:  (optional) 'struct' with the optional fields
%            dD, dw        double, added to Dmin_d and to the weight rule
%                          (default 1, 1);
%            wmin          double, the least weight degree when the support
%                          is read from P itself (default 0: none);
%            wR            double, the weight degree of an R^q space
%                          (default 0, the identity basis);
%            psatz         the terms: 'auto' (default; in one variable the
%                          codes [0 1] with offsets [0 1]; in nv ~= 1
%                          variables the codes [0, 3:2nv+2], the plain term
%                          and the 2nv faces, the faces at offset 1 when
%                          every w_d is at least 2, else at offset 0),
%                          'none' (code 0), 'product' (codes [0 1], offsets
%                          [0 1]), 'faces' (codes [0, 3:2nv+2], offsets 0),
%                          'faces1' (codes [0, 3:2nv+2], offsets [0, 1,
%                          ..., 1]: the faces one weight degree below the
%                          plain term for every w, the set 'auto' takes
%                          when every w_d >= 2; in one variable codes
%                          [0 1], offsets [0 1]),
%                          'corners' (the 'auto' terms and every face
%                          product of two directions, at the offset of the
%                          face terms),
%                          or a double row of term codes (0 plain, 1 box
%                          product, 2d+1 / 2d+2 lower / upper face of
%                          direction d, 2nv+2+d box quadratic of direction
%                          d);
%            psatz_offset  double, with a row of codes: the offset of each
%                          term, or one for all (default 0);
%            psatz_rows, psatz_rows_offset
%                          with a row of codes: R x nv rows of
%                          per-direction factor codes (0 none, 1 face
%                          (s_d-a_d), 2 face (b_d-s_d), 3 box quadratic),
%                          one term per row with the product of its
%                          factors as weight, and their offsets; returned
%                          in TERMS;
%            wcross        double (default []: w_d): the weight degree of
%                          an L2 space in every registry direction it does
%                          not own, in place of w_d;
%            dw_space      1 x M double, M the number of spaces of P
%                          (default zeros): added to the weight degree of
%                          space k in every direction, after the rule and
%                          before wcross; an R^q space ignores it;
%            Dmax          logical (default false): the lift degree D_d =
%                          max(Dmax_d, Dmin_d + dD) in place of
%                          Dmin_d + dD, with w_d from the rule at
%                          Dmin_d + dD;
%            cap_inc       'auto' or 1 x nv double (default []: none): the
%                          cap theta_d + s_d <= cap_inc_d on the monomials
%                          of the lower and upper basis operators, the
%                          field cap_inc of each L2 entry of DEG; 'auto' is
%                          Dmax_d - 1, floored at 0;
%            like          'dopvar' or 'opvar' object, or a cell of them,
%                          with variables among those of P: their pooled
%                          support is read instead of that of P, and DEG is
%                          laid out over the spaces of P.
%
% OUTPUT
% - deg:   1 x M cell, M the number of spaces of P (1 for a block): for an
%          L2 space a 'struct' with fields int (1 x nv double, w) and mult
%          (1 x nv double, D in the directions of the space, 0 in the
%          others), and cap_inc (1 x nv double) when OPTS.cap_inc is set;
%          for an R^q space a 'struct' with field int (1 x nv double, wR);
% - terms: 'struct' with fields codes and offsets, double rows of one
%          length, the term codes and the amounts by which each lowers the
%          weight degree, and rows (R x nv double) and rows_offsets (1 x R),
%          the factor rows and their offsets;
% - info:  'struct' with fields Dmin, Dmax, Mdeg, J, D, w (1 x nv double),
%          Jmax (1 x nv double, the largest total degree of a kernel
%          monomial in the two variables of direction d: i + j over the
%          lower cells, i over multiplier and one-sided cells; -1 where
%          unseen; the per-direction cap theta_d + s_d a lift basis needs
%          to carry the operator read), wR, dD, dw, wmin, nv, M (double),
%          vars (1 x nv cellstr, the registry), isR (1 x M logical, true
%          for an R^q space), Dmax_rule (logical, OPTS.Dmax) and cap_inc
%          (the caps applied, [] for none).
%
% NOTES
% The rule bounds the support of a fixed P. For P linear in a decision
% operator, such as A'*Q + Q'*A with Q free, the support of P is that of
% every value of the decision operator; to size from the positive operator
% of the LPI instead, pass it as OPTS.like.
%
% See also POSLPIVAR_DIRECT, LPI_INEQ_SOP, EQ_OPTS_SOPVAR,
%          GET_LPIVAR_DEGS_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - get_lift_degs
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
% Initial coding MMP, 10/08/2026

% Options; the defaults of dD, dw and wR
if nargin<2 || isempty(opts),   opts = struct();    end
dD = 1;     if isfield(opts,'dD') && ~isempty(opts.dD),     dD = opts.dD;   end
dw = 1;     if isfield(opts,'dw') && ~isempty(opts.dw),     dw = opts.dw;   end
wR = 0;     if isfield(opts,'wR') && ~isempty(opts.wR),     wR = opts.wR;   end
% wmin: the least weight degree read from P itself (default 0, none)
wmin = 0;   if isfield(opts,'wmin') && ~isempty(opts.wmin), wmin = opts.wmin;   end
useDmax = false;    capinc = [];    % the ceiling, the triangular cap
if isfield(opts,'Dmax') && ~isempty(opts.Dmax)
    useDmax = logical(opts.Dmax);
end
if isfield(opts,'cap_inc') && ~isempty(opts.cap_inc)
    capinc = opts.cap_inc;
end

% % % Blocks and spaces
if isa(P,'sdopvar') || isa(P,'sopvar')
    % One block: its sorted variables, one space
    vars = union(P.vars.out,P.vars.in);     vars = reshape(sort(vars),1,[]);
    blocks = {P};   bk = 1;  bkp = 1;
    own = {reshape(find(ismember(vars,P.vars.out)),1,[])};
    if ~isequal(sort(P.vars.out),sort(P.vars.in))
        error("A single block must map a space to itself (vars.out = vars.in).")
    end
elseif isa(P,'dopvar') || isa(P,'opvar')
    if ~isequal(P.space_out,P.space_in) || ~isequal(P.dim_out(:),P.dim_in(:))
        error("P must be square: the same spaces on both sides.")
    end
    % A container: its nonempty blocks with their space indices, and the
    % variables of each space
    vars = reshape(P.vars,1,[]);
    [Mo,Mi] = size(P.C);
    blocks = {};    bk = [];    bkp = [];
    for k = 1:Mo
        for kp = 1:Mi
            if ~isempty(P.C{k,kp}),  blocks{end+1} = P.C{k,kp};  bk(end+1) = k;  bkp(end+1) = kp;  end %#ok<AGROW>
        end
    end
    own = cell(1,Mo);
    for k = 1:Mo,   own{k} = reshape(find(P.space_out(k,:)),1,[]);   end
else
    error("P should be a 'dopvar', 'opvar', 'sdopvar' or 'sopvar'.")
end
% Registry size, number of spaces, the R^q spaces
nv = numel(vars);   M = numel(own);
isR = cellfun(@isempty,own);
% dw_space: a raise of the weight degree per space (default none)
dws = zeros(1,M);
if isfield(opts,'dw_space') && ~isempty(opts.dw_space)
    dws = reshape(opts.dw_space,1,[]);
    if numel(dws)~=M,   error("opts.dw_space needs one entry per space of P (%d).",M);   end
end
% % % opts.like: the support of other containers, over the spaces of P
isRs = isR;
if isfield(opts,'like') && ~isempty(opts.like)
    % One operator or a cell of them: the blocks pooled, so the support is
    % the union; the space indices of each operator are offset past those
    % of the previous ones for the R-space flags isRs
    Ls = opts.like;     if ~iscell(Ls),     Ls = {Ls};  end
    blocks = {};    bk = [];    bkp = [];   isRs = false(1,0);
    for j = 1:numel(Ls)
        L = Ls{j};      o = numel(isRs);
        if ~(isa(L,'dopvar') || isa(L,'opvar'))
            error("opts.like should be a 'dopvar' or 'opvar', or a cell of them.")
        end
        if ~all(ismember(reshape(L.vars,1,[]),vars))
            error("The variables of opts.like must be among those of P.")
        end
        [Mo,Mi] = size(L.C);
        for k = 1:Mo
            for kp = 1:Mi
                if ~isempty(L.C{k,kp}),  blocks{end+1} = L.C{k,kp};  bk(end+1) = o+k;  bkp(end+1) = o+kp;  end %#ok<AGROW>
            end
        end
        isRs = [isRs, reshape(~any(L.space_out,2),1,[])]; %#ok<AGROW>
    end
end

% % % The support, per direction
Dmin = zeros(1,nv);     Dmax = zeros(1,nv);     Mdeg = -ones(1,nv);     J = -ones(1,nv);
% Jmax_d: the largest total degree in (s_d, s'_d) of a kernel monomial,
% i + j over the lower cells, i over multiplier and one-sided cells (-1 if
% unseen): what a per-direction cap theta_d + s_d of a lift basis must
% span to carry the operator.
Jmax = -ones(1,nv);
for b = 1:numel(blocks)
    % Block b: its variables and their registry indices, its monomial
    % grids, dimensions and shared variables S3
    B = blocks{b};
    vout = reshape(B.vars.out,1,[]);    vin = reshape(B.vars.in,1,[]);
    [~,po] = ismember(vout,vars);       [~,pi_] = ismember(vin,vars);
    EL = multiindex_grid(B.ZL,'first_slowest');     NL = size(EL,1);
    ER = multiindex_grid(B.ZR,'first_slowest');     NR = size(ER,1);
    m = B.dims(1);  n = B.dims(2);
    S3 = intersect(vout,vin);   n3 = numel(S3);
    [~,p3o] = ismember(S3,vout);    [~,p3i] = ismember(S3,vin);
    [~,d3] = ismember(S3,vars);
    % True for a block between an R^q space and an L2 space
    otherR = (isRs(bk(b)) && ~isRs(bkp(b))) || (~isRs(bk(b)) && isRs(bkp(b)));
    % The number of parameter cells of the block B: params.A of an
    % 'sdopvar', params of a 'sopvar'
    if isa(B,'sdopvar'),    ncell = numel(B.params.A);     else,   ncell = numel(B.params);  end
    for q = 1:ncell
        % The monomial pairs of the nonzeros of cell q, and its multi-index
        [li,ei] = cell_support(B,q,m,n,NL,NR);
        if isempty(li),     continue,   end
        gam = gamma_of_cell(q,n3);
        % shared directions: lower cells give (i,j), multiplier cells i
        for t = 1:n3
            d = d3(t);  i = EL(li,p3o(t));   j = ER(ei,p3i(t));
            if gam(t)==1
                Mdeg(d) = max([Mdeg(d); i(:)]);
                Jmax(d) = max([Jmax(d); i(:)]);
            elseif gam(t)==2
                Dmin(d) = max([Dmin(d); min(i(:),j(:))]);
                Dmax(d) = max([Dmax(d); max(i(:),j(:))]);
                Jmax(d) = max([Jmax(d); i(:)+j(:)]);
            end
        end
        % one-sided directions: the exponent of the side that has them
        for t = 1:numel(vout)
            if ismember(vout{t},S3),    continue,   end
            d = find(strcmp(vars,vout{t}));  i = EL(li,t);
            if otherR,  J(d) = max([J(d); i(:)]);  else,   Dmax(d) = max([Dmax(d); i(:)]);  end
            Jmax(d) = max([Jmax(d); i(:)]);
        end
        for t = 1:numel(vin)
            if ismember(vin{t},S3),     continue,   end
            d = find(strcmp(vars,vin{t}));   j = ER(ei,t);
            if otherR,  J(d) = max([J(d); j(:)]);  else,   Dmax(d) = max([Dmax(d); j(:)]);  end
            Jmax(d) = max([Jmax(d); j(:)]);
        end
    end
end

% % % The rule
D = Dmin + dD;
w = max([ceil((Dmax-D-1)/2); ceil(Mdeg/2); J-D-wR-1; zeros(1,nv)],[],1) + dw;
% The ceiling after the weight rule: w from the base lift Dmin + dD, so the
% family at D = Dmax contains the base family
if useDmax,     D = max(Dmax,D);     end
% Read from P itself, at least wmin: with dw alone a fixed target whose
% support needs w = 0 gets w = 1
if ~(isfield(opts,'like') && ~isempty(opts.like)),  w = max(w,wmin);  end
% cap_inc: 'auto' = Dmax_d - 1 (floored at 0), else one cap per direction
if ischar(capinc) || isstring(capinc),  capinc = max(Dmax-1,0);    end
if ~isempty(capinc)
    capinc = reshape(capinc,1,[]);
    if isscalar(capinc),    capinc = repmat(capinc,1,nv);    end
end
% Per space: w and the lift D in its own directions, or wR on R^q
deg = cell(1,M);
for k = 1:M
    if isR(k)
        deg{k} = struct('int',wR*ones(1,nv));
    else
        mult = zeros(1,nv);     mult(own{k}) = D(own{k});
        % w plus dw_space(k); then wcross in the directions k does not own,
        % where the space has no lift: wcross 0 drops those slack monomials
        wk = w + dws(k);
        if isfield(opts,'wcross') && ~isempty(opts.wcross)
            wk(~ismember(1:nv,own{k})) = opts.wcross;
        end
        deg{k} = struct('int',wk,'mult',mult);
        if ~isempty(capinc),    deg{k}.cap_inc = capinc;    end
    end
end

% % % The terms
ps = 'auto';
if isfield(opts,'psatz') && ~isempty(opts.psatz),   ps = opts.psatz;    end
% The face products: a face (s_i-a_i) or (b_i-s_i) of direction i times
% one of direction j > i, as rows of per-direction factors, the products
% the codes cannot express; none in one variable
cr = zeros(2*nv*(nv-1),nv);    k = 0;  rows = zeros(0,nv);   roffs = zeros(1,0);
for i = 1:nv
    for j = i+1:nv
        for e = 1:2
            for f = 1:2
                k = k+1;    cr(k,:) = 0;    cr(k,i) = e;    cr(k,j) = f;
            end
        end
    end
end
% A row of codes with its offsets, or a named set of terms
if isnumeric(ps)
    codes = reshape(ps,1,[]);   offs = zeros(size(codes));
    if isfield(opts,'psatz_offset') && ~isempty(opts.psatz_offset)
        offs = reshape(opts.psatz_offset,1,[]);
        if isscalar(offs),  offs = repmat(offs,size(codes));    end
    end
    % Rows given with the codes
    if isfield(opts,'psatz_rows') && ~isempty(opts.psatz_rows)
        rows = opts.psatz_rows;     roffs = zeros(1,size(rows,1));
        if isfield(opts,'psatz_rows_offset') && ~isempty(opts.psatz_rows_offset)
            roffs = reshape(opts.psatz_rows_offset,1,[]);
        end
    end
else
    switch lower(char(ps))
        case 'auto'
            if nv==1,   codes = [0 1];  offs = [0 1];
            else        % faces one below the plain term when every weight allows it, else at w
                codes = [0, 3:2*nv+2];  offs = [0, double(min(w)>=2)*ones(1,2*nv)];
            end
        case 'none',    codes = 0;              offs = 0;
        case 'product', codes = [0 1];          offs = [0 1];
        case 'faces',   codes = [0, 3:2*nv+2];  offs = zeros(1,2*nv+1);
        case 'faces1'   % the faces one weight degree below the plain term, for every w
            if nv==1,   codes = [0 1];  offs = [0 1];   % as 'auto', 'product'
            else,       codes = [0, 3:2*nv+2];  offs = [0, ones(1,2*nv)];   end
        case 'corners'  % the 'auto' terms and every two-direction face product
            % at the offset of the face terms
            codes = [0, 3:2*nv+2];  offs = [0, double(min(w)>=2)*ones(1,2*nv)];
            rows = cr;              roffs = offs(2)*ones(1,size(cr,1));
        otherwise
            error("'psatz' should be 'auto', 'none', 'product', 'faces', 'faces1', 'corners' or a row of term codes.")
    end
end
% The outputs
terms = struct('codes',codes,'offsets',offs,'rows',rows,'rows_offsets',roffs);
info = struct('Dmin',Dmin,'Dmax',Dmax,'Mdeg',Mdeg,'J',J,'Jmax',Jmax,'D',D,'w',w,'wR',wR,...
              'dD',dD,'dw',dw,'wmin',wmin,'nv',nv,'M',M,'vars',{vars},'isR',isR);
info.Dmax_rule = useDmax;   info.cap_inc = capinc;

end


function [li,ei] = cell_support(B,q,m,n,NL,NR)
% [LI,EI] = CELL_SUPPORT(B,Q,M,N,NL,NR) returns the left and right monomial
% indices of every nonzero of parameter cell Q of the block B, over the
% constant part and every decision-variable row.
% INPUT  B:      'sdopvar' or 'sopvar' object;
%        q:      double, the linear index of the parameter cell;
%        m, n:   double, the output and input dimensions of B;
%        NL, NR: double, the numbers of monomials of B.ZL and B.ZR.
% OUTPUT li, ei: double columns, the indices (1..NL and 1..NR) of the
%                left and right monomials of the nonzeros; [] for a zero
%                cell.
if isa(B,'sdopvar')
    % Coefficient positions with a nonzero constant or decision row
    A = B.params.A{q};  Bq = B.params.B{q};
    cols = find(A);     cols = cols(:);
    if ~isempty(Bq),    [~,cB] = find(Bq);   cols = unique([cols; cB(:)]);  end
    if isempty(cols),   li = [];  ei = [];   return,     end
    [ridx,cidx] = ind2sub([m*NL,n*NR],cols);
else
    % The nonzeros of the coefficient matrix
    Cq = B.params{q};
    if isempty(Cq) || nnz(Cq)==0,   li = [];  ei = [];   return,     end
    [ridx,cidx] = find(Cq);
end
% Rows are (matrix row, monomial), columns (matrix column, monomial)
li = mod(ridx-1,NL)+1;  ei = mod(cidx-1,NR)+1;
end
