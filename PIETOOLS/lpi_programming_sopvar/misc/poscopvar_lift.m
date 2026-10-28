function [prog,Pop,M,A,info] = poscopvar_lift(prog,dims,spaces,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP,M,A,INFO] = POSCOPVAR_LIFT(PROG,DIMS,SPACES,DOM,DEG,OPTIONS)
% declares a positive semidefinite, self-adjoint 'cdopvar' decision
% operator on X = X_1 x ... x X_Ms, X_k = R^{m_k} or L_2^{m_k}[s^k], as
%
%       Pop = A^* M A,
%
% A the fixed lift from X to Y = R^{m_R} x L_2^{K_1}[s] x ... x
% L_2^{K_G}[s], s the registry (the sorted union of the variables of the
% s^k), and M a positive decision operator on Y. A places every R^{m_k} in
% R^{m_R} and maps L_2^{m_k}[s^k] to its moments, one per multi-index
% alpha over the variables of s^k,
%
%       (A_{k,alpha} x)(theta) = int I_alpha(theta-s)
%                                    (I_{m_k} kron V^alpha(s)) x(s) ds,
%       I_alpha(r) = prod_d I_{alpha_d}(r_d),     I_1(r) = delta(r),
%       I_2(r) = [r >= 0],     I_3(r) = [r <= 0],     I_4(r) = 1,
%       V^alpha(s) = kron_{d: alpha_d ~= 1} [1; s_d; ...; s_d^{D_d}],
%
% the integral over the domain of s^k, [.] the indicator function, d over
% the variables of s^k in registry order, alpha_d in {1,2,3} ({1,4} in a
% separable direction) and D the lift degree; the moments with equal
% weight degree w form one space L_2^{K_g}[s]. M is
%
%       M = sum_c Lambda_c^* g_c Q_c Lambda_c,     Q_c >= 0,
%
% one Gram matrix Q_c of new decision variables per term c of
% OPTIONS.psatz, g_c the multiplication by its weight g_c(theta), and
% Lambda_c the block-diagonal map y -> (I kron Z_{w_c}(theta)) y over the
% spaces of Y, Z_{w_c}(theta) the column of the monomials of degree at
% most w_c in each theta_d, w_c = max(w - OPTIONS.psatz_offset(c),0) with
% w the weight degree of the space, and y in R^{m_R} constant in theta.
%
% INPUT
% - prog:    'struct', the LPI program to which the Gram matrices Q_c are
%            added;
% - dims:    'double', Ms x 1 positive integers, dims(k) = m_k; a scalar
%            for every space;
% - spaces:  1 x Ms 'cell' of 'cellstr', spaces{k} the variable names of
%            s^k, an empty entry the space R^{m_k}, no name ending in
%            '_int' or '_dum'; a 'cellstr' is one space;
% - dom:     'double', nv x 2, row d the interval [a_d,b_d] of variable d
%            of the registry; or 1 x 2, the interval of every variable; or
%            a 'struct' with fields vars ('cellstr') and dom ('double', one
%            row per name);
% - deg:     the lift degree D and the weight degree w: a 'double' scalar
%            d, D = w = d; or a 'struct' with the fields
%              lift (or mult)    'double', scalar or 1 x nv, D (default w);
%              weight (or int)   'double', scalar or 1 x nv, w (default 1
%                                on L_2^{m_k}[s^k], 0 on R^{m_k});
%              joint             (optional) 'double' scalar, at least
%                                sum(D) + sum(w);
%            or a 1 x Ms 'cell' with one such entry per space, each of
%            which may be a 'cell' with one entry per component of the
%            space, in the order of INFO.basis_list; nonnegative integers.
%            On R^{m_k} only w is read, and all R^{m_k} have the same w;
% - options: (optional) 'struct' with the optional fields
%              sep           'logical', scalar or 1 x nv (default false),
%                            the separable directions;
%              include       1 x Ms 'cell', entry k the components of space
%                            k (default all, in the enumeration of every
%                            alpha, direction 1 fastest): a 'double' array
%                            with one column per variable of s^k, its rows
%                            distinct admissible alpha; a 'logical' mask
%                            over the enumeration; or other 'double'
%                            indices into it;
%              psatz         'double' row of term codes (default 0): 0,
%                            g_c = 1; 1, g_c(theta) = prod_d
%                            (theta_d-a_d)(b_d-theta_d); 2d+1 or 2d+2, d
%                            indexing the registry, g_c(theta) =
%                            (theta_d-a_d)/(b_d-a_d) or
%                            (b_d-theta_d)/(b_d-a_d);
%              psatz_offset  'double', scalar or one entry per code, the
%                            reduction of w in every direction (default 1
%                            for code 1, 0 otherwise);
%              wR            'double', scalar or 1 x nv, the weight degree
%                            of every R^{m_k}, in place of that of DEG;
%              type          'char', 'pos'.
%
% OUTPUT
% - prog:    'struct', PROG with the Gram matrices Q_c and the constraints
%            Q_c >= 0;
% - Pop:     'cdopvar' object on X, Ms x Ms blocks;
% - M:       'cdopvar' object on Y, the spaces INFO.dims_Y, INFO.spaces_Y;
% - A:       'copvar' object from X to Y without decision variables: block
%            row 1 R^{m_R} when X has an R^{m_k} space, then one block row
%            per space L_2^{K_g}[s];
% - info:    'struct' with fields
%              basis_list  'double', N_b x (1+nv), row c = [k, alpha] of
%                          component c, alpha over the registry with 0
%                          where s^k lacks the variable, ordered by space;
%              comp        1 x N_b 'struct' array, per component the fields
%                          space (k), alpha (over the variables of s^k),
%                          dim (its number of rows of Y), lift (1 x nv, D),
%                          weight (1 x nv, w), yspace (its block row of A),
%                          off (its first row in that space, minus 1), ZR
%                          (1 x n_k 'cell', the exponents of V^alpha in
%                          each variable of s^k) and cells (the parameter
%                          cells it fills in its block of A);
%              dims_Y      'double' column, the dimension of each space of
%                          Y;
%              spaces_Y    'cell' row, the variable names of each space of
%                          Y (empty for R^{m_R});
%              isR_Y       'logical' row, isR_Y(1) true when Y has R^{m_R},
%                          every other entry false;
%              vars, dom   1 x nv 'cellstr', the registry, and nv x 2
%                          'double', its intervals;
%              nv, M, nb   'double', the numbers of registry variables, of
%                          spaces of X and of components per space
%                          (1 x Ms);
%              sep         1 x nv 'logical', the separable directions;
%              mult        'struct' of the terms of M, with fields psatz
%                          and offset (1 x nterms 'double', the codes and
%                          the reductions of w), deg (1 x nterms 'cell',
%                          per term a 'cell' of 'struct' with field int,
%                          w_c on each space of Y), Qcell (1 x nterms
%                          'cell', per term the 'cell' of 'cellstr' naming
%                          the entries of the blocks of Q_c), gram
%                          (1 x nterms 'double', the order of each Q_c)
%                          and ndec;
%              ndec        'double', the number of decision variables of M;
%              t_lift, t_mult, t_compose  'double', the seconds taken to
%                          build A, M and A^* M A.
%
% NOTES
% A field subset of DEG and an OPTIONS.type other than 'pos' are errors.
%
% See also LIFT_COPVAR, POSMULT_CDOPVAR, POSCOPVAR, COPQUADVAR,
% GETSOL_LPIVAR_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - poscopvar_lift
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
% Initial coding MMP, 10/07/2026. The separated form A'*M*A proposed from
%                the proof program's representation P = A_D^* M_W A_D:
%                lift and weight as separate inputs, the weight the
%                parameter. Built beside 'poscopvar', which is unchanged,
%                so that the two can be measured against each other.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body and headers
%                of the local functions; code unchanged. The deleted
%                header text was the reading of the degrees in the proof
%                program (with references to notes and a memo), a measured
%                equivalence with 'poscopvar', and remarks on cost and N-D.

% Five inputs at least; OPTIONS a 'struct', type 'pos' only
if nargin<5
    error("Not enough input arguments.")
end
if nargin<6 || isempty(options),    options = struct(); end
if ~isa(options,'struct')
    error("Options should be specified as a 'struct' object.")
end
if isfield(options,'type') && ~isempty(options.type) && ~strcmp(char(options.type),'pos')
    error("'poscopvar_lift' declares a positive operator; the multiplier M is "...
          +"positive by construction, so type '%s' has no meaning here.",char(options.type))
end

% % % Spaces, to read a per-space 'deg'.
[meta,sp] = parse_copvar_spaces(dims,spaces,dom);
Ms = numel(sp);     nv = numel(meta.vars);
isR = cellfun(@isempty,sp);

% % % Degrees: lift and weight per space, or per component (cells).
[dl_spec,w_spec] = parse_deg(deg,Ms,isR,nv,options);

% % % Weight labels: one Y space per distinct weight row among the L2
% components; every R component shares one weight.
wrows = zeros(0,nv);    yg = cell(1,Ms);
wR = [];
for k = 1:Ms
    wk = w_spec{k};
    if ~iscell(wk),     wk = {wk};  end
    % Per component: label 0 on R (its weight kept as wR), else the index
    % of its weight row in wrows
    lab = zeros(1,numel(wk));
    for i = 1:numel(wk)
        if isR(k)
            if isempty(wR),     wR = wk{i};
            elseif ~isequal(wR,wk{i})
                error("Every R^q space must have the same weight degree (one R space of the lift).")
            end
            lab(i) = 0;
        else
            [tf,loc] = ismember(wk{i},wrows,'rows');
            if ~tf,     wrows(end+1,:) = wk{i};     loc = size(wrows,1);     end %#ok<AGROW>
            lab(i) = loc;
        end
    end
    % One label for every component of the space, or one per component
    % (a cell, which 'lift_copvar' reads per component).
    if numel(lab)==1,   yg{k} = lab;    else,   yg{k} = num2cell(lab);  end
end

% % % The lift
% The labels as groups, sep and include as given; timed
lopt = struct('ygroup',{yg});
if isfield(options,'sep'),      lopt.sep = options.sep;         end
if isfield(options,'include'),  lopt.include = options.include; end
t0 = tic;
[A,info] = lift_copvar(dims,spaces,dom,dl_spec,lopt);
t_lift = toc(t0);

% % % The weight degrees of Y: the R space, then the groups in order of
% first appearance, which is the order 'lift_copvar' numbers them in.
nY = numel(info.dims_Y);
wY = cell(1,nY);
if info.isR_Y(1),   wY{1} = wR;     end
for g = 1:size(wrows,1),    wY{g+info.isR_Y(1)} = wrows(g,:);  end
% The weight of each component, from its Y space
for c = 1:numel(info.comp)
    y = info.comp(c).yspace;    info.comp(c).weight = wY{y};
end

% % % The positive multiplier
% On the spaces of Y with their weights; psatz and psatz_offset as given;
% timed
mopt = struct();
if isfield(options,'psatz'),        mopt.psatz = options.psatz;                 end
if isfield(options,'psatz_offset'), mopt.psatz_offset = options.psatz_offset;   end
t0 = tic;
[prog,M,minfo] = posmult_cdopvar(prog,info.dims_Y,info.spaces_Y,info.dom,wY,mopt);
t_mult = toc(t0);

% % % The composition. M*A first: multiplier times lift, no integration;
% then the adjoint lift, which integrates theta out.
t0 = tic;
Pop = A'*(M*A);
t_compose = toc(t0);

% The terms of M, its number of decision variables and the build times
info.mult = minfo;      info.ndec = minfo.ndec;
info.t_lift = t_lift;   info.t_mult = t_mult;   info.t_compose = t_compose;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dl,ws] = parse_deg(deg,Ms,isR,nv,options)
% dl{k}, ws{k}: the lift and weight of space k, each a 1 x nv row or, when
% 'deg' gives one entry per component, a cell of rows.
% INPUT  deg:     the DEG input of POSCOPVAR_LIFT;
%        Ms:      'double', the number of spaces;
%        isR:     1 x Ms 'logical', true for an R^{m_k} space;
%        nv:      'double', the number of registry variables;
%        options: 'struct', the field wR read.
% OUTPUT dl, ws:  1 x Ms 'cell' of 1 x nv 'double' rows, or of 'cell'
%                 arrays of rows, one per component.
% One entry per space
if iscell(deg)
    if numel(deg)~=Ms
        error("A cell 'deg' should have one entry per space.")
    end
    spec = reshape(deg,1,[]);
else
    spec = repmat({deg},1,Ms);
end
% Per space: one entry or one per component; wR replaces the weight on R
dl = cell(1,Ms);    ws = cell(1,Ms);
for k = 1:Ms
    sk = spec{k};
    if iscell(sk)
        dl{k} = cell(1,numel(sk));  ws{k} = cell(1,numel(sk));
        for i = 1:numel(sk)
            [dl{k}{i},ws{k}{i}] = one_deg(sk{i},isR(k),nv,iscell(deg));
        end
    else
        [dl{k},ws{k}] = one_deg(sk,isR(k),nv,iscell(deg));
    end
    if isR(k) && isfield(options,'wR') && ~isempty(options.wR)
        wr = expand_row(options.wR,nv,'wR');
        if iscell(ws{k}),   ws{k} = repmat({wr},size(ws{k}));   else,   ws{k} = wr;     end
    end
end
end


function [dl,w] = one_deg(d,isRk,nv,percell)
% A scalar: lift = weight = scalar (copquadvar's scalar 'deg'); a struct:
% weight (alias int) default 1, lift (alias mult) default weight. An R
% space reads only the weight, default 0 for a struct or a per-space cell.
% INPUT  d:       one entry of DEG, a 'double' scalar or a 'struct';
%        isRk:    'logical', true for an R^{m_k} space;
%        nv:      'double', the number of registry variables;
%        percell: 'logical', true when DEG is a 'cell'.
% OUTPUT dl, w:   1 x nv 'double', the lift and the weight degree.
% Scalar: lift = weight
if isnumeric(d) && isscalar(d)
    dl = expand_row(d,nv,'deg');    w = dl;
    return
end
if ~isa(d,'struct')
    error("Degrees should be specified as a scalar or a 'struct' with fields 'lift' and 'weight'.")
end
% Struct: weight (or int), then lift (or mult), with their defaults
w = [];
if isfield(d,'weight') && ~isempty(d.weight),   w = d.weight;
elseif isfield(d,'int') && ~isempty(d.int),     w = d.int;
end
dl = [];
if isfield(d,'lift') && ~isempty(d.lift),       dl = d.lift;
elseif isfield(d,'mult') && ~isempty(d.mult),   dl = d.mult;
end
if isempty(w)
    if isRk && (percell || isstruct(d)),    w = 0;  else,   w = 1;  end
end
w = expand_row(w,nv,'weight');
if isempty(dl),     dl = w;     end
dl = expand_row(dl,nv,'lift');
% No degree cap: subset refused, joint at least sum(dl) + sum(w)
if isfield(d,'subset') && ~isempty(d.subset)
    error("A 'subset' cap grades the basis of 'copquadvar'; the separated form has no cap. Use 'poscopvar'.")
end
if isfield(d,'joint') && ~isempty(d.joint)
    if ~isscalar(d.joint) || d.joint<sum(dl)+sum(w)
        error("A 'joint' cap below the sum of lift and weight grades the basis of "...
              +"'copquadvar'; the separated form has no cap. Use 'poscopvar'.")
    end
end
end


function v = expand_row(v,nv,name)
% V = EXPAND_ROW(V,NV,NAME) returns V as a row of nonnegative integers,
% one entry per registry variable, a scalar repeated.
% INPUT  v:    'double', a scalar or nv entries;
%        nv:   'double', the number of registry variables;
%        name: 'char', the name used in the error messages.
% OUTPUT v:    1 x max(nv,1) 'double'.
if ~isnumeric(v) || any(v(:)<0) || any(v(:)~=round(v(:)))
    error("'%s' should be nonnegative integers.",name)
end
v = reshape(v,1,[]);
if isscalar(v)
    v = repmat(v,1,max(nv,1));
elseif numel(v)~=nv
    error("'%s' should be a scalar or have one entry per registry variable.",name)
end
end
