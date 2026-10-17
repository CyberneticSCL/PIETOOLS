function [prog,Pop,M,A,info] = poscopvar_lift(prog,dims,spaces,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP,M,A,INFO] = POSCOPVAR_LIFT(PROG,DIMS,SPACES,DOM,DEG,OPTIONS)
% declares a positive semidefinite, self-adjoint 'cdopvar' decision operator
% on a concatenation of mixed spaces in the separated form
%
%       Pop = A' * M * A,
%
% A the FIXED moment lift of 'lift_copvar' and M the positive multiplier of
% 'posmult_cdopvar'. The only decision object is the pointwise PSD
% polynomial matrix M(theta); Pop is its composition with the lift. With one
% Gram term this is the quadratic form of 'poscopvar' at 'int' = weight,
% 'mult' = lift and no joint cap (measured span-identical in
% test_poscopvar_lift). What the separation changes is what is exposed: the
% lift degree, the weight degree and the Positivstellensatz terms are
% separate inputs, and M is returned so that the solved weight - the
% certificate W(theta) of the proof program - can be read with
% 'getsol_lpivar_sop'.
%
% DEGREES, the two knobs of the proof program (proof_roadmap.md Sec. 1; the
% slot map is measured in the memory note degree-map-proof-to-code):
%   lift    D, the degree of the moments int_a^s t^j x(t) dt in the
%           components' own variables ('copquadvar's 'mult', poslpivar's
%           slot 2). For an exact equality Deop - P = 0 it must satisfy
%           D >= max min(i,j) over the monomials s^i t^j of P's lower
%           kernel, for every weight; D at least the kernel degree in each
%           variable is where the proof program's existence results hold;
%   weight  w, the half-degree of M in theta ('int', slot 1), the refinement
%           knob. The proof program bounds it by nothing; a polynomial
%           weight at finite w exists under a margin P >= delta*e_D, and its
%           finiteness in general is the necessity ASSUMPTION the
%           maintainer adopted on 10/07/2026 (1-D; N-D is a further
%           conjecture).
%
% INPUT
% - prog, dims, spaces, dom: as 'poscopvar';
% - deg:     scalar s: lift = weight = s on every space, as a scalar 'deg'
%            in 'copquadvar'; or a struct with fields
%              lift    (alias mult)  scalar or 1 x nv, defaults to weight;
%              weight  (alias int)   scalar or 1 x nv, defaults to 1;
%              joint   accepted only when it imposes no cap, i.e. is at
%                      least the sum of lift and weight over the space's
%                      and the registry's variables; a cap is copquadvar's
%                      graded basis and is refused here;
%              subset  refused;
%            or a cell with one such entry per space, each of which may be
%            a cell with one per component (order of INFO.basis_list). On
%            an R^{m_k} space only the weight is read; it defaults to 0,
%            the identity basis of 'poslpivar', when 'deg' is a struct or
%            a cell, and to s for a scalar 'deg' (copquadvar's reading);
% - options: (optional) struct with fields
%   sep, include   as 'copquadvar': the components of A;
%   psatz          row of term codes for M, as 'posmult_cdopvar'; default
%                  0. NOTE 'poscopvar's psatz = 1 is the single product
%                  term; the plain-plus-product pair of 'lpi_ineq' is
%                  psatz = [0 1], the product term one degree lower;
%   psatz_offset   as 'posmult_cdopvar';
%   wR             weight degree of the R components (scalar or 1 x nv),
%                  overriding the per-space value;
%
% OUTPUT
% - prog:    the program with M's Grams declared and the constraints
%            Q_c >= 0;
% - Pop:     M x M 'cdopvar' on the given spaces, Pop = Pop' >= 0;
% - M:       the positive multiplier, a 'cdopvar' on the lifted space
%            INFO.dims_Y x INFO.spaces_Y;
% - A:       the lift, a 'copvar' from the given spaces to the lifted one;
% - info:    the struct of 'lift_copvar' (basis_list, comp, dims_Y, ...),
%            comp with the added field weight, plus mult (the info of
%            'posmult_cdopvar'), ndec, and the build times t_lift, t_mult,
%            t_compose in seconds.
%
% NOTES
% Y groups the L2 components by their weight degree, so with one weight
% for all components (the usual case) A has at most two block rows and
% Pop = A'*(M*A) is a few compositions of large blocks instead of the
% N_b^2 pair terms of 'copquadvar'. Which is faster at which size is
% measured, not assumed: see test_poscopvar_lift.
%
% N-D: the lift is the tensor product of the 1-D lifts and M a polynomial
% matrix in all theta on the box, the structure of 'copquadvar'. The proof
% program is 1-D; nothing here is asserted beyond that in N-D.
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
for c = 1:numel(info.comp)
    y = info.comp(c).yspace;    info.comp(c).weight = wY{y};
end

% % % The positive multiplier
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

info.mult = minfo;      info.ndec = minfo.ndec;
info.t_lift = t_lift;   info.t_mult = t_mult;   info.t_compose = t_compose;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dl,ws] = parse_deg(deg,Ms,isR,nv,options)
% dl{k}, ws{k}: the lift and weight of space k, each a 1 x nv row or, when
% 'deg' gives one entry per component, a cell of rows.
if iscell(deg)
    if numel(deg)~=Ms
        error("A cell 'deg' should have one entry per space.")
    end
    spec = reshape(deg,1,[]);
else
    spec = repmat({deg},1,Ms);
end
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
if isnumeric(d) && isscalar(d)
    dl = expand_row(d,nv,'deg');    w = dl;
    return
end
if ~isa(d,'struct')
    error("Degrees should be specified as a scalar or a 'struct' with fields 'lift' and 'weight'.")
end
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
