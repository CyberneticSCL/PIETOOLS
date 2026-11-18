function [prog,Pop,Qcell,basis_list] = poscopvar(prog,dims,spaces,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP,QCELL,BASIS_LIST] = POSCOPVAR(PROG,DIMS,SPACES,DOM,DEG,OPTIONS)
% declares in the LPI program PROG a positive semidefinite, self-adjoint
% 'cdopvar' decision operator on a concatenation of spaces, parameterized
% as in Sec. 8.4 of the sopvar document,
%
%       Pop: X -> X,     X = L_2^{m_1}[s^1] x ... x L_2^{m_M}[s^M],
%
%       Pop = sum_{b,b'} Z_b'*(g*Q_{b,b'})*Z_{b'},     Q >= 0,
%
% with ONE Gram matrix Q of decision variables over every pair (b,b') of
% basis operators of all spaces, so that Pop is positive as a whole, and
% g(theta) >= 0 the weight of OPTIONS.psatz. An empty s^k is the space
% R^{m_k}. Basis operator b of space k carries a multi-index alpha over the
% variables of s^k and maps L_2^{m_k}[s^k] to L_2[theta], theta over every
% variable of the registry:
%
%       (Z_b x)(theta) = int I_alpha(theta-s) ...
%                            (I_{m_k} kron Z^alpha(theta,s)) x(s) ds,
%
% I_alpha(r) = prod_d I_{alpha_d}(r_d) over the variables of s^k, alpha_d =
% 1 the multiplier (delta), 2 the lower integral (s_d <= theta_d), 3 the
% upper integral (s_d >= theta_d), 4 the full integral; Z^alpha(theta,s)
% is the vector of monomials of the degrees DEG, and on R^{m_k},
% Z_b x = (I_{m_k} kron Z^alpha(theta)) x. Pop is also appended to the
% cell prog.sop_pos, from which LPI_INEQ_SOP sizes an inequality that
% depends on Pop.
%
% [PROG,POP,QCELL,BASIS_LIST] = POSCOPVAR(PROG,X,DEG,OPTIONS) declares Pop
% on the output spaces of the 'copvar' or 'cdopvar' object X, with the
% domain of X, so that X'*Pop*X is defined: for a PIE, X = T gives the P
% of T'*P*T. For the input spaces of X, pass X'.
%
% INPUT
% - prog:        'struct', the LPI program (LPIPROGRAM);
% - dims:        double, M x 1 positive integers, dims(k) = m_k the number
%                of components of space k; a scalar is used for every
%                space; empty for no space (POP is then 0 x 0, with no
%                decision variable, and PROG is returned unchanged);
% - spaces:      1 x M cell of 'cellstr', spaces{k} the names of the
%                variables of s^k ({} for R^{m_k}); a plain 'cellstr' is
%                one space, so several spaces are nested, as {{},{'s'}};
% - dom:         (optional) double, nv x 2, dom(d,:) the interval of
%                variable d of the registry (the sorted union of all
%                spaces), or one 1 x 2 row for every variable; or a
%                'struct' with fields vars ('cellstr') and dom (double, one
%                row per name). Omitted, as POSCOPVAR(PROG,DIMS,SPACES,DEG,
%                OPTIONS): the domain of PROG (LPIPROGRAM), by name;
% - X:           'copvar' or 'cdopvar' object, in place of DIMS, SPACES
%                and DOM (second form above);
% - deg:         (optional; omitted or [], the storage row of the program   % MMP, 11/17/2026
%                settings: DEG OMITTED below)                               % MMP, 11/17/2026
%                the degrees of Z^alpha(theta,s): a scalar (int and mult),
%                or a 'struct' with the optional fields
%                  int     scalar or 1 x nv, the degree in each theta_d
%                          (default 1);
%                  mult    scalar or 1 x nv, the degree in each s_d of
%                          s^k, 0 where alpha_d = 1 (default int);
%                  joint   scalar, a cap on the total degree (default
%                          none);
%                  subset  2^(2nv) caps, one per subset of [theta, s^k],
%                          for a space with every registry variable only;
%                or a 1 x M cell of these, one per space, each of which
%                may be a cell with one per basis operator of the space;
% - options:     (optional) 'struct' with the optional fields
%                  psatz    double, the weight g: 0 (default, g = 1), 1
%                           (g = prod_d (theta_d-a_d)*(b_d-theta_d)), or
%                           2d+1 / 2d+2, d = 1..nv (g = (theta_d-a_d)/L_d
%                           or (b_d-theta_d)/L_d, L_d = b_d-a_d); any
%                           other value is an error;
%                  psatz_rows, psatz_rows_offset
%                           rows of per-direction factors and their
%                           offsets, passed on to POSCOPVAR_DIRECT (its
%                           help): the products of faces in different
%                           directions;
%                  sep      logical, scalar or 1 x nv: where true, the
%                           lower and upper basis operators of that
%                           direction are replaced by the full integral
%                           (alpha_d = 4);
%                  include  1 x M cell, entry k selecting basis operators
%                           of space k: multi-index rows over the
%                           variables of s^k, a logical mask, or linear
%                           indices (default all).
%                                                                           % MMP, 11/17/2026
% DEG OMITTED (or []): the storage row lpi.pos of the program settings      % MMP, 11/17/2026
% (LPIPROGRAM; settings report 11/16/2026, Sec. 5.1, 5.6, 6.1) gives DEG    % MMP, 11/17/2026
% and the include and sep it implies, which OPTIONS may then not set.       % MMP, 11/17/2026
% Fields other than len_R are 1 x 3 cells, entry e for a space in e         % MMP, 11/17/2026
% variables, the last entry for every higher e; lengths are lengths of      % MMP, 11/17/2026
% monomial vectors, 0 = absent (a pattern with a zero length is dropped):   % MMP, 11/17/2026
%   len_R      scalar, the multiplier length on R^q (default 1: identity);  % MMP, 11/17/2026
%   len_inc    1 x e, d1: incomplete moments, pattern c (lower integral);   % MMP, 11/17/2026
%   len_com    1 x e, d2: complete moments, pattern b;                      % MMP, 11/17/2026
%   len_sos    3 x e, d3: multiplier lengths of the blocks a, b, c;         % MMP, 11/17/2026
%   patterns   'uniform' (a..a, b..b), 'separable' ({a,b}^e), 'all'         % MMP, 11/17/2026
%              ({a,b,c}^e), a cellstr of words ('aab'), or rows (1 a,       % MMP, 11/17/2026
%              2 b, 3 c);                                                   % MMP, 11/17/2026
%   cap        (optional) 1 x e, Inf = none: theta_k+s_k <= cap_k in an     % MMP, 11/17/2026
%              owned direction k; in a direction k the space does not       % MMP, 11/17/2026
%              own, theta_k <= the first cap of the largest owner;          % MMP, 11/17/2026
%   len_cross  (optional) scalar, the monomial length in a direction the    % MMP, 11/17/2026
%              space does not own; absent: the largest int of the basis     % MMP, 11/17/2026
%              operator in its owned directions;                            % MMP, 11/17/2026
%   terms, offsets  'none', []: the storage carries no term (an error       % MMP, 11/17/2026
%              otherwise).                                                  % MMP, 11/17/2026
% Codes: a -> 1; c -> 2; b -> 3 in a direction where some space carries     % MMP, 11/17/2026
% c, else 4 (sep true there). In an owned direction, mult = 0, d2-1,        % MMP, 11/17/2026
% d1-1 and int = d3(a,:)-1, d3(b,:)-1, d3(c,:)-1 for a, b, c. Without       % MMP, 11/17/2026
% len_com the fields mult, int, joint, sep of lpi.pos are read: per         % MMP, 11/17/2026
% space, entry min(e,3): the multiplier basis of degree mult, integral      % MMP, 11/17/2026
% bases of degree int and total joint; a direction owned by a space in      % MMP, 11/17/2026
% two or more variables, or any with sep set, takes the full integral.      % MMP, 11/17/2026
%
% OUTPUT
% - prog:        'struct', the program with Q declared positive
%                semidefinite and Pop appended to prog.sop_pos;
% - Pop:         'cdopvar' object, M x M blocks on X, the operator above;
% - Qcell:       N_b x N_b cell, Qcell{b,b'} the 'cellstr' names of the
%                decision variables of block (b,b') of Q, N_b the number
%                of basis operators of all spaces;
% - basis_list:  double, N_b x (1+nv), row b the basis operator b: its
%                space k, then its multi-index entry alpha_d in each
%                registry direction, 0 where s^k lacks variable d; rows
%                ordered by space, as Qcell.
%
% See also POSCOPVAR_DIRECT, COPQUADVAR, POSSOPVAR, SOPQUADVAR,
%          LPI_EQ_CDOPVAR, LPI_INEQ_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - poscopvar
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
% Initial coding MMP, 09/21/2026
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  mopquadvar -> copquadvar, posmopvar -> poscopvar. File was
%                  'posmopvar.m'.
% MMP, 10/17/2026: Built by 'poscopvar_direct' instead of 'copquadvar': the
%                  same operator and decision variables (test_poscopvar_direct,
%                  24 cases) at 2-4x less time cold, 7-19x warm, and a third
%                  of the 3-D peak memory. The position map of each distinct
%                  shape stays cached for the session (about 16 B per
%                  coefficient nonzero); poscopvar_direct('clear') releases
%                  it, options.cache = false skips it. 'copquadvar' stays the
%                  routine for type 'sym' and the reference of the tests.
%                  options.psatz keeps copquadvar's check (one code, 0, 1 or
%                  a face), since poscopvar_direct accepts more and
%                  test_copquadvar_faces asserts the rejections.
% MMP, 10/21/2026: Pop is recorded in prog.sop_pos, from which lpi_ineq_sop
%                  sizes the certificate of an inequality that depends on
%                  Pop (maintainer: lpi_ineq_sop needs no hand-set options).
% MMP, 10/28/2026: Header to the header protocol, step comments in the body;
%                  code unchanged. Deleted: the wrapper description, the
%                  references to copquadvar and poscopvar_direct for the
%                  inputs and outputs, a test name and the reason
%                  options.type is refused.
% MMP, 10/28/2026: Empty DIMS (no space: an empty channel of a PIE, as the
%                  W of the H2 executives on R^0 when the PIE has no
%                  disturbance or no regulated output) returns the 0 x 0
%                  operator with no decision variable. Before, the space
%                  parser read the empty list as one R^q space and failed
%                  (maintainer: empty channels should not raise an error).
% MMP, 11/04/2026: Help text names LPIPROGRAM as the source of PROG;
%                LPIPROGRAM_SOP is retired.
% MMP, 11/08/2026: Comment pass (maintainer, 11/08/2026): log condensed, See
%                  also without the retired poslpivar; code unchanged.
% MMP, 11/12/2026: DOM may be omitted: the domain of PROG, by name
%                  (maintainer, 11/12/2026: "Settings on the program"; DOM
%                  repeated the domain given to LPIPROGRAM). options.type,
%                  path and cache are no longer documented: no caller, test
%                  included, sets them. 'pos' was the only legal type, and
%                  the check here repeated that of POSCOPVAR_DIRECT, which
%                  also keeps path and cache (defaults 'fast', on).
% MMP, 11/12/2026: Help text: psatz_rows and psatz_rows_offset (POSCOPVAR_DIRECT,
%                  the products of faces in different directions) pass
%                  through in OPTIONS; code unchanged.
% MMP, 11/12/2026: DEG may be omitted: the degrees of the program settings,
%                  lpi.pos, space by space, entry min(d,3) for a space in d
%                  variables (phase 2 of the knob reduction; maintainer, for
%                  two variables: "N-D rule"). SEP being per direction, a
%                  direction owned by a space in two or more variables is
%                  separable for every space. Before, DEG was required.
% MMP, 11/14/2026: Operator form POSCOPVAR(PROG,X,DEG,OPTIONS): DIMS, SPACES
%                  and DOM are the output spaces and the domain of X, read
%                  by COPVAR_SPACE_LIST as in the executives (maintainer,
%                  11/14/2026: the calls in the demos restated what T
%                  carries and were too complicated).
% MMP, 11/17/2026: DEG omitted and lpi.pos carrying len_com: the storage
%                  row of the settings report 11/16/2026 (Sec. 5.1, 5.6,
%                  6.1), read by the local ROW_STORAGE, ported from the
%                  validated scratch reader pf_row_storage (same DEG,
%                  include and sep): patterns, codes, cap, len_cross and
%                  the int of non-owned directions; a term other than
%                  'none' is an error. Without len_com the fields mult,
%                  int, joint, sep of 11/12/2026 are read as before, until
%                  the tier files switch to the row. Help text: the
%                  omitted-DEG note under DEG (11/12/2026) replaced by the
%                  paragraph DEG OMITTED.

% Check the inputs; OPTIONS a struct, of type 'pos' only
% DOM omitted: the 4th input is DEG (a scalar, a 'struct' without field     % MMP, 11/12/2026
% dom, or a 'cell', never n x 2), and the domain is that of PROG by name,   % MMP, 11/12/2026
% per element, since the varname of a polynomial vector is sorted while     % MMP, 11/12/2026
% the rows of prog.dom follow the program order                             % MMP, 11/12/2026
nin = nargin;                                                               % MMP, 11/12/2026
% Operator form POSCOPVAR(PROG,X,DEG,OPTIONS): the output spaces of X and   % MMP, 11/14/2026
% its domain stand for DIMS, SPACES, DOM; DEG and OPTIONS shift by two      % MMP, 11/14/2026
if isa(dims,'copvar') || isa(dims,'cdopvar')                                % MMP, 11/14/2026
    narginchk(2,4);                                                         % MMP, 11/14/2026
    if nin>=4,  options = dom;  end                                         % MMP, 11/14/2026
    if nin>=3,  deg = spaces;   end                                         % MMP, 11/14/2026
    X = dims;   [spaces,dims] = copvar_space_list(X,'out');                 % MMP, 11/14/2026
    dom = struct('vars',{X.vars},'dom',X.dom);  nin = nin+2;                % MMP, 11/14/2026
end                                                                         % MMP, 11/14/2026
% if nin>=4 && ~((isnumeric(dom) && (isempty(dom) || size(dom,2)==2)) ...   % MMP, 11/12/2026 (was)
%                || (isstruct(dom) && isfield(dom,'dom')))                  % MMP, 11/12/2026 (was)
%     if nin>=5,  options = deg;  end                                       % MMP, 11/12/2026 (was)
%     deg = dom;  nin = nin+1;                                              % MMP, 11/12/2026 (was)
% Also DOM and DEG both omitted, POSCOPVAR(PROG,DIMS,SPACES)                % MMP, 11/12/2026
if nin==3 || (nin>=4 && ~((isnumeric(dom) && (isempty(dom) || size(dom,2)==2)) ...
               || (isstruct(dom) && isfield(dom,'dom'))))                   % MMP, 11/12/2026
    if nin>=5,  options = deg;  end                                         % MMP, 11/12/2026
    if nin>=4,  deg = dom;      end                                         % MMP, 11/12/2026
    nin = nin+1;                                                            % MMP, 11/12/2026
    v = cell(1,size(prog.dom,1));                                           % MMP, 11/12/2026
    for k = 1:numel(v),     v{k} = prog.vartable(k).varname{1};     end     % MMP, 11/12/2026
    dom = struct('vars',{v},'dom',prog.dom);                                % MMP, 11/12/2026
end                                                                         % MMP, 11/12/2026
% DEG omitted or empty: the degrees of the program settings (below)         % MMP, 11/12/2026
usedef = nin<5 || isempty(deg);                                             % MMP, 11/12/2026
% % if nargin<5                                                             % MMP, 11/12/2026 (was)
% if nin<5                                                                  % MMP, 11/12/2026 (was)
%     error("Not enough input arguments.")                                  % MMP, 11/12/2026 (was)
% end                                                                       % MMP, 11/12/2026 (was)
% if nargin<6 || isempty(options)                                           % MMP, 11/12/2026 (was)
if nin<6 || isempty(options)                                                % MMP, 11/12/2026
    options = struct();
end
if ~isa(options,'struct')
    error("Options should be specified as a 'struct' object.")
end
% if isfield(options,'type') && ~isempty(options.type) ...                  % MMP, 11/12/2026 (was)
%         && ~strcmp(char(options.type),'pos')                              % MMP, 11/12/2026 (was)
%     error("'poscopvar' declares a positive operator; for type '%s' call "... % MMP, 11/12/2026 (was)
%           +"'copquadvar' directly.",char(options.type))                   % MMP, 11/12/2026 (was)
% end                                                                       % MMP, 11/12/2026 (was)
% options.type = 'pos';                                                     % MMP, 11/12/2026 (was)
% type: POSCOPVAR_DIRECT rejects any but 'pos'                              % MMP, 11/12/2026

% No space (an empty channel of a PIE, as the W of an H2 LPI on R^0): the   % MMP, 10/28/2026
% 0 x 0 operator, with no block and no decision variable; PROG unchanged    % MMP, 10/28/2026
if isempty(dims)                                                            % MMP, 10/28/2026
    meta0 = struct('vars',{cell(1,0)},'dom',zeros(0,2),'space_out',false(0,0),...
                   'space_in',false(0,0),'dim_out',zeros(0,1),'dim_in',zeros(0,1),...
                   'Zd',{cell(0,1)});                                       % MMP, 10/28/2026
    Pop = cdopvar(cell(0,0),meta0);     Qcell = {};     basis_list = zeros(0,1); % MMP, 10/28/2026
    return                                                                  % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026

% [prog,Pop,Qcell,basis_list] = copquadvar(prog,dims,spaces,dom,deg,options); % MMP, 10/17/2026 (was)
meta = parse_copvar_spaces(dims,spaces,dom);    % the registry, for nv      % MMP, 10/17/2026
% DEG omitted: the degrees of prog.settings.lpi.pos (LPIPROGRAM; light for  % MMP, 11/12/2026
% a program without settings) space by space, entry min(d,3) for a space    % MMP, 11/12/2026
% in d variables: on R^n the identity basis; otherwise the multiplier       % MMP, 11/12/2026
% basis of degree mult and integral bases of degree int in each variable    % MMP, 11/12/2026
% and joint in total (int < 0: none), the lower and upper integrals in one  % MMP, 11/12/2026
% variable, one full integral in a separable direction. SEP is per          % MMP, 11/12/2026
% direction: a direction that a space in two or more variables owns is      % MMP, 11/12/2026
% separable for every space, any other when lpi.pos.sep is set.             % MMP, 11/12/2026
if usedef                                                                   % MMP, 11/12/2026
    if any(isfield(options,{'include','sep'}))                              % MMP, 11/12/2026
        error("With DEG omitted the settings give include and sep; OPTIONS may not set them.") % MMP, 11/12/2026
    end                                                                     % MMP, 11/12/2026
    ps = struct('mult',[1 0 0],'int',[1 1 1],'joint',[2 2 2],'sep',0);      % MMP, 11/12/2026
    if isfield(prog,'settings') && isfield(prog.settings,'lpi') && isfield(prog.settings.lpi,'pos') % MMP, 11/12/2026
        ps = prog.settings.lpi.pos;                                         % MMP, 11/12/2026
    end                                                                     % MMP, 11/12/2026
    % lpi.pos with len_com: the storage row (header, DEG OMITTED; settings  % MMP, 11/17/2026
    % report 11/16/2026, Sec. 6.1) by ROW_STORAGE; else the rule above,     % MMP, 11/17/2026
    % kept until the tier files switch to the row                           % MMP, 11/17/2026
    if isfield(ps,'len_com')                                                % MMP, 11/17/2026
        [deg,sepd,incl] = row_storage(ps,meta.space_out);                   % MMP, 11/17/2026
    else                                                                    % MMP, 11/17/2026
    own = meta.space_out;   nd = sum(own,2);    e = @(x,c) x(min(c,end));   % MMP, 11/12/2026
    sepd = logical(ps.sep(1)) | any(own(nd>=2,:),1);                        % MMP, 11/12/2026
    deg = cell(1,numel(nd));    incl = cell(1,numel(nd));                   % MMP, 11/12/2026
    for k = 1:numel(nd)                                                     % MMP, 11/12/2026
        if nd(k)==0,    deg{k} = struct('int',0);   continue;   end         % MMP, 11/12/2026
        c = min(nd(k),3);   m = e(ps.mult,c);   b = e(ps.int,c);            % MMP, 11/12/2026
        Zi = struct('int',b,'mult',b,'joint',e(ps.joint,c));                % MMP, 11/12/2026
        deg{k} = {struct('int',m,'mult',0,'joint',m), Zi, Zi};              % MMP, 11/12/2026
        incl{k} = [1;2;3];                                                  % MMP, 11/12/2026
        if all(sepd(own(k,:)))                                              % MMP, 11/12/2026
            deg{k} = deg{k}(1:2);   incl{k} = [ones(1,nd(k)); 4*ones(1,nd(k))]; % MMP, 11/12/2026
        end                                                                 % MMP, 11/12/2026
        if b<0,     deg{k} = deg{k}(1);     incl{k} = incl{k}(1,:);     end % MMP, 11/12/2026
    end                                                                     % MMP, 11/12/2026
    end                                                                     % MMP, 11/17/2026
    options.sep = sepd;     options.include = incl;                         % MMP, 11/12/2026
end                                                                         % MMP, 11/12/2026
% copquadvar's contract for options.psatz is kept: one code, 0, 1 or a face % MMP, 10/17/2026
% 2d+1 / 2d+2. poscopvar_direct also takes a row of codes and the box       % MMP, 10/17/2026
% quadratics 2nv+2+d (the terms of lpi_ineq_sop); poscopvar rejects them,   % MMP, 10/17/2026
% as before (test_copquadvar_faces part (i)), with copquadvar's message.    % MMP, 10/17/2026
nv = numel(meta.vars);      ps = 0;                                         % MMP, 10/17/2026
if isfield(options,'psatz') && ~isempty(options.psatz)                      % MMP, 10/17/2026
    ps = options.psatz;                                                     % MMP, 10/17/2026
end                                                                         % MMP, 10/17/2026
if ~isscalar(ps) || ~ismember(ps,[0,1,3:2*nv+2])                            % MMP, 10/17/2026
    error("'psatz' should be 0, 1, or a face code 2d+1 / 2d+2 with d "...
          +"in 1..nv, here nv = "+num2str(nv)+"; 2 is "...
          +"poslpivar_2d's ball, not the box.")                             % MMP, 10/17/2026
end                                                                         % MMP, 10/17/2026
% Declare Q and build Pop
[prog,Pop,Qcell,info] = poscopvar_direct(prog,dims,spaces,dom,deg,options); % MMP, 10/17/2026
basis_list = info.basis_list;                   % as copquadvar's           % MMP, 10/17/2026
% The record of lpi_ineq_sop: an inequality that depends on Pop is sized    % MMP, 10/21/2026
% from Pop's support, not from every value of its decision variables.       % MMP, 10/21/2026
% MATLAB copies on write, so the record shares Pop's arrays.                % MMP, 10/21/2026
if ~isfield(prog,'sop_pos'),    prog.sop_pos = {};  end                     % MMP, 10/21/2026
prog.sop_pos{end+1} = Pop;                                                  % MMP, 10/21/2026

end


function [deg,sepv,incl] = row_storage(row,own)                             % MMP, 11/17/2026
% [DEG,SEPV,INCL] = ROW_STORAGE(ROW,OWN) the degrees, sep and include rows  % MMP, 11/17/2026
% of POSCOPVAR_DIRECT for the storage row ROW (header, DEG OMITTED;         % MMP, 11/17/2026
% settings report 11/16/2026, Sec. 5.1, 6.1) on the spaces of OWN.          % MMP, 11/17/2026
% INPUT  row:   'struct', lpi.pos with field len_com;                       % MMP, 11/17/2026
%        own:   M x nv logical, row k the registry variables of space k.    % MMP, 11/17/2026
% OUTPUT deg:   1 x M cell, per space a cell of one degree struct per basis % MMP, 11/17/2026
%               operator (int, mult, cap: 1 x nv), or a struct on R^q;      % MMP, 11/17/2026
%        sepv:  1 x nv logical, true where no owning space carries c;       % MMP, 11/17/2026
%        incl:  1 x M cell, the multi-index rows of each space.             % MMP, 11/17/2026
% Ported from the validated scratch reader pf_row_storage (11/16/2026):     % MMP, 11/17/2026
% the same DEG, INCL and SEPV for the same row. Small arrays only.          % MMP, 11/17/2026
[M,nv] = size(own);     nd = sum(own,2);                                    % MMP, 11/17/2026
e_of = @(k) min(max(nd(k),1),3);    % the entry: the last for e > 3         % MMP, 11/17/2026
if ~isfield(row,'len_R') || isempty(row.len_R),     row.len_R = 1;  end     % MMP, 11/17/2026
% The storage carries no term (settings report, Sec. 6.4)                   % MMP, 11/17/2026
if isfield(row,'terms') && ~isempty(row.terms)                              % MMP, 11/17/2026
    t = row.terms;  if ~iscell(t),  t = {t};    end                         % MMP, 11/17/2026
    if ~all(cellfun(@(x) isempty(x) || isequal(char(x),'none'),t))          % MMP, 11/17/2026
        error("The storage carries no term: lpi.pos.terms should be 'none'.") % MMP, 11/17/2026
    end                                                                     % MMP, 11/17/2026
end                                                                         % MMP, 11/17/2026
% Pattern rows per space (1 a, 2 b, 3 c); a zero length drops the pattern   % MMP, 11/17/2026
pats = cell(1,M);                                                           % MMP, 11/17/2026
for k = 1:M                                                                 % MMP, 11/17/2026
    if nd(k)==0,    continue;   end                                         % MMP, 11/17/2026
    e = e_of(k);    n = nd(k);                                              % MMP, 11/17/2026
    pt = row.patterns{e};                                                   % MMP, 11/17/2026
    if ischar(pt)                                                           % MMP, 11/17/2026
        switch pt                                                           % MMP, 11/17/2026
            case 'uniform',     A = [ones(1,n); 2*ones(1,n)];               % MMP, 11/17/2026
            case 'separable',   A = multiindex_grid(repmat({[1 2]},1,n));   % MMP, 11/17/2026
            case 'all',         A = multiindex_grid(repmat({[1 2 3]},1,n)); % MMP, 11/17/2026
            otherwise,  error("Unknown pattern set '"+pt+"' in lpi.pos.patterns.") % MMP, 11/17/2026
        end                                                                 % MMP, 11/17/2026
    elseif iscell(pt)   % words over a, b, c                                % MMP, 11/17/2026
        A = zeros(numel(pt),n);                                             % MMP, 11/17/2026
        for i = 1:numel(pt),    w = pt{i};    A(i,:) = (w=='b')*1 + (w=='c')*2 + 1;   end % MMP, 11/17/2026
    else                                                                    % MMP, 11/17/2026
        A = pt;                                                             % MMP, 11/17/2026
    end                                                                     % MMP, 11/17/2026
    d1 = fit_len(row.len_inc{e},n);   d2 = fit_len(row.len_com{e},n);   d3 = row.len_sos{e}; % MMP, 11/17/2026
    if size(d3,2)~=n,   d3 = repmat(d3(:,1),1,n);    end                    % MMP, 11/17/2026
    keep = true(size(A,1),1);                                               % MMP, 11/17/2026
    for i = 1:size(A,1)                                                     % MMP, 11/17/2026
        for d = 1:n                                                         % MMP, 11/17/2026
            a = A(i,d);                                                     % MMP, 11/17/2026
            if (a==2 && d2(d)==0) || (a==3 && d1(d)==0) || d3(a,d)==0,   keep(i) = false;  end % MMP, 11/17/2026
        end                                                                 % MMP, 11/17/2026
    end                                                                     % MMP, 11/17/2026
    pats{k} = A(keep,:);                                                    % MMP, 11/17/2026
end                                                                         % MMP, 11/17/2026
% Separable directions: no owning space carries c there (b -> 4, the full   % MMP, 11/17/2026
% integral, spans b and c with equal multiplier lengths; Sec. 6.1)          % MMP, 11/17/2026
sepv = true(1,nv);                                                          % MMP, 11/17/2026
for k = 1:M                                                                 % MMP, 11/17/2026
    if nd(k)==0,    continue;   end                                         % MMP, 11/17/2026
    idx = find(own(k,:));                                                   % MMP, 11/17/2026
    for d = 1:nd(k)                                                         % MMP, 11/17/2026
        if any(pats{k}(:,d)==3),    sepv(idx(d)) = false;    end            % MMP, 11/17/2026
    end                                                                     % MMP, 11/17/2026
end                                                                         % MMP, 11/17/2026
% The number of variables of the largest space owning each direction        % MMP, 11/17/2026
big = zeros(1,nv);                                                          % MMP, 11/17/2026
for d = 1:nv,   o = find(own(:,d));  if ~isempty(o),  big(d) = max(nd(o));  end,  end % MMP, 11/17/2026
% Degrees and include rows                                                  % MMP, 11/17/2026
deg = cell(1,M);    incl = cell(1,M);                                       % MMP, 11/17/2026
for k = 1:M                                                                 % MMP, 11/17/2026
    if nd(k)==0,    deg{k} = struct('int',row.len_R-1);    continue;   end  % MMP, 11/17/2026
    e = e_of(k);    n = nd(k);  idx = find(own(k,:));                       % MMP, 11/17/2026
    d1 = fit_len(row.len_inc{e},n);   d2 = fit_len(row.len_com{e},n);   d3 = row.len_sos{e}; % MMP, 11/17/2026
    if size(d3,2)~=n,   d3 = repmat(d3(:,1),1,n);    end                    % MMP, 11/17/2026
    A = pats{k};    rows = zeros(size(A));   dk = cell(1,size(A,1));        % MMP, 11/17/2026
    for i = 1:size(A,1)                                                     % MMP, 11/17/2026
        mu = zeros(1,nv);   in = zeros(1,nv);                               % MMP, 11/17/2026
        for d = 1:n                                                         % MMP, 11/17/2026
            a = A(i,d);     g = idx(d);                                     % MMP, 11/17/2026
            switch a    % b: 4 (full) where separable, else 3 (upper)       % MMP, 11/17/2026
                case 1,     rows(i,d) = 1;                      mu(g) = 0;  % MMP, 11/17/2026
                case 2,     rows(i,d) = 4 - ~sepv(g);           mu(g) = d2(d)-1; % MMP, 11/17/2026
                case 3,     rows(i,d) = 2;                      mu(g) = d1(d)-1; % MMP, 11/17/2026
            end                                                             % MMP, 11/17/2026
            in(g) = d3(a,d)-1;                                              % MMP, 11/17/2026
        end                                                                 % MMP, 11/17/2026
        % A direction the space does not own: len_cross - 1 when given,     % MMP, 11/17/2026
        % else the operator's largest owned int (Sec. 6.1)                  % MMP, 11/17/2026
        for g = find(~own(k,:))                                             % MMP, 11/17/2026
            in(g) = max(in(idx));                                           % MMP, 11/17/2026
            if isfield(row,'len_cross') && ~isempty(row.len_cross) && ~isempty(row.len_cross{e}) % MMP, 11/17/2026
                in(g) = row.len_cross{e} - 1;                               % MMP, 11/17/2026
            end                                                             % MMP, 11/17/2026
        end                                                                 % MMP, 11/17/2026
        dk{i} = struct('int',in,'mult',mu);                                 % MMP, 11/17/2026
        % Per-direction cap (Sec. 5.6): owned directions from this entry,   % MMP, 11/17/2026
        % a non-owned one from the entry of its largest owner               % MMP, 11/17/2026
        if isfield(row,'cap') && ~isempty(row.cap)                          % MMP, 11/17/2026
            cv = Inf(1,nv);                                                 % MMP, 11/17/2026
            ce = fit_len(row.cap{e},n);                                     % MMP, 11/17/2026
            for d = 1:n,    cv(idx(d)) = ce(d);     end                     % MMP, 11/17/2026
            for g = find(~own(k,:))                                         % MMP, 11/17/2026
                cb = row.cap{min(big(g),3)};    cv(g) = cb(1);              % MMP, 11/17/2026
            end                                                             % MMP, 11/17/2026
            dk{i}.cap = cv;                                                 % MMP, 11/17/2026
        end                                                                 % MMP, 11/17/2026
    end                                                                     % MMP, 11/17/2026
    deg{k} = dk;    incl{k} = rows;                                         % MMP, 11/17/2026
end                                                                         % MMP, 11/17/2026
end                                                                         % MMP, 11/17/2026


function v = fit_len(v,n)                                                   % MMP, 11/17/2026
% V = FIT_LEN(V,N) a row of N lengths from a scalar or a row of N           % MMP, 11/17/2026
v = reshape(v,1,[]);                                                        % MMP, 11/17/2026
if isscalar(v),     v = repmat(v,1,n);                                      % MMP, 11/17/2026
elseif numel(v)~=n                                                          % MMP, 11/17/2026
    error("A length row of lpi.pos has "+numel(v)+" entries for a space in "+n+" variables.") % MMP, 11/17/2026
end                                                                         % MMP, 11/17/2026
end                                                                         % MMP, 11/17/2026
