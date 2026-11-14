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
% - deg:         (optional; omitted or [], the degrees of the program
%                settings, lpi.pos of LPIPROGRAM, space by space, with the
%                include and sep they imply, which OPTIONS may then not set)
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
