function varargout = poslpivar_sop(prog,n,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P] = POSLPIVAR_SOP(PROG,X,D,OPTIONS[,SIDE[,DOM]]) declares a
% positive semidefinite, self-adjoint 'cdopvar' decision operator P on the
% SIDE spaces of the container X, each space R^q or L2^m[s] with s in
% [a,b]:
%   P = Z'*(g*Q)*Z,   Q >= 0 one Gram matrix of new decision variables,
% where Z = blkdiag(Z_1,...,Z_M) stacks the basis operators of the M
% spaces: on R^q the identity, on L2^m[s] those of
%   (Z_k x)(eta) = [Z1(eta)*x(eta); int_a^eta Z2(eta,s)*x(s) ds;
%                   int_eta^b Z3(eta,s)*x(s) ds]
% that OPTIONS.exclude keeps, each Zi a vector of monomials of the degrees
% D times I_m, and g(eta) = (eta-a)*(b-eta) when OPTIONS.psatz = 1, g = 1
% otherwise. With OPTIONS.sep = 1, Z2 integrates over [a,b] and Z3 is
% dropped.
% [PROG,POP,QMAT,ZOP,GS] = POSLPIVAR_SOP(PROG,N[,D[,OPTIONS]]), with N not
% a 'copvar' or 'cdopvar' object, returns POSLPIVAR(PROG,N,D,OPTIONS).
%
% INPUT
% - prog:     'struct', the LPI program to which the decision variables are
%             added;
% - X:        'copvar' or 'cdopvar' object whose SIDE spaces P maps; or
%   N:        'double', the dimensions [n0;n1] of an operator on
%             R^n0 x L2^n1, or the 2x2 .dim of an 'opvar' object;
% - D:        (optional) 'cell' {d1,d2,d3}: d1 the degree of Z1, d2 =
%             [degree in eta, degree in s, joint degree] of Z2, d3 the same
%             for Z3. Gaps are filled: a one-entry d2 or d3 p becomes
%             [p p p], a two-entry [p q] becomes [p q max(p,q)]; a missing
%             d3 is d2; {d1} is {d1,[d1 d1 d1],[d1 d1 d1]}. A 'double' p is
%             {p}, [p q] is {max(p,q),[p q max(p,q)],[p q max(p,q)]},
%             [p q r] is {max(p,q,r),[p q r],[p q r]}. [] or omitted:
%             {1,[1,1,1],[1,1,1]}. With N, D and OPTIONS are passed to
%             POSLPIVAR unchanged;
% - OPTIONS:  (optional) 'struct' with fields, each optional, psatz (g is
%             used only for psatz = 1; default 0), exclude (1x4 logical:
%             entry 1 the R^q basis, entries 2:4 Z1, Z2, Z3; default
%             [0 0 0 0]) and sep (1 or 0; default 0); [] or omitted for the
%             defaults;
% - SIDE:     (optional) 'char', 'out' (the output spaces of X) or 'in' (its
%             input spaces); omitted or '', X must have the same spaces and
%             dimensions in and out, and 'out' is used;
% - DOM:      (optional) the domain: 'double' 1x2 [a b] for every variable
%             or nv x 2 over the sorted variable names, or 'struct' with
%             fields vars ('cellstr') and dom (nv x 2 'double'); default the
%             variables and domain of X.
%
% OUTPUT
% - prog:     'struct', the program with the decision variables of Q and the
%             constraint Q >= 0 added;
% - P:        'cdopvar' object on the SIDE spaces of X.
% With N: prog; POP, 'dopvar' object, POP = ZOP'*(GS*QMAT)*ZOP; QMAT,
% 'dpvar' object, the Gram matrix; ZOP, 'opvar' object, the monomial
% operator; GS, 'polynomial' object, the weight g. For a 2-D program
% (PROG.dom 2x2), POP is a 'dopvar2d' and ZOP an 'opvar2d' object.
%
% NOTES
% Errors: a space of X over two or more variables; exclude(1) set with an
% R^q space present; every L2 basis operator excluded (exclude(2:4), Z3
% counted as excluded under sep) with an L2 space present; more than two
% outputs with a container X; a malformed D; OPTIONS not a 'struct'.
%
% See also POSLPIVAR, POSCOPVAR, POSLPIVAR_SETTINGS_SOP, COPVAR_SPACE_LIST.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - poslpivar_sop
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
% Initial coding MMP, 10/06/2026. One library translation of a 1-D
%                poslpivar call for containers, replacing the test-folder
%                cx_stability_pos, cx_h2_pos, cx_hinf_posdeg and cx_pl2pm
%                (cx_exec, MMP 09/25/2026), which differed on psatz other
%                than 0/1 (error / pass-through / none), on sep (error in
%                cx_h2_pos) and on degree filling (none in cx_hinf_posdeg).
%                This follows poslpivar on all three. On the shipped
%                settings presets every 1-D container program is unchanged.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body; code
%                unchanged. The deleted header text was the comparison with
%                the 1-D branch of poslpivar and a cost note.

% A dimension N: the call is POSLPIVAR's
if ~(isa(n,'copvar') || isa(n,'cdopvar'))
    % Legacy family: poslpivar itself.
    [varargout{1:max(nargout,1)}] = poslpivar(prog,n,varargin{:});
    return
end
% A container X: at most two outputs
X = n;
if nargout>2
    error('poslpivar_sop:nargout',...
        'A container positive variable has two outputs, [prog,P].')
end
% The optional inputs D, OPTIONS, SIDE, DOM
d = [];     options = [];   side = '';  dom = [];
if numel(varargin)>=1,  d = varargin{1};        end
if numel(varargin)>=2,  options = varargin{2};  end
if numel(varargin)>=3,  side = varargin{3};     end
if numel(varargin)>=4,  dom = varargin{4};      end
% SIDE omitted: X must be square, and its output spaces are used
if isempty(side)
    if ~isequal(X.space_out,X.space_in) || ~isequal(X.dim_out(:),X.dim_in(:))
        error('poslpivar_sop:side',...
            "X is not square: name the spaces with SIDE = 'out' or 'in'.")
    end
    side = 'out';
end
% DOM omitted: the variables and domain of X
if isempty(dom),    dom = struct('vars',{X.vars},'dom',X.dom);  end

% Options, as poslpivar.m reads them (defaults; a non-struct is an error).
if isempty(options),    options = struct();     end
if ~isstruct(options)
    error("Options for positive operator should be specified as struct with fields 'psatz', 'exclude', and 'sep'.")
end
if ~isfield(options,'psatz'),   options.psatz = 0;              end
if ~isfield(options,'exclude'), options.exclude = [0 0 0 0];    end
if ~isfield(options,'sep'),     options.sep = 0;                end
% The excluded basis operators; sep excludes Z3
exc = logical(options.exclude(:).');
sep = options.sep==1;
if sep,     exc(4) = true;  end                     % Z3 = Z2 when R1 = R2

% The SIDE spaces of X (variable names, dimensions), each in at most one
% variable
[sp,dm] = copvar_space_list(X,side);
if any(cellfun(@numel,sp)>1)
    error('poslpivar_sop:dim',['poslpivar degrees describe one spatial '...
          'variable per space; a space over %d variables has no such form.'],...
          max(cellfun(@numel,sp)))
end
% Degrees of Z1, Z2, Z3 ('int' caps eta, 'mult' caps s), their
% multi-indices (1 multiplier, 2 lower, 3 upper, 4 full integral), and the
% kept ones
d = fill_degrees(d);
L2 = { struct('int',d{1},'mult',0), ...                         % Z1
       struct('int',d{2}(1),'mult',d{2}(2),'joint',d{2}(3)), ... % Z2
       struct('int',d{3}(1),'mult',d{3}(2),'joint',d{3}(3)) };   % Z3
alph = [1;2;3];     if sep,     alph = [1;4;0];     end         % alpha 4: full integral
keep = ~exc(2:4);
% Per space: the identity basis on R^q, the kept Z1, Z2, Z3 on L2[s]
deg = cell(1,numel(sp));    incl = cell(1,numel(sp));
for k = 1:numel(sp)
    if isempty(sp{k})                               % R^q: the identity basis
        if exc(1)
            error('poslpivar_sop:excludeRn',['exclude(1) drops the R^q basis; '...
                  'poscopvar needs one basis operator per space.'])
        end
        deg{k} = struct('int',0);   incl{k} = [];
    else
        if ~any(keep)
            error('poslpivar_sop:excludeL2','exclude removes every L2 basis operator.')
        end
        deg{k} = L2(keep);          incl{k} = alph(keep);
    end
end
% The positive operator on those spaces, by POSCOPVAR
popt = struct('psatz',double(options.psatz==1),'sep',sep,'include',{incl});
[prog,P] = poscopvar(prog,dm,sp,dom,deg,popt);
varargout = {prog,P};
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = fill_degrees(d)
% D = FILL_DEGREES(D) returns the degrees D as a 1x3 cell {d1,d2,d3} with
% d2, d3 1x3, from the forms the header lists, with the gaps filled.
% INPUT  d: [], 'double' (scalar, 1x2 or 1x3) or 'cell' of one to three
%           entries; any other form is an error.
% OUTPUT d: 'cell' {d1,d2,d3}.
% Default; a 'double' to a cell
if isempty(d),  d = {1,[1,1,1],[1,1,1]};    end
if isnumeric(d)
    d = d(:)';
    if isscalar(d),         d = {d};
    elseif numel(d)==2,     d = {max(d),[d,max(d)],[d,max(d)]};
    elseif numel(d)==3,     d = {max(d),d,d};
    else,                   error('Degrees must be specified as 1x3 cell array.')
    end
elseif ~iscell(d)
    error('Degrees must be specified as 1x3 cell array.')
end
% Gaps: {d1} to three entries; d2, d3 to 1x3; a missing d3 from d2
if numel(d)==1
    d{2} = [d{1},d{1},d{1}];    d{3} = d{2};
else
    if numel(d{2})==1,      d{2} = d{2}*ones(1,3);
    elseif numel(d{2})==2,  d{2}(3) = max(d{2});
    end
    if numel(d)==2
        d{3} = d{2};
    elseif numel(d)==3
        if numel(d{3})==1,      d{3} = d{3}*ones(1,3);
        elseif numel(d{3})==2,  d{3}(3) = max(d{3});
        end
    else
        error('Degrees must be specified as 1x3 cell array.')
    end
end
end
