function prog = lpiprogram_sop(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = LPIPROGRAM_SOP(VARTAB,DUMVARTAB,DOM,DECVARTAB,FREEVARTAB) declares
% an LPI program structure in any number of spatial variables. It is
% 'lpiprogram' without the cap of two spatial variables, and with two more
% ways to name the variables, so that one program serves the legacy
% operator classes and the containers of the 'sopvar' family:
%
%   PROG = LPIPROGRAM_SOP(VARTAB,DUMVARTAB,DOM[,DECVARTAB[,FREEVARTAB]])
%   PROG = LPIPROGRAM_SOP(VARTAB,DOM[,DECVARTAB[,FREEVARTAB]])
%          every form 'lpiprogram' accepts, with the same result;
%   PROG = LPIPROGRAM_SOP(NAMES,DOM,...), NAMES a cellstr (the registry
%          P.vars of a container, with DOM = P.dom);
%   PROG = LPIPROGRAM_SOP(P[,DECVARTAB[,FREEVARTAB]]), P a 'copvar',
%          'cdopvar', 'sopvar' or 'sdopvar': the variables and domains of
%          P, in the order of its sorted registry;
%   PROG = LPIPROGRAM_SOP(S,...), S a struct with fields 'vars' (cellstr)
%          and 'dom' (the struct form of DOM in 'lpivar_cdopvar').
%
% INPUT
% - vartab:     n x 1 'polynomial' array of the spatial variables, or a
%               cellstr of their names;
% - dumvartab:  (optional) n x 1 'polynomial' array (or cellstr) of dummy
%               variables, one per spatial variable; [] gives s_dum for
%               variable s, as 'lpiprogram' does;
% - dom:        n x 2 'double', row i the interval of vartab(i); one 1 x 2
%               row applies to every variable;
% - decvartab:  (optional) q x 1 'dpvar' of decision variables;
% - freevartab: (optional) m x 1 'polynomial' of non-spatial variables.
%
% OUTPUT
% - prog:       SOSTOOLS program structure with the field 'dom' added, the
%               fields and values 'lpiprogram' builds: vartable =
%               [vartab; dumvartab; freevartab] ('polynomial'), dom n x 2,
%               decvartable a cellstr. For n <= 2 it equals the output of
%               'lpiprogram' field by field (tests/test_lpiprogram_sop).
%
% NOTES
% Nothing in 'lpiprogram' depends on the number of variables except the
% cap, so the program for n > 2 is the one 'lpiprogram' would build. The
% container routines (lpivar_cdopvar, poscopvar, lpi_eq_cdopvar,
% lpi_eq_sdopvar, getsol_lpivar_sop) read neither prog.vartable nor
% prog.dom; 'lpidecvar', 'lpisetobj' and 'lpisolve' read neither either.
% The legacy 'lpi_eq' and 'lpi_ineq' read the dummy variables at
% prog.vartable(n+1:2n), which is where they are here, as in 'lpiprogram'.
% The legacy 1-D and 2-D operator routines (lpivar, poslpivar) still
% support at most 2 variables; this lifts only the program's cap.
%
% Cost: O(n) in the number of spatial variables; with DECVARTAB, that of
% 'sosprogram' on it, as in 'lpiprogram'.
%
% See also LPIPROGRAM, LPIDECVAR, LPI_EQ_SOP, LPIGETSOL_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpiprogram_sop
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
% Initial coding MMP, 09/29/2026. 'lpiprogram' (DJ, 10/13/2024; 11/30/2024;
%                01/23/2025) refuses more than 2 spatial variables, and the
%                N-D container scripts (heatNd_lpi, heatNd_poincare,
%                test_copquadvar_faces) rebuilt its output by hand for
%                N > 2. This follows lpiprogram's parsing and construction
%                in full, without the cap, so one code path serves every N
%                and is checked against lpiprogram where lpiprogram runs.
%                Deviations, all on inputs lpiprogram accepts wrongly or
%                not at all: the dummy variables' ispvar check tests the
%                dummy variables (lpiprogram tests vartab again); names
%                may be given as a cellstr or a container; no spatial
%                variable gives dom = zeros(0,2), where lpiprogram errors.

a = varargin;
if isempty(a)
    error("No domain has been specified for the spatial variables.")
end
% % % New input forms: a container, a block, or a struct of vars and dom.
% They carry their own domains, so the forms become (names,[],dom,...).
if is_sop_class(a{1}) || (isstruct(a{1}) && isfield(a{1},'vars') && isfield(a{1},'dom'))
    [nm,dm] = registry_of(a{1});
    a = [{nm,[],dm}, a(2:end)];
end
if numel(a)>5
    error("Too many input arguments.")
end
% Names given as cellstr: the registry of a container ('copvar.vars').
if iscellstr(a{1}) || isstring(a{1}),   a{1} = names2pvar(a{1});    end %#ok<ISCLSTR>

% % % The forms of 'lpiprogram': if the third argument is
% numeric, the second is the dummy variables; else a numeric second
% argument is the domain, and an n x 2 VARTAB carries the dummies.
nin = numel(a);
if nin==1
    error("No domain has been specified for the spatial variables.")
elseif nin>=3 && isnumeric(a{3})
    % (vartab, dumvartab, dom, ...)
elseif nin<=4 && isnumeric(a{2})
    vt = a{1};  dm = a{2};  dvt = [];
    if size(vt,1)==size(dm,1) && size(vt,2)==2
        dvt = vt(:,2);  vt = vt(:,1);           % second column: dummies
    end
    a = [{vt,dvt,dm}, a(3:end)];
    nin = numel(a);
end
vartab = a{1};  dumvartab = a{2};
if nin<3
    error("No domain has been specified for the spatial variables.")
end
dom = a{3};
if iscellstr(dumvartab) || isstring(dumvartab),  dumvartab = names2pvar(dumvartab);  end %#ok<ISCLSTR>

% % % Spatial variables (as lpiprogram, cap removed).
if isempty(vartab)
    vartab = polynomial(zeros(0,1));
elseif ~isa(vartab,'polynomial')
    error("Spatial variables in the LPI optimization program should be specified as nx1 array of type 'polynomial'.")
elseif ~ispvar(vartab)
    error("Each element of the array of spatial variables should correspond to a single polynomial variable.")
end
% lpiprogram's vector test; skipped for no variables, which it would
% refuse (prod of 0x1 is 0) although its own code handles them after.
if ~isempty(vartab) && prod(size(vartab))~=max(size(vartab)) %#ok<PSIZE>
    error("Spatial variables in the LPI optimization program should be specified as nx1 array.")
end
vartab = vartab(:);
n = size(vartab,1);

% % % Dummy variables (as lpiprogram).
if isempty(dumvartab)
    if n==0
        dumvartab = polynomial(zeros(0,1));
    else
        dum_vars = cell(n,1);
        for ii = 1:n
            dum_vars{ii} = [vartab(ii).varname{1},'_dum'];
        end
        dumvartab = polynomial(dum_vars);
    end
else
    if ~isa(dumvartab,'polynomial')
        error("Dummy variables in the LPI optimization program should be specified as nx1 array of type 'polynomial'.")
    elseif ~ispvar(dumvartab)
        error("Each element of the array of dummy variables should correspond to a single polynomial variable.")
    end
    dumvartab = dumvartab(:);
    if length(dumvartab)~=n
        error('Number of dummy variables should match the number of primary spatial variables.')
    end
end

% % % Domain (as lpiprogram). Its [s1,s2]-row branch is not
% reachable there (vartab is a column by then) and is left out here.
if isempty(dom) && n==0
    % No spatial variable (a container on R^n only). lpiprogram sets
    % zeros(0,1), which then fails its own dom(:,2) test below.
    dom = zeros(0,2);
elseif ~isa(dom,'double') && ~(isa(dom,'polynomial') && isdouble(dom))
    error("Spatial domain should be specified as nx2 array of type 'double'.")
elseif size(dom,2)~=2
    error("Spatial domain should be specified as nx2 array.")
elseif size(dom,1)~=n
    if size(dom,1)==1
        dom = repmat(dom,[n,1]);                % one interval for every variable
    else
        error("Number of rows in domain array should match number of spatial variables.")
    end
end
dom = double(dom);
if any(dom(:,2)<=dom(:,1))
    error('Upper boundary of domain should be strictly greater than lower boundary.')
end

% % % Decision and free variables (as lpiprogram).
if nin<=3
    freevartab = polynomial(zeros(0,1));
    decvartab = dpvar(zeros(0,1));
elseif nin==4
    decvartab = a{4};
    if isa(decvartab,'polynomial') && ispvar(decvartab)
        freevartab = decvartab(:);
        decvartab = dpvar([]);
    elseif isa(decvartab,'polynomial')
        error("Each element of the array of additional independent variables should correspond to a single polynomial variable.")
    elseif ~isa(decvartab,'dpvar')
        error("Decision variables in the LPI optimization program should be specified as mx1 array of type 'dpvar'.")
    else
        decvartab = decvartab(:);
        freevartab = polynomial(zeros(0,1));
    end
else
    decvartab = a{4};   freevartab = a{5};
    if isa(freevartab,'dpvar') && isa(decvartab,'polynomial')
        [decvartab,freevartab] = deal(freevartab,decvartab);    % either order
    end
    if ~isa(freevartab,'polynomial')
        error("Additional independent variables in the LPI optimization program should be specified as mx1 array of type 'polynomial'.")
    elseif ~ispvar(freevartab)
        error("Each element of the array of additional independent variables should correspond to a single polynomial variable.")
    end
    freevartab = freevartab(:);
    if ~isa(decvartab,'dpvar')
        error("Decision variables in the optimization program should be specified as qx1 array of type 'dpvar'.")
    end
    decvartab = decvartab(:);
end

% % % The program, as lpiprogram builds it: sosprogram
% sorts independent variables, so they are set by hand in the order
% [primary; dummy; free], which lpi_eq/lpi_ineq index by position.
prog = sosprogram(polynomial([]),decvartab);
prog.vartable = [prog.vartable; vartab(:); dumvartab(:); freevartab(:)];
prog.dom = dom;

end


% ========================================================================
function tf = is_sop_class(x)
tf = isa(x,'copvar') || isa(x,'cdopvar') || isa(x,'sopvar') || isa(x,'sdopvar');
end


% ========================================================================
function [nm,dm] = registry_of(X)
% Names (1 x nv cellstr, sorted) and domains (nv x 2) of the variables of
% X. A container stores both as its registry (vars, dom); a block stores
% them per side (vars.out/in, dom.out/in), which agree on a shared
% variable (sopvar class invariant); a struct gives them directly.
if isa(X,'copvar') || isa(X,'cdopvar')
    nm = X.vars;    dm = X.dom;
elseif isstruct(X)
    nm = cellstr(X.vars);    dm = X.dom;
else
    vo = cellstr(X.vars.out);   vi = cellstr(X.vars.in);
    [nm,io] = unique([vo(:); vi(:)]);
    D = [reshape(X.dom.out,[],2); reshape(X.dom.in,[],2)];
    dm = D(io,:);
end
nm = reshape(cellstr(nm),[],1);
if isempty(nm),     dm = zeros(0,2);    end
end


% ========================================================================
function p = names2pvar(names)
% Column 'polynomial' of the named variables, in the given order.
names = cellstr(names);
if isempty(names),  p = polynomial(zeros(0,1));     return,     end
p = polynomial(names(:));
end
