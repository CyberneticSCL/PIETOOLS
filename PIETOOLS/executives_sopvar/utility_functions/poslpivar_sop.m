function varargout = poslpivar_sop(prog,n,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declares a positive semidefinite PI operator variable, as 'poslpivar'
% does, for either operator family:
%
%   [PROG,POP,QMAT,ZOP,GS] = POSLPIVAR_SOP(PROG,N[,D[,OPTIONS]])
%       N numeric (an opvar's .dim): the legacy routine, poslpivar(...),
%       called unchanged (its 2-D dispatch included);
%   [PROG,P] = POSLPIVAR_SOP(PROG,X,D,OPTIONS[,SIDE[,DOM]])
%       X a 'copvar'/'cdopvar': a psd 'cdopvar' from 'poscopvar' over the
%       SIDE spaces of X, with poslpivar's 1-D degrees D and OPTIONS
%       translated.
%
% INPUT
% - prog:    LPI program structure;
% - n / X:   numeric dimension (legacy), or container whose spaces P maps;
% - d:       poslpivar degrees, {d1,[d2],[d3]} or any form poslpivar fills
%            in ([] is its default {1,[1,1,1],[1,1,1]});
% - options: struct with fields psatz, exclude, sep (each optional), as
%            poslpivar reads them; [] for the defaults;
% - side:    'out' or 'in' (container only); omitted, X must be square and
%            'out' is used;
% - dom:     (container only) the domain as 'poscopvar' takes it; omitted,
%            X's own (vars, dom).
%
% OUTPUT
% - prog:    the program with the Gram matrix declared;
% - P:       'cdopvar' on the listed spaces (container), or poslpivar's
%            outputs (legacy).
%
% NOTES
% The container branch mirrors poslpivar's 1-D branch (lpi_programming/
% poslpivar.m: option defaults; degree gap filling; sep -> exclude(4); the
% psatz weight; the Z1/Z2/Z3 monomial roles):
% - degrees: d is filled exactly as poslpivar fills it, with its errors.
%   An R^q space gets the identity basis (degree 0); an L2 space gets
%   Z1 = {int d1, mult 0}, Z2 = {int d2(1), mult d2(2), joint d2(3)}, Z3
%   likewise from d3. 'int' caps the integration variable (poslpivar's
%   slot 1 after its var2 -> var1 substitution), 'mult' the other.
% - psatz: weight only for psatz == 1, as poslpivar (any other value: none).
% - exclude(2:4) and sep: the kept L2 basis operators become poscopvar's
%   'include' (alpha 1, 2, 3; under sep Z2 is the full integral, alpha 4,
%   and Z3 is dropped, as poslpivar sets exclude(4)).
% Not expressible, so errors: exclude(1) with an R^q space present
% (poscopvar needs one basis operator per space; poslpivar drops the R^q
% term), and exclude(2:4) all set with an L2 space present (poslpivar
% then declares the R^q part only). 1-D only: a space over two or more
% variables is an error (the {d1,[d2],[d3]} vocabulary names three basis
% operators per space).
% Cost: the translation is O(number of spaces); the Gram and its decision
% variables are poscopvar's.
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

if ~(isa(n,'copvar') || isa(n,'cdopvar'))
    % Legacy family: poslpivar itself.
    [varargout{1:max(nargout,1)}] = poslpivar(prog,n,varargin{:});
    return
end
X = n;
if nargout>2
    error('poslpivar_sop:nargout',...
        'A container positive variable has two outputs, [prog,P].')
end
d = [];     options = [];   side = '';  dom = [];
if numel(varargin)>=1,  d = varargin{1};        end
if numel(varargin)>=2,  options = varargin{2};  end
if numel(varargin)>=3,  side = varargin{3};     end
if numel(varargin)>=4,  dom = varargin{4};      end
if isempty(side)
    if ~isequal(X.space_out,X.space_in) || ~isequal(X.dim_out(:),X.dim_in(:))
        error('poslpivar_sop:side',...
            "X is not square: name the spaces with SIDE = 'out' or 'in'.")
    end
    side = 'out';
end
if isempty(dom),    dom = struct('vars',{X.vars},'dom',X.dom);  end

% Options, as poslpivar.m reads them (defaults; a non-struct is an error).
if isempty(options),    options = struct();     end
if ~isstruct(options)
    error("Options for positive operator should be specified as struct with fields 'psatz', 'exclude', and 'sep'.")
end
if ~isfield(options,'psatz'),   options.psatz = 0;              end
if ~isfield(options,'exclude'), options.exclude = [0 0 0 0];    end
if ~isfield(options,'sep'),     options.sep = 0;                end
exc = logical(options.exclude(:).');
sep = options.sep==1;
if sep,     exc(4) = true;  end                     % Z3 = Z2 when R1 = R2

[sp,dm] = copvar_space_list(X,side);
if any(cellfun(@numel,sp)>1)
    error('poslpivar_sop:dim',['poslpivar degrees describe one spatial '...
          'variable per space; a space over %d variables has no such form.'],...
          max(cellfun(@numel,sp)))
end
d = fill_degrees(d);
L2 = { struct('int',d{1},'mult',0), ...                         % Z1
       struct('int',d{2}(1),'mult',d{2}(2),'joint',d{2}(3)), ... % Z2
       struct('int',d{3}(1),'mult',d{3}(2),'joint',d{3}(3)) };   % Z3
alph = [1;2;3];     if sep,     alph = [1;4;0];     end         % alpha 4: full integral
keep = ~exc(2:4);
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
popt = struct('psatz',double(options.psatz==1),'sep',sep,'include',{incl});
[prog,P] = poscopvar(prog,dm,sp,dom,deg,popt);
varargout = {prog,P};
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = fill_degrees(d)
% poslpivar.m's degree default, numeric forms and gap filling, with its
% error messages.
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
