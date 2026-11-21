function prog = lpiprogram(varargin)
% PROG = LPIPROGRAM(VARTAB,DUMVARTAB,DOM,DECVARTAB,FREEVARTAB,SETTINGS)
% declares an LPI program structure in independent variables 'vartab' on
% domain 'dom', with dummy variables 'dumvartab', decision variables
% 'decvartab', additional (non-spatial) independent variables
% 'freevartab', and settings 'settings'.
%
% INPUT
% - vartab:     nx1 array of type 'polynomial', specifying n independent
%               variables for the optimization program, corresponding to
%               spatial variables in the LPI;
% - dumvartab:  (optional) nx1 array of type 'polynomial', specifying dummy
%               variables associated with each of the specified spatial
%               variables. If not specified, or declared as [], default 
%               dummy variables will be generated, so that if vartab(i)=si, 
%               then  dumvartab(i) = si_dum;
% - dom:        nx2 array of type 'double', with each row dom(i,:)=[a,b]
%               specifying the interval on which the independent variable
%               vartab(i) (and corresponding dummy variable) is defined.
% - decvartab:  (optional) qx1 array of type 'dpvar', specifying
%               decision variables to be used in the optimization
%               program.
% - freevartab: (optional) mx1 array of type 'polynomial', specifying
%               additional independent variables in the optimization 
%               program that need not exist on any fixed domain.
% - settings:   (optional, the last input) 'struct' of LPI settings, or a
%               'char' or 'string' naming a tier: 'stripped', 'light',
%               'heavy', 'veryheavy', 'custom' or 'extreme'; default
%               'light'. Stored in prog.settings, from which the LPI
%               routines take the options a call does not give.
%
% OUTPUT
% - prog:       'struct' specifying an optimization program structure as
%               per SOSTOOLS format, but with additional fields 'dom',
%               specifying the domain of the independent variables, and
%               'settings', the settings structure.
%
% NOTES
% The independent variables will be stored in the lpi program structure
% 'prog' under 'prog.vartable', ordered as
%       prog.vartable = [vartab(:); dumvartab(:); freevartab(:)];
% Thus, if size(prog.dom,1)=n, then prog.vartable(1:n) will correspond to
% primary spatial variables, prog.vartable(n+1:2*n) to dummy variables, and
% prog.vartable(2*n+1:end) to non-spatial variables.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpiprogram
%
% Copyright (C)2024  PIETOOLS Team
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
% Initial coding DJ - 10/13/2024
% DJ, 11/30/2024: Introduce separate input for dummy variables;
% DJ, 01/23/2025: Allow dummy variables to be specified as empty array;
% MMP, DJ, 11/21/2026: Changes since PIETOOLS 2025: optional last input
%   SETTINGS (a settings structure or a tier name, default 'light'), stored
%   in prog.settings, the other inputs declared by the new local function
%   declare_program; any number of spatial variables, also none; variable
%   names as cellstr accepted; each dummy variable is checked to be a
%   polynomial variable; help and comments revised.

% % % SETTINGS: a 'struct' or text as the last input (no other input is     % MMP, DJ, 11/21/2026
% either), stored in PROG; the inputs before it declare the program.        % MMP, DJ, 11/21/2026
st = [];                                                                    % MMP, DJ, 11/21/2026
if nargin>=2,   st = varargin{end};     end                                 % MMP, DJ, 11/21/2026
if isstruct(st) || ischar(st) || isstring(st),  varargin(end) = [];         % MMP, DJ, 11/21/2026
else,                                           st = [];                    % MMP, DJ, 11/21/2026
end                                                                         % MMP, DJ, 11/21/2026
prog = lpisetsettings(declare_program(varargin{:}),st);                     % MMP, DJ, 11/21/2026
end                                                                         % MMP, DJ, 11/21/2026


function prog = declare_program(vartab,dumvartab,dom,decvartab,freevartab)  % MMP, DJ, 11/21/2026
% PROG = DECLARE_PROGRAM(VARTAB,DUMVARTAB,DOM,DECVARTAB,FREEVARTAB)         % MMP, DJ, 11/21/2026
% declares the program of LPIPROGRAM from the inputs before SETTINGS.       % MMP, DJ, 11/21/2026
% INPUT, OUTPUT: as LPIPROGRAM, without SETTINGS and prog.settings.         % MMP, DJ, 11/21/2026

% % % Process the inputs
% % First, check if the second argument corresponds to dummy variables or a
% % domain.
if nargin==1
    error("No domain has been specified for the spatial variables.")
elseif nargin>=3 && isnumeric(dom)
    % Third argument is numeric
    % --> second argument must correspond to dummy variables.
elseif nargin>=2 && nargin<=4 && isnumeric(dumvartab)
    % Second argument is numeric, likely specifies the domain.
    if nargin==2 || ~isnumeric(dom)
        % Assume second argument corresponds to domain.
        if size(vartab,1)==size(dumvartab,1) && size(vartab,2)==2
            % Second column of vartab corresponds to dummy vars.
            dumvartab_tmp = vartab(:,2);
            vartab = vartab(:,1);
        else
            dumvartab_tmp = [];
        end
        % Call the function with correct order of input arguments.
        if nargin==4
            freevartab = decvartab;
            decvartab = dom;
            dom = dumvartab;
%           prog = lpiprogram(vartab,dumvartab_tmp,dom,decvartab,freevartab); % MMP, DJ, 11/21/2026 (was)
            prog = declare_program(vartab,dumvartab_tmp,dom,decvartab,freevartab); % MMP, DJ, 11/21/2026
        elseif nargin==3
            decvartab = dom;
            dom = dumvartab;
%           prog = lpiprogram(vartab,dumvartab_tmp,dom,decvartab);          % MMP, DJ, 11/21/2026 (was)
            prog = declare_program(vartab,dumvartab_tmp,dom,decvartab);     % MMP, DJ, 11/21/2026
        else
            dom = dumvartab;
%           prog = lpiprogram(vartab,dumvartab_tmp,dom);                    % MMP, DJ, 11/21/2026 (was)
            prog = declare_program(vartab,dumvartab_tmp,dom);               % MMP, DJ, 11/21/2026
        end
        return
    end
end

% % Check that the spatial variables are properly specified.
if isempty(vartab)
    vartab = polynomial(zeros(0,1));
elseif iscellstr(vartab)                                                    % MMP, DJ, 11/21/2026
    vartab = polynomial(vartab);
elseif ~isa(vartab,'polynomial')
    error("Spatial variables in the LPI optimization program should be specified as nx1 array of type 'polynomial'.")
elseif ~ispvar(vartab)
    error("Each element of the array of spatial variables should correspond to a single polynomial variable.")
end
% if prod(size(vartab))~=max(size(vartab))                                  % MMP, DJ, 11/21/2026 (was)
if ~isempty(vartab) && prod(size(vartab))~=max(size(vartab))  % 0 x 1 is a column % MMP, DJ, 11/21/2026
    error("Spatial variables in the LPI optimization program should be specified as nx1 array.")
end
vartab = vartab(:);
% if size(vartab,1)>2                                                       % MMP, DJ, 11/21/2026 (was)
%     error('LPI programs involving more than 2 spatial variables are currently not supported.') % MMP, DJ, 11/21/2026 (was)
% end                                                                       % MMP, DJ, 11/21/2026 (was)

% % Check that the dummy varaibles are properly specified.
if isempty(dumvartab)                                                       % DJ, 01/23/2025
    % No dummy variables have been specified
    % --> generate default dummy variables.
    if isempty(vartab)
        dumvartab = polynomial(zeros(0,1));
    else
        % For each variable si, generate dummy variable si_dum.
        dum_vars = cell(size(vartab));
        for ii=1:size(vartab,1)
            dum_vars{ii} = [vartab(ii).varname{1},'_dum'];
        end
        dumvartab = polynomial(dum_vars);
    end
else
    if ~isa(dumvartab,'polynomial')
        error("Dummy variables in the LPI optimization program should be specified as nx1 array of type 'polynomial'.")
%   elseif ~ispvar(vartab)                                                  % MMP, DJ, 11/21/2026 (was)
    elseif ~ispvar(dumvartab)                                               % MMP, DJ, 11/21/2026
        error("Each element of the array of dummy variables should correspond to a single polynomial variable.")
    end
    dumvartab = dumvartab(:);
    if length(dumvartab)~=length(vartab)
        error('Number of dummy variables should match the number of primary spatial variables.')
    end
end

% % Check that the spatial domain is properly specified.
if nargin<=2
    error("No domain has been specified for the spatial variables.")
elseif isempty(dom) && isempty(vartab)
%   dom = zeros(0,1);                                                       % MMP, DJ, 11/21/2026 (was)
    dom = zeros(0,2);       % n x 2 for n = 0, so dom(:,2) below is valid   % MMP, DJ, 11/21/2026
elseif ~isa(dom,'double') && ~(isa(dom,'polynomial') && isdouble(dom))
    error("Spatial domain should be specified as nx2 array of type 'double'.")
elseif size(dom,2)~=2
    error("Spatial domain should be specified as nx2 array.")
elseif size(dom,1)~=size(vartab,1)
    % The number of domains should match the number of spatial variables,
    % but we have to be careful that the user may specify [s1,s2] as two
    % spatial variables (on different domains), or [s1,th2] as spatial and
    % dummy variable (on same domain).
    if size(dom,1)==2 && all(size(vartab)==[1,2])
        % Assume the user declares two spatial variables, [s1,s2];
        vartab = vartab';
    elseif size(dom,1)==1
        % Assume all rows of vartable specify different spatial variables, 
        % but set same domain for all spatial variables.
        dom = repmat(dom,[size(vartab,1),1]);
    else
        error("Number of rows in domain array should match number of spatial variables.")
    end
end
dom = double(dom);
% Make sure the domain makes sense.
if any(dom(:,2)<=dom(:,1))
    error('Upper boundary of domain should be strictly greater than lower boundary.')
end



% % Check that the "free" variables and decision variables are 
% % properly specified.
if nargin<=3
    freevartab = polynomial(zeros(0,1));
    decvartab = dpvar(zeros(0,1));
elseif nargin==4
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
elseif nargin==5
    % Allow the user to switch order in which decision variables and
    % independent variables are specified.
    if isa(freevartab,'dpvar') && isa(decvartab,'polynomial')
        new_decvartable = freevartab;
        freevartab = decvartab;
        decvartab = new_decvartable;
    end
    % Check the independent variables.
    if ~isa(freevartab,'polynomial')
        error("Additional independent variables in the LPI optimization program should be specified as mx1 array of type 'polynomial'.")
    elseif ~ispvar(freevartab)
        error("Each element of the array of additional independent variables should correspond to a single polynomial variable.")
    else
        freevartab = freevartab(:);
    end
    % Check the decision variables.
    if ~isa(decvartab,'dpvar')
        error("Decision variables in the optimization program should be specified as qx1 array of type 'dpvar'.")
    else
        decvartab = decvartab(:);
    end
else
    error("Too many input arguments.")
end



% % % Declare the optimization program;

% The LPI program structure is mostly just an SOSTOOLS program structure.
% However, 'sosprogram' reorders the input independent variables in
% alphabetical order. Since the order is important, we declare the
% independent variables manually.
prog = sosprogram(polynomial([]),decvartab);
prog.vartable = [prog.vartable; vartab(:); dumvartab(:); freevartab(:)];

% Also set the domain of the spatial variables, adding a field to the
% optimization program structure. 
prog.dom = dom;

end