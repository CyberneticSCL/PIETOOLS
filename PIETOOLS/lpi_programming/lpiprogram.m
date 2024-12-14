function prog = lpiprogram(vartab,dumvartab,dom,decvartab,freevartab)
% PROG = LPIPROGRAM(VARTAB,DUMVARTAB,DOM,DECVARTAB,FREEVARTAB) declares an
% LPI program structure in independent variables 'vartab' on domain 'dom',
% with dummy variables 'dumvartab', decision variables 'decvartab', and
% additional (non-spatial) independent variables 'freevartab'.
%
% INPUT
% - vartab:     nx1 array of type 'polynomial', specifying n independent
%               variables for the optimization program, corresponding to
%               spatial variables in the LPI;
% - dumvartab:  (optional) nx1 array of type 'polynomial', specifying dummy
%               variables associated with each of the specified spatial
%               variables. If not specified, default dummy variables will
%               be generated, so that if vartab(i)=si, then 
%               dumvartab(i) = si_dum;
% - dom:        nx2 array of type 'double', with each row dom(i,:)=[1,b]
%               specifying the interval on which the independent variable
%               vartab(i) (and corresponding dummy variable) is defined.
% - decvartab:  (optional) qx1 array of type 'dpvar', specifying
%               decision variables to be used in the optimization
%               program.
% - freevartab: (optional) mx1 array of type 'polynomial', specifying
%               additional independent variables in the optimization 
%               program that need not exist on any fixed domain.
%
% OUTPUT
% - prog:       'struct' specifying an optimization program structure as
%               per SOSTOOLS format, but with an additional field 'dom'
%               specifying the domain of the independent variables.
%
% NOTES
% The independent variables will be stored in the lpi program structure
% 'prog' under 'prog.vartable', ordered as
%       prog.vartable = [vartab(:); dumvartab(:); freevartab(:)];
% Thus, if size(prog.dom,1)=n, then prog.vartable(1:n) will correspond to
% primary spatial variables, prog.vartable(n+1:2*n) to dummy variables, and
% prog.vartable(2*n+1:end) to non-spatial variables.

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
% DJ, 11/30/2024: Introduce separate input for dummy variables.


% % % Process the inputs
% % First, check if the second argument corresponds to dummy variables or a
% % domain.
if nargin==1
    error("No domain has been specified for the spatial variables.")
elseif nargin>=2 && nargin<=4 && isnumeric(dumvartab)
    if size(vartab,1)==size(dumvartab,1) && size(vartab,2)==2
        % Second column of vartab corresponds to dummy vars.
        dumvartab_tmp = vartab(:,2);
        vartab = vartab(:,1);
    elseif prod(size(vartab))==max(size(vartab))
        % No dummy variables have been specified
        % --> generate default dummy variables.
        vartab = vartab(:);
        dum_vars = cell(length(vartab),1);
        for ii=1:size(vartab,1)
            dum_vars{ii} = [vartab(ii).varname{1},'_dum'];
        end
        dumvartab_tmp = polynomial(dum_vars);
    else
        error("Spatial variables in the LPI optimization program should be specified as nx1 array.")
    end
    % Call the function with default dummy variables.
    if nargin==4
        freevartab = decvartab;
        decvartab = dom;
        dom = dumvartab;
        prog = lpiprogram(vartab,dumvartab_tmp,dom,decvartab,freevartab);
    elseif nargin==3
        decvartab = dom;
        dom = dumvartab;
        prog = lpiprogram(vartab,dumvartab_tmp,dom,decvartab);
    else
        dom = dumvartab;
        prog = lpiprogram(vartab,dumvartab_tmp,dom);
    end
    return
end

% % Check that the spatial variables are properly specified.
if isempty(vartab)
    vartab = polynomial(zeros(0,2));
elseif ~isa(vartab,'polynomial')
    error("Spatial variables in the LPI optimization program should be specified as nx1 array of type 'polynomial'.")
elseif ~ispvar(vartab)
    error("Each element of the array of spatial variables should correspond to a single polynomial variable.")
end
if prod(size(vartab))~=max(size(vartab))
    error("Spatial variables in the LPI optimization program should be specified as nx1 array.")
end
vartab = vartab(:);
if size(vartab,1)>2
    error('LPI programs involving more than 2 spatial variables are currently not supported.')
end

% % Check that the dummy varaibles are properly specified.
if ~isa(dumvartab,'polynomial')
    error("Dummy variables in the LPI optimization program should be specified as nx1 array of type 'polynomial'.")
elseif ~ispvar(vartab)
    error("Each element of the array of dummy variables should correspond to a single polynomial variable.")
end
dumvartab = dumvartab(:);
if length(dumvartab)~=length(vartab)
    error('Number of dummy variables should match the number of primary spatial variables.')
end

% % Check that the spatial domain is properly specified.
if nargin<=2
    error("No domain has been specified for the spatial variables.")
elseif isempty(dom) && isempty(vartab)
    dom = zeros(0,1); 
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