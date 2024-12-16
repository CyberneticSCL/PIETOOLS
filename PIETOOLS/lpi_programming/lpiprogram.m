function prog = lpiprogram(vartable,dom,freevartable,decvartable)
% PROG = LPIPROGRAM(VARTABLE,DOM,DECVARTABLE) declares an LPI program
% structure in independent variables 'vartable' on domain 'dom', with
% additional independent variables 'freevartable' (non-spatial variables)
% and decision variables 'decvartable'.
%
% INPUT
% - vartable:   nx2 array of type 'polynomial', specifying n independent
%               variables for the optimization program in the first column,
%               along with associated dummy variables used for integration
%               in the second column. If an nx1 array of spatial variables
%               is provided, associated dummy variables will be generated,
%               so that if vartable(i,1)=s1, then vartable(1,2) = s1_dum;
% - dom:        nx2 array of type 'double', with each row dom(i,:)=[1,b]
%               specifying the interval on which the independent variable
%               vartable(i,1) (and corresponding dummy variable) is
%               defined.
% - freevartable:   (optional) mx1 array of type 'polynomial', specifying
%                   additional independent variables in the optimization 
%                   program that need not exist on any fixed domain.
% - decvartable:    (optional) qx1 array of type 'dpvar', specifying
%                   decision variables to be used in the optimization
%                   program.
%
% OUTPUT
% - prog:       'struct' specifying an optimization program structure as
%               per SOSTOOLS format, but with an additional field 'dom'
%               specifying the domain of the independent variables.
%
% NOTES
% The independent variables will be stored in the lpi program structure
% 'prog' under 'prog.vartable', ordered as
%       prog.vartable = [vartable(:,1); vartable(:,2); freevartable];
% Thus, if size(prog.dom,1)=n, then prog.vartable(1:n) will correspond to
% primary spatial variables, prog.vartable(n+1:2*n) to dummy variables, and
% prog.vartable(2*n+1:end) to non-spatial variables.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpiprogram
%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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


% % % Process the inputs
% % Check that independent variables are properly specified
if isempty(vartable)
    vartable = polynomial(zeros(0,2));
elseif ~isa(vartable,'polynomial')
    error("Spatial variables in the LPI optimization program should be specified as nx2 array of type 'polynomial'.")
elseif ~ispvar(vartable)
    error("Each element of the array of spatial variables should correspond to a single polynomial variable.")
elseif size(vartable,2)~=2
    if size(vartable,1)==1
        % Assume a row vector of only primary spatial variables is
        % specified.
        vartable = vartable';
    elseif size(vartable,2)~=1
        error("Spatial variables should be specified as nx2 array, with first column specifying primary variables, and second column dummy variables used for integration.")
    end
end
if size(vartable,1)>2
    error('LPI programs involving more than 2 spatial variables are currently not supported.')
end

% % Check that the spatial domain is properly specified.
if nargin<2
    error("No domain has been specified for the spatial variables.")
end
if isempty(dom) && isempty(vartable)
    dom = zeros(0,1); 
elseif ~isa(dom,'double') && ~(isa(dom,'polynomial') && isdouble(dom))
    error("Spatial domain should be specified as nx2 array of type 'double'.")
elseif size(dom,2)~=2
    error("Spatial domain should be specified as nx2 array.")
elseif size(dom,1)~=size(vartable,1)
    % The number of domains should match the number of spatial variables,
    % but we have to be careful that the user may specify [s1,s2] as two
    % spatial variables (on different domains), or [s1,th2] as spatial and
    % dummy variable (on same domain).
    if size(dom,1)==2 && all(size(vartable)==[1,2])
        % Assume the user declares two spatial variables, [s1,s2];
        vartable = vartable';
    elseif size(dom,1)==1
        % Assume all rows of vartable specify different spatial variables, 
        % but set same domain for all spatial variables.
        dom = repmat(dom,[size(vartable,1),1]);
    else
        error("Number of rows in domain array should match number of spatial variables.")
    end
end
dom = double(dom);
% Make sure the domain makes sense.
if any(dom(:,2)<=dom(:,1))
    error('Upper boundary of domain should be strictly greater than lower boundary.')
end

% % If only spatial variables have been declared, add dummy variables:
if size(vartable,2)==1
    dum_vars = cell(size(vartable,1),1);
    for ii=1:size(vartable,1)
        dum_vars{ii} = [vartable(ii,1).varname{1},'_dum'];
    end
    vartable = [vartable, polynomial(dum_vars)];
end

% % Check that the "free" variables and decision variables are 
% % properly specified.
if nargin<=2
    freevartable = polynomial(zeros(0,1));
    decvartable = dpvar(zeros(0,1));
elseif nargin==3
    if isa(freevartable,'dpvar')
        decvartable = freevartable;
        freevartable = polynomial([]);
    elseif ~isa(freevartable,'polynomial')
        error("Additional independent variables in the LPI optimization program should be specified as mx1 array of type 'polynomial'.")
    elseif ~ispvar(vartable)
        error("Each element of the array of additional independent variables should correspond to a single polynomial variable.")
    else
        freevartable = freevartable(:);
        decvartable = dpvar([]);
    end
elseif nargin==4
    % Allow the user to switch order in which decision variables and
    % independent variables are specified.
    if isa(freevartable,'dpvar') && isa(decvartable,'polynomial')
        new_decvartable = freevartable;
        freevartable = decvartable;
        decvartable = new_decvartable;
    end
    % Check the independent variables.
    if ~isa(freevartable,'polynomial')
        error("Additional independent variables in the LPI optimization program should be specified as mx1 array of type 'polynomial'.")
    elseif ~ispvar(freevartable)
        error("Each element of the array of additional independent variables should correspond to a single polynomial variable.")
    else
        freevartable = freevartable(:);
    end
    % Check the decision variables.
    if ~isa(decvartable,'dpvar')
        error("Decision variables in the optimization program should be specified as qx1 array of type 'dpvar'.")
    else
        decvartable = decvartable(:);
    end
else
    error("Too many input arguments.")
end



% % % Declare the optimization program;

% The LPI program structure is mostly just an SOSTOOLS program structure.
% However, 'sosprogram' reorders the input independent variables in
% alphabetical order. Since the order is important, we declare the
% independent variables manually.
prog = sosprogram(polynomial([]),decvartable);
prog.vartable = [prog.vartable; vartable(:,1); vartable(:,2); freevartable];

% Also set the domain of the spatial variables, adding a field to the
% optimization program structure. 
prog.dom = dom;

end