function prog = lpiprogram(vartable,dom,freevartable,decvartable)
% PROG = LPIPROGRAM(VARTABLE,DOM,DECVARTABLE) declares an LPI program
% structure in independent variables 'vartable' on domain 'dom', with
% decision variables 'decvartable'.
%
% INPUT
% - vartable:   nx2 array of type 'polynomial', specifying n independent
%               variables for the optimization program in the first column,
%               along with associated dummy variables used for integration
%               in the second column.
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



% % % Process the inputs
% % Check that independent variables are properly specified
if ~isa(vartable,'polynomial')
    error("Spatial variables in the LPI optimization program should be specified as nx2 array of type 'polynomial'.")
elseif ~ispvar(vartable)
    error("Each element of the array of spatial variables should correspond to a single polynomial variable.")
elseif size(vartable,2)~=2
    error("Spatial variables should be specified as nx2 array, with first column specifying primary variables, and second column dummy variables used for integration.")
end

% % Check that the spatial domain is properly specified.
if nargin<2
    error("No domain has been specified for the spatial variables.")
end
if ~isa(dom,'double') && ~(isa(dom,'polynomial') && isdouble(dom))
    error("Spatial domain should be specified as nx2 array of type 'double'.")
elseif size(dom,2)~=2
    error("Spatial domain should be specified as nx2 array.")
elseif size(dom,1)~=size(vartable,1)
    if size(dom,1)==1
        % Assume same domain for all spatial variables.
        dom = repmat(dom,[size(vartable,1),1]);
    else
        error("Number of rows in domain array should match number of spatial variables.")
    end
end

% % Check that the "free" variables and decision variables are 
% % properly specified.
if nargin<=2
    freevartable = polynomial([]);
    decvartable = dpvar([]);
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
    elseif ~ispvar(vartable)
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
vartab_1 = vartable(:,1);
prog = sosprogram(vartab_1,decvartable);
prog.vartable = [prog.vartable; vartable(:,2); freevartable];
% Also set the domain of the spatial variables, adding a field to the
% optimization program structure. 
prog.dom = dom;


end