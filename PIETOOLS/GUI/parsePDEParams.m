function PDE_struct = parsePDEParams(varlist, PDE_Text, metadata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parsePDEParams.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function converts the textual PDE of the PDE GUI
% to a matlab PDE object in terms format.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following inputs must be defined passed:
%
%  varlist - variables table in the PDE gui: cell array, one row per
%            variable, its name (x, X, w, u, z or y and a number) in
%            column 1 and, for an X state, its order of differentiability
%            in column 2
% PDE_text - PDE text to be parsed from the GUI: cell array of char, one
%            equation per entry
% metadata - row numbers metadata from PDE gui: 'double' 1x5, the entries
%            of PDE_text that head the PDE, ODE, regulated output, observed
%            output and boundary condition equations
%
% OUTPUT
% PDE_struct - initialized 'pde_struct' object of the PDE (terms format)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding MMP, SS, DJ (PIETOOLS 2021b; in the repository since 12/31/2021)
% MMP, 11/05/2026: The loop over differentiability orders runs to the
%   largest order, max(Xvals(:,1)), instead of max(Xvals): Xvals has two
%   columns, so with two or more PDE states max(Xvals) is a row and R2025b
%   rejects it as a colon operand (33 of the 47 GUI files of the example
%   library failed to parse). Earlier releases used its first element,
%   max(Xvals(:,1)), so the result there is unchanged.
% MMP, 11/08/2026: Each PDE state gets the global variable and domain
%   (PDE_struct.vars, PDE_struct.dom), the values initialize assigns to a
%   state without them. Since c47e093a (11/07/2026) initialize reads the
%   vars of the states in boundary terms before assigning these defaults,
%   so without them 46 of the 47 GUI files of the example library failed.
% MMP, 11/19/2026: Release comment pass (comment-only).
% MMP, 11/20/2026: eval_b: an empty coefficient text gives 1 (index error); a
%   text that fails to evaluate is an error naming it, not 0 (see eval_b).

% remove whitespaces
% PDE_Text = regexprep(PDE_Text,'\s+','');
PDE_text = erase(PDE_Text,'$'); % remove $
PDE_text = strrep(PDE_text, '\theta', 'theta'); %replace latex compatible \theta with pvar theta string

% build a lookup table for indices
xlist = varlist(contains(varlist(:,1),'x'),1);
Xlist = varlist(contains(varlist(:,1),'X'),1);
wlist = varlist(contains(varlist(:,1),'w'),1);
ulist = varlist(contains(varlist(:,1),'u'),1);
zlist = varlist(contains(varlist(:,1),'z'),1);
ylist = varlist(contains(varlist(:,1),'y'),1);

keys = [xlist;wlist;ulist;zlist;ylist];
vals = [1:length(xlist) 1:length(wlist) 1:length(ulist) 1:length(zlist) 1:length(ylist)]; 
if ~isempty(keys)
    PDE.varidx = containers.Map(keys, vals);
else
    PDE.varidx = [];
end
% For each X state: its order of differentiability, and its index among
% the X states of that order
Xvals = cell2mat(varlist(contains(varlist(:,1),'X'),2));
Xvals = [Xvals zeros(size(Xvals))];
% for i=0:max(Xvals)                                                        % MMP, 11/05/2026 (was)
for i=0:max(Xvals(:,1))     % the orders are column 1                       % MMP, 11/05/2026
    Xvals(Xvals(:,1)==i,2) = (1:length(Xvals(Xvals(:,1)==i,2)))';
end

if ~isempty(Xlist)
    PDE.Xvaridx = containers.Map(Xlist,num2cell(Xvals,2));
else
    PDE.Xvaridx = [];
end

% Index of each variable among the variables of its kind, the x and X
% states numbered together (x first)
keys = [xlist;Xlist;wlist;ulist;zlist;ylist];
vals = [1:length(xlist) length(xlist)+1:length(xlist)+length(Xlist) 1:length(wlist) 1:length(ulist) 1:length(zlist) 1:length(ylist)]; 
statenameMap = containers.Map(keys,vals);

% Number of variables of each kind; n_pde(i+1) is the number of X states
% differentiable i times
pvar s theta;
PDE.vars = [s;theta];
PDE.n.nx = sum(contains(varlist(:,1),'x'));
PDE.n.nw = sum(contains(varlist(:,1),'w'));
PDE.n.nu = sum(contains(varlist(:,1),'u'));
np_list = varlist(contains(varlist(:,1),'X'),:); np_list = cell2mat(np_list(:,2));
n_pde = [];
for i=0:max(np_list)
    n_pde(i+1) = sum(np_list==i);
end
PDE.n.n_pde = n_pde;
PDE.n.nz = sum(contains(varlist(:,1),'z'));
PDE.n.ny = sum(contains(varlist(:,1),'y'));


% first set the standard known parameters of pde_struct
PDE_struct.dim = 1;
PDE_struct.dom = [0 1];
PDE_struct.vars = PDE.vars(1);

% extract dimensions of various signals, x X z y w u & BC
no = PDE.n.nx; nw = PDE.n.nw; nu = PDE.n.nu;
n_pde = PDE.n.n_pde; nz = PDE.n.nz; ny = PDE.n.ny;
N = length(n_pde)-1; %max derivative
nBC = sum(n_pde.*(0:N));

% now initialize cell structure for the dynamics, outputs and BCs
PDE_struct.x = cell(no+sum(n_pde),1);
PDE_struct.w = cell(nw,1);
PDE_struct.u = cell(nu,1);
PDE_struct.z = cell(nz,1);
PDE_struct.y = cell(ny,1);
PDE_struct.BC = cell(nBC,1);

% specify first no cells as ODE states
for i=1:no
    PDE_struct.x{i}.vars = [];
end
% specify order of differentiation for PDE states
diffvals = varlist(contains(varlist(:,1),'X'),2);
for i=no+1:no+sum(n_pde)
    PDE_struct.x{i}.diff = diffvals{i-no};
    PDE_struct.x{i}.vars = PDE_struct.vars;     % defaults of initialize,   % MMP, 11/08/2026
    PDE_struct.x{i}.dom = PDE_struct.dom;       % set before it reads them  % MMP, 11/08/2026
end


% get PDE params
iLoc = no+1;
for i=metadata(1)+1:metadata(2)-1
    eqn = PDE_text{i}; % get line
    Lstate = regexp(eqn,'[\w]+ =','match'); % get Lstate
    % separate LHS and RHS, also remove partial derivative in time
    eqn = regexprep(eqn,['\\partial_t ',Lstate],'');
    Lstate = Lstate{1}(1:end-2); % remove equal sign
    PDE_struct.x{iLoc}.term = [];
    
    % separate terms on RHS
    eqn = split(eqn);
    ode_coeffs = eqn(contains(eqn,'x')); % ode terms 
    dis_coeffs = eqn(contains(eqn,'w')); % disturbance terms
    inp_coeffs = eqn(contains(eqn,'u')); % input terms
    pde_coeffs = eqn(endsWith(eqn,'(s)')); %pde multiplier terms
    pde_bcoeffs = eqn(endsWith(eqn,'(0)')|endsWith(eqn,'(1)')); %pde boundary terms
    pde_icoeffs = eqn(endsWith(eqn,'dtheta')); %pde integral terms
    
    % place terms in appropriate structure
    for k=1:length(ode_coeffs)
        term =ode_coeffs{k};
        Rstate = regexp(term,'x[0-9]+','match');
        coeff = regexprep(term,'x[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.x = statenameMap(keyval);
        PDE_struct.x{iLoc}.term{k} = tmp;
        clear tmp;
    end
    % Disturbance terms w<i>: coefficient and input index
    for k=1:length(dis_coeffs)
        term =dis_coeffs{k};
        Rstate = regexp(term,'w[0-9]+','match');
        coeff = regexprep(term,'w[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.w = statenameMap(keyval);
        PDE_struct.x{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Control-input terms u<i>: coefficient and input index
    for k=1:length(inp_coeffs)
        term =inp_coeffs{k};
        Rstate = regexp(term,'u[0-9]+','match');
        coeff = regexprep(term,'u[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.u = statenameMap(keyval);
        PDE_struct.x{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % PDE-state terms X<i>(s): coefficient, state index, derivative order D
    for k=1:length(pde_coeffs)
        term =pde_coeffs{k};
        Rstate = regexp(term,'X[0-9]+\(s\)','match');
        coeff = regexprep(term,'X[0-9]+\(s\)','');
        % Order of the s-derivative: 0, 1, or k for \partial_s^k
        der_order = regexp(coeff,'(\\partial_s)+(\^.)?','match');
        if isempty(der_order)
            der_order = 0;
        elseif strcmp(der_order{1},'\partial_s')
            der_order = 1;
        else
            der_order = str2double(regexprep(der_order{1},'\\partial_s\^',''));
        end
        coeff = regexprep(coeff,'(\\partial_s)+(\^.)?','');
        tmp.C = eval_b(coeff);
        keyval = regexprep(Rstate{1},'\(s\)','');
        tmp.x = statenameMap(keyval);
        tmp.D = der_order;
        PDE_struct.x{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Boundary values X<i>(0), X<i>(1): coefficient, location loc, order D
    for k=1:length(pde_bcoeffs)
        term =pde_bcoeffs{k};
        Rstate = regexp(term,'X[0-9]+\(.\)','match');
        coeff = regexprep(term,'X[0-9]+\(.\)','');
        % Order of the s-derivative: 0, 1, or k for \partial_s^k
        der_order = regexp(coeff,'(\\partial_s)+(\^.)?','match');
        if isempty(der_order)
            der_order = 0;
        elseif strcmp(der_order{1},'\partial_s')
            der_order = 1;
        else
            der_order = str2double(regexprep(der_order{1},'\\partial_s\^',''));
        end
        coeff = regexprep(coeff,'(\\partial_s)+(\^.)?','');
        tmp.C = eval_b(coeff);
        tmp.loc = regexprep(Rstate{1},'X[0-9]+\(',''); tmp.loc = str2double(tmp.loc(1:end-1));
        keyval = regexprep(Rstate{1},'\(.\)','');
        tmp.x = statenameMap(keyval);
        tmp.D = der_order;
        PDE_struct.x{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Integral terms of X<i>: coefficient, limits I{1}, derivative order D
    for k=1:length(pde_icoeffs)
        term =pde_icoeffs{k};
        Rstate = regexp(term,'X[0-9]+','match');
        intlimits = regexp(term,'\\int_.\^.','match');
        b = regexp(intlimits{1},'\^.','match');
        a = regexp(intlimits{1},'_.\^','match');
        coeff = regexprep(term,'(\\int_.\^.)+|(*X[0-9]+(\(theta\)dtheta)','');
        % Order of the s-derivative: 0, 1, or k for \partial_s^k
        der_order = regexp(coeff,'(\\partial_s)+(\^.)?','match');
        if isempty(der_order)
            der_order = 0;
        elseif strcmp(der_order{1},'\partial_s')
            der_order = 1;
        else
            der_order = str2double(regexprep(der_order{1},'\\partial_s\^',''));
        end
        coeff = regexprep(coeff,'(\\partial_s)+(\^.)?','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.x = statenameMap(keyval);
        tmp.D = der_order;
        tmp.I{1} = [eval_b(a{1}(2:end-1)) eval_b(b{1}(2:end))];
        PDE_struct.x{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    iLoc = iLoc+1;
end
% get ODE params
iLoc = 1;
for i=metadata(2)+1:metadata(3)-1
    eqn = PDE_text{i}; % get line
    Lstate = regexp(eqn,'[\w]+ =','match'); % get Lstate
    % separate LHS and RHS, also remove partial derivative in time
    eqn = regexprep(eqn,['\\partial_t ',Lstate],'');
    Lstate = Lstate{1}(1:end-2); % remove equal sign
    PDE_struct.x{iLoc}.term = [];
    
    % separate terms on RHS
    eqn = split(eqn);
    ode_coeffs = eqn(contains(eqn,'x')); % ode terms 
    dis_coeffs = eqn(contains(eqn,'w')); % disturbance terms
    inp_coeffs = eqn(contains(eqn,'u')); % input terms
    pde_bcoeffs = eqn(endsWith(eqn,'(0)')|endsWith(eqn,'(1)')); %pde boundary terms
    pde_icoeffs = eqn(endsWith(eqn,'ds')); %pde integral terms
    
    % place terms in appropriate structure
    for k=1:length(ode_coeffs)
        term =ode_coeffs{k};
        Rstate = regexp(term,'x[0-9]+','match');
        coeff = regexprep(term,'x[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.x = statenameMap(keyval);
        PDE_struct.x{iLoc}.term{k} = tmp;
        clear tmp;
    end
    % Disturbance terms w<i>: coefficient and input index
    for k=1:length(dis_coeffs)
        term =dis_coeffs{k};
        Rstate = regexp(term,'w[0-9]+','match');
        coeff = regexprep(term,'w[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.w = statenameMap(keyval);
        PDE_struct.x{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Control-input terms u<i>: coefficient and input index
    for k=1:length(inp_coeffs)
        term =inp_coeffs{k};
        Rstate = regexp(term,'u[0-9]+','match');
        coeff = regexprep(term,'u[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.u = statenameMap(keyval);
        PDE_struct.x{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Boundary values X<i>(0), X<i>(1): coefficient, location loc, order D
    for k=1:length(pde_bcoeffs)
        term =pde_bcoeffs{k};
        Rstate = regexp(term,'X[0-9]+\(.\)','match');
        coeff = regexprep(term,'X[0-9]+\(.\)','');
        % Order of the s-derivative: 0, 1, or k for \partial_s^k
        der_order = regexp(coeff,'(\\partial_s)+(\^.)?','match');
        if isempty(der_order)
            der_order = 0;
        elseif strcmp(der_order{1},'\partial_s')
            der_order = 1;
        else
            der_order = str2double(regexprep(der_order{1},'\\partial_s\^',''));
        end
        coeff = regexprep(coeff,'(\\partial_s)+(\^.)?','');
        tmp.C = eval_b(coeff);
        tmp.loc = regexprep(Rstate{1},'X[0-9]+\(',''); tmp.loc = str2double(tmp.loc(1:end-1));
        keyval = regexprep(Rstate{1},'\(.\)','');
        tmp.x = statenameMap(keyval);
        tmp.D = der_order;
        PDE_struct.x{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Integral terms of X<i>: coefficient, limits I{1}, derivative order D
    for k=1:length(pde_icoeffs)
        term =pde_icoeffs{k};
        Rstate = regexp(term,'X[0-9]+','match');
        coeff = regexprep(term,'(\\int_.\^.)+|(*X[0-9]+(\(s\)ds)','');
        % Order of the s-derivative: 0, 1, or k for \partial_s^k
        der_order = regexp(coeff,'(\\partial_s)+(\^.)?','match');
        if isempty(der_order)
            der_order = 0;
        elseif strcmp(der_order{1},'\partial_s')
            der_order = 1;
        else
            der_order = str2double(regexprep(der_order{1},'\\partial_s\^',''));
        end
        coeff = regexprep(coeff,'(\\partial_s)+(\^.)?','');
        tmp.C = eval_b(coeff);
        tmp.I{1} = [0,1];
        keyval = Rstate{1};
        tmp.x = statenameMap(keyval);
        tmp.D = der_order;
        PDE_struct.x{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    iLoc = iLoc+1;
end
% get regulated output params
iLoc = 1;
for i=metadata(3)+1:metadata(4)-1
    eqn = PDE_text{i}; % get line
    Loutput = regexp(eqn,'[\w]+ =','match'); % get Loutput
    
    % separate LHS and RHS, also remove partial derivative in time
    eqn = regexprep(eqn,Loutput,'');
    Lstate = Loutput{1}(1:end-2); % remove equal sign
    PDE_struct.z{iLoc}.term = [];
    % separate terms on RHS
    eqn = split(eqn);
    ode_coeffs = eqn(contains(eqn,'x')); % ode terms 
    dis_coeffs = eqn(contains(eqn,'w')); % disturbance terms
    inp_coeffs = eqn(contains(eqn,'u')); % input terms
    pde_bcoeffs = eqn(endsWith(eqn,'(0)')|endsWith(eqn,'(1)')); %pde boundary terms
    pde_icoeffs = eqn(endsWith(eqn,'ds')); %pde integral terms
    
    % place terms in appropriate structure
    for k=1:length(ode_coeffs)
        term =ode_coeffs{k};
        Rstate = regexp(term,'x[0-9]+','match');
        coeff = regexprep(term,'x[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.x = statenameMap(keyval);
        PDE_struct.z{iLoc}.term{k} = tmp;
        clear tmp;
    end
    % Disturbance terms w<i>: coefficient and input index
    for k=1:length(dis_coeffs)
        term =dis_coeffs{k};
        Rstate = regexp(term,'w[0-9]+','match');
        coeff = regexprep(term,'w[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.w = statenameMap(keyval);
        PDE_struct.z{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Control-input terms u<i>: coefficient and input index
    for k=1:length(inp_coeffs)
        term =inp_coeffs{k};
        Rstate = regexp(term,'u[0-9]+','match');
        coeff = regexprep(term,'u[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.u = statenameMap(keyval);
        PDE_struct.z{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Boundary values X<i>(0), X<i>(1): coefficient, location loc, order D
    for k=1:length(pde_bcoeffs)
        term =pde_bcoeffs{k};
        Rstate = regexp(term,'X[0-9]+\(.\)','match');
        coeff = regexprep(term,'X[0-9]+\(.\)','');
        % Order of the s-derivative: 0, 1, or k for \partial_s^k
        der_order = regexp(coeff,'(\\partial_s)+(\^.)?','match');
        if isempty(der_order)
            der_order = 0;
        elseif strcmp(der_order{1},'\partial_s')
            der_order = 1;
        else
            der_order = str2double(regexprep(der_order{1},'\\partial_s\^',''));
        end
        coeff = regexprep(coeff,'(\\partial_s)+(\^.)?','');
        tmp.C = eval_b(coeff);
        tmp.loc = regexprep(Rstate{1},'X[0-9]+\(',''); tmp.loc = str2double(tmp.loc(1:end-1));
        keyval = regexprep(Rstate{1},'\(.\)','');
        tmp.x = statenameMap(keyval);
        tmp.D = der_order;
        PDE_struct.z{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Integral terms of X<i>: coefficient, limits I{1}, derivative order D
    for k=1:length(pde_icoeffs)
        term =pde_icoeffs{k};
        Rstate = regexp(term,'X[0-9]+','match');
        coeff = regexprep(term,'(\\int_.\^.)+|(*X[0-9]+(\(s\)ds)','');
        % Order of the s-derivative: 0, 1, or k for \partial_s^k
        der_order = regexp(coeff,'(\\partial_s)+(\^.)?','match');
        if isempty(der_order)
            der_order = 0;
        elseif strcmp(der_order{1},'\partial_s')
            der_order = 1;
        else
            der_order = str2double(regexprep(der_order{1},'\\partial_s\^',''));
        end
        coeff = regexprep(coeff,'(\\partial_s)+(\^.)?','');
        tmp.C = eval_b(coeff);
        tmp.I{1} = [0,1];
        keyval = Rstate{1};
        tmp.x = statenameMap(keyval);
        tmp.D = der_order;
        PDE_struct.z{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    iLoc = iLoc+1;
end
% get observed output params
iLoc = 1;
for i=metadata(4)+1:metadata(5)-1
    eqn = PDE_text{i}; % get line
    Loutput = regexp(eqn,'[\w]+ =','match'); % get Loutput
    
    % separate LHS and RHS, also remove partial derivative in time
    eqn = regexprep(eqn,Loutput,'');
    Lstate = Loutput{1}(1:end-2); % remove equal sign
    PDE_struct.y{iLoc}.term = [];
    % separate terms on RHS
    eqn = split(eqn);
    ode_coeffs = eqn(contains(eqn,'x')); % ode terms 
    dis_coeffs = eqn(contains(eqn,'w')); % disturbance terms
    inp_coeffs = eqn(contains(eqn,'u')); % input terms
    pde_bcoeffs = eqn(endsWith(eqn,'(0)')|endsWith(eqn,'(1)')); %pde boundary terms
    pde_icoeffs = eqn(endsWith(eqn,'ds')); %pde integral terms
    
    % place terms in appropriate structure
    for k=1:length(ode_coeffs)
        term =ode_coeffs{k};
        Rstate = regexp(term,'x[0-9]+','match');
        coeff = regexprep(term,'x[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.x = statenameMap(keyval);
        PDE_struct.y{iLoc}.term{k} = tmp;
        clear tmp;
    end
    % Disturbance terms w<i>: coefficient and input index
    for k=1:length(dis_coeffs)
        term =dis_coeffs{k};
        Rstate = regexp(term,'w[0-9]+','match');
        coeff = regexprep(term,'w[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.w = statenameMap(keyval);
        PDE_struct.y{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Control-input terms u<i>: coefficient and input index
    for k=1:length(inp_coeffs)
        term =inp_coeffs{k};
        Rstate = regexp(term,'u[0-9]+','match');
        coeff = regexprep(term,'u[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.u = statenameMap(keyval);
        PDE_struct.y{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Boundary values X<i>(0), X<i>(1): coefficient, location loc, order D
    for k=1:length(pde_bcoeffs)
        term =pde_bcoeffs{k};
        Rstate = regexp(term,'X[0-9]+\(.\)','match');
        coeff = regexprep(term,'X[0-9]+\(.\)','');
        % Order of the s-derivative: 0, 1, or k for \partial_s^k
        der_order = regexp(coeff,'(\\partial_s)+(\^.)?','match');
        if isempty(der_order)
            der_order = 0;
        elseif strcmp(der_order{1},'\partial_s')
            der_order = 1;
        else
            der_order = str2double(regexprep(der_order{1},'\\partial_s\^',''));
        end
        coeff = regexprep(coeff,'(\\partial_s)+(\^.)?','');
        tmp.C = eval_b(coeff);
        tmp.loc = regexprep(Rstate{1},'X[0-9]+\(',''); tmp.loc = str2double(tmp.loc(1:end-1));
        keyval = regexprep(Rstate{1},'\(.\)','');
        tmp.x = statenameMap(keyval);
        tmp.D = der_order;
        PDE_struct.y{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    clear tmp;
    % Integral terms of X<i>: coefficient, limits I{1}, derivative order D
    for k=1:length(pde_icoeffs)
        term =pde_icoeffs{k};
        Rstate = regexp(term,'X[0-9]+','match');
        coeff = regexprep(term,'(\\int_.\^.)+|(*X[0-9]+(\(s\)ds)','');
        % Order of the s-derivative: 0, 1, or k for \partial_s^k
        der_order = regexp(coeff,'(\\partial_s)+(\^.)?','match');
        if isempty(der_order)
            der_order = 0;
        elseif strcmp(der_order{1},'\partial_s')
            der_order = 1;
        else
            der_order = str2double(regexprep(der_order{1},'\\partial_s\^',''));
        end
        coeff = regexprep(coeff,'(\\partial_s)+(\^.)?','');
        tmp.C = eval_b(coeff);
        tmp.I{1} = [0,1];
        keyval = Rstate{1};
        tmp.x = statenameMap(keyval);
        tmp.D = der_order;
        PDE_struct.y{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    iLoc = iLoc+1;
end
% get BC params
iLoc = 1;
for i=metadata(5)+1:length(PDE_Text)
    eqn = PDE_text{i}; % get line
    
    % separate LHS and RHS, 
    eqn = regexprep(eqn,'0 = ','');
    
    % identify BC number
    BCn = i-metadata(5); %find BC number to arrange the matrices
    % separate terms on RHS
    eqn = split(eqn);
    ode_coeffs = eqn(contains(eqn,'x')); % ode terms 
    dis_coeffs = eqn(contains(eqn,'w')); % disturbance terms
    inp_coeffs = eqn(contains(eqn,'u')); % input terms
    pde_bcoeffs = eqn(endsWith(eqn,'(0)')|endsWith(eqn,'(1)')); %pde boundary terms
    pde_icoeffs = eqn(endsWith(eqn,'ds')); %pde integral terms
    PDE_struct.BC{iLoc}.term = [];
    % place terms in appropriate structure
    for k=1:length(ode_coeffs)
        term =ode_coeffs{k};
        Rstate = regexp(term,'x[0-9]+','match');
        coeff = regexprep(term,'x[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.x = statenameMap(keyval);
        PDE_struct.BC{iLoc}.term{k} = tmp;
        clear tmp;
    end
    % Disturbance terms w<i>: coefficient and input index
    for k=1:length(dis_coeffs)
        term =dis_coeffs{k};
        Rstate = regexp(term,'w[0-9]+','match');
        coeff = regexprep(term,'w[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.w = statenameMap(keyval);
        PDE_struct.BC{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Control-input terms u<i>: coefficient and input index
    for k=1:length(inp_coeffs)
        term =inp_coeffs{k};
        Rstate = regexp(term,'u[0-9]+','match');
        coeff = regexprep(term,'u[0-9]+','');
        tmp.C = eval_b(coeff);
        keyval = Rstate{1};
        tmp.u = statenameMap(keyval);
        PDE_struct.BC{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Boundary values X<i>(0), X<i>(1): coefficient, location loc, order D
    for k=1:length(pde_bcoeffs)
        term =pde_bcoeffs{k};
        Rstate = regexp(term,'X[0-9]+\(.\)','match');
        coeff = regexprep(term,'X[0-9]+\(.\)','');
        % Order of the s-derivative: 0, 1, or k for \partial_s^k
        der_order = regexp(coeff,'(\\partial_s)+(\^.)?','match');
        if isempty(der_order)
            der_order = 0;
        elseif strcmp(der_order{1},'\partial_s')
            der_order = 1;
        else
            der_order = str2double(regexprep(der_order{1},'\\partial_s\^',''));
        end
        coeff = regexprep(coeff,'(\\partial_s)+(\^.)?','');
        tmp.C = eval_b(coeff);
        tmp.loc = regexprep(Rstate{1},'X[0-9]+\(',''); tmp.loc = str2double(tmp.loc(1:end-1));
        keyval = regexprep(Rstate{1},'\(.\)','');
        tmp.x = statenameMap(keyval);
        tmp.D = der_order;
        PDE_struct.BC{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    % Integral terms of X<i>: coefficient, limits I{1}, derivative order D
    for k=1:length(pde_icoeffs)
        term =pde_icoeffs{k};
        Rstate = regexp(term,'X[0-9]+','match');
        coeff = regexprep(term,'(\\int_.\^.)+|(*X[0-9]+(\(s\)ds)','');
        % Order of the s-derivative: 0, 1, or k for \partial_s^k
        der_order = regexp(coeff,'(\\partial_s)+(\^.)?','match');
        if isempty(der_order)
            der_order = 0;
        elseif strcmp(der_order{1},'\partial_s')
            der_order = 1;
        else
            der_order = str2double(regexprep(der_order{1},'\\partial_s\^',''));
        end
        coeff = regexprep(coeff,'(\\partial_s)+(\^.)?','');
        tmp.C = eval_b(coeff);
        tmp.I{1} = [0,1];
        keyval = Rstate{1};
        tmp.x = statenameMap(keyval);
        tmp.D = der_order;
        PDE_struct.BC{iLoc}.term{end+1} = tmp;
        clear tmp;
    end
    iLoc = iLoc+1;
end


% Initialize and check the PDE
PDE_struct = initialize_PIETOOLS_PDE(PDE_struct);
end

function val = eval_b(exp)
% VAL = EVAL_B(EXP) evaluates the coefficient text EXP of a term in the
% polynomial variables s and theta. A trailing '*' is dropped; an empty
% text, '+' or '*' gives 1; a text that fails is evaluated with a leading
% '(' (a text with an unmatched ')'); if that fails too, an error names
% the text.
%
% INPUT
% - exp:    char, the coefficient text;
%
% OUTPUT
% - val:    'double' or 'polynomial', the coefficient.

% Variables the text may contain
pvar s theta;
% if strcmp(exp(end),'*')                                                   % MMP, 11/20/2026 (was)
if ~isempty(exp) && strcmp(exp(end),'*')    % '' errored at exp(end)        % MMP, 11/20/2026
    exp = exp(1:end-1);
end
if strcmp(exp,'+')||strcmp(exp,'*')||isempty(exp)
    val = 1;
else
    try 
        val = eval(exp);
    catch err
        % BEGIN MMP, 11/20/2026: the retry with a leading '(' was taken only
        % for one exact message text (mismatched delimiters; R2025b words
        % the error of an unmatched ')' differently), and every other
        % failure gave 0 without notice. Now every failure is retried,
        % which succeeds only for an unmatched ')', and a failed retry is
        % an error naming the coefficient.
%       if strcmp(err.message, 'Error: Invalid expression. When calling a function or indexing a variable, use parentheses. Otherwise, check for mismatched delimiters.') % MMP, 11/20/2026 (was)
%           val = eval(['(',exp]);                                          % MMP, 11/20/2026 (was)
%       else                                                                % MMP, 11/20/2026 (was)
%           val=0;                                                          % MMP, 11/20/2026 (was)
%       end                                                                 % MMP, 11/20/2026 (was)
        try                                                                 % MMP, 11/20/2026
            val = eval(['(',exp]);                                          % MMP, 11/20/2026
        catch                                                               % MMP, 11/20/2026
            error('The coefficient ''%s'' of a term could not be evaluated: %s',exp,err.message) % MMP, 11/20/2026
        end                                                                 % MMP, 11/20/2026
        % END MMP, 11/20/2026
    end
end
end
