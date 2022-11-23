function PDE_struct = parsePDEParams(varlist, PDE_Text, metadata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parsePDEParams.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is exclusively used by the PDE_GUI to convert textual PDE
% to a matlab PDE object in terms format.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following inputs must be defined passed:
%
%  varlist - variables table in the PDE gui
% PDE_text - PDE text to be parsed from the GUI
% metadata - row numbers metadata from PDE gui

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
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
Xvals = cell2mat(varlist(contains(varlist(:,1),'X'),2));
Xvals = [Xvals zeros(size(Xvals))];
for i=0:max(Xvals)
    Xvals(Xvals(:,1)==i,2) = (1:length(Xvals(Xvals(:,1)==i,2)))';
end

if ~isempty(Xlist)
    PDE.Xvaridx = containers.Map(Xlist,num2cell(Xvals,2));
else
    PDE.Xvaridx = [];
end

keys = [xlist;Xlist;wlist;ulist;zlist;ylist];
vals = [1:length(xlist) length(xlist)+1:length(xlist)+length(Xlist) 1:length(wlist) 1:length(ulist) 1:length(zlist) 1:length(ylist)]; 
statenameMap = containers.Map(keys,vals);

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
    for k=1:length(pde_coeffs)
        term =pde_coeffs{k};
        Rstate = regexp(term,'X[0-9]+\(s\)','match');
        coeff = regexprep(term,'X[0-9]+\(s\)','');
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
    for k=1:length(pde_bcoeffs)
        term =pde_bcoeffs{k};
        Rstate = regexp(term,'X[0-9]+\(.\)','match');
        coeff = regexprep(term,'X[0-9]+\(.\)','');
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
    for k=1:length(pde_icoeffs)
        term =pde_icoeffs{k};
        Rstate = regexp(term,'X[0-9]+','match');
        intlimits = regexp(term,'\\int_.\^.','match');
        b = regexp(intlimits{1},'\^.','match');
        a = regexp(intlimits{1},'_.\^','match');
        coeff = regexprep(term,'(\\int_.\^.)+|(*X[0-9]+(\(theta\)dtheta)','');
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
    for k=1:length(pde_bcoeffs)
        term =pde_bcoeffs{k};
        Rstate = regexp(term,'X[0-9]+\(.\)','match');
        coeff = regexprep(term,'X[0-9]+\(.\)','');
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
    for k=1:length(pde_icoeffs)
        term =pde_icoeffs{k};
        Rstate = regexp(term,'X[0-9]+','match');
        coeff = regexprep(term,'(\\int_.\^.)+|(*X[0-9]+(\(s\)ds)','');
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
    for k=1:length(pde_bcoeffs)
        term =pde_bcoeffs{k};
        Rstate = regexp(term,'X[0-9]+\(.\)','match');
        coeff = regexprep(term,'X[0-9]+\(.\)','');
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
    for k=1:length(pde_icoeffs)
        term =pde_icoeffs{k};
        Rstate = regexp(term,'X[0-9]+','match');
        coeff = regexprep(term,'(\\int_.\^.)+|(*X[0-9]+(\(s\)ds)','');
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
    for k=1:length(pde_bcoeffs)
        term =pde_bcoeffs{k};
        Rstate = regexp(term,'X[0-9]+\(.\)','match');
        coeff = regexprep(term,'X[0-9]+\(.\)','');
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
    for k=1:length(pde_icoeffs)
        term =pde_icoeffs{k};
        Rstate = regexp(term,'X[0-9]+','match');
        coeff = regexprep(term,'(\\int_.\^.)+|(*X[0-9]+(\(s\)ds)','');
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
    for k=1:length(pde_bcoeffs)
        term =pde_bcoeffs{k};
        Rstate = regexp(term,'X[0-9]+\(.\)','match');
        coeff = regexprep(term,'X[0-9]+\(.\)','');
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
    for k=1:length(pde_icoeffs)
        term =pde_icoeffs{k};
        Rstate = regexp(term,'X[0-9]+','match');
        coeff = regexprep(term,'(\\int_.\^.)+|(*X[0-9]+(\(s\)ds)','');
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


PDE_struct = initialize_PIETOOLS_PDE(PDE_struct);
end

function val = eval_b(exp)
pvar s theta;
if strcmp(exp(end),'*')
    exp = exp(1:end-1);
end
if strcmp(exp,'+')||strcmp(exp,'*')||isempty(exp)
    val = 1;
else
    try 
        val = eval(exp);
    catch err
        if strcmp(err.message, 'Error: Invalid expression. When calling a function or indexing a variable, use parentheses. Otherwise, check for mismatched delimiters.')
            val = eval(['(',exp]);
        else
            val=0;
        end
    end
end
end
