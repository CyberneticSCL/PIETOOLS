function PDE_params = parsePDEParams_old(varlist, PDE_Text, metadata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parsePDEParams.m     PIETOOLS 2022
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

PDE.ODE = struct(); PDE.PDE = struct(); PDE.BC = struct();
%initialize PDE parameters
% init_PDE_params;
pvar s theta;
PDE.vars = [s;theta];
% get PDE params
PDE.PDE.Bpx ={};PDE.PDE.Bpw ={};PDE.PDE.Bpu ={};PDE.PDE.App ={};PDE.PDE.Bpb ={};PDE.PDE.Bpi ={};
for i=metadata(1)+1:metadata(2)-1
    eqn = PDE_text{i}; % get line
    Lstate = regexp(eqn,'[\w]+ =','match'); % get Lstate
    % separate LHS and RHS, also remove partial derivative in time
    eqn = regexprep(eqn,['\\partial_t ',Lstate],'');
    Lstate = Lstate{1}(1:end-2); % remove equal sign
    
    % identify Lstate number
    Xn = str2double(Lstate(2:end)); %find output number to arrange the matrices
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        PDE.PDE.Bpx{end+1} = tmp;
    end
    for k=1:length(dis_coeffs)
        term =dis_coeffs{k};
        Rstate = regexp(term,'w[0-9]+','match');
        coeff = regexprep(term,'w[0-9]+','');
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        PDE.PDE.Bpw{end+1} = tmp;
    end
    for k=1:length(inp_coeffs)
        term =inp_coeffs{k};
        Rstate = regexp(term,'u[0-9]+','match');
        coeff = regexprep(term,'u[0-9]+','');
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        PDE.PDE.Bpu{end+1} = tmp;        
    end
    for k=1:length(pde_coeffs)
        term =pde_coeffs{k};
        Rstate = regexp(term,'X[0-9]+\(.\)','match');
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        tmp.D = der_order;
        PDE.PDE.App{end+1} = tmp;
    end
    clear tmp;
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        tmp.D = der_order;
        PDE.PDE.Bpb{end+1} = tmp;
    end
    clear tmp;
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        tmp.D = der_order;
        tmp.I = [eval_b(a{1}(2:end-1)) eval_b(b{1}(2:end))];
        PDE.PDE.Bpi{end+1} = tmp;
    end
    clear tmp;
end
% get ODE params
PDE.ODE.A={}; PDE.ODE.Bxw={};PDE.ODE.Bxu={};PDE.ODE.Bxb={};PDE.ODE.Bxi={};
for i=metadata(2)+1:metadata(3)-1
    eqn = PDE_text{i}; % get line
    Lstate = regexp(eqn,'[\w]+ =','match'); % get Lstate
    % separate LHS and RHS, also remove partial derivative in time
    eqn = regexprep(eqn,['\\partial_t ',Lstate],'');
    Lstate = Lstate{1}(1:end-2); % remove equal sign
    
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        PDE.ODE.A{end+1} = tmp;
    end
    for k=1:length(dis_coeffs)
        term =dis_coeffs{k};
        Rstate = regexp(term,'w[0-9]+','match');
        coeff = regexprep(term,'w[0-9]+','');
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        PDE.ODE.Bxw{end+1} = tmp;
    end
    for k=1:length(inp_coeffs)
        term =inp_coeffs{k};
        Rstate = regexp(term,'u[0-9]+','match');
        coeff = regexprep(term,'u[0-9]+','');
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        PDE.ODE.Bxu{end+1} = tmp;
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        tmp.D = der_order;
        PDE.ODE.Bxb{end+1} = tmp;
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        tmp.D = der_order;
        PDE.ODE.Bxi{end+1} = tmp;
        clear tmp;
    end
end
% get regulated output params
PDE.ODE.Cz={}; PDE.ODE.Dzw={};PDE.ODE.Dzu={};PDE.ODE.Dzb={};PDE.ODE.Dzi={};
for i=metadata(3)+1:metadata(4)-1
    eqn = PDE_text{i}; % get line
    Loutput = regexp(eqn,'[\w]+ =','match'); % get Loutput
    
    % separate LHS and RHS, also remove partial derivative in time
    eqn = regexprep(eqn,Loutput,'');
    Lstate = Loutput{1}(1:end-2); % remove equal sign
    
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        PDE.ODE.Cz{end+1} = tmp;
    end
    for k=1:length(dis_coeffs)
        term =dis_coeffs{k};
        Rstate = regexp(term,'w[0-9]+','match');
        coeff = regexprep(term,'w[0-9]+','');
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        PDE.ODE.Dzw{end+1} = tmp;
    end
    for k=1:length(inp_coeffs)
        term =inp_coeffs{k};
        Rstate = regexp(term,'u[0-9]+','match');
        coeff = regexprep(term,'u[0-9]+','');
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        PDE.ODE.Dzu{end+1} = tmp;
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        tmp.D = der_order;
        PDE.ODE.Dzb{end+1} = tmp;
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        tmp.D = der_order;
        PDE.ODE.Dzi{end+1} = tmp;
    end
    clear tmp;
end
% get observed output params
PDE.ODE.Cy={}; PDE.ODE.Dyw={};PDE.ODE.Dyu={};PDE.ODE.Dyb={};PDE.ODE.Dyi={};
for i=metadata(4)+1:metadata(5)-1
    eqn = PDE_text{i}; % get line
    Loutput = regexp(eqn,'[\w]+ =','match'); % get Loutput
    
    % separate LHS and RHS, also remove partial derivative in time
    eqn = regexprep(eqn,Loutput,'');
    Lstate = Loutput{1}(1:end-2); % remove equal sign
    
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        PDE.ODE.Cy{end+1} = tmp;
    end
    for k=1:length(dis_coeffs)
        term =dis_coeffs{k};
        Rstate = regexp(term,'w[0-9]+','match');
        coeff = regexprep(term,'w[0-9]+','');
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        PDE.ODE.Dyw{end+1} = tmp;
    end
    for k=1:length(inp_coeffs)
        term =inp_coeffs{k};
        Rstate = regexp(term,'u[0-9]+','match');
        coeff = regexprep(term,'u[0-9]+','');
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        PDE.ODE.Dyu{end+1} = tmp;
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        tmp.D = der_order;
        PDE.ODE.Dyb{end+1} = tmp;
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = Lstate;
        tmp.Rstate = Rstate{1};
        tmp.D = der_order;
        PDE.ODE.Dyi{end+1} = tmp;
    end
    clear tmp;
end
% get BC params
PDE.BC.Ex = {};PDE.BC.Ew = {};PDE.BC.Eu = {};PDE.BC.Ebb = {};PDE.BC.Ebp = {};
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
    
    % place terms in appropriate structure
    for k=1:length(ode_coeffs)
        term =ode_coeffs{k};
        Rstate = regexp(term,'x[0-9]+','match');
        coeff = regexprep(term,'x[0-9]+','');
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = BCn;
        tmp.Rstate = Rstate{1};
        PDE.BC.Ex{end+1} = tmp;
    end
    for k=1:length(dis_coeffs)
        term =dis_coeffs{k};
        Rstate = regexp(term,'w[0-9]+','match');
        coeff = regexprep(term,'w[0-9]+','');
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = BCn;
        tmp.Rstate = Rstate{1};
        PDE.BC.Ew{end+1} = tmp;
    end
    for k=1:length(inp_coeffs)
        term =inp_coeffs{k};
        Rstate = regexp(term,'u[0-9]+','match');
        coeff = regexprep(term,'u[0-9]+','');
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = BCn;
        tmp.Rstate = Rstate{1};
        PDE.BC.Eu{end+1} = tmp;
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = BCn;
        tmp.Rstate = Rstate{1};
        tmp.D = der_order;
        PDE.BC.Ebb{end+1} = tmp;
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
        tmp.coeff = eval_b(coeff);
        tmp.Lstate = BCn;
        tmp.Rstate = Rstate{1};
        tmp.D = der_order;
        PDE.BC.Ebp{end+1} = tmp;
    end
    clear tmp;
end

PDE_params = getPDEstructure(PDE);
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
