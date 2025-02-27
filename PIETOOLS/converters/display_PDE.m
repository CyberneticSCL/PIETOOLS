function [] = display_PDE(PDE,name)
% display_PDE(PDE,name) displays the PDE defined by the struct or
% pde_struct object "PDE", with name "name".
%
% INPUTS:
% PDE:  A struct or pde_struct class object, defining a PDE in a manner as
%       outlined in "initialize_PIETOOLS_PDE".
% name: str, specifying the name of the PDE. (not currently used...)
% 
% OUTPUTS:
% Command line display of PDE, as well as the output equations and BCs.
%
% NOTES:
% - Although spatial variable names may be specified in the PDE structure,
%   these will not be used in the display. Instead, s_i will be used to
%   denote primary variables, and phi_i to denote dummy variables.
% - Coefficient matrices term{j}.C will not be shown, unless they are a
%   constant scalar value. If they are matrix-valued and/or polynomial, "C"
%   will be used to denote the factor, to avoid the expression from getting
%   too messy.
% - Display is only possible for initialized PDEs, as the fields x_tab
%   through BC_tab constructed in "initialize_PIETOOLS_PDE" are used to
%   establish variable dependence of the different states, inputs, and
%   outputs. An unfinished PDE will likely produce issues in the
%   initialization function, so only complete systems can be displayed.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 07/08/2022: Initial coding;
% DJ, 12/30/2024: Bugfix for display of equation number i in C_{i,j};

% Initialize the system
if isa(PDE,'sys')
    try PDE = PDE.params;
    catch
        disp(PDE);
        return
    end
elseif isa(PDE,'struct')
    PDE = pde_struct(PDE);
end
try PDE = initialize(PDE,true);
catch
    error(['The presented PDE is not finished, or not appropriate.',...
            ' Please check that the PDE is properly specified, and run "initialize_PIETOOLS_PDE" to check for any errors.'])
end

%fprintf([name,' = \n']);

% Maximal number of characters we allow on a single line in the command 
% window.
width_max = 150;


% % % Set up a "library" of UNICODE characters
% % Greek letters
%omega = '\x3C9';
%theta = '\x03B8';
%rho = ' ';
phi = '\x03D5';
tau = '\x03C4';

% % Subscripts
sub_num = mat2cell([repmat('\x208',[10,1]),num2str((0:9)')],ones(10,1),6);
% sub_min = '\x208B';
% sub_a = '\x2090';
sub_i = '\x1D62';
%sub_j = '\x2C7C';
sub_k = '\x2096';
sub_s = '\x209B';
sub_t = '\x209C';
% sub_rho = '\x1D68';
sub_phi = '\x1D69';

% % Superscripts
sup_num = cell(10,1);       
sup_num{1} = '\x2070';      % Superscript 0
sup_num{2} = '\xB9';        % Superscript 1
sup_num{3} = '\xB2';        % etc.
sup_num{4} = '\xB3';
sup_num{5} = '\x2074';
sup_num{6} = '\x2075';
sup_num{7} = '\x2076';
sup_num{8} = '\x2077';
sup_num{9} = '\x2078';
sup_num{10} = '\x2079';
%sup_min = '\x207B';

%sup_min = '\x207B';
%sup_a = '\x0363';
% https://rupertshepherd.info/resource_pages/superscript-letters-in-unicode
% sup_a = '\x1D43';
% sup_b = '\x1D47';
% sup_c = '\x1D9C';
% sup_d = '\x1D48';
% sup_e = '\x1D49';
% sup_f = '\x1DA0';
% sup_g = '\x1D4D';
% sup_h = '\x02B0';
% sup_i = '\x2071';
% sup_j = '\x02B2';
% sup_k = '\x1D4F';
% sup_l = '\x02E1';
% sup_m = '\x1D50';
% sup_n = '\x207F';
% sup_o = '\x1D52';
% sup_p = '\x1D56';
% sup_q = 'q';        % No superscript for q?
% sup_r = '\x02B3';
sup_s = '\x02E2';
% sup_t = '\x1D57';
% sup_u = '\x1D58';
% sup_v = '\x1D5B';
% sup_w = '\x02B7';
% sup_x = '\x02E3';
% sup_y = '\x02B8';
% sup_z = '\x1DBB';

% % Other symbols
partial = '\x2202';
%int = '\x222B';
%xdot = '\x1E8B';
xdot = [partial,sub_t,' x'];
%xddot = '\x1E8D';
xddot = [partial,sub_t,sup_num{3},' x'];


% Establish cell with (UNICODE) characters for the primary variables
% {'s_1';'s_2';...}, the dummy variables {'phi_1';'phi_2';...}, the
% superscript primary variables {'^s1';'^s2';...}, the subscript primary
% variables {'_s1';'_s2';...}, and the subscript dummy variables 
% {'_phi1';'_phi2';...}
nvars = size(PDE.vars,1);
if nvars==1
    % If we have a single spatial variable, there is no need to use
    % subscripts on the variables.
    var1_list = {'s'};
    var2_list = {phi};
    sub_var1_list = {sub_s};
    sub_var2_list = {sub_phi};
    sup_var1_list = {sup_s};
else
    var1_list = cell(nvars,1);
    var2_list = cell(nvars,1);
    sup_var1_list = cell(nvars,1);
    sub_var1_list = cell(nvars,1);
    sub_var2_list = cell(nvars,1);
    for vv = 1:nvars
        sub_var_num = cell2mat(sub_num(str2num(num2str(vv)')+1)');
        var1_list{vv} = ['s',sub_var_num];
        var2_list{vv} = [phi,sub_var_num];
        sub_var1_list{vv} = [sub_s,sub_var_num];
        sup_var1_list{vv} = [sup_s,cell2mat(sup_num(str2num(num2str(vv)')+1)')];
        sub_var2_list{vv} = [sub_phi,sub_var_num];
    end
end
ndelays = size(PDE.tau,1);
if ndelays==1
    tau_list = {tau};
else
    tau_list = cell(ndelays,1);
    for vv = 1:ndelays
        sub_var_num = cell2mat(sub_num(str2num(num2str(vv)')+1)');
        tau_list{vv} = [tau,sub_var_num];
    end
end

var_list = [var1_list, var2_list, sub_var1_list, sub_var2_list, sup_var1_list];

% If the coefficients in the terms are too big, they will be displayed as
% "C_{ij}" instead. We use "use_Cij" to keep track of whether or not this
% is done at all.
use_Cij = false;

% % % Display the equations of each type, starting with state equations,
% % % observed outputs, regulated outputs, and then BCs
eq_types = {'x';'y';'z';'BC'};
objs_LHS = {'x';'y';'z';'0'};
eq_num_strt = 0;                                                            % DJ, 12/30/2024
for kk=1:numel(eq_types)
    eq_type_kk = eq_types{kk};
    obj_LHS_kk = objs_LHS{kk};
    if numel(PDE.(eq_type_kk))==0
        % No equations to display.
        continue
    end
    for eq_num=1:numel(PDE.(eq_type_kk))
        fprintf("\n  ");
           
        % % First construct the LHS of the PDE
        LHS_comp = PDE.(eq_type_kk){eq_num}; % Extract the info for the considered state
        LHS_length = 2;         % Keep track of the size of the expression on the LHS
        eq_info = PDE.([eq_type_kk,'_tab'])(eq_num,:);
        eq_ID = eq_info(1);
        
        % Set the order for the temporal derivative of the state.
        if isfield(LHS_comp,'tdiff')
            tdiff = LHS_comp.tdiff;
        elseif kk==1
            tdiff = 1;
        else
            tdiff = 0;
        end
        
        if kk<=3
            % Establish the set of variables on which the state component or
            % output depends.
            Lvar_indcs = eq_info(1,3:3+nvars-1);
            if ~any(Lvar_indcs)
                LHS_length = LHS_length + 3;
                Lvar_str = '(t)';
            else
                Lvars = var1_list(logical(Lvar_indcs));
                LHS_length = LHS_length + 3;
                Lvar_str = '(t';
                for idx=1:length(Lvars)
                    LHS_length = LHS_length + 2 + length(num2str(idx));
                    Lvar_str = [Lvar_str,',',Lvars{idx}];
                end
                Lvar_str = [Lvar_str,')'];
            end
            % Establish the (subscript) index for the state component.
            if size(PDE.([eq_type_kk,'_tab']),1)==1
                % There is only one state component --> no need to give index
                Lstate_idx = '';
            elseif max(PDE.([eq_type_kk,'_tab'])(:,1))<=9
                % The state number consists of a single decimal
                LHS_length = LHS_length + 1;
                Lstate_idx = sub_num{eq_ID+1};
            else
                % The state number consists of multiple decimals
                Lstate_idx = cell2mat(sub_num(str2num(num2str(eq_ID)')+1)');
                LHS_length = LHS_length + max(PDE.([eq_type_kk,'_tab'])(:,1));
            end
        else
            % For boundary conditions, don't display variables or index.
            Lvar_str = '';
            Lstate_idx = '';
        end
        
        % % Construct the left-hand side of the PDE
        if kk==1
            if tdiff==1
                LHS = [xdot,Lstate_idx,Lvar_str];
                LHS_length = LHS_length + 1;
            elseif tdiff==2
                LHS = [xddot,Lstate_idx,Lvar_str];
                LHS_length = LHS_length + 1;
            elseif tdiff<=9
                sup_tdiff = sup_num{tdiff+1};
                LHS = [partial,sub_t,sup_tdiff,' x',Lstate_idx,Lvar_str];
                LHS_length = LHS_length + 3;
            else
                % The order of the temporal derivative consists of muliple decimals
                sup_tdiff = cell2mat(sub_num(str2num(num2str(tdiff)')+1)');
                LHS = [partial,sub_t,sup_tdiff,' x',Lstate_idx,Lvar_str];
                LHS_length = LHS_length + 2 + length(num2str(tdiff));
            end
        else
            % For output equations or boundary conditions, there is no
            % temporal derivative taken.
            LHS = [obj_LHS_kk,Lstate_idx,Lvar_str];
            LHS_length = LHS_length + 1;
        end
        % Initialize a str for the PDE
        eq_kk = LHS;
        LHS_length = LHS_length + 2;    % Account for size of ' ='.
        
        
        % % Next, add the right-hand side
        % If no terms are present, set RHS equal to 0
        if ~isfield(LHS_comp,'term') || isempty(LHS_comp.term)
            eq_kk = [eq_kk , ' = 0'];
            fprintf(eq_kk);
            continue
        end
        % Otherwise, loop over the terms, adding each to the equation.
        for trm = 1:numel(LHS_comp.term)
            C_row_idx = eq_num_strt+eq_num;
            [term_str,use_Cij_ii] = construct_term(PDE,LHS_comp,C_row_idx,eq_info,trm,var_list,tau_list);
            use_Cij = use_Cij || use_Cij_ii;
            eq_kk_new = [eq_kk, term_str];
            
            % If the size of the equation becomes too large, start a new line.
            if trm==numel(LHS_comp.term)
                fprintf([eq_kk_new,';'])
            elseif length(eq_kk_new) > 4*width_max
                fprintf(eq_kk_new)
                fprintf('\n');
                eq_kk = repmat(' ',1,LHS_length+5);
            else
                eq_kk = eq_kk_new;
            end
        end
    end
    % Keep track of how many equations we have already had for display of
    % coefficients;
    eq_num_strt = eq_num_strt + numel(PDE.(eq_type_kk));                    % DJ, 12/30/2024
    % Add some space between different equation types.
    fprintf('\n')
end


if use_Cij
    fprintf(['\n Call "PDE.C{i,k}" to see the value of coefficients C',sub_i,sub_k,' as in the displayed equations.\n'])
end

fprintf('\n')

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %






%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [term_str,use_Cij] = construct_term(PDE,eq,eq_num,eq_info,term_num,var_list,tau_list)
% Build a str to represent term number "trm" in the equation "eq" of the
% PDE structure "PDE".

PDE_term = eq.term{term_num};

var1_list = var_list(:,1);
var2_list = var_list(:,2);
sub_var1_list = var_list(:,3);
sub_var2_list = var_list(:,4);
sup_var1_list = var_list(:,5);

nvars = length(var1_list);
ndelays = length(tau_list);
% Extract varnames one at a time, since "polynomial" class stores varnames
% alphabetically.
var1_name = cell(nvars,1);
var2_name = cell(nvars,1);
for vv=1:nvars
    var1_name{vv} = PDE.vars(vv,1).varname{1};
    var2_name{vv} = PDE.vars(vv,2).varname{1};
end
tau_name = cell(ndelays,1);
for vv=1:ndelays
    tau_name{vv} = PDE.tau(vv,1).varname{1};
end
%has_vars_eq = eq_info(3:2+nvars);

% % % Define a number of symbols we'll need to display the system
% Subscripts
sub_num = mat2cell([repmat('\x208',[10,1]),num2str((0:9)')],ones(10,1),6);
sub_min = '\x208B';
sub_a = '\x2090';
sup_b = '\x1D47';
sub_T = '\x1D1B';
sub_lp = '\x208D';
sub_rp = '\x208E';

% Superscripts
sup_num = cell(10,1);    % Superscripts
sup_num{1} = '\x2070';  % Superscript 0
sup_num{2} = '\xB9';
sup_num{3} = '\xB2';
sup_num{4} = '\xB3';
sup_num{5} = '\x2074';
sup_num{6} = '\x2075';
sup_num{7} = '\x2076';
sup_num{8} = '\x2077';
sup_num{9} = '\x2078';
sup_num{10} = '\x2079';
sup_min = '\x207B';

% Partial derivative and integral symbol
partial = '\x2202';
int_symbol = '\x222B';
%tau = '\x1D70F';
        
% % % Establish which state or input the term involves
% % First check what type of object (state or input) is involved.
if isfield(PDE_term,'x')
    obj_RHS = 'x';
    Rstate_trm = 'x';
elseif isfield(PDE_term,'w')
    obj_RHS = 'w';
    Rstate_trm = 'w';
elseif isfield(PDE_term,'u')
    obj_RHS = 'u';
    Rstate_trm = 'u';
end
% Then determine the index and ID associated to the object.
tab_row = PDE_term.(obj_RHS);
rr_ID = PDE.([obj_RHS,'_tab'])(tab_row,1);
%rr = PDE_term.(obj_RHS);
if numel(PDE.(obj_RHS))==1
    % There is only one component of the considered type --> no need to give index
    Rstate_idx = '';
elseif rr_ID<=9
    % The state index consists of a single decimal
    Rstate_idx = sub_num{rr_ID+1};
else
    % The state number consists of multiple decimals
    Rstate_idx = cell2mat(sub_num(str2num(num2str(rr_ID)')+1)');
end
%tab_row = find(PDE.x_tab(:,1)==rr);
%tab_row = rr;

% % Establish which variables the state or input depends on
Rvar_indcs = logical(PDE.([obj_RHS,'_tab'])(tab_row,3:3+nvars-1));
%     if isa(PDE.x{tab_row}.vars,'polynomial')
%         Rvar1_name = PDE.x{rr}.vars(:,1).varname;
%         Rvar2_name = PDE.x{rr}.vars(:,2).varname;
%     else
%         Rvar1_name = cell(0,1);
%         Rvar2_name = cell(0,1);
%     end
Rstate_trm = [Rstate_trm,Rstate_idx];
Rvar1_list = var1_list(Rvar_indcs);
Rvar2_list = var2_list(Rvar_indcs);
sub_Rvar1_list = sub_var1_list(Rvar_indcs);
sub_Rvar2_list = sub_var2_list(Rvar_indcs);
%sup_Rvar1_list = sup_var1_list(Rvar_indcs);

%has_vars_eq = has_vars_eq(Rvar_indcs);


% % % Display the delay
int_trm = '';
dtheta_trm = '';
Rvar1_str = '(t';
if isfield(PDE_term,'delay')
    delay = PDE_term.delay;
    if (isa(delay,'double') || isdouble(delay)) && double(delay)~=0
        % % Delay is a real value a, just display x(t-a)
        Rvar1_str = [Rvar1_str,'-',num2str(double(PDE_term.delay))];
    elseif isa(delay,'polynomial') && ~isdouble(delay)
        % % Delay is distributed, we have to build an integral.
        % First determine which delay variable is involved, and display
        % x(t-tau)
        delay = PDE_term.delay;
        tau_indx = ismember(tau_name,delay.varname{1});
        if sum(tau_indx)>1
            error('The term involves multiple delays...')
        end
        tau_indx = find(tau_indx,1,'first');
        tau1 = tau_list(tau_indx);
        Rvar1_str = [Rvar1_str,'+',tau1{1}];
        
        % % Now for the integral
        % Add lower limit
        L_kk = abs(double(PDE.tau(tau_indx,2)));
        if (round(L_kk)-L_kk)~=0
            % The lower limit is not integer --> indicate
            % as just a_kk
            int_trm = [int_trm,sub_min,sub_T,sub_num{tau_indx+1}];
        else
            int_trm = [int_trm,sub_min];
            if L_kk<=9
                % The value consists of a single decimal
                int_trm = [int_trm,sub_num{L_kk+1}];
            else
                % The value consists of multiple decimals
                int_trm = [int_trm, cell2mat(sub_num(str2num(num2str(L_kk)')+1)')];
            end
        end
        int_trm = [int_trm,int_symbol,sup_num{1}];
        % Add dtau at end of integral
        dtheta_trm = ['d',tau1{1},' ',dtheta_trm,];
    end
end


% % % Display the integrals
use_theta_trm = false(length(Rvar1_list),1);
if isfield(PDE_term,'I') && ~isempty(PDE_term.I)
    for kk=1:numel(PDE_term.I)
        if isempty(PDE_term.I{kk})
            continue
        end
        
        % % Add lower limit of integral.
        L_kk = PDE_term.I{kk}(1);
        if isa(L_kk,'double') || (isa(L_kk,'polynomial') && isdouble(L_kk))
            L_kk = double(L_kk);
            if (round(L_kk)-L_kk)~=0
                % The lower limit is not integer --> indicate
                % as just a_kk
                int_trm = [int_trm,sub_a,sub_num{kk+1}];
            else
                if L_kk < 0
                    % Negative value: add a minus
                    int_trm = [int_trm,sub_min];
                    L_kk = abs(L_kk);
                end
                if L_kk<=9
                    % The value consists of a single decimal
                    int_trm = [int_trm,sub_num{L_kk+1}];
                else
                    % The value consists of multiple decimals
                    int_trm = [int_trm, cell2mat(sub_num(str2num(num2str(L_kk)')+1)')];
                end
            end
        else
            % The lower limit is a variable, set to s_kk
            Lkk_var1_indx = ismember(var1_name,L_kk.varname{1});
            int_trm = [int_trm, sub_var1_list{Lkk_var1_indx}];
        end
        
        % % Add integral symbol
        int_trm = [int_trm,int_symbol];
        
        % % Add upper limit of integral.
        U_kk = PDE_term.I{kk}(2);
        if isa(U_kk,'double') || (isa(U_kk,'polynomial') && isdouble(U_kk))
            U_kk = double(U_kk);
            if (round(U_kk)-U_kk)~=0
                % The upper limit is not integer --> indicate
                % as just n_kk
                int_trm = [int_trm,sup_b,sup_num{kk+1}];
            else
                if U_kk < 0
                    % Negative value: add a minus
                    int_trm = [int_trm,sup_min];
                    U_kk = abs(U_kk);
                end
                if U_kk<=9
                    % The value consists of a single decimal
                    int_trm = [int_trm,sup_num{U_kk+1}];
                else
                    % The value consists of multiple decimals
                    int_trm = [int_trm, cell2mat(sup_num(str2num(num2str(U_kk)')+1)')];
                end
            end
        else
            % The upper limit is a variable, set to s_kk
            Ukk_var1_indx = ismember(var1_name,U_kk.varname{1});
            int_trm = [int_trm, sup_var1_list{Ukk_var1_indx}];
        end
        
        % Keep track of whether primary or dummy variable is used in
        % integral.
        use_theta_trm(kk) = true;
        if ~isfield(PDE_term,'loc') || (ispvar(PDE_term.loc(kk)) && ismember(PDE_term.loc(kk).varname{1},var1_name))
            % Use standard spatial variable for integration
            dtheta_trm = ['d',Rvar1_list{kk},' ',dtheta_trm,];
        else
            % Use dummy variable for integration.
            Rvar1_list{kk} = Rvar2_list{kk}; % Evaluate state at dummy var
            sub_Rvar1_list{kk} = sub_Rvar2_list{kk};
            dtheta_trm = ['d',Rvar2_list{kk},' ',dtheta_trm,];
        end      
    end
end
if ~isempty(int_trm)
    int_trm = [int_trm,'['];
    %int_trm = [int_trm,'['];
    dtheta_trm = [']',strtrim(dtheta_trm)];
end


% % % Add the coefficient
C_trm = '';
if term_num==1
    sign_trm = ' = ';
else
    sign_trm = ' + ';
end
if isfield(PDE_term,'C') && ~isempty(PDE_term.C)
    Cval = PDE_term.C;
    % Convert constant polynomial to double
    if isa(Cval,'polynomial') && isdouble(Cval)
        Cval = double(Cval);
    end
    if isa(Cval,'double') && all(size(Cval,1)==size(Cval,2)) && ~any(any(Cval-diag(diag(Cval)))) && ...
        ~any(any(Cval-Cval(1,1)*eye(size(Cval))))
        % % Constant scalar factor
        use_Cij = false;
        Cval = Cval(1,1);
        if Cval<0 && term_num==1
            sign_trm = ' = - ';
            Cval = abs(Cval);
        elseif Cval<0
            sign_trm = ' - ';
            Cval = abs(Cval);
        end
        if Cval~=1
            C_trm = [C_trm,num2str(Cval),' * '];
        end
    elseif isa(Cval,'double')
        % % Constant matrix-valued factor   
        use_Cij = true;
        % Add subscripts indicating the equation and term number.
        if eq_num<=9
            % The equation number consists of a single decimal
            eq_indx = sub_num{eq_num+1};
        else
            % The equation number consists of multiple decimals
            eq_indx = cell2mat(sub_num(str2num(num2str(eq_num)')+1)');
        end
        if term_num<=9
            % The equation number consists of a single decimal
            trm_indx = sub_num{term_num+1};
        else
            % The equation number consists of multiple decimals
            trm_indx = cell2mat(sub_num(str2num(num2str(term_num)')+1)');
        end
        if eq_num>9 || term_num>9
            C_trm = [C_trm,'C',sub_lp,eq_indx,',',trm_indx,sub_rp,' * '];
        else
            C_trm = [C_trm,'C',eq_indx,'',trm_indx,' * '];
        end
    else
        % % Polynomial function.
        use_Cij = true;
        Ctau_indcs = ismember(tau_name,Cval.varname);
        Ctau_list = tau_list(Ctau_indcs);
        Cvar1_indcs = ismember(var1_name,Cval.varname);
        Cvar1_list = var1_list(Cvar1_indcs);
        Cvar2_indcs = ismember(var2_name,Cval.varname);
        Cvar2_list = var2_list(Cvar2_indcs);
        % Account for delay.
        if any(Ctau_indcs)
            Cvar_str = ['(',Ctau_list{1},','];
        elseif any(Cvar1_indcs) || any(Cvar2_indcs)
            Cvar_str = '(';
        else
            Cvar_str = '';
        end
        % Account for spatial dependence.
        if any(Cvar1_indcs)
            Cvar_str = [Cvar_str,Cvar1_list{1}];
            for idx=2:length(Cvar1_list)
                Cvar_str = [Cvar_str,',',Cvar1_list{idx}];
            end
            for idx=1:length(Cvar2_list)
                Cvar_str = [Cvar_str,',',Cvar2_list{idx}];
            end
            Cvar_str = [Cvar_str,')'];
        elseif any(Cvar2_indcs)
            Cvar_str = [Cvar_str,Cvar2_list{1}];
            for idx=2:length(Cvar2_list)
                Cvar_str = [Cvar_str,',',Cvar2_list{idx}];
            end
            Cvar_str = [Cvar_str,')'];
        end
        % Add subscripts indicating the equation and term number.
        if eq_num<=9
            % The equation number consists of a single decimal
            eq_indx = sub_num{eq_num+1};
        else
            % The equation number consists of multiple decimals
            eq_indx = cell2mat(sub_num(str2num(num2str(eq_num)')+1)');
        end
        if term_num<=9
            % The equation number consists of a single decimal
            trm_indx = sub_num{term_num+1};
        else
            % The equation number consists of multiple decimals
            trm_indx = cell2mat(sub_num(str2num(num2str(term_num)')+1)');
        end
        if eq_num>9 || term_num>9
            C_trm = [C_trm,'C',sub_lp,eq_indx,',',trm_indx,sub_rp,Cvar_str,' * '];
        else
            C_trm = [C_trm,'C',eq_indx,'',trm_indx,Cvar_str,' * '];
        end
    end
%     if ~isempty(C_trm)
%         C_trm = [strtrim(C_trm),' '];
%     end
end
%         if ~isempty(C_trm)
%              C_trm = [' ',C_trm];
%         end

% % % Display the partial derivatives
D_trm = '';
if isfield(PDE_term,'D') && ~isempty(PDE_term.D)
    for kk=1:size(PDE_term.D,2)
        D_kk = PDE_term.D(kk);
        if D_kk==0
            % Don't display partial symbol to 0th power
            continue
        elseif D_kk==1
            % Don't display superscript 1 in derivative
            D_trm = [D_trm,' ',partial,sub_Rvar1_list{kk}];
        elseif D_kk<=9
            D_trm = [D_trm,' ',partial,sup_num{D_kk+1},sub_Rvar1_list{kk}];
        else
            D_sup = cell2mat(sup_num(str2num(num2str(D_kk)')+1)');
            D_trm = [D_trm,' ',partial,D_sup,sub_Rvar1_list{kk}];
        end
    end
    if ~isempty(D_trm)
        D_trm = [strtrim(D_trm),' '];
    end
end


% % % Display the location at which to evaluate the state
if (isfield(PDE_term,'loc') && ~isempty(PDE_term.loc))
    for kk=1:size(PDE_term.loc,2)
        loc_kk = PDE_term.loc(kk);
        if isa(loc_kk,'double') || (isa(loc_kk,'polynomial') && isdouble(loc_kk))
            loc_kk = double(loc_kk);
            Rvar1_list{kk} = num2str(loc_kk);
        end
    end
end
for idx=1:size(Rvar1_list,1)
    Rvar1_str = [Rvar1_str,',',Rvar1_list{idx}];
end
Rvar1_str = [Rvar1_str,')'];

    % % % Finally, add the term
    term_str = [sign_trm,int_trm,C_trm,D_trm,Rstate_trm,Rvar1_str,dtheta_trm];

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %