function [term_str,term_sign,use_Cij] = display_term(PDE,PDE_term,C_idcs)
% Build a 'char' object to display the term "PDE_term" in the 'pde_struct'
% object "PDE".
%
% INPUT
% - PDE:        'pde_struct' object representing the PDE equations (or a
%               set of loose PDE terms) to which the input term
%               "PDE_term" belongs.
% - PDE_term:   struct specifying a single term in the PDE, as per the
%               'pde_struct' format.
% - C_idcs:     1x2 array of type 'double' specifying indices (i,j) for the
%               coefficients C_{i,j} in the term if these coefficients
%               cannot be displayed explicitly.
%
% OUTPUT
% - term_str:   1xn char object representing the term
% - term_sign:  scalar integer set to 1 or -1, to specify the sign of the
%               term, when the absolute value of the coefficients is taken.
% - use_Cij:    boolean object set to true if the coefficients multiplying
%               the term are printed as "C_{ik}" rather than the actual
%               value the coefficients assume. This is necessary if C is
%               matrix-valued, or polynomial.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 06/23/2024

% % % Set up a "library" of UNICODE characters
% % Greek letters
%phi = '\x03D5';
%tau = '\x03C4';

% % Subscripts
% sub_i = '\x1D62';
% sub_k = '\x2096';
% sub_s = '\x209B';
sub_t = '\x209C';
%sub_phi = '\x1D69';


% https://rupertshepherd.info/resource_pages/superscript-letters-in-unicode
%sup_s = '\x02E2';
%sup_t = '\x1D57';
% sup_u = '\x1D58';
% sup_v = '\x1D5B';
% sup_w = '\x02B7';
% sup_x = '\x02E3';
% sup_y = '\x02B8';
% sup_z = '\x1DBB';



% % % Define a number of symbols we'll need to display the system
% Subscripts
sub_num = mat2cell([repmat('\x208',[10,1]),num2str((0:9)')],ones(10,1),6);
sub_min = '\x208B';
sub_a = '\x2090';
sup_b = '\x1D47';
sub_lp = '\x208D';      % Left parenthesis
sub_rp = '\x208E';      % Right parenthesis

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

if nargin<=2
    eq_num = 1;
    term_num = 1;
    vars = polynomial(zeros(0,1));
elseif nargin<=3
    eq_num = C_idcs(1);
    term_num = C_idcs(2);
    vars = polynomial(zeros(0,1));
end

        
% % % Establish which state or input the term involves
if isfield(PDE_term,'x')
    Robj = 'x';
    tab_row = PDE_term.x;
    % Display potential temporal derivative of state.
    if isfield(PDE_term,'tdiff') && PDE_term.tdiff>0
        tdiff_Robj = PDE_term.tdiff;
        if tdiff_Robj==1
            Rstate_trm = [partial,sub_t,' ',Robj];
        elseif tdiff_Robj<=9
            Rstate_trm = [partial,sub_t,sup_num{tdiff_Robj+1},Robj];
        else
            tD_sup = cell2mat(sup_num(str2num(num2str(tdiff_Robj)')+1)');
            Rstate_trm = [partial,sub_t,tD_sup,Robj];
        end
    else
        Rstate_trm = Robj;
    end
else
    if isfield(PDE_term,'y')
        Robj = 'y';
    elseif isfield(PDE_term,'z')
        Robj = 'z';
    elseif isfield(PDE_term,'u')
        Robj = 'u';
    elseif isfield(PDE_term,'w')
        Robj = 'w';
    else
        error("The term to display is not properly specified: it does not involve any state, input, or output variable.")
    end
    tab_row = PDE_term.(Robj);
    Rstate_trm = Robj;
end
    
% Assign an index to the state, input, or output
R_ID = PDE.(Robj){tab_row}.ID;
if R_ID<=9
    % The state index consists of a single decimal
    Rstate_idx = sub_num{R_ID+1};
else
    % The state number consists of multiple decimals
    Rstate_idx = cell2mat(sub_num(str2num(num2str(R_ID)')+1)');
end
Rstate_trm = [Rstate_trm,Rstate_idx];

% Extract spatial variables on which the component depends.
Rvars = PDE.(Robj){tab_row}.vars;
Rvarname = cell(size(Rvars,1),1);
for kk=1:length(Rvarname)
    Rvarname{kk} = Rvars(kk).varname{1};
end



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
        % First build display for x(t-tau)
        %delay = PDE_term.delay;
        error('Display of distributed delays is currently not supported.')
        % Rvar1_str = [Rvar1_str,'+',tau];
        % 
        % % % Now for the integral
        % % Add lower limit
        % L_kk = abs(double(PDE.tau(tau_indx,2)));
        % if (round(L_kk)-L_kk)~=0
        %     % The lower limit is not integer --> indicate
        %     % as just a_kk
        %     int_trm = [int_trm,sub_min,sub_T,sub_num{tau_indx+1}];
        % else
        %     int_trm = [int_trm,sub_min];
        %     if L_kk<=9
        %         % The value consists of a single decimal
        %         int_trm = [int_trm,sub_num{L_kk+1}];
        %     else
        %         % The value consists of multiple decimals
        %         int_trm = [int_trm, cell2mat(sub_num(str2num(num2str(L_kk)')+1)')];
        %     end
        % end
        % int_trm = [int_trm,int,sup_num{1}];
        % % Add dtau at end of integral
        % dtheta_trm = ['d',tau1{1},' ',dtheta_trm,];
    end
end



% % % Display the integrals
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
            % % Add integral symbol
            int_trm = [int_trm,int_symbol];
        else
            % The lower limit is a variable, set to s_kk
            int_trm = [int_trm,'{',L_kk.varname{1},'}_',int_symbol];
            %Lkk_var1_indx = ismember(var1_name,L_kk.varname{1});
            %int_trm = [int_trm, sub_var1_list{Lkk_var1_indx}];
        end

        
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
            int_trm = [int_trm,'^{',U_kk.varname{1},'}'];
            %Ukk_var1_indx = ismember(var1_name,U_kk.varname{1});
            %int_trm = [int_trm, sup_var1_list{Ukk_var1_indx}];
        end
        
        % Add the dummy variables used for integration
        varname_kk = PDE_term.loc(kk).varname{1};
        dtheta_trm = ['d',varname_kk,' ',dtheta_trm,]; 
    end
end
% Place brackets around the integrand.
if ~isempty(int_trm)
    int_trm = [int_trm,'['];
    %int_trm = [int_trm,'['];
    dtheta_trm = [']',strtrim(dtheta_trm)];
end


% % % Add the coefficient
C_trm = '';
term_sign = 1;   % Assume positive sign of coefficients.
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
        if Cval<0
            term_sign = -1;
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

        % Set display of the variables on which the function depends
        Cvar_str = ['(',Cval.varname{1}];
        for kk=2:length(Cval.varname)
            Cvar_str = [Cvar_str,',',Cval.varname{kk}];
        end
        Cvar_str = [Cvar_str,')'];

        % Add subscripts indicating the equation and term number.
        if eq_num<=9
            % The equation number consists of a single decimal
            eq_indx = sub_num{eq_num+1};
        else
            % The equation number consists of multiple decimals
            eq_indx = cell2mat(sub_num(str2num(num2str(eq_num)')+1)');
        end
        if term_num<=9
            % The term number consists of a single decimal
            trm_indx = sub_num{term_num+1};
        else
            % The term number consists of multiple decimals
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
            D_trm = [D_trm,'(',partial,'/',partial,Rvarname{kk},')'];
        elseif D_kk<=9
            D_trm = [D_trm,'(',partial,'/',partial,Rvarname{kk},')',sup_num{D_kk+1}];
        else
            D_sup = cell2mat(sup_num(str2num(num2str(D_kk)')+1)');
            D_trm = [D_trm,'(',partial,'/',partial,Rvarname{kk},')',D_sup];
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
            loc_kk = num2str(double(loc_kk));
        else
            loc_kk = loc_kk.varname{1};
        end
        Rvar1_str = [Rvar1_str,',',loc_kk];
    end
elseif ~isfield(PDE_term,'loc')
    % If no location is specified (e.g. for inputs or outputs), assume the
    % object is not evaluated at any position.
    if isfield(PDE.(Robj){tab_row},'vars')
        vars = PDE.(Robj){tab_row}.vars;
        for kk=1:size(vars,1)
            loc_kk = vars(kk,1).varname{1};
            Rvar1_str = [Rvar1_str,',',loc_kk];
        end
    end
end
Rvar1_str = [Rvar1_str,')'];



% % % Finally, combine to construct the full term
term_str = [int_trm,C_trm,D_trm,Rstate_trm,Rvar1_str,dtheta_trm];

end