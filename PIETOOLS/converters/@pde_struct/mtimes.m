function PDE_out = mtimes(C,PDE_in)
% PDE_OUT = MTIMES(C,PDE_IN) declares a new PDE term corresponding to
% the left-product of a PDE term "PDE_in" with coefficients "C".
%
% INPUT
% - C:          mxn object of type 'polynomial' or 'double', represting a
%               matrix of coefficients with which to multiply terms in the
%               PDE.
% - PDE_in:     'pde_struct' object corresponding to either a single PDE
%               variable (state, input, or output) of size n or
%               corresponding to one or several terms in a PDE, belonging
%               to n equations.
%
% OUTPUT
% - PDE_out:    'pde_struct' object representing one or several terms in a
%               PDE, corresponding to the product C*PDE_in;
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

% % % Process the inputs

% % Check that the first input makes sense.
if ~isa(C,'double') && ~isa(C,'polynomial') && ~isa(C,'opvar')
    error("The first factor in the product should be of type 'double' or 'polynomial'.")
elseif ~isa(C,'opvar')
    C = polynomial(C);
end

% % Check that the second input makes sense.
% Make sure the input is a pde_struct object.
if ~isa(PDE_in,'pde_struct')
    error("The second factor in the product should correspond to a single term in the PDE.")
end
% Make sure the input corresponds to a single state/input/output
% or set of terms.
[is_pde_var_in,obj] = is_pde_var(PDE_in);
if ~is_pde_var_in && ~is_pde_term(PDE_in)
    error("The second factor in the product should correspond to a single term in the PDE.")
end

% Initialize the output PDE.
PDE_out = PDE_in;


% % % Deal with the case that the coefficients are specified as 'opvar'
% % % object, requiring (partial) integrals to be taken.
if isa(C,'opvar')
    % % Convert separate PDE variables to free terms.
    if is_pde_var_in
        PDE_in = var2term(PDE_in);
    end
    % % Check that the size makes sense.
    n_eqs = size(PDE_in,'free','vec_size_tot');
    if sum(C.dim(:,2))~=n_eqs
        error("Column dimension of the factor should match row dimension of PDE object.")
    elseif n_eqs==0
        PDE_out = PDE_in;
        return
    elseif sum(C.dim(:,1))==0
        % Empty output.
        PDE_out = PDE_in;
        PDE_out.free = cell(0,1);
        return
    end
    
    % % Now, the opvar object has a structure
    % %     C = [P, Q1; Q2, R]
    % % Extract information from the object.
    n_eqs_1 = C.dim(1,2);       n_eqs_2 = C.dim(2,2);
    n_eqs_out1 = C.dim(1,1);    n_eqs_out2 = C.dim(2,1);
    Cvar1 = C.var1;      Cvar2 = C.var2;
    Cdom = C.I;
    if n_eqs_1~=0
        % % Determine which equations to multiply with [P;Q2], and which to
        % % multiply with [Q1;R]
        eq_num = 0;     eq_tot = 0;
        while eq_tot<n_eqs_1
            % Keep adding equations until we have "n_eqs_1" many.
            eq_num = eq_num + 1;
            eq_tot = eq_tot + PDE_in.free{eq_num}.size;
        end
        if eq_tot~=n_eqs_1
            error("Multiplication with opvar is not supported: some vector-valued state or input components appear to have both finite- and infinite-dimensional elements")
        end
        % Extract equations to multiply with [P;Q2].
        PDE_in1 = PDE_in;
        PDE_in1.free = PDE_in1.free(1:eq_num);
        if n_eqs_2~=0
            % Extract equations to multiply with [Q1;R];
            PDE_in2 = PDE_in;
            PDE_in2.free = PDE_in2.free(eq_num+1:end);
        end
    else
        % No equations to multiply with [P;Q2].
        PDE_in2 = PDE_in;
    end

    % % Apply the finite-dimensional components of C to the PDE;
    if n_eqs_1~=0
        % Apply P to the PDE;
        PDE_out11 = C.P*PDE_in1;
        % Apply Q2 to the PDE;
        PDE_out21 = C.Q2*PDE_in1;
    else
        PDE_out11 = 0;
        PDE_out21 = 0;
    end
    % % Apply the infinite-dimensional components of C to the PDE;
    if n_eqs_2~=0
        % Apply Q1 to the PDE;
        PDE_out12 = int(C.Q1*PDE_in2,Cvar1,Cdom);
        % % Apply R to the PDE
        % First apply multiplier operator.
        PDE_out22 = C.R.R0*PDE_in2;
        % Check if any integration is performed.
        R1 = cleanpoly(polynomial(C.R.R1),1e-12);
        R2 = cleanpoly(polynomial(C.R.R2),1e-12);
        if any(any(abs(R1.C)))
            % Partial integral a^s is taken.
            if any(any(abs(R2.C)))
                % Partial integral s^b is also taken
                % --> re-organize as full integral and Volterra operator.
                PDE_out22 = PDE_out22 +int(R2*PDE_in2,Cvar1,Cdom);
                R1 = cleanpoly(R1-R2,1e-12);
                if any(any(abs(R1.C)))
                    PDE_subs = subs(PDE_in2,Cvar1,Cvar2);
                    PDE_out22 = PDE_out22 +int(R1*PDE_subs,Cvar2,[Cdom(1),Cvar1]);
                end
            else
                % Only partial integral a^s is taken;
                PDE_subs = subs(PDE_in2,Cvar1,Cvar2);
                PDE_out22 = PDE_out22 +int(R1*PDE_subs,Cvar2,[Cdom(1),Cvar1]);
            end
        elseif any(any(abs(R2.C)))
            % Only partial integral s^b is taken;
            PDE_subs = subs(PDE_in2,Cvar1,Cvar2);
            PDE_out22 = PDE_out22 +int(R2*PDE_subs,Cvar2,[Cvar1,Cdom(2)]);
        end
    else
        PDE_out12 = 0;
        PDE_out22 = 0;
    end
    % % Add the terms, and concatenate
    if n_eqs_out1>0 && n_eqs_out2>0
        PDE_out = [PDE_out11+PDE_out12; PDE_out21+PDE_out22];
    elseif n_eqs_out1>0
        PDE_out = PDE_out11+PDE_out12;
    elseif n_eqs_out2>0
        PDE_out = PDE_out21+PDE_out22;
    end
    return
end

% % % Deal with the case that the input PDE corresponds to just a PDE 
% % % variable (state, input, output)
% % % --> Initialize a PDE term which to multiply with C
if is_pde_var_in
    % Check that the number of columns in C matches the number of rows of
    % the PDE variable.
    if prod(size(C))~=1 && size(C,2)~=PDE_in.([obj,'_tab'])(1,2)
        error("Column dimension of the factor should match row dimension of PDE object.")
    end
    
    % Extract spatial variables on which the PDE object depends
    vars = PDE_in.(obj){1}.vars;    nvars = size(vars,1);

    % Initialize a PDE term with no derivative or integral
    PDE_out.free{1}.term{1}.(obj) = 1; % or should it be PDE_in.([obj,'_tab'])(1,1);
    if strcmp(obj,'x')
        PDE_out.free{1}.term{1}.loc = vars';
        PDE_out.free{1}.term{1}.D = zeros(1,nvars);
    end
    PDE_out.free{1}.term{1}.I = cell(nvars,1);

    % Set the coefficients and size of the term.
    if all(size(C)==1)
        % C acts as scalar --> size of term matches that of "obj".
        PDE_out.free{1}.term{1}.C = C*eye(PDE_in.([obj,'_tab'])(1,2));
        PDE_out.free{1}.size = PDE_in.([obj,'_tab'])(1,2);
    else
        % C determines the size.
        PDE_out.free{1}.term{1}.C = C;
        PDE_out.free{1}.size = size(C,1);        
    end

    % Keep track of which variables the newly started equation depends on.
    new_varnames = unique([vars.varname;C.varname]);
    PDE_out.free{1}.vars = polynomial(new_varnames);
    %PDE_out.free{1}.vars = addvars(vars,C);
    %PDE_out.free{1}.term{1}.vars = addvars(vars,C);
    %PDE_out.free{1}.vars = PDE_out.free{1}.term{1}.vars;

    % Set the variables and size of the equation to which the term is to be
    % added.    
    %PDE_out.free{1}.size = PDE_out.free{1}.term{1}.size;
    return
end

% % % At this point, assume the PDE correponds to one or multiple terms.

% % Check that the number of equations matches the column dimension of the 
% % coefficients.
eq_size_list = size(PDE_in,'free','vec_size');
if size(C,1)==1 && size(C,2)==1
    % Multiplication with a scalar: build an identity matrix;
    C = C*eye(sum(eq_size_list));
    % Convert to a cell, with each row corresponding to an equation in the
    % output structure, and each column to an equation in the input.
    C_cell = mat2cell_poly(C,eq_size_list,eq_size_list);
    % Keep track of which elements of C_cell are nonzero.
    include_eqs = eye(numel(eq_size_list));
elseif sum(eq_size_list)==size(C,2)
    % % Convert C to a cell, with each column corresponding to an input
    % % equation, and each column to an output equation.
    if size(C,1)==sum(eq_size_list)
        % Assume the number of input equations matches number of output
        % equations.
        C_cell = mat2cell_poly(C,eq_size_list,eq_size_list);
    else
        % Assume only a single output equation is desired.
        C_cell = mat2cell_poly(C,size(C,1),eq_size_list);
    end
    % Assume all elements of C_cell are nonzero.
    include_eqs = ones(size(C_cell));
elseif numel(eq_size_list)==size(C,2)
    % The size of C matches the number of equations.
    % Presume each element of C corresponds to a scalar factor which to
    % multiply the equation with.
    if length(unique(eq_size_list))~=1
        error("For addition of equations with a scalar factor, the size of the terms in the equation should match.")
    end
    if size(C,1)==1 && length(unique(eq_size_list))==1
        % Single output equation, assume the same size as all inputs.
        eq_size_cell_out = num2cell(sum(eq_size_list(1))*ones(1,size(C,2)));
    elseif size(C,1)==size(C,2)
        % Same amount of input and output equations.
        eq_size_cell_out = num2cell(repmat(eq_size_list,1,size(C,2)));
    elseif all(eq_size_list==1)
        % For all scalar equations, assume also scalar output equations.
        eq_size_cell_out = num2cell(ones(size(C)));
    else
        error("Multiplication is not supported: number of output equations is ambiguous.")
    end
    eq_size_cell_in = num2cell(repmat(eq_size_list',size(C,1),1));
    C_cell = mat2cell_poly(C);
    C_cell = cellfun(@(a,nr,nc) a*eye(nr,nc),C_cell,eq_size_cell_out,eq_size_cell_in,'UniformOutput',false);
    %C = cell2mat(C_cell);
    % Keep track of which element of C_cell are nonzero
    include_eqs = ~isequal(C,0);
else
    error("Number of columns of the coefficients should match number of rows of equations.")
end

% % % Check that the number of output equations matches the number of input
% % % equations.
% if numel(eq_size_list)>1 && size(C,1)~=size(C,2)
%     error("For multiplication with multiple equations, the number of rows of C must match the number of columns, to avoid ambiguity.")
% end

% % % For now, assume only a single equation
% if numel(PDE_in.free)>1 && any(size(C)>1)
%     error('Multiplication with multiple equations is currently not supported.')
% end

% % Loop over all the output equations, adding terms from each of the
% % input equations, scaled with the appropriate factor.
PDE_out.free = cell(size(C_cell,1),1);
for ii=1:numel(PDE_out.free)
    % Set the size of the output equation.
    PDE_out.free{ii}.size = size(C_cell{ii,1},1);
    PDE_out.free{ii}.eq = {};

    % Keep track of which variables the term depends on.
    vars_ii = polynomial(zeros(0,1));

    % Set the terms of the output equations.
    use_eqs = find(include_eqs(ii,:));
    term_num = 1;
    for kk=use_eqs
        % % Add terms from equation kk.
        Cik = C_cell{ii,kk};        eq_kk = PDE_in.free{kk};
        if all(all(isequal(Cik,0))) || ~isfield(eq_kk,'term') || isempty(eq_kk.term)
            % If all coefficients are zero, this equation contributes no
            % terms.
            continue
        end
        for jj=1:numel(eq_kk.term)
            % Add the term to the output equation.
            PDE_out.free{ii}.term{term_num} = eq_kk.term{jj};
            % Check first if the term corresponds to an output object, or a
            % temporal derivative of a state.
            if is_LHS_term(eq_kk.term{jj})
                % The term corresponds to a temporal derivative of a state, or to
                % an output signal. Only multiplication with +1 or -1 is supported.
                if isdouble(Cik) && isequal(double(Cik),eye(size(Cik,2)))
                    PDE_out.free{ii}.term{term_num}.C = 1;
                elseif isdouble(Cik) && isequal(double(Cik),-eye(size(Cik,2)))
                    PDE_out.free{ii}.term{term_num}.C = -1;
                else
                    error("Multiplication with output signals or temporal derivatives of state variables is not supported.")
                end
            else
                % Add a standard term, with coefficients multiplied with
                % Cik;
                if isfield(eq_kk.term{jj},'C')
                    PDE_out.free{ii}.term{term_num}.C = Cik*eq_kk.term{jj}.C;
                else
                    PDE_out.free{ii}.term{term_num}.C = Cik;
                end
            end
            % Move to the next term;
            term_num = term_num+1;
        end
        % Add variables from equation kk to the list;
        varnames_kk = unique([eq_kk.vars.varname;Cik.varname]);
        varnames_ii = unique([vars_ii.varname;varnames_kk]);
        vars_ii = polynomial(varnames_ii);
        % vars_kk = addvars(eq_kk.vars,Cik);
        % vars_ii = [vars_ii;vars_kk];
        % vars_ii = polynomial(vars_ii.varname); % remove duplicate variables.
    end

    % Set the variables on which the output equation depends.
    PDE_out.free{ii}.vars = vars_ii;
end

% for ii=1:numel(PDE_out.free)
% 
%     % Check that the size of the coefficients is appropriate.
%     if all(size(C)==1)
%         % C is scalar-values --> adjust to appropriate size.
%         Cik = C*eye(PDE_in.free{ii}.size);
%     elseif size(C,2)~=PDE_in.free{ii}.size
%         error("Column dimension of the factor should match row dimension of PDE object.")   
%     else
%         Cik = C;
%     end
%     PDE_out.free{ii}.size = size(Cik,1);
% 
%     % % Set the new variables of the equation.
%     %PDE_out.free{1}.vars = addvars(PDE_in.free{1}.vars,C);
% 
% 
%     % For multiplication with 0, just get rid of all the terms
%     if all(all(isequal(Cik,0)))
%         PDE_out.free{ii}.term = {};
%         continue
%     end
%     % % Otherwise, loop over all terms, multiplying the coefficients with C
%     for jj=1:numel(PDE_in.free{ii}.term)
%         % Check first if the term corresponds to an output object, or a
%         % temporal derivative of a state.
%         if is_LHS_term(PDE_in.free{ii}.term{jj})
%             % The term corresponds to a temporal derivative of a state, or to
%             % an output signal. Only multiplication with +1 or -1 is supported.
%             if isdouble(Cik) && isequal(double(Cik),eye(size(Cik,2)))
%                 PDE_out.free{ii}.term{jj}.C = 1;
%             elseif isdouble(Cik) && isequal(double(Cik),-eye(size(Cik,2)))
%                 PDE_out.free{ii}.term{jj}.C = -1;
%             else
%                 error("Multiplication with output signals or temporal derivatives of state variables is not supported.")
%             end
%             continue
%         end
%         PDE_out.free{ii}.term{jj}.C = Cik*PDE_in.free{ii}.term{jj}.C;
%         %PDE_out.free{1}.term{jj}.vars = addvars(PDE_in.free{1}.term{jj}.vars,C);
%         %PDE_out.free{1}.term{jj}.vars = addvars(PDE_in.free{1}.term{jj}.vars,C);
%     end
%     % Keep track of which variables the equation depends on
%     PDE_out.free{ii}.vars = addvars(PDE_in.free{ii}.vars,Cik);
% 
% end


% % % Check that the number of columns in C matches the number of equations.
% % Determine the number of equations.
% RHS_size = 0;
% for ii=1:numel(PDE_in.free)
%     RHS_size = RHS_size + PDE_in.free.size;
% end
% 
% if all(size(C)==1)
%     % C is scalar-valued --> adjust to appropriate size
%     C = C*eye(RHS_size);
% elseif size(C,2)~=RHS_size && size(C,2)==numel(PDE_in.free)
%     % Each column acts on one vector-sized equation --> adjust to
%     % appropriate size
%     if size(C,1)~=size(C,2)
%         error('For multiplication with multiple terms, the cooefficient matrix should be diagonal')
%     end
% 
% end




end

% function newvars = addvars(oldvars,C)
% % Add variables that appear in polynomial object C to the list of variables
% % oldvars.
% 
% if isdouble(C)
%     % C introduces no new variables.
%     newvars = oldvars;
% else
%     % Combine spatial variabels of the PDE variable and the factor C
%     % Has to be done somewhat akwardly, as the polynomial structure
%     % rearranges the variables alphabetically
%     nvars = size(oldvars,1);
%     vnames = cell(nvars+length(C.varname),1);
%     for jj=1:nvars
%         vnames{jj} = oldvars(jj,1).varname{1};
%     end
%     vnames(nvars+1:end) = C.varname;
%     vnames = unique(vnames,'stable');
%     newvars = polynomial(vnames);
% end
% 
% end

function P_cell = mat2cell_poly(P,nr_list,nc_list)
% Function mat2cell for 'polynomial' classs objects.

% If no numbers of rows and columns in each cell element are specified,
% assume each element contains a scalar element of the polynomial.
if nargin<=2
    nc_list = ones(size(P,2),1);
end
if nargin==1
    nr_list = ones(size(P,1),1);
end

% Set start and end indices for which rows/columns in P will end up in each
% element in P_cell;
nnr_list = cumsum([0;nr_list(:)]);
nnc_list = cumsum([0;nc_list(:)]);

% Build the output cell, filling it in one element at a time.
P_cell = cell(numel(nr_list),numel(nc_list));
for kk=1:numel(P_cell)
    [ii,jj] = ind2sub(size(P_cell),kk);
    r_idcs = nnr_list(ii)+1:nnr_list(ii+1);
    c_idcs = nnc_list(jj)+1:nnc_list(jj+1);
    P_cell{ii,jj} = P(r_idcs,c_idcs);
end

end