function b = subsasgn(a,L,RHS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b = subsasgn(a,L,RHS) assigns values RHS to field/slice L of an nopvar 
% object a
% 
% INPUT
% a:    nopvar class object
% L:    a struct specifying the component to be adjusted/assigned
% RHS:  a value to be assigned to the component
%
% OUTPUT
% b:    nopvar object with desired component value
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2025 PIETOOLS Team
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
% Initial coding AT, 01/30/2026;

switch L(1).type
    case '.'
        b = builtin('subsasgn',a,L,RHS); 
    case '()'

        % Check if the indices are supported.
        if isscalar(L(1).subs)
            error('Type ''()'' subsasgn using linear indices is not supported.')
        elseif length(L(1).subs)>=3
            error('opvar objects take at most two subscripts.')
        end
        % Make sure we can work with the inputs.
        if any(any(a.dom~=RHS.dom))|| any(~strcmp(a.vars(:, 1).varname,RHS.vars(:, 1).varname)) || any(~strcmp(a.vars(:, 2).varname,RHS.vars(:, 2).varname))
            error('The spatial domain and variables of the right-hand side operator should match those of the old operator.')
        end
        indr = L(1).subs{1}; %output
        indc = L(1).subs{2}; %input
        dim_old=a.dim;

        nr1_old=dim_old(1); % output dimension
        nr2_old=dim_old(2); % input dimension

        % Allow indices to be specified as e.g. (i,:);
        if strcmp(indr,':') && strcmp(indc,':')
            out = obj;
            return
        end
        if strcmp(indr,':')
            indr = (1:nr1_old);
        end
        if strcmp(indc,':')
            indc = (1:nr2_old);
        end


        if size(indr, 1) == 1
            indr = indr';
        end
        if size(indc, 1) == 1
            indc = indc';
        end

        % Check that the output dimensions match
        if (length(indr)~=RHS.dim(1))
            error('Output dimensions of ndopvar/nopvar objects do not match')
        end
        if (length(indc)~=RHS.dim(2))
            error('Input dimensions of ndopvar/nopvar objects do not match')
        end
        % Check if the RHS is of appropriate type. 
        if ~isa(RHS, 'nopvar') && ~isa(RHS, 'ndopvar')
            error('Currently implemented only for nopvar/ndopvar')
        end
        % if the degrees are different convert to the same
        if any(a.deg(:)~=RHS.deg(:))
            max_degree = max(a.deg, RHS.deg);
            a = change_degree(a, max_degree);
            RHS = change_degree(RHS, max_degree);
        end


        % if operators have different decvar, convert to the same decvar 
        if isa(RHS, 'ndopvar') % if RHS is nopvar, convert it to ndopvar
            RHS_dvarname = RHS.dvarname;
            a = change_dec_var(a, RHS_dvarname); % convert to the same dec var
            % a(L(1).subs{1}, L(1).subs{2}) = RHS; % use subsasign for ndopvar

            b = a;
            b(L(1).subs{1}, L(1).subs{2}) = RHS;
            return
        end

                % Construct the sliced opvar
        b = nopvar();
        b.dom = a.dom;
        b.deg = a.deg;
        b.vars = a.vars; 
        N = size(a.dom,1);
        b.C = cell([3*ones(1,N),1]);

        dec_var_NEW = []; % nopvar no dec var
        binStr = dec2base(0:(3^N-1), 3);% create indexing for C{i}
        
        size_of_monom = prod(a.deg+1); % size of monomimials 
        
        Loc_L_R = size_of_monom*(indr - 1) + 1; % starting location of C rows
        Loc_L_R= [kron(Loc_L_R, ones(size_of_monom, 1)), kron( ones(length(Loc_L_R), 1), (0:(size_of_monom-1))' )];
        Loc_L_R = sum(Loc_L_R, 2);
        for iter_idx = 1:size(binStr, 1)
        
            substr = binStr(iter_idx, :);
            Qindx = str2num(substr')'+1; % array of .C indeces
            Qindx_cell = num2cell(Qindx); % cell array of .C indeces
        
            index_for_monom_in_theta = Qindx == 1; % include monomials for theta
        
            subdegree1 = a.deg;
            subdegree1(index_for_monom_in_theta) = 0;
        
            size_of_monom_c = prod(subdegree1 + 1);
            Loc_L_C = size_of_monom_c*(indc - 1) + 1; % starting location of C columns
            Loc_L_C= [kron(Loc_L_C, ones(size_of_monom_c, 1)), kron( ones(length(Loc_L_C), 1), (0:(size_of_monom_c-1))' )];
            Loc_L_C = sum(Loc_L_C, 2);
        
            % change Cop values
            full_matrix = a.C{Qindx_cell{:}};
            full_matrix(Loc_L_R, Loc_L_C) = RHS.C{Qindx_cell{:}};
            b.C{Qindx_cell{:}} = full_matrix;
        end

        % Extract the dimensions of the operators a and RHS.
        % end
    case '{}'
        error('{}- like subsassign is not supported for nopvar objects.');
end

end