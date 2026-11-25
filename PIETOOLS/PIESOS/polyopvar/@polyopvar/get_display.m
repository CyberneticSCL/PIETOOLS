function row_disp = get_display(P,expand)
% ROW_DISP = GET_DISPLAY(P,EXPAND) returns, for each row-block i of P.C,
% the terms of the polynomial
%   C_i1 x + ... + C_ij (x^{od} o y^{or}),
% with the operators C_ij acting on the monomials displayed by name, or
% as the actual operators if they are small enough, for the 'display'
% method of 'polyopvar'.
%
% INPUTS
% - P:          'polyopvar' object with nr row-blocks;
% - expand:     (optional) logical, false to always display the operators
%               as C_ij. Defaults to true;
%
% OUTPUTS
% - row_disp:   nr x 1 'struct' array, with for row-block i the fields
%     terms:    1 x n 'cell' of char, the terms of the polynomial without
%               their sign, in the order of the monomials: C_ij followed
%               by the monomial, or the terms of the display of C_ij
%               acting on the monomial;
%     signs:    1 x n array of 1 and -1, the sign of each term;
%     expanded: logical, true if the operators are displayed as the
%               actual operators, false if as C_ij;
%     states:   'cellstr', the display names of the states in the
%               monomials on which a nonempty block of the row acts;
%     out_args: 'cellstr', the display names of the output variables,
%               with a subscript for the copy of a variable if it has
%               several copies;
%
% NOTES
% - The operators of a row are displayed as the actual operators only if
%   all of them can be displayed by 'get_display' of 'tensopvar' or
%   'intop' (which requires them to be scalar-valued, or identity
%   operators, displayed as the monomial), and the row has at most 6
%   terms, each of at most 100 characters. Otherwise, all operators of
%   the row are displayed as C_ij.
% - A monomial is displayed as the tensor product of its factors, in the
%   order of P.Z, with d consecutive copies of a state x as x^{od}, and a
%   factor that is a multiplier from R as 1. A functional ('intop') acts
%   on the monomial as a function of the variables of all its factors,
%   e.g. (x o y)(s_1,s_2).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - get_display
%
% Copyright (C) 2026 PIETOOLS Team
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 11/24/2026: Initial coding;

% Default settings
max_terms = 6;          % maximal number of monomials and terms for which the operators
                        % are displayed as the actual operators

if nargin<2
    expand = true;
end
if ~P.chkval
    P = initialize(P);
end
state_names = cellfun(@varname2display,reshape(P.varname,1,[]),'UniformOutput',false);
[nr,nZ] = size(P.C.ops);

otimes_sym = [' ',varname2display('otimes'),' '];
[~,sub_dig] = varname2display('0123456789');                 % subscript digits
sub_str = @(n) sub_dig(sprintf('%d',n)-'0'+1);

% Display each monomial as the tensor product of its factors: consecutive
% copies of a state x as a tensor power x^{od}, and a multiplier from R as 1
nfac = full(num_factors(P));
mono_strs = cell(1,nZ);
mono_states = cell(1,nZ);                % states in each monomial
mono_fctrs = cell(1,nZ);                 % the state of each factor of each monomial
for j=1:nZ
    fctr_states_j = full(P.Z(j,1:min(nfac(j),size(P.Z,2))));
    fctr_states_j = [fctr_states_j,zeros(1,nfac(j)-numel(fctr_states_j))];           % zeros for factors beyond Z
    mono_fctrs{j} = fctr_states_j;
    mono_states{j} = unique(fctr_states_j(fctr_states_j>0));
    fctr_groups = cell(1,0);
    k = 1;
    while k<=numel(fctr_states_j)
        % Number of consecutive copies of the same state
        r = 1;
        while fctr_states_j(k)>0 && k+r<=numel(fctr_states_j) && fctr_states_j(k+r)==fctr_states_j(k)
            r = r+1;
        end
        if fctr_states_j(k)==0
            fctr_groups{end+1} = '1';                    %#ok<AGROW>
        else
            fctr_groups{end+1} = state_names{fctr_states_j(k)};           %#ok<AGROW>
            if r>1
                [~,~,pow_str] = varname2display(sprintf('otimes%d',r));      % tensor power
                fctr_groups{end} = [fctr_groups{end},pow_str];
            end
        end
        k = k+r;
    end
    % Tensor product of the groups, in brackets if there are several
    mono_strs{j} = strjoin(fctr_groups,otimes_sym);
    if numel(fctr_groups)>1
        mono_strs{j} = ['(',mono_strs{j},')'];
    end
end

% Per row: the states in its monomials, the output variables, and the terms
nvars_out = numel(P.vars.out);
row_disp = repmat(struct('terms',{cell(1,0)},'signs',zeros(1,0),'expanded',false,'states',{cell(1,0)},'out_args',{cell(1,0)}),nr,1);
for i=1:nr
    mono_idcs = find(~cellfun('isempty',P.C.ops(i,:)));
    state_idcs = unique([mono_states{mono_idcs}]);
    row_disp(i).states = state_names(state_idcs);
    % Output variables, numbered by copy if a variable has several copies
    out_args = cell(1,0);
    for k=1:nvars_out
        ncopies = P.depmat.out(i,k);
        for c=1:ncopies
            if ncopies>1
                var_disp = varname2display(P.vars.out{k},c);
            else
                var_disp = varname2display(P.vars.out{k});
            end
            out_args{end+1} = var_disp;                       %#ok<AGROW>
        end
    end
    row_disp(i).out_args = out_args;

    expanded = expand && numel(mono_idcs)<=max_terms;
    terms = cell(1,0);      signs = zeros(1,0);
    for ii=1:numel(mono_idcs)
        if ~expanded
            break
        end
        j = mono_idcs(ii);
        blk = P.C.ops{i,j};
        % The names of the functions on which the factors of the block act
        fctr_names = cell(1,nfac(j));
        for fctr_idx=1:nfac(j)
            if mono_fctrs{j}(fctr_idx)>0
                fctr_names{fctr_idx} = state_names{mono_fctrs{j}(fctr_idx)};
            else
                fctr_names{fctr_idx} = '1';
            end
        end
        if isa(blk,'tensopvar') && size(blk.ops,2)==nfac(j)
            blk_disp = get_display(blk,fctr_names);
            is_ok = blk_disp.expanded;
            % The display of an identity operator on several factors names the function it acts
            % on by the generic symbol x, not by the factors: use the monomial instead, with the
            % output variables, e.g. x^{o2}(s_1,s_2)
            if blk_disp.is_id && nfac(j)>1
                id_str = mono_strs{j};
                if ~isempty(blk_disp.out_args)
                    id_str = [id_str,'(',strjoin(blk_disp.out_args,','),')'];
                end
                blk_disp.terms = {id_str};
            end
        elseif isa(blk,'intop')
            % A functional acts on the monomial as a function of the variables
            % of all its factors, e.g. (x o y)(s_1,s_2)
            x_fctr_str = mono_strs{j};
            if nfac(j)>1 && x_fctr_str(1)~='('
                x_fctr_str = ['(',x_fctr_str,')'];
            end
            blk_disp = get_display(blk,x_fctr_str);
            is_ok = blk_disp.ok;
        elseif isempty(mono_states{j}) && isequal(size(blk),[1,1]) && (isnumeric(blk) || ((isa(blk,'dpvar') || isa(blk,'polynomial')) && isempty(blk.varname)))
            % A constant term c0 + sum_k c_k*d_k, with decision variables d_k
            if isa(blk,'dpvar')
                cf = full(blk.C(:));        nms = [{''};blk.dvarname(:)];
            else
                cf = full(double(blk));     nms = {''};
            end
            blk_disp = struct('terms',{cell(1,0)},'signs',zeros(1,0));
            for k=reshape(find(cf),1,[])
                if k>1 && abs(cf(k))==1
                    blk_disp.terms{end+1} = nms{k};
                else
                    blk_disp.terms{end+1} = strtrim([num2str(abs(cf(k)),'%g'),' ',nms{k}]);
                end
                blk_disp.signs(end+1) = sign(cf(k));
            end
            is_ok = true;
        else
            is_ok = false;
        end
        if ~is_ok
            expanded = false;
            break
        end
        terms = [terms,blk_disp.terms];                   %#ok<AGROW>
        signs = [signs,blk_disp.signs];                   %#ok<AGROW>
    end
    % Too many terms, or a block not displayed explicitly: display the
    % operators by name, C_ij acting on monomial j
    if expanded && numel(terms)>max_terms
        expanded = false;
    end
    if ~expanded
        terms = cell(1,numel(mono_idcs));
        for ii=1:numel(mono_idcs)
            j = mono_idcs(ii);
            % Subscripts i (several rows) and j (several monomials), with a comma
            % if an index exceeds 9
            op_name = 'C';
            if nr>1
                op_name = [op_name,sub_str(i)];
            end
            if nr>1 && nZ>1 && max(i,j)>9
                op_name = [op_name,','];
            end
            if nZ>1
                op_name = [op_name,sub_str(j)];
            end
            terms{ii} = [op_name,' ',mono_strs{j}];
        end
        signs = ones(1,numel(terms));
    end
    % Store the terms of row i
    row_disp(i).terms = terms;
    row_disp(i).signs = signs;
    row_disp(i).expanded = expanded;
end

end
