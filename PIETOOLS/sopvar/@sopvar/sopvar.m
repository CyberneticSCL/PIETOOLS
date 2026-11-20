classdef (InferiorClasses={?polynomial,?dpvar})sopvar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - sopvar class
    % Represents PI maps from L_2^p[S1,S3] to L_2^q[S2,S3]
    %
    % This defines PI operators from one L2 space to another.
    %
    %   Pop: L_2^p[S1,S3] to L_2^q[S2,S3]
    %
    %   The operator has the form
    %         y(S2,S3)=sum_alpha int_S1 int_S3dum  I_alpha(S_3-S_3dum) Z_d(S_2,S_3)
    %         C{\alpha} Z_d(S_3dum,S_1) x(S_3dum,S_1)
    %                         alpha \in \{0,1,-1\}^{n_3}
    %
    %   where S_1 is the variables in x(S1,S3) (input) not in S3 (common)
    %         S_2 is the variables in y(S2,S3) (output) not in S3 (common)
    %         S_3 is the variables common to x (input) and y (output)
    %         S_3dum are dummy versions of the variables in S_3
    %   and the cell index gamma in {1,2,3}^n_3 of params stores alpha with
    %   gamma_k = 1, 2, 3 for alpha_k = 0 (multiplier, delta(s_k-s_k')),
    %   1 (integral from the lower limit to s_k), -1 (from s_k to the upper
    %   limit), direction 1 varying fastest.
    %
    % CLASS properties
    % - P.params:   3^n3 (stored as 3x3...x3) Cell array of sparse matrices, where n_3 is the number of variables in S_3;
    % - P.dims:   1x2 vector specifying the rows and columns of the operator i.e. [q,p];
    % - P.vars: a structure with fields 'in' and 'out'
    % - P.vars.in: a cell array specifying the unique names of the spatial
    % variables in the input function space, stored in the canonical order
    % [S_3,S_1]. S_1 and S_3 are each sorted, but the concatenation is not,
    % so do not assume vars.in is sorted; look up a variable by name
    % - P.vars.out: a cell array specifying the unique names of the spatial
    % variables in the output function space, canonical order [S_2,S_3]
    % - P.dom.: struct object with fields 'in' and 'out'
    % - P.dom.in: ordered array (2 columns) specifying the spatial domain of each
    % variable in P.vars.in, so that dom(i,:) = [ai,bi] -> vars.in(i) \in [ai,bi]
    % - P.dom.out: same as P.dom.in but for out variables stored in
    % P.vars.out
    % - ZL: cell array of column vectors representing exponents of
    % monomials in P.vars.out. That is monomials are constructed by as
    % follows: ZL(P.vars.out) = P.vars.out(1).^ZL{1}\otimes ... P.vars.out(i).^ZL{i}... 
    % - ZR: cell array of column vectors representing exponents of
    % monomials in P.vars.in. That is monomials are constructed by as
    % follows: ZR(P.vars.in) = P.vars.in(1).^ZR{1}\otimes ... P.vars.in(i).^ZR{i}...
    % - P.kernels: (dependent, hidden) cell of the size of P.params, the
    % dims(1) x dims(2) kernel of each parameter as a 'polynomial' object.
    %
    % CONSTRUCTOR P = sopvar(params,vars,ZL,ZR,dom,dims)
    %
    % INPUT
    % - params: cell array with one cell per multi-index gamma in
    %           {1,2,3}^n3 (3 x ... x 3 with n3 directions; 3 x 1 for
    %           n3 = 1; a single cell for n3 = 0), cell gamma the sparse
    %           coefficient matrix C_gamma of the kernel of that cell, of
    %           size dims(1)*nL x dims(2)*nR with nL, nR the numbers of
    %           monomials of ZL and ZR (rows and columns are component
    %           outer, monomial inner); [] or scalar 0 for a zero cell;
    % - vars:   struct with fields 'in' and 'out', cellstr of the spatial
    %           variable names of the input and output spaces, in any
    %           order;
    % - ZL, ZR: cell arrays with one column of monomial exponents per
    %           variable of vars.out and of vars.in, in the same order;
    % - dom:    struct with fields 'in' and 'out', one row [a,b] per
    %           variable of vars.in and of vars.out;
    % - dims:   1 x 2 array [q,p], the numbers of output and input
    %           components.
    %
    % OUTPUT
    % - P:      'sopvar' object, with the variables in the canonical order
    %           above (the coefficients permuted to match) and every
    %           multiplier cell in canonical form; a warning is issued when
    %           a parameter had to be rewritten.
    %
    % CANONICAL MULTIPLIER FORM
    % A parameter P.params{alpha} with alpha(k)==1 carries a factor
    % delta(s_k-s_k'), which identifies s_k with s_k'. A coefficient at left
    % degree a and right degree b in direction k then contributes the same
    % monomial s_k^(a+b) as one at left degree a+b and right degree 0, so
    % the stored coefficients are not determined by the operator. Objects of
    % this class therefore satisfy
    %
    %   alpha(k)==1  ==>  P.params{alpha} has no content outside the columns
    %                     whose ZR monomial has degree 0 in direction k,
    %
    % that is, a multiplier direction carries no dummy variable degree. The
    % constructor enforces this by rewriting anything that does not already
    % satisfy it, enlarging ZL where the folded degrees require it. The
    % representation is then unique: equal operators have equal
    % coefficients. Assigning to P.params directly bypasses the constructor
    % and so bypasses the invariant; build a new object instead.
    %
    % For support, contact M. Peet, Arizona State University at mpeet@asu.edu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    % SS,MP, 01/15/2026: Initial coding
    properties
        vars = struct('in',{},'out',{});
        dom = struct('in',zeros(0,2),'out',zeros(0,2));
        dims = [1,1];
        ZL = {};
        ZR = {};
        params = {0};
    end
    properties(Dependent,Hidden)
        kernels;    % cell of the kernels defined by params, as 'polynomial' objects
    end
    properties(Access=protected)
        vars_S1;   % these variables are integrated out (full integral only)
        vars_S2;   % these variables are introduced (multiplier only)
        vars_S3;   % these variables pass through (3 PI operators)
        dom_1;   % domain of s1
        dom_2;   % domain of s2
        dom_3;   % domain of s3
    end

    methods
        function P = sopvar(params,vars,ZL,ZR,dom,dims)
            % Canonical variable order first, vars.out = [S2,S3] and
            % vars.in = [S3,S1]: ZL{i} is bound to vars.out(i), so a
            % reordering permutes the monomial index of every coefficient
            % (rowIdx, colIdx, lifted over the matrix dimension; empty when
            % nothing moves). A pass-through variable given two domains is
            % rejected there.
            [vars,dom,ZL,ZR,rowIdx,colIdx,reordered] = ...
                canonical_var_order(vars,dom,ZL,ZR,dims);
            if reordered
                for ii = 1:numel(params)
                    % [] and scalar 0 (zero-block shorthand) cannot be
                    % indexed, and permuting a zero block is a no-op.
                    if numel(params{ii})>1
                        params{ii} = params{ii}(rowIdx,colIdx);
                    end
                end
            end
            % Then the canonical multiplier form, the class invariant of
            % the help text: a non-canonical parameter is rewritten, a
            % canonical one is not touched.
            [params,ZL,ZR,was_rewritten] = ...
                canonicalize_multiplier(params,vars,ZL,ZR,dims);
            if was_rewritten
                warning('sopvar:noncanonicalMultiplier',...
                    ['A multiplier parameter depended on a dummy variable, so it was rewritten in '...
                     'canonical form and ZL was enlarged to hold the folded degrees. The operator '...
                     'is unchanged. Build multiplier parameters with the degree on ZL to avoid this.']);
            end
            % Store the canonical data
            P.params = params;
            P.vars = vars;
            P.dom = dom;
            P.dims = dims;
            P.ZR = ZR;
            P.ZL = ZL;

            % The variable partition S1 (input only), S2 (output only), S3
            % (shared), each sorted, with their domains
            P.vars_S1 = setdiff(vars.in,vars.out);
            P.vars_S2 = setdiff(vars.out,vars.in);
            P.vars_S3 = intersect(vars.in,vars.out);
            [~,idx] = ismember(P.vars_S1,vars.in);
            P.dom_1 = dom.in(idx,:);
            [~,idx] = ismember(P.vars_S2,vars.out);
            P.dom_2 = dom.out(idx,:);
            [~,idx] = ismember(P.vars_S3,vars.in);
            P.dom_3 = dom.in(idx,:);
        end

        function K = get.kernels(P)
            % K = GET.KERNELS(P) returns the kernels defining the operator P,
            % as a cell of the same size as P.params, with element K{i}
            % the dims(1) x dims(2) kernel K_i associated with P.params{i},
            % as a 'polynomial' object in the variables vars.out, and in
            % the dummy variables, named as the variables in vars.in with
            % _dum appended. Kernel K_i enters the operator with the
            % multiplier or integrals of the multi-index of cell i (class
            % help).
            %
            % INPUT
            % - P:  'sopvar' object;
            %
            % OUTPUT
            % - K:  cell of the size of P.params of 'polynomial' objects.
            K = cell(size(P.params));
            % One kernel per parameter cell
            for k=1:numel(K)
                K{k} = params2kernel(P,k);
            end
        end
    end
    methods(Static)
        [Minv, R1Minv, R2Minv, info] = inv_1D(M, K1, K2, interval, opts);
    end

end
