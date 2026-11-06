classdef(InferiorClasses={?polynomial,?sopvar,?dpvar}) sdopvar              % MMP, 09/29/2026
    % Represents PI maps from L_2^p[S1,S3] to L_2^q[S2,S3]
    %
    % This defines PI operators from one L2 space to another.
    %
    %   Pop: L_2^p[S1,S3] to L_2^q[S2,S3]
    %
    % Elements of this class are intended to be included in a cdopvar container object.
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
    %         C(alpha) = A(alpha) + B(alpha)^T*d, where d is column vector
    %         of decision variables
    %
    % CLASS properties
    % - P.dims:   2x1 vector specifying the rows and columns of the operator i.e. [p,q];
    % - P.vars: a structure with fields 'in' and 'out'
    % - P.vars.in: a cell array specifying the unique names of the spatial  % MMP, 09/09/2026
    % variables in the input function space, stored in the canonical order  % MMP, 09/09/2026
    % [S_3,S_1]. S_1 and S_3 are each sorted, but the concatenation is not, % MMP, 09/09/2026
    % so do not assume vars.in is sorted; look up a variable by name        % MMP, 09/09/2026
    % - P.vars.out: a cell array specifying the unique names of the spatial % MMP, 09/09/2026
    % variables in the output function space, canonical order [S_2,S_3]     % MMP, 09/09/2026
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
    %
    % - P.params:   struct with fields 'A' and 'B', each a 3^n3 (stored as
    % 3x3...x3) cell array, n_3 the number of variables in S_3, with
    % params.A{alpha} the vec of the constant coefficient matrix of cell
    % alpha and params.B{alpha} the nd x numel(params.A{alpha}) sparse
    % matrix of the decision variable coefficients;
    % - P.Zd:   nd x 1 cell array of the names of the decision variables d.
    %
    % Each kernel of the PI operator is represented as
    % K_alpha =
    % (I_dims(1)\otimes ZL(s2,s3)')
    %          unvec(params.A(alpha) + params.B(alpha)^T*Zd)
    %                (I_dims(2)\otimes ZR(t1,t3))
    % where K_alpha represents kernel for multiplier/integral operator
    % specified by alpha index
    % params.A has column vectors
    % Zd is a vector of decision variables, sorted (all monomial degrees 1)
    % ZL, ZR, vars, dom, dims, are same as in sopvar
    %
    % CONSTRUCTOR P = sdopvar(params,vars,Zd,ZL,ZR,dom,dims)
    %
    % INPUTS
    % - params: struct with fields 'A' and 'B', each a 3^n3 cell array, one
    %           cell per multi-index alpha in {1,2,3}^n3; params.A{i} is the
    %           vectorized constant part of the coefficient, params.B{i} the
    %           nd x numel(params.A{i}) matrix of decision variable
    %           coefficients (see CLASS properties);
    % - vars:   struct with fields 'in' and 'out', cell arrays of the names of
    %           the spatial variables of the input and output spaces;
    % - Zd:     list of decision variables, with nd entries;
    % - ZL, ZR: cell arrays of column vectors of monomial exponents in
    %           vars.out and vars.in, respectively;
    % - dom:    struct with fields 'in' and 'out', each an n x 2 array with the
    %           interval of each variable in vars.in and vars.out;
    % - dims:   1 x 2 vector, the size [p,q] of the operator;
    %
    % OUTPUTS
    % - P:      'sdopvar' object, with variables put in the canonical order
    %           and multipliers in canonical form (with a warning if rewritten);
    %
    %
    % CANONICAL MULTIPLIER FORM                                             % MMP, 08/29/2026
    % A parameter with alpha(k)==1 carries a factor delta(s_k-s_k'), which  % MMP, 08/29/2026
    % identifies s_k with s_k'. A coefficient at left degree a and right    % MMP, 08/29/2026
    % degree b in direction k then contributes the same monomial s_k^(a+b)  % MMP, 08/29/2026
    % as one at left degree a+b and right degree 0, so the stored           % MMP, 08/29/2026
    % coefficients are not determined by the operator. Objects of this class % MMP, 08/29/2026
    % therefore satisfy                                                     % MMP, 08/29/2026
    %                                                                       % MMP, 08/29/2026
    %   alpha(k)==1  ==>  the parameter has no content outside the columns  % MMP, 08/29/2026
    %                     whose ZR monomial has degree 0 in direction k,    % MMP, 08/29/2026
    %                                                                       % MMP, 08/29/2026
    % that is, a multiplier direction carries no dummy variable degree. The % MMP, 08/29/2026
    % constructor enforces this by rewriting anything that does not already % MMP, 08/29/2026
    % satisfy it, enlarging ZL where the folded degrees require it; see     % MMP, 08/29/2026
    % 'canonicalize_multiplier'. The representation is then unique, which is % MMP, 08/29/2026
    % what lets 'eq' compare coefficients and 'lpi_eq_sdopvar' constrain    % MMP, 08/29/2026
    % them to zero. Assigning to the properties directly bypasses the       % MMP, 08/29/2026
    % constructor and so bypasses the invariant; build a new object instead. % MMP, 08/29/2026
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - sdopvar
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
    % MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
    %                  mdopvar -> cdopvar, with every file and function named
    %                  after them. Mechanical rename, no functional change. The
    %                  mention above names the container that holds sdopvar
    %                  blocks, which is cdopvar: copvar admits sopvar blocks
    %                  only.
    % MMP, 09/29/2026: ?dpvar added to InferiorClasses, so that a dpvar
    %                  operand (gam*P, P - gam, [gam, P]) dispatches to this
    %                  class's methods, which route it to 'dpvar_op_copvar'
    %                  (legacy @dopvar semantics); before, the dpvar's own
    %                  method was called and failed. The classdef line
    %                  before is the (was) line at the end of this log.
    % MMP, 09/29/2026: Dead files of this class (no live caller), all by
    %                  other developers, moved to the dead-code folder (now
    %                  cleaning/dead_code/, whose README lists origin and
    %                  author): mtimes_AT, plus_decparam, plus_decparam_batch,
    %                  minus_decparam, lrmultiply_batch, rand_sdopvar,
    %                  test_script, private/MatrixMultiply and
    %                  private/CombineDecisionBasis.
    % MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): a blank line
    %                  that cut the help text short is now a comment line
    %                  (params and Zd described there), the 09/29/2026
    %                  entries condensed, the (was) classdef line moved to
    %                  the end of the log, constructor comments tightened;
    %                  code unchanged.
    % DJ, 11/06/2026: Add the dependent property 'kernels'.
% classdef(InferiorClasses={?polynomial,?sopvar}) sdopvar                   % MMP, 09/07/2026 % MMP, 09/29/2026 (was)
    %

    properties
        vars = struct('in',{},'out',{});
        dom = struct('in',zeros(0,2),'out',zeros(0,2));
        dims = [1,1];
        ZL = {};
        ZR = {};
        Zd = {};
        params = struct('A', {0},'B', {});  % represents C(d) = unvec(A + B'*d)
        % for an m-by-n decision opvar, A is size m*n-by-1, B is size nd-by-m*n.
    end
    properties(Dependent,Hidden)
        kernels;    % cell of the kernels defined by params, see get.kernels
    end

    methods
        function P = sdopvar(params,vars,Zd,ZL,ZR,dom,dims)
            % Canonical variable order first, vars.out = [S2,S3] and        % MMP, 09/09/2026
            % vars.in = [S3,S1]: ZL{i} is bound to vars.out(i), so a        % MMP, 09/09/2026
            % reordering permutes the vec index of every coefficient; the   % MMP, 09/09/2026
            % decision variables (rows of B) are carried through, so this   % MMP, 09/09/2026
            % does not scale with their number. A pass-through variable     % MMP, 09/09/2026
            % given two domains is rejected there.                          % MMP, 09/09/2026
            [vars,dom,ZL,ZR,rowIdx,colIdx,reordered] = ...                  % MMP, 09/09/2026
                canonical_var_order(vars,dom,ZL,ZR,dims);                   % MMP, 09/09/2026
            if reordered                                                    % MMP, 09/09/2026
                nrow = dims(1)*prod([cellfun(@numel,ZL),1]);                % MMP, 09/09/2026
                vecIdx = reshape(rowIdx(:)+(colIdx(:).'-1)*nrow,[],1);      % MMP, 09/09/2026
                for ii = 1:numel(params.A)                                  % MMP, 09/09/2026
                    % [] and scalar 0 (zero-block shorthand) cannot be      % MMP, 09/09/2026
                    % indexed, and permuting a zero block is a no-op.       % MMP, 09/09/2026
                    if numel(params.A{ii})>1                                % MMP, 09/09/2026
                        params.A{ii} = params.A{ii}(vecIdx);                % MMP, 09/09/2026
                    end                                                     % MMP, 09/09/2026
                    if size(params.B{ii},2)>1                               % MMP, 09/09/2026
                        params.B{ii} = params.B{ii}(:,vecIdx);              % MMP, 09/09/2026
                    end                                                     % MMP, 09/09/2026
                end                                                         % MMP, 09/09/2026
            end                                                             % MMP, 09/09/2026
            % Then the canonical multiplier form, the class invariant of    % MMP, 08/29/2026
            % the help text; the check reads only the forbidden positions,  % MMP, 08/29/2026
            % so a canonical object costs numel(C) per parameter and nothing % MMP, 08/29/2026
            % in the decision variables.                                    % MMP, 08/29/2026
            [params,ZL,ZR,was_rewritten] = ...                              % MMP, 08/29/2026
                canonicalize_multiplier(params,vars,ZL,ZR,dims);            % MMP, 08/29/2026
            if was_rewritten                                                % MMP, 08/29/2026
                warning('sdopvar:noncanonicalMultiplier',...                % MMP, 08/29/2026
                    ['A multiplier parameter depended on a dummy variable, so it was rewritten in '...% MMP, 08/29/2026
                     'canonical form and ZL was enlarged to hold the folded degrees. The operator '...% MMP, 08/29/2026
                     'is unchanged. Build multiplier parameters with the degree on ZL to avoid this.']);% MMP, 08/29/2026
            end                                                             % MMP, 08/29/2026
            % Store the canonical data
            P.params = params;
            P.vars = vars;
            P.dom = dom;
            P.dims = dims;
            P.ZR = ZR; 
            P.ZL = ZL;
            P.Zd = Zd;
        end

        function K = get.kernels(P)
            % K = GET.KERNELS(P) returns the kernels defining the operator P,
            % as a cell of the same size as P.params.A, with element K{i}
            % the dims(1) x dims(2) kernel K_i associated with the parameter
            % P.params.A{i}, P.params.B{i}, as a 'dpvar' object in the
            % variables vars.out, in the dummy variables, named as the
            % variables in vars.in with _dum appended, and in the decision
            % variables Zd. See 'display' for the integrals corresponding to
            % each kernel.
            K = cell(size(P.params.A));
            for k=1:numel(K)
                K{k} = params2kernel(P,k);
            end
        end
    end
end
