function [prog, V3, dV3] = V3_Liediff(prog, PIE, d, opdeg)
    % [prog,V3,dV3] = V3_Liediff(...) constructs the Lyapunov
    % functional V3 from Sec. 7.5 of the Automatica paper and returns its
    % Lie derivative along the polynomial PIE, as stated in Lem. 13.
    %
    % The PIE is assumed to have the form
    %
    %                 T*x_t = f(x) = sum_l C_l*x^(otimes l),
    %
    % and V3 is constructed in the pseudo-quadratic form
    %
    %       V3(x) = sum_i sum_j <U^i*x^i,Q_ij*(U*T)^j*x^j>_L2.
    %
    % The symmetry condition Qhat_ij^* = Qhat_ji from Lem. 13 is imposed
    % explicitly in this function. As this needs to be imposed block-wise,
    % this is not naturally enforced in LocalStability.m. The resulting 
    % Lie derivative is
    %
    %  dV3 = 2*sum_i sum_j sum_k sum_l
    %       <U^i*x^i,Q_ij*[ (U*T*x)^(k-1) otimes (U*C_l*x^l)
    %                              otimes (U*T*x)^(j-k) ]>_L2.
    %
    % where V3 and dV3 are built together so that they use and constrain the 
    % same free Gram blocks Q_ij.
    %
    % INPUTS
    % - prog    Current PIESOS program.
    % - PIE     Polynomial PIE structure with fields:
    %             PIE.T   inverse-map PI operator T;
    %             PIE.f   polyopvar right-hand side f(x);
    %             PIE.dom spatial domain [a,b].
    % - d       Maximum degree in Z_d. V3 has degree 2*d.
    % - opdeg   Degree of the spatial monomial basis defining U.
    %
    % OUTPUTS
    % - prog    Updated PIESOS program containing the free Q_ij variables
    %           and the adjoint-symmetry equality constraints from Lem. 13.
    % - V3      polyopvar representation of the Sec. 7.5 Lyapunov
    %           functional.
    % - dV      polyopvar representation of L_{T,C}V3 from Lem. 13.
    %
    % NOTES
    % - Lem. 13 is stated for a scalar fundamental state x in L2[a,b].
    %   This implementation therefore currently supports scalar PIEs.
    % - The factor of two uses the symmetry hypothesis in Lem. 13. This
    %   function imposes that hypothesis directly on every diagonal block
    %   and every distinct off-diagonal block pair.
    % - Constant terms in PIE.f are excluded because Lem. 13 assumes the
    %   polynomial PIE is a sum over degrees l >= 1.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - V3_Liediff
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
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License
    % along with this program; if not, write to the Free Software
    % Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % If you modify this code, document all changes carefully and include date
    % authorship, and a brief description of modifications
    %
    % CR, 10/05/2026: Initial coding.
    % DJ, 10/05/2026: Modify construction of V3 and dV3;
    %                   account for lower-diagonal (j<i) terms in V3
    %                   avoid "combine_terms" in computation of dV3 until
    %                   after substitution

    narginchk(4,4);

    %% Check the PIE and degree data.

    if ~isstruct(PIE) || ~isfield(PIE,'T') || ~isfield(PIE,'f') || ...
            ~isfield(PIE,'dom')
        error('V3_Liediff:InvalidPIE', ...
            'PIE must be a structure with fields T, f, and dom.');
    end

    if ~isa(PIE.f,'polyopvar')
        error('V3_Liediff:InvalidRHS', ...
            'PIE.f must be a polyopvar representation of the PIE RHS.');
    end

    if ~isnumeric(d) || ~isscalar(d) || d < 1 || d ~= floor(d)
        error('V3_Liediff:InvalidDegree', ...
            'd must be a positive integer.');
    end

    if ~isnumeric(opdeg) || ~isscalar(opdeg) || opdeg < 0 || ...
            opdeg ~= floor(opdeg)
        error('V3_Liediff:InvalidOperatorDegree', ...
            'opdeg must be a nonnegative integer.');
    end
    
    % Unpack PIE.
    Top = PIE.T; % Inverse map.
    f   = PIE.f; % Polynomial PIE.
    dom = PIE.dom; % Spatial domain.
    x   = f.vartab; % Fundamental state.

    % Lem. 13 and the current DMB implementation concern one scalar
    % distributed state. Extending this check requires a corresponding
    % multivariable version of Def. 4.
    if numel(f.varname) ~= 1 || size(f.degmat,2) ~= 1 || ...
            f.matdim(1) ~= 1
        error('V3_Liediff:ScalarPIERequired', ...
            'V3_Liediff currently supports a scalar fundamental state only.');
    end

    f_degs = sum(f.degmat,2);
    if any(f_degs==0)
        error('V3_Liediff:ConstantRHS', ...
            ['Lemma 13 assumes PIE.f contains only positive-degree ', ...
             'distributed monomials.']);
    end
    

    %% Separate the terms C_l*x^(otimes l) in the polynomial PIE into f_terms{l}.

    f_deg   = size(f.degmat,1);
    f_terms = cell(f_deg,1);
    for l = 1:f_deg
        f_terms{l}           = f;
        f_terms{l}.degmat    = f.degmat(l,:);
        f_terms{l}.C.ops     = f.C.ops(:,l);
        f_terms{l}.C.depmat2 = f.C.depmat2(l,:);
    end


    %% Build the monomial basis operator U for constructing V3 and dV3.

    pvar s s_dum
    Zmon     = monomials([s,s_dum],0:opdeg);
    Zop      = opvar();
    Zop.var1 = s;
    Zop.var2 = s_dum;
    Zop.I    = dom;

    % The degree-one case admits the multiplier component of the 3-PI basis. 
    % The current higher tensor conversion does not support multiplier factors, 
    % so only the two integral components are used for d > 1.
    if d == 1
        Zmon0    = monomials(s,0:opdeg);
        Zop.R.R0 = [Zmon0;0*Zmon;0*Zmon];
        Zop.R.R1 = [0*Zmon0;Zmon;0*Zmon];
        Zop.R.R2 = [0*Zmon0;0*Zmon;Zmon];
    else
        Zop.R.R0 = [0*Zmon;0*Zmon];
        Zop.R.R1 = [Zmon;0*Zmon];
        Zop.R.R2 = [0*Zmon;Zmon];
    end

    Z    = dopvar2ndopvar(Zop);
    ZTop = Z*Top;
    Zx   = Z*x;
    ZTx  = ZTop*x;

    %% Construct the T-Pi operators Zs1{i}=U^i*x^i and Zs2{j}=(U*T)^j*x^j.

    Zs1 = cell(d,1);
    Zs2 = cell(d,1);
    for i = 1:d
        if i == 1
            Zs1{i} = Zx;
            Zs2{i} = ZTx;
        else
            Zs1{i} = DMB(Zs1{i-1},Zx);
            Zs2{i} = DMB(Zs2{i-1},ZTx);
        end
    end


    %% Form the right factors present in dV3.

    % right_factors{j,k} = (U*T*x)^(k-1) otimes (U*x) otimes (U*T*x)^(j-k).
    %
    % After forming the associated FDP, the placeholder U*x is replaced by
    % U*C_l*x^(otimes l) using polyopvar/subs. This invokes the Lem. 7-8
    % implementation already present in PIETOOLS.
    right_factors = cell(d,d);
    for j = 1:d
        for k = 1:j
            Ztmp = Zx; % This is the placeholder which will be substituted.
            if k > 1
                Ztmp = DMB(Zs2{k-1},Ztmp);
            end
            if k < j
                Ztmp = DMB(Ztmp,Zs2{j-k});
            end
            right_factors{j,k} = Ztmp;
        end
    end


    %% Declare the unconstrained block Gram operator, common to V3 and dV3.

    % The same Q_ij variables must define both V3 and dV3. This is why the
    % two FDPs are constructed together.
    [prog, Pcell] = polyopvar_sosquadvar(prog,Zs1,Zs2,'free');
    

    %% Enforce the adjoint-symmetry condition and construct V3, dV3.

    % Equality of V3 to a scalar SOS-FDP constrains only the symmetric part
    % visible when every tensor factor is the same state x. It does not, by
    % itself, enforce Qhat_ij^* = Qhat_ji on the underlying operators.
    %
    % Introduce two independent dummy distributed states xL and xR solely
    % for these constraints. For each i <= j, impose
    %
    %  <(U*xL)^i,P_ij*(U*T*xR)^j>
    %       = <(U*T*xL)^i,P_ji'*(U*xR)^j>.
    %
    % Since xL and xR are independent, piesos_eq equates the complete operator
    % kernels rather than only their action on repeated copies of one
    % state. Including i=j therefore enforces self-adjointness of every
    % diagonal block. Restricting to i<=j avoids duplicate constraints.
    
    % Dummy distributed states.
    state_name = f.varname{1};
    left_name  = [state_name,'_left'];
    right_name = [state_name,'_right'];
    Zs1_left  = Zs1;
    Zs1_right = Zs1;
    Zs2_left  = Zs2;
    Zs2_right = Zs2;
    
    
    % Construct the LF and derivative
    V3  = 0;
    dV3 = 0;
    for i = 1:d
        Zs1_left{i}.varname  = {left_name};
        Zs1_right{i}.varname = {right_name};
        Zs2_left{i}.varname  = {left_name};
        Zs2_right{i}.varname = {right_name};
        for j = i:d
            % Symmetry constraints.
            Qhat_ij = innerprod_v2(Zs1_left{i}, Zs2_right{j}, Pcell{i,j});
            Qhat_ji_adj = innerprod_v2(Zs1_right{j}, Zs2_left{i}, Pcell{j,i});
            prog = piesos_eq(prog,Qhat_ij-Qhat_ji_adj);
            % Construct V3.
            if i==j                                                         % DJ, 10/05/2026 
                V3 = V3 + innerprod_v2(Zs1{i},Zs2{j},Pcell{i,j});
            else
                V3 = V3 + 2*innerprod_v2(Zs1{i},Zs2{j},Pcell{i,j});
            end
            % Construct the derivative
            block = 0;
            for k = 1:j
                % R_ij: substitute factor k of right_factors{j,k}, which
                % sits at global position i+k after Zs1{i}'s own i factors.
                Rijk = innerprod_v2(Zs1{i}, right_factors{j,k}, Pcell{i,j}, 'skip_combine'); % DJ, 10/05/2026
                for l = 1:f_deg
                    block = block + subs(Rijk,i+k,f_terms{l});
                end
            end
            for m = 1:i
                % L_ij = R_ji: substitute factor m of right_factors{i,m},
                % at global position j+m after Zs1{j}'s own j factors.
                Ljim = innerprod_v2(Zs1{j}, right_factors{i,m}, Pcell{j,i}, 'skip_combine'); % DJ, 10/05/2026
                for l = 1:f_deg
                    block = block + subs(Ljim,j+m,f_terms{l});
                end
            end
            if i==j
                % Add the derivative of diagonal block to full derivative
                dV3 = dV3 + block;
            else
                % Double the contribution of block (i,j) to account for
                % symmetry
                dV3 = dV3 + 2*block;
            end
        end
    end

end