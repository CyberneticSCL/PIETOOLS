classdef(InferiorClasses={?polynomial,?opvar}) sdopvar
    % Represents PI maps from L_2^p[S1,S3] to L_2^q[S2,S3]
    %
    % This defines PI operators from one L2 space to another.
    %
    %   Pop: L_2^p[S1,S3] to L_2^q[S2,S3]
    %
    % Elements of this class are intended to be included in a container object.
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
    % CANONICAL MULTIPLIER FORM                                              % MMP, 08/29/2026
    % A parameter with alpha(k)==1 carries a factor delta(s_k-s_k'), which   % MMP, 08/29/2026
    % identifies s_k with s_k'. A coefficient at left degree a and right     % MMP, 08/29/2026
    % degree b in direction k then contributes the same monomial s_k^(a+b)   % MMP, 08/29/2026
    % as one at left degree a+b and right degree 0, so the stored            % MMP, 08/29/2026
    % coefficients are not determined by the operator. Objects of this class % MMP, 08/29/2026
    % therefore satisfy                                                     % MMP, 08/29/2026
    %                                                                        % MMP, 08/29/2026
    %   alpha(k)==1  ==>  the parameter has no content outside the columns   % MMP, 08/29/2026
    %                     whose ZR monomial has degree 0 in direction k,     % MMP, 08/29/2026
    %                                                                        % MMP, 08/29/2026
    % that is, a multiplier direction carries no dummy variable degree. The  % MMP, 08/29/2026
    % constructor enforces this by rewriting anything that does not already  % MMP, 08/29/2026
    % satisfy it, enlarging ZL where the folded degrees require it; see      % MMP, 08/29/2026
    % 'canonicalize_multiplier'. The representation is then unique, which is % MMP, 08/29/2026
    % what lets 'eq' compare coefficients and 'lpi_eq_sdopvar' constrain     % MMP, 08/29/2026
    % them to zero. Assigning to the properties directly bypasses the        % MMP, 08/29/2026
    % constructor and so bypasses the invariant; build a new object instead. % MMP, 08/29/2026

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

    methods
        function P = sdopvar(params,vars,Zd,ZL,ZR,dom,dims)
            % The canonical multiplier form is a class invariant; see the    % MMP, 08/29/2026
            % CANONICAL MULTIPLIER FORM note above. Anything handed in that  % MMP, 08/29/2026
            % is not canonical is rewritten here, so no 'sdopvar' object can % MMP, 08/29/2026
            % violate it. The check inspects only the positions the form     % MMP, 08/29/2026
            % forbids, so for a canonical object it costs numel(C) per       % MMP, 08/29/2026
            % parameter and does not scale with the decision variables.      % MMP, 08/29/2026
            [params,ZL,ZR,was_rewritten] = ...                               % MMP, 08/29/2026
                canonicalize_multiplier(params,vars,ZL,ZR,dims);             % MMP, 08/29/2026
            if was_rewritten                                                 % MMP, 08/29/2026
                warning('sdopvar:noncanonicalMultiplier',...                 % MMP, 08/29/2026
                    ['A multiplier parameter depended on a dummy variable, so it was rewritten in '...% MMP, 08/29/2026
                     'canonical form and ZL was enlarged to hold the folded degrees. The operator '...% MMP, 08/29/2026
                     'is unchanged. Build multiplier parameters with the degree on ZL to avoid this.']);% MMP, 08/29/2026
            end                                                              % MMP, 08/29/2026
            P.params = params;
            P.vars = vars;
            P.dom = dom;
            P.dims = dims;
            P.ZR = ZR; 
            P.ZL = ZL;
            P.Zd = Zd;
        end
    end
end