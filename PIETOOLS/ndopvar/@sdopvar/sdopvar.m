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
    %
    % CLASS properties
    % Each kernel represents
    % K_alpha = 
    % (I_dims(1)\otimes ZL(s2,s3)') 
    %          unvec(params.A(alpha) + params.Bt(alpha)*Zd)
    %                (I_dims(2)\otimes ZR(t1,t3))
    % where K_alpha represents kernel for multiplier/integral operator
    % specified by alpha index
    % params.A has column vectors
    % Zd is a vector of decision variables, sorted (all monomial degrees 1)
    % ZL, ZR, vars, dom, dims, are same as in sopvar

    properties
        vars = struct('in',{},'out',{});
        dom = struct('in',zeros(0,2),'out',zeros(0,2));
        dims = [1,1];
        ZL = {};
        ZR = {};
        Zd = {};
        params = struct('A', {0},'Bt', {});
    end

    methods
        function P = sopvar(params,vars,Zd,ZL,ZR,dom,dims)
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