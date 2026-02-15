classdef quadPoly
    % This class represents polynomials in a quadratic form described
    % below.
    % For F(s,t) \in R^{m\times n}, in variables
    % s = {s1, ..., sk} and t = {t1, ..., tl}, we have the form
    % F(s,t) = (I_m ⊗ Z(s)^T) * C * (I_n ⊗ Z(t))
    % here Z(s) = Zs1(s1) ⊗ ... ⊗ Zsk(sk)
    %      Z(t) = Zt1(t1) ⊗ ... ⊗ Ztl(sl)
    % where Zsi, Zti are column vectors of monomial degrees.
    %
    % DJ, 01/21/2026: allow conversion from matrix to quadpoly
    properties
        C; % sparse coefficent matrix,  size: dim(1)*d_s x dim(2)*d_t
        Zs; % cell of exponent vectors, size: 1 x length(ns)
        % each element of cell is a SORTED column vector of size n_si
        Zt; % cell of exponent vectors, size: 1 x length(nt)
        % each element of cell is a SORTED column vector of size n_ti
        dim; % matrix dimension [m,n]
        ns; % cell of strings with variable names, size: 1 x ns
        nt; % cell of strings with variable names, size: 1 x nt
        % note, ns and nt MUST NOT have repeated varnames and must be
        % sorted
    end

    methods
        function obj = quadPoly(C, Zs, Zt, dim, ns, nt, presorted)
            % A simple constructor using all properties
            if nargin==1
                % Convert matrix to quadPoly
                Zs = {};
                Zt = {};
                dim = size(C);
                ns = {};
                nt = {};
                presorted = 1;
            elseif nargin<6 || nargin>7
                error("quadPoly must have 1 or 6 inputs, nothing else.");
            elseif nargin==6
                presorted = 0;
            end
            obj.C = sparse(C);
            obj.Zs = Zs;
            obj.Zt = Zt;
            obj.dim = double(dim(:).');
            obj.ns = ns;
            obj.nt = nt;

            if ~presorted
                obj = clean(obj);
                obj = combine(obj);
            end
        end
    end

    methods (Static, Access=?sopvar)
        % creating an auxiliary constructor for sopvar, since in sopvar we
        % will ensure that the inputs to quadPoly are in the correct
        % structure and no sorting is needed
        function obj = create(C,Zs,Zt,dim,ns,nt)
            obj = quadPoly(C,Zs,Zt,dim,ns,nt,1);
        end
    end
    methods(Static)
        % this creates a randquadPoly of size `dim', with `nmons' number of monomials,
        % in variables `var_s, var_t'. Monomials are restricted to degree
        % `deg' or lower and coefficient matrix C has sparsity `density'.
        obj = randquadPoly(dim,nmons,var_s,var_t,deg,density);
        obj = Sym2quadPoly(Fsym, sVars, tVars);
        obj = polynomial2quadPoly(Fpoly,sVars,tVars);
    end
end
