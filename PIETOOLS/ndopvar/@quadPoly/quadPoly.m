classdef quadPoly
    % Polynomial of the form
    % F(s,t) = (I \otimes Z(s)^T) * C * (I \otimes Z(t)) \in R^{m\times n}
    properties
        C; % sparse coefficent size: dim(1)*d_s x dim(2)*d_t
        Zs; % exponent matrix, size: d_s x length(ns)
        Zt; % exponent matrix, size: d_t x length(nt)
        dim; % matrix dimension [m,n]
        ns; % cell of strings with variable names
        nt; % cell of strings with variable names
    end

    methods
        function obj = quadPoly(C, Zs, Zt, dim, ns, nt)
            obj.C = sparse(C);
            obj.Zs = double(Zs);
            obj.Zt = double(Zt);
            obj.dim = double(dim(:).'); 
            obj.ns = ns;
            obj.nt = nt;
        end
    end

    methods(Static)
        obj = randquadPoly(dim,nmons,var_s,var_t,deg,density);
    end
end
