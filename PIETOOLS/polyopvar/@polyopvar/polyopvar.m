classdef (InferiorClasses={?polynomial,?dpvar})polyopvar
    properties
        C = dictionary([],[]);
        pwr = zeros(0,2);
        deg = zeros(0,2);
        dom = zeros(0,2);
        vars = polynomial(zeros(0,2));
        dim = [0,0];
        dvarname = cell(0,1);
    end

    methods
        function P = polyopvar(C,pwr,deg,dom,vars)
            P.C = C;
            P.pwr = pwr;
            P.deg = deg;
            P.dom = dom;
            P.vars = vars;
        end

        function [dim] = get.dim(obj)
        % % Determine the dimensions of the operator, m x n, from the
        % % dimensions of the coefficient matrices,
        % %     m*prod(deg+1)*(q+1) x n*prod(deg+1).

        % Get the number of monomials and decision variables
        degs = obj.deg;
        pwrs = obj.pwr;
        if isempty(degs)
            degs = 0;
        end
        nZL = prod((degs(:,1)+1).^pwrs(:,1));
        N = numel(obj.deg);
        q = numel(obj.dvarname);

        % Check the dimensions of the coefficient matrices
        m_min = inf;    n_min = inf;
        m_max = 0;      n_max = 0;
        for ii=1:numel(obj.C)
            [m,n] = size(obj.C{ii});
            m_min = min(m_min,m/(nZL*(q+1)));
            m_max = max(m_min,m/(nZL*(q+1)));
            % Get the number of monomials in the dummy variables
            alpha_ii = key2inttype(pwrs,ii);
            nZR = 1;
            for jj=1:N
                is_int = all(alpha_ii{jj},1);
                nZR = nZR*(degs(jj,2)+1)^(sum(is_int));
            end
            n_min = min(n_min,n/(nZR));
            n_max = max(n_min,n/(nZR));
        end
        
        % Set the dimensions
        dim = [nan,nan];
        if m_min==m_max && round(m_min)==m_min
            dim(1) = m_min;
        end
        if n_min==n_max && round(n_min)==n_min
            dim(2) = n_min;
        end
    end

    end

end