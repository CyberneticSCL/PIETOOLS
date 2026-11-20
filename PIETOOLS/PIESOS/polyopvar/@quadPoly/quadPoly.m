classdef quadPoly
    % This class represents polynomials in a quadratic form described
    % below.
    % For F(s,t) \in R^{m\times n}, in variables
    % s = {s1, ..., sk} and t = {t1, ..., tl}, we have the form
    % F(s,t) = (I_m ⊗ Z(s)^T) * C * (I_n ⊗ Z(t))
    % here Z(s) = Zs1(s1) ⊗ ... ⊗ Zsk(sk)
    %      Z(t) = Zt1(t1) ⊗ ... ⊗ Ztl(sl)
    % where Zsi, Zti are column vectors of monomial degrees.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - quadPoly
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
    % Initial coding SS - 01/16/2026

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
            % OBJ = QUADPOLY(C,ZS,ZT,DIM,NS,NT,PRESORTED) sets the properties;
            % OBJ = QUADPOLY(C) converts a matrix C to a constant quadPoly
            % (no variables), or returns a 'quadPoly' C unchanged.
            %
            % INPUT
            % - C:          dim(1)*d_s x dim(2)*d_t coefficient matrix
            %               (stored sparse), or an m x n matrix (one input);
            % - Zs, Zt:     1 x numel(ns), 1 x numel(nt) cells of sorted
            %               column vectors of exponents;
            % - dim:        1 x 2 array [m,n];
            % - ns, nt:     cellstr of sorted, distinct variable names;
            % - presorted:  (optional) logical, unused (no sorting is
            %               performed);
            %
            % OUTPUT
            % - obj:        'quadPoly' object.
            if nargin==1
                if isa(C,'quadPoly')
                    obj = C;
                    return
                else
                    % Convert matrix to quadPoly
                    Zs = {};
                    Zt = {};
                    dim = size(C);
                    ns = {};
                    nt = {};
                    presorted = 1;
                end
            % Otherwise 6 inputs, or 7 with presorted
            elseif nargin<6 || nargin>7
                error("quadPoly must have 1 or 6 inputs, nothing else.");
            elseif nargin==6
                presorted = 0;
            end
            % Set the properties, the coefficients stored sparse
            obj.C = sparse(C);
            obj.Zs = Zs;
            obj.Zt = Zt;
            obj.dim = double(dim(:).');
            obj.ns = ns;
            obj.nt = nt;

            % Sorting and combining are disabled: the condition is always false
            if ~presorted && 0
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
        % OBJ = CREATE(C,ZS,ZT,DIM,NS,NT) declares a 'quadPoly' object with the
        % properties given, without sorting; the inputs must already be sorted.
        %
        % INPUT
        % - C, Zs, Zt, dim, ns, nt:     the properties, as for the constructor;
        %
        % OUTPUT
        % - obj:    'quadPoly' object.
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
        obj = quadPoly2polynomial(Fquadpoly);
    end
end
