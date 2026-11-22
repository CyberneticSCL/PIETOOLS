function [prog,g,Qmat,Zop,Z] = piesos_ineq(prog,f,opts)
% [PROG,G,QMAT,ZOP,Z] = PIESOS_INEQ(PROG,F,OPTS) adds the distributed SOS
% constraint F>=0, for a distributed polynomial F, to the PIESOS program
% PROG.
%
% INPUTS
% - prog:   'struct' specifying a PIESOS program;
% - f:      'polyopvar' object representing a distributed polynomial.
%           The object must be a scalar-valued functional: its
%           coefficients must be 'intop' objects (or constants for
%           degree-0 monomials), with f.matdim=[1,1] and no output
%           variables;
% - opts:   (optional) struct, passed on to `piesos_sosvar`, with fields
%               'psatz':  1 x N array of orders of the psatz multipliers
%                         along each spatial variable: order q gives the
%                         multipliers (s-a)^j*(b-s)^(q-j), j=0,...,q;
%               'degs', 'exclude':  set here, based on the kernels of f
%                         (see NOTES), only if not specified. A field
%                         'deg' is renamed to 'degs', with a warning;
%               'tol':    tolerance below which coefficients of f are
%                         taken to be zero. Defaults to 1e-14;
%               'otype':  output of the operator Zop of g: pointwise (1) or
%                         tensor (2) products, or both (3), see
%                         `piesos_quadopvar`, or 'auto' for the smallest
%                         type whose terms have the structure of those of
%                         f (see `get_ineq_opts`). Defaults to 'auto';
%
% OUTPUTS
% - prog:   'struct' specifying the same PIESOS program as the input, but
%           with the added constraint F>=0. This constraint is enforced
%           by declaring a distributed SOS polynomial G using
%           `piesos_sosvar`, and enforcing F==G using `piesos_eq`;
% - g:      'polyopvar' object representing the distributed SOS variable
%           for which F==G is enforced to impose positivity of F;
% - Qmat:   m x m cell of 'dpvar' blocks representing the symmetric
%           positive semidefinite (SOS in the primary variables)
%           multiplier operator used to
%           parameterize g in the quadratic form
%               g = <Zop*Z(x), Qmat*Zop*Z(x)>;
% - Zop:    n x m 'tensopmat' object representing the monomial basis
%           operator used to parameterize g in the quadratic form;
% - Z:      m x 1 'polyopvar' object representing the distributed monomial
%           basis used to define g in the quadratic form;
%
% NOTES
% - If F is a 'dpvar', 'polynomial', or 'double' object, F>=0 is
%   enforced using `sosineq`, and if F is an 'opvar', 'opvar2d', 'dopvar',
%   or 'dopvar2d' object, using `lpi_ineq`. In these cases, only the
%   output PROG is returned.
%
% - A monomial in the basis for g (Z) cannot mix states that depend on
%   different spatial variables, e.g. x1(s1)*x2(s1,s2). Such monomials
%   are discarded with a warning (id
%   PIESOS:ineq:mixedDependence). This makes the constraint F>=0
%   conservative, in particular for F with terms such as (x1*x2)^2 that
%   require them. If no monomials remain, F=0 is imposed, with g=0 and
%   with Qmat, Zop and Z empty.
%
% - The degrees d (of the dummy variables in Zop) and p (of the primary
%   variables, in Qmat), and the types of terms in Zop, follow from the
%   exponents of the monomials in the dummy variables of the kernels of
%   f, per spatial variable. For a monomial, let e(1)>=e(2)>=... be the
%   sorted exponents of its free dummy variables, where dummy variables
%   evaluated at the same point (a delta term) count as one variable.
%   The terms of g, with dummy degree d and primary degree p, have all
%   but one exponent <=d, and that one <=d+2p+1 (from int_{s_i}^b s^k
%   ds). So, to recover a monomial,
%       d>=e(2),    d+2p+1>=e(1),    2p>=e_delta,
%   where e_delta is the exponent of a delta term.
%
% - The multiplier term in Zop (exclude(1)) is included only along
%   variables in which f has delta terms: otherwise its Gram blocks are
%   zero, as the delta part of g is a square. The integral terms
%   (exclude(2:3)) are included only along variables in which f has free
%   dummy variables: otherwise the multiplier block reproduces f.
%
% - If f has no decision variables, d and p are the smallest degrees
%   recovering all monomials. The total degree of the kernels does not
%   give d: for f=(int s^3*x ds)^2 it is 6, but d=3 is needed.
%
% - If f has decision variables, only the monomials with a numeric
%   coefficient need to be recovered: the others may be set to zero by
%   the solver. However, a certificate may need other monomials than
%   those of f: for -dV of a Lyapunov functional V, d=5 and p=2,
%   recovering all monomials (e(1)=9), is infeasible, whereas d=p=3 is
%   feasible, so that the smallest d and p depend on the data, not on
%   the monomials. Therefore, d and p are at least those for the numeric
%   monomials, and at least the smaller of ceil(D/4), for D the maximal
%   total degree of the kernels, and the d=p for which all monomials
%   fit: d+2p+1>=e(1), 2p>=e_delta.
%
% - A psatz multiplier of order q limits the degree of Qmat in s to
%   q+2(p-ceil(q/2)).
%
% - The same d is used for each factor of a distributed monomial, as the
%   order of the factors in the kernels of f is not unique. If F>=0 is
%   infeasible for the d and p chosen, increase them using opts.degs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_ineq
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
% DJ, 10/15/2026: Initial coding;
% DJ, 11/22/2026: Option 'otype' ('auto' by default), and use 
%                   `get_ineq_opts` to set basis, degrees and options;

% % Check the inputs
% Check the program structure
if ~isa(prog,'struct')
    error("PIESOS program structure must be specified as 'struct'.")
end
% Check the polynomial functional
if isa(f,'dpvar') || isa(f,'polynomial') || isa(f,'double')
    % Enforce positivity of finite-dimensional polynomials using sos_ineq;
    if nargout>1
        error("Multiple output arguments is not supported in this case.")
    end
    prog = sosineq(prog,f);
    return
elseif isa(f,'dopvar') || isa(f,'dopvar2d') || isa(f,'opvar') || isa(f,'opvar2d')
    % Enforce positivity of operator variables using lpi_ineq
    if nargout>1
        error("Multiple output arguments is not supported in this case.")
    end
    prog = lpi_ineq(prog,f);
    return
elseif ~isa(f,'polyopvar')
    error("Distributed polynomial functional must be specified as object of type 'polyopvar'.")
elseif f.is_adjoint
    error("An adjoint (u') is not a polynomial in the states: use u'*v instead.")
end
% Make sure f is properly declared, and a scalar-valued functional, as
% piesos_eq(prog,f-g) compares it with the scalar functional g
if ~f.chkval
    f = initialize(f); 
end
if ~isequal(f.matdim,[1,1]) || ~isempty(f.vars.out)
    error("Distributed polynomial must be a scalar-valued functional (no output variables).")
end
% Check the options
if nargin==2
    opts = struct();
elseif nargin==3
    if ~isa(opts,'struct')
        error("Options for enforcing positivity must be specified as struct.")
    end
else
    error("Too many input arguments.")
end
% Tolerance below which coefficients are taken to be zero
if isfield(opts,'tol') && ~isempty(opts.tol)
    tol = opts.tol;
else
    tol = 1e-14;
end
% Merge duplicate rows, drop zero rows, and remove unused variables
f = clean_poly(f,tol);
% Type of output of Zop: by default from the structure of f ('auto',
% resolved by get_ineq_opts). An otype in opts.L or opts.R still applies
% to that operator (see piesos_quadopvar)
if ~isfield(opts,'otype')
    opts.otype = 'auto';
end

% The distributed monomial basis Z of g (half the Newton polytope of f),
% and the options of piesos_sosvar: otype for 'auto', and the degrees and
% term types of Zop from the kernels of f unless specified (see NOTES)
[Z,opts] = get_ineq_opts(f,opts);
if isempty(Z)
    % The only way for the polynomial to be positive semidefinite, is for
    % it to be zero
    g = constant2polyopvar(0);
    prog = enforce_eq(prog,f,g,tol);
    Qmat = cell(0,0);
    Zop = tensopmat();
    Z = polyopvar();
    return
end

% Declare distributed SOS functional
[prog,g,Qmat,Zop] = piesos_sosvar(prog,Z,opts);

% Enforce f==g>=0

prog = enforce_eq(prog,f,g,tol);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function prog = enforce_eq(prog,f,g,tol)
% PROG = ENFORCE_EQ(PROG,F,G,TOL) adds the constraint F=G to the program, for F
% and G scalar 'polyopvar' functionals. The constant (degree-0) terms are
% equated directly, as they are stored as 'tensopvar' in F but as 'intop'
% in G, and 'tensopvar' + 'intop' is not supported.
% The remaining terms are equated using piesos_eq.
%
% INPUT
% - prog:   'struct', the PIESOS program;
% - f, g:   scalar 'polyopvar' functionals;
% - tol:    scalar, tolerance below which coefficients are taken zero;
%
% OUTPUT
% - prog:   the program with the constraint F=G added.

% Equate the constant terms (degree-0 monomials) as scalars
[f,cf] = split_constant(f);
[g,cg] = split_constant(g);
if ~isempty(cf) || ~isempty(cg)
    if isempty(cf)
        cf = 0;
    end
    if isempty(cg)
        cg = 0;
    end
    c = cf-cg;
    if isnumeric(c)
        c = polynomial(c);
    end
    prog = soseq(prog,c);
end

% Equate the remaining terms
has_f = size(f.Z,1)>0;
has_g = size(g.Z,1)>0;
if has_f && has_g
    prog = piesos_eq(prog,f-g,tol);
elseif has_f
    prog = piesos_eq(prog,f,tol);
elseif has_g
    prog = piesos_eq(prog,-g,tol);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [f,c] = split_constant(f)
% [F,C] = SPLIT_CONSTANT(F) removes the degree-0 monomial from the scalar
% 'polyopvar' F, and returns its coefficient C (empty if F has none), a
% 'double' (for a 'tensopvar' coefficient), or a 'dpvar' (for an 'intop'
% coefficient).
%
% INPUT
% - f:  scalar 'polyopvar' functional;
%
% OUTPUT
% - f:  the functional without its degree-0 monomial;
% - c:  'double' or 'dpvar', the sum of the constant coefficients, or []
%       if f has no degree-0 monomial.

% Sum the scalar coefficients of the degree-0 monomials: the kernel of an
% 'intop' without integrals, the parameter of a single 'sopvar' factor
c = [];
is_const = ~any(f.Z,2);
for i=reshape(find(is_const),1,[])
    Cop = f.C.ops{1,i};
    if isa(Cop,'intop') && isscalar(Cop.params) && isempty(Cop.intmat)
        Kc = kernels(Cop);
        ci = Kc{1};
    elseif isa(Cop,'tensopvar') && isequal(size(Cop.ops),[1,1]) && isa(Cop.ops{1,1},'sopvar') && isscalar(Cop.ops{1,1}.params)
        ci = full(Cop.ops{1,1}.params{1});
        if isempty(ci)
            ci = 0;
        end
    elseif isnumeric(Cop)
        ci = Cop;
    else
        error("Constant coefficients must be scalar 'tensopvar' or 'intop' objects without variables.")
    end
    % Add to the sum
    if isempty(c)
        c = ci;
    else
        c = c+ci;
    end
end
% Remove the constant monomial and its coefficient
if any(is_const)
    keep = ~is_const;
    f.C.ops = f.C.ops(:,keep);
    f.C.dims.in = f.C.dims.in(keep);
    f.C.depmat.in = f.C.depmat.in(keep,:);
    f.C.chkval = true;
    f.Z = f.Z(keep,:);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
