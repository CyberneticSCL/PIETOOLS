function prog = lpi_ineq(prog,P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = LPI_INEQ(PROG,P,OPTS) takes an LPI program structure 'prog'
% and a (PI operator) decision variable P, and adds inequality constraints 
% enforcing P>=0.
% 
% INPUT
% - prog:   'struct' specifying LPI program structure to modify;
% - P:      Object of type 'dpvar' or 'dopvar', specifying a variable 
%           or PI operator variable of which to enforce positivity;
% - opts:   (optional) determines if the operator needs to be pure or 
%           to include psatz term. If opts.pure=1, then R_0 term is 
%           excluded. If opts.psatz=1, then a psatz term is added.
% 
% OUTPUT 
% - prog:   Same LPI program structure as input, but with the added
%           constraints enforcing P==0.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - lpi_ineq
%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding MMP, SS, DJ  - 09/26/2021
% 07/07/2021, DJ: Adjusted (slightly) for dpvar (dopvar) implementation;
% 10/19/2024, DJ: Add support for non-opvar decision variables;
% 11/29/2024, DJ: Make sure dummy variables match those in LPI program;


% Extract the inputs
switch nargin
    case 1
        error('Not enough inputs!')
    case 2
        opts.psatz = 0;
        opts.pure = 0;
        opts.sep = 0;
    case 3
        if ~isfield(opts,'psatz')
            opts.psatz=0;
        end
        if ~isfield(opts,'pure')
            opts.pure=0;
        end
        if ~isfield(opts,'sep')
            opts.sep=0;
        end
end

% Check that the specified (operator) variable is of appropriate type.
if isa(P,'polynomial') || isa(P,'double')
    error('Enforcing equality constraints on fixed values or polynomials is not supported.')
elseif isa(P,'dpvar')
    % If P is not an opvar, we can enforce inequality using just sosineq.
    prog = sosineq(prog,P);
    return
elseif ~isa(P,'dopvar') && ~isa(P,'dopvar2d')
    error('Variable of which to enforce positivity should be specified as object of type ''dpvar'', ''dopvar'', or ''dopvar2d''.')
end
% Make sure the opvar dummy variables match those in the LPI program.       (11/29/2024, DJ)
n = size(prog.dom,1);
dumvars = prog.vartable(n+1:2*n);
if any(~ismember(P.var2.varname,dumvars.varname))
    P = subs_vars_op(P,P.var2,dumvars);
end
% Pass 2D operator to associated lpi_ineq function.
if isa(P,'dopvar2d')
    prog = lpi_ineq_2d(prog,P,opts);
    return
end

% Establish dimension of new operator.
dim = P.dim;
if dim(:,1)~=dim(:,2)
    error('Non-symmetric operators cannot be sign definite. Unable to set the inequality');
end

% Establish degrees of new operator to match those of P.
d2 = degbalance(P);

% Check if certain parameters can be excluded.
tol = 1e-14;
if opts.pure == 1
    options2.exclude = [0,1,0,0];
    options3.exclude = [0,1,0,0];
else
    options2.exclude = [0,0,0,0];
    options3.exclude = [0,0,0,0];
end
if (isa(P.P,'double') && max(max(abs(P.P)))<tol) || all(max(max(abs(P.P.C)))<tol)
    options2.exclude(1) = 1;
    options3.exclude(1) = 1;
end
if (isa(P.R.R0,'double') && max(max(abs(P.R.R0)))<tol) || all(max(max(abs(P.R.R0.C)))<tol)
    options2.exclude(2) = 1;
    options3.exclude(2) = 1;
end
if (isa(P.R.R1,'double') && max(max(abs(P.R.R1)))<tol) || all(max(max(abs(P.R.R1.C)))<tol)
    options2.exclude(3) = 1;
    options3.exclude(3) = 1;
end
if (isa(P.R.R2,'double') && max(max(abs(P.R.R2)))<tol) || all(max(max(abs(P.R.R2.C)))<tol)
    options2.exclude(4) = 1;
    options3.exclude(4) = 1;
end


% Declare a positive operator Deop and enforce Dop-P==0.
if opts.psatz == 1
    options3.psatz=1;
    [prog, Deop] = poslpivar(prog,dim,d2,options2);
    [prog, De2op] = poslpivar(prog,dim,d2,options3);
    prog = lpi_eq(prog,Deop+De2op-P,'symmetric'); %Dop=Deop+De2op
else
    [prog, Deop] = poslpivar(prog,dim,d2,options2);
    prog = lpi_eq(prog,Deop-P,'symmetric'); %Dop=Deop
end
end