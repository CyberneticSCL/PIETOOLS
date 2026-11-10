function display(P,name,width_max,expand,call_name)
% DISPLAY(P,NAME,WIDTH_MAX,EXPAND,CALL_NAME) displays the distributed
% polynomial P in the command window, with each row of P.C displayed as a
% separate polynomial,
%   NAME_i(x,y)(s) = C_i1 x + ... + C_ij (x^{od} o y^{or}),
% for x, y the states, o the tensor product, and x^{od} the d-fold tensor
% product of x, with the operators C_ij displayed as the actual operators
% acting on the monomials if they are small enough. Called automatically
% when displaying a 'polyopvar' object without semicolon.
%
% INPUTS
% - P:      'polyopvar' object;
% - name:   (optional) char, the name to display for the polynomial.
%           Defaults to the name of the variable passed to display, or to
%           'p' if there is none. The name 'ans' is replaced by 'p';
% - width_max:  (optional) maximal number of characters of a line with
%           multiple terms. Defaults to the width of the command window, or
%           to 100 if this cannot be determined;
% - expand: (optional) logical, false to always display the operators as
%           C_ij, instead of as the actual operators if they are small
%           enough. Defaults to true;
% - call_name:  (optional) char, the name of the variable that the user can
%           call to extract the operators, if this is not NAME. Used to
%           display the polynomial that an adjoint is the adjoint of;
%
% OUTPUTS
% Display of the polynomials in the command window, with the terms
% constructed by 'get_display'.
%
% NOTES
% - If each row of P is a state variable (an identity operator acting on a
%   single state, see 'ispvar'), as for P.vartab, only the states are
%   displayed, one per row, with the spatial variables that they depend on,
%   e.g. x(s_1,s_2).
% - Each polynomial is displayed as a sum over the monomials with a nonzero
%   operator C_ij = P.C.ops{i,j}, acting on the monomial j, defined by the
%   state variables in the jth row of P.Z. The operators are displayed as
%   the actual operators if they are small enough, see 'get_display'.
%   Otherwise, they are displayed as C_ij, and can be extracted by calling
%   P.C{i,j}.
% - A monomial is displayed as the tensor product of its factors, in the
%   order of P.Z, with the tensor product of d consecutive copies of a
%   state x displayed as x^{od}. A factor of a monomial that is a known
%   function (a multiplier from R, see 'polyopvar') is displayed as 1.
% - The polynomials are named NAME_i, dropping i if P has only one row.
%   They are displayed as functions of the states in their monomials,
%   followed by the output variables, those in P.vars.out, with a
%   subscript to distinguish copies of the same variable. The output
%   variables are only displayed once, not after each of the terms. The
%   operators are named C_ij, dropping i if there is only one row, and j
%   if there is only one monomial.
% - A number at the end of the name of a state or a variable is displayed
%   as a subscript, followed by the subscript of the copy of a variable,
%   separated by a comma, and a state ending in _pie is displayed with a
%   hat over its name, e.g. x1_pie as x_1 with a hat on the x, see
%   'varname2display'.
% - All terms are displayed for each polynomial. The terms are placed on
%   the same line, unless the line would have more characters than fit in
%   the command window, in which case the whole next term is moved to the
%   next line, so that no term is split by the command window, with the +
%   aligned under the equality sign.
% - An adjoint (P.is_adjoint, see 'ctranspose') is displayed as
%       NAME = <p(x_1,...), \cdot >, where:
%   or, for a weight W other than the identity, as
%       NAME = <p(x_1,...), W, \cdot >, where:
%   with the states x_1,... in the monomials of the polynomial p that P is
%   the adjoint of, followed by the display of p, named p (or q if NAME is
%   p) and, for a weight W other than the identity, the display of W.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - display
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
% DJ, 11/06/2026: Initial coding;
% DJ, 11/06/2026: Get the subscripts from 'varname2display';
% DJ, 11/07/2026: Display adjoints, see 'ctranspose'; add input call_name;
% DJ, 11/07/2026: Test for the identity weight with is_eye;

if nargin<2 || isempty(name)
    name = inputname(1);
    if isempty(name)
        name = 'p';
    end
end
if nargin<5 || isempty(call_name)
    call_name = name;                   % the name that the user can call
end
if strcmp(name,'ans')
    % Display an unnamed object under a default name
    name = 'p';
    % An adjoint is displayed as ans = <p(x...), \cdot>, where:
    if ~P.is_adjoint
        fprintf('\n  ans = %s, where\n',name);
    end
end
if nargin<3
    width_max = [];
end
if nargin<4
    expand = true;
end

if isempty(P.C.ops)
    fprintf('\n  %s is an empty distributed polynomial, with input dimension %d and output dimension %d.\n\n',name,sum(P.dims.in),sum(P.dims.out));
    return
end
if ~P.chkval
    try
        P = initialize(P);
    catch
        builtin('disp',P);
        return
    end
end
[nr,nZ] = size(P.C.ops);

% An adjoint stores the polynomial p that it is the adjoint of, which is
% displayed after the line NAME = <p(x...), W, \cdot >, where:, followed by
% the weight W if it is not 1. The polynomial p is named p, or q if the
% adjoint is named p
if P.is_adjoint
    Pp = P;                     Pp.is_adjoint = false;      Pp.weight = 1;
    pname = 'p';
    if strcmp(call_name,'p')
        pname = 'q';
    end
    % States in the monomials of p, displayed as the arguments of p
    xnames_p = P.varname(unique(nonzeros(full(P.Z))));
    args_p = strjoin(cellfun(@varname2display,reshape(xnames_p,1,[]),'UniformOutput',false),',');
    if ~isempty(args_p)
        args_p = ['(',args_p,')'];
    end
    % The weight is not displayed if it is the identity
    has_weight = ~is_eye(P.weight,size(P.C,1));
    if has_weight
        fprintf('\n  %s = <%s%s, W, %s >, where:\n',call_name,pname,args_p,char(183));
    else
        fprintf('\n  %s = <%s%s, %s >, where:\n',call_name,pname,args_p,char(183));
    end
    display(Pp,pname,width_max,expand,call_name);
    if has_weight
        W = P.weight;
        display(W);
    end
    return
end

% A state variable (acted upon by an identity operator) is displayed as the
% state itself, with the spatial variables that it depends on, e.g. x(s_1,s_2)
[is_state,xnames,~,xvars] = ispvar(P);
if is_state
    fprintf('\n');
    for i=1:nr
        state_str = varname2display(xnames{i});
        if ~isempty(xvars{i})
            state_str = [state_str,'(',strjoin(cellfun(@varname2display,xvars{i},'UniformOutput',false),','),')'];
        end
        fprintf('  %s\n',state_str);
    end
    fprintf('\n');
    return
end

row_disp = get_display(P,expand);
fprintf('\n');
use_C = false;
for i=1:nr
    % Name of the polynomial, with the states and the output variables that it
    % depends on
    poly_name = name;
    if nr>1
        [~,sub_i] = varname2display(sprintf('%d',i));
        poly_name = [poly_name,sub_i];
    end
    if ~isempty(row_disp(i).states)
        poly_name = [poly_name,'(',strjoin(row_disp(i).states,','),')'];
    end
    if ~isempty(row_disp(i).out_args)
        poly_name = [poly_name,'(',strjoin(row_disp(i).out_args,','),')'];
    end
    print_terms([poly_name,' ='],row_disp(i).terms,row_disp(i).signs,width_max);
    use_C = use_C || (~row_disp(i).expanded && ~isempty(row_disp(i).terms));
end
fprintf('\n');

if use_C
    op_name = 'C';
    if nr>1,    [~,sub_i] = varname2display('i');   op_name = [op_name,sub_i];  end         % subscript i
    if nZ>1,    [~,sub_j] = varname2display('j');   op_name = [op_name,sub_j];  end         % subscript j
    if nr>1,    iname = 'i';    else,   iname = '1';    end
    if nZ>1,    jname = 'j';    else,   jname = '1';    end
    fprintf('  Call "%s.C{%s,%s}" to extract operator %s acting on monomial %s.\n\n',call_name,iname,jname,op_name,jname);
end

end
