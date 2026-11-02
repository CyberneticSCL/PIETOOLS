function Q = subs(P,invar,inval,opts)
% Q = SUBS(P,INVAR,INVAL,OPTS) evaluates the output of the PI operator P (a
% 'sopvar' object) at a boundary value of one of its output variables,
%       (P*u)(invar=inval) = Q*u,
% with Q again a PI operator, from the space of P to the space of the output
% variables of P without invar. Writing the input variables as S1 (integrated
% out), S3 (shared with the output), and the output variables as S2
% (introduced), S3, with the parameters of P indexed by S3 as in 'sopvar':
% - invar in S2 (the output is a polynomial in invar, e.g. x --> s^2*x):
%       s^p is replaced by inval^p;
% - invar in S3 and the PI parameter acting along invar is
%     * int_{a}^{s} K(s,theta) u(theta) dtheta
%       --> int_{a}^{b} K(b,theta) u(theta) dtheta (s=b) or 0 (s=a);
%     * int_{s}^{b} K(s,theta) u(theta) dtheta
%       --> int_{a}^{b} K(a,theta) u(theta) dtheta (s=a) or 0 (s=b);
%     * multiplier K(s) u(s): UNSUPPORTED unless K(inval)=0, as
%           K(a)*u(a) (s=a) or K(b)*u(b) (s=b) cannnot be expressed in
%           terms of a PI operator.
% After the evaluation, invar is no longer an output variable. If it was in S3,
% it becomes an input variable that is integrated out (S1).
%
% INPUT
% - P:      'sopvar' object;
% - invar:  char, string, cell of char, or 'polynomial' object (pvar),
%           specifying names of output variables of P;
% - inval:  numeric array, with one value for each variable in invar, each
%           of which should be a boundary value of the interval on which 
%           the corresponding variable is defined, P.dom.out(i,:);
% - opts:   (optional) numeric scalar, or struct with field 'tol': the
%           tolerance, relative to the largest coefficient involved, to
%           which a multiplier must vanish. Defaults to 1e-10;
%
% OUTPUT
% - Q:      'sopvar' object, with the same dimensions and input variables 
%           as P, and the output variables of P without invar, such that
%           for all u, Q*u is the output of P*u evaluated at invar=inval.
%
% NOTES
% - A variable in S2 could also be evaluated at an interior point, but, as
%   for variables in S3, only boundary values are supported here.
% - Several variables are evaluated one at a time, in the order of invar.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subs (sopvar)
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
% DJ, 11/01/2026: Initial coding;

% % % Check the inputs
if ~isa(P,'sopvar')
    error("Operator must be specified as object of type 'sopvar'.")
end
if ischar(invar) || isstring(invar)
    invar = cellstr(invar);
elseif isa(invar,'polynomial')
    invar = invar.varname;
elseif ~iscell(invar)
    error("Variables must be specified as char, cell of char, or 'polynomial'.")
end
invar = invar(:);
if isa(inval,'cell')
    inval = cellfun(@(x) x,inval(:));
end
if isa(inval,'polynomial')
    % A variable evaluated at itself (x=x) is left unchanged
    same = false(numel(invar),1);
    for i=1:numel(invar)
        same(i) = ispvar(inval(i)) && strcmp(inval(i).varname{1},invar{i});
    end
    invar = invar(~same);   inval = double(inval(~same));
end
inval = inval(:);
if numel(invar)~=numel(inval)
    error("The number of variables must match the number of values at which to evaluate them.")
end
if numel(unique(invar))~=numel(invar)
    error("Variable names must be unique.")
end
tol = 1e-10;
if nargin>=4 && ~isempty(opts)
    if isstruct(opts) && isfield(opts,'tol')
        tol = opts.tol;
    elseif isnumeric(opts)
        tol = opts;
    end
end

% % % Evaluate along each variable in turn, the output of Q at name=val
Q = P;
for ii=1:numel(invar)
    P = Q;      name = invar{ii};       val = inval(ii);
    vars_in = P.vars.in(:)';            vars_out = P.vars.out(:)';
    kL = find(strcmp(vars_out,name));
    if isempty(kL)
        error("'%s' is not an output variable of the operator.",name)
    end
    % Snap val to the boundary of the domain of 'name'; atb=true for s=b.
    % Interior points are rejected: a PI operator cannot represent them for
    % S3 variables.
    dk = P.dom.out(kL,:);
    tolb = 1e-12*max(1,max(abs(dk)));
    if abs(val-dk(1))<=tolb        
        val = dk(1);    atb = false;
    elseif abs(val-dk(2))<=tolb
        val = dk(2);    atb = true;
    else
        error("Evaluation of '%s' at %g is not supported: only the boundary values %g and %g are.",name,val,dk(1),dk(2))
    end

    % Dimensions of the monomial bases. The row index of a parameter is
    % (matrix row outer, monomial inner), with the first variable the
    % outermost in the monomial index, and likewise for the columns
    nZL_arr = cellfun(@numel,P.ZL(:)');
    nnZL_arr = prod(nZL_arr);           
    nnZR_arr = prod(cellfun(@numel,P.ZR));
    p = P.dims(1);                      q = P.dims(2);
    strideL = prod(nZL_arr(kL+1:end));     % stride of the digit of name in the monomial index
    nZL_k = nZL_arr(kL);
    nL_new = nnZL_arr/nZL_k;
    degk = P.ZL{kL}(:);

    % Pass-through variables S3 (in both input and output), in the order in
    % which they index P.params. Each parameter slot along a variable of S3
    % is 1 (multiplier), 2 (int_a^s) or 3 (int_s^b). If name is in S3, all
    % parameters along that variable collapse into just an integral along
    % that variable (S3 variables becomes S1 variable).
    S3 = intersect(vars_in,vars_out);
    n3 = numel(S3);
    S3_idx = find(strcmp(S3,name));
    n3_new = n3-~isempty(S3_idx);
    % Collect the nonzero values of the new parameters, along with the
    % associated row indices and column indices.
    params_new = cell([3*ones(1,n3_new),1]);
    vals_cell = cell(size(params_new));
    rridcs_cell = vals_cell;
    ccidcs_cell = vals_cell;

    for lidx=1:numel(P.params)
        if isempty(P.params{lidx}) || nnz(P.params{lidx})==0
            continue
        end
        [rridcs,ccidcs,vals] = find(P.params{lidx});
        rridcs = rridcs(:)-1;     ccidcs = ccidcs(:);       vals = vals(:);
        % Split the row index into the matrix row (ridcs) and the monomial
        % index ZL_idcs
        ridcs = floor(rridcs/nnZL_arr);             
        ZL_idcs = mod(rridcs,nnZL_arr);
        ZLk_idcs = mod(floor(ZL_idcs/strideL),nZL_k);
        % Replace the output monomial by val^degk        
        ZL_idcs_new = floor(ZL_idcs/(strideL*nZL_k))*strideL + mod(ZL_idcs,strideL);
        rridcs_new = ridcs*nL_new + ZL_idcs_new + 1;
        vals = vals.*val.^degk(ZLk_idcs+1);
        lin_new = lidx;
        if ~isempty(S3_idx)
            % The input variable gets integrated out: collapse parameters
            % along this dimension
            slot = cell(1,n3);      
            [slot{:}] = ind2sub([3*ones(1,n3),1],lidx);
            slot = cell2mat(slot);
            lin_new = 1;
            if n3_new>0
                sl = num2cell(slot([1:S3_idx-1,S3_idx+1:n3]));
                lin_new = sub2ind([3*ones(1,n3_new),1],sl{:});
            end
            if slot(S3_idx)==1
                % Multiplier K(s)u(s), at s=val: only allowed if K(val)=0
                Kval = sparse(rridcs_new,ccidcs,vals,p*nL_new,q*nnZR_arr);
                if nnz(Kval)>0 && max(abs(nonzeros(Kval)))>tol*max(abs(vals))
                    error('sopvar:subs:multiplier',['Evaluation of ''',name,''' at ',num2str(val),' is not supported: the multiplier along ',name,' does not vanish there, so that the result cannot be expressed as a PI operator acting on the input.'])
                end
                continue
            elseif (slot(S3_idx)==2 && ~atb) || (slot(S3_idx)==3 && atb)
                continue                    % int_a^s at s=a, or int_s^b at s=b
            end
            % int_a^s at s=b, or int_s^b at s=a: full integral
        end
        rridcs_cell{lin_new} = [rridcs_cell{lin_new};rridcs_new];
        ccidcs_cell{lin_new} = [ccidcs_cell{lin_new};ccidcs];
        vals_cell{lin_new} = [vals_cell{lin_new};vals];
    end
    for lidx=1:numel(params_new)
        params_new{lidx} = sparse(rridcs_cell{lidx},ccidcs_cell{lidx},vals_cell{lidx},p*nL_new,q*nnZR_arr);
    end

    % Remove the variable from the output, and update the bases (ZR, and the
    % input variables, are unchanged: only the output side was evaluated)
    vars = P.vars;      dom = P.dom;        ZL = P.ZL;
    vars.out(kL) = [];  dom.out(kL,:) = []; ZL(kL) = [];
    Q = sopvar(params_new,vars,ZL,P.ZR,dom,P.dims);
end

end
