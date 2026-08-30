function Px = apply_sopvar(Pop,x)
% PX = APPLY_SOPVAR(POP,X) returns the output of the operator POP
% (specified as 'sopvar' object) to the polynomial function x
%
% INPUTS
% - Pop:    m x n 'sopvar' object representing a PI operator
%               Pop: L2^{n}[dom_in] --> L2^{m}[dom_out];
% - x:      n x q 'polynomial' object representing a polynomial function
%               x in L2^{m}[dom_in]
%
% OUTPUTS
% - Px:     m x q 'polynomial' object representing the output of Pop when
%           applied to x;
%
% NOTES
% - Integration will be performed only along variables Pop.vars.in. The
%   function x may also depend on other variables, in which case Pop will
%   act as multiplier along those directions.
% - The current implementation uses the 'polynomial' class rather than the
%   'quadPoly' class to represent the function and the operator kernels, as
%   there are issues with the 'quadPoly' class functions. Speed could
%   likely be improved using 'quadPoly'.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026  PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 06/05/2026: Initial coding;


if ~isa(Pop,'sopvar')
    error("Operator must be specified as object of type 'sopvar'.")
end
% Determine the primary variable names
vars = Pop.vars;    dom = Pop.dom;
vars1 = vars.out;   M = numel(vars1);
vars2 = vars.in;    N = numel(vars2);
dom2 = dom.in;


% Determine along which variables P3 acts as 3-PI operator
[varsS3,S3_idcs] = intersect(vars2,vars.out,'stable');
S2_idcs = setdiff(1:N,S3_idcs);
N3 = numel(varsS3);
domS3 = [dom2(S3_idcs,1),polynomial(varsS3(:)),dom2(S3_idcs,2)];   % set the domain, [a,s,b]

% Check that the function is properly specified
if isa(x,'quadPoly')
    x = quadPoly2polynomial(x);
elseif isa(x,'double')
    x = polynomial(x);
elseif ~isa(x,'polynomial')
    error("Function must be specified as object of type 'polynomial'.")
end

% Make sure the dimensions match
[m,n] = size(Pop);
if size(Pop,2)~=size(x,1)
    error("Inner dimensions of operator and function must match.")
end

% Express x in terms of dummy variables for integration
%x = subs(x,vars.in,vars2);

% Compute the left-monomial vector
if m==0
    ZL = polynomial(zeros(0,0));
else
    ZLmons = zeros(1,0);
    nZL = 1;
    for i=1:M
        nZL_i = numel(Pop.ZL{i});
        ZLmons = [repelem(ZLmons,nZL_i,1),repmat(Pop.ZL{i},nZL,1)];
        nZL = nZL_i*nZL;
    end
    ridcs = repmat((1:nZL)',m,1);
    cidcs = (1:m:m*nZL)' + (0:m-1)*(m*nZL+1);
    ZL_C = sparse(ridcs(:),cidcs(:),1);
    ZL = polynomial(ZL_C,ZLmons,vars1,[m,m*nZL]);
end

% Compute the right-monomial vector
ZRmons = zeros(1,0);
if n==0
    ZR = polynomial(zeros(0,0));
else
    nZR = 1;
    for i=1:N
        nZR_i = numel(Pop.ZR{i});
        ZRmons = [repelem(ZRmons,nZR_i,1),repmat(Pop.ZR{i},nZR,1)];
        nZR = nZR_i*nZR;
    end
    ridcs = repmat((1:nZR)',n,1);
    cidcs = (1:nZR)' + (0:n-1)*nZR + (0:n:n^2-1)*nZR;
    ZR_C = sparse(ridcs(:),cidcs(:),1);
    ZR = polynomial(ZR_C,ZRmons,vars2,[n*nZR,n]);
end

% Apply right-monomials to the function and perform full integrals
ZRx = ZR*x;
if ~isempty(S2_idcs)
    for vnum=S2_idcs
        ZRx = int(ZRx,polynomial(vars2(vnum)),dom2(vnum,1),dom2(vnum,2));
    end
end

% For each term in the operator, compute the appropriate partial integral,
% and multiply with coefficients
if N3>0
    Px = 0;
    for lidx=1:numel(Pop.params)
        % Determine what integral must be taken
        idcs = cell(1,N3+1);
        [idcs{:}] = ind2sub([3*ones(1,N3),1],lidx);
        idcs = cell2mat(idcs);
        is_int = idcs~=1;
        % Take the desired integrals
        trm = ZRx;
        for j=find(is_int)
            var = domS3(j,2);
            L = domS3(j,idcs(j)-1);
            U = domS3(j,idcs(j));
            trm = int(trm,var,L,U);
        end
        % Multiply with the appropriate coefficients
        C = Pop.params{lidx};
        trm = C*trm;
        Px = Px+trm;
    end
else
    % There are no 3-PI operators along any dimension
    Px = Pop.params{1}*ZRx;
end

% Multiply with the left monomials
Px = ZL*Px;

end