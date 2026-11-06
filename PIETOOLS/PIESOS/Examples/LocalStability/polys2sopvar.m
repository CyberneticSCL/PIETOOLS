function Z = polys2sopvar(K,var,dom)
% Z = POLYS2SOPVAR(K,VAR,DOM) declares the 'sopvar' object Z mapping
% L2[DOM] to L2^m[DOM], with the three kernels K{1:3} (multiplier,
% integral from a to s, and integral from s to b) given as polynomials.
% Replaces building an 'opvar' with R.R0, R.R1, R.R2 and converting it.
%
% INPUTS
% - K:      3 x 1 cell of m x 1 'polynomial' objects, representing the
%           kernels in the variables VAR and its dummy VAR_dum, so that
%               (Zx)(s) = K{1}(s)*x(s) + int_a^s K{2}(s,th)*x(th) dth
%                           + int_s^b K{3}(s,th)*x(th) dth;
%           K{1} may not depend on VAR_dum;
% - var:    'polynomial' object representing a single variable s;
% - dom:    1 x 2 array [a,b] specifying the domain of s;
%
% OUTPUTS
% - Z:      'sopvar' object of dimension m x 1, in the variable s, with
%           the same kernels as K;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - polys2sopvar
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
% If you modify this code, document all changes carefully and
% include date authorship, and a brief description of modifications
%
% DJ (CC), 11/06/2026: Initial coding;

vname = var.varname{1};
dname = [vname,'_dum'];
m = size(K{1},1);

% Exponents of s (rows) and s_dum (columns) of every term of every kernel;
% term t of kernel i has coefficient row cf{i}(t,:) over the m entries
ea = cell(3,1);     eb = cell(3,1);     cf = cell(3,1);
for i=1:3
    p = polynomial(K{i});
    ids = find(strcmp(p.varname(:),vname));
    idd = find(strcmp(p.varname(:),dname));
    if numel(ids)+numel(idd)~=numel(p.varname) || ~isequal(size(p),[m,1])
        error("Kernels must be "+m+" x 1 polynomials in '"+vname+"' and '"+dname+"' only.")
    end
    nt = size(p.coefficient,1);
    ea{i} = zeros(nt,1);    eb{i} = zeros(nt,1);
    if ~isempty(ids),   ea{i} = full(p.degmat(:,ids));  end
    if ~isempty(idd),   eb{i} = full(p.degmat(:,idd));  end
    cf{i} = p.coefficient;
end
if ~isempty(eb{1}) && any(eb{1})
    error("The multiplier kernel K{1} may not depend on '"+dname+"'.")
end

% Monomial bases ZL in s, ZR in s_dum: kernel i = (I_m o ZL)'*C{i}*ZR, so
% C{i}((r-1)*nZL+a,b) is the coefficient of entry r at s^ZL(a)*s_dum^ZR(b)
Zs = unique(vertcat(ea{:}));        if isempty(Zs), Zs = 0; end
Zt = unique(vertcat(eb{:}));        if isempty(Zt), Zt = 0; end
params = cell(3,1);
for i=1:3
    [tt,rr,vv] = find(cf{i});
    [~,ia] = ismember(ea{i}(tt),Zs);
    [~,ib] = ismember(eb{i}(tt),Zt);
    params{i} = sparse((rr-1)*numel(Zs)+ia,ib,vv,m*numel(Zs),numel(Zt));
end

vars = struct('in',{{vname}},'out',{{vname}});
dom = struct('in',dom(:)','out',dom(:)');
Z = sopvar(params,vars,{Zs},{Zt},dom,[m,1]);

end
