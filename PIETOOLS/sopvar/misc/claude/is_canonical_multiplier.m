function [tf,info] = is_canonical_multiplier(params,vars,ZL,ZR,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TF,INFO] = IS_CANONICAL_MULTIPLIER(PARAMS,VARS,ZL,ZR,DIMS) tests whether
% the parameters of a 'sopvar' or 'sdopvar' operator are in canonical
% multiplier form.
%
% A parameter whose multi-index GAM has GAM(k)==1 carries a factor
% delta(s_k-s_k'), which identifies s_k with s_k'. Monomial degrees may
% therefore be moved freely between ZL and ZR in direction k without changing
% the operator, so the stored coefficients are not determined by the operator
% alone. The canonical form removes that freedom by putting all of the degree
% on the left:
%
%   GAM(k)==1  ==>  the parameter has no content outside the columns whose
%                   ZR monomial has degree 0 in direction k.
%
% Every routine that consumes the stored coefficients directly, rather than
% the operator they represent, relies on this. In particular 'eq' can compare
% coefficients only because they are canonical, and 'lpi_eq_sdopvar' can
% constrain coefficients to zero only because that is equivalent to
% constraining the operator to zero.
%
% INPUT
% - params: parameters of the operator, either a cell array of coefficient
%           matrices ('sopvar') or a struct with fields 'A' and 'B' holding
%           cell arrays ('sdopvar');
% - vars:   'struct' with fields 'in' and 'out', each a cellstr;
% - ZL:     cell array of exponent vectors, one per variable in vars.out;
% - ZR:     cell array of exponent vectors, one per variable in vars.in;
% - dims:   1x2 array [m,n];
%
% OUTPUT
% - tf:     'logical', true if every multiplier parameter is canonical.
%           A declaration this routine cannot interpret, or a parameter whose
%           size disagrees with the declared bases, returns FALSE: callers use
%           the result to decide whether it is safe to work on the stored
%           coefficients, so an object that cannot be assessed must not pass.
%           An operator with no variable common to its input and output space
%           has no multiplier parameter and satisfies the form vacuously;
% - info:   'struct' with fields
%             .bad        indices of the offending parameters;
%             .dirs       for each, the multiplier directions violated;
%             .unchecked  indices of parameters that could not be assessed;
%             .message    a description of the first problem, or ''.
%
% See also CANONICAL_ADJOINT_MAP, SOPVAR, SDOPVAR.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - is_canonical_multiplier
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
% MMP, 08/29/2026: Initial coding

tf = true;
info = struct('bad',[],'dirs',{{}},'unchecked',[],'message','');

if ~isstruct(vars) || ~isfield(vars,'in') || ~isfield(vars,'out')
    [tf,info] = cannot_assess(info,"'vars' is not a struct with fields 'in' and 'out'");
    return
end
vin  = reshape(vars.in,1,[]);
vout = reshape(vars.out,1,[]);
vS3  = intersect(vin,vout);
n3   = numel(vS3);
if n3==0
    return              % no common direction, so no multiplier parameter,
end                     % and the form holds vacuously

ZL = reshape(ZL,1,[]);      ZR = reshape(ZR,1,[]);
if numel(ZL)~=numel(vout) || numel(ZR)~=numel(vin) || numel(dims)<2
    [tf,info] = cannot_assess(info,"the monomial bases or dimensions do not match "...
                                   +"the declared variables");
    return
end
[~,posR] = ismember(vS3,vin);

nR = cellfun(@numel,ZR);        NR = prod([nR,1]);
nL = cellfun(@numel,ZL);        NL = prod([nL,1]);
m  = dims(1);       n = dims(2);
nrow = m*NL;        nC = nrow*(n*NR);


is_sdop = isstruct(params);
if is_sdop
    if ~isfield(params,'A') || ~isfield(params,'B')
        [tf,info] = cannot_assess(info,"the parameter struct has no 'A' and 'B' fields");
        return
    end
    ncell = numel(params.A);
else
    ncell = numel(params);
end
if ncell~=3^n3
    [tf,info] = cannot_assess(info,"the parameter array holds "+num2str(ncell)+" cells, not "...
                                   +"the 3^n3 = "+num2str(3^n3)+" the declared variables require");
    return
end

sz_C = [3*ones(1,n3),1];
for k = 1:ncell
    idcs = cell(1,n3);
    [idcs{:}] = ind2sub(sz_C,k);
    gam = cell2mat(idcs);
    mult_dirs = find(gam==1);
    if isempty(mult_dirs)
        continue
    end

    if is_sdop
        Ak = params.A{k};       Bk = params.B{k};
        if (~isempty(Ak) && numel(Ak)~=nC) || (~isempty(Bk) && size(Bk,2)~=nC)
            [tf,info] = mark_unchecked(info,k);
            continue
        end
        if nnz(Ak)==0 && nnz(Bk)==0
            continue
        end
    else
        Ak = params{k};         Bk = [];
        if isempty(Ak)
            continue
        end
        if numel(Ak)~=nC
            [tf,info] = mark_unchecked(info,k);
            continue
        end
        if nnz(Ak)==0
            continue
        end
    end

    % % % In a multiplier direction the ZR degree has to be zero. Test each
    % % % direction by inspecting only the positions it forbids: that costs
    % % % numel(C) and, for a canonical parameter, touches none of its
    % % % nonzeros, so it does not scale with the decision variables.
    viol = false(1,numel(mult_dirs));
    for ii = 1:numel(mult_dirs)
        p = posR(mult_dirs(ii));
        v = (ZR{p}(:)==0);
        kk = true(1,1);
        for i = 1:numel(nR)
            if i==p
                kk = kron(kk,v);
            else
                kk = kron(kk,true(nR(i),1));
            end
        end
        bBad = find(~kk)-1;
        if isempty(bBad)
            continue
        end
        [bG,jG,rG] = ndgrid(bBad,(0:n-1)',(0:nrow-1)');
        linBad = (jG(:)*NR + bG(:))*nrow + rG(:) + 1;
        hit = ~isempty(Ak) && nnz(Ak(linBad))>0;
        if ~hit && ~isempty(Bk)
            hit = nnz(Bk(:,linBad))>0;
        end
        viol(ii) = hit;
    end
    if any(viol)
        tf = false;
        info.bad(end+1) = k;
        info.dirs{end+1} = mult_dirs(viol);
        if isempty(info.message)
            info.message = char("Parameter "+num2str(k)+" carries a multiplier in "...
                +"direction "+strjoin(cellstr(vS3(mult_dirs(viol))),', ')...
                +" but depends on the corresponding dummy variable. Degrees in a "...
                +"multiplier direction belong on the left basis ZL; see "...
                +"'is_canonical_multiplier'.");
        end
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tf,info] = cannot_assess(info,why)
% A declaration this routine cannot interpret is reported as NOT canonical.
% Callers use the result to decide whether it is safe to work on the stored
% coefficients, so an unassessable object has to fail, not pass by default.

tf = false;
if isempty(info.message)
    info.message = char("The canonical multiplier form could not be assessed because "...
                        +why+".");
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tf,info] = mark_unchecked(info,k)
% A parameter whose size disagrees with the declared bases cannot be tested,
% so it fails for the same reason a whole declaration that cannot be read does.

tf = false;
info.unchecked(end+1) = k;
if isempty(info.message)
    info.message = char("Parameter "+num2str(k)+" does not have the number of entries "...
                        +"the declared dimensions and monomial bases require, so its "...
                        +"canonical multiplier form could not be assessed.");
end

end
