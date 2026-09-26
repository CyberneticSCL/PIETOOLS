function [prog,Pm] = cx_stability_pos(prog,X,side,dom,dd,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,PM] = CX_STABILITY_POS(PROG,X,SIDE,DOM,DD,OPT) is the container
% counterpart of the 1-D stock call  [prog,P] = poslpivar(prog,X.dim,DD,OPT)
% used by every 1-D stability executive: a psd 'cdopvar' (poscopvar) over the
% SIDE ('out'|'in') spaces of container X, with poslpivar's degrees DD =
% {d1,[d2],[d3]} and options psatz / exclude / sep translated.
%
% Translation (poslpivar.m:220-310, 349-379; copquadvar.m:339-404):
%   degrees  DD normalized as poslpivar.m:220-261 does, then cx_pl2pm:
%            R^n -> identity basis; L2: Z1 (alpha 1), Z2 (alpha 2), Z3 (alpha 3);
%   psatz    0/1 passed through; copquadvar puts g at theta for every pair,
%            as poslpivar puts gs/gth/geta (poslpivar.m:302-310);
%   sep      poslpivar sep=1 drops Z3 and makes Z2 a full integral
%            (poslpivar.m:271-273, 440-441) = poscopvar sep=true, alpha {1,4};
%   exclude  entries 2:4 become options.include on each L2 space; entry 1 (the
%            R^n term) is NOT expressible (copquadvar.m:396-398 needs one
%            basis per space) and errors when an R^n space is present.
% Absent spaces (n0=0 or n1=0) need no handling: X carries none.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6 || isempty(opt),    opt = struct();     end
psatz = 0;      exc = [0 0 0 0];    sep = 0;        % poslpivar.m:179-187
if isfield(opt,'psatz') && ~isempty(opt.psatz),       psatz = opt.psatz;     end
if isfield(opt,'exclude') && ~isempty(opt.exclude),   exc = opt.exclude(:)'; end
if isfield(opt,'sep') && ~isempty(opt.sep),           sep = opt.sep;         end
if ~ismember(psatz,[0 1])
    error('cx_stability:psatz','poslpivar psatz=%g has no poscopvar form (0/1 only).',psatz)
end
if sep==1,  exc(4) = 1;     end                     % poslpivar.m:271-273

% Degree normalization, poslpivar.m:220-261.
if isempty(dd),     dd = {1,[1,1,1],[1,1,1]};   end
if isnumeric(dd)
    dd = dd(:)';
    if isscalar(dd),            dd = {dd};
    elseif numel(dd)==2,        dd = {max(dd),[dd,max(dd)],[dd,max(dd)]};
    else,                       dd = {max(dd),dd,dd};
    end
end
if numel(dd)==1,    dd{2} = dd{1}*[1 1 1];  dd{3} = dd{2};  end
for k = 2:3
    if numel(dd)<k,             dd{k} = dd{2};  end
    if numel(dd{k})==1,         dd{k} = dd{k}*[1 1 1];
    elseif numel(dd{k})==2,     dd{k}(3) = max(dd{k});
    end
end

[sp,dm] = cx_space_list(X,side);
deg = cx_pl2pm(dd,sp);
alph = [1;2;3];     if sep==1,  alph = [1;4;0];   end   % Z1, Z2 (full if sep), Z3
incl = cell(1,numel(sp));
for k = 1:numel(sp)
    if isempty(sp{k})                               % R^n space: identity basis
        if exc(1)
            error('cx_stability:excludeRn',['exclude(1)=1 with an R^n state: '...
                  'poscopvar needs one basis operator per space (copquadvar.m:396-398).'])
        end
        incl{k} = [];                               % [] = all: the one identity basis
    else
        keep = ~logical(exc(2:4));
        if ~any(keep)
            error('cx_stability:excludeL2','exclude(2:4) removes every L2 basis operator.')
        end
        incl{k} = alph(keep);                       % multi-indices (one variable)
        deg{k} = deg{k}(keep);                      % per-basis degrees, include order
    end
end
popt = struct('psatz',psatz,'sep',logical(sep),'include',{incl});
[prog,Pm] = poscopvar(prog,dm,sp,dom,deg,popt);
end
