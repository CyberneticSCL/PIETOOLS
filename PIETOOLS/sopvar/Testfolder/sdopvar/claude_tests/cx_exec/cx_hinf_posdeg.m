function [deg,copts] = cx_hinf_posdeg(dd,opts,sp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DEG,COPTS] = CX_HINF_POSDEG(DD,OPTS,SP) translates one 1-D 'poslpivar'
% call, poslpivar(prog,dim,DD,OPTS), into the degree and option arguments of
% 'poscopvar' over the container spaces SP (a cx_space_list cell).
%
% Degrees: cx_pl2pm (d{1} -> Z1, d{2} -> Z2, d{3} -> Z3; R^n is the
% identity). Options, as poslpivar reads them (poslpivar.m, nargin==4
% branch and the excludeL block):
%   psatz    0/1, passed through; copquadvar's box g(theta) is poslpivar's
%            g at the integration variable, on every pair.
%   exclude  1x4 [R^n Z1 Z2 Z3]. Entries 2:4 become copquadvar 'include'
%            on each L_2 space. Entry 1 (drop the R^n basis) has no
%            poscopvar form - a space needs one basis operator - so it is an
%            error when an R^n space is present.
%   sep      1: poslpivar sets exclude(4)=1 and makes Z2 a full-domain
%            integral; copquadvar sep=1 does the same with alpha = 4, which
%            keeps Z2's degrees.
% With the shipped light/heavy options (exclude 0, sep 0) this is
% cx_pl2pm plus psatz, the mapping test_poscopvar_vs_poslpivar verifies.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psatz = getopt(opts,'psatz',0);
exc   = getopt(opts,'exclude',[0 0 0 0]);
sep   = getopt(opts,'sep',0);
deg = cx_pl2pm(dd,sp);
copts = struct('psatz',double(psatz==1),'sep',logical(sep));
isR = cellfun(@isempty,sp);
if exc(1) && any(isR)
    error('cx_hinf_posdeg:exclude1',['poslpivar exclude(1) drops the R^n basis; '...
          'poscopvar needs one basis operator per space, so this is not expressible.'])
end
if ~any(exc(2:4)) && ~sep,  return,     end
% Included L_2 basis operators, poslpivar's order Z1, Z2, Z3 (alpha 1,2,3).
a = find(~exc(2:4));
if sep,     a(a==3) = [];   end                 % sep: Z3 is Z2's duplicate
spec = deg{find(~isR,1)};                       % the L_2 spec {Z1,Z2,Z3}
incl = a(:);    incl(incl==2 & logical(sep)) = 4;   % alpha 4 = full domain
copts.include = cell(1,numel(sp));
for k = find(~isR)
    deg{k} = spec(a);                           % one spec per kept operator
    copts.include{k} = incl;
end
end

function v = getopt(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)),  v = s.(f);  else,   v = d;  end
end
