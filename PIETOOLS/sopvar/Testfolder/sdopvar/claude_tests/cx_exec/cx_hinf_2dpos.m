function [deg,copts,lossy] = cx_hinf_2dpos(pdeg,popts,sp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DEG,COPTS,LOSSY] = CX_HINF_2DPOS(PDEG,POPTS,SP) translates one 2-D
% 'poslpivar_2d(prog,dim,PDEG,POPTS)' call into poscopvar degrees and options
% over the container spaces SP (cx_space_list cell, registry {s1,s2}).
%
%   R^n space      identity basis, struct('int',[0 0]) - poslpivar_2d's u0
%   L_2[s1,s2]     settings2possopvar's 2-D translation, which carries the
%                  whole build_monoms subset-cap array (deg.subset), maps
%                  exclude(8:16) to include and sep(3:6) to a per-direction
%                  sep (EXACT by its header; checked here by Gram size
%                  against the stock program)
%   L_2[s1], L_2[s2]  NOT transcribed: no plant in scope has them, and
%                  their dx/dy + cross-space cap mapping is unverified. Error.
%
% copquadvar's sep is per registry variable and so shared by every space;
% with only R^n and L_2[s1,s2] spaces nothing else reads it. psatz 0/1 is
% passed through; psatz = 2 (ball) or 3-6 is rejected by copquadvar and is
% reported in LOSSY by settings2possopvar before that.
% exclude(1) (drop the R^n basis) is not expressible and is an error when an
% R^n space is present.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = struct('is2D',1);   s.LF_deg = pdeg;    s.LF_opts = popts;
o = settings2possopvar(s);
lossy = o.report.lossy;
lossy = lossy(startsWith(lossy,'LF'));  % only the block passed in; the 'eq'
                                        % block of s is absent on purpose
if ~isempty(lossy)          % a cone the stock builds and this one cannot
    warning('cx_hinf_2dpos:lossy','%s',strjoin(lossy,' '));
end
exc = zeros(1,16);
if isfield(popts,'exclude') && ~isempty(popts.exclude)
    exc(1:numel(popts.exclude)) = popts.exclude;
end
deg = cell(1,numel(sp));    incl = cell(1,numel(sp));
for k = 1:numel(sp)
    v = sort(sp{k});
    if isempty(v)
        if exc(1)
            error('cx_hinf_2dpos:exclude1',['exclude(1) drops the R^n basis; '...
                  'poscopvar needs one basis operator per space.'])
        end
        deg{k} = struct('int',[0 0]);
    elseif isequal(v,{'s1','s2'})
        deg{k} = o.LF.deg;      incl{k} = o.LF.opts.include;
    else
        error('cx_hinf_2dpos:space',['Space {%s}: only R^n and L2[s1,s2] spaces '...
              'are transcribed.'],strjoin(v,','))
    end
end
% psatz taken from POPTS itself: settings2possopvar maps any nonzero value to
% 1, which would turn a ball (2) into a box silently; copquadvar rejects it.
ps = 0;     if isfield(popts,'psatz') && ~isempty(popts.psatz),  ps = popts.psatz;  end
copts = struct('psatz',ps,'sep',o.LF.opts.sep);
copts.include = incl;           % assigned apart: a cell value in struct()
end                             % would build a struct array
