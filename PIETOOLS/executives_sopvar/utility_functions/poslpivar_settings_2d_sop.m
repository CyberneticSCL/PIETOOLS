function [prog,P] = poslpivar_settings_2d_sop(prog,X,st,role,side,PIE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P] = POSLPIVAR_SETTINGS_2D_SOP(PROG,X,ST,ROLE,SIDE,PIE) the 2-D
% counterpart of POSLPIVAR_SETTINGS_SOP: a psd 'cdopvar' over the SIDE
% spaces of the container X, declared through 'poscopvar' with the degrees
% and options the stock 2-D executives pass to 'poslpivar_2d', translated by
% 'settings2possopvar':
%   ROLE 'lf'     the storage: LF_deg, LF_opts and the LF psatz terms;
%   ROLE 'slack'  the stock slack of an equality-imposed inequality on X:
%                 eq_deg, eq_opts with the exclusions and separations
%                 'get_eq_opts_2D' derives from the structure of X (the
%                 zero cells, the one-sided cell differences), and the eq
%                 psatz terms (exclude / sep or-ed with those of eq_opts).
% The settings are read by SETTINGS_2D_SOP from ST (an lpisettings struct
% or a 2-D settings struct). Only R^n and L2[x,y] spaces are handled: the
% stock settings name degrees for those, and the structure rules for L2[x]
% and L2[y] spaces are not transcribed (an error, not a wrong operator).
% A cone the stock declares that poscopvar cannot (settings2possopvar's
% 'lossy' report) is a warning for 'lf' and for 'slack'.
%
% The logic is that of the 10/06/2026 cx_exec translators cx_hinf_lf2d,
% cx_hinf_2dpos, cx_hinf_slack2d and cx_hinf_eqopts2d (sopvar/Testfolder/
% sdopvar/claude_tests/cx_exec), which produced programs identical to the
% stock 2-D executives on the test plants; those files are left as they
% are.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = settings_2d_sop(st);
dom = struct();     dom.vars = reshape(PIE.vars(:,1).varname,1,[]);    dom.dom = PIE.dom;
[sp,dm] = copvar_space_list(X,side);
switch role
    case 'lf'
        [deg,co] = pos2d(S.LF_deg,S.LF_opts,sp);
        [prog,P] = poscopvar(prog,dm,sp,dom,deg,co);
        for j = 1:numel(S.LF_use_psatz)
            if S.LF_use_psatz(j)~=0
                [deg,co] = pos2d(S.LF_deg_psatz{j},S.LF_opts_psatz{j},sp);
                [prog,P2] = poscopvar(prog,dm,sp,dom,deg,co);
                P = P + P2;
            end
        end
    case 'slack'
        if S.use_sosineq
            error('poslpivar_settings_2d_sop:sosineq','use_sosineq selects lpi_ineq_2d, which has no container counterpart.');
        end
        eq_opts = eqopts2d(X,S.eq_opts,1e-12);
        [deg,co] = pos2d(S.eq_deg,eq_opts,sp);
        [prog,P] = poscopvar(prog,dm,sp,dom,deg,co);
        for j = 1:numel(S.eq_use_psatz)
            if S.eq_use_psatz(j)~=0
                o = S.eq_opts_psatz{j};
                o.exclude = o.exclude | eq_opts.exclude;    o.sep = o.sep | eq_opts.sep;
                [deg,co] = pos2d(S.eq_deg_psatz{j},o,sp);
                [prog,P2] = poscopvar(prog,dm,sp,dom,deg,co);
                P = P + P2;
            end
        end
    otherwise
        error('poslpivar_settings_2d_sop:role','ROLE should be ''lf'' or ''slack''.');
end
end


function [deg,copts] = pos2d(pdeg,popts,sp)
% poslpivar_2d degrees and options -> poscopvar degrees and options over the
% spaces sp (R^n: identity basis; L2[x,y]: the translation)
s = struct('is2D',1);   s.LF_deg = pdeg;    s.LF_opts = popts;
o = settings2possopvar(s);
lossy = o.report.lossy;
lossy = lossy(startsWith(lossy,'LF'));
if ~isempty(lossy)
    warning('poslpivar_settings_2d_sop:lossy','%s',strjoin(lossy,' '));
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
            error('poslpivar_settings_2d_sop:exclude1',['exclude(1) drops the R^n basis; '...
                  'poscopvar needs one basis operator per space.']);
        end
        deg{k} = struct('int',[0 0]);
    elseif numel(v)==2
        deg{k} = o.LF.deg;      incl{k} = o.LF.opts.include;
    else
        error('poslpivar_settings_2d_sop:space',['Space {%s}: only R^n and L2[x,y] spaces '...
              'are handled.'],strjoin(v,','));
    end
end
ps = 0;     if isfield(popts,'psatz') && ~isempty(popts.psatz),  ps = popts.psatz;  end
copts = struct('psatz',ps,'sep',o.LF.opts.sep);
copts.include = incl;           % assigned apart: a cell value in struct() builds a struct array
end


function eq_opts = eqopts2d(Km,eq_opts,ztol)
% the exclusions and separations get_eq_opts_2D derives from the structure
% of Km, for R^n and L2[x,y] spaces
exc = eq_opts.exclude;      sep = eq_opts.sep;
nvs = sum(Km.space_out,2);
if any(nvs==1) || numel(Km.vars)>2
    error('poslpivar_settings_2d_sop:space','Only R^n and L2[x,y] spaces are handled.');
end
iR = find(nvs==0);      iX = find(nvs==2);
if ~isempty(iR) && all_zero(Km,iR,1)
    exc(1) = 1;
end
if ~isempty(iX)
    pidx = [8,9,10,11,13,14,12,15,16];
    for ii = 1:9
        if all_zero(Km,iX,ii),  exc(pidx(ii)) = 1;  end
    end
    c = @(g1,g2) sub2ind([3,3],g1,g2);
    if ~exc(9) && ~exc(10) && diff_le(Km,iX,c(2,1),c(3,1),ztol)
        sep(3) = 1;
    end
    if ~exc(11) && ~exc(12) && diff_le(Km,iX,c(1,2),c(1,3),ztol)
        sep(4) = 1;
    end
    if ((~exc(13) && ~exc(14)) || (~exc(15) && ~exc(16))) && ...
            diff_le(Km,iX,c(2,2),c(3,2),ztol) && diff_le(Km,iX,c(2,3),c(3,3),ztol)
        sep(5) = 1;
    end
    if ((~exc(13) && ~exc(15)) || (~exc(14) && ~exc(16))) && ...
            diff_le(Km,iX,c(2,2),c(2,3),ztol) && diff_le(Km,iX,c(3,2),c(3,3),ztol)
        sep(6) = 1;
    end
end
eq_opts.exclude = exc;      eq_opts.sep = sep;
end


function tf = all_zero(Km,idx,q)
tf = true;
for i = idx(:)',    for j = idx(:)'
    [A,Bq] = cellAB(Km.C{i,j},q,1e-12);
    if any(A) || nnz(Bq),   tf = false;     return,     end
end,                end
end


function tf = diff_le(Km,idx,q1,q2,ztol)
tf = true;
for i = idx(:)',    for j = idx(:)'
    [A1,B1] = cellAB(Km.C{i,j},q1,ztol);    [A2,B2] = cellAB(Km.C{i,j},q2,ztol);
    dB = B1 - B2;
    if any(A1-A2 > ztol) || any(nonzeros(dB) > ztol)    % implicit zeros pass
        tf = false;     return
    end
end,                end
end


function [A,Bq] = cellAB(B,q,ztol)
% cell q of a block as (constant part, decision part), zero shorthands expanded
A = zeros(0,1);     Bq = sparse(0,0);
if isempty(B),  return,     end
nL = prod([cellfun(@numel,B.ZL),1]);    nR = prod([cellfun(@numel,B.ZR),1]);
mn = B.dims(1)*nL*B.dims(2)*nR;
if isa(B,'sdopvar')
    A = B.params.A{q};      Bq = B.params.B{q};
    if numel(A)~=mn,        A = zeros(mn,1);    end         % 0 / [] shorthand
    if size(Bq,2)~=mn,      Bq = sparse(numel(B.Zd),mn);    end
else                                                        % fixed sopvar block
    A = B.params{q};        Bq = sparse(0,mn);
    if numel(A)~=mn,        A = zeros(mn,1);    end
end
A = full(A(:));     A(abs(A)<=ztol) = 0;
Bq = Bq.*(abs(Bq)>ztol);
end
