function tl = exec_tools_sop(PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TL = EXEC_TOOLS_SOP(PIE,ST) the dimension-dependent pieces the builders of
% the container executives share, for a 1-D or a 2-D PIE:
%   nv             PIE.dim;
%   dom            the domain argument of lpivar_cdopvar / poscopvar
%                  (PIE.dom in 1-D; struct vars/dom in 2-D);
%   eppos, epneg   the strictness constants the stock executives read
%                  (1-D: ST.eppos, ST.eppos2, ST.epneg; 2-D: the settings_2d
%                  fields, SETTINGS_2D_SOP), eppos one weight per space;
%   eppos_default  true when the settings carried no eppos;
%   storage(prog,X,side)   the stock storage over the SIDE spaces of X
%                  (POSLPIVAR_SETTINGS_SOP 'lf' / POSLPIVAR_SETTINGS_2D_SOP);
%   slack(prog,X)  the stock slack for an inequality on X ('slack' role);
%   eye(side,e)    diag(e_k I) on the state space (the local eye_state);
%   ident(Bop,side) the identity on the SIDE space of the opvar Bop, as a
%                  container (I_w = ident(B1,2), I_z = ident(C1,1));
%   toop(C)        a solved container to an opvar (1-D) / opvar2d (2-D);
%   zdeg           the degrees of a free operator variable (1-D: ST.ddZ;
%                  2-D: one cap, the largest entry of settings_2d.Zop_deg,
%                  on every role of lpivar_cdopvar: a superset of the stock
%                  lpivar_2d monomials);
%   S              the 2-D settings struct ([] in 1-D).
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nv = PIE.dim;
vars = PIE.vars;
vnames = reshape(vars(:,1).varname,1,[]);
tl = struct('nv',nv,'S',[]);
if nv==1
    tl.dom = PIE.dom;
    tl.eppos = [getf(st,'eppos',1e-4); getf(st,'eppos2',1e-6)];
    tl.epneg = getf(st,'epneg',0);
    tl.eppos_default = ~isfield(st,'eppos') || isempty(st.eppos);
    tl.storage = @(prog,X,side) poslpivar_settings_sop(prog,X,st,'lf',side,PIE.dom);
    tl.slack   = @(prog,X) poslpivar_settings_sop(prog,X,st,'slack','out',PIE.dom);
    tl.toop    = @(C) copvar2opvar(C);
    tl.zdeg    = getf(st,'ddZ',1);
elseif nv==2
    S = settings_2d_sop(st);
    tl.S = S;
    tl.dom = struct('vars',{vnames},'dom',PIE.dom);
    tl.eppos = S.eppos;     tl.epneg = S.epneg;     tl.eppos_default = S.eppos_default;
    tl.storage = @(prog,X,side) poslpivar_settings_2d_sop(prog,X,st,'lf',side,PIE);
    tl.slack   = @(prog,X) poslpivar_settings_2d_sop(prog,X,st,'slack','out',PIE);
    tl.toop    = @(C) copvar2opvar2d(C,vnames);
    tl.zdeg    = zcaps(S.Zop_deg);
else
    error('exec_tools_sop:dim','The stock executives exist for 1-D and 2-D PIEs only.');
end
tl.eye   = @(side,e) eye_state(PIE,side,e);
tl.ident = @(Bop,side) op2copvar_sop(mat2opvar(eye(size(Bop,side)),Bop.dim(:,side),vars,PIE.dom));
end


function Im = eye_state(PIE,side,e)
% the weighted identity diag(e_k I_{n_k}) on the state space of the PIE as a
% fixed container, n_k the component counts of PIE.T on SIDE (1: output,
% the PDE state; 2: input, the fundamental state), e_k the weight of space k
% (two entries in 1-D, four in 2-D, or one scalar for every space): the
% eppos I the stock executives add (blkdiag on mat2opvar in 1-D, the
% four-entry eppos on opvar2d in 2-D). Distinct from the library's
% eye_copvar_sop(dims,spaces,dom), the plain identity on given spaces.
n = PIE.T.dim(:,side);
if isscalar(e),     e = e*ones(numel(n),1);     end
if numel(e)~=numel(n)
    error('exec_tools_sop:e','E should have one weight per space (%d) or be a scalar.',numel(n));
end
blocks = cell(1,numel(n));
for k = 1:numel(n),     blocks{k} = e(k)*eye(n(k));    end
Imat = blkdiag(blocks{:});
if PIE.dim==1
    Pop = mat2opvar(Imat,n,PIE.vars,PIE.dom);
else
    Pop = opvar2d(Imat,[n,n],PIE.dom,PIE.vars);
end
Im = op2copvar_sop(Pop);
end


function v = getf(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)),  v = s.(f);  else,   v = d;  end
end


function d = zcaps(Z)
% one degree cap from the stock 2-D free-variable degrees (fields dx, dy, d2
% of nested cells of arrays): their largest entry
m = max([cellmax(Z.dx), cellmax(Z.dy), cellmax(Z.d2), 0]);
d = struct('mult',m,'int',[m m],'out',m,'in',m);
end


function m = cellmax(c)
m = 0;
if iscell(c)
    for k = 1:numel(c),     m = max(m,cellmax(c{k}));   end
elseif isnumeric(c) && ~isempty(c)
    m = max(c(:));
end
end
