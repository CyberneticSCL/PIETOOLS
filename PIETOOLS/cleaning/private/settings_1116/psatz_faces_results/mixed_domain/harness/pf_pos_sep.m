function [prog,P] = pf_pos_sep(prog,X,side,dom,sepv)
% [PROG,P] = PF_POS_SEP(PROG,X,SIDE,DOM,SEPV) declares the positive operator
% on the SIDE spaces of X at the light degrees of prog.settings.lpi.pos
% (entry min(d,3) for a space in d variables: the multiplier basis of degree
% mult, integral bases of degree int per variable and joint in total), with
% EVERY basis operator of each space: in a registry direction d with SEPV(d)
% false the multiplier, the lower and the upper integral (alpha_d in
% {1,2,3}); with SEPV(d) true the multiplier and the full integral (alpha_d
% in {1,4}). The default storage of POSCOPVAR with DEG omitted keeps only
% the pure multiplier and the pure full integral of a 2-D space, and forces
% every direction of a 2-D space separable; this builder tests the split
% kernels (a transport or delay certificate has a kink on the diagonal).
% SEPV empty: no direction separable. SCRATCH, MMP, 11/16/2026.
[sp,dm] = copvar_space_list(X,side);
ps = struct('mult',[1 0 0],'int',[1 1 1],'joint',[2 2 2],'sep',0);
if isfield(prog,'settings') && isfield(prog.settings,'lpi') && isfield(prog.settings.lpi,'pos')
    ps = prog.settings.lpi.pos;
end
if isappdata(0,'PF_POS') && ~isempty(getappdata(0,'PF_POS')),   ps = getappdata(0,'PF_POS');   end    % degree override (M6)
vars = dom.vars;    nv = numel(vars);
if nargin<5 || isempty(sepv),   sepv = false(1,nv);     end
sepv = logical(sepv(:).');      assert(numel(sepv)==nv,'pf_pos_sep: SEPV should have one entry per registry variable');
e = @(x,c) x(min(c,end));
M = numel(sp);  deg = cell(1,M);    incl = cell(1,M);
for k = 1:M
    nd = numel(sp{k});
    if nd==0,   deg{k} = struct('int',0);   continue;   end     % R^n: the identity basis
    c = min(nd,3);  m = e(ps.mult,c);   b = e(ps.int,c);    j = e(ps.joint,c);
    [tf,idx] = ismember(sp{k},vars);    assert(all(tf),'pf_pos_sep: space variables not in the registry');
    sets = cell(1,nd);
    for d = 1:nd
        if sepv(idx(d)),    sets{d} = [1 4];    else,   sets{d} = [1 2 3];  end
    end
    g = cell(1,nd);     [g{:}] = ndgrid(sets{:});
    rows = zeros(numel(g{1}),nd);
    for d = 1:nd,   rows(:,d) = g{d}(:);    end
    rows = sortrows(rows);                      % the multiplier (all ones) first
    Zm = struct('int',m,'mult',0,'joint',m);    Zi = struct('int',b,'mult',b,'joint',j);
    dk = cell(1,size(rows,1));
    for r = 1:size(rows,1)
        if all(rows(r,:)==1),   dk{r} = Zm;     else,   dk{r} = Zi;     end
    end
    deg{k} = dk;    incl{k} = rows;
end
opts = struct('sep',sepv);  opts.include = incl;
[prog,P] = poscopvar(prog,dm,sp,dom,deg,opts);
end
