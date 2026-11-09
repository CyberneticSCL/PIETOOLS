% Test_polyopvar_remove_unused_states.m -- removal of the state variables
% that do not appear in any monomial of a distributed polynomial. The
% polynomial must keep its value at all states (evaluated by integrating
% the kernels, with apply_op), and the remaining variables must be those
% that appear in some monomial, in the original order.
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-10;
pvar s
dom = [0,1];
ncases = 0;  nremoved = 0;
rng(11);

% Functionals in states a,b,c,d, in which only some of the states appear
names = {'a';'b';'c';'d'};
x = polyopvar(names,[1;1;1;1],{'s'},dom,true(4,1));
st = struct();
for k=1:4
    st.(names{k}) = rand_poly([1,1],s,2,5);
end
used_sets = {[2],[1,3],[2,4],[1,2,3,4],[4]};
for iu = 1:numel(used_sets)
    used = used_sets{iu};
    % f = sum of products <xi,xj>, xi,xj in the used states
    f = [];
    for i=used
        for j=used(used>=i)
            xi = polyopvar(names(i),1,{'s'},dom,true);
            xj = polyopvar(names(j),1,{'s'},dom,true);
            term = (1+0.5*i+j)*innerprod(xi,xj);
            if isempty(f),  f = term;   else,   f = f+term;     end
        end
    end
    % Declare the unused states as variables of f, at random positions: f has the states
    % of f0 (sorted), the unused states are inserted among them (as piesos_LF returns)
    f0 = f;
    unused = setdiff(1:4,used);
    names0 = f0.varname(:);
    [~,idcs0] = ismember(names0,names);                % states in f0, as index in names
    all_idcs = [reshape(idcs0,1,[]),unused];
    perm = randperm(numel(all_idcs));
    new_idcs = all_idcs(perm);                          % order of the states in f
    old2new = zeros(numel(all_idcs),1);                 % index in f0 -> index in f
    for k=1:numel(idcs0),   old2new(k) = find(perm==k);     end
    f = f0;
    f.varname = names(new_idcs);
    [ri,ci,vals] = find(f0.Z);
    f.Z = sparse(ri,ci,old2new(vals),size(f0.Z,1),size(f0.Z,2));
    f.dims.in = ones(numel(new_idcs),1);
    f.depmat.in = true(numel(new_idcs),1);
    f_vars = f.varname;
    g = remove_unused_states(f);
    % The variables of g are those in f that appear in a monomial of f
    appears = false(numel(f_vars),1);
    Zf = full(f.Z);
    for k=1:numel(f_vars),  appears(k) = any(Zf(:)==k);     end
    if ~isequal(reshape(g.varname,[],1),reshape(f_vars(appears),[],1))
        error("Test_polyopvar_remove_unused_states failed (case %d): remaining variables incorrect.",iu)
    end
    % The same value
    vf = eval_poly(f,st);       vg = eval_poly(g,st);
    if abs(vf)<1e-8
        error("Test_polyopvar_remove_unused_states (case %d): value is zero, the test is vacuous.",iu)
    elseif abs(vf-vg)>tol*max(1,abs(vf))
        error("Test_polyopvar_remove_unused_states failed (case %d): value %g became %g.",iu,vf,vg)
    end
    if size(g.depmat.in,1)~=numel(g.varname) || numel(g.dims.in)~=numel(g.varname)
        error("Test_polyopvar_remove_unused_states failed (case %d): dims/depmat do not match the variables.",iu)
    end
    nremoved = nremoved+(numel(g.varname)<numel(f.varname));
    ncases = ncases+1;
end
if nremoved<2
    error("Test_polyopvar_remove_unused_states: the variables were removed in only %d cases, the test is vacuous.",nremoved)
end

% A polynomial in which all states appear is unchanged
f = innerprod(x,x);
g = remove_unused_states(f);
if ~isequal(g.varname,f.varname) || ~isequal(g.Z,f.Z)
    error("Test_polyopvar_remove_unused_states failed: polynomial in which all states appear changed.")
end
ncases = ncases+1;

fprintf('Test_polyopvar_remove_unused_states passed (%d cases).\n',ncases)


function Y = eval_poly(f,st)
% Value of the functional f at the states st.(name), polynomials in s, by integrating
% the kernels, with the dummy variables numbered over the factors of each monomial
Y = 0;
Zf = full(f.Z);
for k=1:size(Zf,1)
    op = f.C.ops{1,k};
    if isempty(op),     continue,   end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    v = 1;
    for j=1:d
        xj = st.(f.varname{idcs(j)});
        if d>1
            xj.varname(strcmp(xj.varname,'s')) = {['s_',num2str(j)]};
        end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
Y = double(polynomial(Y));
end
