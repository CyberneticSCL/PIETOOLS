function X = on_registry_sop(X,R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X = ON_REGISTRY_SOP(X,R) restates the container X on the variable
% registry (vars, dom) of the container R, which must contain X's. The
% blocks are untouched; only the space masks are widened. Needed before
% adding operators built on different registries, e.g. -gam*I_w on R -> R
% and Tw'*(P*Bw) on the state's variable (@cdopvar/plus demands identical
% registries).
%
% Initial coding MMP, 10/08/2026 (the logic of the frozen test helper
% cx_on_registry, in the library so that the executives depend on the
% library only).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(X.vars,R.vars) && isequal(X.dom,R.dom),  return,     end
meta = metadata(X);
[tf,loc] = ismember(X.vars,R.vars);
if ~all(tf)
    error('on_registry_sop:notSubset','X''s registry is not contained in R''s.')
end
so = false(size(X.C,1),numel(R.vars));  so(:,loc) = X.space_out;
si = false(size(X.C,2),numel(R.vars));  si(:,loc) = X.space_in;
meta.vars = R.vars;     meta.dom = R.dom;
meta.space_out = so;    meta.space_in = si;
if isa(X,'cdopvar'),    X = cdopvar(X.C,meta);
else,                   X = copvar(X.C,meta);
end
end
