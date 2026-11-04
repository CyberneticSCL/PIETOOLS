function p0_f2cmp(a, b)
% P0_F2CMP(A,B) compares the F2 operators of passes A and B independently
% of their block partition, and writes ../compare_F2_pass<A>_pass<B>.txt.
% Per PIE and operator present in both passes: both containers are applied,
% block by block (apply_sopvar), to the same random polynomial test
% functions drawn per scalar input component, and every scalar output
% component is evaluated at the same random points; the line gives the
% largest relative difference max|Ya-Yb| / max(1,max|Ya|), or the reason
% the operators cannot be compared (input or output components, their
% variables, or the domain differ).
%
% INPUT
% - a, b:   double, the pass numbers (F2/pass<A>, F2/pass<B>).
%
% NOTES
% - P0_COMPARE compares the stored outputs of F2, which P0_APPLY draws per
%   block: a container split into more blocks (one per component, as the
%   converter of b4256de6 builds them) gets other test functions there, so
%   equal operators can differ. This routine does not depend on the split.
% - Needs the classes on the path (P0_PATH): pass 1 objects load into the
%   classes of the current export.
% - 3 test functions per input component, 5 points per output component,
%   total degree <= 2 (as P0_APPLY); the seed depends on the label and the
%   operator name only.
% - A PIE with an operator above 1e-8, or not comparable, is compared once
%   more up to a reordering of its state components (P0_F2PERM): the
%   converter of b4256de6 orders the state components differently.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
p0_path();
q = 3;  P = 5;
rep = fullfile(root,sprintf('compare_F2_pass%d_pass%d.txt',a,b));
fid = fopen(rep,'w');
fprintf(fid,'Phase 0, F2 by components: pass %d vs pass %d, %s\n', a, b, char(datetime('now')));
fprintf(fid,'(label operator: relative difference | reason not compared; [blocks a / blocks b])\n');
cnt = struct('eq12',0,'eq8',0,'gt8',0,'struct',0,'only',0,'perm',0,'permfail',0);
worst = 0;  wlab = '';
for f = dir(fullfile(root,'F2',sprintf('pass%d',a),'*.mat'))'
    g = fullfile(root,'F2',sprintf('pass%d',b),f.name);
    if ~exist(g,'file'),    fprintf(fid,'%s missing in pass %d\n',f.name,b);  continue,   end
    A = load(fullfile(f.folder,f.name));    B = load(g);
    na = {A.fix.name};  nb = {B.fix.name};
    for o = setxor(na, nb)
        in = 'b';  if ismember(o{1},na),   in = 'a';   end
        fprintf(fid,'%s %s: in pass %s only\n', A.label, o{1}, in);
        cnt.only = cnt.only + 1;
    end
    bad = false;
    for o = intersect(na, nb, 'stable')
        Xa = A.fix(strcmp(na,o{1})).X;  Xb = B.fix(strcmp(nb,o{1})).X;
        tag = sprintf('[%s / %s]', mat2str(size(Xa.C)), mat2str(size(Xb.C)));
        [d, why] = p0_cmpop(Xa, Xb, q, P, sum(double([A.label,o{1}])));
        if isempty(why)
            fprintf(fid,'%s %s: %.3g %s\n', A.label, o{1}, d, tag);
            if d <= 1e-12,     cnt.eq12 = cnt.eq12 + 1;
            elseif d <= 1e-8,  cnt.eq8 = cnt.eq8 + 1;
            else,              cnt.gt8 = cnt.gt8 + 1;   bad = true;
            end
            if d > worst,   worst = d;  wlab = [A.label,' ',o{1}];  end
        else
            fprintf(fid,'%s %s: %s %s\n', A.label, o{1}, why, tag);
            cnt.struct = cnt.struct + 1;    bad = true;
        end
    end
    % A PIE with an operator above 1e-8 or not comparable: retry up to a
    % reordering of its state components (one permutation for all operators)
    if bad
        [e, px, det] = p0_f2perm(A, B);
        fprintf(fid,'%s: up to the state permutation %s, largest relative error %.3g (%s)\n', ...
            A.label, mat2str(px), e, det);
        if e <= 1e-10,  cnt.perm = cnt.perm + 1;    else,   cnt.permfail = cnt.permfail + 1;    end
    end
end
fprintf(fid,['F2 by components: %d operators <= 1e-12, %d in (1e-12, 1e-8], %d > 1e-8 ', ...
    '(largest %.3g, %s); %d not comparable; %d in one pass only\n'], ...
    cnt.eq12, cnt.eq8, cnt.gt8, worst, wlab, cnt.struct, cnt.only);
fprintf(fid,['F2 up to a state permutation: of the PIEs with an operator above 1e-8 or not ', ...
    'comparable, %d agree to <= 1e-10, %d do not\n'], cnt.perm, cnt.permfail);
fclose(fid);
fprintf('\nDONE F2 by components %d %d: %s\n', a, b, rep);
end

function [d, why] = p0_cmpop(Xa, Xb, q, P, seed)
% Relative difference of two containers on the same per-component inputs,
% or WHY they cannot be compared ('' if they can).
d = NaN;    why = '';
% Variables of each scalar component, by name (the order of vars may differ)
[ia, oa] = p0_comps(Xa);    [ib, ob] = p0_comps(Xb);
if numel(ia)~=numel(ib),    why = sprintf('input components %d / %d',numel(ia),numel(ib));    return,     end
if numel(oa)~=numel(ob),    why = sprintf('output components %d / %d',numel(oa),numel(ob));  return,     end
if ~isequal(ia,ib),     why = 'variables of the input components differ';    return,     end
if ~isequal(oa,ob),     why = 'variables of the output components differ';   return,     end
va = cellstr(string(Xa.vars(:)));   vb = cellstr(string(Xb.vars(:)));
[~, kb] = ismember(va, vb);
if any(kb==0) || ~isequal(Xa.dom, Xb.dom(kb,:)),  why = 'domains differ';   return,     end
rng(seed,'twister');
xin = cell(numel(ia),1);
for c = 1:numel(ia),    xin{c} = p0_randpoly(1, q, ia{c});      end
pts = cell(numel(oa),1);
for c = 1:numel(oa)
    di = Xa.dom(ismember(va, oa{c}),:);
    pts{c} = di(:,1)' + rand(P, numel(oa{c})).*(di(:,2)-di(:,1))';
end
Ya = p0_applyc(Xa, xin, oa, pts);   Yb = p0_applyc(Xb, xin, ob, pts);
d = max(abs(Ya(:)-Yb(:))) / max(1, max(abs(Ya(:))));
if isempty(d),  d = 0;  end
end

function [ic, oc] = p0_comps(X)
% Variable names of each scalar input (IC) and output (OC) component, in
% the order of X.vars, as cells of cellstr.
v = cellstr(string(X.vars(:)));
ic = {};    oc = {};
for j = 1:numel(X.dim_in),  ic = [ic; repmat({v(logical(X.space_in(j,:)))'},X.dim_in(j),1)];      end %#ok<AGROW>
for i = 1:numel(X.dim_out), oc = [oc; repmat({v(logical(X.space_out(i,:)))'},X.dim_out(i),1)];    end %#ok<AGROW>
end

function Y = p0_applyc(X, xin, oc, pts)
% Outputs of container X on the per-component inputs XIN, every output
% component evaluated at its points: numel(oc) x q x P.
N = numel(X.dim_in);    M = numel(X.dim_out);
q = size(xin{1},2);     P = size(pts{1},1);
ej = [0; cumsum(X.dim_in(:))];  ei = [0; cumsum(X.dim_out(:))];
Y = zeros(numel(oc), q, P);
for i = 1:M
    acc = polynomial(zeros(X.dim_out(i), q));
    for j = 1:N
        B = X.C{i,j};
        if isempty(B),  continue,   end
        acc = acc + apply_sopvar(B, vertcat(xin{ej(j)+1:ej(j+1)}));
    end
    for r = 1:X.dim_out(i)
        c = ei(i) + r;
        Y(c,:,:) = p0_evalpoly(acc(r,:), oc{c}, pts{c});
    end
end
end
