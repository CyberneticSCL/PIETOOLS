function D = lb_settings_check()
% D = LB_SETTINGS_CHECK() compares the settings structs that lpisettings
% builds for the six tiers, separately for the 1-D part and for
% settings_2d, and prints which tiers coincide and in which fields the
% 'custom' tier differs from the others.
%
% OUTPUT
%   D:  struct with fields 'tiers' (cellstr), 'same1' and 'same2'
%       (logical 6x6: isequal of the 1-D parts / of settings_2d)

tiers = {'extreme','stripped','light','heavy','veryheavy','custom'};
S = cell(1,numel(tiers));
for k = 1:numel(tiers), S{k} = lpisettings(tiers{k},0,'','mosek'); end
n = numel(tiers);
same1 = false(n);  same2 = false(n);
for i = 1:n
    for j = 1:n
        a = rmfield(S{i},'settings_2d');  b = rmfield(S{j},'settings_2d');
        same1(i,j) = isequaln(a,b);
        same2(i,j) = isequaln(S{i}.settings_2d,S{j}.settings_2d);
    end
end
D = struct('tiers',{tiers},'same1',same1,'same2',same2);
fprintf('\nTier pairs with identical 1-D settings: ');
lb_pairs(tiers,same1);
fprintf('Tier pairs with identical 2-D settings: ');
lb_pairs(tiers,same2);
c = S{end};
for k = 1:n-1
    f1 = lb_diff(rmfield(c,'settings_2d'), rmfield(S{k},'settings_2d'));
    f2 = lb_diff(c.settings_2d, S{k}.settings_2d);
    fprintf('custom vs %-9s: 1-D fields differing: %s | 2-D: %s\n', tiers{k}, ...
        strjoin(f1,' '), strjoin(f2,' '));
end
end

function lb_pairs(tiers,same)
p = {};
for i = 1:numel(tiers)
    for j = i+1:numel(tiers)
        if same(i,j), p{end+1} = [tiers{i},'=',tiers{j}]; end %#ok<AGROW>
    end
end
if isempty(p), fprintf('none\n'); else, fprintf('%s\n',strjoin(p,', ')); end
end

function f = lb_diff(a,b)
% Top-level fields whose values differ (or exist in one struct only).
fa = fieldnames(a);  fb = fieldnames(b);
f = setdiff(union(fa,fb),intersect(fa,fb))';
for k = 1:numel(fa)
    if isfield(b,fa{k}) && ~isequaln(a.(fa{k}),b.(fa{k})), f{end+1} = fa{k}; end %#ok<AGROW>
end
if isempty(f), f = {'(none)'}; end
end
