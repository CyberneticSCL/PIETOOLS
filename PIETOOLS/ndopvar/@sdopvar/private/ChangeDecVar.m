function P = ChangeDecVar(P,Zd)
% P = ChangeDecVar(P,Zd)
% Change the decision-variable order of an sdopvar object.
% this routine does not change the operator P. It only rewrites it so that B rows match a new decision-variable ordering Zd,
% inserting zero rows for newly introduced decision variables.
%% First, convert the existing decision-variable list P.Zd and the desired decision-variable list Zd 
% into column string vectors for comparison and indexing.
Zd_old = P.Zd(:).';
Zd_new = Zd(:).';
% Checks the case where the decision variables already match exactly, including order.
if isequal(Zd_old,Zd_new)
    P.Zd = Zd_new;
    return
end
% Each old decision variable need to appear in the new list.
[tf,loc] = ismember(Zd_old,Zd_new);
if any(~tf)
    error('New decision-variable list must contain all existing variables.');
end
% old number of rows in each B{i}.
n_old = numel(Zd_old);
% new number of rows in each B{i}.
n_new = numel(Zd_new);
params = P.params;
for ii=1:numel(params.A)
    Aii = params.A{ii};
    Bii = params.B{ii};
    % case where no decision coefficient matrix was stored.
    if isempty(Bii)
        % Treats an empty B as a sparse matrix with zero decision-variable rows
        % and one column per scalar entry of Aii.
        Bii = sparse(0,numel(Aii));
    end
    block = numel(Aii);
    % Creates a new zero sparse B matrix with one row per new decision variable.
    Bnew = sparse(n_new,block);
    % Only remap old rows if there were old decision variables.
    if n_old>0
        % Moves the old B rows into their new locations.
        %Example: if old d2 is now row 2 and old d1 is now row 1 this line reorders the rows accordingly. Any newly added decision variables get zero rows.
        Bnew(loc,:) = reshape(Bii,n_old,block);
    end
    params.B{ii} = Bnew;
end
P.params = params;
P.Zd = Zd_new;
end
