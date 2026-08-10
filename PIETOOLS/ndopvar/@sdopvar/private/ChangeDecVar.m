function P = ChangeDecVar(P,Zd)
% P = ChangeDecVar(P,Zd)
% Change the decision-variable order of an sdopvar object.
% this routine does not change the operator P. It only rewrites it so that Bt colummns match a new decision-variable ordering Zd,
% inserting zero colummns for newly introduced decision variables.
%% First, convert the existing decision-variable list P.Zd and the desired decision-variable list Zd 
% into column cell arrays of character vectors.
Zd_old = cellstr(string(P.Zd(:)));
Zd_new = cellstr(string(Zd(:)));
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
% old number of columms in each Bt{i}.
n_old = numel(Zd_old);
% new number of columms in each Bt{i}.
n_new = numel(Zd_new);
params = P.params;
for ii=1:numel(params.A)
    Aii = params.A{ii};
    Bii = params.Bt{ii};
    % case where no decision coefficient matrix was stored.
    if isempty(Bii)
        %Treats an empty Bt as a sparse matrix with zero decision-variable collumns
        % and one row per scalar entry of Aii.
        Bii = sparse(numel(Aii),0);
    end
    block = numel(Aii);
    % Creates a new zero sparse Bt matrix with one column per new decision variable.
    Bnew = sparse(block,n_new);
    % Only remap old rows if there were old decision variables.
    if n_old>0
        % Moves the old Bt rows into their new locations.
        %Example: if old d2 is now row 2 and old d1 is now row 1 this line reorders the rows accordingly. Any newly added decision variables get zero rows.
        Bnew(:,loc) = reshape(Bii,block,n_old);
    end
    params.Bt{ii} = Bnew;
end
P.params = params;
P.Zd = Zd_new;
end
