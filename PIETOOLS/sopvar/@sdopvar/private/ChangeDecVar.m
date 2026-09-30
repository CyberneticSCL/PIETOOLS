function P = ChangeDecVar(P,Zd,loc)                                         % MMP, 09/25/2026
% function P = ChangeDecVar(P,Zd)                                           % MMP, 09/25/2026 (was)
% P = ChangeDecVar(P,Zd)
% Change the decision-variable order of an sdopvar object.
% this routine does not change the operator P. It only rewrites it so that B rows match a new decision-variable ordering Zd,
% inserting zero rows for newly introduced decision variables.
%
% MMP, 09/25/2026: Optional third input LOC, the position in Zd of each
%                  entry of P.Zd, for a caller that moves many blocks sharing
%                  one list and has the map already (the container
%                  concatenations get it free from one 'unique'). It skips
%                  the per-block 'ismember', which was 4.4 of 5.3 s of a
%                  two-list [D, E] at q = 4e5. LOC's size and range are
%                  checked, a scan of LOC, and three entries by name; a
%                  permutation within the right length and range is NOT
%                  caught, so LOC must come from the list itself. Without
%                  LOC the behaviour is unchanged.
% MMP, 09/29/2026: The row remap is now one product per cell, B_new = Pi*B,
%                  with Pi = T1' for the list change d_old = T1*d_new (spec
%                  sec. 8.1.4), i.e. Pi(loc(j),j) = 1, built once per call.
%                  It replaces a scatter into a preallocated sparse, the
%                  costliest statement of the class in a 3-D heat build.
%                  Measured on that build's two largest calls (n_old 1.07e6
%                  and 9.2e5, 27 cells, nnz 2.4e7 and 1.9e7), best of 3,
%                  with LOC: 0.87-0.94 -> 0.24-0.27 s and 0.73-0.75 ->
%                  0.21-0.23 s; allocated 373 -> 398 MB, the difference Pi.
%                  A supplied LOC must have distinct entries, which the probe
%                  does not check. A name repeated in P.Zd (the constructor
%                  does not forbid it) gives a repeated LOC entry on the
%                  name path too; Pi*B then sums its rows, which is C(d)
%                  read by name, where the scatter kept only the last row
%                  (wrong). B is made sparse before the product, since
%                  Pi*full(B) is full while the scatter into a sparse stayed
%                  sparse; LOC is made double, since sparse() refuses mixed
%                  index classes and the scatter accepted any numeric LOC.
%                  Memory: Pi, 24 B per old decision variable (240 MB at
%                  q = 1e7); peak memory inside the call rises by Pi when a
%                  B cell is smaller than Pi (34 -> 69 MB at q = 3e6). The
%                  no-LOC branch is re-indented (whitespace).
%% First, convert the existing decision-variable list P.Zd and the desired decision-variable list Zd 
% into column string vectors for comparison and indexing.
Zd_old = P.Zd(:).';
Zd_new = Zd(:).';
if nargin<3 || isempty(loc)                                                 % MMP, 09/25/2026
    % Checks the case where the decision variables already match exactly, including order. % MMP, 09/29/2026
    if isequal(Zd_old,Zd_new)                                               % MMP, 09/29/2026
        P.Zd = Zd_new;                                                      % MMP, 09/29/2026
        return                                                              % MMP, 09/29/2026
    end                                                                     % MMP, 09/29/2026
    % Each old decision variable need to appear in the new list.            % MMP, 09/29/2026
    [tf,loc] = ismember(Zd_old,Zd_new);                                     % MMP, 09/29/2026
    if any(~tf)                                                             % MMP, 09/29/2026
        error('New decision-variable list must contain all existing variables.'); % MMP, 09/29/2026
    end                                                                     % MMP, 09/29/2026
else                                                                        % MMP, 09/25/2026
    % A supplied map is trusted for the names - that is the saving - but a  % MMP, 09/25/2026
    % map built for another list is caught by its size, its range, and a    % MMP, 09/25/2026
    % name check at three positions.                                        % MMP, 09/25/2026
    loc = loc(:).';                                                         % MMP, 09/25/2026
    nl = numel(Zd_old);                                                     % MMP, 09/25/2026
    probe = unique([1,ceil(nl/2),nl]);                                      % MMP, 09/25/2026
    if numel(loc)~=nl || any(loc<1) || any(loc>numel(Zd_new)) ||...
            ~isequal(Zd_new(loc(probe)),Zd_old(probe))                      % MMP, 09/25/2026
        error('ChangeDecVar:badLoc',['LOC must give, for each entry of '...
            'P.Zd, its position in Zd.']);                                  % MMP, 09/25/2026
    end                                                                     % MMP, 09/25/2026
end                                                                         % MMP, 09/25/2026
% old number of rows in each B{i}.
n_old = numel(Zd_old);
% new number of rows in each B{i}.
n_new = numel(Zd_new);
params = P.params;
% Pi = T1' (spec sec. 8.1.4): row loc(j) of Pi*B is row j of B, the other   % MMP, 09/29/2026
% rows are zero. Exact: each output entry is 1 times one input entry.       % MMP, 09/29/2026
if n_old>0,     Pi = sparse(double(loc(:)),(1:n_old).',1,n_new,n_old); end  % MMP, 09/29/2026
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
%   % Creates a new zero sparse B matrix with one row per new decision variable. % MMP, 09/29/2026 (was)
%   Bnew = sparse(n_new,block);                                             % MMP, 09/29/2026 (was)
    % Only remap old rows if there were old decision variables.
    if n_old>0
        % Moves the old B rows into their new locations.
        %Example: if old d2 is now row 2 and old d1 is now row 1 this line reorders the rows accordingly. Any newly added decision variables get zero rows.
%       Bnew(loc,:) = reshape(Bii,n_old,block);                             % MMP, 09/29/2026 (was)
        Bii = reshape(Bii,n_old,block);                                     % MMP, 09/29/2026
        % Pi*full(B) would be full; the scatter it replaces stayed sparse.  % MMP, 09/29/2026
        if ~issparse(Bii),  Bii = sparse(Bii);  end                         % MMP, 09/29/2026
        params.B{ii} = Pi*Bii;                                              % MMP, 09/29/2026
    else                                                                    % MMP, 09/29/2026
        params.B{ii} = sparse(n_new,block);                                 % MMP, 09/29/2026
    end
%   params.B{ii} = Bnew;                                                    % MMP, 09/29/2026 (was)
end
P.params = params;
P.Zd = Zd_new;
end
