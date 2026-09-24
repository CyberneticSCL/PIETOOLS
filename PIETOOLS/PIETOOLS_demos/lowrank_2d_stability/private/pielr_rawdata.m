function D = pielr_rawdata(sos)                                             % CC, 09/23/2026
% PIELR_RAWDATA  Extract the SeDuMi data of an LPI program: the EQUALITY rows,
% the objective, and the x-vector layout.  Generalises raw_data, which is
% correct only for the stability programs.
%
% TWO DEFECTS IN raw_data THIS FIXES, both measured on the 1-D l2gain program
% (PIETOOLS_Hinf_gain lines, Ex_Transport_Eq_with_Disturbance, 'light'):
%
% 1. IT CONCATENATES INEQUALITY EXPRESSIONS AS EQUALITIES.  raw_data loops
%    i = 1:sos.expr.num and appends every expr.At{i}/expr.b{i}.  On that
%    program expr 1 has sos.expr.type{1} == 'ineq' -- it is the gam >= 0
%    constraint from lpi_ineq(prog,gam).  Treating it as an equality row PINS
%    gamma to a value instead of bounding it, which silently changes the
%    problem.  Here 'eq' rows build A(x) = b and 'ineq' rows are returned
%    separately in D.Gineq / D.hineq for the caller to impose as bounds.
%
% 2. IT ASSUMES A CONTIGUOUS BLOCK LAYOUT.  See pielr_layout's header: with
%    lpivar in the program the free coordinates are interleaved between the
%    Gram blocks, and raw_data loses them.  Layout comes from pielr_layout.
%
% THE OBJECTIVE.  raw_data returns no objective at all, so an optimisation
% executive handed to it degrades silently into a feasibility problem.
% sossolve.m:217-222 shows the objective vector is sos.objective, added into
% c over 1:sos.var.idx{end}-1; it is returned here as D.c.
%
% OUTPUT (struct)
%   D.At      Ntot x m, the equality rows (SeDuMi's At; constraint At'*x = b)
%   D.b       m x 1
%   D.c       Ntot x 1 objective (all zero when the program has none)
%   D.Gineq   Ntot x mi, inequality rows (At'*x >= b, SOSTOOLS' convention)
%   D.hineq   mi x 1
%   D.L       pielr_layout(sos)
%   D.has_obj true iff any(D.c ~= 0)

L = pielr_layout(sos);

Ae = {};  be = {};  Ai = {};  bi = {};
for i = 1:sos.expr.num
    if isfield(sos.expr,'type') && numel(sos.expr.type)>=i ...
            && strcmp(sos.expr.type{i},'ineq')
        Ai{end+1} = sos.expr.At{i};   bi{end+1} = sos.expr.b{i}; %#ok<AGROW>
    else
        Ae{end+1} = sos.expr.At{i};   be{end+1} = sos.expr.b{i}; %#ok<AGROW>
    end
end
% one horzcat/vertcat of the whole list, not a growing concatenation: the row
% dimension here is the decision-variable count and a loop would repeatedly
% copy it (see CLAUDE.md section 2)
D.At = [Ae{:}];
D.b  = full(vertcat(be{:}));
if isempty(Ai)
    D.Gineq = sparse(L.Ntot,0);   D.hineq = zeros(0,1);
else
    D.Gineq = [Ai{:}];            D.hineq = full(vertcat(bi{:}));
end

c = zeros(L.Ntot,1);
if isfield(sos,'objective') && ~isempty(sos.objective)
    no = numel(sos.objective);
    if no > L.Ntot
        error('pielr_rawdata:obj','objective longer (%d) than x (%d)',no,L.Ntot);
    end
    c(1:no) = full(sos.objective(:));
end
D.c = c;
D.has_obj = any(c~=0);
D.L = L;

% The layout must account for every coordinate; a silent mismatch is exactly
% the failure mode this routine exists to remove, so it is a hard error.
nacc = numel(L.free) + sum(cellfun(@numel,L.rows));
if nacc ~= L.Ntot
    error('pielr_rawdata:layout', ...
          'layout accounts for %d of %d coordinates',nacc,L.Ntot);
end
end
