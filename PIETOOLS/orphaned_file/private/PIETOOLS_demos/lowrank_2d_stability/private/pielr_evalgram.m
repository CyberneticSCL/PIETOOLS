function P = pielr_evalgram(prog,P)
% PIELR_EVALGRAM -- substitute a solved decision vector into a dopvar2d,
% directly, with sparse algebra only.
%
%   P = pielr_evalgram(prog,P)
%
% Returns the opvar2d obtained by evaluating every dpvar parameter cell of
% the dopvar2d P at the decision values prog.solinfo.RRx (ordered as
% prog.decvartable) -- the same operator lpigetsol/getsol_lpivar_2d returns,
% without going through sosgetsol.
%
% WHY THIS EXISTS.  On PIETOOLS master (998a68cc) sosgetsol handles a dpvar
% by dpvar2poly + combine, which expands every decision variable into a
% polynomial variable: on this package's negativity operator (~277k dvars,
% heat n=1) getsol_lpivar_2d did not finish in 35 minutes there, vs ~1 s on
% the development branch.  opcheck_2d / opcheck_pos2d were the package's
% only getsol dependence, so verification substitutes here instead and the
% package runs on both trees.
%
% LAYOUT (verified against the live class 2026-08-28; the dpvar classdef
% header claims the opposite and is WRONG): dpvar(C,degmat,varname,dvarname,
% [m,n]) stores C with ROWS matrix-row outer, decision-variable inner -- one
% block of nd+1 consecutive rows per matrix row, constant term FIRST, i.e. a
% left factor kron(I_m,[1 d']) -- and COLUMNS matrix-column outer, monomial
% inner.  Evaluation at d is therefore ONE sparse product per cell,
%     Cnew = kron(speye(m),sparse([1;d]')) * C          (m x n*nZ)
% followed by an index remap into polynomial's T x (m*n) column-major
% coefficient layout.  Nothing is ever full()'d on a decision-variable
% dimension (CLAUDE.md performance rules; rows of C are the ~q axis).
%
% Each cell's dvarname (a subset of the program's dvars, in its own order)
% is mapped to the global vector BY NAME against prog.decvartable, whose
% positions index RRx directly -- the same mapping sosgetsol's own dpvar
% path uses.  A dvar missing from decvartable is an ERROR, never a symbolic
% leftover: in this package every dvar comes from the program, so a miss
% means the wrong program/operator pairing, and hiding it would corrupt the
% operator-level verification this feeds.
%
% INPUT
%   prog - LPI/SOS program struct; only .decvartable (cellstr) and
%          .solinfo.RRx (decision vector, decvartable order) are used
%   P    - dopvar2d to evaluate; an opvar2d passes through unchanged
% OUTPUT
%   P    - opvar2d with every cell evaluated (polynomial or double)
%
% The 36-leaf enumeration below is getsol_lpivar_2d's own; if the class ever
% grows a parameter field, opcheck's 16-field/36-leaf assert on this
% function's output catches it.
%
% CC, 09/19/2026: initial version.
% CC, 09/23/2026: the local function evaldp moved out to private/pielr_evaldp.m
%   unchanged, because pielr_evalop needs the identical per-cell arithmetic for
%   the 6 leaves of a 1-D dopvar.  Two real callers, and a second copy would be
%   free to drift; nothing else about this file changed.

if isa(P,'opvar2d')
    return                              % nothing to substitute
end
if ~isa(P,'dopvar2d')
    error('pielr_evalgram:input','input must be a dopvar2d (or opvar2d) object.');
end
dtab = prog.decvartable;
RRx  = prog.solinfo.RRx;

for f = {'R00','R0x','R0y','R02','Rx0','Rxy','Ry0','Ryx','R20'}
    P.(f{:}) = pielr_evaldp(P.(f{:}),dtab,RRx);
end
for i = 1:3
    P.Rxx{i,1} = pielr_evaldp(P.Rxx{i,1},dtab,RRx);
    P.Rx2{i,1} = pielr_evaldp(P.Rx2{i,1},dtab,RRx);
    P.R2x{i,1} = pielr_evaldp(P.R2x{i,1},dtab,RRx);
    P.Ryy{1,i} = pielr_evaldp(P.Ryy{1,i},dtab,RRx);
    P.Ry2{1,i} = pielr_evaldp(P.Ry2{1,i},dtab,RRx);
    P.R2y{1,i} = pielr_evaldp(P.R2y{1,i},dtab,RRx);
    for j = 1:3
        P.R22{i,j} = pielr_evaldp(P.R22{i,j},dtab,RRx);
    end
end
P = opvar2d(P);                         % dopvar2d with polynomial cells -> opvar2d
end
% CC, 09/23/2026: the local function evaldp that stood here was moved verbatim
% to private/pielr_evaldp.m and the six call sites above renamed to match.  It
% is not reproduced here because it is unchanged -- see that file's header.
