function [prog,P] = poslpivar_settings_2d_sop(prog,X,st,role,side,PIE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P] = POSLPIVAR_SETTINGS_2D_SOP(PROG,X,ST,ROLE,SIDE,PIE) declares
% the positive operator P = P0 + P1 + ... that the 2-D LPI settings in ST
% define for ROLE on the SIDE spaces of the container X, each term a
% positive semidefinite 'cdopvar' Z'*(g*Q)*Z, Q >= 0, declared by
% POSCOPVAR: Z the identity on an R^n space and, on an L2[x,y] space, the
% basis operators Z2^{ij} of the settings that are kept, with their degree
% caps, as SETTINGS2POSSOPVAR translates them; g the psatz weight of the
% term's options:
%   ROLE 'lf'     P0 of LF_deg and LF_opts; one Pj of LF_deg_psatz{j} and
%                 LF_opts_psatz{j} per nonzero LF_use_psatz(j);
%   ROLE 'slack'  P0 of eq_deg and eq_opts; one Pj of eq_deg_psatz{j} and
%                 eq_opts_psatz{j} per nonzero eq_use_psatz(j). The exclude
%                 and sep entries of every term are or-ed with those the
%                 structure of X allows: exclude(1) when every R^n -> R^n
%                 block of X is zero, the exclude entry of Z2^{ij} when
%                 parameter cell (i,j) of every L2[x,y] -> L2[x,y] block of
%                 X is zero, and sep(3:6) when the lower- and upper-integral
%                 cells that the entry identifies agree, to 1e-12, in those
%                 blocks.
% No constraint other than the positivity of each term is added.
%
% INPUT
% - prog:   'struct', the LPI program to which the decision variables are
%           added;
% - X:      'copvar' or 'cdopvar' object whose SIDE spaces P maps, each R^n
%           or L2[x,y] (in the two variables of the PIE);
% - st:     'struct' of LPI settings: an LPISETTINGS struct with field
%           settings_2d, or a 2-D settings struct (settings_PIETOOLS_*_2D,
%           field is2D = 1). Fields read: LF_deg, LF_opts, LF_use_psatz,
%           LF_deg_psatz, LF_opts_psatz ('lf'); eq_deg, eq_opts,
%           eq_use_psatz, eq_deg_psatz, eq_opts_psatz, use_sosineq
%           ('slack'). Each *_deg is a 'struct' with field d2, the 3x3
%           'cell' whose entry {i,j} is the array of degree caps of Z2^{ij}
%           (one cap per subset of its variables; i, j = 1 multiplier, 2
%           lower, 3 upper integral in x and in y); each *_opts a 'struct'
%           with fields exclude (1x16), sep (1x6) and psatz; each
%           *_use_psatz a 'double' array, nonzero for a psatz term in use;
%           *_deg_psatz and *_opts_psatz a 'cell' with one entry per entry
%           of *_use_psatz, one entry (or a 'struct') for every term, or a
%           longer 'cell' whose entry *_use_psatz(j) serves term j;
% - role:   'char', 'lf' or 'slack';
% - side:   'char', 'out' (the output spaces of X) or 'in' (its input
%           spaces);
% - PIE:    'pie_struct' object of a 2-D PIE: PIE.vars(:,1) ('polynomial',
%           the two spatial variables) and PIE.dom (2x2 'double') give the
%           domain of P.
%
% OUTPUT
% - prog:   'struct', the program with the decision variables and Gram
%           constraints of P0, P1, ... added in that order;
% - P:      'cdopvar' object, the sum of the terms.
%
% NOTES
% Errors: ROLE 'slack' with use_sosineq set; a ROLE other than 'lf' and
% 'slack'; a space of X other than R^n and L2[x,y]; exclude(1) set (for
% 'slack' also from a zero R^n block of X) with an R^n space present; a
% term with psatz = 2, which POSCOPVAR rejects. A part of the settings the
% translation cannot carry (sep(2+k) ~= sep(4+k) in direction k, taken as
% not separable; a missing d2 cell, whose basis operator is omitted) is the
% warning 'poslpivar_settings_2d_sop:lossy'.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - poslpivar_settings_2d_sop
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 10/08/2026
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body; code unchanged. The deleted header text was
%                the provenance of the code (the test-folder translators it
%                transcribes) and comparisons with the stock 2-D routines.

% The 2-D settings, the domain of P (variable names of the PIE), and the
% SIDE spaces of X
S = settings_2d_sop(st);
dom = struct();     dom.vars = reshape(PIE.vars(:,1).varname,1,[]);    dom.dom = PIE.dom;
[sp,dm] = copvar_space_list(X,side);
switch role
    case 'lf'
        % P0, then one term per psatz entry in use
        [deg,co] = pos2d(S.LF_deg,S.LF_opts,sp);
        [prog,P] = poscopvar(prog,dm,sp,dom,deg,co);
        for j = 1:numel(S.LF_use_psatz)
            if S.LF_use_psatz(j)~=0
                [deg,co] = pos2d(S.LF_deg_psatz{j},S.LF_opts_psatz{j},sp);
                [prog,P2] = poscopvar(prog,dm,sp,dom,deg,co);
                P = P + P2;
            end
        end
    case 'slack'
        % use_sosineq (lpi_ineq_2d in place of the equality) is not
        % supported
        if S.use_sosineq
            error('poslpivar_settings_2d_sop:sosineq','use_sosineq selects lpi_ineq_2d, which has no container counterpart.');
        end
        % The exclusions and separations the structure of X allows; P0,
        % then one term per psatz entry in use, with the same exclusions and
        % separations added to its own
        eq_opts = eqopts2d(X,S.eq_opts,1e-12);
        [deg,co] = pos2d(S.eq_deg,eq_opts,sp);
        [prog,P] = poscopvar(prog,dm,sp,dom,deg,co);
        for j = 1:numel(S.eq_use_psatz)
            if S.eq_use_psatz(j)~=0
                o = S.eq_opts_psatz{j};
                o.exclude = o.exclude | eq_opts.exclude;    o.sep = o.sep | eq_opts.sep;
                [deg,co] = pos2d(S.eq_deg_psatz{j},o,sp);
                [prog,P2] = poscopvar(prog,dm,sp,dom,deg,co);
                P = P + P2;
            end
        end
    otherwise
        error('poslpivar_settings_2d_sop:role','ROLE should be ''lf'' or ''slack''.');
end
end


function [deg,copts] = pos2d(pdeg,popts,sp)
% [DEG,COPTS] = POS2D(PDEG,POPTS,SP) returns the POSCOPVAR degrees and
% options of one term with the 2-D settings degrees PDEG and options POPTS
% over the spaces SP: the identity basis on R^n, the basis operators of
% SETTINGS2POSSOPVAR on L2[x,y].
% INPUT  pdeg:   'struct' with field d2 (3x3 'cell' of degree-cap arrays);
%        popts:  'struct' with fields exclude (1x16), sep (1x6), psatz;
%        sp:     1xM 'cell' of 'cellstr', the variable names of each space.
% OUTPUT deg:    1xM 'cell', the degree specification of each space;
%        copts:  'struct' with fields psatz (POPTS.psatz, default 0), sep
%                (1x2 logical) and include (1xM 'cell' of multi-indices).
% The translation (as the LF block of a 2-D settings struct); what it
% cannot carry for that block is a warning
s = struct('is2D',1);   s.LF_deg = pdeg;    s.LF_opts = popts;
o = settings2possopvar(s);
lossy = o.report.lossy;
lossy = lossy(startsWith(lossy,'LF'));
if ~isempty(lossy)
    warning('poslpivar_settings_2d_sop:lossy','%s',strjoin(lossy,' '));
end
% exclude, padded to 16 entries
exc = zeros(1,16);
if isfield(popts,'exclude') && ~isempty(popts.exclude)
    exc(1:numel(popts.exclude)) = popts.exclude;
end
% Per space: the identity basis on R^n (an error under exclude(1)), the
% translated basis operators on L2[x,y], an error otherwise
deg = cell(1,numel(sp));    incl = cell(1,numel(sp));
for k = 1:numel(sp)
    v = sort(sp{k});
    if isempty(v)
        if exc(1)
            error('poslpivar_settings_2d_sop:exclude1',['exclude(1) drops the R^n basis; '...
                  'poscopvar needs one basis operator per space.']);
        end
        deg{k} = struct('int',[0 0]);
    elseif numel(v)==2
        deg{k} = o.LF.deg;      incl{k} = o.LF.opts.include;
    else
        error('poslpivar_settings_2d_sop:space',['Space {%s}: only R^n and L2[x,y] spaces '...
              'are handled.'],strjoin(v,','));
    end
end
% psatz of POPTS unchanged, sep of the translation
ps = 0;     if isfield(popts,'psatz') && ~isempty(popts.psatz),  ps = popts.psatz;  end
copts = struct('psatz',ps,'sep',o.LF.opts.sep);
copts.include = incl;           % assigned apart: a cell value in struct() builds a struct array
end


function eq_opts = eqopts2d(Km,eq_opts,ztol)
% EQ_OPTS = EQOPTS2D(KM,EQ_OPTS,ZTOL) sets the exclude and sep entries of
% EQ_OPTS that the structure of KM allows, for R^n and L2[x,y] spaces.
% With (i,j) the parameter cells of an L2[x,y] -> L2[x,y] block (i, j = 1
% multiplier, 2 lower, 3 upper integral in x, in y), and each condition
% holding in every such block of KM:
%   exclude(1)     every R^n -> R^n block zero;
%   exclude(8:16)  cell (1,1), (2,1), (3,1), (1,2), (1,3), (2,2), (3,2),
%                  (2,3), (3,3) zero;
%   sep(3)         (2,1) - (3,1) <= ZTOL, exclude(9:10) unset;
%   sep(4)         (1,2) - (1,3) <= ZTOL, exclude(11:12) unset;
%   sep(5)         (2,2) - (3,2) <= ZTOL and (2,3) - (3,3) <= ZTOL,
%                  exclude(13:14) or exclude(15:16) unset;
%   sep(6)         (2,2) - (2,3) <= ZTOL and (3,2) - (3,3) <= ZTOL,
%                  exclude([13 15]) or exclude([14 16]) unset;
% zero meaning no coefficient above 1e-12 in magnitude, and a difference
% <= ZTOL holding in every coefficient.
% INPUT  Km:      'copvar' or 'cdopvar' object with the same spaces in and
%                 out (the output spaces are read);
%        eq_opts: 'struct' with fields exclude (1x16) and sep (1x6);
%        ztol:    'double', the tolerance of the differences.
% OUTPUT eq_opts: EQ_OPTS with those entries set, the others unchanged.
exc = eq_opts.exclude;      sep = eq_opts.sep;
% Spaces: R^n (no variable) and L2[x,y] (two variables) only
nvs = sum(Km.space_out,2);
if any(nvs==1) || numel(Km.vars)>2
    error('poslpivar_settings_2d_sop:space','Only R^n and L2[x,y] spaces are handled.');
end
iR = find(nvs==0);      iX = find(nvs==2);
% exclude(1): every R^n -> R^n block zero
if ~isempty(iR) && all_zero(Km,iR,1)
    exc(1) = 1;
end
if ~isempty(iX)
    % exclude(8:16): pidx(ii) is the exclude entry of the cell of linear
    % index ii in the 3x3 grid of cells
    pidx = [8,9,10,11,13,14,12,15,16];
    for ii = 1:9
        if all_zero(Km,iX,ii),  exc(pidx(ii)) = 1;  end
    end
    % sep(3:6): the cell differences of the header
    c = @(g1,g2) sub2ind([3,3],g1,g2);
    if ~exc(9) && ~exc(10) && diff_le(Km,iX,c(2,1),c(3,1),ztol)
        sep(3) = 1;
    end
    if ~exc(11) && ~exc(12) && diff_le(Km,iX,c(1,2),c(1,3),ztol)
        sep(4) = 1;
    end
    if ((~exc(13) && ~exc(14)) || (~exc(15) && ~exc(16))) && ...
            diff_le(Km,iX,c(2,2),c(3,2),ztol) && diff_le(Km,iX,c(2,3),c(3,3),ztol)
        sep(5) = 1;
    end
    if ((~exc(13) && ~exc(15)) || (~exc(14) && ~exc(16))) && ...
            diff_le(Km,iX,c(2,2),c(2,3),ztol) && diff_le(Km,iX,c(3,2),c(3,3),ztol)
        sep(6) = 1;
    end
end
eq_opts.exclude = exc;      eq_opts.sep = sep;
end


function tf = all_zero(Km,idx,q)
% TF = ALL_ZERO(KM,IDX,Q) is true when parameter cell Q of every block
% (i,j) of KM, i and j in IDX, has no coefficient above 1e-12 in magnitude
% (an empty block counts as zero).
% INPUT  Km:  'copvar' or 'cdopvar' object;
%        idx: 'double' array, space indices;
%        q:   'double', the parameter cell.
% OUTPUT tf:  'logical'.
% Over every block (i,j); the first nonzero cell ends the search
tf = true;
for i = idx(:)',    for j = idx(:)'
    [A,Bq] = cellAB(Km.C{i,j},q,1e-12);
    if any(A) || nnz(Bq),   tf = false;     return,     end
end,                end
end


function tf = diff_le(Km,idx,q1,q2,ztol)
% TF = DIFF_LE(KM,IDX,Q1,Q2,ZTOL) is true when, in every block (i,j) of KM
% with i and j in IDX, parameter cell Q1 minus cell Q2 is at most ZTOL in
% every coefficient, constant and decision part (in each cell the
% coefficients up to ZTOL in magnitude taken as zero).
% INPUT  Km:     'copvar' or 'cdopvar' object;
%        idx:    'double' array, space indices;
%        q1, q2: 'double', parameter cells;
%        ztol:   'double', the tolerance.
% OUTPUT tf:     'logical'.
% Over every block (i,j); the first difference above ZTOL ends the search
tf = true;
for i = idx(:)',    for j = idx(:)'
    [A1,B1] = cellAB(Km.C{i,j},q1,ztol);    [A2,B2] = cellAB(Km.C{i,j},q2,ztol);
    dB = B1 - B2;
    if any(A1-A2 > ztol) || any(nonzeros(dB) > ztol)    % implicit zeros pass
        tf = false;     return
    end
end,                end
end


function [A,Bq] = cellAB(B,q,ztol)
% [A,BQ] = CELLAB(B,Q,ZTOL) returns parameter cell Q of the block B as its
% constant part A and decision part BQ, a zero shorthand (0 or []) of the
% cell expanded to full size and coefficients up to ZTOL in magnitude set
% to zero.
% INPUT  B:    'sopvar' or 'sdopvar' object, or [] (a zero block);
%        q:    'double', the parameter cell;
%        ztol: 'double', the zero tolerance.
% OUTPUT A:    'double' column, the constant coefficients (0x1 for []);
%        Bq:   sparse 'double', the decision coefficients, one row per
%              decision variable (no rows for a 'sopvar' block).
% An empty block: empty parts
A = zeros(0,1);     Bq = sparse(0,0);
if isempty(B),  return,     end
% Coefficient count of a cell: (dims(1) x left monomials) x (dims(2) x
% right monomials)
nL = prod([cellfun(@numel,B.ZL),1]);    nR = prod([cellfun(@numel,B.ZR),1]);
mn = B.dims(1)*nL*B.dims(2)*nR;
if isa(B,'sdopvar')
    A = B.params.A{q};      Bq = B.params.B{q};
    if numel(A)~=mn,        A = zeros(mn,1);    end         % 0 / [] shorthand
    if size(Bq,2)~=mn,      Bq = sparse(numel(B.Zd),mn);    end
else                                                        % fixed sopvar block
    A = B.params{q};        Bq = sparse(0,mn);
    if numel(A)~=mn,        A = zeros(mn,1);    end
end
% Coefficients up to ztol in magnitude to zero
A = full(A(:));     A(abs(A)<=ztol) = 0;
Bq = Bq.*(abs(Bq)>ztol);
end
