function [prog,P] = poslpivar_settings_2d_sop(prog,X,st,role,side,PIE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P] = POSLPIVAR_SETTINGS_2D_SOP(PROG,X,ST,ROLE,SIDE,PIE) declares
% the positive operator P = P0 + P1 + ... that the 2-D LPI settings in ST
% define for ROLE on the SIDE spaces of the container X, each term a
% positive semidefinite 'cdopvar' Z'*(g*Q)*Z, Q >= 0, declared by
% POSCOPVAR: Z the identity on an R^n space and, on an L2[x,y] space, the
% basis operators Z2^{ij} of the settings that are kept, with their degree
% caps, as SETTINGS2POSSOPVAR translates them; g the psatz weight of the
% term's options. On an L2[x] (L2[y]) space Z is the bases Zx^o, Zx^a,
% Zx^b (Zy^o, Zy^a, Zy^b) of the degrees dx (dy) that are kept, by the 1-D
% translation of SETTINGS2POSSOPVAR, with exclude(2:4) and sep(1)
% (exclude(5:7) and sep(2)):
%   ROLE 'lf'     P0 of LF_deg and LF_opts; one Pj of LF_deg_psatz{j} and
%                 LF_opts_psatz{j} per nonzero LF_use_psatz(j);
%   ROLE 'slack'  P0 of eq_deg and eq_opts; one Pj of eq_deg_psatz{j} and
%                 eq_opts_psatz{j} per nonzero eq_use_psatz(j). The exclude
%                 and sep entries of every term are or-ed with those the
%                 structure of X allows: exclude(1) when every R^n -> R^n
%                 block of X is zero; exclude(2:4) (exclude(5:7)) when the
%                 multiplier, lower- or upper-integral cell of every L2[x]
%                 -> L2[x] (L2[y] -> L2[y]) block of X is zero, and sep(1)
%                 (sep(2)) when its upper-integral cell is nonzero and its
%                 lower- and upper-integral cells agree, to 1e-12; the
%                 exclude entry of Z2^{ij} when parameter cell (i,j) of
%                 every L2[x,y] -> L2[x,y] block of X is zero, and sep(3:6)
%                 when the lower- and upper-integral cells that the entry
%                 identifies agree, to 1e-12, in those blocks.
% No constraint other than the positivity of each term is added.
%
% INPUT
% - prog:   'struct', the LPI program to which the decision variables are
%           added;
% - X:      'copvar' or 'cdopvar' object whose SIDE spaces P maps, each R^n,
%           L2[x], L2[y] or L2[x,y] (x, y the two variables of the PIE);
% - st:     'struct' of LPI settings: an LPISETTINGS struct with field
%           settings_2d, or a 2-D settings struct (settings_PIETOOLS_*_2D,
%           field is2D = 1). Fields read: LF_deg, LF_opts, LF_use_psatz,
%           LF_deg_psatz, LF_opts_psatz ('lf'); eq_deg, eq_opts,
%           eq_use_psatz, eq_deg_psatz, eq_opts_psatz, use_sosineq
%           ('slack'). Each *_deg is a 'struct' with fields dx (3x1
%           'cell': the degree in x of Zx^o, and the caps [x; theta; joint]
%           of Zx^a and Zx^b), dy (1x3 'cell', the same in y), read when
%           an L2[x] or L2[y] space is present, and d2, the 3x3
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
% 'slack'; a space of X in a variable other than the two of the PIE;
% exclude(1) set (for 'slack' also from a zero R^n block of X) with an R^n
% space present; every basis of an L2[x] or L2[y] space excluded; a term
% with psatz = 2, which POSCOPVAR rejects. A part of the settings the
% translation cannot carry is the warning 'poslpivar_settings_2d_sop:lossy':
% sep(2+k) ~= sep(4+k) in direction k, taken as not separable; sep(k)
% different from sep(2+k) and sep(4+k) with both an L2[x,y] space and a
% space in direction k alone, direction k then not separable (POSCOPVAR
% takes one flag per direction for every space); a missing d2 cell, whose
% basis operator is omitted.
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
% MMP, 10/28/2026: L2[x] and L2[y] spaces (one-variable states, as in the
%                PIEs of coupled ODE-PDE and delay systems). Before, such a
%                space was an error in pos2d and eqopts2d, so the container
%                executives could not run those 2-D PIEs (example library
%                PDE 42, 43, 50). pos2d declares their bases from dx, dy by
%                the 1-D translation of settings2possopvar (new local
%                function pos1d), and settles the one separation flag per
%                direction that poscopvar takes; eqopts2d adds the
%                exclude(2:7), sep(1:2) rules of get_eq_opts_2D. The help
%                text describes the new spaces; R^n and L2[x,y] spaces are
%                translated as before.

% The 2-D settings, the domain of P (variable names of the PIE), and the
% SIDE spaces of X
S = settings_2d_sop(st);
dom = struct();     dom.vars = reshape(PIE.vars(:,1).varname,1,[]);    dom.dom = PIE.dom;
[sp,dm] = copvar_space_list(X,side);
switch role
    case 'lf'
        % P0, then one term per psatz entry in use; pos2d reads the PIE's   % MMP, 10/28/2026
        % variable names to tell an L2[x] space from an L2[y] space         % MMP, 10/28/2026
%       [deg,co] = pos2d(S.LF_deg,S.LF_opts,sp);                            % MMP, 10/28/2026 (was)
        [deg,co] = pos2d(S.LF_deg,S.LF_opts,sp,dom.vars);                   % MMP, 10/28/2026
        [prog,P] = poscopvar(prog,dm,sp,dom,deg,co);
        for j = 1:numel(S.LF_use_psatz)
            if S.LF_use_psatz(j)~=0
%               [deg,co] = pos2d(S.LF_deg_psatz{j},S.LF_opts_psatz{j},sp);  % MMP, 10/28/2026 (was)
                [deg,co] = pos2d(S.LF_deg_psatz{j},S.LF_opts_psatz{j},sp,dom.vars); % MMP, 10/28/2026
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
%       eq_opts = eqopts2d(X,S.eq_opts,1e-12);                              % MMP, 10/28/2026 (was)
%       [deg,co] = pos2d(S.eq_deg,eq_opts,sp);                              % MMP, 10/28/2026 (was)
        eq_opts = eqopts2d(X,S.eq_opts,1e-12,dom.vars);                     % MMP, 10/28/2026
        [deg,co] = pos2d(S.eq_deg,eq_opts,sp,dom.vars);                     % MMP, 10/28/2026
        [prog,P] = poscopvar(prog,dm,sp,dom,deg,co);
        for j = 1:numel(S.eq_use_psatz)
            if S.eq_use_psatz(j)~=0
                o = S.eq_opts_psatz{j};
                o.exclude = o.exclude | eq_opts.exclude;    o.sep = o.sep | eq_opts.sep;
%               [deg,co] = pos2d(S.eq_deg_psatz{j},o,sp);                   % MMP, 10/28/2026 (was)
                [deg,co] = pos2d(S.eq_deg_psatz{j},o,sp,dom.vars);          % MMP, 10/28/2026
                [prog,P2] = poscopvar(prog,dm,sp,dom,deg,co);
                P = P + P2;
            end
        end
    otherwise
        error('poslpivar_settings_2d_sop:role','ROLE should be ''lf'' or ''slack''.');
end
end


% function [deg,copts] = pos2d(pdeg,popts,sp)                               % MMP, 10/28/2026 (was)
function [deg,copts] = pos2d(pdeg,popts,sp,xy)                              % MMP, 10/28/2026
% [DEG,COPTS] = POS2D(PDEG,POPTS,SP,XY) returns the POSCOPVAR degrees and
% options of one term with the 2-D settings degrees PDEG and options POPTS
% over the spaces SP: the identity basis on R^n, the bases of POS1D on
% L2[x] and L2[y], the basis operators of SETTINGS2POSSOPVAR on L2[x,y].
% INPUT  pdeg:   'struct' with fields dx, dy (read for an L2[x], L2[y]
%                space) and d2 (3x3 'cell' of degree-cap arrays);
%        popts:  'struct' with fields exclude (1x16), sep (1x6), psatz;
%        sp:     1xM 'cell' of 'cellstr', the variable names of each space;
%        xy:     1x2 'cellstr', the names of x and y.
% OUTPUT deg:    1xM 'cell', the degree specification of each space;
%        copts:  'struct' with fields psatz (POPTS.psatz, default 0), sep
%                (1x2 logical) and include (1xM 'cell' of multi-indices).
% The direction (1 x, 2 y) of each one-variable space; a variable other     % MMP, 10/28/2026
% than x and y is an error                                                  % MMP, 10/28/2026
dir1 = zeros(1,numel(sp));                                                  % MMP, 10/28/2026
for k = 1:numel(sp)                                                         % MMP, 10/28/2026
    if numel(sp{k})==1                                                      % MMP, 10/28/2026
        d = find(strcmp(xy,sp{k}{1}));                                      % MMP, 10/28/2026
        if isempty(d)                                                       % MMP, 10/28/2026
            error('poslpivar_settings_2d_sop:space','Space {%s}: not a variable of the PIE.',sp{k}{1}); % MMP, 10/28/2026
        end                                                                 % MMP, 10/28/2026
        dir1(k) = d;                                                        % MMP, 10/28/2026
    end                                                                     % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026
% One separation flag per direction serves every space in POSCOPVAR. In a   % MMP, 10/28/2026
% direction held by a one-variable space and by the L2[x,y] space, flags    % MMP, 10/28/2026
% that differ (sep(d) against sep(2+d) and sep(4+d)) leave the direction    % MMP, 10/28/2026
% not separated: more basis operators than the stock, never fewer           % MMP, 10/28/2026
s6 = zeros(1,6);                                                            % MMP, 10/28/2026
if isfield(popts,'sep') && ~isempty(popts.sep),  s6(1:numel(popts.sep)) = popts.sep(:)';  end % MMP, 10/28/2026
has2 = any(cellfun(@numel,sp)==2);                                          % MMP, 10/28/2026
sep1 = s6(1:2);                 % the flags of the one-variable spaces      % MMP, 10/28/2026
for d = 1:2                                                                 % MMP, 10/28/2026
    if any(dir1==d) && has2 && (s6(2+d) && s6(4+d))~=s6(d)                  % MMP, 10/28/2026
        warning('poslpivar_settings_2d_sop:lossy',['sep(%d) = %d against sep(%d), sep(%d) = %d, %d: '...
                'direction %d not separated.'],d,s6(d),2+d,4+d,s6(2+d),s6(4+d),d); % MMP, 10/28/2026
        popts.sep(2+d) = 0;     popts.sep(4+d) = 0;     sep1(d) = 0;        % MMP, 10/28/2026
    end                                                                     % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026
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
% bases of POS1D on L2[x] and L2[y], the translated basis operators on
% L2[x,y], an error otherwise
deg = cell(1,numel(sp));    incl = cell(1,numel(sp));
for k = 1:numel(sp)
    v = sort(sp{k});
    if isempty(v)
        if exc(1)
            error('poslpivar_settings_2d_sop:exclude1',['exclude(1) drops the R^n basis; '...
                  'poscopvar needs one basis operator per space.']);
        end
        deg{k} = struct('int',[0 0]);
    elseif numel(v)==1                                                      % MMP, 10/28/2026
        [deg{k},incl{k}] = pos1d(pdeg,exc,dir1(k),sep1(dir1(k)),numel(xy)); % MMP, 10/28/2026
    elseif numel(v)==2
        deg{k} = o.LF.deg;      incl{k} = o.LF.opts.include;
    else
%       error('poslpivar_settings_2d_sop:space',['Space {%s}: only R^n and L2[x,y] spaces '...
%             'are handled.'],strjoin(v,','));                              % MMP, 10/28/2026 (was)
        error('poslpivar_settings_2d_sop:space',['Space {%s}: only R^n, L2[x], L2[y] and '...
              'L2[x,y] spaces are handled.'],strjoin(v,','));               % MMP, 10/28/2026
    end
end
% psatz of POPTS unchanged, sep of the translation
ps = 0;     if isfield(popts,'psatz') && ~isempty(popts.psatz),  ps = popts.psatz;  end
copts = struct('psatz',ps,'sep',o.LF.opts.sep);
copts.include = incl;           % assigned apart: a cell value in struct() builds a struct array
% In a direction held by a one-variable space, the flag of that space,      % MMP, 10/28/2026
% equal to the translation's when an L2[x,y] space is present (see above)   % MMP, 10/28/2026
for d = unique(dir1(dir1>0))                                                % MMP, 10/28/2026
    copts.sep(d) = logical(sep1(d));                                        % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026
end


function [deg,incl] = pos1d(pdeg,exc,d,sepd,nv)                             % MMP, 10/28/2026
% [DEG,INCL] = POS1D(PDEG,EXC,D,SEPD,NV) returns the POSCOPVAR degrees and  % MMP, 10/28/2026
% included basis operators of an L2[x] (D = 1) or L2[y] (D = 2) space: the  % MMP, 10/28/2026
% bases Zx^o, Zx^a, Zx^b (Zy^o, Zy^a, Zy^b) of the degrees dx (dy) that     % MMP, 10/28/2026
% exclude(2:4) (exclude(5:7)) keeps, separated when SEPD is set, by the     % MMP, 10/28/2026
% 1-D translation of SETTINGS2POSSOPVAR. That translation reads the caps    % MMP, 10/28/2026
% of Zx^a in the order [middle variable; outer variable; joint], as         % MMP, 10/28/2026
% POSLPIVAR_2D applies dx{2} (Zxa(ss1,tt1), ss1 becoming the integrated     % MMP, 10/28/2026
% variable of P); its scalar int, mult become entry D of 1 x NV rows,       % MMP, 10/28/2026
% 0 in the direction the space lacks.                                       % MMP, 10/28/2026
% INPUT  pdeg:  'struct' with field dx (D = 1) or dy (D = 2), a 3-entry     % MMP, 10/28/2026
%               'cell' {d1; [d2(1) d2(2) d2(3)]; [d3(1) d3(2) d3(3)]};      % MMP, 10/28/2026
%        exc:   1x16 'double', the exclude entries of the term;             % MMP, 10/28/2026
%        d:     'double', 1 or 2;                                           % MMP, 10/28/2026
%        sepd:  'double' or 'logical', the separation flag of direction D;  % MMP, 10/28/2026
%        nv:    'double', the number of registry variables (2).             % MMP, 10/28/2026
% OUTPUT deg:   1xN 'cell' of 'struct' (fields int, mult: 1 x NV; joint),   % MMP, 10/28/2026
%               one per kept basis operator;                                % MMP, 10/28/2026
%        incl:  Nx1 'double', their multi-indices (1 multiplier, 2 lower,   % MMP, 10/28/2026
%               3 upper, 4 full integral).                                  % MMP, 10/28/2026
% The degrees of the direction; none is an error                            % MMP, 10/28/2026
f = {'dx','dy'};                                                            % MMP, 10/28/2026
if ~isfield(pdeg,f{d}) || isempty(pdeg.(f{d}))                              % MMP, 10/28/2026
    error('poslpivar_settings_2d_sop:deg1d','The settings give no %s degrees for the L2 space in that variable alone.',f{d}); % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026
c = pdeg.(f{d});                                                            % MMP, 10/28/2026
% The 1-D settings of the direction: exclude [R^n (not in this space),      % MMP, 10/28/2026
% multiplier, lower, upper]; no psatz term (the caller adds those)          % MMP, 10/28/2026
s = struct();                                                               % MMP, 10/28/2026
s.dd1 = {c{1}, reshape(c{2},1,[]), reshape(c{3},1,[])};                     % MMP, 10/28/2026
s.options1 = struct('exclude',[1, exc(3*d-1:3*d+1)],'sep',double(sepd));    % MMP, 10/28/2026
s.override1 = 1;                                                            % MMP, 10/28/2026
o = settings2possopvar(s);                                                  % MMP, 10/28/2026
lossy = o.report.lossy;     lossy = lossy(startsWith(lossy,'LF'));          % MMP, 10/28/2026
if ~isempty(lossy)                                                          % MMP, 10/28/2026
    warning('poslpivar_settings_2d_sop:lossy','%s',strjoin(lossy,' '));     % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026
% No basis left: POSCOPVAR needs one per space                              % MMP, 10/28/2026
if isempty(o.LF.deg)                                                        % MMP, 10/28/2026
    error('poslpivar_settings_2d_sop:exclude1d',['exclude(%d:%d) drops every basis of the '...
          'L2 space in %s alone; poscopvar needs one basis operator per space.'],3*d-1,3*d+1,f{d}(2)); % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026
% The scalar caps of each basis operator laid over the registry             % MMP, 10/28/2026
deg = cell(1,numel(o.LF.deg));                                              % MMP, 10/28/2026
for b = 1:numel(o.LF.deg)                                                   % MMP, 10/28/2026
    g = o.LF.deg{b};                                                        % MMP, 10/28/2026
    t = struct('int',zeros(1,nv),'mult',zeros(1,nv),'joint',g.joint);       % MMP, 10/28/2026
    t.int(d) = g.int;   t.mult(d) = g.mult;     deg{b} = t;                 % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026
incl = o.LF.opts.include;                                                   % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026


% function eq_opts = eqopts2d(Km,eq_opts,ztol)                              % MMP, 10/28/2026 (was)
function eq_opts = eqopts2d(Km,eq_opts,ztol,xy)                             % MMP, 10/28/2026
% EQ_OPTS = EQOPTS2D(KM,EQ_OPTS,ZTOL,XY) sets the exclude and sep entries of
% EQ_OPTS that the structure of KM allows, for R^n, L2[x], L2[y] and
% L2[x,y] spaces, as get_eq_opts_2D does for an 'opvar2d'. With (i,j) the
% parameter cells of an L2[x,y] -> L2[x,y] block (i, j = 1 multiplier, 2
% lower, 3 upper integral in x, in y), q the cells of an L2[x] -> L2[x]
% block (L2[y] -> L2[y]; q = 1 multiplier, 2 lower, 3 upper integral), and
% each condition holding in every such block of KM:
%   exclude(1)     every R^n -> R^n block zero;
%   exclude(2:4)   L2[x] cell q = 1, 2, 3 zero (exclude(5:7): L2[y]);
%   sep(1)         L2[x] cell 3 nonzero, exclude(3) unset, and cell 2 -
%                  cell 3 <= ZTOL (sep(2): L2[y], exclude(6));
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
%        ztol:    'double', the tolerance of the differences;
%        xy:      1x2 'cellstr', the names of x and y.
% OUTPUT eq_opts: EQ_OPTS with those entries set, the others unchanged.
exc = eq_opts.exclude;      sep = eq_opts.sep;
% Spaces: R^n (no variable), L2[x] and L2[y] (one), L2[x,y] (two); a
% registry of more than two variables is an error
nvs = sum(Km.space_out,2);
% if any(nvs==1) || numel(Km.vars)>2                                        % MMP, 10/28/2026 (was)
%     error('poslpivar_settings_2d_sop:space','Only R^n and L2[x,y] spaces are handled.'); % MMP, 10/28/2026 (was)
if numel(Km.vars)>2                                                         % MMP, 10/28/2026
    error('poslpivar_settings_2d_sop:space','Only spaces in the two variables of the PIE are handled.'); % MMP, 10/28/2026
end
iR = find(nvs==0);      iX = find(nvs==2);
% exclude(1): every R^n -> R^n block zero
if ~isempty(iR) && all_zero(Km,iR,1)
    exc(1) = 1;
end
% exclude(2:7), sep(1:2): the spaces in x alone, then in y alone, by the    % MMP, 10/28/2026
% rules of get_eq_opts_2D for Rxx (Ryy): a zero cell q excludes its basis;  % MMP, 10/28/2026
% with cell 3 nonzero and the lower-integral basis kept, cells 2 and 3      % MMP, 10/28/2026
% that agree separate the direction                                         % MMP, 10/28/2026
for d = 1:2                                                                 % MMP, 10/28/2026
    i1 = find(nvs==1 & any(Km.space_out(:,strcmp(Km.vars,xy{d})),2));       % MMP, 10/28/2026
    if isempty(i1),     continue,   end                                     % MMP, 10/28/2026
    e0 = 3*d - 1;       % exclude(2) in x, exclude(5) in y                  % MMP, 10/28/2026
    if all_zero(Km,i1,1),   exc(e0) = 1;    end                             % MMP, 10/28/2026
    if all_zero(Km,i1,2),   exc(e0+1) = 1;  end                             % MMP, 10/28/2026
    if all_zero(Km,i1,3)                                                    % MMP, 10/28/2026
        exc(e0+2) = 1;                                                      % MMP, 10/28/2026
    elseif ~exc(e0+1) && diff_le(Km,i1,2,3,ztol)                            % MMP, 10/28/2026
        sep(d) = 1;                                                         % MMP, 10/28/2026
    end                                                                     % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026
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
