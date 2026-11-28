function test_settings_prog
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_SETTINGS_PROG tests the settings of an LPI program (LPIPROGRAM,
% LPISETSETTINGS) and the domain that POSLPIVAR and LPIVAR read
% from the program, and errors at the first failed check:
%   (1) LPIPROGRAM stores the settings given as its last input, a tier name
%       or a structure, in each input form (vars and dom; vars, dummies and
%       dom; with free variables; with decision and free variables), and
%       'light' without them; the program is otherwise that of the same
%       form without settings;
%   (2) LPISETSETTINGS sets a tier, a structure, and 'light' for []; an
%       unknown tier name and a 'double' are refused;
%   (3) on a program whose variables are declared out of sorted order, s2
%       on [0,2] before s1 on [0,1], POSLPIVAR with DOM omitted builds the
%       operator of the call with DOM given, on R^2 x L2[s1] x L2[s1,s2]
%       (an ODE, 1-D and 2-D coupling), and not that of the swapped domain;
%   (4) LPIVAR with DOM omitted, and with DOM and DEG omitted,
%       builds the operator of the explicit call (DEG: the caps of 'light'
%       block by block, [2 2 2] for L2[s1], 1 for L2[s1,s2]);
%   (5) LPI_INEQ takes its degree defaults from prog.settings: on the
%       target Pt of TEST_LPI_INEQ the weight degree is 2 at 'light'
%       and 1 at 'stripped' (wmin 0); a settings structure passed to the
%       call wins over prog.settings, and so does an explicit dw;
%   (6) LPISOLVE without options passes prog.settings.sos_opts, and with a
%       settings structure its sos_opts (an unknown solver name is reported
%       by SOSSOLVE in both cases); with the tier options it solves;
%   (7) POSLPIVAR without DEG builds, on R^2 x L2[s1] with lpi.pos in the
%       stock-field form of 'light' and of 'heavy' on 11/12/2026, the
%       operator of the explicit specification of lpi.pos
%       (the multiplier, lower and upper bases), and on R^2 x L2[s1] x
%       L2[s1,s2] the multiplier and full-integral bases in every space
%       (s1 separable, being owned by the 2-D space); with DEG omitted an
%       OPTIONS.include is refused; LPIVAR without DEG on a 1-D
%       program at 'heavy' builds the operator of lpi.var.deg{1} = [4 2 2];
%   (8) EXEC_SOP_SETTINGS refuses a retired stock-format field (dd1) and
%       settings_nd.LF_int with the identifier exec_sop_settings:retired;
%       LPISETTINGS_STOCK('heavy') sets sop.slack = 'stock' and returns
%       sop.stock with the former dd2 of the heavy level, {3,[2 2 3],[2 2 3]};
%   (9) LPIVAR without DEG at 'heavy' (lpi.var.deg {[4 2 2],1,1})
%       on R x L2[s1] x L2[s1,s2]: the largest exponent of each variable
%       in the left basis of a block is the cap of the entry for its output
%       space, in the right basis that of its input space (L2[s1]: 4 left,
%       2 right in s1 shared, 4 right in s1 input-only; L2[s1,s2]: 1); on
%       R x L2[s1,s2] the operator of entry 2;
%  (10) POSLPIVAR without DEG reads the storage row lpi.pos (len_com) of
%       every tier: on R^2 x L2[s1], L2[s1,s2], L2[s1,s2,s3] and R^2 x
%       L2[s1] x L2[s1,s2], the basis operators (space, multi-index) and
%       Gram sizes derived by hand from the tier files and the header of
%       POSLPIVAR (6 tiers x 4 configurations);
%  (11) on a row with lengths distinct per letter, direction and entry,
%       the operator of each configuration equals that of the explicit
%       DEG of the header, and a b with its theta- and s-lengths exchanged
%       does not; the Gram sizes with len_cross, owned and non-owned caps,
%       len_R 2 and absent, letter words (equal to code rows) and a zero
%       length that drops a letter; with both lpi.pos layouts the row is
%       read; DOM given with DEG omitted or empty (the call of POS_ROW)
%       gives the operator of both omitted; five inputs refused by message;
%  (12) PIETOOLS_stability at 'light' certifies x_t = x_ss + pi^2/2 x,
%       x(0) = x(1) = 0, and its storage, the first entry of prog.sop_pos,
%       is the operator of the hand specification of entry 1 on L2[s];
%  (13) a space in 4 variables reads entry 3, the last entry: 'uniform',
%       a of length 1 and b of 2 x 2 per direction, from the scalar entry
%       3 of (11) and from 'light', whose entry 3 holds rows of 3 lengths.
%
% INPUT:  none.
% OUTPUT: none; one line per passed part, then a summary line.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_settings_prog
%
% Copyright (C)2026  PIETOOLS Team
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
% Initial coding MMP - 11/12/2026: settings on the program (maintainer,
%   11/12/2026: "Settings on the program"). Part (3) guards the domain by
%   name: the varname of a polynomial vector is sorted, so a domain read
%   from it would follow the sorted names, not the rows of prog.dom.
% MMP, 11/12/2026: Parts (7) and (8): the degrees of POSCOPVAR and
%   LPIVAR_CDOPVAR from the program settings when DEG is omitted, the
%   retired stock-format fields and LPISETTINGS_STOCK_SOP (phase 2 of the
%   knob reduction).
% MMP, 11/12/2026: Part (4) compares DEG omitted with the caps of each block
%   from its two spaces, part (9) reads those caps off the bases (phase 3:
%   LPIVAR_CDOPVAR defaults per block, no longer the largest space).
% MMP, 11/23/2026: Part (7) runs on lpi.pos in the stock-field form (mult,
%   int, joint, sep; the 'light' and 'heavy' entries of 11/12/2026), which
%   POSCOPVAR reads without len_com; the tiers hold the storage row since
%   11/18/2026, and part (7) errored at ps.int. No part checks the
%   storage-row branch of POSCOPVAR. Help to match. Part (9) runs on the
%   'heavy' lpi.var of 11/12/2026, {[4 2 2],1,1}, which its expected
%   exponents assume (the tier has entry 2 = 2 since 11/18/2026); its
%   failure message passed a 1 x 2 block index to one %d.
% MMP, 11/24/2026: Parts (10)-(13) check the storage-row branch of POSCOPVAR
%   (lpi.pos with len_com: every tier since 11/18/2026, and every executive
%   storage through POS_ROW), which no part checked: the tier tables, each
%   rule on a synthetic row against the explicit DEG, the storage of an
%   executive, a space in 4 variables. Local GRAM_SIZES; help to match.

warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
npass = 0;

%% (1) LPIPROGRAM: the settings as the last input, in every input form
vars = polynomial({'s2';'s1'});     dum = polynomial({'s2_dum';'s1_dum'});
dom = [0 2; 0 1];                   free = polynomial({'t'});
light = lpisettings();          heavy = lpisettings('heavy');
p0 = lpiprogram(vars,dom);
pairs = {p0, lpiprogram(vars,dom,'heavy'); ...
         lpiprogram(vars,dum,dom), lpiprogram(vars,dum,dom,heavy); ...
         lpiprogram(vars,dum,dom,free), lpiprogram(vars,dum,dom,free,"heavy"); ...
         lpiprogram(vars,dum,dom,dpvar(zeros(0,1)),free), ...
         lpiprogram(vars,dum,dom,dpvar(zeros(0,1)),free,'heavy')};
for k = 1:size(pairs,1)
    a = pairs{k,1};     b = pairs{k,2};
    % isequaln: the tiers hold NaN entries (sop.cap_delta)
    if ~isequaln(a.settings,light) || ~isequaln(b.settings,heavy)
        error('test_settings_prog: (1) form %d stores the wrong settings.',k);
    end
    if ~isequaln(rmfield(a,'settings'),rmfield(b,'settings'))
        error('test_settings_prog: (1) form %d: the settings changed the program.',k);
    end
end
fprintf('  passed: (1) settings stored in %d input forms, programs otherwise equal\n',size(pairs,1));
npass = npass+1;

%% (2) LPISETSETTINGS
custom = lpisettings('custom');
ok = isequaln(lpisetsettings(p0,'veryheavy').settings,lpisettings('veryheavy')) ...
     && isequaln(lpisetsettings(p0,custom).settings,custom) ...
     && isequaln(lpisetsettings(p0,[]).settings,light);
if ~ok
    error('test_settings_prog: (2) LPISETSETTINGS stores the wrong settings.');
end
refused = 0;
try,    lpisetsettings(p0,'nosuchtier');    catch,  refused = refused+1;    end
try,    lpisetsettings(p0,5);               catch,  refused = refused+1;    end
if refused~=2
    error('test_settings_prog: (2) an unknown tier or a double was accepted.');
end
fprintf('  passed: (2) tier, structure and [] set; unknown tier and double refused\n');
npass = npass+1;

%% (3) POSLPIVAR: DOM from the program, by name
sp = {{},{'s1'},{'s1','s2'}};       dm = [2;1;1];
[pA,PA] = poslpivar(p0,dm,sp,[0 1; 0 2],1);     % rows in sorted order s1, s2
[pB,PB] = poslpivar(p0,dm,sp,1);
[~,PC] = poslpivar(p0,dm,sp,[0 2; 0 1],1);      % the swapped domain
if ~isequal(PA,PB) || ~isequal(pA.decvartable,pB.decvartable)
    error('test_settings_prog: (3) DOM omitted builds another operator.');
end
if isequal(PB,PC)
    error('test_settings_prog: (3) the swapped domain gives the same operator; the check is void.');
end
fprintf('  passed: (3) poslpivar on R^2 x L2[s1] x L2[s1,s2]: DOM omitted = given (%d names), swapped differs\n',...
        numel(pB.decvartable));
npass = npass+1;

%% (4) LPIVAR: DOM, and DOM and DEG, from the program
dz = struct('out',[1;1],'in',[1;1]);    sz = {{'s1'},{'s1','s2'}};
[~,ZA] = lpivar(p0,dz,sz,[0 1; 0 2],1);
[~,ZB] = lpivar(p0,dz,sz,1);
[~,ZC] = lpivar(p0,dz,sz);
% DEG omitted at light: entry 1 ([2 2 2]) for L2[s1], entry 2 (1) for       % MMP, 11/12/2026
% L2[s1,s2]; left caps from the output space, right from the input space    % MMP, 11/12/2026
dgc = {[2 2 2], struct('mult',2,'int',[2 1],'out',2,'in',1); ...
       struct('mult',1,'int',[1 2],'out',1,'in',2), 1};                     % MMP, 11/12/2026
[~,ZD] = lpivar(p0,dz,sz,[0 1; 0 2],dgc);                                   % MMP, 11/12/2026
% if ~isequal(ZA,ZB) || ~isequal(ZA,ZC)                                     % MMP, 11/12/2026 (was)
if ~isequal(ZA,ZB) || ~isequal(ZD,ZC)                                       % MMP, 11/12/2026
    error('test_settings_prog: (4) DOM or DEG omitted builds another operator.');
end
fprintf('  passed: (4) lpivar: DOM omitted, and DOM and DEG omitted, = explicit\n');
npass = npass+1;

%% (5) LPI_INEQ: degree defaults from prog.settings
% Pt of TEST_LPI_INEQ on L2[0,1]: lower kernel 2 - s, upper 2 - t
v1 = {'s1'};    d1 = [0 1];     ZL = {[0;1]};   ZR = {[0;1]};
params = reshape({sparse(2,2),sparse([2 0;-1 0]),sparse([2 -1;0 0])},[3 1 1]);
Pt = opvar({sopvar(params,struct('out',{v1},'in',{v1}),ZL,ZR,struct('out',d1,'in',d1),[1 1])});
pL = lpiprogram(polynomial({'s1'}),[],d1);
pS = lpiprogram(polynomial({'s1'}),[],d1,'stripped');
[~,~,iL] = lpi_ineq(pL,Pt);
[~,~,iS] = lpi_ineq(pS,Pt);
[~,~,iX] = lpi_ineq(pS,Pt,heavy);                   % the call wins
[~,~,iD] = lpi_ineq(pL,Pt,struct('dw',1));          % explicit dw: no floor
w = [iL.degrees.w, iS.degrees.w, iX.degrees.w, iD.degrees.w];
if ~isequal(w,[2 1 2 1])
    error('test_settings_prog: (5) weights %s, expected [2 1 2 1] (light, stripped, call heavy, dw 1).',mat2str(w));
end
fprintf('  passed: (5) lpi_ineq weight: light %d, stripped %d, call heavy %d, dw 1 %d\n',w);
npass = npass+1;

%% (6) LPISOLVE: the solver options of the program, or of a settings structure
[pE,g] = lpidecvar(pL,'g');
pE = lpi_ineq(pE,g);    pE = lpisetobj(pE,g);
bad = pE;   bad.settings.sos_opts.solver = 'nosuchsolver';
stb = light;    stb.sos_opts.solver = 'nosuchsolver';
seen = 0;
try,    lpisolve(bad);      catch ME,   seen = seen + contains(ME.message,'nosuchsolver');    end
try,    lpisolve(pE,stb);   catch ME,   seen = seen + contains(ME.message,'nosuchsolver');    end
if seen~=2
    error('test_settings_prog: (6) LPISOLVE did not pass the solver options of the settings.');
end
evalc('ps = lpisolve(pE);');
gv = double(lpigetsol(ps,g));
if ~(abs(gv) < 1e-6)
    error('test_settings_prog: (6) min g s.t. g >= 0 returned g = %g.',gv);
end
fprintf('  passed: (6) lpisolve: prog.settings.sos_opts and a settings structure passed; min g = %.1e\n',gv);
npass = npass+1;

%% (7) POSLPIVAR and LPIVAR without DEG: the program settings
% lpi.pos in the stock-field form of 11/12/2026, which POSLPIVAR reads      % MMP, 11/23/2026
% without len_com (the tiers hold the storage row since 11/18/2026)         % MMP, 11/23/2026
pf = @(m) struct('mult',[m 0 0],'int',[1 1 1],'joint',[2 2 2],'sep',0);     % MMP, 11/23/2026
o1 = light;     o1.lpi.pos = pf(1);         % 'light' on 11/12/2026         % MMP, 11/23/2026
o2 = heavy;     o2.lpi.pos = pf(2);         % 'heavy' on 11/12/2026         % MMP, 11/23/2026
% q1 = lpiprogram(polynomial({'s1'}),[0 1]);          % light               % MMP, 11/23/2026 (was)
% q2 = lpiprogram(polynomial({'s1'}),[0 1],'heavy');                        % MMP, 11/23/2026 (was)
q1 = lpiprogram(polynomial({'s1'}),[0 1],o1);       % light                 % MMP, 11/23/2026
q2 = lpiprogram(polynomial({'s1'}),[0 1],o2);       % heavy                 % MMP, 11/23/2026
sp1 = {{},{'s1'}};      dm1 = [2;1];
for q = {q1, q2}
    ps = q{1}.settings.lpi.pos;
    Zi = struct('int',ps.int(1),'mult',ps.int(1),'joint',ps.joint(1));
    dx = {struct('int',0), {struct('int',ps.mult(1),'mult',0,'joint',ps.mult(1)),Zi,Zi}};
    ox = struct('sep',false,'include',{{[],[1;2;3]}});
    [qa,Pa] = poslpivar(q{1},dm1,sp1);
    [qb,Pb] = poslpivar(q{1},dm1,sp1,dx,ox);
    if ~isequaln(Pa,Pb) || ~isequaln(qa.decvartable,qb.decvartable)
        error('test_settings_prog: (7) poslpivar without DEG differs from the explicit lpi.pos (mult %d).',ps.mult(1));
    end
end
% the coupled spaces: s1 separable for the 1-D space too
% st2 = p0.settings.lpi.pos;                                                % MMP, 11/23/2026 (was)
p0o = lpiprogram(vars,dom,o1);              % p0 with the stock fields      % MMP, 11/23/2026
st2 = p0o.settings.lpi.pos;                                                 % MMP, 11/23/2026
Zi1 = struct('int',st2.int(1),'mult',st2.int(1),'joint',st2.joint(1));
Zi2 = struct('int',st2.int(2),'mult',st2.int(2),'joint',st2.joint(2));
dx = {struct('int',0), {struct('int',st2.mult(1),'mult',0,'joint',st2.mult(1)),Zi1}, ...
      {struct('int',st2.mult(2),'mult',0,'joint',st2.mult(2)),Zi2}};
ox = struct('sep',true(1,2));   ox.include = {[],[1;4],[1 1;4 4]};
% [~,Pc] = poscopvar(p0,dm,sp);                                             % MMP, 11/23/2026 (was)
% [~,Pd] = poscopvar(p0,dm,sp,dx,ox);                                       % MMP, 11/23/2026 (was)
[~,Pc] = poslpivar(p0o,dm,sp);                                              % MMP, 11/23/2026
[~,Pd] = poslpivar(p0o,dm,sp,dx,ox);                                        % MMP, 11/23/2026
if ~isequaln(Pc,Pd)
    error('test_settings_prog: (7) the coupled default differs from multiplier and full integral in every space.');
end
refused = 0;
try,    poslpivar(q1,dm1,sp1,[],struct('include',{{[],[1;2]}}));   catch,  refused = 1;    end
if ~refused
    error('test_settings_prog: (7) OPTIONS.include with DEG omitted was accepted.');
end
dz1 = struct('out',[1;1],'in',[1;1]);   sz1 = {{},{'s1'}};
[~,Za] = lpivar(q2,dz1,sz1);
[~,Zb] = lpivar(q2,dz1,sz1,[4 2 2]);
if ~isequaln(Za,Zb)
    error('test_settings_prog: (7) lpivar without DEG differs from lpi.var.deg{1}.');
end
fprintf('  passed: (7) poslpivar and lpivar without DEG: the settings at light and heavy, coupled spaces separable\n');
npass = npass+1;

%% (8) Retired fields refused; the stock-slack utility
bad = {setfield(lpisettings(),'dd1',{1,[1 1 2],[1 1 2]}), ...
       setfield(lpisettings(),'settings_nd',struct('LF_int',1))};
nref = 0;
for k = 1:numel(bad)
    try,    exec_sop_settings(bad{k});
    catch ME,   nref = nref + strcmp(ME.identifier,'exec_sop_settings:retired');
    end
end
[sk,stock] = lpisettings_stock('heavy');
if nref~=2 || ~strcmp(sk.sop.slack,'stock') || ~isequal(stock,sk.sop.stock) || ...
        ~isequal(stock.dd2,{3,[2 2 3],[2 2 3]})
    error('test_settings_prog: (8) retired fields accepted (%d of 2 refused) or the stock slack of heavy wrong.',nref);
end
fprintf('  passed: (8) retired fields refused (2 of 2); the stock slack of heavy from LPISETTINGS_STOCK\n');
npass = npass+1;

%% (9) LPIVAR without DEG: the caps of each block from its two spaces
% heavy: entry 1 [4 2 2] (mult = out = in = 4, int = [2 2]), entry 2 every
% cap 1. Rows: block, side, its variables, their largest exponents.
% h = lpiprogram(polynomial({'s1';'s2'}),[0 1; 0 2],'heavy');               % MMP, 11/23/2026 (was)
hv = heavy;     hv.lpi.var.deg = {[4 2 2],1,1};     % 'heavy' on 11/12/2026 % MMP, 11/23/2026
h = lpiprogram(polynomial({'s1';'s2'}),[0 1; 0 2],hv);                      % MMP, 11/23/2026
sp3 = {{},{'s1'},{'s1','s2'}};      d3 = struct('out',[1;1;1],'in',[1;1;1]);
[~,Zh] = lpivar(h,d3,sp3);
want = {[2 2],'ZL',{'s1'},4;            [2 2],'ZR',{'s1'},2;
        [3 3],'ZL',{'s1','s2'},[1 1];   [3 3],'ZR',{'s1','s2'},[1 1];
        [2 3],'ZL',{'s1'},4;            [2 3],'ZR',{'s1','s2'},[1 1];
        [3 2],'ZL',{'s1','s2'},[1 1];   [3 2],'ZR',{'s1'},2;
        [1 2],'ZR',{'s1'},4;            [1 3],'ZR',{'s1','s2'},[1 1];
        [2 1],'ZL',{'s1'},4;            [3 1],'ZL',{'s1','s2'},[1 1]};
for k = 1:size(want,1)
    B = Zh.C{want{k,1}(1),want{k,1}(2)};
    if strcmp(want{k,2},'ZL'),  v = B.vars.out;     Zb = B.ZL;
    else,                       v = B.vars.in;      Zb = B.ZR;
    end
    [in,loc] = ismember(want{k,3},v);
    if ~all(in) || numel(v)~=numel(want{k,3}) || ~isequal(cellfun(@max,Zb(loc)),want{k,4})
%             want{k,1},want{k,2},strjoin(v,','),mat2str(cellfun(@max,Zb)),mat2str(want{k,4})); % MMP, 11/23/2026 (was)
        error('test_settings_prog: (9) block (%d,%d) %s over %s: largest exponents %s, expected %s.', ...
              want{k,1}(1),want{k,1}(2),want{k,2},strjoin(v,','),mat2str(cellfun(@max,Zb)),mat2str(want{k,4})); % MMP, 11/23/2026
    end
end
% One dimension besides R^n: every block from that entry, as before
dz2 = struct('out',[1;1],'in',[1;1]);   sz2 = {{},{'s1','s2'}};
[~,Ze] = lpivar(h,dz2,sz2);
[~,Zf] = lpivar(h,dz2,sz2,1);
if ~isequal(Ze,Zf)
    error('test_settings_prog: (9) on R x L2[s1,s2] the default differs from entry 2.');
end
fprintf('  passed: (9) lpivar without DEG: %d block sides capped by their own space; R x L2[s1,s2] = entry 2\n',size(want,1));
npass = npass+1;

% ---- MMP, 11/24/2026: begin of parts (10)-(13), the storage row ----
%% (10) POSLPIVAR without DEG: the storage row of every tier
% Per tier and configuration, the basis operators and Gram sizes that the
% header of POSLPIVAR (DEG OMITTED) gives for the row of the tier file,
% derived by hand. Rows: space; alpha_d per registry direction (1
% multiplier, 2 lower, 3 upper, 4 full integral, 0 not owned); Gram size
% m_k*|Z_b|, |Z_b| the product over owned directions of theta-length times
% s-length (1 for a), times len_cross = 1 in a direction not owned. In the
% coupled configuration the 1-D space carries c in s1, so b is the upper
% integral in s1 in both spaces, the full integral in s2.
cf = {polynomial({'s1'}),           [0 1],          [2;1],   {{},{'s1'}};
      polynomial({'s1';'s2'}),      [0 1;0 2],      1,       {{'s1','s2'}};
      polynomial({'s1';'s2';'s3'}), [0 1;0 2;0 1],  1,       {{'s1','s2','s3'}};
      polynomial({'s1';'s2'}),      [0 1;0 2],      [2;1;1], {{},{'s1'},{'s1','s2'}}};
% stripped: entry 1 'all' (a 2; b, c 2x2); 2, 3 'uniform' (b 2x2, 1x1)
Est = {[1 0 2; 2 1 2; 2 3 4; 2 2 4],  [1 1 1 1; 1 4 4 16],  [1 1 1 1 1; 1 4 4 4 1], ...
       [1 0 0 2; 2 1 0 2; 2 3 0 4; 2 2 0 4; 3 1 1 1; 3 3 4 16]};
% light, custom: entry 2 'separable' (a 1; b 2x2), entry 3 'uniform'
Eli = {[1 0 2; 2 1 2; 2 3 4; 2 2 4],  [1 1 1 1; 1 4 1 4; 1 1 4 4; 1 4 4 16], ...
       [1 1 1 1 1; 1 4 4 4 64], ...
       [1 0 0 2; 2 1 0 2; 2 3 0 4; 2 2 0 4; 3 1 1 1; 3 3 1 4; 3 1 4 4; 3 3 4 16]};
% heavy: entry 1 a 3; entry 2 a 2, b 2x2
Ehe = {[1 0 2; 2 1 3; 2 3 4; 2 2 4],  [1 1 1 4; 1 4 1 8; 1 1 4 8; 1 4 4 16], ...
       [1 1 1 1 1; 1 4 4 4 64], ...
       [1 0 0 2; 2 1 0 3; 2 3 0 4; 2 2 0 4; 3 1 1 4; 3 3 1 8; 3 1 4 8; 3 3 4 16]};
% veryheavy: entry 1 a 4, b and c 4x4; entry 2 a 2, b 3x3
Evh = {[1 0 2; 2 1 4; 2 3 16; 2 2 16],  [1 1 1 4; 1 4 1 18; 1 1 4 18; 1 4 4 81], ...
       [1 1 1 1 1; 1 4 4 4 64], ...
       [1 0 0 2; 2 1 0 4; 2 3 0 16; 2 2 0 16; 3 1 1 4; 3 3 1 18; 3 1 4 18; 3 3 4 81]};
% extreme: the multiplier alone, a of length 2, 1, 1
Eex = {[1 0 2; 2 1 2],  [1 1 1 1],  [1 1 1 1 1],  [1 0 0 2; 2 1 0 2; 3 1 1 1]};
tiers = {'stripped',Est; 'light',Eli; 'custom',Eli; 'heavy',Ehe; 'veryheavy',Evh; 'extreme',Eex};
for t = 1:size(tiers,1)
    for c = 1:size(cf,1)
        [~,~,Qc,bl] = poslpivar(lpiprogram(cf{c,1},cf{c,2},tiers{t,1}),cf{c,3},cf{c,4});
        G = [bl,gram_sizes(Qc)];
        if ~isequal(G,tiers{t,2}{c})
            error('test_settings_prog: (10) %s, configuration %d: rows %s, expected %s.', ...
                  tiers{t,1},c,mat2str(G),mat2str(tiers{t,2}{c}));
        end
    end
end
fprintf('  passed: (10) poslpivar without DEG: the storage rows of %d tiers, %d configurations\n', ...
        size(tiers,1),size(cf,1));
npass = npass+1;

%% (11) The rules of the row, on lengths distinct per letter, direction, entry
% The explicit DEG of the header, by hand, is the oracle: equal operator,
% decision variables, Gram names and basis rows; Gram sizes as in (10).
% No len_cross: in a direction not owned, the largest owned int.
S = struct('len_R',1,'len_inc',{{3,[2 0],1}},'len_com',{{2,[3 2],2}}, ...
           'len_sos',{{[2;3;4],[1 2;2 3;3 1],[1;2;0]}},'patterns',{{'all','all','uniform'}}, ...
           'cap',{{Inf,[Inf Inf],Inf}},'terms',{{'none','none','none'}},'offsets',{{[],[],[]}});
sS = light;     sS.lpi.pos = S;
mk = @(i,m) struct('int',i,'mult',m);
% R x L2[s1]: a; b upper (s1 carries c), int 2, mult 1; c lower, int 3, mult 2
X1 = {{struct('int',0),{mk(1,0),mk(2,1),mk(3,2)}}, struct('sep',false,'include',{{[],[1;3;2]}})};
% L2[s1,s2]: 'all' less c in s2 (len_inc 0 there): b upper in s1, full in s2
D2 = {mk([0 1],[0 0]),mk([1 1],[2 0]),mk([2 1],[1 0]),mk([0 2],[0 1]),mk([1 2],[2 1]),mk([2 2],[1 1])};
I2 = [1 1; 3 1; 2 1; 1 4; 3 4; 2 4];
X2 = {{D2}, struct('sep',[false true],'include',{{I2}})};
% L2[s1,s2,s3]: 'uniform', the column of len_sos in every direction
X3 = {{{mk([0 0 0],[0 0 0]),mk([1 1 1],[1 1 1])}}, struct('sep',true(1,3),'include',{{[1 1 1; 4 4 4]}})};
% R x L2[s1] x L2[s1,s2]: the 1-D space has in s2 the int of s1
X4 = {{struct('int',0),{mk([1 1],[0 0]),mk([2 2],[1 0]),mk([3 3],[2 0])},D2}, ...
      struct('sep',[false true],'include',{{[],[1;3;2],I2}})};
cs = {[1;1],   X1, [1 0 1; 2 1 2; 2 3 6; 2 2 12];
      1,       X2, [1 1 1 2; 1 3 1 12; 1 2 1 12; 1 1 4 6; 1 3 4 36; 1 2 4 36];
      1,       X3, [1 1 1 1 1; 1 4 4 4 64];
      [1;1;1], X4, [1 0 0 1; 2 1 0 4; 2 3 0 18; 2 2 0 48; 3 1 1 2; 3 3 1 12; 3 2 1 12; ...
                    3 1 4 6; 3 3 4 36; 3 2 4 36]};
pr0 = cell(1,4);    P0 = cell(1,4);
for c = 1:4
    pr0{c} = lpiprogram(cf{c,1},cf{c,2},sS);
    [qa,P0{c},Qa,ba] = poslpivar(pr0{c},cs{c,1},cf{c,4});
    [qb,Pb,Qb,bb] = poslpivar(pr0{c},cs{c,1},cf{c,4},cs{c,2}{:});
    if ~isequaln(P0{c},Pb) || ~isequaln(qa.decvartable,qb.decvartable) || ~isequal(Qa,Qb) || ~isequal(ba,bb)
        error('test_settings_prog: (11) configuration %d: the row differs from the explicit DEG.',c);
    end
    if ~isequal([ba,gram_sizes(Qa)],cs{c,3})
        error('test_settings_prog: (11) configuration %d: rows %s, expected %s.',c, ...
              mat2str([ba,gram_sizes(Qa)]),mat2str(cs{c,3}));
    end
end
% Control: b of R x L2[s1] with its theta- and s-lengths exchanged keeps the
% Gram size and changes the operator, so the comparison separates int, mult
Xs = X1;    Xs{1}{2}{2} = mk(1,2);
[~,Ps] = poslpivar(pr0{1},[1;1],cf{1,4},Xs{:});
if isequaln(P0{1},Ps)
    error('test_settings_prog: (11) b with theta- and s-lengths exchanged builds the same operator.');
end
% Variants, by Gram size: len_cross 1 (the 1-D space has no theta-monomial
% in s2); theta_1+s_1 <= 2 (b 5 of 3x2, c 6 of 4x3); cap{2} = [1 Inf],
% read by the 2-D space in s1 and by the 1-D space in s2 (theta_2 <= 1: 2
% of 2, 3, 4); len_R 2 (degree 1 in theta_1 on R^2) and absent (1); the
% words {'ab','ba','cb'} ('cb' puts c in s1: b upper there); a zero
% theta-length of b in s2 drops b there
vr = {@(p) setfield(p,'len_cross',{1,1,1}), 4, [1;1;1], ...
         [1 0 0 1; 2 1 0 2; 2 3 0 6; 2 2 0 12; 3 1 1 2; 3 3 1 12; 3 2 1 12; 3 1 4 6; 3 3 4 36; 3 2 4 36];
      @(p) setfield(p,'cap',{2,[Inf Inf],Inf}), 1, [1;1], [1 0 1; 2 1 2; 2 3 5; 2 2 6];
      @(p) setfield(p,'cap',{Inf,[1 Inf],Inf}), 4, [1;1;1], ...
         [1 0 0 1; 2 1 0 4; 2 3 0 12; 2 2 0 24; 3 1 1 2; 3 3 1 6; 3 2 1 6; 3 1 4 6; 3 3 4 18; 3 2 4 18];
      @(p) setfield(p,'len_R',2), 1, [2;1], [1 0 4; 2 1 2; 2 3 6; 2 2 12];
      @(p) rmfield(p,'len_R'), 1, [2;1], [1 0 2; 2 1 2; 2 3 6; 2 2 12];
      @(p) setfield(p,'patterns',{'all',{'ab','ba','cb'},'uniform'}), 2, 1, [1 1 4 6; 1 3 1 12; 1 2 4 36];
      @(p) setfield(p,'len_sos',{[2;3;4],[1 2;2 0;3 1],[1;2;0]}), 2, 1, [1 1 1 2; 1 3 1 12; 1 2 1 12]};
Pv = cell(size(vr,1),1);
for r = 1:size(vr,1)
    v = sS;     v.lpi.pos = vr{r,1}(S);     c = vr{r,2};
    [~,Pv{r},Qv,bv] = poslpivar(lpiprogram(cf{c,1},cf{c,2},v),vr{r,3},cf{c,4});
    if ~isequal([bv,gram_sizes(Qv)],vr{r,4})
        error('test_settings_prog: (11) variant %d: rows %s, expected %s.',r, ...
              mat2str([bv,gram_sizes(Qv)]),mat2str(vr{r,4}));
    end
end
% The code rows [1 2;2 1;3 2] spell the patterns of the words above
v = sS;     v.lpi.pos.patterns = {'all',[1 2;2 1;3 2],'uniform'};
[~,Pr] = poslpivar(lpiprogram(cf{2,1},cf{2,2},v),1,cf{2,4});
% Both layouts in lpi.pos: len_com decides, the stock fields are not read
v = sS;     v.lpi.pos.mult = [3 0 0];   v.lpi.pos.int = [2 2 2];
v.lpi.pos.joint = [4 4 4];  v.lpi.pos.sep = 1;
[~,Pt] = poslpivar(lpiprogram(cf{1,1},cf{1,2},v),[1;1],cf{1,4});
% DOM given, DEG omitted (the call of POS_ROW in EXEC_TOOLS) or empty
dm4 = struct('vars',{{'s1','s2'}},'dom',[0 1;0 2]);
[~,Pf1] = poslpivar(pr0{4},[1;1;1],cf{4,4},dm4);
[~,Pf2] = poslpivar(pr0{4},[1;1;1],cf{4,4},dm4,[]);
eq = [isequaln(Pv{6},Pr),isequaln(Pt,P0{1}),isequaln(Pf1,P0{4}),isequaln(Pf2,P0{4})];
if ~all(eq)
    error('test_settings_prog: (11) words = rows, both layouts, DOM given: %s (1 = equal).',mat2str(eq));
end
% Refused, by message: a term in the storage, an unknown pattern set, a
% length row of 3 entries for a space in 2 variables; include or sep in
% OPTIONS with DEG omitted
bad = {@(p) setfield(p,'terms',{'psatz','none','none'}),      'carries no term';
       @(p) setfield(p,'patterns',{'bogus','all','uniform'}), 'Unknown pattern set';
       @(p) setfield(p,'len_com',{2,[3 2 1],2}),              'entries for a space in 2 variables'};
nref = 0;
for r = 1:size(bad,1)
    v = sS;     v.lpi.pos = bad{r,1}(S);
    try,    poslpivar(lpiprogram(cf{4,1},cf{4,2},v),[1;1;1],cf{4,4});
    catch ME,   nref = nref + contains(ME.message,bad{r,2});
    end
end
for o = {struct('include',{{[],[1;2],[1 1]}}), struct('sep',true)}
    try,    poslpivar(pr0{4},[1;1;1],cf{4,4},dm4,[],o{1});
    catch ME,   nref = nref + contains(ME.message,'OPTIONS may not set them');
    end
end
if nref~=5
    error('test_settings_prog: (11) %d of 5 inputs refused with their message.',nref);
end
fprintf('  passed: (11) the row rules in 4 configurations, %d variants, 5 refusals\n',size(vr,1));
npass = npass+1;

%% (12) The storage of an executive: PIETOOLS_stability at 'light'
% x_t = x_ss + pi^2/2 x, x(0) = x(1) = 0 (decay rate pi^2/2): the storage
% that POS_ROW of EXEC_TOOLS declares from the row, the first entry of
% prog.sop_pos, is the hand specification of entry 1 on L2[s]
pvar s t
x = pde_var(1,s,[0,1]);
PDE = [diff(x,t,1)==diff(x,s,2)+0.5*pi^2*x;  subs(x,s,0)==0;  subs(x,s,1)==0]; %#ok<NASGU> read in evalc
PIE = [];   pe = [];    ie = [];    % assigned in evalc
evalc('PIE = initialize(convert(PDE));');
sl = light;     sl.sop.verbose = false;
evalc('[pe,~,ie] = PIETOOLS_stability(PIE,sl,''P'');');
pf = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom,sl);
[~,Ph] = poslpivar(pf,1,{{'s'}},{{mk(1,0),mk(1,1),mk(1,1)}},struct('sep',false,'include',{{[1;3;2]}}));
if ~ie.certified || ~isequaln(pe.sop_pos{1},Ph)
    error('test_settings_prog: (12) stability at light: certified %d, storage of the row %d.', ...
          ie.certified,isequaln(pe.sop_pos{1},Ph));
end
fprintf('  passed: (12) PIETOOLS_stability at light: certified, its storage the row of entry 1\n');
npass = npass+1;

%% (13) A space in 4 variables reads entry 3, the last entry
% 'uniform': a of length 1, b of 2 x 2 per direction; from the scalar entry
% 3 of (11), and from 'light', whose entry 3 holds rows of 3 lengths
v4 = polynomial({'s1';'s2';'s3';'s4'});     d4 = repmat([0 1],4,1);
s4 = {{'s1','s2','s3','s4'}};   E4 = [1 1 1 1 1 1; 1 4 4 4 4 256];
for o = {'the row of (11)','light'; sS,light}
    try
        [~,~,Q4,b4] = poslpivar(lpiprogram(v4,d4,o{2}),1,s4);
    catch ME
        error('test_settings_prog: (13) %s on L2[s1,s2,s3,s4]: %s',o{1},ME.message);
    end
    if ~isequal([b4,gram_sizes(Q4)],E4)
        error('test_settings_prog: (13) %s on L2[s1,s2,s3,s4]: rows %s, expected %s.',o{1}, ...
              mat2str([b4,gram_sizes(Q4)]),mat2str(E4));
    end
end
fprintf('  passed: (13) a space in 4 variables reads entry 3 (the row of (11), light)\n');
npass = npass+1;
% ---- MMP, 11/24/2026: end of parts (10)-(13) ----

% The summary line
fprintf('\ntest_settings_prog_sop passed (%d checks).\n',npass);
end


function g = gram_sizes(Qc)                                                 % MMP, 11/24/2026
% G = GRAM_SIZES(QC) the size of each diagonal block of the Gram matrix.    % MMP, 11/24/2026
% INPUT  Qc: N_b x N_b cell, the Gram names of POSLPIVAR.                   % MMP, 11/24/2026
% OUTPUT g:  N_b x 1 double, the number of rows of block (b,b).             % MMP, 11/24/2026
g = zeros(size(Qc,1),1);                                                    % MMP, 11/24/2026
for b = 1:numel(g),     g(b) = size(Qc{b,b},1);     end                     % MMP, 11/24/2026
end                                                                         % MMP, 11/24/2026
