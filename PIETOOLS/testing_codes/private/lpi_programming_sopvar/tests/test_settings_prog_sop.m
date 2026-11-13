function test_settings_prog_sop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_SETTINGS_PROG_SOP tests the settings of an LPI program (LPIPROGRAM,
% LPISETSETTINGS) and the domain that POSCOPVAR and LPIVAR_CDOPVAR read
% from the program, and errors at the first failed check:
%   (1) LPIPROGRAM stores the settings given as its last input, a tier name
%       or a structure, in each input form (vars and dom; vars, dummies and
%       dom; with free variables; with decision and free variables), and
%       'light' without them; the program is otherwise that of the same
%       form without settings;
%   (2) LPISETSETTINGS sets a tier, a structure, and 'light' for []; an
%       unknown tier name and a 'double' are refused;
%   (3) on a program whose variables are declared out of sorted order, s2
%       on [0,2] before s1 on [0,1], POSCOPVAR with DOM omitted builds the
%       operator of the call with DOM given, on R^2 x L2[s1] x L2[s1,s2]
%       (an ODE, 1-D and 2-D coupling), and not that of the swapped domain;
%   (4) LPIVAR_CDOPVAR with DOM omitted, and with DOM and DEG omitted,
%       builds the operator of the explicit call (DEG: the caps of 'light'
%       block by block, [2 2 2] for L2[s1], 1 for L2[s1,s2]);
%   (5) LPI_INEQ_SOP takes its degree defaults from prog.settings: on the
%       target Pt of TEST_LPI_INEQ_SOP the weight degree is 2 at 'light'
%       and 1 at 'stripped' (wmin 0); a settings structure passed to the
%       call wins over prog.settings, and so does an explicit dw;
%   (6) LPISOLVE without options passes prog.settings.sos_opts, and with a
%       settings structure its sos_opts (an unknown solver name is reported
%       by SOSSOLVE in both cases); with the tier options it solves;
%   (7) POSCOPVAR without DEG builds, on R^2 x L2[s1] at 'light' and at
%       'heavy', the operator of the explicit specification of lpi.pos
%       (the multiplier, lower and upper bases), and on R^2 x L2[s1] x
%       L2[s1,s2] the multiplier and full-integral bases in every space
%       (s1 separable, being owned by the 2-D space); with DEG omitted an
%       OPTIONS.include is refused; LPIVAR_CDOPVAR without DEG on a 1-D
%       program at 'heavy' builds the operator of lpi.var.deg{1} = [4 2 2];
%   (8) EXEC_SOP_SETTINGS refuses a retired stock-format field (dd1) and
%       settings_nd.LF_int with the identifier exec_sop_settings:retired;
%       LPISETTINGS_STOCK_SOP('heavy') sets sop.slack = 'stock' and returns
%       sop.stock with the former dd2 of the heavy level, {3,[2 2 3],[2 2 3]};
%   (9) LPIVAR_CDOPVAR without DEG at 'heavy' (lpi.var.deg {[4 2 2],1,1})
%       on R x L2[s1] x L2[s1,s2]: the largest exponent of each variable
%       in the left basis of a block is the cap of the entry for its output
%       space, in the right basis that of its input space (L2[s1]: 4 left,
%       2 right in s1 shared, 4 right in s1 input-only; L2[s1,s2]: 1); on
%       R x L2[s1,s2] the operator of entry 2.
%
% INPUT:  none.
% OUTPUT: none; one line per passed part, then a summary line.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_settings_prog_sop
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

warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
npass = 0;

%% (1) LPIPROGRAM: the settings as the last input, in every input form
vars = polynomial({'s2';'s1'});     dum = polynomial({'s2_dum';'s1_dum'});
dom = [0 2; 0 1];                   free = polynomial({'t'});
light = lpisettings_sop();          heavy = lpisettings_sop('heavy');
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
        error('test_settings_prog_sop: (1) form %d stores the wrong settings.',k);
    end
    if ~isequaln(rmfield(a,'settings'),rmfield(b,'settings'))
        error('test_settings_prog_sop: (1) form %d: the settings changed the program.',k);
    end
end
fprintf('  passed: (1) settings stored in %d input forms, programs otherwise equal\n',size(pairs,1));
npass = npass+1;

%% (2) LPISETSETTINGS
custom = lpisettings_sop('custom');
ok = isequaln(lpisetsettings(p0,'veryheavy').settings,lpisettings_sop('veryheavy')) ...
     && isequaln(lpisetsettings(p0,custom).settings,custom) ...
     && isequaln(lpisetsettings(p0,[]).settings,light);
if ~ok
    error('test_settings_prog_sop: (2) LPISETSETTINGS stores the wrong settings.');
end
refused = 0;
try,    lpisetsettings(p0,'nosuchtier');    catch,  refused = refused+1;    end
try,    lpisetsettings(p0,5);               catch,  refused = refused+1;    end
if refused~=2
    error('test_settings_prog_sop: (2) an unknown tier or a double was accepted.');
end
fprintf('  passed: (2) tier, structure and [] set; unknown tier and double refused\n');
npass = npass+1;

%% (3) POSCOPVAR: DOM from the program, by name
sp = {{},{'s1'},{'s1','s2'}};       dm = [2;1;1];
[pA,PA] = poscopvar(p0,dm,sp,[0 1; 0 2],1);     % rows in sorted order s1, s2
[pB,PB] = poscopvar(p0,dm,sp,1);
[~,PC] = poscopvar(p0,dm,sp,[0 2; 0 1],1);      % the swapped domain
if ~isequal(PA,PB) || ~isequal(pA.decvartable,pB.decvartable)
    error('test_settings_prog_sop: (3) DOM omitted builds another operator.');
end
if isequal(PB,PC)
    error('test_settings_prog_sop: (3) the swapped domain gives the same operator; the check is void.');
end
fprintf('  passed: (3) poscopvar on R^2 x L2[s1] x L2[s1,s2]: DOM omitted = given (%d names), swapped differs\n',...
        numel(pB.decvartable));
npass = npass+1;

%% (4) LPIVAR_CDOPVAR: DOM, and DOM and DEG, from the program
dz = struct('out',[1;1],'in',[1;1]);    sz = {{'s1'},{'s1','s2'}};
[~,ZA] = lpivar_cdopvar(p0,dz,sz,[0 1; 0 2],1);
[~,ZB] = lpivar_cdopvar(p0,dz,sz,1);
[~,ZC] = lpivar_cdopvar(p0,dz,sz);
% DEG omitted at light: entry 1 ([2 2 2]) for L2[s1], entry 2 (1) for       % MMP, 11/12/2026
% L2[s1,s2]; left caps from the output space, right from the input space    % MMP, 11/12/2026
dgc = {[2 2 2], struct('mult',2,'int',[2 1],'out',2,'in',1); ...
       struct('mult',1,'int',[1 2],'out',1,'in',2), 1};                     % MMP, 11/12/2026
[~,ZD] = lpivar_cdopvar(p0,dz,sz,[0 1; 0 2],dgc);                           % MMP, 11/12/2026
% if ~isequal(ZA,ZB) || ~isequal(ZA,ZC)                                     % MMP, 11/12/2026 (was)
if ~isequal(ZA,ZB) || ~isequal(ZD,ZC)                                       % MMP, 11/12/2026
    error('test_settings_prog_sop: (4) DOM or DEG omitted builds another operator.');
end
fprintf('  passed: (4) lpivar_cdopvar: DOM omitted, and DOM and DEG omitted, = explicit\n');
npass = npass+1;

%% (5) LPI_INEQ_SOP: degree defaults from prog.settings
% Pt of TEST_LPI_INEQ_SOP on L2[0,1]: lower kernel 2 - s, upper 2 - t
v1 = {'s1'};    d1 = [0 1];     ZL = {[0;1]};   ZR = {[0;1]};
params = reshape({sparse(2,2),sparse([2 0;-1 0]),sparse([2 -1;0 0])},[3 1 1]);
Pt = copvar({sopvar(params,struct('out',{v1},'in',{v1}),ZL,ZR,struct('out',d1,'in',d1),[1 1])});
pL = lpiprogram(polynomial({'s1'}),[],d1);
pS = lpiprogram(polynomial({'s1'}),[],d1,'stripped');
[~,~,iL] = lpi_ineq_sop(pL,Pt);
[~,~,iS] = lpi_ineq_sop(pS,Pt);
[~,~,iX] = lpi_ineq_sop(pS,Pt,heavy);               % the call wins
[~,~,iD] = lpi_ineq_sop(pL,Pt,struct('dw',1));      % explicit dw: no floor
w = [iL.degrees.w, iS.degrees.w, iX.degrees.w, iD.degrees.w];
if ~isequal(w,[2 1 2 1])
    error('test_settings_prog_sop: (5) weights %s, expected [2 1 2 1] (light, stripped, call heavy, dw 1).',mat2str(w));
end
fprintf('  passed: (5) lpi_ineq_sop weight: light %d, stripped %d, call heavy %d, dw 1 %d\n',w);
npass = npass+1;

%% (6) LPISOLVE: the solver options of the program, or of a settings structure
[pE,g] = lpidecvar(pL,'g');
pE = lpi_ineq_sop(pE,g);    pE = lpisetobj(pE,g);
bad = pE;   bad.settings.sos_opts.solver = 'nosuchsolver';
stb = light;    stb.sos_opts.solver = 'nosuchsolver';
seen = 0;
try,    lpisolve(bad);      catch ME,   seen = seen + contains(ME.message,'nosuchsolver');    end
try,    lpisolve(pE,stb);   catch ME,   seen = seen + contains(ME.message,'nosuchsolver');    end
if seen~=2
    error('test_settings_prog_sop: (6) LPISOLVE did not pass the solver options of the settings.');
end
evalc('ps = lpisolve(pE);');
gv = double(lpigetsol_sop(ps,g));
if ~(abs(gv) < 1e-6)
    error('test_settings_prog_sop: (6) min g s.t. g >= 0 returned g = %g.',gv);
end
fprintf('  passed: (6) lpisolve: prog.settings.sos_opts and a settings structure passed; min g = %.1e\n',gv);
npass = npass+1;

%% (7) POSCOPVAR and LPIVAR_CDOPVAR without DEG: the program settings
q1 = lpiprogram(polynomial({'s1'}),[0 1]);          % light
q2 = lpiprogram(polynomial({'s1'}),[0 1],'heavy');
sp1 = {{},{'s1'}};      dm1 = [2;1];
for q = {q1, q2}
    ps = q{1}.settings.lpi.pos;
    Zi = struct('int',ps.int(1),'mult',ps.int(1),'joint',ps.joint(1));
    dx = {struct('int',0), {struct('int',ps.mult(1),'mult',0,'joint',ps.mult(1)),Zi,Zi}};
    ox = struct('sep',false,'include',{{[],[1;2;3]}});
    [qa,Pa] = poscopvar(q{1},dm1,sp1);
    [qb,Pb] = poscopvar(q{1},dm1,sp1,dx,ox);
    if ~isequaln(Pa,Pb) || ~isequaln(qa.decvartable,qb.decvartable)
        error('test_settings_prog_sop: (7) poscopvar without DEG differs from the explicit lpi.pos (mult %d).',ps.mult(1));
    end
end
% the coupled spaces: s1 separable for the 1-D space too
st2 = p0.settings.lpi.pos;
Zi1 = struct('int',st2.int(1),'mult',st2.int(1),'joint',st2.joint(1));
Zi2 = struct('int',st2.int(2),'mult',st2.int(2),'joint',st2.joint(2));
dx = {struct('int',0), {struct('int',st2.mult(1),'mult',0,'joint',st2.mult(1)),Zi1}, ...
      {struct('int',st2.mult(2),'mult',0,'joint',st2.mult(2)),Zi2}};
ox = struct('sep',true(1,2));   ox.include = {[],[1;4],[1 1;4 4]};
[~,Pc] = poscopvar(p0,dm,sp);
[~,Pd] = poscopvar(p0,dm,sp,dx,ox);
if ~isequaln(Pc,Pd)
    error('test_settings_prog_sop: (7) the coupled default differs from multiplier and full integral in every space.');
end
refused = 0;
try,    poscopvar(q1,dm1,sp1,[],struct('include',{{[],[1;2]}}));   catch,  refused = 1;    end
if ~refused
    error('test_settings_prog_sop: (7) OPTIONS.include with DEG omitted was accepted.');
end
dz1 = struct('out',[1;1],'in',[1;1]);   sz1 = {{},{'s1'}};
[~,Za] = lpivar_cdopvar(q2,dz1,sz1);
[~,Zb] = lpivar_cdopvar(q2,dz1,sz1,[4 2 2]);
if ~isequaln(Za,Zb)
    error('test_settings_prog_sop: (7) lpivar_cdopvar without DEG differs from lpi.var.deg{1}.');
end
fprintf('  passed: (7) poscopvar and lpivar_cdopvar without DEG: the settings at light and heavy, coupled spaces separable\n');
npass = npass+1;

%% (8) Retired fields refused; the stock-slack utility
bad = {setfield(lpisettings_sop(),'dd1',{1,[1 1 2],[1 1 2]}), ...
       setfield(lpisettings_sop(),'settings_nd',struct('LF_int',1))};
nref = 0;
for k = 1:numel(bad)
    try,    exec_sop_settings(bad{k});
    catch ME,   nref = nref + strcmp(ME.identifier,'exec_sop_settings:retired');
    end
end
[sk,stock] = lpisettings_stock_sop('heavy');
if nref~=2 || ~strcmp(sk.sop.slack,'stock') || ~isequal(stock,sk.sop.stock) || ...
        ~isequal(stock.dd2,{3,[2 2 3],[2 2 3]})
    error('test_settings_prog_sop: (8) retired fields accepted (%d of 2 refused) or the stock slack of heavy wrong.',nref);
end
fprintf('  passed: (8) retired fields refused (2 of 2); the stock slack of heavy from LPISETTINGS_STOCK_SOP\n');
npass = npass+1;

%% (9) LPIVAR_CDOPVAR without DEG: the caps of each block from its two spaces
% heavy: entry 1 [4 2 2] (mult = out = in = 4, int = [2 2]), entry 2 every
% cap 1. Rows: block, side, its variables, their largest exponents.
h = lpiprogram(polynomial({'s1';'s2'}),[0 1; 0 2],'heavy');
sp3 = {{},{'s1'},{'s1','s2'}};      d3 = struct('out',[1;1;1],'in',[1;1;1]);
[~,Zh] = lpivar_cdopvar(h,d3,sp3);
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
        error('test_settings_prog_sop: (9) block (%d,%d) %s over %s: largest exponents %s, expected %s.', ...
              want{k,1},want{k,2},strjoin(v,','),mat2str(cellfun(@max,Zb)),mat2str(want{k,4}));
    end
end
% One dimension besides R^n: every block from that entry, as before
dz2 = struct('out',[1;1],'in',[1;1]);   sz2 = {{},{'s1','s2'}};
[~,Ze] = lpivar_cdopvar(h,dz2,sz2);
[~,Zf] = lpivar_cdopvar(h,dz2,sz2,1);
if ~isequal(Ze,Zf)
    error('test_settings_prog_sop: (9) on R x L2[s1,s2] the default differs from entry 2.');
end
fprintf('  passed: (9) lpivar_cdopvar without DEG: %d block sides capped by their own space; R x L2[s1,s2] = entry 2\n',size(want,1));
npass = npass+1;

% The summary line
fprintf('\ntest_settings_prog_sop passed (%d checks).\n',npass);
end
