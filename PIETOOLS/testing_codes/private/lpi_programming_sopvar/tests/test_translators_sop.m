function R = test_translators_sop()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_TRANSLATORS_SOP() checks poslpivar_sop, poslpivar_settings_sop,
% get_lpivar_degs_sop, trace_rn_sop and copvar_space_list in one spatial
% variable s1 on [0,1], against their definitions, and errors at the first
% failed check:
%   (i)   none since 11/08/2026: it compared the translators with the
%         stock poslpivar and get_lpivar_degs on legacy inputs, retired;
%   (ii)  each of 12 inputs built to fail raises an error, with the
%         identifier expected where one is listed (poslpivar_sop:side,
%         :dim, :excludeRn, :excludeL2, :nargout;
%         poslpivar_settings_sop:stock, :role; get_lpivar_degs_sop:class;
%         trace_rn_sop:square; copvar_space_list:side);
%   (iii) on the container R x L2[s1], with lpisettings_stock_sop('custom')
%         and override2 = 0, poslpivar_settings_sop returns for the role
%         'lf' the program and operator of the former stock-format
%         dd1 = {1,[1 1 2],[1 1 2]}, sep 0 (poslpivar_sop), and for the
%         role 'slack' those of its two poslpivar_sop terms summed;
%   (iv)  with the integral bases excluded (exclude [0 0 1 1]),
%         get_lpivar_degs_sop returns 0 as third degree;
%   (v)   at random decision values, trace_rn_sop equals the sum of the
%         traces of the R^q diagonal blocks of the substituted operator, to
%         1e-12 (1 + |trace|), for a free operator on R^3, a free operator
%         W on R^2 x L2[s1], Mz*W*Mz' (Mz a random fixed operator), the
%         identity on R^3 x L2[s1] and the identity on L2^2[s1] (trace 0);
%   (vi)  none since 11/08/2026: it built its operators as 'opvar'
%         objects (opvar2copvar), retired.
%
% INPUT
% none.
%
% OUTPUT
% - R:      'struct' with field checks ('double', the number of checks
%           passed); one printed line with that number.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_translators_sop
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
% Initial coding MMP, 10/06/2026
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and headers of the local functions; code
%                unchanged. The deleted header text was a history note on
%                the copies these routines replace and a reference to
%                another test.
% MMP, 11/04/2026: Programs declared by LPIPROGRAM; LPIPROGRAM_SOP, whose
%                extra input forms (names, containers) only filled fields that
%                container routines do not read, is retired.
% MMP, 11/08/2026: Settings from lpisettings_sop (was the stock lpisettings,
%                whose folder leaves the path); the preset loop of (iv)
%                covers the same six names, 'extreme' now a tier there.
% MMP, 11/08/2026: (iv) Z2, Z3 excluded: GET_LPIVAR_DEGS_SOP now returns
%                max(deg2-1,0), so the check expects 0 and compares with the
%                stock degrees floored at 0 (the stock rule gives -1).
% MMP, 11/08/2026: The comparisons with the stock routines (10/06/2026)
%                deleted: (i), the two message checks of (ii), the stock
%                loop of (iv) and its stock comparison, and (vi), which
%                built its operators as opvars. POSLPIVAR, GET_LPIVAR_DEGS,
%                LPIVAR and the opvar classes are retired (cleaning/
%                private). (iv) keeps the container check of the third
%                degree; the settings of (iii) are set before (ii). Help to
%                match.
% MMP, 11/08/2026: (ii) case 10 passes a double to GET_LPIVAR_DEGS_SOP and
%                expects :class. It passed the container X2 and expected
%                :dim, which the N-D reading (10/08/2026) replaced by a
%                returned struct. Help to match.
% MMP, 11/08/2026: Comment pass (maintainer, 11/08/2026): three (was)
%                comment lines removed; code unchanged.
% MMP, 11/12/2026: Native settings (phase 2 of the knob reduction,
%                maintainer 11/12/2026): (ii) checks the error of a 'slack'
%                without sop.stock in place of sosineq_on, retired; (iii)
%                compares 'lf' with the former dd1 of 'custom', so that the
%                native lpi.pos is checked against what the tier held, and
%                'slack' with the two terms of sop.stock
%                (LPISETTINGS_STOCK_SOP). override1 and the 'lf' psatz term
%                are gone.

% No warnings for non-canonical multipliers; seed; s1 on [0,1]; check
% counter
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
rng(20261006);
pvar s1 s1_dum
dom = [0 1];
n = 0;

% Deleted: (i) (MMP, 10/06/2026), poslpivar_sop, poslpivar_settings_sop,    % MMP, 11/08/2026
% get_lpivar_degs_sop and trace_rn_sop against POSLPIVAR, GET_LPIVAR_DEGS   % MMP, 11/08/2026
% and LPIVAR on legacy inputs; all retired. Kept: the settings, used by     % MMP, 11/08/2026
% (ii) and (iii).                                                           % MMP, 11/08/2026
% st = lpisettings('custom');                                               % MMP, 11/08/2026 (was)
% st = lpisettings_sop('custom');                                           % MMP, 11/12/2026 (was); MMP, 11/08/2026
% st.override1 = 0;   st.override2 = 0;           % both second terms declared % MMP, 11/12/2026 (was)
st = lpisettings_stock_sop('custom');               % with the stock slack  % MMP, 11/12/2026
st.sop.stock.override2 = 0;                         % its second term       % MMP, 11/12/2026

% % % (ii) errors
% Containers: R^1 x L2[s1], L2[s1] -> R^2, one space over two variables
Xm = eye_copvar_sop([1;1],{{},{'s1'}},dom);         % R^1 x L2[s1]
Xr = rand_copvar(struct('out',{{{}}},'in',{{{'s1'}}}),struct('out',2,'in',1),dom,1,0.8);  % L2 -> R^2
X2 = eye_copvar_sop(1,{{'s1','s2'}},[0 1;0 1]);     % one space over two variables
pc = lpiprogram(s1,s1_dum,dom);
E = { @() poslpivar_sop(pc,Xr,1,[]),                                   'poslpivar_sop:side'
%       @() poslpivar_sop(lpiprogram_sop({'s1','s2'},[0 1;0 1]),X2,1,[]), 'poslpivar_sop:dim' % MMP, 11/04/2026 (was)
      @() poslpivar_sop(lpiprogram(polynomial({'s1';'s2'}),[0 1;0 1]),X2,1,[]), 'poslpivar_sop:dim' % MMP, 11/04/2026
      @() poslpivar_sop(pc,Xm,1,struct('exclude',[1 0 0 0])),         'poslpivar_sop:excludeRn'
      @() poslpivar_sop(pc,Xm,1,struct('exclude',[0 1 1 1])),         'poslpivar_sop:excludeL2'
      @() nout3(pc,Xm),                                                'poslpivar_sop:nargout'
      @() poslpivar_sop(pc,Xm,[1 2 3 4],[]),                           ''
      @() poslpivar_sop(pc,Xm,1,7),                                    ''
      @() poslpivar_settings_sop(pc,Xm,lpisettings_sop('custom'),'slack'), 'poslpivar_settings_sop:stock' % MMP, 11/12/2026
      @() poslpivar_settings_sop(pc,Xm,st,'other'),                    'poslpivar_settings_sop:role'
%       @() get_lpivar_degs_sop(X2),                                     'get_lpivar_degs_sop:dim' % MMP, 11/08/2026 (was)
      @() get_lpivar_degs_sop(eye(2)),                                 'get_lpivar_degs_sop:class' % MMP, 11/08/2026
      @() trace_rn_sop(rand_copvar(struct('out',{{{}}},'in',{{{}}}),struct('out',3,'in',2),dom,0,0.8)), 'trace_rn_sop:square'
      @() copvar_space_list(Xm,'both'),                                'copvar_space_list:side' };
% Each case raises an error, with the listed identifier where one is given
for k = 1:size(E,1)
    id = errid(E{k,1});
    assert(~isempty(id) || ~isempty(errmsg(E{k,1})),'(ii) case %d raised no error',k);
    if ~isempty(E{k,2})
        assert(strcmp(id,E{k,2}),'(ii) case %d: got ''%s'', expected ''%s''',k,id,E{k,2});
    end
    n = n+1;
end
% Deleted: cases 6 and 7 against the messages of POSLPIVAR (MMP,            % MMP, 11/08/2026
% 10/06/2026), retired; the loop above still requires an error for each.    % MMP, 11/08/2026

% % % (iii) poslpivar_settings_sop composition, container
% % For each role, its two settings terms by poslpivar_sop, summed          % MMP, 11/12/2026 (was)
% 'lf': the one term of the former dd1 of 'custom'; 'slack': the two        % MMP, 11/12/2026
% terms of sop.stock by poslpivar_sop, summed                               % MMP, 11/12/2026
for role = {'lf','slack'}
%   if strcmp(role{1},'lf'),  f = {'dd1','options1','dd12','options12'};    % MMP, 11/12/2026 (was)
%   else,                     f = {'dd2','options2','dd3','options3'};      % MMP, 11/12/2026 (was)
%   end                                                                     % MMP, 11/12/2026 (was)
    [pa,Pa] = poslpivar_settings_sop(pc,Xm,st,role{1},'out',dom);
%   [pb,P1] = poslpivar_sop(pc,Xm,st.(f{1}),st.(f{2}),'out',dom);           % MMP, 11/12/2026 (was)
%   [pb,P2] = poslpivar_sop(pb,Xm,st.(f{3}),st.(f{4}),'out',dom);           % MMP, 11/12/2026 (was)
%   [tf,where] = same_val_sop({pa,Pa},{pb,P1+P2},role{1});                  % MMP, 11/12/2026 (was)
    if strcmp(role{1},'lf')                                                 % MMP, 11/12/2026
        [pb,Pb] = poslpivar_sop(pc,Xm,{1,[1 1 2],[1 1 2]}, ...
                                struct('sep',0,'exclude',[0 0 0 0]),'out',dom); % MMP, 11/12/2026
    else                                                                    % MMP, 11/12/2026
        sk = st.sop.stock;                                                  % MMP, 11/12/2026
        [pb,P1] = poslpivar_sop(pc,Xm,sk.dd2,sk.options2,'out',dom);        % MMP, 11/12/2026
        [pb,P2] = poslpivar_sop(pb,Xm,sk.dd3,sk.options3,'out',dom);        % MMP, 11/12/2026
        Pb = P1+P2;                                                         % MMP, 11/12/2026
    end                                                                     % MMP, 11/12/2026
    [tf,where] = same_val_sop({pa,Pa},{pb,Pb},role{1});                     % MMP, 11/12/2026
    assert(tf,'(iii) %s pair differs at %s',role{1},where);    n = n+1;
end

% % % (iv) the third Q degree without integral bases                        % MMP, 11/08/2026
% Deleted: the loop over two plants and six presets comparing               % MMP, 11/08/2026
% get_lpivar_degs_sop with GET_LPIVAR_DEGS on the stock operators, and the  % MMP, 11/08/2026
% stock side of the check below (MMP, 10/06/2026; GET_LPIVAR_DEGS and the   % MMP, 11/08/2026
% opvar classes are retired).                                               % MMP, 11/08/2026
[~,R0] = poslpivar_sop(pc,Xm,1,struct('exclude',[0 0 1 1]));      % no R1, R2: deg2 = 0
q0 = get_lpivar_degs_sop(R0);
assert(q0(3)==0,'(iv) Z2, Z3 excluded: container %s (expected third degree 0)',mat2str(q0)) % MMP, 11/08/2026
n = n+1;

% % % (v) trace against its definition
% The five operators: free on R^3, free W on R^2 x L2[s1], Mz*W*Mz', and
% two identities
[pw,W3] = lpivar_cdopvar(pc,struct('out',3,'in',3),struct('out',{{{}}},'in',{{{}}}),dom,[1 1 1]);
[pw,Wmx] = lpivar_cdopvar(pw,struct('out',[2;1],'in',[2;1]), ...
                          struct('out',{{{},{'s1'}}},'in',{{{},{'s1'}}}),dom,[1 1 1]);
Mz = rand_copvar(struct('out',{{{}}},'in',{{{},{'s1'}}}),struct('out',1,'in',[2;1]),dom,1,0.8);
Xs = {W3, Wmx, (Mz*Wmx)*Mz', eye_copvar_sop([3;1],{{},{'s1'}},dom), eye_copvar_sop(2,{{'s1'}},dom)};
for k = 1:numel(Xs)
    % trace_rn_sop, and the operator and its trace at random values
    X = Xs{k};
    tr = trace_rn_sop(X);
    if isa(X,'cdopvar')
        v = randn(numel(X.Zd),1);
        Xv = subs_dvar_sop(X,X.Zd,v);    trv = double(subs_dvar_sop(tr,X.Zd,v));
    else
        Xv = X;     trv = double(tr.C(1));
    end
    % The definition: traces of the diagonal blocks on spaces without a
    % spatial variable
    ref = 0;
    for kk = find(~any(Xv.space_out,2)).'
        if ~isempty(Xv.C{kk,kk}),   ref = ref + trace(full(Xv.C{kk,kk}.params{1}));   end
    end
    assert(abs(trv-ref) <= 1e-12*(1+abs(ref)),'(v) case %d: trace %g, definition %g',k,trv,ref);
    n = n+1;
end

% Deleted: (vi) (MMP, 10/06/2026), copvar_space_list on the 'copvar' of     % MMP, 11/08/2026
% an 'opvar' (opvar2copvar) for each pattern of zero dimensions; the opvar  % MMP, 11/08/2026
% classes are retired.                                                      % MMP, 11/08/2026

% Summary and output
fprintf('test_translators_sop: %d checks passed.\n',n);
R = struct('checks',n);
end


function nout3(p,X)
% NOUT3(P,X) calls poslpivar_sop(P,X,1,[]) with three outputs.
% INPUT  p: 'struct', an LPI program; X: 'copvar' object.
% OUTPUT none.
[~,~,~] = poslpivar_sop(p,X,1,[]);
end

function id = errid(f)
% ID = ERRID(F) returns the identifier of the error that F() raises (''
% when it raises none).
% INPUT  f: 'function_handle', called with no argument.
% OUTPUT id: 'char'.
id = '';
try,        f();
catch ME,   id = ME.identifier;
end
end

function m = errmsg(f)
% M = ERRMSG(F) returns the message of the error that F() raises ('' when
% it raises none).
% INPUT  f: 'function_handle', called with no argument.
% OUTPUT m: 'char'.
m = '';
try,        f();
catch ME,   m = ME.message;
end
end
