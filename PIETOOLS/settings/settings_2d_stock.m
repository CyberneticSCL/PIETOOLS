function s2 = settings_2d_stock(level)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S2 = SETTINGS_2D_STOCK(LEVEL) returns the structure settings_2d of the
% level LEVEL: the degrees dx, dy, d2 and the options of the construction
% of the positive operator P in two variables from 16 bases, the degrees
% of the free operator, and the slack of the 'stock' slack mode. As the
% field settings_2d of LPI settings it replaces, in two variables, the
% storage row lpi.pos. The values are the 2-D degree levels of M. Peet,
% S. Shivakumar and D. Jagt (2021-2022).
%
% INPUT
% - level:  (optional) char, 'light' (default), 'stripped', 'heavy',
%           'veryheavy', 'extreme' or 'custom' (= 'light');
%
% OUTPUT
% - s2:     'struct' with fields is2D, LF_deg, LF_opts, LF_use_psatz,
%           LF_opts_psatz, LF_deg_psatz (P); Zop_deg, Zop_opts (free
%           operator); use_sosineq and ineq_opts or eq_deg, eq_opts,
%           eq_use_psatz, eq_opts_psatz, eq_deg_psatz ('stock' slack).
%           The body gives what each field does.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - settings_2d_stock
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
% Initial coding MMP, 11/17/2026

% LEVEL 'light' when omitted
if nargin<1 || isempty(level),  level = 'light';    end

% P = P0 + P1 + ..., each term Z'*Q*Z with Q >= 0 (role 'lf'), Z built
% from the 16 bases
% [R^n; Zx^o, Zx^a, Zx^b; Zy^o, Zy^a, Zy^b; Z2^oo, Z2^ao, Z2^bo, Z2^oa,
% Z2^ob, Z2^aa, Z2^ba, Z2^ab, Z2^bb]: o the multiplier, a the lower and b
% the upper integral, in x for the first index and in y for the second;
% tt and nu the integration variables in x and y.
% Degrees: dx on L2[x], dy on L2[y], d2 on L2[x,y]. The caps of Z2^{ij}
% have rows [-; x] (i = o) or [-; x; tt; x*tt] (i = a, b) and columns
% [-, y] (j = o) or [-, y, nu, y*nu] (j = a, b), each entry the cap on the
% joint degree in its row and column variables.
% Per level: sep, exclude (LF_opts); psatz terms with zero sep/exclude
% (pz0) or those of LF_opts; Zop_opts isscalar, ismultiplier; use_sosineq
% default and ineq_opts.psatz; raises [x y 2] of eq_deg and eq_deg_psatz;
% eq_use_psatz.
d2 = cell(3,3);
switch level
    case {'light','stripped','custom'}      % one level for the three
        dx = {1; [1;1;2]; [1;1;2]};         % Zx^o: degree in x; Zx^a, Zx^b: caps [x; tt; x*tt]
        dy = {1, [1 1 2], [1 1 2]};         % the same in y and nu
        d2{1,1} = [0 0;0 0];                                % Z2^oo
        d2{1,2} = [0 1 1 2;1 1 1 2];    d2{1,3} = d2{1,2};  % Z2^oa, Z2^ob
        d2{2,1} = [0 1;1 1;1 1;2 2];    d2{3,1} = d2{2,1};  % Z2^ao, Z2^bo
        d2{2,2} = [0 1 1 2;1 1 1 2;1 1 1 2;2 2 2 2];        % Z2^aa, ab, ba, bb
        sep = [1;1;1;1;1;1];    excl = [0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0];
        pz0 = true;     zsc = 0;    zmu = 0;    usi = 1;    ips = 2;
        up_eq = [3 3 3];    up_ps = [2 2 2];    eup = [0;2];
    case 'heavy'                            % the heavy level
        dx = {1; [1;1;1]; [1;1;1]};
        dy = {1, [1 1 1], [1 1 1]};
        d2{1,1} = [0 1;1 2];
        d2{1,2} = [0 1 1 1;1 2 2 2];    d2{1,3} = d2{1,2};
        d2{2,1} = [0 1;1 2;1 2;1 2];    d2{3,1} = d2{2,1};
        d2{2,2} = [0 1 1 1;1 2 2 2;1 2 2 2;1 2 2 2];
        sep = [1;1;1;1;1;1];    excl = [0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0];
        pz0 = false;    zsc = 0;    zmu = 1;    usi = 0;    ips = 0;
        up_eq = [2 2 3];    up_ps = [1 1 2];    eup = [0;0];
    case 'veryheavy'                        % the veryheavy level
        dx = {1; [2;2;3]; [2;2;3]};
        dy = {1, [2 2 3], [2 2 3]};
        d2{1,1} = [0 1;1 2];
        d2{1,2} = [0 2 2 3;1 1 2 4];    d2{1,3} = d2{1,2};
        d2{2,1} = [0 1;2 3;1 2;2 4];    d2{3,1} = d2{2,1};
        d2{2,2} = [0 2 2 4;2 2 2 4;2 2 2 4;4 4 4 4];
        sep = [0;0;1;1;1;1];    excl = [0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0];
        pz0 = true;     zsc = 0;    zmu = 0;    usi = 0;    ips = 2;
        up_eq = [3 3 3];    up_ps = [2 2 2];    eup = [0;2];
    case 'extreme'                          % the extreme level
        dx = {0; [1;1;1]; [1;1;1]};
        dy = {0, [1 1 1], [1 1 1]};
        d2{1,1} = [0 0;0 0];
        d2{1,2} = [0 1 1 1;0 1 1 1];    d2{1,3} = d2{1,2};
        d2{2,1} = [0 0;1 1;1 1;1 1];    d2{3,1} = d2{2,1};
        d2{2,2} = [0 1 1 1;1 2 2 2;1 2 2 2;1 2 2 2];
        sep = [1;1;1;1;1;1];    excl = [0 0 1 1 0 1 1 0 1 1 1 1 1 1 1 1];
        pz0 = false;    zsc = 1;    zmu = 1;    usi = 0;    ips = 0;
        up_eq = [2 2 3];    up_ps = [1 1 2];    eup = [0;0];
    otherwise
        error('settings_2d_stock:level',['LEVEL should be ''light'', ''stripped'', ' ...
              '''heavy'', ''veryheavy'', ''extreme'' or ''custom''.']);
end
d2{2,3} = d2{2,2};  d2{3,2} = d2{2,2};  d2{3,3} = d2{2,2};  % equal at every level

% Options: exclude (1x16) drops the bases in the order above; sep (6x1)
% equates the lower- and upper-integral parameters of P: Rxx{2} = Rxx{3}
% (1), Ryy{2} = Ryy{3} (2), R22{2,1} = R22{3,1} (3), R22{1,2} = R22{1,3}
% (4), R22{2,j} = R22{3,j} for j = 2,3 (5), R22{i,2} = R22{i,3} for
% i = 2,3 (6), with R22{i,j} the parameter of Z2^{ij}; psatz 1 weights Q
% by the domain, 2 by a disk (rejected by the declaration of P). Each nonzero entry
% of LF_use_psatz adds a term of the options LF_opts_psatz{j} and the
% degrees LF_deg_psatz{j}.
s2.is2D = 1;
s2.LF_deg  = struct('dx',{dx},'dy',{dy},'d2',{d2});
s2.LF_opts = struct('sep',sep,'exclude',excl,'psatz',0);
s2.LF_use_psatz  = 0;                               % no added term: the two fields below unread
if pz0,     sp = zeros(1,6);    ex = zeros(1,16);
else,       sp = sep;           ex = excl;      end
s2.LF_opts_psatz = {struct('psatz',1,'sep',sp,'exclude',ex), ...
                    struct('psatz',2,'sep',sp,'exclude',ex)};
s2.LF_deg_psatz  = {s2.LF_deg, s2.LF_deg};
% Free operator of the synthesis LPIs: twice the degrees of P, of which only
% the largest entry is read, and its options
twice = @(c) cellfun(@(d) 2*d, c, 'UniformOutput', false);
s2.Zop_deg  = struct('dx',{twice(dx)},'dy',{twice(dy)},'d2',{twice(d2)});
s2.Zop_opts = struct('sep',[0;0;0;0;0],'isscalar',zsc,'ismultiplier',zmu);
% Slack of the 'stock' slack mode (role 'slack'), read only with sop.slack
% = 'stock'. use_sosineq is read from the base workspace (sosineq_on, else
% use_sosineq) when it exists there: 1 selects an SOS inequality, which
% the executives reject, and sets only ineq_opts. With 0, the slack has
% the degrees of P raised by up_eq in dx, dy and d2, and one term per
% nonzero entry of eq_use_psatz with those raised by up_ps.
s2.use_sosineq = usi;
if evalin('base','exist(''sosineq_on'',''var'')')
    s2.use_sosineq = evalin('base','sosineq_on');
elseif evalin('base','exist(''use_sosineq'',''var'')')
    s2.use_sosineq = evalin('base','use_sosineq');
end
if s2.use_sosineq
    s2.ineq_opts = struct('psatz',ips);
else
    up = @(c,k) cellfun(@(d) d+k, c, 'UniformOutput', false);
    s2.eq_deg  = struct('dx',{up(dx,up_eq(1))},'dy',{up(dy,up_eq(2))},'d2',{up(d2,up_eq(3))});
    s2.eq_opts = struct('psatz',0,'sep',zeros(1,6),'exclude',zeros(1,16));
    s2.eq_use_psatz = eup;
    for j = 1:numel(s2.eq_use_psatz)
        s2.eq_opts_psatz{j} = struct('psatz',s2.eq_use_psatz(j),'sep',zeros(1,6),'exclude',zeros(1,16));
        s2.eq_deg_psatz{j}  = struct('dx',{up(dx,up_ps(1))},'dy',{up(dy,up_ps(2))},'d2',{up(d2,up_ps(3))});
    end
end
end
