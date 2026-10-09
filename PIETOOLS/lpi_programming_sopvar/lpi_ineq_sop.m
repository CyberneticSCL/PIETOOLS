function [prog,Pop,info] = lpi_ineq_sop(prog,P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP,INFO] = LPI_INEQ_SOP(PROG,P,OPTS) adds to the LPI program PROG
% the constraint P >= 0 for a self-adjoint PI operator P of either family.
% It is 'lpi_ineq' with the container family added:
%   class of P                              handled by
%   'cdopvar', 'copvar'                     this file: P = Pop, Pop >= 0
%   anything else ('dopvar', 'dpvar', ...)  lpi_ineq, unchanged
%
% CONTAINER PATH. The positive operator is sized from P by the lift/weight
% rules of the separated form, GET_LIFT_DEGS, and declared in ONE call of
% POSCOPVAR_DIRECT with every Positivstellensatz term as a Gram-side shift
% on a shared map (the multipliers are added before the lift is applied),
% then P - Pop == 0 is imposed by LPI_EQ_SOP ('symmetric'). Measured on
% the heat benchmark (heatNd_tailor, 10/08/2026): the defaults reach the
% exact 1-D rate with 2.25 times fewer SDP variables than the stock sizing,
% and in 2-D the operator-wise degrees save 31% of the variables at the
% same rate while the face terms stay necessary.
%
% INPUT
% - prog:  'struct', the LPI program ('lpiprogram_sop' or 'lpiprogram');
% - P:     the operator to make positive: 'cdopvar' (decision), or 'copvar'
%          (fixed; then the constraint certifies P >= 0); square, assumed
%          self-adjoint, as 'lpi_ineq' assumes;
% - opts:  (optional) struct with fields
%   dD, dw, wR, psatz, psatz_offset, like   as GET_LIFT_DEGS ('like': size
%                from the fixed positive operator of the LPI instead of from
%                P, for a P linear in a free operator such as a KYP slack:
%                the Q forms want dw 1, the coercive form dw 2, measured
%                10/08/2026, README);
%   deg          a degree specification for 'poscopvar_direct' used as
%                given instead of the reader's (then psatz, psatz_offset
%                as given, default [0 1] with [0 1] in 1-D, the faces at
%                w-1 in N-D when every weight is at least 2, at w
%                otherwise: GET_LIFT_DEGS 'auto', changed 10/08/2026);
%   prune        logical (default true): restrict the basis operators of
%                each L2 space to those that can reach the support of the
%                diagonal block of P ('eq_opts_sopvar', lossless by its
%                argument);
%   path, cache  as 'poscopvar_direct';
%   margin       eps >= 0 (default 0): certify P + eps*e_D >= 0 instead,
%                e_D = A'A the energy of the plain term's lift, the
%                regularizer of the proof program (a solver margin; the
%                certified statement is the perturbed one: INFO.margin, eD).
%   For a legacy P, OPTS is passed to 'lpi_ineq' unchanged.
%
% OUTPUT
% - prog:  the program with the Grams declared and P - Pop == 0 imposed;
% - Pop:   the positive 'cdopvar' (container path; [] for a legacy P);
% - info:  struct with fields deg (the specification used), terms (codes,
%          offsets), degrees (GET_LIFT_DEGS info), include, direct (the
%          info of 'poscopvar_direct'), tag (the row record of the
%          equality in prog.sopeq, for LPIGETDUAL_SOP), margin and eD.
%
% NOTES
% Why one call: the SDP is the same whether the terms are declared
% together or one 'poscopvar' at a time and summed, but the sum re-merges a
% q-length decision list per term ('plus_batch': 6 s of the 56 s Q stage
% of the 3-D heat benchmark) and rebuilds a map per term; the one-call
% route builds each map once, shares it between terms on the same index
% set and the operator once (POSCOPVAR_DIRECT, 10/08/2026).
%
% See also LPI_INEQ, LPI_EQ_SOP, POSCOPVAR_DIRECT, GET_LIFT_DEGS,
%          EQ_OPTS_SOPVAR, LPIPROGRAM_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpi_ineq_sop
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
% Initial coding MMP, 10/08/2026. Tier 2 of the container parity map: the
%                negativity constraint on the direct map, degrees from the
%                target by GET_LIFT_DEGS.
% MMP, 10/08/2026: opts.margin (the e_D regularizer) and info.tag (the row
%                record of the equality, for the dual read-back), for the
%                container executives.
% MMP, 10/08/2026: The pruned 'include' expands the full-integral index 4
%                (eq_opts_sopvar marks a direction separable when the
%                diagonal block's lower and upper cells agree, as the stock
%                2-D storage's do under options.sep) to the indices 2 and 3
%                of the unseparated lift; the 2-D H2 Q forms, whose
%                positivity block carries the storage itself, failed in
%                alpha_select ("inadmissible multi-index") before this.

if nargin<2
    error("Not enough input arguments.")
end
if nargin<3,    opts = [];  end

% % % Legacy classes: 'lpi_ineq', unchanged.
if ~(isa(P,'cdopvar') || isa(P,'copvar'))
    if isempty(opts),   prog = lpi_ineq(prog,P);
    else,               prog = lpi_ineq(prog,P,opts);
    end
    Pop = [];   info = struct('legacy',true);
    return
end
if isempty(opts),   opts = struct();    end
if ~isstruct(opts)
    error("Options should be specified as a 'struct' object.")
end
if ~isequal(P.space_out,P.space_in) || ~isequal(P.dim_out(:),P.dim_in(:))
    error('lpi_ineq_sop:square',"P must be square: the same spaces on both sides.")
end

% % % Spaces, as poscopvar takes them.
[sp,dm] = copvar_space_list(P,'out');
dom = struct('vars',{P.vars},'dom',P.dom);
M = numel(sp);

% % % Degrees and terms: the reader, or the given specification.
if isfield(opts,'deg') && ~isempty(opts.deg)
    deg = opts.deg;     dinfo = struct('given',true);
    nv = numel(P.vars);
    if isfield(opts,'psatz') && ~isempty(opts.psatz) && isnumeric(opts.psatz)
        codes = reshape(opts.psatz,1,[]);   offs = zeros(size(codes));
        if isfield(opts,'psatz_offset') && ~isempty(opts.psatz_offset)
            offs = reshape(opts.psatz_offset,1,[]);
            if isscalar(offs),  offs = repmat(offs,size(codes));    end
        end
    else
        ro = opts;  ro.deg = [];
        [~,t0,~] = get_lift_degs(P,ro);     codes = t0.codes;   offs = t0.offsets;
    end
    terms = struct('codes',codes,'offsets',offs);
else
    [deg,terms,dinfo] = get_lift_degs(P,opts);
end

% % % The basis operators each L2 space needs: those whose diagonal term can
% reach the support of P's diagonal block.
prune = true;
if isfield(opts,'prune') && ~isempty(opts.prune),   prune = logical(opts.prune);  end
include = cell(1,M);
if prune
    for k = 1:M
        if isempty(sp{k}) || isempty(P.C{k,k}),    continue,   end
        eo = eq_opts_sopvar(P.C{k,k});
%       include{k} = eo.include;                                            % MMP, 10/08/2026 (was)
        inc = eo.include;                                                   % MMP, 10/08/2026
        % eq_opts_sopvar marks a direction separable (its lower and upper   % MMP, 10/08/2026
        % cells agree, as the stock storage's do under options.sep) and     % MMP, 10/08/2026
        % lists the full-integral index 4 there; the lift here is not       % MMP, 10/08/2026
        % separated, so 4 is expanded to the lower and upper indices 2 and  % MMP, 10/08/2026
        % 3 (a superset of the unseparated reading: lossless).              % MMP, 10/08/2026
        for t = 1:size(inc,2)                                               % MMP, 10/08/2026
            r = inc(:,t)==4;                                                % MMP, 10/08/2026
            if any(r)                                                       % MMP, 10/08/2026
                lo = inc(r,:);  up = lo;    lo(:,t) = 2;    up(:,t) = 3;    % MMP, 10/08/2026
                inc = unique([inc(~r,:); lo; up],'rows');                   % MMP, 10/08/2026
            end                                                             % MMP, 10/08/2026
        end                                                                 % MMP, 10/08/2026
        include{k} = inc;                                                   % MMP, 10/08/2026
    end
end

% % % One positive operator, all terms, then the equality.
popts = struct('psatz',terms.codes,'psatz_offset',terms.offsets);
if prune && any(~cellfun(@isempty,include)),   popts.include = include;    end
for f = {'path','cache'}
    if isfield(opts,f{1}) && ~isempty(opts.(f{1})),  popts.(f{1}) = opts.(f{1});    end
end
[prog,Pop,~,pinfo] = poscopvar_direct(prog,dm,sp,dom,deg,popts);
% % % The e_D margin (opts.margin = eps > 0): certify P + eps*e_D >= 0 with  % MMP, 10/08/2026
% e_D = A'A the energy of the plain term's own lift (lift_copvar), the      % MMP, 10/08/2026
% perturbation under which the proof program guarantees a polynomial weight % MMP, 10/08/2026
% at a fixed lift (sopvar_lift_notes Sec. 7.3). A solver regularisation:    % MMP, 10/08/2026
% the certified statement is the perturbed one, and the caller reports it.  % MMP, 10/08/2026
eD = [];    margin = 0;                                                     % MMP, 10/08/2026
if isfield(opts,'margin') && ~isempty(opts.margin) && opts.margin>0         % MMP, 10/08/2026
    margin = opts.margin;                                                   % MMP, 10/08/2026
    lopt = struct();                                                        % MMP, 10/08/2026
    if isfield(popts,'include'),    lopt.include = popts.include;   end     % MMP, 10/08/2026
    dl = 1;                         % the plain term's lift degree          % MMP, 10/08/2026
    dc = deg;   if ~iscell(dc),     dc = {dc};  end                         % MMP, 10/08/2026
    for k = 1:numel(dc)                                                     % MMP, 10/08/2026
        if isstruct(dc{k}) && isfield(dc{k},'mult') && ~isempty(dc{k}.mult) % MMP, 10/08/2026
            dl = max([dl, reshape(dc{k}.mult,1,[])]);                       % MMP, 10/08/2026
        elseif isnumeric(dc{k}) && isscalar(dc{k}),  dl = max(dl,dc{k});   % MMP, 10/08/2026
        end                                                                 % MMP, 10/08/2026
    end                                                                     % MMP, 10/08/2026
    A_D = lift_copvar(dm,sp,dom,dl,lopt);                                   % MMP, 10/08/2026
    eD = A_D'*A_D;                                                          % MMP, 10/08/2026
    P = P + margin*eD;                                                      % MMP, 10/08/2026
end                                                                         % MMP, 10/08/2026
% prog = lpi_eq_sop(prog,P-Pop,'symmetric');                                % MMP, 10/08/2026 (was)
[prog,tag] = lpi_eq_sop(prog,P-Pop,'symmetric');                            % MMP, 10/08/2026
info = struct('deg',{deg},'terms',terms,'degrees',dinfo,'include',{include},'direct',pinfo, ...
              'tag',tag,'margin',margin,'eD',eD);                           % MMP, 10/08/2026

end
