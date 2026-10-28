function [prog,Pop,info] = lpi_ineq_sop(prog,P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP,INFO] = LPI_INEQ_SOP(PROG,P,OPTS) adds to the LPI program PROG
% the constraint P >= 0 for a self-adjoint PI operator P. For a 'cdopvar'
% or 'copvar' P on X = L_2^{m_1}[s^1] x ... x L_2^{m_M}[s^M] (an empty s^k
% is R^{m_k}), it declares the positive operator
%
%       Pop = sum_c Z_c'*(g_c*Q_c)*Z_c,     Q_c >= 0,
%
% one term c per weight g_c(theta) >= 0 on the domain, Z_c a column of
% multiplier and integral operators from X to L_2[theta] with monomial
% degrees read from the support of P or of OPTS.like (GET_LIFT_DEGS), Q_c
% a Gram matrix of decision variables (POSCOPVAR_DIRECT), and imposes
% P - Pop == 0 (LPI_EQ_SOP, 'symmetric'). Before declaring Pop, it prints
% the number of decision variables Pop adds and the program total. For P
% of any other class it calls LPI_INEQ(PROG,P,OPTS).
%
% INPUT
% - prog:  'struct', the LPI program (LPIPROGRAM_SOP or LPIPROGRAM); read
%          are decvartable and, if present, sop_pos (cell of the 'cdopvar'
%          objects POSCOPVAR declared) and sopeq (cell of the equality
%          records read by LPIGETDUAL_SOP);
% - P:     'cdopvar' or 'copvar' object, square (the same output and input
%          spaces) and self-adjoint (not checked), the operator to make
%          positive; for a 'copvar' P the constraint certifies P >= 0. Any
%          other class ('dopvar', 'dopvar2d', 'opvar', 'opvar2d', 'dpvar',
%          'polynomial') is passed to LPI_INEQ;
% - opts:  (optional) 'struct' with the optional fields
%            dD, dw        double, added to the lift degree and to the
%                          weight degree read from the support (default
%                          1, 1);
%            wR            double, the weight degree on an R^q space
%                          (default 0);
%            psatz         the terms of Pop: 'auto' (default; in one
%                          variable the plain term and the box product one
%                          weight degree lower, in nv ~= 1 variables the
%                          plain term and the 2nv faces of the box, one
%                          degree lower when every weight degree is at
%                          least 2), 'none' (the plain term), 'product' (as
%                          'auto' in one variable), 'faces' (the plain term
%                          and the 2nv faces), or a double row of term
%                          codes: 0 plain, 1 box product, 2d+1 / 2d+2 lower
%                          / upper face of direction d, 2nv+2+d box
%                          quadratic of direction d;
%            psatz_offset  double, with a row of codes: the weight degree
%                          by which each term is lowered, one per code or
%                          one for all (default 0);
%            like          'cdopvar' or 'copvar' object, or a cell of them,
%                          with variables among those of P: their pooled
%                          support is read instead of that of P. Default
%                          for a 'cdopvar' P without deg: the operators of
%                          prog.sop_pos with variables among those of P
%                          whose decision variables P carries (tested on
%                          up to five of them); if there is none, P;
%            deg           a degree specification used instead of the one
%                          read from the support: a scalar, a 'struct' with
%                          fields int, mult, joint, subset, or a 1 x M cell
%                          of them, one per space; the terms are then psatz
%                          and psatz_offset for a numeric psatz, else those
%                          that psatz names;
%            prune         logical (default true): on each space with
%                          spatial variables, keep only the basis operators
%                          whose diagonal term reaches a nonzero parameter
%                          cell of the diagonal block of P;
%            path          char, 'fast' (default) or 'general', the route
%                          of POSCOPVAR_DIRECT;
%            cache         logical (default true): keep the position maps
%                          of POSCOPVAR_DIRECT for the session;
%            margin        double >= 0 (default 0): impose
%                          P + margin*eD - Pop == 0 instead, eD = A'*A with
%                          A the lift of LIFT_COPVAR at the largest 'mult'
%                          degree of the specification (at least 1);
%          for P of another class, OPTS is passed to LPI_INEQ (fields
%          psatz, pure, sep).
%
% OUTPUT
% - prog:  'struct', the program with the Gram matrices of Pop declared
%          positive semidefinite and P - Pop == 0 imposed; for P of another
%          class, the program of LPI_INEQ;
% - Pop:   'cdopvar' object on the spaces of P, the positive operator Pop
%          ([] for P of another class);
% - info:  'struct' with fields
%            deg         the degree specification used: a 1 x M cell of
%                        'struct' objects read from the support, or
%                        opts.deg;
%            terms       'struct' with fields codes and offsets (double
%                        rows), the terms of Pop;
%            degrees     'struct' of the degrees read from the support,
%                        with fields Dmin, Dmax, Mdeg, J, D, w, wR, dD, dw,
%                        nv, M, vars, isR (GET_LIFT_DEGS), or
%                        struct('given',true) with opts.deg;
%            include     1 x M cell, entry k the multi-index rows of the
%                        basis operators kept on space k ([] when not
%                        pruned);
%            direct      'struct', the INFO output of POSCOPVAR_DIRECT;
%            tag         double, the index in prog.sopeq of the record of
%                        P - Pop == 0, for LPIGETDUAL_SOP (0 when PROG has
%                        no field sopeq);
%            margin      double, the margin used;
%            eD          'copvar' object, eD ([] for margin 0);
%            N           double row, the Gram dimension of each term;
%            ndec        double, sum(N.^2), the entries Pop adds to
%                        prog.decvartable (Pop has N.*(N+1)/2 distinct
%                        decision variables per term);
%            sized_from  char, the source of the degrees: 'the operator',
%                        'opts.like', 'opts.deg', or, for the default
%                        like, the number of operators found followed by
%                        ' operator(s) of poscopvar';
%          for P of another class, struct('legacy',true).
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
% MMP, 10/21/2026: Maintainer: lpi_ineq_sop always estimates the decision
%                variables needed, and needs no hand-set options. (1) The
%                decision variables of the certificate are counted before
%                anything is declared (POSCOPVAR_DIRECT 'estimate') and
%                printed with the program total; INFO.ndec, INFO.N. (2)
%                Without opts.like and opts.deg, the support is read from
%                the positive operators of 'poscopvar' that P depends on
%                (prog.sop_pos), as the executives pass their storage: sized
%                from P itself, the 2-D stability LPI of examples/
%                stability_2d_sop declared 12036880 decision variables and
%                was unsolved after 900 s and 31 GB, against 197920 from
%                its operator (README Sec. 4).
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body, INPUT and
%                OUTPUT lines for the local function; code unchanged. The
%                deleted header text was the comparison with lpi_ineq and
%                its dispatch table, measurement results, the rationale
%                for one call (NOTES), dates and README references. The
%                help text of the 10/21/2026 entry is in the paragraph and
%                in the entries of opts.like and info, without markers.

% Check the inputs; OPTS is optional
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
% OPTS as a struct; P square
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

% % % Without opts.like and opts.deg: the support of the positive           % MMP, 10/21/2026
% operators that 'poscopvar' recorded in prog.sop_pos and P depends on, as  % MMP, 10/21/2026
% the executives pass their storage. Read from P itself, a P linear in a    % MMP, 10/21/2026
% decision operator gets a certificate for every value of that operator     % MMP, 10/21/2026
% (GET_LIFT_DEGS NOTES). P depends on L when it carries L's decision        % MMP, 10/21/2026
% names: up to five sampled names, each one O(q) strcmp over P's list.      % MMP, 10/21/2026
src = 'the operator';                                                       % MMP, 10/21/2026
hasopt = @(f) isfield(opts,f) && ~isempty(opts.(f));                        % MMP, 10/21/2026
if hasopt('like'),  src = 'opts.like';  end                                 % MMP, 10/21/2026
if hasopt('deg'),   src = 'opts.deg';   end                                 % MMP, 10/21/2026
if ~hasopt('like') && ~hasopt('deg') && isfield(prog,'sop_pos') && isa(P,'cdopvar') % MMP, 10/21/2026
    zP = dvar_names(P);     like = {};                                      % MMP, 10/21/2026
    for i = 1:numel(prog.sop_pos)                                           % MMP, 10/21/2026
        L = prog.sop_pos{i};    zL = dvar_names(L);                         % MMP, 10/21/2026
        if isempty(zL) || ~all(ismember(L.vars,P.vars)),   continue,   end  % MMP, 10/21/2026
        smp = zL(unique(round(linspace(1,numel(zL),min(5,numel(zL))))));    % MMP, 10/21/2026
        if any(cellfun(@(z) any(strcmp(zP,z)),smp)),  like{end+1} = L;  end %#ok<AGROW> % MMP, 10/21/2026
    end                                                                     % MMP, 10/21/2026
    if ~isempty(like)                                                       % MMP, 10/21/2026
        opts.like = like;                                                   % MMP, 10/21/2026
        src = sprintf('%d operator(s) of poscopvar',numel(like));           % MMP, 10/21/2026
    end                                                                     % MMP, 10/21/2026
end                                                                         % MMP, 10/21/2026

% % % Degrees and terms: the reader, or the given specification.
if isfield(opts,'deg') && ~isempty(opts.deg)
    deg = opts.deg;     dinfo = struct('given',true);
    nv = numel(P.vars);
    % Terms: a numeric opts.psatz and its offsets, else those of the reader
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
    % Degrees and terms read from the support
    [deg,terms,dinfo] = get_lift_degs(P,opts);
end

% % % The basis operators each L2 space needs: those whose diagonal term can
% reach the support of P's diagonal block.
prune = true;
if isfield(opts,'prune') && ~isempty(opts.prune),   prune = logical(opts.prune);  end
include = cell(1,M);
if prune
    for k = 1:M
        % Spaces with variables and a stored diagonal block only
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
% path and cache, when given
for f = {'path','cache'}
    if isfield(opts,f{1}) && ~isempty(opts.(f{1})),  popts.(f{1}) = opts.(f{1});    end
end
% % % The decision variables of the certificate, counted before anything    % MMP, 10/21/2026
% is declared (POSCOPVAR_DIRECT 'estimate': the grids of the build below,   % MMP, 10/21/2026
% no Gram, no map) and printed, so an oversized program shows before the    % MMP, 10/21/2026
% build and the solve.                                                      % MMP, 10/21/2026
eopts = popts;  eopts.estimate = true;                                      % MMP, 10/21/2026
[~,~,~,einfo] = poscopvar_direct(prog,dm,sp,dom,deg,eopts);                 % MMP, 10/21/2026
% Counted as prog.decvartable counts them: N^2 per Gram block of dimension  % MMP, 10/21/2026
% N (the N(N+1)/2 of einfo.ndec are the distinct names of Pop).             % MMP, 10/21/2026
nd = sum(einfo.N.^2);                                                       % MMP, 10/21/2026
fprintf(' - lpi_ineq_sop: %d decision variables in the certificate (Gram blocks %s, sized from %s); %d in the program.\n',...
        nd,mat2str(einfo.N),src,numel(prog.decvartable)+nd);                % MMP, 10/21/2026
% Declare Pop, every term in one call
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
% Impose P - Pop == 0, one block per adjoint pair; TAG indexes its record
% prog = lpi_eq_sop(prog,P-Pop,'symmetric');                                % MMP, 10/08/2026 (was)
[prog,tag] = lpi_eq_sop(prog,P-Pop,'symmetric');                            % MMP, 10/08/2026
% The info output
% info = struct('deg',{deg},'terms',terms,'degrees',dinfo,'include',{include},'direct',pinfo, ...
%               'tag',tag,'margin',margin,'eD',eD);                         % MMP, 10/21/2026 (was)
info = struct('deg',{deg},'terms',terms,'degrees',dinfo,'include',{include},'direct',pinfo, ...
              'tag',tag,'margin',margin,'eD',eD,'ndec',nd,'N',einfo.N,'sized_from',src); % MMP, 10/21/2026

end


function z = dvar_names(X)                                                  % MMP, 10/21/2026
% The decision names of a 'cdopvar': every 'sdopvar' block holds the one    % MMP, 10/21/2026
% list of the container (cdopvar unify_dvars); {} for a fixed operator.     % MMP, 10/21/2026
% INPUT  X: 'cdopvar' or 'copvar' object.
% OUTPUT z: cellstr, the decision variable names ({} for a 'copvar').
z = {};                                                                     % MMP, 10/21/2026
for i = 1:numel(X.C)                                                        % MMP, 10/21/2026
    if isa(X.C{i},'sdopvar'),   z = X.C{i}.Zd;  return,     end             % MMP, 10/21/2026
end                                                                         % MMP, 10/21/2026
end                                                                         % MMP, 10/21/2026
