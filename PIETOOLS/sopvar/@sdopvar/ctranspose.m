function At = ctranspose(A)
% Pt = ctranspose(P) transposes a sdopvar operator P: L_2^p[S1,S3] to L_2^q[S2,S3]
% Date: 8/17/26
%
% MMP, 08/29/2026: Preserve the canonical multiplier form. Swapping the bases
%                  in a multiplier direction moves the monomial degrees onto
%                  the dummy variable, which represents the same operator but
%                  leaves the coefficients no longer determined by it, so
%                  P==P' failed for self-adjoint P and constraining the
%                  coefficients of P-P' to zero was stronger than
%                  constraining the operator. The adjoint is now assembled
%                  through 'canonical_adjoint_map', which swaps only the
%                  integral directions.
% Compute the number of variables S3 that appear in both the input and output spaces.
nvars3 = numel(intersect(A.vars.in,A.vars.out));
% each common variable has three kernel types: multiplier, lower integral, upper integral.
%1 = multiplier
%2 = lower integral
%3 = upper integral
%Under adjoint/transposition, lower and upper integral kernels swap. Multiplier terms stay multiplier
idx = repmat({[1 3 2]},1,nvars3);
params = A.params;
if numel(params.A)~=3^nvars3                                                % MMP, 08/29/2026
    error("An 'sdopvar' object should hold 3^n3 parameters, one per "...    % MMP, 08/29/2026
          +"multi-index in {1,2,3}^n3.")                                    % MMP, 08/29/2026
end                                                                         % MMP, 08/29/2026
if nvars3>0                                                                 % MMP, 08/29/2026
    lin = reshape(1:3^nvars3,[3*ones(1,nvars3),1]);                         % MMP, 08/29/2026
    adjidx = lin(idx{:});    % where the adjoint of parameter k belongs     % MMP, 08/29/2026
else                                                                        % MMP, 08/29/2026
    adjidx = 1;                                                             % MMP, 08/29/2026
end                                                                         % MMP, 08/29/2026

% The coefficients are rearranged, not merely transposed. In a direction     % MMP, 08/29/2026
% carrying an integral the primary and dummy variables exchange roles, so    % MMP, 08/29/2026
% the monomial bases swap. In a direction carrying a multiplier they do      % MMP, 08/29/2026
% not: delta(s-s') collapses the kernel onto the diagonal, so the adjoint    % MMP, 08/29/2026
% is the same polynomial in s with the matrix indices transposed, and the    % MMP, 08/29/2026
% degree stays on the left basis. That is what keeps the canonical           % MMP, 08/29/2026
% multiplier form invariant under taking adjoints.                           % MMP, 08/29/2026
[Tcell,ZLt,ZRt] = canonical_adjoint_map(A.vars,A.ZL,A.ZR,A.dims);           % MMP, 08/29/2026
nC_t = A.dims(2)*prod([cellfun(@numel,ZLt),1]) ...                          % MMP, 08/29/2026
       * A.dims(1)*prod([cellfun(@numel,ZRt),1]);                           % MMP, 08/29/2026
q = size_decvars(A);                                                        % MMP, 08/29/2026

pA = cell(size(params.A));      pB = cell(size(params.B));                  % MMP, 08/29/2026
for ii=1:numel(params.A)                                                    % MMP, 08/29/2026
    Aii = params.A{ii};         Bii = params.B{ii};                         % MMP, 08/29/2026
    if isempty(Aii)                                                         % MMP, 08/29/2026
        pA{adjidx(ii)} = sparse(nC_t,1);                                    % MMP, 08/29/2026
    else                                                                    % MMP, 08/29/2026
        pA{adjidx(ii)} = Tcell{ii}*sparse(Aii(:));                          % MMP, 08/29/2026
    end                                                                     % MMP, 08/29/2026
    if isempty(Bii)                                                         % MMP, 08/29/2026
        pB{adjidx(ii)} = sparse(q,nC_t);                                    % MMP, 08/29/2026
    else                                                                    % MMP, 08/29/2026
        pB{adjidx(ii)} = Bii*Tcell{ii}.';                                   % MMP, 08/29/2026
    end                                                                     % MMP, 08/29/2026
end                                                                         % MMP, 08/29/2026
params.A = pA;      params.B = pB;                                          % MMP, 08/29/2026

At_vars = struct('in',{A.vars.out},'out',{A.vars.in});
At_dom = struct('in',A.dom.out,'out',A.dom.in);
At = sdopvar(params,At_vars,A.Zd,ZLt,ZRt,At_dom,[A.dims(2),A.dims(1)]);     % MMP, 08/29/2026
end

function q = size_decvars(A)                                                % MMP, 08/29/2026
% Number of decision variables, taken from a stored parameter where possible % MMP, 08/29/2026
% so that an operator with none is handled the same way.                     % MMP, 08/29/2026

q = numel(A.Zd);                                                            % MMP, 08/29/2026
for ii = 1:numel(A.params.B)                                                % MMP, 08/29/2026
    if ~isempty(A.params.B{ii})                                             % MMP, 08/29/2026
        q = size(A.params.B{ii},1);                                         % MMP, 08/29/2026
        return                                                              % MMP, 08/29/2026
    end                                                                     % MMP, 08/29/2026
end                                                                         % MMP, 08/29/2026
end                                                                         % MMP, 08/29/2026

