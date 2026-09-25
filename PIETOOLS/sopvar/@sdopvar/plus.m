function C = plus(A,B)
% Add two sdopvar operators, preserving affine decision representation.
%
% MMP, 09/07/2026: Accept a fixed 'sopvar' summand, promoting it with
% 'sopvar2sdopvar'. Adding a known block to a decision operator is normal
% usage, and since 'sdopvar' now lists ?sopvar among its InferiorClasses
% both A+Pop and Pop+A are dispatched here.
% MMP, 09/21/2026: Compare the variable lists with 'isequal' rather than
% 'any(~strcmp(...))'. '@sopvar/plus' took the same fix on 09/07/2026 and
% this copy was missed. Reached by adding the strict-positivity identity to a
% 'poscopvar' container, whose R^n block has an empty variable list: '{}' and
% a 1 x 0 cellstr name the same space, and strcmp of two different-sized
% cellstr does not compare them.
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
% mdopvar -> cdopvar, with every file and function named after them. Mechanical
% rename, no functional change. Renamed here: posmopvar -> poscopvar.

% A fixed 'sopvar' operand is promoted to a decision operator with a        % MMP, 09/07/2026
% zero B, so that A+Pop and Pop+A work; mixing fixed and decision blocks    % MMP, 09/07/2026
% is normal usage. The decision variable list is taken from the sdopvar     % MMP, 09/07/2026
% operand so 'CombineDecisionBasis' below takes its fast path. Placed        % MMP, 09/07/2026
% before the checks because an 'sopvar' has no 'params.A'.                  % MMP, 09/07/2026
if isa(A,'sopvar') || isa(B,'sopvar')                                       % MMP, 09/07/2026
    if      isa(B,'sdopvar'),   Zd_p = B.Zd;                                % MMP, 09/07/2026
    elseif  isa(A,'sdopvar'),   Zd_p = A.Zd;                                % MMP, 09/07/2026
    else,                       Zd_p = cell(0,1);                           % MMP, 09/07/2026
    end                                                                     % MMP, 09/07/2026
    if isa(A,'sopvar'),     A = sopvar2sdopvar(A,Zd_p);     end             % MMP, 09/07/2026
    if isa(B,'sopvar'),     B = sopvar2sdopvar(B,Zd_p);     end             % MMP, 09/07/2026
end                                                                         % MMP, 09/07/2026

% Error handling: Checks to ensure A and B are compatible
if any(A.dims~=B.dims)
    error('Dimensions of summands A and B do not match');
end
% 'isequal' rather than 'any(~strcmp(...))', the same fix '@sopvar/plus'     % MMP, 09/21/2026
% took on 09/07/2026 and this line was missed by: strcmp of two cellstr of   % MMP, 09/21/2026
% different size does not compare them, so an operand on {} against one on   % MMP, 09/21/2026
% {'s1'} either throws a size error or silently passes the check. The empty   % MMP, 09/21/2026
% case is reached as soon as a container holds an R^n block, since '{}' and   % MMP, 09/21/2026
% a 1 x 0 cellstr name the same space but are not the same size.             % MMP, 09/21/2026
if ~isequal(A.vars.in(:),B.vars.in(:)) || ~isequal(A.vars.out(:),B.vars.out(:)) % MMP, 09/21/2026
%if any(~strcmp(A.vars.in,B.vars.in)) || any(~strcmp(A.vars.out,B.vars.out)) % MMP, 09/21/2026 (was)
    error('Summands A and B map between different spaces');
end
if any(any(A.dom.in~=B.dom.in)) || any(any(A.dom.out~=B.dom.out))
    error('Input or output variables in summands A and B have different domains');
end
if numel(A.params.A)~=numel(B.params.A)
    error('number of terms in summands is not equal -- one of them is probably malformed');
end
% Put A and B on a common decision variable list and common monomial        % MMP, 09/07/2026
% bases. 'sync_basis' is N-ary and replaces the prologue that plus, eq,     % MMP, 09/07/2026
% horzcat and vertcat each carried inline; T{k} is empty when operand k     % MMP, 09/07/2026
% needs no remapping, and 'apply_basis_map' then skips the multiply.        % MMP, 09/07/2026
[ops,T,Zd,ZL,ZR] = sync_basis({A,B});                                       % MMP, 09/07/2026
A = ops{1};     B = ops{2};                                                 % MMP, 09/07/2026
%[A,B,Zd] = CombineDecisionBasis(A,B);                                     % MMP, 09/07/2026 (was)
%[ZR,C1R,C2R] = UnionBasisMonomials(A.ZR,B.ZR);                            % MMP, 09/07/2026 (was)
%[ZL,C1L,C2L] = UnionBasisMonomials(A.ZL,B.ZL);                            % MMP, 09/07/2026 (was)
%C1L = kron(eye(A.dims(1)),C1L);                                           % MMP, 09/07/2026 (was)
%C2L = kron(eye(B.dims(1)),C2L);                                           % MMP, 09/07/2026 (was)
%C1R = kron(eye(A.dims(2)),C1R);                                           % MMP, 09/07/2026 (was)
%C2R = kron(eye(B.dims(2)),C2R);                                           % MMP, 09/07/2026 (was)
%T1 = kron(C1R', C1L');                                                    % MMP, 09/07/2026 (was)
%T2 = kron(C2R', C2L');                                                    % MMP, 09/07/2026 (was)
% Once C_1, C_2 have compatible sizes and variables,
% C_1(d) + C_2(d) = unvec(A_1 + A_2 + (Bt_1 + Bt_2)*d)
% with the change of monomial basis,
% C_1L C_1(d) C_1R + C_2L C_2(d) C_2R =unvec(T_1 A_1 +T_2 A_2 +(T_1 Bt_1 +T_2 Bt_2)*d)
params_new.A = cell(size(A.params.A));
params_new.B = cell(size(A.params.B));
for i=1:numel(A.params.A)
   % params.A{ii} = params.A{ii}+pb.A{ii};
  %  params.Bt{ii} = params.Bt{ii}+pb.Bt{ii};
    [A1,B1] = apply_basis_map(T{1},A.params.A{i},A.params.B{i});            % MMP, 09/07/2026
    [A2,B2] = apply_basis_map(T{2},B.params.A{i},B.params.B{i});            % MMP, 09/07/2026
    params_new.A{i} = A1 + A2;                                              % MMP, 09/07/2026
    params_new.B{i} = B1 + B2;                                              % MMP, 09/07/2026
%   params_new.A{i} = T1*A.params.A{i} + T2*B.params.A{i};                  % MMP, 09/07/2026 (was)
%   params_new.B{i} = A.params.B{i}*T1.' + B.params.B{i}*T2.';              % MMP, 09/07/2026 (was)
end
C = sdopvar(params_new,A.vars,Zd,ZL,ZR,A.dom,A.dims);
end
