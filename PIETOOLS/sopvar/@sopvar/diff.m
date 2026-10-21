function [dPop,is_zero] = diff(Pop,vars,d,chk_override)
% [DPOP,IS_ZERO] = DIFF(POP,VARS,D) takes the composition of the PI
% operator POP with a Dth order differential operator in variables VARS
%   DPOP = d^D/(d vars^D) o POP
% VARS and D can be row vectors to simultaneously differentiate with
% respect to multiple variables. D can also be a column vector or matrix to
% take multiple derivatives, concatenated vertically.
%
% INPUTS
% - Pop:    m x n 'sopvar' object representing a PI operator;
% - vars:   1 x r 'cellstr' or 'polynomial' (pvar) array specifying the
%           names of the spatial variables with respect to which to take
%           the derivative;
% - d:      (optional) p x r array of nonnegative integers specifying the
%           order of the derivatives to be taken with respect to each of
%           the spatial variables. Defaults to 1;
% - chk_override:   (optional) boolean value. Set to true to avoid
%                   checking whether Pop, vars and d are appropriately
%                   specified, in which case vars must be a 1 x r cellstr
%                   of unique names, and d a p x r array of nonnegative
%                   integers;
%
% OUTPUTS
% - dPop:   m*p x n 'sopvar' representing the composition of the
%           requested differential operator with the PI operator Pop;
% - is_zero:    boolean value set to true if the output operator is known
%               to be zero (false does not imply the operator is not zero);
%
% NOTES
% - The derivative is only computed if the composition is again a PI
%   operator. For example, an error will be thrown when differentiating a
%   pure multiplier operator, (Pop*u)(s) = P(s)*u(s), since
%       d/ds P(s)*u(s) = (d/ds P)(s)*u(s) + P(s)*(d/ds u)(s)
%   which cannot be expressed as dPop*u for a PI operator dPop;
% - The derivative is computed using the Leibniz integral rule,
%   d/ds int_{a(s)}^{b(s)} P(s,t) u(t) dt
%       = (d/ds b)(s)*P(s,b(s))*u(b(s)) - (d/ds a)(s)*P(s,a(s))*u(a(s))
%           + int_{a(s)}^{b(s)} (d/ds P)(s,t) u(t) dt
%   The boundary terms are multipliers, so a repeated derivative in the
%   same variable is only supported if they cancel, K_2(s,s)=K_3(s,s) (up
%   to roundoff, relative 1e-12);
% - A variable that is not an output variable of Pop does not affect
%   Pop*u, so differentiating with respect to it returns the zero operator;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - diff
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
% DJ, 10/20/2026: Initial coding;

if nargin<=3 || ~chk_override
    % Check that the operator is properly specified
    if ~isa(Pop,'sopvar')
        error("PI operator must be specified as object of type 'sopvar'.")
    end

    % Check that the variables and orders of differentiation are properly
    % specified
    if nargin<=2
        [vars,d] = parse_diff_inputs(vars);
    else
        [vars,d] = parse_diff_inputs(vars,d);
    end
end
p = size(d,1);

% Return quickly for zero derivative or zero outputs
is_zero = false;
if isempty(vars)
    % Differentiation with no variables returns the original operator
    dPop = repmat(Pop,p,1);
    return
end

% Each row of d is a separate derivative, applied one variable and one
% order at a time (the derivatives commute); rows are stacked vertically.
is_out = ismember(vars,Pop.vars.out);
rows = cell(p,1);
row_zero = false(p,1);
for i=1:p
    di = d(i,:);
    if any(di(~is_out)>0)
        % Pop*u depends only on the output variables, so d/dv = 0 for
        % any other variable v.
        rows{i} = 0*Pop;
        row_zero(i) = true;
        continue
    end
    rows{i} = Pop;
    % Consecutively apply derivative with respect to each variable
    for k=find(di>0)
        for j=1:di(k)
            [rows{i},row_zero(i)] = diff_step(rows{i},vars{k});
            if row_zero(i)
                break
            end
        end
        if row_zero(i)
            break
        end
    end
end
if p==1
    dPop = rows{1};
else
    dPop = vertcat(rows{:});
end
is_zero = all(row_zero);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function [dP,is_zero] = diff_step(P,v)
% [DP,IS_ZERO] = DIFF_STEP(P,V) returns d/dV o P for the single output
% variable V (char), by the Leibniz rule applied to the kernel form of
% sopvar, 
%   y = sum_gam int K_gam(s,s') I_gam(s-s') u(s') ds', 
% with
%   K_gam(s,s') = (I_m o ZL(s))' unvec(C_gam) (I_n o ZR(s')).
% Cell gam has gam_k = 1 (multiplier), 2 (int_a^s) or 3 (int_s^b) in the
% k-th shared variable. For V = s_k, per cell:
%   gam_k = 1:      d/ds (M u) contains u_s, so M must be zero (error);
%   gam_k = 2,3:    d/ds int_{a|s}^{s|b} K u dt = +-K(s,s) u(s)
%                   + int_{a|s}^{s|b} (d/ds K) u dt, the first term a
%                   multiplier of cell gam with gam_k=1, K_2(s,s)-K_3(s,s),
%   V an output 
%   variable (S2):  no integral, d/ds K only.
%
% d/ds K lowers the degree a of s in ZL by one with factor a. K(s,s) moves
% the degree b of the dummy variable onto ZL, (a,b) -> (a+b,0), which is
% the canonical multiplier form; ZR is enlarged by degree 0 if it lacks it.

% Extract variables and monomials of operator P 
vin = P.vars.in(:)';    vout = P.vars.out(:)';
M = numel(vout);        N = numel(vin);
m = P.dims(1);          n = P.dims(2);
ZL = P.ZL;              ZR = P.ZR;
% Determine which input/output variable in P corresponds to var
vl_idx = find(strcmp(vout,v),1);        % position of v in ZL
vr_idx = find(strcmp(vin,v),1);         % position of v in ZR (empty for S2)
vS3 = intersect(vin,vout);          n3 = numel(vS3);
v3_idx = find(strcmp(vS3,v),1);         % shared direction of v (empty for S2)
nL = cellfun(@numel,ZL(:)');        nR = cellfun(@numel,ZR(:)');
NL = prod([nL,1]);                  NR = prod([nR,1]);
sL = kron_strides(nL);              sR = kron_strides(nR);
ncell = 3^n3;

% % Build the parameters and monomials of the derivative manually:
% % For each parameter in P, determine to which parameters in dP it
% % contributes, determining the row indices, column indices, and left- and
% % right-monomials for each coefficient
param_idcs = cell(2*ncell,1);                   % parameter in dP
ridcs = param_idcs;        cidcs = param_idcs;  % row and column indices  
degsL = param_idcs;        degsR = param_idcs;  % left and right monomials
vals = param_idcs;                              % coefficient values
% Loop over all the parameters
cntr = 0;
need_zero_ZR = false;   % keep track of whether we need to add a zero monomial to ZR
for k=1:ncell
    gam = gamma_of_cell(k,n3);
    % Do not support differentiation of multipliers
    if ~isempty(v3_idx) && gam(v3_idx)==1
        if nnz(P.params{k})>0
            error("Desired composition with derivative cannot be expressed as a PI operator");
        end
        continue
    end
    % Take the derivative of the kernel within the integral
    [ridcs_k,cidcs_k,vals_k] = find(P.params{k});
    ridcs_k = ridcs_k(:);     
    cidcs_k = cidcs_k(:);     
    vals_k = vals_k(:);
    if ~isempty(vals_k)
        % Determine left- and right-monomial degrees associated with 
        % coefficients (ridcs_k,cidcs_k)
        [i1,dl,j1,dr] = term_degrees(ridcs_k,cidcs_k,ZL,ZR,NL,NR,sL,sR);
        % Retain only terms which depend on the differentiation variables
        rtn_idcs = dl(:,vl_idx)>0;
        cntr = cntr+1;
        % Derivative of kernel contributes to the same parameter
        param_idcs{cntr} = k*ones(nnz(rtn_idcs),1);
        ridcs{cntr} = i1(rtn_idcs);          cidcs{cntr} = j1(rtn_idcs);
        % Coefficients get scaled by monomial degree when differentiating
        vals{cntr} = vals_k(rtn_idcs).*dl(rtn_idcs,vl_idx);
        % Monomial degree gets reduced by 1
        dl(:,vl_idx) = dl(:,vl_idx)-1;
        degsR{cntr} = dr(rtn_idcs,:);        degsL{cntr} = dl(rtn_idcs,:);
    end
    % Compute the boundary terms K_2(s,s)-K_3(s,s)
    % --> becomes a multiplier in direction t
    if ~isempty(v3_idx) && gam(v3_idx)==2
        % Determine index of parameter for upper integral (k3) and
        % multiplier (k1) associated with parameter k
        k3 = k+3^(v3_idx-1);             k1 = k-3^(v3_idx-1);
        % Take the difference B(s,t) = K_2(s,t) - K_3(s,t);
        Pk = P.params{k};           Pk3 = P.params{k3};
        if nnz(Pk)==0 && nnz(Pk3)==0
            B = sparse(m*NL,n*NR);
        elseif nnz(Pk)==0
            B = -Pk3;
        elseif nnz(Pk3)==0
            B = Pk;
        else
            B = Pk - Pk3;
        end
        % Evaluate B(s,t) at t=s
        [ridcs_k,cidcs_k,vals_k] = find(B);
        ridcs_k = ridcs_k(:);     cidcs_k = cidcs_k(:);     vals_k = vals_k(:);
        if ~isempty(vals_k)
            % Determine left- and right-monomial degrees associated with 
            % coefficients (ridcs_k,cidcs_k)
            [i1,dl,j1,dr] = term_degrees(ridcs_k,cidcs_k,ZL,ZR,NL,NR,sL,sR);
            cntr = cntr+1;
            param_idcs{cntr} = k1*ones(numel(vals_k),1);
            ridcs{cntr} = i1;            cidcs{cntr} = j1;            
            vals{cntr} = vals_k;
            dl(:,vl_idx) = dl(:,vl_idx)+dr(:,vr_idx);         % B(s,s): s^a s'^b -> s^(a+b)
            dr(:,vr_idx) = 0;                       % B(s,s) does not depend on t
            degsL{cntr} = dl;            degsR{cntr} = dr;
            need_zero_ZR = true;
        end
    end
end
% cntr==0 implies the output does not depend on one of the variables
if cntr==0
    dP = 0*P;
    is_zero = true;
    return
end
% Concatenate
is_zero = false;
param_idcs = vertcat(param_idcs{1:cntr});     
ridcs = vertcat(ridcs{1:cntr});     degsL = vertcat(degsL{1:cntr});
cidcs = vertcat(cidcs{1:cntr});     degsR = vertcat(degsR{1:cntr});     
vals = vertcat(vals{1:cntr});

% % Build the new monomial bases
ZLn = ZL;       ZRn = ZR;
ZLn{vl_idx} = unique(degsL(:,vl_idx));
if need_zero_ZR && ~any(ZR{vr_idx}==0)
    % Make sure ZR has constant monomial
    ZRn{vr_idx} = union(ZR{vr_idx}(:),0);
end
nLn = cellfun(@numel,ZLn(:)');      nRn = cellfun(@numel,ZRn(:)'); 
NLn = prod([nLn,1]);                NRn = prod([nRn,1]);    
sLn = kron_strides(nLn);            sRn = kron_strides(nRn);

% % Build the new parameters
rowL = zeros(numel(vals),1);      colR = zeros(numel(vals),1);
for q=1:M
    rowL = rowL + sLn(q)*monomial_position(monomial_position_map(ZLn{q}),degsL(:,q));
end
for q=1:N
    colR = colR + sRn(q)*monomial_position(monomial_position_map(ZRn{q}),degsR(:,q));
end
rowL = ridcs*NLn + rowL + 1;       colR = cidcs*NRn + colR + 1;
params = cell(size(P.params));
for k=1:ncell
    idx = (param_idcs==k);
    params{k} = sparse(rowL(idx),colR(idx),vals(idx),m*NLn,n*NRn);
    if ~isempty(v3_idx) && any(idx)
        % The multiplier in the differentiation variable only contains the
        % boundary terms K_2(s,s)-K_3(s,s), which cancel for kernels that
        % are continuous at t=s. Remove coefficients below tolerance
        gam = gamma_of_cell(k,n3);
        if gam(v3_idx)==1
            [ri,ci,vi] = find(params{k});
            keep = abs(vi)>1e-12*max(abs(vals(idx)));
            params{k} = sparse(ri(keep),ci(keep),vi(keep),m*NLn,n*NRn);
        end
    end
end

% % Build the new operator
dP = sopvar(params,P.vars,ZLn,ZRn,P.dom,P.dims);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function [i1,dl,j1,dr] = term_degrees(ri,ci,ZL,ZR,NL,NR,sL,sR)
% Matrix row/column and monomial degrees of the entries (ri,ci) of a
% parameter; rows are ordered (matrix row, left monomial), columns (matrix
% column, right monomial), monomials as kron(Z{1},...,Z{end}).
i1 = floor((ri-1)/NL);      a = kron_split(mod(ri-1,NL),sL);
j1 = floor((ci-1)/NR);      b = kron_split(mod(ci-1,NR),sR);
dl = zeros(numel(ri),numel(ZL));
for q=1:numel(ZL)
    dq = ZL{q}(:);      dl(:,q) = dq(a(:,q)+1);
end
dr = zeros(numel(ci),numel(ZR));
for q=1:numel(ZR)
    dq = ZR{q}(:);      dr(:,q) = dq(b(:,q)+1);
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

