% Test_tensopmat_native_cost.m -- measures the cost of the operations on a
% 'tensopmat' with a block of type 'dpvar' acting on R, as a function of
% the number of spatial variables K=1,...,4 and of the number q of
% decision variables (1e3, 1e5, 1e6), with a fixed number of nonzero
% coefficients. The coefficients of a 'dpvar' have a row for every
% decision variable, so the cost of the operations on the blocks must not
% grow with q (but with the number of nonzero coefficients, and with the
% number of monomials, 2^K here). The times are printed as a table, and
% compared between q=1e3 and q=1e6 for each K: the operations below,
% except the sum, must not be more than a factor 20 slower (plus 0.2 s).
% The sum of two 'dpvar' merges their lists of decision variables, which
% is linear in q, and is only reported. The peak memory of the session is
% reported and must stay below 4 GB.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - The peak-memory check uses System.Diagnostics (.NET, Windows); where
%   that is unavailable the check is skipped (pws=0).

Ks = 1:4;
qs = [1e3,1e5,1e6];
opnames = {'construct','isvalid','clean_op','uminus','scalar*C','D*C','C*K','[C,C]','C(1,:)','apply_op','C+C'};
nops = numel(opnames);
nflat = 10;                                 % the first nflat operations must not depend on q
T = zeros(numel(Ks),numel(qs),nops);
nz = 40;                                    % nonzero coefficients, independent of q and K
for iq=1:numel(qs)
    q = qs(iq);
    dvarnames = cellstr(compose('d%d',(1:q)'));
    for iK=1:numel(Ks)
        K = Ks(iK);
        names = arrayfun(@(k) sprintf('x%d',k),1:K,'UniformOutput',false);
        vars = struct('in',{cell(1,0)},'out',{names});
        dm = struct('in',zeros(0,2),'out',[zeros(K,1),ones(K,1)]);
        nmon = 2^K;
        degmat = dec2bin(0:nmon-1)-'0';                         % all monomials of degree <=1 in each variable
        degmat = reshape(degmat,nmon,K);
        % coefficients: rows (matrix row outer, constant term first, then the decision variables), columns
        % (matrix column outer, monomial inner), nz of them spread over all q decision variables
        rows = round(linspace(1,(q+1)*2,nz))';
        cols = mod((1:nz)'*7,nmon*2)+1;
        Cq = sparse(rows,cols,1+(1:nz)'/nz,(q+1)*2,nmon*2);
        Dq = dpvar(Cq,degmat,names,dvarnames,[2,2]);
        Sc = polynomial(1,[1,zeros(1,K-1)],names,[1,1]) + 2;    % 2+x1
        Dm = [1,2];     Kc = tensopmat({[1,2;3,4]});    X = {[1;2]};

        t0 = tic;   C = tensopmat({Dq},vars,dm);                T(iK,iq,1) = toc(t0);
        t0 = tic;   lv = isvalid(C);                            T(iK,iq,2) = toc(t0);
        t0 = tic;   Cc = clean_op(C,1e-12);                     T(iK,iq,3) = toc(t0);
        t0 = tic;   Cm = -C;                                    T(iK,iq,4) = toc(t0);
        t0 = tic;   Cs = Sc*C;                                  T(iK,iq,5) = toc(t0);
        t0 = tic;   Cd = Dm*C;                                  T(iK,iq,6) = toc(t0);
        t0 = tic;   Ck = C*Kc;                                  T(iK,iq,7) = toc(t0);
        t0 = tic;   Ch = [C,C];                                 T(iK,iq,8) = toc(t0);
        t0 = tic;   Cr = C(1,:);                                T(iK,iq,9) = toc(t0);
        t0 = tic;   Y = apply_op(C,X);                          T(iK,iq,10) = toc(t0);
        t0 = tic;   Cp = C+C;                                   T(iK,iq,11) = toc(t0);
        if ~lv || ~isa(Cs.ops{1},'dpvar') || ~isequal(size(Cd),[1,2]) || ~isequal(size(Ck),[2,2]) ...
                || ~isequal(size(Ch),[2,4]) || ~isequal(size(Cr),[1,2]) || ~isequal(size(Y),[2,1]) ...
                || ~isequal(size(Cp),[2,2]) || ~isequal(size(Cc),[2,2]) || ~isequal(size(Cm),[2,2])
            error("Test_tensopmat_native_cost failed: object for K=%d, q=%d.",K,q)
        end
    end
    clear dvarnames Dq Cq C Cc Cm Cs Cd Ck Ch Cr Cp Y
end

for iK=1:numel(Ks)
    fprintf('K=%d: times [%s] for q=%s\n',Ks(iK),strjoin(opnames,', '),mat2str(qs))
    for iq=1:numel(qs)
        fprintf('     q=%-8g %s\n',qs(iq),sprintf('%9.4f',squeeze(T(iK,iq,:))))
    end
end
try
    pws = System.Diagnostics.Process.GetCurrentProcess().PeakWorkingSet64/2^30;
    fprintf('Peak working set of the session: %.2f GB\n',pws)
catch
    pws = 0;
end
if pws>4
    error("Test_tensopmat_native_cost failed: peak memory %.2f GB.",pws)
end
bad = T(:,end,1:nflat) > 20*T(:,1,1:nflat)+0.2;
if any(bad(:))
    [iK,~,k] = ind2sub(size(bad),find(bad,1));
    error("Test_tensopmat_native_cost failed: the cost of '%s' for K=%d grows with the number of decision variables (%.3f s for q=%g, %.3f s for q=%g).", ...
        opnames{k},Ks(iK),T(iK,end,k),qs(end),T(iK,1,k),qs(1))
end
fprintf('Test_tensopmat_native_cost passed.\n')
