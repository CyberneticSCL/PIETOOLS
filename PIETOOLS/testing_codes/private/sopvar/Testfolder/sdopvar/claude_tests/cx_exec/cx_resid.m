function [rel_b,psd_min,psd_relmin,trivial] = cx_resid(sol)
% [REL_B,PSD_MIN,PSD_RELMIN,TRIVIAL] = CX_RESID(SOL): solution quality of a
% solved LPI program, for the stock and container paths alike.
%
% Solution quality, for solves the numerr flag does not settle (2-D):
% rel_b = ||At'*x - b|| / ||b|| in decvartable order (expr.At rows and
% solinfo.RRx are both decvartable order), and the PSD margin of the Gram
% blocks read from solinfo.x, which is SeDuMi's CONE order - RRx sliced by
% cone blocks gives a wrong margin (baseline harness trap 1). psd_relmin
% divides by the largest eigenvalue over ALL blocks: a numerically zero
% block's own-max ratio is rounding noise. trivial flags X = 0, which an
% infeasible solve returns and which makes every cone measure look perfect.
%
% Initial coding MMP, 09/25/2026
S = cx_shape(sol);
rel_b = NaN;    psd_min = NaN;  psd_relmin = NaN;   trivial = true;
if ~isfield(sol,'solinfo') || ~isfield(sol.solinfo,'RRx') || isempty(sol.solinfo.RRx)
    return
end
Atf = [];   bf = [];
for i = 1:sol.expr.num
    Atf = [Atf, sol.expr.At{i}];    bf = [bf; sol.expr.b{i}];               %#ok<AGROW>
end
xr = sol.solinfo.RRx(:);
nb = norm(full(bf));
rel_b = norm(full(Atf'*xr(1:size(Atf,1)) - bf))/max(nb,eps);
x = sol.solinfo.x(:);
off = S.Kf;     emin = inf;     emax = -inf;
Ks = cone_order_blocks(sol);
for k = 1:numel(Ks)
    N = Ks(k);
    if off+N^2 > numel(x),  break,  end
    X = reshape(x(off+(1:N^2)),N,N);    ev = eig((X+X')/2);
    emin = min(emin,min(ev));   emax = max(emax,max(ev));
    off = off + N^2;
end
psd_min = emin;     psd_relmin = emin/max(emax,eps);
trivial = norm(x)<=1e-12 || abs(rel_b-1)<=1e-6;
end

function Ks = cone_order_blocks(sol)
% PSD block sizes in the order sossolve lays them out (var 'sos', then
% extravar), unlike cx_shape's sorted multiset.
Ks = [];
for i = 1:sol.var.num
    if strcmp(sol.var.type{i},'sos')
        Ks(end+1) = round(sqrt(sol.var.idx{i+1}-sol.var.idx{i}));           %#ok<AGROW>
    end
end
if isfield(sol,'extravar') && isfield(sol.extravar,'num')
    for i = 1:sol.extravar.num
        Ks(end+1) = round(sqrt(sol.extravar.idx{i+1}-sol.extravar.idx{i})); %#ok<AGROW>
    end
end
end
