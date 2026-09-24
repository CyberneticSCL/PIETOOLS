function L = pielr_layout(sos)                                              % CC, 09/23/2026
% PIELR_LAYOUT  Map the SeDuMi x-vector of an LPI program into free
% coordinates and PSD (Gram) blocks, making NO assumption about their order.
%
% WHY THIS EXISTS.  raw_data models the x-vector as
%     [ Kf free coordinates ][ Gram block 1 ][ Gram block 2 ] ...
% with Kf = sos.var.idx{1}-1 and the Gram blocks contiguous from there.  That
% is true for the stability LPIs and FALSE for any program built with lpivar:
% a free operator variable is stored as 'poly' entries in sos.var, INTERLEAVED
% between the 'sos' entries.  MEASURED on the 1-D l2gain program of
% PIETOOLS_Hinf_gain (Ex_Transport_Eq_with_Disturbance, 'light'):
%     x(1)        gam                     free
%     x(2..101)   sos   N=10              Gram  (Rop)
%     x(102..144) poly  43 coordinates    free  (Qop, from lpivar)
%     x(145..433) sos   N=17              Gram  (De1op)
%     x(434..497) sos   N=8               Gram  (De2op)
% raw_data returns Kf=1, Ns=[10 17 8], so Kf+sum(Ns.^2) = 454 against 497
% actual coordinates: the 43 Qop coordinates are silently dropped, and every
% Gram block after the first is then indexed at the wrong offset.
%
% This routine walks sos.var.idx / sos.var.type and sos.extravar instead, so
% the map is read from the program rather than assumed.
%
% OUTPUT (struct)
%   L.free   column vector of free coordinate indices (any type but 'sos')
%   L.rows   1 x B cell, L.rows{i} the coordinate indices of Gram block i
%   L.N      1 x B, side length of each Gram block   (numel(rows{i}) == N^2)
%   L.Ntot   total number of coordinates
%   L.Kf     numel(L.free), kept for callers that only need the count
%
% NEUTRALITY.  On a program whose Gram blocks ARE contiguous after a free
% prefix -- every stability program this package builds -- L.free is exactly
% 1:Kf, L.N is exactly Ns, and L.rows{i} is exactly the offset range raw_data
% computes.  Asserted in the T-layout regression test.

ntot = size(sos.expr.At{1},1);
mask = false(ntot,1);                 % true where a coordinate is a Gram entry
N    = zeros(1,0);
rows = {};
for i = 1:sos.var.num
    lo = sos.var.idx{i};   hi = sos.var.idx{i+1}-1;
    if hi < lo, continue; end         % zero-length blocks occur and are not errors
    if strcmp(sos.var.type{i},'sos')
        n = hi-lo+1;   Ni = round(sqrt(n));
        if Ni*Ni ~= n
            error('pielr_layout:notSquare', ...
                  'sos block %d has length %d, not a perfect square',i,n);
        end
        N(end+1)    = Ni;             %#ok<AGROW>
        rows{end+1} = (lo:hi)';       %#ok<AGROW>
        mask(lo:hi) = true;
    end
end
for i = 1:sos.extravar.num
    lo = sos.extravar.idx{i};   hi = sos.extravar.idx{i+1}-1;
    if hi < lo, continue; end
    n = hi-lo+1;   Ni = round(sqrt(n));
    if Ni*Ni ~= n
        error('pielr_layout:notSquare', ...
              'extravar block %d has length %d, not a perfect square',i,n);
    end
    N(end+1)    = Ni;                 %#ok<AGROW>
    rows{end+1} = (lo:hi)';           %#ok<AGROW>
    mask(lo:hi) = true;
end
L.free = find(~mask);
L.rows = rows;
L.N    = N;
L.Ntot = ntot;
L.Kf   = numel(L.free);
end
