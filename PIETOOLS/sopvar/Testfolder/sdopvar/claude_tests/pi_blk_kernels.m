function Kc = pi_blk_kernels(B,dval)
% KC = PI_BLK_KERNELS(B,DVAL) returns the kernels of ONE block of a 'copvar'
% or 'cdopvar' container, as matrix-valued 'polynomial' objects, whichever
% class the block happens to be.
%
% A container's blocks need not all be decision operators: adding a fixed
% operator to a decision one leaves a 'sopvar' block among 'sdopvar' ones,
% and a 'copvar' has only fixed blocks. 'pi_sdopvar_kernels' covers the
% decision case; this wraps it and handles the fixed case, which stores its
% coefficients directly and has no decision variables to substitute.
%
% INPUTS
% - B:      'sopvar' or 'sdopvar' object, one block of a container;
% - dval:   numel(B.Zd) x 1 array of decision variable values, in the order
%           of B.Zd. Ignored, and may be anything, for a fixed 'sopvar'
%           block;
%
% OUTPUTS
% - Kc:     cell array of the same size as the block's parameter cell,
%           holding each kernel as a B.dims(1) x B.dims(2) 'polynomial' in
%           the output variables B.vars.out and the dummy input variables,
%           the latter named by appending '_dum' to B.vars.in. The cell is
%           indexed over the block's SHARED variables only, so it has
%           3^numel(intersect(B.vars.in,B.vars.out)) entries -- one for a
%           cross block whose two spaces share nothing.
%
% See also PI_SDOPVAR_KERNELS, PI_MONOM_VECTOR, PI_GAMMA_INDEX.
%
% MMP, 09/21/2026: Initial coding
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change.

if isa(B,'sdopvar')
    Kc = pi_sdopvar_kernels(B,dval);
    return
end
if ~isa(B,'sopvar')
    error("A container block must be an 'sopvar' or an 'sdopvar'; got '"...
          +class(B)+"'.")
end

m1 = B.dims(1);     m2 = B.dims(2);
NL = prod([cellfun(@numel,B.ZL),1]);
NR = prod([cellfun(@numel,B.ZR),1]);
ZLp = pi_monom_vector(B.ZL,B.vars.out);
ZRp = pi_monom_vector(B.ZR,strcat(B.vars.in,'_dum'));

Kc = cell(size(B.params));
for k = 1:numel(B.params)
    Ck = reshape(full(B.params{k}),m1*NL,m2*NR);
    Kk = polynomial(zeros(m1,m2));
    for p = 1:m1
        for q = 1:m2
            Kk(p,q) = ZLp.'*(Ck((p-1)*NL+(1:NL),(q-1)*NR+(1:NR))*ZRp);
        end
    end
    Kc{k} = Kk;
end

end
