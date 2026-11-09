function s = p0_infonote(info)
% S = P0_INFONOTE(INFO) returns the part of a container executive's info
% struct that Phase 0 records per row, as text: the classifier status, the
% SDP shape from lpi_shape_sop (equality rows m, PSD block sides Ks, primal
% length nx, decision variables ndv), and the numerical witness (maxre for
% stability, gain or h2 otherwise, NaN when absent).
s = '';
if ~isstruct(info),     return,     end
if isfield(info,'status'),      s = [s,'status=',char(info.status),' '];    end
if isfield(info,'shape') && isstruct(info.shape)
    sh = info.shape;
    s = [s,sprintf('m=%d nx=%d ndv=%d Ks=%s ', sh.m, sh.nx, sh.ndv, mat2str(sh.Ks))];
end
w = struct('maxre',NaN,'gain',NaN,'h2',NaN,'ok',NaN);
if isfield(info,'maxre'),   w.maxre = double(info.maxre);   end
if isfield(info,'witness') && isstruct(info.witness)
    W = info.witness;
    for f = {'gain','h2','ok'}
        if isfield(W,f{1}) && isnumeric(W.(f{1})) || (isfield(W,f{1}) && islogical(W.(f{1})))
            w.(f{1}) = double(W.(f{1}));
        end
    end
end
s = [s,sprintf('witness_ok=%g maxre=%.10g wgain=%.10g wh2=%.10g',w.ok,w.maxre,w.gain,w.h2)];
end
