function code = lb_blocks(m, libdir, outdir)
% CODE = LB_BLOCKS(M, LIBDIR, OUTDIR) uncomments one example block of the
% DDE or NDS/DDF library and writes it to OUTDIR/<id>.m.
%
% INPUT
%   m:      one entry of lb_meta (lib 'DDE' or 'NDSDDF'), with the block
%           line range m.lines and the single-comment prose lines m.prose
%   libdir: folder of the library files (PIETOOLS_examples)
%   outdir: folder for the uncommented block files
%
% OUTPUT
%   code:   char, the uncommented block (newline separated)
%
% NOTES
%   - Exactly one comment level is removed from every line in the range,
%     except the prose lines listed in m.prose: a line '% code' becomes
%     ' code', a line '% % note' becomes ' % note' (still a comment).
%   - Lines without a leading '%' are copied unchanged. A block left
%     active in the library (DDE26, DDF03) lists its whole range in
%     m.prose, so its own one-level comments and the commented
%     alternatives inside it stay commented.
%   - m.prose was set by reading each block; a prose line with one '%'
%     would otherwise run as command syntax.

txt = splitlines(fileread(fullfile(libdir,m.fn)));
a = m.lines(1);  b = m.lines(2);
out = cell(b-a+1,1);
for n = a:b
    L = txt{n};
    if ~ismember(n,m.prose)
        L = regexprep(L,'^(\s*)%','$1','once');
    end
    out{n-a+1} = L;
end
code = strjoin([{sprintf('%% %s lines %d-%d, one comment level removed (lb_blocks)',m.fn,a,b)}; out],newline);
if ~exist(outdir,'dir'), mkdir(outdir); end
fid = fopen(fullfile(outdir,[m.id,'.m']),'w');
fprintf(fid,'%s\n',code);
fclose(fid);
end
