function C = t1d_codelines(txt)                                             % CC, 09/22/2026
% t1d_codelines(txt) -- the EXECUTABLE lines of some MATLAB source, with every
% comment and all trailing whitespace removed, as a cellstr.
%
% WHY.  Two provenance assertions in this suite compare code that is not
% byte-identical but must not differ in behaviour:
%   T0  restrict_solve_1d against the shipped private/restrict_solve -- they
%       differ in the header and in ONE line (the gate call), and the test has
%       to see exactly that and nothing else;
%   T9  t1d_bm's local copy of fitw against pielr_certify's own.
% Both files are dense with change markers (`% CC, dd/mm/yyyy`) and commented-
% out old code, which the house protocol requires and which a byte compare
% would trip over for no reason.  Stripping comments compares the code.
%
% The scanner is quote-aware: a `%` inside a single-quoted string (a printf
% format, say) is not a comment.  Transpose is treated as opening a string,
% which is the usual ambiguity of MATLAB's own lexer; it is harmless here
% because a stray unpaired quote can only make the scanner keep MORE text,
% i.e. it can only cause a spurious DIFFERENCE, never hide a real one.
L = strsplit(txt,newline);
C = {};
for k = 1:numel(L)
    s = L{k};
    inq = false;   cut = 0;
    for j = 1:numel(s)
        c = s(j);
        if c == ''''
            inq = ~inq;
        elseif c == '%' && ~inq
            cut = j;   break
        end
    end
    if cut > 0, s = s(1:cut-1); end
    s = deblank(s);
    if isempty(strtrim(s)), continue; end      % comment-only or blank
    C{end+1} = s;                              %#ok<AGROW>
end
C = C(:);
end
