function [h,txt,file] = t1d_subfun(src,name,outdir)                         % CC, 09/22/2026
% t1d_subfun(src,name) -- lift the subfunction `name` out of the .m file `src`,
% write it to its own file in a temporary directory, and return a HANDLE to it.
%
% WHY.  Two of the behaviours this suite must guard live in subfunctions of
% pielr_certify.m -- padw, whose zero fill is defect B2, and fitw, whose random
% fill is the repair the .w0 hook depends on (B7).  MATLAB gives no way to call
% a subfunction from outside its file, and pielr_certify itself is 2-D only and
% far too expensive to run here.  Reading the source and asserting on the text
% would test a comment; extracting and EXECUTING it tests the shipped code.
%
% The slice runs from the `function ... name(` line to the LAST `end` at column
% 1 before the next `function` line (or end of file).  Taking the FIRST such
% `end` is wrong and was measured to be: both fitw and padw close their `for`
% loop with an unindented `end`, so a first-match slice silently dropped the
% function's own terminator -- harmless for those two only because MATLAB lets
% a single-function file omit it, and a silent truncation waiting to happen the
% moment either gains a statement after the loop.  If the file is reformatted
% past this rule the extraction errors rather than grabbing the wrong lines,
% which is the correct outcome: the suite would otherwise be testing a function
% it did not actually find.
%
% OUTPUT  h    handle to the extracted function
%         txt  its source text, for the copy-equality assertion in T9
%         file where it was written (caller may add its folder to the path)
if nargin<3 || isempty(outdir), outdir = tempname; end
if ~exist(outdir,'dir'), mkdir(outdir); end
L = strsplit(fileread(src),newline);
% the declaration: `function <out> = <name>(` anywhere in the file
pat = ['^\s*function\s+.*\<' name '\s*\('];
i0 = find(~cellfun('isempty',regexp(L,pat,'once')),1);
if isempty(i0)
    error('t1d_subfun: no subfunction %s in %s',name,src);
end
% the terminator: the LAST unindented `end` before the next function header.
% `\s` and not `\b` after `function`: MATLAB's regexp does NOT implement `\b`
% as a word boundary (it has `\<` and `\>` instead), so `^\s*function\b`
% matches NOTHING -- measured, and it silently ran the slice to end of file and
% swept up padw, tern and getf with it.
nx = i0 + find(~cellfun('isempty',regexp(L(i0+1:end),'^\s*function\s','once')),1);
if isempty(nx), nx = numel(L)+1; end
ce = i0 + find(~cellfun('isempty',regexp(L(i0+1:nx-1),'^end\s*$','once')));
if isempty(ce)
    error('t1d_subfun: no unindented `end` terminates %s in %s',name,src);
end
i1 = ce(end);
txt  = strjoin(L(i0:i1),newline);
file = fullfile(outdir,[name '.m']);
fid  = fopen(file,'w');
if fid<0, error('t1d_subfun: cannot write %s',file); end
fprintf(fid,'%s\n',txt);   fclose(fid);
addpath(outdir);
h = str2func(name);
% Confirm the handle resolves to the file just written, not to a same-named
% function already on the path -- the whole point is to exercise the SHIPPED
% body, and silently calling something else would invert the test's meaning.
w = which(name);
if ~strcmp(w,file)
    error('t1d_subfun: %s resolves to %s, not the extracted %s',name,w,file);
end
end
