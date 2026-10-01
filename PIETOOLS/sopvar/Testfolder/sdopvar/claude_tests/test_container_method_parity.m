function test_container_method_parity(mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_CONTAINER_METHOD_PARITY checks that the parallel method files of
% @copvar and @cdopvar stay in step. The two classes duplicate their
% methods by design, as @sopvar and @sdopvar do ('derive_copvar_meta'), so
% a fix applied to one file has no mechanism carrying it to the other; this
% test is that mechanism.
%
% TEST_CONTAINER_METHOD_PARITY('print') prints each pair's actual diff as
% the MATLAB literal of EXPECTED below, to update it after a deliberate
% change to a decision-variable branch.
%
% CHECKS
% (1) The two class folders hold exactly the method files listed here: the
%     16 shared names, plus the class-specific files. A method added to one
%     class only fails the test.
% (2) For each shared name, both files are NORMALIZED - comments and block
%     comments removed (string-aware: a '%' inside a quoted literal is
%     kept), leading and trailing blanks removed, internal runs of blanks
%     collapsed, blank lines dropped, and the words 'copvar' and 'cdopvar'
%     replaced by one token - and compared line by line through a longest
%     common subsequence.
%     - 10 pairs must be IDENTICAL: cat, ctranspose, eq, minus,
%       numArgumentsFromSubscript, size, subsasgn, subsref, transpose,
%       uminus.
%     - 6 pairs must differ by EXACTLY the lines in EXPECTED, which are the
%       decision-variable branches: blkdiag, horzcat, vertcat (the sdopvar
%       promotion in copvar; the Zd reconciliation in cdopvar), mtimes
%       (the decision x decision guard, promotion, the Zd reconciliation,
%       'plus_batch'), plus (promotion, the Zd reconciliation and sharing)
%       and verify (the admitted block classes).
% The oracle for (2) is the source text of each pair, not either routine's
% behaviour, so it cannot be fooled by two copies that err alike: it
% detects that they stopped being copies.
%
% Block class names ('sopvar', 'sdopvar') are NOT mapped, so a branch on
% the block class shows as a difference.
%
% Initial coding MMP, 09/30/2026. Audit item C10: nothing checked that the
%                10 code-identical pairs stay identical or that the other 6
%                differ only in their decision-variable branches.
% MMP, 09/30/2026: EXPECTED.verify follows the verify input rename Mop -> P.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1,    mode = 'test';  end
printing = strcmpi(mode,'print');

dirA = fileparts(which('copvar'));      % the @copvar folder
dirB = fileparts(which('cdopvar'));     % the @cdopvar folder
assert(endsWith(dirA,[filesep '@copvar']) && endsWith(dirB,[filesep '@cdopvar']),...
    'copvar -> %s, cdopvar -> %s',dirA,dirB);

IDENT = {'cat','ctranspose','eq','minus','numArgumentsFromSubscript','size',...
         'subsasgn','subsref','transpose','uminus'};
DIFF  = {'blkdiag','horzcat','vertcat','mtimes','plus','verify'};
ONLY_A = {'copvar','canonicalize','copvar2nopvar','copvar2opvar','copvar2opvar2d'};
ONLY_B = {'cdopvar'};
nchk = 0;

% % % (1) The method sets.
fa = mfiles(dirA);      fb = mfiles(dirB);
assert(isequal(sort(fa),sort([IDENT,DIFF,ONLY_A])),...
    '@copvar method files differ from the list: extra {%s}, missing {%s}',...
    strjoin(setdiff(fa,[IDENT,DIFF,ONLY_A]),','),strjoin(setdiff([IDENT,DIFF,ONLY_A],fa),','));
assert(isequal(sort(fb),sort([IDENT,DIFF,ONLY_B])),...
    '@cdopvar method files differ from the list: extra {%s}, missing {%s}',...
    strjoin(setdiff(fb,[IDENT,DIFF,ONLY_B]),','),strjoin(setdiff([IDENT,DIFF,ONLY_B],fb),','));
nchk = nchk+2;

% % % (2) The pairs.
E = expected();
for name = [IDENT,DIFF]
    nm = name{1};
    a = normalize(fullfile(dirA,[nm '.m']));
    b = normalize(fullfile(dirB,[nm '.m']));
    [onlyA,onlyB] = lcs_diff(a,b);
    if printing
        fprintf('E.%s = {{ ...\n',nm);     print_lines(onlyA);
        fprintf('    }, { ...\n');          print_lines(onlyB);
        fprintf('    }};\n');
        continue
    end
    if any(strcmp(nm,IDENT))
        want = {cell(0,1),cell(0,1)};
    else
        want = E.(nm);
    end
    ok = isequal(onlyA(:),want{1}(:)) && isequal(onlyB(:),want{2}(:));
    assert(ok,['@copvar/%s.m and @cdopvar/%s.m differ from the expected '...
        'difference.\n  only in copvar:\n    %s\n  only in cdopvar:\n    %s\n'...
        'Run test_container_method_parity(''print'') after a deliberate change.'],...
        nm,nm,strjoin(onlyA,'\n    '),strjoin(onlyB,'\n    '));
    nchk = nchk+1;
end
if printing,    return,     end

% % % Controls: the comparison sees a one-token drift in an identical pair
% % % and an extra line in a differing one; the comment stripper keeps a
% % % '%' inside a literal and reads x' as a transpose.
a = normalize(fullfile(dirA,'uminus.m'));   b = a;      k = ceil(numel(a)/2);
b{k} = [b{k} ' X'];
[oa,ob] = lcs_diff(a,b);
assert(isequal(oa,a(k)) && isequal(ob,b(k)),'control: a one-token drift was not isolated');
a = normalize(fullfile(dirA,'plus.m'));     b = normalize(fullfile(dirB,'plus.m'));
b = [b(1:end-1); {'x = 1;'}; b(end)];
[oa,ob] = lcs_diff(a,b);
assert(isequal(oa(:),E.plus{1}(:)) && isequal(ob(:),[E.plus{2}(:); {'x = 1;'}]),...
    'control: an extra line in a differing pair was not isolated');
assert(strcmp(strip_comment('x = ''a%b''; % c'),'x = ''a%b''; ') && ...
       strcmp(strip_comment('y = x''; % c'),'y = x''; ') && ...
       strcmp(strip_comment('s = "p%q"; % c'),'s = "p%q"; ') && ...
       strcmp(strip_comment('z = [a'' ''%'']; % c'),'z = [a'' ''%'']; '),...
    'control: comment stripping');
nchk = nchk+3;

fprintf('test_container_method_parity passed (%d checks).\n',nchk);
end


% ========================================================================
function E = expected()
% The lines, after normalization, that one file of the pair has and the
% other lacks, in file order. CONT stands for 'copvar'/'cdopvar'.
E = struct();
E.blkdiag = {{ ...
    'isdec = cellfun(@(a) isa(a,''sdopvar''),varargin);'
    'if any(isdec)'
    'varargin(isdec) = cellfun(@(a) CONT({a}),varargin(isdec),''uni'',0);'
    'P = blkdiag(varargin{:});'
    'return'
    'end'
    '[C,meta] = cat_copvar_grid(''d'',varargin,''CONT'');'
    }, { ...
    '[C,meta,Zds,src] = cat_copvar_grid(''d'',varargin,''CONT'');'
    '[C,meta.Zd] = merge_dvar_lists(C,Zds,src);'
    }};
E.horzcat = {{ ...
    'isdec = cellfun(@(a) isa(a,''sdopvar''),varargin);'
    'if any(isdec)'
    'varargin(isdec) = cellfun(@(a) CONT({a}),varargin(isdec),''uni'',0);'
    'P = horzcat(varargin{:});'
    'return'
    'end'
    '[C,meta] = cat_copvar_grid(''h'',varargin,''CONT'');'
    }, { ...
    '[C,meta,Zds,src] = cat_copvar_grid(''h'',varargin,''CONT'');'
    '[C,meta.Zd] = merge_dvar_lists(C,Zds,src);'
    }};
E.vertcat = {{ ...
    'isdec = cellfun(@(a) isa(a,''sdopvar''),varargin);'
    'if any(isdec)'
    'varargin(isdec) = cellfun(@(a) CONT({a}),varargin(isdec),''uni'',0);'
    'P = vertcat(varargin{:});'
    'return'
    'end'
    '[C,meta] = cat_copvar_grid(''v'',varargin,''CONT'');'
    }, { ...
    '[C,meta,Zds,src] = cat_copvar_grid(''v'',varargin,''CONT'');'
    '[C,meta.Zd] = merge_dvar_lists(C,Zds,src);'
    }};
E.mtimes = {{ ...
    '''multiplication by a matrix changes the spaces and must be a CONT.''])'
    'error(''mtimes:badInput'',''Factors must be CONT objects, or one a scalar.'')'
    'CA = A.C; CB = B.C;'
    'if nt>=1'
    'meta.space_in = B.space_in; meta.dim_in = B.dim_in;'
    }, { ...
    '''multiplication by a matrix changes the spaces and must be a container.''])'
    'if isa(A,''CONT'') && isa(B,''CONT'')'
    'error(''CONT:decisionTimesDecision'',...'
    '[''Both factors are CONT, so the product would be quadratic in the ''...'
    '''decision variables and is not representable. Compose a decision ''...'
    '''operator only against a CONT, as in T''''*P*T.''])'
    'end'
    'if isa(A,''CONT''), A = CONT(A); end'
    'if isa(B,''CONT''), B = CONT(B); end'
    'error(''mtimes:badInput'',...'
    '''Factors must be CONT or CONT objects, or one a scalar.'')'
    '[CA,CB,Zd] = merge_dvar_pair(A.C,B.C,A.Zd,B.Zd);'
    'if nt==1'
    'elseif nt>1 && any(cellfun(@(t) isa(t,''sdopvar''),terms(1:nt)))'
    'Cc{i,j} = plus_batch(terms{1:nt});'
    'elseif nt>1'
    'Cc{i,j} = terms{1};'
    'meta.space_in = B.space_in; meta.dim_in = B.dim_in; meta.Zd = Zd;'
    }};
E.plus = {{ ...
    'error(''plus:badInput'',''Both summands must be CONT objects.'')'
    'CA = A.C; CB = B.C;'
    'meta = metadata(A);'
    }, { ...
    'if isa(A,''CONT''), A = CONT(A); end'
    'if isa(B,''CONT''), B = CONT(B); end'
    'error(''plus:badInput'',''Summands must be CONT or CONT objects.'')'
    '[CA,CB,Zd] = merge_dvar_pair(A.C,B.C,A.Zd,B.Zd);'
    'if isa(Cc{ii},''sdopvar''), Cc{ii}.Zd = Zd; end'
    'meta = metadata(A); meta.Zd = Zd;'
    }};
E.verify = {{ ...
    'info = verify_copvar_meta(P,''CONT'',{''sopvar''});'
    }, { ...
    'info = verify_copvar_meta(P,''CONT'',{''sopvar'',''sdopvar''});'
    }};
end


function f = mfiles(d)
% Method file names of a class folder, without '.m'.
s = dir(fullfile(d,'*.m'));
f = cellfun(@(n) n(1:end-2),{s.name},'uni',0);
end


function L = normalize(file)
% Code lines of FILE, comments removed, blanks trimmed and collapsed, the
% container class names mapped to CONT; a column cell of char.
txt = fileread(file);
raw = regexp(txt,'\r?\n','split');
L = cell(0,1);
inblock = false;
for k = 1:numel(raw)
    s = strtrim(raw{k});
    if strcmp(s,'%{'),  inblock = true;     continue,   end
    if inblock
        if strcmp(s,'%}'),  inblock = false;    end
        continue
    end
    s = strtrim(regexprep(strip_comment(s),'\s+',' '));
    if isempty(s),  continue,   end
    L{end+1,1} = regexprep(s,'\<c(d)?opvar\>','CONT');                     %#ok<AGROW>
end
end


function s = strip_comment(s)
% s up to its first '%' outside a quoted literal. A single quote opens a
% char literal unless it follows an identifier character, a closing
% bracket, a dot or another quote, where it is a transpose.
inS = false;    inD = false;    prev = ' ';     k = 1;
while k<=numel(s)
    c = s(k);
    if inS
        if c==''''
            if k<numel(s) && s(k+1)==''''
                k = k+1;                    % '' inside a char literal
            else
                inS = false;
            end
        end
    elseif inD
        if c=='"'
            if k<numel(s) && s(k+1)=='"'
                k = k+1;                    % "" inside a string literal
            else
                inD = false;
            end
        end
    elseif c=='%'
        s = s(1:k-1);
        return
    elseif c=='"'
        inD = true;
    elseif c==''''
        if ~(isletter(prev) || (prev>='0' && prev<='9') || any(prev=='_)]}.'''))
            inS = true;
        end
    end
    prev = c;
    k = k+1;
end
end


function [onlyA,onlyB] = lcs_diff(a,b)
% Lines of a not in a longest common subsequence with b, and vice versa,
% each in file order.
na = numel(a);      nb = numel(b);
T = zeros(na+1,nb+1);
for i = na:-1:1
    for j = nb:-1:1
        if strcmp(a{i},b{j})
            T(i,j) = T(i+1,j+1)+1;
        else
            T(i,j) = max(T(i+1,j),T(i,j+1));
        end
    end
end
onlyA = cell(0,1);      onlyB = cell(0,1);
i = 1;      j = 1;
while i<=na && j<=nb
    if strcmp(a{i},b{j})
        i = i+1;    j = j+1;
    elseif T(i+1,j)>=T(i,j+1)
        onlyA{end+1,1} = a{i};  i = i+1;                                    %#ok<AGROW>
    else
        onlyB{end+1,1} = b{j};  j = j+1;                                    %#ok<AGROW>
    end
end
onlyA = [onlyA; a(i:na)];      onlyB = [onlyB; b(j:nb)];
end


function print_lines(L)
for k = 1:numel(L)
    fprintf('    ''%s''\n',strrep(L{k},'''',''''''));
end
end
