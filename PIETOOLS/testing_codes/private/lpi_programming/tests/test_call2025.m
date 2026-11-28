function R = test_call2025(record)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_CALL2025(RECORD) checks that a call of poslpivar, lpivar or
% mat2opvar in a PIETOOLS 2025 form stops with the identifier
% poslpivar:call2025, lpivar:call2025 or mat2opvar:call2025, and that the
% current forms of their headers run; errors at the first failed check. On
% the programs lpiprogram of s on [0,1] and of (s1,s2) on [0,1] x [0,2]:
%   (1) 31 calls in the forms of the 2025 headers, N the dimensions of a
%       2025 'opvar' and D its degrees:
%       - poslpivar(prog,N), (prog,N,D), (prog,N,D,options): N = [1;2] and
%         [1 1;2 2]; D = 2, [1 2 3], [] and {1,[1 1 2],[1 1 2]}; in 2-D,
%         N = [1;2;2;1] and D the 'struct' with fields dx, dy, d2 (11);
%       - lpivar(prog,N), (prog,N,D), (prog,N,D,opts): N = [1 1;2 2], D = 2,
%         [1 2 3] and []; in 2-D, N = [1 1;2 2;0 0;1 1] and the 'struct' D
%         (8);
%       - mat2opvar(Pmat,Pdim), (Pmat,Pdim,vars), (Pmat,Pdim,dom),
%         (Pmat,Pdim,vars,dom), (Pmat,Pdim,dom,vars): Pmat 'double',
%         'polynomial' and 'dpvar'; Pdim [1 1;2 2], [1;2] and, in 2-D,
%         [1 1;2 2;0 0;0 0]; vars [s,s_dum] as 'polynomial' and 'cellstr'
%         (12);
%       each raises the identifier of its function;
%   (2) 23 calls in the current forms of the headers: counts and spaces
%       with and without DOM and DEG; spaces given as [], {}, a char, a
%       'polynomial' and a 'struct'; the operator form; no space (empty
%       counts); 2-D spaces; the two examples of mat2opvar. Each raises
%       nothing and returns the class and the size (components) stated.
%
% INPUT
% - record:  (optional) logical, default false. True skips every check and
%            records what each call of (1) raises and what each call of
%            (2) returns, for a comparison with an older tree.
%
% OUTPUT
% - R:       'struct' with fields
%              n     'double', the number of checks passed;
%              old   31 x 1 cell, the labels of (1);
%              errs  31 x 1 'struct', fields identifier and message ('char'),
%                    the error each call of (1) raised ('' for none);
%              cur   23 x 1 cell, the labels of (2);
%              out   23 x 1 cell, each {prog, Pop} (poslpivar, lpivar) or
%                    {Pop} (mat2opvar), the outputs of (2); with RECORD,
%                    {identifier, message} where a call of (2) raised one.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_call2025
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
% Initial coding MMP, 11/27/2026

if nargin<1,    record = false;     end
% Programs in s on [0,1] and in (s1,s2) on [0,1] x [0,2]; a scalar dpvar g
% for the decision forms of mat2opvar
s = polynomial({'s'});      sd = polynomial({'s_dum'});
p1 = lpiprogram(s,[],[0 1]);
p2 = lpiprogram(polynomial({'s1';'s2'}),[],[0 1;0 2]);
[~,g] = lpidecvar(p1,'g');
% 2025 inputs: the degrees of poslpivar_2d and lpivar_2d (a 'struct' of
% cells), the options of poslpivar, a 1-D degree cell, a 1-D dimension
% array, the variable pairs of mat2opvar
d2 = struct('dx',{{1,[1 1 1],[1 1 1]}},'dy',{{1,[1 1 1],[1 1 1]}},'d2',{cell(3,3)});
o25 = struct('psatz',0,'exclude',[0 0 0 0],'sep',0);
dc = {1,[1 1 2],[1 1 2]};   N2 = [1 1;2 2];     N4 = [1 1;2 2;0 0;1 1];     M = eye(3);
v25 = [s,sd];   c25 = {'s','s_dum'};
v2d = [polynomial({'s1'}),polynomial({'s1_dum'});polynomial({'s2'}),polynomial({'s2_dum'})];

% (1) The 2025 forms: label, function, inputs
old = {
 "poslpivar(prog,[1;2])",                       @poslpivar, {p1,[1;2]};
 "poslpivar(prog,[1 1;2 2])",                   @poslpivar, {p1,N2};
 "poslpivar(prog,[1;2],2)",                     @poslpivar, {p1,[1;2],2};
 "poslpivar(prog,[1;2],[1 2 3])",               @poslpivar, {p1,[1;2],[1 2 3]};
 "poslpivar(prog,[1;2],[])",                    @poslpivar, {p1,[1;2],[]};
 "poslpivar(prog,[1;2],{1,[1 1 2],[1 1 2]})",   @poslpivar, {p1,[1;2],dc};
 "poslpivar(prog,[1 1;2 2],{...},options)",     @poslpivar, {p1,N2,dc,o25};
 "poslpivar(prog,[1;2],2,options)",             @poslpivar, {p1,[1;2],2,o25};
 "poslpivar(prog2,[1;2;2;1])",                  @poslpivar, {p2,[1;2;2;1]};
 "poslpivar(prog2,[1;2;2;1],d2d)",              @poslpivar, {p2,[1;2;2;1],d2};
 "poslpivar(prog2,[1;2;2;1],d2d,options)",      @poslpivar, {p2,[1;2;2;1],d2,o25};
 "lpivar(prog,[1 1;2 2])",                      @lpivar,    {p1,N2};
 "lpivar(prog,[1 1;2 2],2)",                    @lpivar,    {p1,N2,2};
 "lpivar(prog,[1 1;2 2],[1 2 3])",              @lpivar,    {p1,N2,[1 2 3]};
 "lpivar(prog,[1 1;2 2],[])",                   @lpivar,    {p1,N2,[]};
 "lpivar(prog,[1 1;2 2],[1 2 3],opts)",         @lpivar,    {p1,N2,[1 2 3],struct()};
 "lpivar(prog2,[1 1;2 2;0 0;1 1])",             @lpivar,    {p2,N4};
 "lpivar(prog2,[1 1;2 2;0 0;1 1],d2d)",         @lpivar,    {p2,N4,d2};
 "lpivar(prog2,[1 1;2 2;0 0;1 1],d2d,opts)",    @lpivar,    {p2,N4,d2,struct()};
 "mat2opvar(M,[1 1;2 2])",                      @mat2opvar, {M,N2};
 "mat2opvar(M,[1;2])",                          @mat2opvar, {M,[1;2]};
 "mat2opvar(M,[1 1;2 2],[s,s_dum])",            @mat2opvar, {M,N2,v25};
 "mat2opvar(M,[1 1;2 2],{'s','s_dum'})",        @mat2opvar, {M,N2,c25};
 "mat2opvar(M,[1 1;2 2],[0 1])",                @mat2opvar, {M,N2,[0 1]};
 "mat2opvar(M,[1 1;2 2],[s,s_dum],[0 1])",      @mat2opvar, {M,N2,v25,[0 1]};
 "mat2opvar(M,[1 1;2 2],{'s','s_dum'},[0 1])",  @mat2opvar, {M,N2,c25,[0 1]};
 "mat2opvar(M,[1 1;2 2],[0 1],[s,s_dum])",      @mat2opvar, {M,N2,[0 1],v25};
 "mat2opvar(M,[1;2],[s,s_dum],[0 1])",          @mat2opvar, {M,[1;2],v25,[0 1]};
 "mat2opvar(s*M,[1 1;2 2],[s,s_dum],[0 1])",    @mat2opvar, {s*M,N2,v25,[0 1]};
 "mat2opvar(g*M,[1 1;2 2],[s,s_dum],[0 1])",    @mat2opvar, {g*M,N2,v25,[0 1]};
 "mat2opvar(M,[1 1;2 2;0 0;0 0],vars2d,dom2d)", @mat2opvar, {M,[N2;0 0;0 0],v2d,[0 1;0 2]}};
n = 0;
errs = repmat(struct('identifier','','message',''),size(old,1),1);
for k = 1:size(old,1)
    f = old{k,2};   a = old{k,3};
    errs(k) = err_of(@() f(a{:}));
    if record,  continue;   end
    id = [func2str(f),':call2025'];
    assert(strcmp(errs(k).identifier,id),"Case %d, %s: raised '%s' (%s), not %s.",...
           k,old{k,1},errs(k).identifier,errs(k).message,id);
    n = n+1;
end

% (2) The current forms: label, function, inputs, class and size of Pop.
% T: a fixed operator on R x L2^2[s] for the operator forms
T = mat2opvar(eye(3),[1;2],{{},{'s'}},[0 1]);
sp = {{},{'s'}};    sp2 = {{'s1'},{'s1','s2'}};     dom2 = [0 1;0 2];
cur = {
 "poslpivar(prog,[1;2],sp,[0 1],2)",            @poslpivar, {p1,[1;2],sp,[0 1],2},  'dopvar', [3 3];
 "poslpivar(prog,[1;2],sp,2)",                  @poslpivar, {p1,[1;2],sp,2},        'dopvar', [3 3];
 "poslpivar(prog,[1;2],sp)",                    @poslpivar, {p1,[1;2],sp},          'dopvar', [3 3];
 "poslpivar(prog,3,[])",                        @poslpivar, {p1,3,[]},              'dopvar', [3 3];
 "poslpivar(prog,2,{})",                        @poslpivar, {p1,2,{}},              'dopvar', [2 2];
 "poslpivar(prog,2,'s',[0 1],1)",               @poslpivar, {p1,2,'s',[0 1],1},     'dopvar', [2 2];
 "poslpivar(prog,2,s,[0 1],1)",                 @poslpivar, {p1,2,s,[0 1],1},       'dopvar', [2 2];
 "poslpivar(prog,T,1)",                         @poslpivar, {p1,T,1},               'dopvar', [3 3];
 "poslpivar(prog,[],{})",                       @poslpivar, {p1,[],{}},             'dopvar', [0 0];
 "poslpivar(prog2,[1;1],sp2,dom2,1)",           @poslpivar, {p2,[1;1],sp2,dom2,1},  'dopvar', [2 2];
 "lpivar(prog,[1;2],sp,[0 1],2)",               @lpivar,    {p1,[1;2],sp,[0 1],2},  'dopvar', [3 3];
 "lpivar(prog,struct(out,in),sp,[0 1],[1 2 3])",@lpivar,    {p1,struct('out',[1;2],'in',[2;1]),sp,[0 1],[1 2 3]}, 'dopvar', [3 3];
 "lpivar(prog,struct(out,in),struct(out,in),[0 1],1)", @lpivar, {p1,struct('out',2,'in',1), ...
                    struct('out',{{{'s'}}},'in',{{{}}}),[0 1],1},           'dopvar', [2 1];
 "lpivar(prog,T,2)",                            @lpivar,    {p1,T,2},               'dopvar', [3 3];
 "lpivar(prog,[1;2],sp)",                       @lpivar,    {p1,[1;2],sp},          'dopvar', [3 3];
 "lpivar(prog,1,{})",                           @lpivar,    {p1,1,{}},              'dopvar', [1 1];
 "lpivar(prog2,[1;1],sp2,dom2,1)",              @lpivar,    {p2,[1;1],sp2,dom2,1},  'dopvar', [2 2];
 "mat2opvar(eye(3),[1;2],sp,[0 1])",            @mat2opvar, {eye(3),[1;2],sp,[0 1]}, 'opvar', [3 3];
 "mat2opvar(eye(2),2,{{}},[])",                 @mat2opvar, {eye(2),2,{{}},[]},     'opvar',  [2 2];
 "mat2opvar(-g*eye(2),2,{{}},[])",              @mat2opvar, {-g*eye(2),2,{{}},[]},  'dopvar', [2 2];
 "mat2opvar(ones(2,3),struct,struct,[0 1])",    @mat2opvar, {ones(2,3),struct('out',[1;1],'in',[1;2]), ...
                    struct('out',{sp},'in',{sp}),[0 1]},                    'opvar',  [2 3];
 "mat2opvar(eye(2),2,'s',[0 1])",               @mat2opvar, {eye(2),2,'s',[0 1]},   'opvar',  [2 2];
 "mat2opvar(eye(2),2,{'s'},[0 1],mult_only)",   @mat2opvar, {eye(2),2,{'s'},[0 1],struct('mult_only',true)}, 'opvar', [2 2]};
out = cell(size(cur,1),1);
for k = 1:size(cur,1)
    f = cur{k,2};   a = cur{k,3};
    try
        if strcmp(func2str(f),'mat2opvar'),     Pop = f(a{:});      out{k} = {Pop};
        else,   [pr,Pop] = f(a{:});     out{k} = {pr,Pop};
        end
    catch ME
        if record,  out{k} = {ME.identifier,ME.message};    continue;   end
        error('test_call2025:current',"Case %d, %s raised '%s': %s",...
              k,cur{k,1},ME.identifier,ME.message);
    end
    if record,  continue;   end
    assert(isa(Pop,cur{k,4}) && isequal(size(Pop),cur{k,5}),...
           "Case %d, %s: a %s of size %s, not a %s of size %s.",k,cur{k,1},...
           class(Pop),mat2str(size(Pop)),cur{k,4},mat2str(cur{k,5}));
    n = n+1;
end
R = struct('n',n,'old',{old(:,1)},'errs',errs,'cur',{cur(:,1)},'out',{out});
if ~record
    fprintf('test_call2025: %d checks passed (%d calls in 2025 forms, %d in current forms).\n',...
            n,size(old,1),size(cur,1));
end
end


function e = err_of(f)
% E = ERR_OF(F) returns the message and identifier of the error that F()
% raises (both '' when it raises none).
% INPUT  f: 'function_handle', called with no argument.
% OUTPUT e: 'struct' with fields identifier and message ('char').
e = struct('identifier','','message','');
try,    f();    catch ex,   e = struct('identifier',ex.identifier,'message',ex.message);   end
end
