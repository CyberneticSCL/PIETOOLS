% Test_polyopvar_varname_column.m -- verify that polyopvar.varname is a
% column for all routes.
%
% The names of the state variables, polyopvar.varname, are stored as a
% p x 1 column (see the class documentation), whatever the route
% (constructor, vertcat/plus via common_states, subs, clean_poly, ...)
% and the orientation of the input, e.g. PIE.varname_pie, a 1 x p row.
%
% Each check uses at least 2 states, as for 1 state a row and a column
% are the same array. An object whose varname was set as a row
% (assignment to the property, as in several tests) is used to check that
% every function that sets or returns varname gives a column regardless
% of the orientation it receives.
%
% DJ (CC), 10/15/2026: Initial coding;

pvar s1 s2

%%%% Constructor: row, column, declaration syntax, defaults
x = polyopvar({'b','a'},[s1,s2],[0,1]);
chk_col(x,'constructor: row of names');
if ~isequal(x.varname,{'a';'b'})
    error("constructor: names must be sorted and stored as column, got %s.",mat2str(size(x.varname)))
end
chk_col(polyopvar({'b';'a'},s1,[0,1]),'constructor: column of names');
chk_col(polyopvar({'c','b','a'},{'s1','s2'},[0,1],true(3,2)),'constructor: 3 names, vars given as cellstr');
chk_col(polyopvar({'a','b'},[1;2],s1,[0,1],[true;false]),'constructor: dims and depmat given');
chk_col(polyopvar(),'default object');
if ~isequal(size(polyopvar().varname),[0,1])
    error("The default varname must be 0x1.")
end
polyopvar y1(s1) y2(s1,s2)
chk_col([y1;y2],'declaration syntax, vertcat');

%%%% Functions combining states (common_states): different and equal names,
%%%% with names set as row (property assignment) or column
a = polyopvar('a',s1,[0,1]);     b = polyopvar('b',s1,[0,1]);     c = polyopvar('c',s1,[0,1]);
chk_col([a;b;c],'vertcat');
chk_col(a+b,'plus');
chk_col(a.*b,'times');
chk_col(a.*b+c.*c,'plus of products');
chk_col(a-b,'minus');
chk_col(-(a+b),'uminus');
ab = a+b;       ab_row = ab;        ab_row.varname = ab.varname';       % names set as row
chk_col(ab_row+ab,'plus, one object with row names (equal names)');
chk_col(ab_row+ab_row,'plus, both objects with row names (equal names)');
chk_col(ab_row+c,'plus, row names, different names');
[A2,B2] = common_states(ab_row,ab_row);
chk_col(A2,'common_states A, equal row names');         chk_col(B2,'common_states B, equal row names');
[A2,B2] = common_states(ab_row,c);
chk_col(A2,'common_states A, different names');        chk_col(B2,'common_states B, different names');
if ~isequal(A2.varname,{'a';'b';'c'})
    error("common_states must return the sorted union of the names as column.")
end

%%%% initialize: names set as row
f = ab;     f.varname = f.varname';     f.chkval = false;
f = initialize(f);
chk_col(f,'initialize');

%%%% clean_poly and remove_unused_states: names set as row
[A2,~] = common_states(a.*a,b);                % states a,b; only a appears in the monomials
A2.varname = A2.varname';        A2.chkval = true;      % (chkval=true: initialize would set the column)
chk_col(remove_unused_states(A2),'remove_unused_states');
chk_col(clean_poly(A2),'clean_poly');
if ~isequal(remove_unused_states(A2).varname,{'a'}) || ~isequal(clean_poly(A2).varname,{'a'})
    error("The unused state b must be removed.")
end
% ... with 2 states left
[A3,~] = common_states(a.*b,c);                % states a,b,c; only a,b appear
A3.varname = A3.varname';        A3.chkval = true;
r1 = remove_unused_states(A3);      r2 = clean_poly(A3);
chk_col(r1,'remove_unused_states, 2 states left');
chk_col(r2,'clean_poly, 2 states left');
if ~isequal(r1.varname,{'a';'b'}) || ~isequal(r2.varname,{'a';'b'})
    error("The states a,b must remain, as column.")
end

%%%% subs: renaming states, names set as row
g = a.*b;       g.varname = g.varname';
g1 = subs(g,{'a'},{'d'});
chk_col(g1,'subs, renaming');
if ~isequal(g1.varname,{'b';'d'})
    error("subs must return the sorted names as column, got %s.",strjoin(g1.varname,','))
end
g2 = subs(g,{'a','b'},{'u','v'});
chk_col(g2,'subs, renaming all states');

%%%% ispvar
ab_p = [a;b];
[is_p,names] = ispvar(ab_p);
if ~is_p || ~isequal(names,{'a';'b'})
    error("ispvar must return the names as column.")
end
abr = ab_p;     abr.varname = abr.varname';
[~,names] = ispvar(abr);
if ~isequal(names,{'a';'b'})
    error("ispvar must return the names as column also for row names.")
end

%%%% dmonomials, innerprod, constant2polyopvar
chk_col(dmonomials([a;b],[1;2]),'dmonomials');
chk_col(innerprod(a,a)+innerprod(b,b),'innerprod');
chk_col(innerprod(a,a)+innerprod(b,b)+1,'functional with a constant');
chk_col(constant2polyopvar(1)+constant2polyopvar(2),'constants (no states)');

%%%% The state names of a PIE (PIE.x, varname_pie, varname_pde)
clear stateNameGenerator
pvar s t
u1 = pde_var(s,[0,1]);      u2 = pde_var(s,[0,1]);
PDE = [diff(u1,t)==diff(u1,s,2)+u1*u2;
       diff(u2,t)==diff(u2,s,2)-u1*u2;
       subs(u1,s,0)==0;  subs(u1,s,1)==0;  subs(u2,s,0)==0;  subs(u2,s,1)==0];
out = evalc('PIE = convert(PDE);');         %#ok<NASGU>
if numel(PIE.varname_pie)<2
    error("Test setup: the PIE must have at least 2 states.")
end
chk_col(PIE.x,'PIE.x');
if size(PIE.varname_pie,2)~=1 || size(PIE.varname_pde,2)~=1
    error("PIE.varname_pie and PIE.varname_pde must be columns.")
end
if ~isequal(PIE.x.varname,sort(PIE.varname_pie))
    error("PIE.x.varname must be the sorted PIE.varname_pie.")
end

disp('Test_polyopvar_varname_column passed: varname is a column for all routes.')


function chk_col(f,label)
% f.varname must be a p x 1 cell (0 x 1 for no states) of char arrays
v = f.varname;
if ~iscell(v) || size(v,2)~=1 || ~ismatrix(v) || ~iscellstr(v)
    error("Test_polyopvar_varname_column failed (%s): varname is %s, not a column cellstr.",label,mat2str(size(v)))
end
end
