% Test_polyopvar_subsasgn.m -- verify @polyopvar/subsasgn.m.
%
% Assigning a field of a 'polyopvar' outside of the methods of the class,
% e.g.
%   f.C = Cnew;        f.dims.in(2) = 1;        f.varname = names;
% sets f.chkval = false (subsasgn.m), so that methods which check
%   if ~f.chkval, f = initialize(f); end
% validate the new value. Assigning chkval itself just sets the value.
% Assignment to elements, f(i,j)=... or f{i,j}=..., is not supported.
%
% DJ (CC), 10/15/2026: Initial coding;

pvar s
x = polyopvar({'a','b'},s,[0,1]);
if ~x.chkval
    error("Test setup: a polyopvar declared by the constructor must have chkval=true.")
end

%%%% Every field: chkval becomes false, and the value is stored
vals = struct('C',x.C,'Z',x.Z,'varname',{x.varname},'vars',x.vars,'dom',x.dom,'dims',x.dims,'depmat',x.depmat);
fields = fieldnames(vals);
for k=1:numel(fields)
    fn = fields{k};
    y = x;
    y.(fn) = vals.(fn);                         % dynamic field
    if y.chkval
        error("f.%s = value must set chkval to false.",fn)
    end
    z = x;
    switch fn
        case 'C',       z.C = vals.C;
        case 'Z',       z.Z = vals.Z;
        case 'varname', z.varname = vals.varname;
        case 'vars',    z.vars = vals.vars;
        case 'dom',     z.dom = vals.dom;
        case 'dims',    z.dims = vals.dims;
        case 'depmat',  z.depmat = vals.depmat;
    end
    if z.chkval
        error("f.%s = value (static field name) must set chkval to false.",fn)
    elseif ~isequal(size(z.(fn)),size(vals.(fn)))
        error("f.%s = value did not store the value.",fn)
    end
end

%%%% Elements and nested fields
y = x;      y.Z(2,1) = 1;
if y.chkval || full(y.Z(2,1))~=1
    error("f.Z(2,1) = 1 must store the value and set chkval=false.")
end
y = x;      y.dims.in(2) = 3;
if y.chkval || y.dims.in(2)~=3
    error("f.dims.in(2) = 3 must store the value and set chkval=false.")
end
y = x;      y.dims.out = [1;1];
if y.chkval
    error("f.dims.out = ... must set chkval=false.")
end
y = x;      y.vars.in = {'t'};
if y.chkval || ~isequal(y.vars.in,{'t'})
    error("f.vars.in = ... must store the value and set chkval=false.")
end
y = x;      y.depmat.in(1,1) = false;
if y.chkval || y.depmat.in(1,1)
    error("f.depmat.in(1,1) = false must store the value and set chkval=false.")
end
y = x;      y.varname{1} = 'c';
if y.chkval || ~strcmp(y.varname{1},'c')
    error("f.varname{1} = 'c' must store the value and set chkval=false.")
end
y = x;      y.C.ops = y.C.ops;
if y.chkval
    error("f.C.ops = ... must set chkval=false.")
end

%%%% The value is the same as before for the other fields
y = x;      y.varname = {'b';'a'};
if ~isequal(y.Z,x.Z) || ~isequal(y.dims,x.dims) || ~isequal(y.C.dims,x.C.dims)
    error("Assigning f.varname must not change other fields.")
end

%%%% Assigning chkval: just sets the value
y = x;      y.Z = y.Z;
y.chkval = true;
if ~y.chkval
    error("f.chkval = true must set chkval to true.")
end
y.chkval = false;
if y.chkval
    error("f.chkval = false must set chkval to false.")
end

%%%% initialize checks the object and sets chkval=true
y = x;      y.Z = y.Z;
y = initialize(y);
if ~y.chkval
    error("initialize must set chkval to true.")
end

%%%% An inconsistent value is now caught by the methods that check chkval
% (dims.in(1)=5 does not match the identity operator of size 1 in C: initialize errors)
y = x;      y.dims.in(1) = 5;
try
    combine_terms(y);
    error("Test_polyopvar_subsasgn:noerror","combine_terms must check an object with changed fields.")
catch ME
    if strcmp(ME.identifier,"Test_polyopvar_subsasgn:noerror")
        rethrow(ME)
    elseif ~contains(ME.message,'does not match') && ~contains(ME.message,'conflicts')
        error("Unexpected error for inconsistent dims: %s",ME.message)
    end
end
% ... while a consistent change is accepted after initialize
y = x;      y.depmat.in = y.depmat.in;
y = combine_terms(y);
if ~y.chkval
    error("combine_terms must return a checked object.")
end

%%%% Element assignment is not supported
try
    y = x;      y(1) = x;
    error("Test_polyopvar_subsasgn:noerror","f(1) = x must error.")
catch ME
    if strcmp(ME.identifier,"Test_polyopvar_subsasgn:noerror") || ~contains(ME.message,'subsasgn')
        error("Unexpected result of f(1) = x: %s",ME.message)
    end
end
try
    y = x;      y{1,1} = x;
    error("Test_polyopvar_subsasgn:noerror","f{1,1} = x must error.")
catch ME
    if strcmp(ME.identifier,"Test_polyopvar_subsasgn:noerror") || ~contains(ME.message,'subsasgn')
        error("Unexpected result of f{1,1} = x: %s",ME.message)
    end
end

%%%% Dependent properties still cannot be assigned
try
    y = x;      y.matdim = [1,1];
    error("Test_polyopvar_subsasgn:noerror","Assigning the dependent property matdim must error.")
catch ME
    if strcmp(ME.identifier,"Test_polyopvar_subsasgn:noerror")
        rethrow(ME)
    end
end

disp('Test_polyopvar_subsasgn passed: field assignment sets chkval=false.')
