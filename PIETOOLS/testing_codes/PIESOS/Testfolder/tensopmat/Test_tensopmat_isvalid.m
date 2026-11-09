% Test_tensopmat_isvalid.m -- verify @tensopmat/isvalid.m.
% Checks (1) a properly built object is reported valid, and (2) several
% corrupted copies are reported INVALID via the [logval,errmsg] return
% contract, i.e. isvalid never throws (callers such as common_vars.m/
% plus.m/mtimes.m use '~A.chkval && ~isvalid(A)' guards to probe objects).

% DJ (CC), 10/15/2026: Initial coding;

p = randi([1,3]);  d = randi([1,3]);
N = randi([2,3]);  M = randi([2,3]);
dom = [0,1];

[Cop,~] = rand_tensopmat(p,d,N,M,dom,0.7);
Cop.chkval = false;

% A properly constructed object must be valid
try
    [logval,errmsg] = isvalid(Cop);
catch ME
    error("Test_tensopmat_isvalid failed: isvalid threw on a valid object: %s",ME.message)
end
if ~logval
    error("Test_tensopmat_isvalid failed: properly built object reported invalid (%s).",errmsg)
end

check_invalid = @(Cbad,desc) local_check_invalid(Cbad,desc);

% Corruption 1: 'ops' not a cell
Cbad = Cop;     Cbad.ops = 5;
check_invalid(Cbad,'ops not a cell');

% Corruption 2: dims.in wrong shape (row instead of column) -- only a
% genuine corruption when d>1 (a 1x1 vector is its own transpose).
if d>1
    Cbad = Cop;     Cbad.dims.in = Cop.dims.in';
    check_invalid(Cbad,'dims.in wrong shape');
end

% Corruption 3: dims.in value mismatched with the block's column size
% (must return logval=false, not throw)
if size(Cop.ops,2)>=1
    Cbad = Cop;     Cbad.dims.in(1) = Cbad.dims.in(1)+1;
    check_invalid(Cbad,'dims.in(1) mismatched with block size');
end

% Corruption 4: dims.out value mismatched with the actual block's row size
if size(Cop.ops,1)>=1
    Cbad = Cop;     Cbad.dims.out(1) = Cbad.dims.out(1)+1;
    check_invalid(Cbad,'dims.out(1) mismatched with block size');
end

% Corruption 5: depmat.in wrong size
Cbad = Cop;     Cbad.depmat.in = [Cop.depmat.in, Cop.depmat.in];
check_invalid(Cbad,'depmat.in wrong size');

% Corruption 6: block of wrong type
Cbad = Cop;
idx_nonempty = find(~cellfun(@isempty,Cop.ops),1);
Cbad.ops{idx_nonempty} = 5;
check_invalid(Cbad,'block of wrong type');

disp('Test_tensopmat_isvalid passed.')

function local_check_invalid(Cbad,desc)
Cbad.chkval = false;
try
    [logval,errmsg] = isvalid(Cbad);
catch ME
    error("Test_tensopmat_isvalid failed: isvalid threw on corruption '%s': %s",desc,ME.message)
end
if logval
    error("Test_tensopmat_isvalid failed: corruption '%s' not detected (reported valid).",desc)
end
if ~(ischar(errmsg) || isstring(errmsg)) || strlength(string(errmsg))==0
    error("Test_tensopmat_isvalid failed: corruption '%s' returned empty/invalid errmsg.",desc)
end
end
