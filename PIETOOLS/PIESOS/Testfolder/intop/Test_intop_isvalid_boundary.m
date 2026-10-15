% Test_intop_isvalid_boundary.m -- Tests that isvalid.m rejects any intmat
% row with A(v)==B(v)==0 (evaluation at the left domain endpoint) or
% A(v)==B(v)==di+1 (right endpoint), even when another variable in the
% same row spans the full domain. The expected outcome is an error from the
% constructor.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - It is an open question whether this restriction is intended, or should
%   be relaxed to match the class definition chi_{[x,y]}(z):=delta(x,z),
%   which does not exclude x==y==boundary. This test records what
%   isvalid.m does, not what it should do.

dom = [0,1];
di = 2;

% Left-boundary delta (A2=B2=0), var1 spans full domain [0,di+1]=[0,3]
intmat_left = {[0,3, 0,0]};
threw_error = false;
try
    Kop = intop({polynomial(1)},intmat_left,{'s'},dom); %#ok<NASGU>
catch
    threw_error = true;
end
if ~threw_error
    error("isvalid unexpectedly ACCEPTED a left-boundary delta (A=B=0); if this is now intended to be legal, this test (and its assumption) needs updating, not silently left failing.")
end

% Right-boundary delta (A2=B2=3=di+1), var1 spans full domain
intmat_right = {[0,3, 3,3]};
threw_error = false;
try
    Kop = intop({polynomial(1)},intmat_right,{'s'},dom); %#ok<NASGU>
catch
    threw_error = true;
end
if ~threw_error
    error("isvalid unexpectedly ACCEPTED a right-boundary delta (A=B=di+1); if this is now intended to be legal, this test (and its assumption) needs updating, not silently left failing.")
end

disp('Test_intop_isvalid_boundary passed: isvalid.m currently rejects boundary-point deltas, as expected.')
