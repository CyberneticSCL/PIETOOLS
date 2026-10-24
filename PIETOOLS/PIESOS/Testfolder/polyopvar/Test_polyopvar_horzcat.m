% Test_polyopvar_horzcat.m -- structural check of horzcat.m (vertcat is
% covered by Test_polyopvar_vertcat.m). Matrix-valued (matdim(2)>1)
% 'polyopvar' objects are not supported, so horzcat(x,y) must error
% immediately ('Matrix-valued') rather than build an object that cannot
% be evaluated.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - Deliberate restriction, not a bug: the test fails once horzcat.m
%   supports matrix-valued states, and must then be rewritten.

for trial = 1:10
    rng(9000+trial);

    pvar s1 s2
    x = polyopvar('x',[s1,s2],[0,1]);
    y = polyopvar('y',[s1,s2],[0,1]);

    try
        h = horzcat(x,y);
        error("Test_polyopvar_horzcat failed: expected horzcat to reject matdim(2)>1 construction, but it succeeded -- update this test (the maintainer may have implemented matrix-valued support).")
    catch ME
        if ~contains(ME.message,'Matrix-valued')
            rethrow(ME);
        end
    end
end

disp('Test_polyopvar_horzcat passed (confirms horzcat.m''s fail-fast matdim(2)>1 restriction).')
