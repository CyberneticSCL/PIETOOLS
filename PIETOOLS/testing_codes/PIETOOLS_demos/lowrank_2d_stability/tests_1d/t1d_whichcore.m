function W = t1d_whichcore()                                                % CC, 09/22/2026
% t1d_whichcore -- where the seven shared core names resolve FROM INSIDE this
% directory.  Returns a 7x2 cellstr {name, resolved path}.
%
% WHY A FUNCTION IN tests_1d RATHER THAN A `which` IN THE TEST SCRIPT.  The
% shipped core lives in lowrank_2d_stability/private, which is visible only to
% functions IN lowrank_2d_stability -- not to this directory, and not to a
% script run from anywhere else.  `which('bm_setup')` therefore answers a
% different question depending on where it is asked from, and the question
% that matters is the one asked here: which file executes when t1d_bm calls
% bm_setup.  That must be this directory's copy, which T0 then proves is
% byte-identical to the shipped one.  Asking from the wrong scope is how a
% suite convinces itself it is testing shipped code when it is not.
nm = t1d_corefiles();
W = cell(numel(nm),2);
for k = 1:numel(nm)
    W{k,1} = nm{k};
    W{k,2} = which(nm{k});
end
end
