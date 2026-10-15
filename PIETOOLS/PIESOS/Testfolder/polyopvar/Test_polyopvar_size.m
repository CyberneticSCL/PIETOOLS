% Test_polyopvar_size.m -- verify size.m against known construction/
% concatenation dimensions, in both 1- and 2-output-argument form.
%
% Oracle: the number p of stacked states.
%
% DJ (CC), 10/15/2026: Initial coding;

for trial = 1:20
    rng(6000+trial);

    pvar s1 s2
    p = randi([1,4]);
    xs = cell(p,1);
    for i=1:p
        xs{i} = polyopvar(['x',num2str(i)],[s1,s2],[0,1]);
    end
    v = vertcat(xs{:});

    if ~isequal(size(v),[p,1])
        error("Test_polyopvar_size failed: size(v) mismatch (trial=%d).",trial)
    end
    if size(v,1)~=p
        error("Test_polyopvar_size failed: size(v,1) mismatch (trial=%d).",trial)
    end
    if size(v,2)~=1
        error("Test_polyopvar_size failed: size(v,2) mismatch (trial=%d).",trial)
    end
    [nr,nc] = size(v);
    if nr~=p || nc~=1
        error("Test_polyopvar_size failed: [nr,nc]=size(v) mismatch (trial=%d).",trial)
    end
    try
        size(v,3);
        error("Test_polyopvar_size failed: size(v,3) should have errored (trial=%d).",trial)
    catch ME
        if contains(ME.message,'Test_polyopvar_size')
            rethrow(ME);
        end
    end

    % size must track through arithmetic that changes matdim(1) but not (2)
    v2 = v.*v;
    if ~isequal(size(v2),[p,1])
        error("Test_polyopvar_size failed: size after .* mismatch (trial=%d).",trial)
    end
end

disp('Test_polyopvar_size passed.')
