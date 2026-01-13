function C = mtimes2(A,B)
Aparams = A.params;  Bparams = B.params;
keyA = keys(Aparams); keyB = keys(Bparams);

% Pull values once (cells of polynomial)
valA = Aparams(keyA);
valB = Bparams(keyB);

Avars = vars(A);
Bvars = vars(B);

nA = numel(Avars);
nB = numel(Bvars);

Cvars = union(A.vars_out, B.vars_in);
nC = numel(Cvars);
maxkeyC = 4^nC;

paramsC = cell(1, maxkeyC);
touched = false(1, maxkeyC);     % faster pruning than cellfun('isempty',...)

% Precompute base-4 powers for key generation (same as index2keys)
pow4C = 4.^(0:nC-1).';           % column

% Precompute key digit tables
keyAidx = nDopvar.keys2index(keyA, nA);
keyBidx = nDopvar.keys2index(keyB, nB);

% ----- Precompute variable position maps (NO ismember/find in the inner loops)

midvars = A.vars_in;            % your "midvars"
nmid = numel(midvars);

% Map each Cvar -> position in Avars/Bvars (0 if absent)
[tfC_in_A, posC_in_A] = ismember(Cvars, Avars);
[tfC_in_B, posC_in_B] = ismember(Cvars, Bvars);

% Map each midvar -> position in Avars/Bvars/Cvars
[~, posMidA] = ismember(midvars, Avars);
[~, posMidB] = ismember(midvars, Bvars);
[~, posMidC] = ismember(midvars, Cvars);

% Cache pvar objects for midvars ONCE per mtimes call
varMid = cell(1, nmid);
dumMid = cell(1, nmid);
symMid = cell(1, nmid);
for t = 1:nmid
    nm = midvars{t};
    varMid{t} = pvar(nm);
    dumMid{t} = pvar([nm,'_dum']);
    symMid{t} = pvar([nm,'_mid']);
end

% ----- Main loops
for i = 1:numel(keyA)
    Ai    = valA{i};
    AiIdx = keyAidx(i,:);

    for j = 1:numel(keyB)
        Bj    = valB{j};
        BjIdx = keyBidx(j,:);

        % Build the base keyC digits from A and B once per (i,j)
        keyCdig = zeros(1, nC);

        if any(tfC_in_A)
            keyCdig(tfC_in_A) = AiIdx(posC_in_A(tfC_in_A));
        end
        if any(tfC_in_B)
            keyCdig(tfC_in_B) = BjIdx(posC_in_B(tfC_in_B));
        end

        % Apply your midvar logic; accumulate directly into paramsC
        for t = 1:nmid
            ia = posMidA(t);
            ib = posMidB(t);
            ic = posMidC(t);

            keyinA = AiIdx(ia);
            keyinB = BjIdx(ib);

            var    = varMid{t};
            dumvar = dumMid{t};
            midvar = symMid{t};

            % Helper: set digit, compute key, add polynomial
            % (nested function would be slightly cleaner, but keep inline for speed)
            switch keyinA
                case 0
                    % keyinB in {0,1,2}
                    keyCdig(ic) = keyinB;
                    kk = keyCdig * pow4C + 1;
                    addTerm(kk, Ai * Bj);

                case 1
                    if keyinB == 0
                        keyCdig(ic) = 1;
                        kk = keyCdig * pow4C + 1;
                        addTerm(kk, Ai * subs(Bj, var, dumvar));

                    elseif keyinB == 1
                        keyCdig(ic) = 1;
                        kk = keyCdig * pow4C + 1;
                        addTerm(kk, int(subs(Ai, dumvar, midvar) * subs(Bj, var, midvar), midvar, dumvar, var));

                    else % keyinB == 2
                        AB = subs(Ai, dumvar, midvar) * subs(Bj, var, midvar);

                        keyCdig(ic) = 1;
                        kk = keyCdig * pow4C + 1;
                        addTerm(kk, int(AB, midvar, 0, dumvar));

                        keyCdig(ic) = 2;
                        kk = keyCdig * pow4C + 1;
                        addTerm(kk, int(AB, midvar, 0, var));
                    end

                case 2
                    if keyinB == 0
                        keyCdig(ic) = 2;
                        kk = keyCdig * pow4C + 1;
                        addTerm(kk, Ai * subs(Bj, var, dumvar));

                    elseif keyinB == 1
                        AB = subs(Ai, dumvar, midvar) * subs(Bj, var, midvar);

                        keyCdig(ic) = 1;
                        kk = keyCdig * pow4C + 1;
                        addTerm(kk, int(AB, midvar, var, 1));

                        keyCdig(ic) = 2;
                        kk = keyCdig * pow4C + 1;
                        addTerm(kk, int(AB, midvar, dumvar, 1));

                    else % keyinB == 2
                        keyCdig(ic) = 2;
                        kk = keyCdig * pow4C + 1;
                        addTerm(kk, int(subs(Ai, dumvar, midvar) * subs(Bj, var, midvar), midvar, var, dumvar));
                    end

                case 3
                    % keyinB in {0,1,2}
                    keyCdig(ic) = 3;
                    kk = keyCdig * pow4C + 1;

                    if keyinB == 0
                        addTerm(kk, Ai * Bj);
                    elseif keyinB == 1
                        addTerm(kk, int(subs(Ai * Bj, var, midvar), midvar, var, 1));
                    else % keyinB == 2
                        addTerm(kk, int(subs(Ai * Bj, var, midvar), midvar, 0, var));
                    end
            end
        end
    end
end

% Prune using touched mask (cheaper than cellfun on a big cell array)
keyC = find(touched);
paramsOut = paramsC(keyC);

C = nDopvar(B.vars_in, A.vars_out, [A.dims(1), B.dims(2)], dictionary(keyC, paramsOut));

    function addTerm(k, term)
        % k is a scalar key in 1..maxkeyC
        if ~touched(k)
            paramsC{k} = term;
            touched(k) = true;
        else
            paramsC{k} = paramsC{k} + term;
        end
    end
end
