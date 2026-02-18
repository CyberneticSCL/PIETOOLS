
% This script is to test that composition of 'nopvar' or 'ndopvar' objects
% in 2D produces the same result as using the 'opvar2d' and 'dopvar2d'
% composition functions.
clear
pvar s1 s1_dum s2 s2_dum
dom = [0,1;-1,1];
d = 3;
m = 3;
q = 4;
n = 3;

use_dpvar = true;

% Generate random opvar/dopvar objects
if use_dpvar
    %opts.exclude = [1,0,0,0];
    prog = lpiprogram([s1,s1_dum; s2,s2_dum],dom);
    [prog,Qop] = lpivar(prog,[0,0;0,0;0,0;m,q],d);
    [prog,Qop2] = lpivar(prog,[0,0;0,0;0,0;m,q],d);
else
    Qop = rand_opvar2d([m,q],d,dom,[s1;s2],[s1_dum;s2_dum]);
    Qop2 = rand_opvar2d([m,q],d,dom,[s1;s2],[s1_dum;s2_dum]);
end
Rop = rand_opvar2d([q,n],d,dom,[s1;s2],[s1_dum;s2_dum]);
Rop2 = rand_opvar2d([q,n],d,dom,[s1;s2],[s1_dum;s2_dum]);

%nopvar+ndopvar test 
Rop_nd = dopvar2d2ndopvar(Rop,d);
[prog,Qop3] = lpivar(prog,[0,0;0,0;0,0;q,n],d);
Qop_nd3 = dopvar2d2ndopvar(Qop3,d+1);
QRop =  Qop3 + Rop;
QRop2_nd   = Rop_nd + Qop_nd3;
QRop_nd_exact = dopvar2d2ndopvar(QRop,d);
if QRop_nd_exact == QRop2_nd
    fprintf('plus nopvar test passed\n')
else
    fprintf('plus nopvar test failed\n')
end


%horzcat nopvar
Rop_nd = dopvar2d2ndopvar(Rop,d);
Rop_nd2 = dopvar2d2ndopvar(Rop2,d+1);
RRop = [Rop Rop2];
RRop_nd = [Rop_nd, Rop_nd2];
RRop_nd_exact = dopvar2d2ndopvar(RRop,d);
if RRop_nd_exact == RRop_nd
    fprintf('horzcat nopvar test passed\n')
else
    fprintf('horzcat nopvar test failed\n')
end

%horzcat nopvar ndopvar
RRop = [Rop Qop3];
RRop_nd = [Rop_nd, Qop_nd3];
RRop_nd_exact = dopvar2d2ndopvar(RRop,d);
if RRop_nd_exact == RRop_nd
    fprintf('horzcat [nopvar ndopvar] test passed\n')
else
    fprintf('horzcat [nopvar ndopvar] test failed\n')
end

%vertcat nopvar
RRop = [Rop; Rop2];
RRop_nd = [Rop_nd; Rop_nd2];
RRop_nd_exact = dopvar2d2ndopvar(RRop,d);
if RRop_nd_exact == RRop_nd
    fprintf('vertcat nopvar test passed\n')
else
    fprintf('vertcat nopvar test failed\n')
end


%blkdiag nopvar ndopvar
RRop = blkdiag(Rop, Rop2);
RRop_nd = blkdiag(Rop_nd, Rop_nd2);
RRop_nd_exact = dopvar2d2ndopvar(RRop,d);
if RRop_nd_exact == RRop_nd
    fprintf('blkdiag nopvar test passed\n')
else
    fprintf('blkdiag nopvar test failed\n')
end


%vertcat nopvar ndopvar
RRop = [Rop; Qop3];
RRop_nd = [Rop_nd; Qop_nd3];
RRop_nd_exact = dopvar2d2ndopvar(RRop,d);
if RRop_nd_exact == RRop_nd
    fprintf('horzcat [nopvar; ndopvar] test passed\n')
else
    fprintf('horzcat [nopvar; ndopvar] test failed\n')
end

%blkdiag nopvar ndopvar
RRop = blkdiag(Rop, Qop3);
RRop_nd = blkdiag(Rop_nd, Qop_nd3);
RRop_nd_exact = dopvar2d2ndopvar(RRop,d);
if RRop_nd_exact == RRop_nd
    fprintf('blkdiag [nopvar; ndopvar] test passed\n')
else
    fprintf('blkdiag [nopvar; ndopvar] test failed\n')
end

%blkdiag ndopvar nopvar
RRop = blkdiag(Qop3, Rop);
RRop_nd = blkdiag(Qop_nd3, Rop_nd);
RRop_nd_exact = dopvar2d2ndopvar(RRop,d);
if RRop_nd_exact == RRop_nd
    fprintf('blkdiag [ndopvar; nopvar] test passed\n')
else
    fprintf('blkdiag [ndopvar; nopvar] test failed\n')
end

%subsref nopvar
RRop = Rop(1:2, 2:3); 
RRop_nd = dopvar2d2ndopvar(RRop, d+1);
Rop_nd_sliced = Rop_nd(1:2, 2:3);
if Rop_nd_sliced == RRop_nd
    fprintf('subsref nopvar test passed\n')
else
    fprintf('subsref nopvar test failed\n')
end

% subsasgn nopvar
RRop = Rop;
RRop(1:2, 2:3) = Rop2(2:3, 1:2);
RRop_nd = dopvar2d2ndopvar(RRop, d+1);
Rop_nd_subsas = Rop_nd;
Rop_nd_subsas(1:2, 2:3) = Rop_nd2(2:3, 1:2);
if Rop_nd_subsas == RRop_nd
    fprintf('subsasgn nopvar test passed\n')
else
    fprintf('subsasgn nopvar test failed\n')
end

%horzcat ndopvar
Qop_nd = dopvar2d2ndopvar(Qop,d);
Qop_nd2 = dopvar2d2ndopvar(Qop2,d+1);
QQop = [Qop Qop2];
QQop_nd = [Qop_nd, Qop_nd2];
QQop_nd_exact = dopvar2d2ndopvar(QQop,d);
if QQop_nd_exact == QQop_nd
    fprintf('horzcat ndopvar test passed\n')
else
    fprintf('horzcat ndopvar test failed\n')
end

%vertcat ndopvar
QQop = [Qop; Qop2];
QQop_nd = [Qop_nd; Qop_nd2];
QQop_nd_exact = dopvar2d2ndopvar(QQop,d);
if QQop_nd_exact == QQop_nd
    fprintf('vertcat ndopvar test passed\n')
else
    fprintf('vertcat ndopvar test failed\n')
end

%blkdiag ndopvar 
QQop = blkdiag(Qop, Qop2);
QQop_nd = blkdiag(Qop_nd, Qop_nd2);
QQop_nd_exact = dopvar2d2ndopvar(QQop,d);
if QQop_nd_exact == QQop_nd
    fprintf('blkdiag ndopvar test passed\n')
else
    fprintf('blkdiag ndopvar test failed\n')
end




%subsref ndopvar
QQop = Qop(1:2, 2:3); 
QQop_nd = dopvar2d2ndopvar(QQop, d+1);
Qop_nd_sliced = Qop_nd(1:2, 2:3);
Qop_nd_sliced = clean_ndopvar(Qop_nd_sliced);
if Qop_nd_sliced == QQop_nd
    fprintf('subsref ndopvar test passed\n')
else
    fprintf('subsref ndopvar test failed\n')
end

% subsasgn ndopvar
% test subsasign ndopvar := ndopvar
QQop = Qop;
% QQop(1:2, 2:3) = Qop2(2:3, 1:2);
QQop_nd = dopvar2d2ndopvar(QQop, d+1);
Qop_nd_subsas = Qop_nd;
% for ind= 1:4
    
incr = 1:2;
incc = 2:3;
Qop_nd_subsas(incr, incc) = Qop_nd2(incr, incc);
t = true; 
for i = 1:3
    for j = 1:4
        arg1 = Qop_nd_subsas(i,j); % arg1 result of ndopvar subsasign
        % arg2 what should be produced
        if (sum(i == incr) > 0) && (sum(j == incc) > 0)
            arg2 = Qop2(i, j);
        else
            arg2 = Qop(i, j);
        end
        arg2_nd = dopvar2d2ndopvar(arg2, d);
        t = t*(arg2_nd == arg1);
    end
end

if t == 1
    fprintf('subsasgn ndopvar test1 passed\n')
else
    fprintf('subsasgn ndopvar test1 failed\n')
end
% test subsasign ndopvar := nopvar
t = true; 
Qop_nd_subsas = Qop_nd;
incr = 1:2;
incc = 2:3;
Qop_nd_subsas(incr, incc) = Rop_nd2(incr, incc);
for i = 1:3
    for j = 1:4
        arg1 = Qop_nd_subsas(i,j); % arg1 result of ndopvar subsasign
        % arg2 what should be produced
        if (sum(i == incr) > 0) && (sum(j == incc) > 0)
            arg2 = Rop2(i, j);
        else
            arg2 = Qop(i, j);
        end
        arg2_nd = dopvar2d2ndopvar(arg2, d);
        t = t*(arg2_nd == arg1);
    end
end

if t == 1
    fprintf('subsasgn ndopvar test2 passed\n')
else
    fprintf('subsasgn ndopvar test2 failed\n')
end



% test subsasign ndopvar := nopvar
t = true; 
Qop_nd_subsas = Rop_nd;
incr = 1:2;
incc = 2:3;
Qop_nd_subsas(incr, incc) = Qop_nd2(incr, incc);
for i = 1:4
    for j = 1:3
        arg1 = Qop_nd_subsas(i,j); % arg1 result of ndopvar subsasign
        % arg2 what should be produced
        if (sum(i == incr) > 0) && (sum(j == incc) > 0)
            arg2 = Qop2(i, j);
        else
            arg2 = Rop(i, j);
        end
        arg2_nd = dopvar2d2ndopvar(arg2, d);
        t = t*(arg2_nd == arg1);
    end
end

if t == 1
    fprintf('subsasgn ndopvar test3 passed\n')
else
    fprintf('subsasgn ndopvar test3 failed\n')
end

 

