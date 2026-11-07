% t3_galerkin.m -- validate opgram against operators with KNOWN induced norms
% before it is used to judge any solver output.  The Volterra constant 2/pi is
% the sharp test: non-trivial, exact, and the Galerkin estimate must converge
% to it FROM BELOW as the basis degree D grows.
cuadmm_path;
pvar s th
one = polynomial(1);  zero = polynomial(0);

fprintf('T3 --- identity: G must equal I exactly at every D\n');
for D = [4 8 16 32]
    G = cuadmm_private('opgram',one,[],[],s,th,[0,1],D);
    fprintf('T3 D=%-3d ||G-I||=%.3e  norm=%.12f\n', D, norm(G-eye(D+1)), norm(G));
end

fprintf('T3 --- Volterra  int_0^s : exact ||V|| = 2/pi = %.12f\n', 2/pi);
for D = [4 8 16 32 48]
    G = cuadmm_private('opgram',zero,one,[],s,th,[0,1],D);
    fprintf('T3 D=%-3d ||G||=%.12f  err=%+.3e\n', D, norm(G), norm(G)-2/pi);
end

fprintf('T3 --- adjoint  int_s^1 : must match 2/pi\n');
for D = [16 48]
    G = cuadmm_private('opgram',zero,[],one,s,th,[0,1],D);
    fprintf('T3 D=%-3d ||G||=%.12f  err=%+.3e\n', D, norm(G), norm(G)-2/pi);
end

fprintf('T3 --- rank one f(s)g(th)=s*1 over all th: ||K||=1/sqrt(3)=%.12f\n',1/sqrt(3));
for D = [4 12 24]
    G = cuadmm_private('opgram',zero,s,s,s,th,[0,1],D);
    fprintf('T3 D=%-3d ||G||=%.12f  err=%+.3e\n', D, norm(G), norm(G)-1/sqrt(3));
end

fprintf('T3 --- non-unit domain [0,2], Volterra scales linearly: 2*2/pi=%.12f\n',2*2/pi);
G = cuadmm_private('opgram',zero,one,[],s,th,[0,2],32);
fprintf('T3 dom[0,2] ||G||=%.12f err=%+.3e\n', norm(G), norm(G)-2*2/pi);
fprintf('T3DONE\n');
