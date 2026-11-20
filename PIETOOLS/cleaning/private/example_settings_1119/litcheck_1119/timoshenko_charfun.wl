(* Characteristic functions of library case 16 (one equation, alp = bet = 1, gam = 2) and of the
   clamped-free Timoshenko beam with all parameters 1 (two equations, as in case 14 and the restored
   case 15). Run 11/19/2026; outputs in PHASE_C.md section 10.
   w = Exp[lam t] W[s]:  W'''' - 2 lam^2 W'' + (lam^4 + lam^2) W = 0,  Y = {W, W', W'', W'''}. *)
Mmat[mu_] := {{0,1,0,0},{0,0,1,0},{0,0,0,1},{-(mu^2+mu),0,2 mu,0}};
(* coded BCs: W(0) = 0, W'(0) = 0, W''(1) = W(1), W'''(1) = W'(1) *)
Fcoded[lam_] := Module[{mu = lam^2, E1}, E1 = MatrixExp[Mmat[mu]];
  Det[Join[{{1,0,0,0},{0,1,0,0}}, {{-1,0,1,0},{0,-1,0,1}}.E1]]];
(* BCs of command_line_input_Examples/stability_16_3.m: W(0) = 0, W'''(0) + (1 - lam^2) W'(0) = 0,
   W''(1) = lam^2 W(1), W'''(1) = 2 lam^2 W'(1) *)
Fcorr[lam_] := Module[{mu = lam^2, E1}, E1 = MatrixExp[Mmat[mu]];
  Det[Join[{{1,0,0,0},{0,1-mu,0,1}}, {{-mu,0,1,0},{0,-2 mu,0,1}}.E1]]];
(* two-equation beam: lam^2 W = W'' - Phi', lam^2 Phi = Phi'' + W' - Phi; Y = {W, W', Phi, Phi'};
   BCs W(0) = Phi(0) = 0, Phi'(1) = 0, W'(1) - Phi(1) = 0 *)
Mt[lam_] := {{0,1,0,0},{lam^2,0,0,1},{0,0,0,1},{0,-1,lam^2+1,0}};
Ftimo[lam_] := Module[{E1 = MatrixExp[Mt[lam]]},
  Det[Join[{{1,0,0,0},{0,0,1,0}}, {{0,0,0,1},{0,1,-1,0}}.E1]]];
(* number of zeros in the box a <= Re lam <= b, |Im lam| <= c (argument principle, n points per side) *)
count[F_, a_, b_, c_, n_:4000] := Module[{pts, vals},
  pts = Join[Table[a - I c + (b - a) t, {t, 0, 1, 1/n}], Table[b - I c + 2 I c t, {t, 0, 1, 1/n}],
             Table[b + I c - (b - a) t, {t, 0, 1, 1/n}], Table[a + I c - 2 I c t, {t, 0, 1, 1/n}]];
  vals = F /@ N[pts, 30];
  Round[Total[Arg[Rest[vals]/Most[vals]]]/(2 Pi)]];
(* zeros i*omega, 0 < omega <= 12 (F(i omega) is real) *)
roots[F_] := Module[{g, grid, vals, br},
  g[w_?NumericQ] := Re[F[I w]];
  grid = Range[1/100, 12, 1/200]; vals = g /@ N[grid, 30];
  br = Pick[Most[grid], Sign[Most[vals]] Sign[Rest[vals]], -1];
  Table[w /. FindRoot[g[w], {w, b, b + 1/200}, WorkingPrecision -> 30, Method -> "Brent"], {b, br}]];
{count[Fcoded, 1/1000, 6, 12], count[Fcorr, 1/1000, 6, 25, 8000], count[Fcorr, 1/1000, 50, 50, 8000],
 count[Ftimo, 1/1000, 6, 25, 8000], Precision[Fcorr[N[50 + 50 I, 30]]]}
(* -> {8, 0, 0, 0, 27.07} *)
{N[roots[Fcorr], 8], N[roots[Ftimo], 8]}
(* -> both {1.2739649, 1.8939628, 4.6441092, 4.8304462, 7.8136751, 7.9254182, 10.966909, 11.046721} *)
