function M = lb_meta()
% M = LB_META() returns the hand-checked metadata of the three PIETOOLS
% example libraries, one entry per example (primary executive only; the
% executives themselves are read from the toggles by lb_build_catalog).
%
% OUTPUT
%   M:  struct array with fields
%       lib     - 'PDE', 'DDE' or 'NDSDDF'
%       id      - example id: PDE switch index ('6', '2.2'), DDE/NDS/DDF
%                 block name ('DDE05', 'NDS01', 'DDF03')
%       name    - example name (PDE: example function name)
%       fn      - PDE: example function; DDE/NDSDDF: library file name
%       idxarg  - PDE: true if the example function takes the index first
%       lines   - DDE/NDSDDF: [first last] line of the block in the library
%       prose   - DDE/NDSDDF: lines with one comment level that are prose,
%                 left commented when the block is uncommented
%       param   - bisection parameter ('' if none)
%       setter  - sprintf template applied for a parameter value: a PDE
%                 params entry, or code evaluated after the block
%       doc     - documented value of the primary executive, as written
%                 in the library; 'tier=value;...' for tier-specific values
%       docdual - documented value for the dual executive ('' if none)
%       docnum  - documented critical value of the bisection (NaN if none)
%       lo, hi  - bisection bracket (lo: a value the block or notes call
%                 stable; hi: about 1.25 x the documented critical value)
%       bisect  - true for a bisection on param
%       run     - false for entries listed but not run (disabled blocks,
%                 duplicates, the random generator)
%       note    - remarks from the inspection of the block
%
% NOTES
%   - Values in 'doc' are quoted from the library comments; they are not
%     results. Notes marked (inferred) are derived here, not quoted.

P = 'PIETOOLS_PDE_Ex_';
% {id, function, idxarg, param, setter, doc, docdual, docnum, lo, hi, bisect, note}
pde = {
 '1',   [P,'Transport_Eq'],0,'','','stable (v=-1)','',NaN,NaN,NaN,0,''
 '2',   [P,'Lamare'],1,'','','stable','',NaN,NaN,NaN,0,'index 2 = 2.1'
 '2.2', [P,'Lamare'],1,'','','stable','',NaN,NaN,NaN,0,''
 '2.3', [P,'Lamare'],1,'','','stable','',NaN,NaN,NaN,0,''
 '2.4', [P,'Lamare'],1,'','','stable','',NaN,NaN,NaN,0,''
 '3',   [P,'Diagne'],0,'','','stable','',NaN,NaN,NaN,0,''
 '4',   [P,'Saba'],1,'','','stable','',NaN,NaN,NaN,0,'index 4 = 4.1'
 '4.2', [P,'Saba'],1,'','','stable','',NaN,NaN,NaN,0,''
 '4.3', [P,'Saba'],1,'','','stable','',NaN,NaN,NaN,0,''
 '4.4', [P,'Saba'],1,'','','stable','',NaN,NaN,NaN,0,''
 '4.5', [P,'Saba'],1,'','','stable','',NaN,NaN,NaN,0,''
 '5',   [P,'Reaction_Diffusion_Eq'],0,'lam','lam=%.10g;','9.8696 (stable for lam<pi^2)','',pi^2,0.75*pi^2,1.25*pi^2,1,''
 '6',   [P,'Reaction_Diffusion_Gahlawat'],0,'lam','lam=%.10g;','2.467 (stable for lam<pi^2/4)','',pi^2/4,0.75*pi^2/4,1.25*pi^2/4,1,''
 '7',   [P,'Parabolic_Eq_Gahlawat'],0,'lam','lam=%.10g; c=-0.5*s^3+1.3*s^2-1.5*s+0.7+lam;','4.66 (unstable for lam>4.66)','',4.66,0.75*4.66,1.25*4.66,1, ...
        'c is built from lam before params are evaluated: a lam override alone has no effect; the setter rebuilds c'
 '8',   [P,'Heat_Eq_with_ODE_Das'],0,'lam','lam=%.10g;','9.8696 (stable for lam<pi^2)','',pi^2,0.75*pi^2,1.25*pi^2,1,''
 '9',   [P,'Heat_Eq_with_ODE_in_BC'],0,'k','k=%.10g;','0 (stable for k<0)','',0,-1,0.25,1,'tolerance relative to the bracket width'
 '10',  [P,'Heat_Eq_with_ODE_Tang'],0,'','','stable (k=-2)','',NaN,NaN,NaN,0,''
 '11',  [P,'Reaction_Diffusion_Ahmadi'],0,'R','R=%.10g;','2.7 (stable for R<2.7; not tight)','',2.7,0.75*2.7,1.25*2.9,1, ...
        'default R=2.9 exceeds the documented 2.7; n=1 Dirichlet mode gives R*=pi^2/max(eig(Cm))=2.9306 (inferred)'
 '12',  [P,'Reaction_Diffusion_Ahmadi_2'],0,'R','R=%.10g;','21 (stable for R<21; not tight)','',21,0.75*21,1.25*30,1, ...
        'default R=30; Cm strictly lower triangular: cascade of heat equations, stable for every R>0 (inferred)'
 '13',  [P,'Euler_Bernoulli_Beam_Eq'],1,'','','stable','',NaN,NaN,NaN,0,'index 13 = 13.1'
 '13.2',[P,'Euler_Bernoulli_Beam_Eq'],1,'','','stable','',NaN,NaN,NaN,0,''
 '14',  [P,'Timoshenko_Beam_1'],0,'','','stable (TOC: PDE stable without damping)','',NaN,NaN,NaN,0,''
 '15',  [P,'Timoshenko_Beam_2'],0,'','','TOC: (unstable); PDE stable with damping c>0','',NaN,NaN,NaN,0,''
 '16',  [P,'Timoshenko_Beam_3'],0,'','','','',NaN,NaN,NaN,0,'no documented verdict'
 '17',  [P,'Wave_Eq_Boundary_Damped'],0,'','','stable for k>0 (k=1)','',NaN,NaN,NaN,0,''
 '18',  [P,'Wave_Eq_Datko_Boundary_Damped_1'],0,'','','stable (TOC: dual stable)','',NaN,NaN,NaN,0,''
 '19',  [P,'Wave_Eq_Datko_Boundary_Damped_2'],0,'','','stable','',NaN,NaN,NaN,0,''
 '20',  [P,'Transport_Eq_with_Disturbance'],0,'','','veryheavy=0.5000;heavy=0.5011;light=0.5162','veryheavy=0.5000;heavy=0.5058;light=0.5204',NaN,NaN,NaN,0,''
 '21',  [P,'Wave_Eq_Tip_Damped'],0,'','','veryheavy=0.8419;heavy=0.8419;light=0.8424','veryheavy=0.8419;heavy=0.8552;light=0.8722',NaN,NaN,NaN,0,'k=0.5'
 '22',  [P,'Heat_Eq_with_Distributed_Disturbance'],0,'','','veryheavy=0.2083;heavy=0.2083;light=0.2085','',NaN,NaN,NaN,0,''
 '23',  [P,'Parabolic_Eq_with_Disturbance'],0,'','','light=25.249;heavy=23.706;veryheavy=23.622','',NaN,NaN,NaN,0,'lam=4.6'
 '24',  [P,'Reaction_Diffusion_Eq_with_Disturbance'],0,'','','8.2136','',NaN,NaN,NaN,0,'lam=(1-1e-2)*pi^2'
 '25',  [P,'Diffusive_Eq_with_Distributed_Disturbance'],0,'','','0.8102','',NaN,NaN,NaN,0,'R=2.6'
 '26',  [P,'Diffusive_Eq_with_Distributed_Disturbance_2'],0,'','','2.1457','',NaN,NaN,NaN,0,'R=21-1e-3'
 '27',  [P,'Transport_Eq_with_Observed_Output'],0,'','','1.0004','',NaN,NaN,NaN,0,''
 '28',  [P,'Heat_Eq_with_Observed_Output'],0,'','','1.0005','',NaN,NaN,NaN,0,''
 '29',  [P,'Reaction_Diffusion_Eq_with_Controlled_Input'],0,'','','1.0000','',NaN,NaN,NaN,0,'lam=10'
 '30',  [P,'Reaction_Diffusion_Boundary_ODE_Control'],0,'','','heavy=0.6594','',NaN,NaN,NaN,0,'TOC: MOSEK numerr=2; lam=5'
 '31',  [P,'Euler_Bernoulli_Beam_with_Control'],0,'','','light=2.7407e-04','',NaN,NaN,NaN,0,'TOC: MOSEK numerr=2; c=0.1'
 '32',  [P,'Wave_Eq_Boundary_ODE_Control'],0,'','','heavy=0.5971','',NaN,NaN,NaN,0,'TOC: MOSEK numerr=2; c=1'
 '33',  [P,'2D_Transport_Eq'],0,'','','','',NaN,NaN,NaN,0,''
 '34',  [P,'2D_KISS_Model'],0,'','','stable for lam<=2*pi^2 (lam=19)','',NaN,NaN,NaN,0,''
 '35',  [P,'2D_Reaction_Diffusion_Eq'],0,'','','','',NaN,NaN,NaN,0,''
 '36',  [P,'2D_Telegraph_Eq'],0,'','','','',NaN,NaN,NaN,0,''
 '37',  [P,'2D_Parabolic_Eq'],0,'','','stable for a<=pi^2/2 (a=4)','',NaN,NaN,NaN,0,''
 '38',  [P,'2D_Wave_Eq'],0,'','','','',NaN,NaN,NaN,0,''
 '39',  [P,'2D_NS_Antonelli'],0,'','','','',NaN,NaN,NaN,0,''
 '40',  [P,'2D_Heat_Eq_with_ODE'],0,'','','','',NaN,NaN,NaN,0,''
 '41',  [P,'Heat_Eq_w_Delayed_Boundary_Input'],0,'','','','',NaN,NaN,NaN,0,''
 '42',  [P,'Coupled_ODE_Heat_Eq_w_Delay'],0,'','','','',NaN,NaN,NaN,0,''
 '43',  [P,'Heat_Eq_w_Interior_Delay'],0,'','','stable for tau<1.0347 (tau=1)','',NaN,NaN,NaN,0,''
 '44',  [P,'Wave_Eq_w_Boundary_Delay'],0,'mu','mu=%.10g;','0.5 (stable for mu<0.5)','',0.5,0.75*0.5,1.25*0.5,1, ...
        'listed among the 2-D examples but converts to a 1-D PIE (size [1 4])'
 '45',  [P,'2D_ReactionDiffusion_DDDN'],0,'','','exact=0.11101','',NaN,NaN,NaN,0,''
 '46',  [P,'2D_ReactionDiffusion_DDDD'],0,'','','exact=0.17110','',NaN,NaN,NaN,0,''
 '47',  [P,'2D_ReactionDiffusion_DNDN'],0,'','','exact=0.20668','',NaN,NaN,NaN,0,''
 '48',  [P,'2D_ReactionDiffusion_DNDN_EdgeObserve'],0,'','','','',NaN,NaN,NaN,0,''
 '49',  [P,'Transport_Eq_with_Disturbance2'],0,'','','veryheavy=0.1016;heavy=0.1017;light=0.1043','veryheavy=0.1016;heavy=0.1162;light=0.1362',NaN,NaN,NaN,0,'TOC: 0.1016 by numerical integration'
 '50',  'PIETOOLS_PDDE_Ex',0,'','','custom+psatz=0.8502','custom+psatz=0.8398',NaN,NaN,NaN,0,'TOC values used custom settings with settings_2d.ineq_opts.psatz=1 set by hand'
 '51',  [P,'Unstable_ReactionDiffusion'],0,'','','veryheavy=1.2683','',NaN,NaN,NaN,0,'closed-loop H2 norm'
 '52',  [P,'Unstable_ReactionDiffusion_Control'],0,'','','','',NaN,NaN,NaN,0,''
};

% {id, name, lines, prose, param, setter, doc, docdual, docnum, lo, hi, bisect, run, note}
T = 'DDE.tau(1)=%.10g;';
TN = 'tauN=%.10g; DDE.tau(1)=tauN/2; DDE.tau(2)=tauN;';
dde = {
 'DDE01','Distributed delay example',[47 52],52,'tau(1)',T,'2.04 (analytical range [0.2, 2.04])','',2.04,0.75*2.04,1.25*2.04,1,1,'A0 unstable: delay-stabilized; lo inside the documented range'
 'DDE02','Distributed delay, [1] p.270',[54 61],[],'tauN',TN,'2 ([1]: tauNmax<2)','',2,1.5,2.5,1,1,''
 'DDE03','Example B, [1] p.193',[64 71],[],'tau(1)',T,'1.71785 (taumin=.100174)','',1.71785,0.75*1.71785,1.25*1.71785,1,1,''
 'DDE04','Example, [1] p.192',[74 80],[],'tau(1)',T,'6.17258','',6.17258,0.75*6.17258,1.25*6.17258,1,1,''
 'DDE05','Example A, [2]',[82 90],[],'tau(1)',T,'1.5707','',1.5707,0.75*1.5707,1.25*1.5707,1,1,''
 'DDE06','Example C, [2],[3]',[92 99],[],'b=Ai{1}','DDE.Ai{1}=%.10g;','3 (bmax)','',3,2.25,3.75,1,1,''
 'DDE07','Example D, [2]',[101 109],[],'tauN',TN,'1.3722','',1.3722,1.3,1.25*1.3722,1,1,'A0 unstable: lo is the block value'
 'DDE08','Example E, [2] from [9]',[111 120],[],'tau(1)',T,'3.95 (tauNmax in [3.9,3.95], requires veryheavy)','',3.95,1.1,1.25*3.95,1,1,'block sets tauN=1.37 but uses DDE.tau(1)=1.1'
 'DDE09','Example 1 (unsourced), Hinf gain',[127 139],[],'','','25.907 (tau=1.5)','25.907 (tau=1.5)',NaN,NaN,NaN,0,1,'also documented: 2 at tau=.5'
 'DDE10','Example A.1, [3]',[141 150],[],'','','0.2364','0.2364',NaN,NaN,NaN,0,1,''
 'DDE11','Example A.2, [3]',[153 161],[],'','','2.87','2.87',NaN,NaN,NaN,0,1,''
 'DDE12','Example D derived from [2], Hinf gain',[164 174],[],'','','4.9784','4.9784',NaN,NaN,NaN,0,1,''
 'DDE13','Example 3, [4] from [8],[7]',[181 193],[],'','','1.339 (tau=2)','',NaN,NaN,NaN,0,1,'also documented: .1000 at tau=.99'
 'DDE14','Example B.2, [4] from [7]',[196 206],[],'','','0.3953','',NaN,NaN,NaN,0,1,''
 'DDE15','Example B.3, [4] from [6]',[209 218],[],'','','0.6104','',NaN,NaN,NaN,0,1,''
 'DDE16','Scalability example, [4], n=3, K=3',[220 241],241,'','','0.9630','',NaN,NaN,NaN,0,1,''
 'DDE17','Spring-mass network delay, n=3',[243 276],[],'h','h=%.10g; DDE.tau(1)=h; for i=2:n, DDE.tau(i)=h*i; end','stable at h=.08, unstable at h=1','',NaN,0.05,1.0,1,1, ...
         'block assigns B2(n+1,1)=1 to a local B2, so DDE.B2=0 (library defect); tolerance relative to the bracket width'
 'DDE18','Multiple showering people, [2], n=3',[278 318],[316 317 318],'','','0.4092 (n=3)','',NaN,NaN,NaN,0,1,''
 'DDE19','Wind tunnel, [4] from [5]',[321 335],[],'','','1.964 (degrees 2,4 and 4,8); 1.9369 (6,12)','',NaN,NaN,NaN,0,1,'block uses fields C0, C{1}, D1, D2'
 'DDE20','Example 1, [12] from [11]',[341 350],[],'','','0.2357','',NaN,NaN,NaN,0,1,''
 'DDE21','Example 2, [12] from [13]',[353 362],[],'','','2.327','',NaN,NaN,NaN,0,1,''
 'DDE22','Example 3, [12]',[365 375],[],'','','1.3501','',NaN,NaN,NaN,0,1,''
 'DDE23','Example 1, [10]',[377 387],[],'','','1.808','',NaN,NaN,NaN,0,1,''
 'DDE24','Example 2, [10]',[389 401],[],'','','','',NaN,NaN,NaN,0,1,'no documented value'
 'DDE25','Example 1, [13], H2 norm',[407 418],[],'','','3.6724','3.6724',NaN,NaN,NaN,0,1,''
 'DDE26','Example 2b, [13], H2 norm',[420 438],420:438,'','','3.8307','3.8307',NaN,NaN,NaN,0,1,'uncommented in the library file'
 'DDE_RAND','Random example generator',[444 498],[],'','','','',NaN,NaN,NaN,0,0,'not run: random system, no documented value'
};

TS = 'NDS.tau(1)=%.10g;';
nds = {
 'NDS01','Neutral, version of Han 2005 [1]',[48 57],[],'tau(1)',TS,'2.2255 (analytic max, uncited)','',2.2255,0.75*2.2255,1.25*2.2255,1,1,''
 'NDS02','Neutral, Park and Won 2000 Ex.3 [4]',[61 69],69,'','','stable at tau=20 (literature: tau<.7516)','',NaN,NaN,NaN,0,1,''
 'NDS03','Neutral, Park and Won 1999a Ex.2.2 [2]',[73 84],84,'','','stable (tau=[.3 .6])','',NaN,NaN,NaN,0,1,''
 'NDS04','Neutral, Park and Won 1999a Ex.2.3 [2]',[88 99],99,'tau(2)','NDS.tau(2)=%.10g;','0.552 (stable for any h1>0, h2<.552)','',0.552,0.75*0.552,1.25*0.552,1,1,''
 'NDS05','Neutral, Zhao 2018 [7], tau=[h 3h]',[103 112],[],'h','h=%.10g; NDS.tau=[h 3*h];','1.5804? (stability verified for h<2.04)','',1.5804,0.75*1.5804,1.25*2.04,1,1,''
 'NDS06','Neutral, Zhao 2018 [7] (?????)',[115 123],[],'','','0.7722? (max)','',NaN,NaN,NaN,0,0,'not run: code commented twice in the library (disabled)'
 'NDS07','Neutral, Zhao 2018 [7], scalar, tau=[h 2h]',[126 133],[],'h','h=%.10g; NDS.tau=[h 2*h];','0.5? (stability verified for h<.603)','',0.5,0.75*0.5,1.25*0.603,1,1,''
 'NDS08','Neutral, Olgac and Sipahi 2004',[138 154],154,'tau(1)',TS,'0.02755','',0.02755,0.75*0.02755,1.25*0.02755,1,1,''
 'NDS09','Neutral, Ramirez and Sipahi 2020',[157 176],176,'tau(1)',TS,'0.603','',0.603,0.75*0.603,1.25*0.603,1,1,''
 'NDS10','Neutral, source uncertain',[179 185],[],'','','','',NaN,NaN,NaN,0,1,'no documented value'
 'NDS11','Neutral, PI-controlled (unsourced)',[190 203],[],'','','','',NaN,NaN,NaN,0,1,'no documented value'
 'NDS12','Neutral, 2 delays, scalar, tau=[h 2h]',[206 215],[],'h','h=%.10g; NDS.tau=[h 2*h];','0.4448?','',0.4448,0.75*0.4448,1.25*0.4448,1,1,''
 'NDS_FRAG','Fragment (ndim=size(An,1))',[217 219],[],'','','','',NaN,NaN,NaN,0,0,'not run: fragment, An undefined'
 'DDF01','Pure difference equation',[239 243],[],'','','','',NaN,NaN,NaN,0,1,'no documented value'
 'DDF02','Pure difference equation (repeat)',[248 252],[],'','','','',NaN,NaN,NaN,0,0,'not run: identical to DDF01'
 'DDF03','Multiple showering people, DDF, n=2',[259 314],259:314,'','','','',NaN,NaN,NaN,0,1,'uncommented in the library file; note "min in [.32 .38]" is a leftover gam_guess comment'
};

M = struct('lib',{},'id',{},'name',{},'fn',{},'idxarg',{},'lines',{},'prose',{}, ...
           'param',{},'setter',{},'doc',{},'docdual',{},'docnum',{},'lo',{},'hi',{}, ...
           'bisect',{},'run',{},'note',{});
for k = 1:size(pde,1)
    r = pde(k,:);
    nm = erase(r{2},'PIETOOLS_PDE_Ex_');
    M(end+1) = struct('lib','PDE','id',r{1},'name',nm,'fn',r{2},'idxarg',logical(r{3}), ...
        'lines',[],'prose',[],'param',r{4},'setter',r{5},'doc',r{6},'docdual',r{7}, ...
        'docnum',r{8},'lo',r{9},'hi',r{10},'bisect',logical(r{11}),'run',true,'note',r{12}); %#ok<AGROW>
end
libs = {'DDE',dde,'examples_DDE_library_PIETOOLS.m'; 'NDSDDF',nds,'examples_NDSDDF_library_PIETOOLS.m'};
for L = 1:2
    C = libs{L,2};
    for k = 1:size(C,1)
        r = C(k,:);
        M(end+1) = struct('lib',libs{L,1},'id',r{1},'name',r{2},'fn',libs{L,3},'idxarg',false, ...
            'lines',r{3},'prose',r{4},'param',r{5},'setter',r{6},'doc',r{7},'docdual',r{8}, ...
            'docnum',r{9},'lo',r{10},'hi',r{11},'bisect',logical(r{12}),'run',logical(r{13}),'note',r{14}); %#ok<AGROW>
    end
end
end
