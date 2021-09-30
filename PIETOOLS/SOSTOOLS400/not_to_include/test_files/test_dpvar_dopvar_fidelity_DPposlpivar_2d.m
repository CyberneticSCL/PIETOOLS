% Code testing proper functioning of DPposlpivar_2dQ
% The code compares the results to those attain building the positive 
% dopvar directly using the quadratic product: Zop'*LLL*Zop
% At this point, the comparison does not work for cases with separability,
% since the positive matrix LLL would have to be adjusted.

% Note: comparison with poly version works only under certain conditions,
% and requires DPsosquadvar and DPsosposmatr to be adjusted to use
% int2coeff instead of fastint2str for dvar declaration

clear

tol = 1e-14;
ntests = 1;

% maximal degree of the monomials
dx0_max = 1;     dx1_max = 1;     dx2_max = 2;
dy0_max = 1;     dy1_max = 1;     dy2_max = 2;
d200_max = 1;    d201_max = 1;    d202_max = 2;
d210_max = 1;    d211_max = 1;    d212_max = 2;
d220_max = 2;    d221_max = 2;    d222_max = 4;

% dx0_max = 1;     dx1_max = 1;     dx2_max = 1;
% dy0_max = 1;     dy1_max = 1;     dy2_max = 1;
% d200_max = 1;    d201_max = 1;    d202_max = 1;
% d210_max = 1;    d211_max = 1;    d212_max = 1;
% d220_max = 1;    d221_max = 1;    d222_max = 1;
% maximal size of the dopvar
n0_max = 1;     nx_max = 1;     ny_max = 1;     n2_max = 2;
% options in calling the function
test_p = 0;     test_s = 0;     test_e = 0;

if test_p==0
    psatz = 0;
elseif test_p==1
    psatz = -1;
elseif test_p==2
    psatz = -2;
end

sep = zeros(1,6);
if test_s>=1
    sep(test_s) = 1;
end

excludeL = zeros(1,16);
if min(test_e)>=1
    excludeL(test_e) = 1;
end

% % % Run a number of (random) tests % % %

for testnumber = 1:ntests

if test_p>=3
    psatz = -(randi(3)-1);
end
if test_s>=6
    sep = randi(2)-1;
end
if max(test_e)>=16
    excludeL = randi(2,1,16) - 1;
    if all(excludeL)
        excludeL = zeros(1,16);
        excludeL(randi(16)) = 1;
    end
end

% % % Set up a random example % % %

I = sort(2*rand(2,2) - 1,2,'ascend');    % domain of the variables s,theta

n0 = randi(n0_max+1)-1;
nx = randi(nx_max+1)-1;
ny = randi(ny_max+1)-1;
n2 = randi(n2_max+1)-1;
n = [n0,nx,ny,n2];        % size of problem
excluden = [excludeL(1),all(excludeL(2:4)),all(excludeL(5:7)),all(excludeL(8:16))];
if all(n==0 | excluden)
    n0 = randi(n0_max);
    nx = randi(nx_max);
    ny = randi(ny_max);
    n2 = randi(n2_max);
    n = [n0,nx,ny,n2];
end

dx0 = randi(dx0_max);
dx1 = randi(dx1_max,1,3);     dx1(end) = max(dx1);
dx2 = randi(dx2_max,1,3);     dx2(end) = max(dx2);
dx = {dx0;dx1;dx2};     % degrees of monomials Zxo_s, Zxa_st, Zxb_st

dy0 = randi(dy0_max);
dy1 = randi(dy1_max,1,3);     dy1(end) = max(dy1);
dy2 = randi(dy2_max,1,3);     dy2(end) = max(dy2);
dy = {dy0,dy1,dy2};     % degrees of monomials Zyo_s, Zya_st, Zyb_st

d200 = randi(d200_max,1,3);   d201 = sort(randi(d201_max,1,5),'ascend');   d202 = sort(randi(d202_max,1,5),'ascend'); 
d210 = randi(d210_max,1,5);   d211 = sort(randi(d211_max,3,3),'ascend');   d212 = sort(randi(d212_max,3,3),'ascend');
d220 = randi(d220_max,1,5);   d221 = sort(randi(d221_max,3,3),'ascend');   d222 = sort(randi(d222_max,3,3),'ascend');
d2 = {d200,d201,d202;d210,d211,d212;d220,d221,d222};     % degrees of monomials Z2oo_ss, ..., Z2bb_sstt

options.psatz = psatz;          % toggle to use psatz or not
options.exclude = excludeL;     % binary array indicating monomials to exclude
options.sep = sep;              % toggle to set R2=R1

pvar ss1 ss2 tt1 tt2
var11 = ss1;   var12 = tt1;
var21 = ss2;   var22 = tt2;
prog = sosprogram([ss1,ss2,tt2,tt2]);


% % % Build the poslpivars % % %

tic
[~,P1] = DPposlpivar_2dQ(prog,n,I,dx,dy,d2,options);
toc

% if ~any(excludeL) && ~any(sep)
%     options2 = options;
%     options2.sep = zeros(1,5);  % poly implementation assumes sep is of size 5
%     tic
%     [~,P2] = poslpivar_2d(prog,n,I,dx,dy,d2,options2);
%     toc
% end

% % % Check correctness of the results % % %

%if all(dx1==dx2) || excludeL(4) || sep % in these cases, we can compare with poly version
%    [~,P0] = poslpivar_2d(prog,n,I,dx,options);
%
if 1 % ~psatz
% If psatz=0, we can directly construct the positive dopvar using the
% monomials and LLL
% doesn't work yet....

tic
    Zop = Zop_construct(n,I,dx,dy,d2,options);
    
    % % % Build the positive matrix variable LLL % % %  
    
    dimLLL =  size(Zop,1);
    [~,LLL] = DPsosposmatr(prog,dimLLL);
    %LLLop = array2dopvar2d(LLL,Zop.I,[0,0;0,0;0,0;dimLLL,dimLLL]);
    LLLop = dopvar2d([0,0;0,0;0,0;dimLLL,dimLLL]);
    LLLop.I = Zop.I;
    LLLop.R22{1,1} = LLL;

    
    % % % Build the positive PI variable % % %
    P0 = (Zop' * LLLop) * Zop;
toc
    
end

% % % Compare the dopvars % % %

if ~eq(P0,P1,tol)
    error('DPposlpivar_2dQ does not produce the correct object!')
end
%dif = P0-P1;    
%max(dif.R00.C,[],'all')

% if ~any(excludeL) && ~any(sep)
%     if ~eq(P1,dopvar2d(P2),tol)
%         error('DPposlpivar_2dQ does not produce the same object as poslpivar_2d!')
%     end
% end

end

disp('No errors were encountered testing DPposlpivar_2dQ!')
