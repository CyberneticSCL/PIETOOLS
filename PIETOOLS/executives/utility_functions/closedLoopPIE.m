function PIEcl = closedLoopPIE(PIE, K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% closedLoopPIE.m     PIETOOLS 2021a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a PIE of the form
% T \dot x +Tw w + Tu u= A x + B1 w + B2 u
% z = C1 x + D11 w + D12 u
% y = C2 x + D21 w + D22 u
% where the PI operators are stored as fields of PIE object
% this function returns the closed loop PIE object that corresponds to the 
% PIE 
% 
% (T + Tu*K) \dot x +Tw w = (A + B2*K) x + B1 w 
% z = (C1 + D12*K)x + D11 w 
% y = (C2 + D22*K) x + D21 w
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE: a structure with fields T Tw Tu A B1 B2 C1 C2 D11 D12 D21 D22
% K : an opvar object with dimension equal to [B1.dim(:,2) T.dim(:,1)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIEcl: a structure with fields Tcl Tw Acl B1 C1cl D11 C2cl D21

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding SS - 5/30/2021

if nargin~=2
    error('Closed loop PIE construction requires a PIE structure and controllers gains K as inputs.')
end

if ~isa(K, 'opvar')
    error('Input K must be an opvar class object');
end

fieldList = fieldnames(PIE);

requiredList = {'T'; 'Tw'; 'Tu'; 'A'; 'B1'; 'B2'; 'C1'; 'D11'; 'D12'; 'C2'; 'D21'; 'D22'};

if ~ismember(fieldList, requiredList)
    error('PIE object must have fields T Tw Tu A B1 B2 C1 C2 D11 D12 D21 D22');
end

if ~isa(PIE.T,'opvar')
    error('PIE.T must be an opvar class object');
end
if ~isa(PIE.Tw,'opvar')
    error('PIE.Tw must be an opvar class object');
end
if ~isa(PIE.Tu,'opvar')
    error('PIE.Tu must be an opvar class object');
end
if ~isa(PIE.A,'opvar')
    error('PIE.A must be an opvar class object');
end
if ~isa(PIE.B1,'opvar')
    error('PIE.B1 must be an opvar class object');
end
if ~isa(PIE.B2,'opvar')
    error('PIE.B2 must be an opvar class object');
end
if ~isa(PIE.C1,'opvar')
    error('PIE.C1 must be an opvar class object');
end
if ~isa(PIE.C2,'opvar')
    error('PIE.C2 must be an opvar class object');
end
if ~isa(PIE.D11,'opvar')
    error('PIE.D11 must be an opvar class object');
end
if ~isa(PIE.D12,'opvar')
    error('PIE.D12 must be an opvar class object');
end
if ~isa(PIE.D21,'opvar')
    error('PIE.D21 must be an opvar class object');
end
if ~isa(PIE.D22,'opvar')
    error('PIE.D22 must be an opvar class object');
end

PIEcl.T = PIE.T+PIE.Tu*K;
PIEcl.Tw = PIE.Tw;
PIEcl.A = PIE.A+PIE.B2*K; 
PIEcl.B1 = PIE.B1;
PIEcl.C1 = PIE.C1+PIE.D12*K;
PIEcl.D11 = PIE.D11;
PIEcl.C2 = PIE.C2+PIE.D22*K;
PIEcl.D21 = PIE.D21;

opvar dummy;
dummy.I = PIE.T.I; dummy.dim = [PIE.T.dim(:,1),[0;0]];
dummy.var1 = PIE.T.var1; dummy.var2 = PIE.T.var2;

PIEcl.Tu = dummy;
PIEcl.B2 = dummy;
dummy.dim = [PIE.C1.dim(:,1),[0;0]];
PIEcl.D12 = dummy;
dummy.dim = [PIE.C2.dim(:,1),[0;0]];
PIEcl.D22 = dummy;

end