
clear

n=5; % order of descriptor plant
nu=4; % rank of its E matrix
mf=2; % number of fault inputs
md=2; % number of disturbance inputs
mu=2; % number of controls
p=3; % number of measurements
pz=1; % number of controlled outputs

E=[eye(4),zeros(4,1); zeros(1,5)];

A=[2 -3 -2 -2 2;
3 3 2 -3 -1;
-1 0 2 -4 -2;
-4 1 -3 -4 -1;
-4 3 1 -1 0]-3*[eye(4) zeros(4,1); zeros(1,5)]; % E(5,5)=A(5,5)=0 hence the plant is impulsive

Bd=[0.5 -0.4;
0.5 0.2;
-0.3 -0.5;
-0.5 -1;
-0.1 -0.2];

Bu=[2 1; 0 -1; -1 0; 0 -2; 1 0];

Cz=[1, 0, -2, -1, 2];
Dzf=[0,1];
Dzd=[0,0];
Dzu=[3, 0];

C=[2 -2 0 1 1;
-2 -1 -2 1 1;
0 -2 2 1 2];
Dd=[0 0; 
0 0;
0 0]; 
Du=[-2 1; 1 0 ; 2 -1];

% Definition of matrices Bf and Bd so that s0 is an invariant zero of subsystem (A,Bf,C,Df)
s0=1; % required invariant zero (unstable)
pom=[A-s0*eye(n); C];
pomK=[0 1.1; -1.4 0; 0.8 0.2; -0.2 -0.2; 0 0]; % auxiliary matrix (arbitrary)
pom=pom*pomK;
Bf=pom(1:n,:);
Df=pom(n+1:n+p,:);

% Preliminary control u=K0*y+u1 so that the resulting system is a state-space system
% with balanced realization denoted by
% (tA,tBf,tBd,tBu,tCz,tDzf,tDzd,tDzu,tC,tDf,tDd,tDu)

K0=[ 0.2, -1, -2; 0, 0.5, -0.5]; % arbitrary matrix such that matrix A22 is nonsingular
pom=A+Bu*K0*C;
A11=pom(1:4,1:4);
A12=pom(1:4,5);
A21=pom(5,1:4);
A22=pom(5,5);

pom=Bf+Bu*K0*Df;
B1f=pom(1:4,:);
B2f=pom(5,:);

pom=Bd+Bu*K0*Dd;
B1d=pom(1:4,:);
B2d=pom(5,:);

B1u=Bu(1:4,:);
B2u=Bu(5,:);

pom=Cz+Dzu*K0*C;
Cz1=pom(:,1:4);
Cz2=pom(:,5);
bDzf=Dzf+Dzu*K0*Df;
bDzd=Dzd+Dzu*K0*Dd;

C1=C(:,1:4);
C2=C(:,5);

tA=A11-A12*A22^(-1)*A21;
tBf=B1f-A12*A22^(-1)*B2f;
tBd=B1d-A12*A22^(-1)*B2d;
tBu=B1u-A12*A22^(-1)*B2u;

tCz=Cz1-Cz2*A22^(-1)*A21;
tDzf=bDzf-Cz2*A22^(-1)*B2f;
tDzd=bDzd-Cz2*A22^(-1)*B2d;
tDzu=Dzu-Cz2*A22^(-1)*B2u;

tC=C1-C2*A22^(-1)*A21;
tDf=Df-C2*A22^(-1)*B2f;
tDd=Dd-C2*A22^(-1)*B2d;
tDu=-C2*A22^(-1)*B2u;

sys=ss(tA,[tBf,tBd,tBu],[tCz;tC],[tDzf,tDzd,tDzu;tDf,tDd,tDu]);
sysb=balreal(sys);
[AA,BB,CC,DD]=ssdata(sysb);
tA=AA;
tBf=BB(:,1:2);
tBd=BB(:,3:4);
tBu=BB(:,5:6);
tCz=CC(1,:);
tC=CC(2:4,:)
tDzf=DD(1,1:2);
tDzd=DD(1,3:4);
tDzu=DD(1,5:6);
tDf=DD(2:4,1:2);
tDd=DD(2:4,3:4);
tDu=DD(2:4,5:6);

% Definition of transfer matrices of the plant subsystems

tGzf=ss(tA,tBf,tCz,tDzf);
tGzd=ss(tA,tBd,tCz,tDzd);
tGzu=ss(tA,tBu,tCz,tDzu);
tzero(tGzu) % tGzu has no zeros in C_+ (it is a wide 1\times 2 matrix) 

tGf=ss(tA,tBf,tC,tDf);
tzero(tGf) % tGf has a zero in s0=1
tGd=ss(tA,tBd,tC,tDd);
tGu=ss(tA,tBu,tC,tDu);

% FE filter/observer
[HH,Y0]=FE(tA,tBf,tBd,tBu,tCz,tDzf,tDzd,tDzu,tC,tDf,tDd,tDu)

% transfer matrix of filter F
F=ss(tA-HH*tC,HH,-Y0*tC,Y0);
norm(eye(2)-minreal(F*tGf),'inf') % computed 3.0039 = \alpha_F
norm(minreal(F*tGd),'inf') % computed 0.6611 = \beta_F

% Nominal controller
CoK=nom_contr(tA,tBf,tBd,tBu,tCz,tDzf,tDzd,tDzu,tC,tDf,tDd,tDu)

PP1.A=tA;
PP1.B1=tBd;
PP1.B2=tBu;
PP1.C1=tCz;
PP1.D11=tDzd;
PP1.D12=tDzu;
PP1.C2=tC;
PP1.D21=tDd;
PP1.D22=tDu;

PP2.A=tA;
PP2.B1=tBf;
PP2.B2=tBu;
PP2.C1=tCz;
PP2.D11=tDzf;
PP2.D12=tDzu;
PP2.C2=tC;
PP2.D21=tDf;
PP2.D22=tDu;

sysd=ss(PP1.A, [PP1.B1, PP1.B2], [PP1.C1; PP1.C2], [PP1.D11 PP1.D12; PP1.D21 PP1.D22]);
sysf=ss(PP2.A, [PP2.B1, PP2.B2], [PP2.C1; PP2.C2], [PP2.D11 PP2.D12; PP2.D21 PP2.D22]);

% checking
Tzd=lft(sysd,CoK);
pole(Tzd)
norm(Tzd,'inf') % computed  0.0075 = \beta

Tzf=lft(sysf,CoK);
pole(Tzf)
norm(Tzf,'inf') % computed 1.4968 = \alpha
% %%%%%%%%%%%%%%%

% FFF compensator
[AK,BK,CK,DK]=FFF_comp(tA,tBf,tBd,tBu,tCz,tDzf,tDzd,tDzu,tC,tDf,tDd,tDu)

[muu,muu1]=size(AK);
BK1=BK(:,1:mf);
BK2=BK(:,mf+1:mf+pz);
CK1=CK(1:p,:);
CK2=CK(p+1:p+mu,:);
DK11=DK(1:p,1:mf);
DK12=DK(1:p,mf+1:mf+pz);
DK21=DK(p+1:p+mu,1:mf);
DK22=DK(p+1:p+mu,mf+1:mf+pz);

CK=[CK1; CK2];
DK1=[DK11; DK21];
DK2=[DK12; DK22];

% checking
pomX=ss(AK,BK1,CK2,DK21);
pomG=ss(AK,BK1,CK1,DK11);
norm(tGzf+tGzu*pomX,'inf')
norm(tGf+tGu*pomX-pomG,'inf') % should be zero
% %%%%%%%%%%%%%%%%%%%%


% Numerical simulation of the CLS tracking
% including simulation of FE

dt=0.0001; % time increment
N=1000000; 
om=5;

% definition of time vector t, and disturbance input vector d 
magd=1; %1;
for i=1:1.1*N
    t(i)=i*dt;
    d(:,i)=magd*[sin(om*t(i)); cos(om*t(i))];
end

ef=N/500;
% definition of fault input vector f
magf=5; %5;
for i=1:N/2;
    f(:,i)=magf*[0;0];
end
for i=N/2+1:N/2+ef
f(:,i)=magf*(i-N/2)/(ef)*[-1,2];
end
for i=N/2+ef+1:1.1*N
    f(:,i)=magf*[-1; 2];
end

% definition of to-be-tracked signal r(scalar)
magr=10; %10;
for i=1:N/4
    r(i)=magr*i/(N/4);
end
for i=N/4+1:3*N/4
    r(i)=magr;
end
for i=3*N/4+1:1.1*N
    r(i)=magr*(N-i)/(N/4);
end

w=[d;f]; % join inputs d and f

Dzw=[tDzd,tDzf];
Dw=[tDd,tDf];
Bw=[tBd,tBf];

[ACK,BCK,CCK,DCK]=ssdata(CoK); % controller matrices
L=(eye(2)-DCK*tDu)^(-1);

% Intial values
x0=[-3,2,-2,1]';
u(:,1)=[0;0]; % control
x(:,1)=x0; % plant state
hx(:,1)=zeros(4,1);  % observer/FE estimator state
xi(:,1)=zeros(3,1); % state of compensator
eta(:,1)=zeros(2,1); % state of controller

eA=expm(tA*dt);
eAK=expm(AK*dt);
eACK=expm(ACK*dt);
pt=N/75; % N/500, for passive, N/75 for active
for i=1:1.1*N
    sign=0;
    if i>N/2+pt
        sign=1;
    end
    
    % plant
    x(:,i+1)=eA*x(:,i)+tA^(-1)*(eA-eye(4))*(Bw*w(:,i)+tBu*u(:,i));
    z(i)=tCz*x(:,i)+Dzw*w(:,i)+tDzu*u(:,i);
    y(:,i)=tC*x(:,i)+Dw*w(:,i)+tDu*u(:,i);

    % FE observer
    hx(:,i+1)=eA*hx(:,i)+tA^(-1)*(eA-eye(4))*(tBu*u(:,i)+HH*(y(:,i)-tC*hx(:,i)-tDu*u(:,i)));
    hf(:,i)=Y0*(y(:,i)-tC*hx(:,i)-tDu*u(:,i));
    
    % FD filter
    th=y(:,i)-tC*hx(:,i)-tDu*u(:,i);
    nth(i)=norm(th);
    tht(i)=10.0613;
    
    % FFF compensator
    xi(:,i+1)=eAK*xi(:,i)+AK^(-1)*(eAK-eye(3))*(BK1*sign*hf(:,i)+BK2*r(i));
    h(:,i)=CK*xi(:,i)+DK1*sign*hf(:,i)+DK2*r(i);
    
    % control law
    eta(:,i+1)=eACK*eta(:,i)+ACK^(-1)*(eACK-eye(2))*( BCK*y(:,i)+[-BCK,zeros(2,2)]*h(:,i) );
    %u(:,i+1)=CCK*eta(:,i)+DCK*y(:,i)+[-DCK,eye(2)]*h(:,i);
    u(:,i+1)=L*(CCK*eta(:,i)+DCK*tC*x(:,i)+DCK*Dw*w(:,i)+[-DCK,eye(2)]*h(:,i));
    
end

plot(t,z,t,r) % controlled output and to-be-tracked signal
%plot(t,hf(1,:),t,hf(2,:)) % fault estimates
%plot(t,nth,t,tht) % fault detection signal and threshold level


function [HH,Y0]=FE(tA,tBf,tBd,tBu,tCz,tDzf,tDzd,tDzu,tC,tDf,tDd,tDu)
% finding FE filter

md=2;
mf=2;
p=3;
betF=3; % min 0.66
Dw=[tDd,tDf];
Bw=[tBd,tBf];
Jg=[eye(md) zeros(md,mf); zeros(mf,md), -betF^2*eye(mf)];
[XX,LL,cL]=care(tA',tC',Bw*Jg*Bw',Dw*Jg*Dw',Bw*Jg*Dw');
cL=cL';
tAL=tA-cL*tC;
tBfL=tBf-cL*tDf;

% checking the eig
eig(tAL)
% %%%%%%%%%%%%

[V,D]=eig(Dw*Jg*Dw');
[VVV,DDD]=ordmax(V,D);

f1=sqrt(DDD(1,1));
f2=sqrt(-DDD(2,2));
f3=sqrt(-DDD(3,3));
fd=diag([f1,f2,f3]);
cM=VVV*fd;

% checking
J=diag([1,-1,-1]);
Dw*Jg*Dw'-cM*J*cM'
% %%%%%%%%%%%%%%

cN=[zeros(mf,p-mf),eye(mf)]*cM^(-1);
sst=0; % the poles of thefilter are required to be in \Re(s)<-sst
[XXX,LL,K1]=care((tAL+sst*eye(4))',(cN*tC)',1*eye(4),1*eye(2));
K1=-K1';

% checking the eig
eig(tAL+K1*cN*tC)
% %%%%%%%%%%%%

HH=cL-K1*cN; % matrix of observer-based FE filter
%FNf=ss(tAL+K1*cN*tC,tBfL+K1*cN*tDf,cN*tC,cN*tDf);
Y=cN*tDf-cN*tC*(tAL+K1*cN*tC)^(-1)*(tBfL+K1*cN*tDf);
Y=Y^(-1); % matrix of observer-based FE filter
Y0=Y*cN;
end

function CoK=nom_contr(tA,tBf,tBd,tBu,tCz,tDzf,tDzd,tDzu,tC,tDf,tDd,tDu)
% Finding nominal controller

aa=0.01; % weighting coefficient
sysj=ss(tA,[[tBd,aa*tBf],tBu],[tCz;tC], [[tDzd,aa*tDzf],tDzu; [tDd,aa*tDf],tDu]);

[CoK,CL,GAM] = hinfsyn(sysj,3,2);
CoK=reduce(CoK,2); % reducing to 2 the controller order

end

function [AK,BK,CK,DK]=FFF_comp(tA,tBf,tBd,tBu,tCz,tDzf,tDzd,tDzu,tC,tDf,tDd,tDu)
% construction of FFF compensator

mf=2;
md=2;
p=3;
pz=1;
mu=2;

bD=[2 2.5]; % arbitrary
pom1=[tDzu; bD]^(-1);
U=pom1(:,1:pz);
bU=pom1(:,pz+1:mu);
es=0; % % the poles of the compensator are required to be in \Re(s)< -es
[Xpom,bC]=icare(tA-tBu*U*tCz+es*eye(4),tBu*bU,1*eye(4),1);

% checking
eig(tA-tBu*U*tCz-tBu*bU*bC)
% %%%%%%%%%%%%%%%%%%%%

Vf=[0 0];
V=0;
pom=minreal(dss([tA,tBu; tCz,tDzu; bC, bD],[tBf,zeros(4,pz); tDzf, -eye(pz); Vf V], [tC,tDu; zeros(mu,4), eye(mu)], [tDf,zeros(p,pz); zeros(mu,mf+pz)], [eye(4),zeros(4,mu); zeros(mu,4+mu)]));
e=pole(pom) % available poles for canceling. We choose e(2)

ttA=tA-tBu*pom1*[tCz; bC];
pom2=null((e(2)*eye(4)-ttA)');
pom2=pom2';
VfV=(pom2*tBu*bU)^(-1)*pom2*[tBf-tBu*U*tDzf,tBu*U];
pom3=balreal(minreal(dss([tA,tBu; tCz,tDzu; bC, bD],[tBf,zeros(4,pz); tDzf, -eye(pz); VfV], [tC,tDu; zeros(mu,4), eye(mu)], [tDf,zeros(p,pz); zeros(mu,mf+pz)], [eye(4),zeros(4,mu); zeros(mu,4+mu)])));

% FFF compensator
[AK,BK,CK,DK]=ssdata(pom3); 

end

function HH1=FD(tA,tBf,tBd,tBu,tCz,tDzf,tDzd,tDzu,tC,tDf,tDd,tDu)
% finding FD filter

md=2;
mf=2;
p=3;
betF=0.66; % min 0.66
Dw=[tDd,tDf];
Bw=[tBd,tBf];
Jg=[eye(md) zeros(md,mf); zeros(mf,md), -betF^2*eye(mf)];
[XX,LL,cL]=care(tA',tC',Bw*Jg*Bw',Dw*Jg*Dw',Bw*Jg*Dw');
cL=cL';
tAL=tA-cL*tC;
tBfL=tBf-cL*tDf;

% checking the eig
eig(tAL)
% %%%%%%%%%%%%

[V,D]=eig(Dw*Jg*Dw');
[VVV,DDD]=ordmax(V,D);

f1=sqrt(DDD(1,1));
f2=sqrt(-DDD(2,2));
f3=sqrt(-DDD(3,3));
fd=diag([f1,f2,f3]);
cM=VVV*fd;

% checking
J=diag([1,-1,-1]);
Dw*Jg*Dw'-cM*J*cM'
% %%%%%%%%%%%%%%

cN=[zeros(mf,p-mf),eye(mf)]*cM^(-1);
sst=0; % the poles of thefilter are required to be in \Re(s)<-sst
[XXX,LL,K1]=care((tAL+sst*eye(4))',(cN*tC)',1*eye(4),1*eye(2));
K1=-K1';

% checking the eig
eig(tAL+K1*cN*tC)
% %%%%%%%%%%%%

HH1=cL-K1*cN; % matrix of observer-based FD filter
end
