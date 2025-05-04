
clear
addpath('C:\Users\jovans\Desktop\hanso2_0');
addpath('C:\Users\jovans\Desktop\HIFOO3.501');

n=5;
nu=4;
mf=2;
md=2;
p=3;
pz=1;
mu=2;

E=[eye(4),zeros(4,1); zeros(1,5)];
A=[2 -3 -2 -2 2;
3 3 2 -3 -1;
-1 0 2 -4 -2;
-4 1 -3 -4 -1;
-4 3 1 -1 0]-3*[eye(4) zeros(4,1); zeros(1,5)]; % E(5,5)=A(5,5)=0 sledi infinite zero

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

s0=1; % Nula na Gf
pom=[A-s0*eye(n); C];
pomK=[0 1.1; -1.4 0; 0.8 0.2; -0.2 -0.2; 0 0]; % definiranje Bf i Bd takvi da se dobie nestabilna nula s0 na Gf
pom=pom*pomK;
Bf=pom(1:n,:);
Df=pom(n+1:n+p,:);


%K0=[-0.5 0.15 0.2; 21 -7.3 -8.6]; 
K0=[ 0.2, -1, -2; 0, 0.5, -0.5]; 
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

%PP.A=tA;
%PP.B=tBu;
%PP.C=tC;
%[K,val]=hifoo(PP,'s');
%CoK=ss(K.a,K.b,K.c,K.d);

% Gzu treba da nema nuli vo C+
%Gzf=dss(A,Bf,Cz,Dzf,E);
%Gzd=dss(A,Bd,Cz,Dzd,E);
%Gzu=dss(A,Bu,Cz,Dzu,E);
%tzero(Gzu)

% Gf treba da ima nula vo s0=1
%Gf=dss(A,Bf,C,Df,E);
%tzero(Gf)
%Gd=dss(A,Bd,C,Dd,E);
%Gu=dss(A,Bu,C,Du,E);

% tGzu treba da nema nuli vo C+
tGzf=ss(tA,tBf,tCz,tDzf);
tGzd=ss(tA,tBd,tCz,tDzd);
tGzu=ss(tA,tBu,tCz,tDzu);
tzero(tGzu)

% tGf treba da ima nula vo s0=1
tGf=ss(tA,tBf,tC,tDf);
tzero(tGf)
tGd=ss(tA,tBd,tC,tDd);
tGu=ss(tA,tBu,tC,tDu);


% ponatamu, FE filter
betF=3; % min 0.66
Dw=[tDd,tDf];
Bw=[tBd,tBf];
Jg=[eye(md) zeros(md,mf); zeros(mf,md), -betF^2*eye(mf)];
[XX,LL,cL]=care(tA',tC',Bw*Jg*Bw',Dw*Jg*Dw',Bw*Jg*Dw');
cL=cL';
tAL=tA-cL*tC;
tBfL=tBf-cL*tDf;

%proverka
eig(tAL)
% %%%%%%%%%%%%

[V,D]=eig(Dw*Jg*Dw');
[VVV,DDD]=ordmax(V,D);

f1=sqrt(DDD(1,1));
f2=sqrt(-DDD(2,2));
f3=sqrt(-DDD(3,3));
fd=diag([f1,f2,f3]);
cM=VVV*fd;

% proverka
J=diag([1,-1,-1]);
Dw*Jg*Dw'-cM*J*cM'
% %%%

cN=[zeros(mf,p-mf),eye(mf)]*cM^(-1);
sst=0; % polovi na filterot vo \Re(s)<-sst
[XXX,LL,K1]=care((tAL+sst*eye(4))',(cN*tC)',1*eye(4),1*eye(2));
K1=-K1';

%proverka
eig(tAL+K1*cN*tC)
% %%%%%%%%%%%%

HH=cL-K1*cN; % matrica na observer-based FE filter
%FNf=ss(tAL+K1*cN*tC,tBfL+K1*cN*tDf,cN*tC,cN*tDf);
Y=cN*tDf-cN*tC*(tAL+K1*cN*tC)^(-1)*(tBfL+K1*cN*tDf);
Y=Y^(-1); % matrica na observer-based FE filter



% Nominal controller

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


load('sa24e.mat','CoK');

sysd=ss(PP1.A, [PP1.B1, PP1.B2], [PP1.C1; PP1.C2], [PP1.D11 PP1.D12; PP1.D21 PP1.D22]);
sysf=ss(PP2.A, [PP2.B1, PP2.B2], [PP2.C1; PP2.C2], [PP2.D11 PP2.D12; PP2.D21 PP2.D22]);
% proverka
Tzd=lft(sysd,CoK);
pole(Tzd)
norm(Tzd,'inf')
Tzf=fb(sysf,CoK);
pole(Tzf)
norm(Tzf,'inf')
% %%%%%%%%%%%%%%%

% construction of FFF compensator
bD=[2 2.5]; % arbitrary
pom1=[tDzu; bD]^(-1);
U=pom1(:,1:pz);
bU=pom1(:,pz+1:mu);
es=0;
[Xpom,bC]=icare(tA-tBu*U*tCz+es*eye(4),tBu*bU,1*eye(4),1);

% proverka
eig(tA-tBu*U*tCz-tBu*bU*bC)
% %%%%%%%%%%%%%%%%%%%%

Vf=[0 0];
V=0;
pom=minreal(dss([tA,tBu; tCz,tDzu; bC, bD],[tBf,zeros(4,pz); tDzf, -eye(pz); Vf V], [tC,tDu; zeros(mu,4), eye(mu)], [tDf,zeros(p,pz); zeros(mu,mf+pz)], [eye(4),zeros(4,mu); zeros(mu,4+mu)]));
e=pole(pom) % polovi koi mozhat da se skratat. Go birame e(2)

ttA=tA-tBu*pom1*[tCz; bC];
pom2=null((e(2)*eye(4)-ttA)');
pom2=pom2';
VfV=(pom2*tBu*bU)^(-1)*pom2*[tBf-tBu*U*tDzf,tBu*U];
pom3=balreal(minreal(dss([tA,tBu; tCz,tDzu; bC, bD],[tBf,zeros(4,pz); tDzf, -eye(pz); VfV], [tC,tDu; zeros(mu,4), eye(mu)], [tDf,zeros(p,pz); zeros(mu,mf+pz)], [eye(4),zeros(4,mu); zeros(mu,4+mu)])));

% FFF compensator
[AK,BK,CK,DK]=ssdata(pom3); 
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

% proverka
pomX=ss(AK,BK1,CK2,DK21);
pomG=ss(AK,BK1,CK1,DK11);
norm(tGzf+tGzu*pomX,'inf')
norm(tGf+tGu*pomX-pomG,'inf')
% %%%%%%%%%%%%%%%%%%%%

% simulacija na CLS
% sose simulacija na FE

dt=0.0001;
N=1000000;
om=5;

% def. na t i na d
magd=1; %1;
for i=1:1.1*N
    t(i)=i*dt;
    d(:,i)=magd*[sin(om*t(i)); cos(om*t(i))];
end

ef=N/500;
% def. na f
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

% def. na r
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

w=[d;f];

Dzw=[tDzd,tDzf];
pom=balreal(CoK);
[ACK,BCK,CCK,DCK]=ssdata(pom);
L=(eye(2)-DCK*tDu)^(-1);
% iitial values
x0=[-3,2,-2,1]';
u(:,1)=[0;0]; % control
x(:,1)=x0; % plant state
hx(:,1)=zeros(4,1);  % observer/FE estimator state
xi(:,1)=zeros(3,1); % state of compensator
eta(:,1)=zeros(2,1); % state of controller

eA=expm(tA*dt);
eAK=expm(AK*dt);
eACK=expm(ACK*dt);
pt=N/75;; % N/500, za passive, N/100 za active
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
    hf(:,i)=Y*cN*(y(:,i)-tC*hx(:,i)-tDu*u(:,i));
    
    % FFF compensator
    xi(:,i+1)=eAK*xi(:,i)+AK^(-1)*(eAK-eye(3))*(BK1*sign*hf(:,i)+BK2*r(i));
    h(:,i)=CK*xi(:,i)+DK1*sign*hf(:,i)+DK2*r(i);
    
    % control law
    eta(:,i+1)=eACK*eta(:,i)+ACK^(-1)*(eACK-eye(2))*( BCK*y(:,i)+[-BCK,zeros(2,2)]*h(:,i) );
    %u(:,i+1)=CCK*eta(:,i)+DCK*y(:,i)+[-DCK,eye(2)]*h(:,i);
    u(:,i+1)=L*(CCK*eta(:,i)+DCK*tC*x(:,i)+DCK*Dw*w(:,i)+[-DCK,eye(2)]*h(:,i));
    
 
end


plot(t,z,t,r)

% filter
F=ss(tA-HH*tC,HH,-Y*cN*tC,Y*cN);
pom=minreal(Tzf*(eye(2)-minreal(F*tGf)));


%[AA BB CC DD]=ssdata(Tzf);
%eAA=expm(AA*dt);
%xx(:,i)=[x0; 0;0]';
%for i=1:1.1*N
%    xx(:,i+1)=eAA*xx(:,i)+AA^(-1)*(eAA-eye(6))*BB*f(:,i);
%    zz(:,i)=CC*xx(:,i)+DD*f(:,i);
%end
%plot(t,zz)
