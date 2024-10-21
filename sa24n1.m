
clear
E=[eye(4),zeros(4,1); zeros(1,5)];
A=[2 -3 -2 -2 2;
3 3 2 -3 -1;
-1 0 2 -4 -2;
-4 1 -3 -4 -1;
-4 3 1 -1 0]-4*[eye(4) zeros(4,1); zeros(1,5)]; % E(5,5)=A(5,5)=0 sledi infinite zero
Bf=[ -4 -2;
-1 2;
-4 -1;
-2 3;
-2 -0];
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
Df=[1 3;
2 -2;
1 -4];
Dd=[0 0; 
0 0;
0 0]; 
Du=[-2 1; 1 0 ; 2 -1];

n=5;
nu=4;
mf=2;
md=2;
p=3;
pz=1;
mu=2;

% Gzu i Gf treba da nemaat nuli vo C+
Gzf=dss(A,Bf,Cz,Dzf,E);
Gzd=dss(A,Bd,Cz,Dzd,E);
Gzu=dss(A,Bu,Cz,Dzu,E);
tzero(Gzu)

Gf=dss(A,Bf,C,Df,E);
tzero(Gf)
Gd=dss(A,Bd,C,Dd,E);
Gu=dss(A,Bu,C,Du,E);

% sledno, kontroller so Tzd=0 (nominalen kontr.)
UUU=[0 0 -1; -1 -1 0]; % se voveduva zaradi stable inverse
pom=UUU*Gd;

R=minreal(-Gzu(:,1)^(-1)*Gzd*pom^(-1));
pole(R)
R=[R*UUU; [ 0, 0, 0]];

Knm=balreal(minreal(R*(eye(p)+Gu*R)^(-1))); % ne se koristi 

% ponatamu, FE filter
A11=A(1:nu,1:nu);
A12=A(1:nu,nu+1:n);
A21=A(nu+1:n,1:nu);
A22=A(nu+1:n,nu+1:n);
B1f=Bf(1:nu,:);
B2f=Bf(nu+1:n,:);
C1=C(:,1:nu);
C2=C(:,nu+1:n);

bDf=[1; -1; -1]; % random taka shto [A22 B2f bB2f; C2 D2f bDf] e nesing.
bB2f=-2; % random

S=[A22 B2f bB2f; C2 Df bDf];
UU=S^(-1);
U=UU(1:3,:);
bU=UU(4,:);

AA=A11-[A12,B1f]*U*[A21; C1];
CC=bU*[A21; C1];

[XX,LL,GG]=care(AA',CC',eye(nu),eye(n-nu));
GG=GG';
%proverka
eig(AA-GG*CC)
% %%%%%%%%%%%%

bB1f=GG;
bBf=[bB1f; bB2f];

cE=[E, zeros(n,p); zeros(p,n+p)];
VV=[2;1]; % random
Fpr=dss([A, Bf, bBf; C, Df, bDf], [zeros(n,p); eye(p)], [zeros(mf,n), -eye(mf), VV], zeros(2,3), cE);

% proverka
tA=A11-[A12,B1f,bB1f]*S^(-1)*[A21;C1];
pom1=ss(tA, [eye(4),-[A12,B1f,bB1f]*S^(-1)], [eye(4); -S^(-1)*[A21; C1]], [zeros(4,8); zeros(4,4), -S^(-1)]);
pom2=dss([A, Bf, bBf; C, Df, bDf], eye(8), eye(8), zeros(8,8), cE);
minreal(pom1-pom2) % nula e

pom3=dss(A,eye(5),eye(5),zeros(5,5),E);
minreal(pom2*[eye(n); zeros(3,n)]-[eye(n);zeros(3,n)]*pom3-pom2*[zeros(5,n);C]*pom3) % prov. na formuata: nula e

prob=dss([A, Bf, bBf; C, Df, bDf], [Bf; Df], [zeros(mf,n), -eye(mf), VV], zeros(2,2), cE)
minreal(prob) % eye(2)

Fpr1=[zeros(mf,n), -eye(mf), VV]*pom1*[zeros(n,p); eye(p)];

pole(Fpr1)
minreal(Fpr1*Gf-eye(mf))
% %%%%%%%%%%%%%%%%%%%%%

FGd=dss([A, Bf, bBf; C, Df, bDf], [Bd; Dd], [zeros(mf,n), eye(mf), zeros(2,1); zeros(1,7), 1], zeros(3,2), cE);
[AR,BR,CR,DR]=ssdata(FGd)
FGd1=ss(AR,BR,CR(1:2,:),DR(1:2,:));
FGd2=ss(AR,BR,CR(3,:),DR(3,:));
D0=null(DR(3,:));
Dpl=pinv(DR(3,:));
F0=ss(AR-BR*Dpl*CR(3,:),BR*D0,-Dpl*CR(3,:),D0);
% prov
%minreal(FGd2*F0)
% %%%
pom=[FGd1*F0;F0];
pom=minreal(pom,1e-7) % mozi i vaka a i kako sto sledi

AAR=AR-BR*Dpl*CR(3,:);
BBR=BR*D0;
CCR1=CR(1:2,:)-DR(1:2,:)*Dpl*CR(3,:);
CCR2=-Dpl*CR(3,:);
CCR=[CCR1;CCR2];
DDR1=DR(1:2,:)*D0;
DDR2=D0;
DDR=[DDR1;DDR2];

% ottuka (mesto hinfsyn) mojata metodata od SICON
sA=AAR';
sB=CCR'
sB1=sB(:,1:2);
sB2=sB(:,3:4);
sC=BBR';
sD=DDR';
sD1=sD(:,1:2);
sD2=sD(:,3:4);

sn=4;
sm1=2;
sm2=2;
sp=1;

gamm1=0.05465234291; %%% optimalno
J1=[eye(sm1) zeros(sm1,sm2); zeros(sm2,sm1), -gamm1^2*eye(sm2)];
ssG=(sD*J1*sD')^(-1);
AAA=sA-sB*J1*sD'*ssG*sC;
SSS=sC'*ssG*sC;
QQQ=sB*J1*sB'-sB*J1*sD'*ssG*sD*J1*sB';
Ham=[AAA', -SSS; -QQQ, -AAA];
eH=eig(Ham)

v1=null(Ham-eH(3)*eye(8));
v2=null(Ham-eH(4)*eye(8));
v3=null(Ham-eH(7)*eye(8));
v4=null(Ham-eH(8)*eye(8));
V=[v1,v2,v3,v4];
V1=V(1:4,:);
V2=V(5:8,:);
XXX=V2*V1^(-1);
XXX=real(XXX);
eX=eig(XXX); % isto e so X podolu

[X,sL,sG] = care(-sA',-sC',-sB*J1*sB',-sD*J1*sD',-sB*J1*sD');
e1=eig(X);
e2=eig(sD*J1*sD');
[vpa(e1(1)),vpa(e2)]

hA=sA-sG'*sC;
hB=sB-sG'*sD;
h=chol(-sD*J1*sD');
LR=h';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JJ=[zeros(sm1,sm2) eye(sm1);  gamm1*eye(sm2) zeros(sm2,sm1)];
cA=-hA';
cB=sC'*(LR^(-1))';
cC=-JJ'*hB';
cD=JJ'*sD'*(LR^(-1))';
JJJ=[eye(sm2) zeros(sm2,sm1); zeros(sm1,sm2), -eye(sm1)];

% proverka
Th=ss(cA,cB,cC,cD);
minreal(Th'*JJJ*Th)
X1=lyap(cA',cC'*JJJ*cC);
X1-X
cB'*X1+cD'*JJJ*cC
cD'*JJJ*cD
% %%%%%%%%%%%%%%%%%%

bD=null(cD'*JJJ);
pom=bD'*JJJ*bD;
[VV,DD]=eig(pom);
pom1=[sqrtm(-DD(1:2,1:2)) zeros(2,1); zeros(1,2) eye(1)];
VV=VV*pom1;
V1=VV(:,1:2);
V2=VV(:,3);
V=[V2,V1];
U=V';

% prov
(U')^(-1)*pom*U^(-1)
%%% 

ccD=bD*U^(-1);
%E11=[cD(1:4,1:2) ccD(1:4,1:2)]; % ne treba
%E21=[cD(5:6,1:2) ccD(5:6,1:2)]; % ne treba
E12=ccD(1:2,2:3);
E22=ccD(3:4,2:3);
cC1=cC(1:2,:);
cC2=cC(3:4,:);

XX=dss(X*cA+cC1'*E12*E22^(-1)*cC2-cC2'*cC2, cC1'*E12*E22^(-1)-cC2',-(cC1-E12*E22^(-1)*cC2), E12*E22^(-1), X);

% proverka
T1=ss(sA,sB1,sC,sD1);
T2=ss(sA,sB2,sC,sD2);
pom5=minreal(T2*XX*gamm1-T1); % inf norma 10^(-9)
norm(pom5,'inf')
% %%

[uA,uB,uC,uD,uE]=dssdata(XX);
[VVV,DDD]=eig(uE);
EE22=DDD(2:4,2:4);
pom=VVV'*uA*VVV;
AA11=pom(1,1);
AA12=pom(1,2:4);
AA21=pom(2:4,1);
AA22=pom(2:4,2:4);
pom=uC*VVV;
CC1=pom(:,1);
CC2=pom(:,2:4);
pom=VVV'*uB;
BB1=pom(1,:);
BB2=pom(2:4,:);
hAA=AA22-AA21*AA11^(-1)*AA12;
hBB=BB2-AA21*AA11^(-1)*BB1;
hCC=CC2-CC1*AA11^(-1)*AA12;
hDD=uD-CC1*AA11^(-1)*BB1;
XX1=gamm1*balreal(dss(hAA,hBB,hCC,hDD,EE22));
[aa, ab, ac ad]=ssdata(XX1);

%proverka
pom6=minreal(T2*XX1-T1); % inf norma 10^(-9)
norm(pom6,'inf')
% 
XXm=minreal(XX); % ne se kratat polovite

% ponatamu se naogja RM V

ba=aa';
bb=ac';
bc=ab';
bd=ad';
tzd=ss(ba,bb,bc,bd);

G21pl=ss(AR-BR*Dpl*CR(3,:),BR*Dpl,-Dpl*CR(3,:),Dpl);
V=minreal(-(tzd-FGd1)*G21pl)

% dotuka bashka. ponatamu isto so sa24b.m

FF=[zeros(mf,n), -eye(mf), V]*pom2*[zeros(n,p); eye(p)];
[AF,BF,CF,DF]=ssdata(FF);
F=minreal(FF);
F=balreal(F);
[AF,BF,CF,DF]=ssdata(F)

% proverka
%Gf=dss(A,Bf,C,Df,E);
pom11=F*Gf;
pom12=minreal(pom11);
pole(pom12) % vekje e pom12 priblizhno nula
norm(pom12-eye(mf),'inf')
%Gd=dss(A,Bd,C,Dd,E);
pom21=F*Gd;
pom22=minreal(pom21);
pole(pom22)
norm(pom22,'inf')

% frekf. karakt. na F
do=0.01;
No=10000;
for i=1:No
    om(i)=i*do;
    pom=DF+CF*(eye(3)*complex(0,1)*om(i)-AF)^(-1)*BF;
    ff(i)=norm(pom);
end
plot(om,ff);

bC=[1 0 4 0 2];
bD=[2 2.5];
Vf=[0 0];
V=0;
pom=minreal(dss([A,Bu; Cz,Dzu; bC, bD],[Bf,zeros(n,pz); Dzf, -eye(pz); Vf V], [C,Du; zeros(mu,n), eye(mu)], [Df,zeros(p,pz); zeros(mu,mf+pz)], [E,zeros(n,mu); zeros(mu,n+mu)]));
e=pole(pom)

pom1=[Dzu; bD]^(-1);
U=pom1(:,1:pz);
bU=pom1(:,pz+1:mu);
tA=A-Bu*pom1*[Cz; bC]
pom2=null((e(3)*eye(n)*E-tA)');
pom2=pom2';
VfV=(pom2*Bu*bU)^(-1)*pom2*[Bf-Bu*U*Dzf,Bu*U];
pom3=balreal(minreal(dss([A,Bu; Cz,Dzu; bC, bD],[Bf,zeros(n,pz); Dzf, -eye(pz); VfV], [C,Du; zeros(mu,n), eye(mu)], [Df,zeros(p,pz); zeros(mu,mf+pz)], [E,zeros(n,mu); zeros(mu,n+mu)])));

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

% nov nominalen kontr.

% Preliminarno u=-Ky+u1, za da se kurtulime od impulsi, t.e. kontroler da naprai impulse-free
K=[ 0, 0, 0.3; 0, 0.1, 0]; 

K1=(eye(2)+K*Du)^(-1)*K;

% ravenkata so descr.system
AAA=A-Bu*K1*C;

A11=AAA(1:nu,1:nu);
A12=AAA(1:nu,nu+1:n);
A21=AAA(nu+1:n,1:nu);
A22=AAA(nu+1:n,nu+1:n);

BBf=Bf-Bu*K1*Df;
B1f=BBf(1:nu,:);
B2f=BBf(nu+1:n,:);

BBd=Bd-Bu*K1*Dd;
B1d=BBd(1:nu,:);
B2d=BBd(nu+1:n,:);

B1u=Bu(1:nu,:);
B2u=Bu(nu+1:n,:);

AA1=A11-A12*A22^(-1)*A21;
BBf1=B1f-A12*A22^(-1)*B2f;
BBd1=B1d-A12*A22^(-1)*B2d;
BBu1=B1u-A12*A22^(-1)*B2u;

%ravenkata so izlez y
CCC=C-Du*K1*C;
CCC1=CCC(:,1:nu);
CCC2=CCC(:,nu+1:n);

CC1=CCC1-CCC2*A22^(-1)*A21;

Daf=Df-Du*K1*Df;
Dad=Dd-Du*K1*Dd;
DDf=Daf-CCC2*A22^(-1)*B2f;
DDd=Dad-CCC2*A22^(-1)*B2d;
DDu=Du-CCC2*A22^(-1)*B2u;

% ravenkata so izlez z
CCC=Cz-Dzu*K1*C;
CCC1=CCC(:,1:nu);
CCC2=CCC(:,nu+1:n);

CCz=CCC1-CCC2*A22^(-1)*A21;

Dad=Dzd-Dzu*K1*Dd;
DDzd=Dad-CCC2*A22^(-1)*B2d;
DDzu=Dzu-CCC2*A22^(-1)*B2u;

K3=(eye(mu)+K*Du)^(-1);
sys=ss(AA1, [BBd1,BBu1*K3], [CCz; CC1], [DDzd,DDzu*K3; DDd DDu*K3]);
options=hinfsynOptions;
options.METHOD='ric';
%options.RelTol=0.01;
%options.AbsTol=0.01;
[Knm,CL,GAML,INFO] = hinfsyn(sys,p,mu,options);

% proverka
GGu=ss(AA1,BBu1*K3,CC1,DDu*K3);
RRR=minreal(Knm*(eye(p)-GGu*Knm)^(-1));
pole(RRR) % istite poloi so CL
% %%%%

Knm=balreal(Knm-K);
R=minreal(Knm*(eye(p)-Gu*Knm)^(-1));

% simulacija na FE. Treba prvo da gi obezbedam f,d,y. Za simultaneous kje trebaat
% i x1, x2. Preliminarno u=-Ky, za da se kurtulime od impulsi
dt=0.0001;
N=200000;
om=10;
magd=1;
for i=1:N
    t(i)=i*dt;
    d(:,i)=magd*[sin(om*t(i)); cos(om*t(i))];
end
magf=5;
for i=1:N/2;
    f(:,i)=magf*[0;0];
    vz(:,i)=BBf1*f(:,i)+BBd1*d(:,i);
end
for i=N/2+1:N
    f(:,i)=magf*[0.5*(t(i)-(N/2+1)*dt); 2];
    %f(:,i)=magf*[-1; 2];
    vz(:,i)=BBf1*f(:,i)+BBd1*d(:,i);
end

magr=10;
for i=1:N/4
    r(i)=magr*i/(N/4);
end
for i=N/4+1:3*N/4
    r(i)=magr;
end
for i=3*N/4+1:N
    r(i)=magr*(N-i)/(N/4);
end

plot(t,r);

eA=expm(AA1*dt);
x1(:,1)=[3; 2; -2; 1];
x2(:,1)=-3; % ne mora
for i=1:N-1
    x1(:,i+1)=eA*x1(:,i)+AA1^(-1)*(eA-eye(nu))*vz(:,i);
end    
for i=1:N
    y(:,i)=CC1*x1(:,i)+DDf*f(:,i)+DDd*d(:,i);
    x2(:,i)=-A22^(-1)*(A21*x1(:,i)+B2f*f(:,i)+B2d*d(:,i));
    u(:,i)=-K*y(:,i);
end

% sega filtriranje
Gu=dss(A,Bu,C,Du,E);
F1F2=minreal([F,F*Gu]);
[AF,BF12,CF,DF12]=ssdata(F1F2);

for i=1:N
    y1(:,i)=[y(:,i);-u(:,i)];
end

eA=expm(AF*dt);
hx(:,1)=zeros(3,1);
for i=1:N-1
    hx(:,i+1)=eA*hx(:,i)+AF^(-1)*(eA-eye(3))*BF12*y1(:,i);
end    
for i=1:N
    hf(:,i)=CF*hx(:,i)+DF12*y1(:,i);
end

plot(t,hf(1,:),t,hf(2,:))

% filtriranje na estimate of f so Butterworth
fc = 0.00011*N;
fs = N;

[b,a] = butter(4,fc/(fs/2));

hhf(1,:) = filter(b,a,hf(1,:));
hhf(2,:) = filter(b,a,hf(2,:));

plot(t,hhf(1,:),t,hhf(2,:))

% filtriranje na estimate of f so Bessel
n4=6;
[Aa,Ba,Ca,Da]=besself(n4,8);
eA=expm(Aa*dt);
hhhx1(:,1)=zeros(n4,1);
hhhx2(:,1)=zeros(n4,1);
for i=1:N-1
    hhhx1(:,i+1)=eA*hhhx1(:,i)+Aa^(-1)*(eA-eye(n4))*Ba*hf(1,i);
    hhhx2(:,i+1)=eA*hhhx2(:,i)+Aa^(-1)*(eA-eye(n4))*Ba*hf(2,i);
end  
for i=1:N
    hhhf1(i)=Ca*hhhx1(:,i)+Da*hf(1,i);
    hhhf2(i)=Ca*hhhx2(:,i)+Da*hf(2,i);
end

plot(t,hhhf1,t,hhhf2)
% %%

% crtanje na grafik na greskata
Gd=dss(A,Bd,C,Dd,E);
Tefd=minreal(-F*Gd);
pom=minreal(Tefd);
pom=balreal(pom);
[Ae,Be,Ce,De]=ssdata(pom);

eA=expm(Ae*dt);
hx(:,1)=zeros(3,1);
for i=1:N-1
    hx(:,i+1)=eA*hx(:,i)+Ae^(-1)*(eA-eye(3))*Be*d(:,i);
end    
for i=1:N
    ef(:,i)=Ce*hx(:,i)+De*d(:,i);
end

%plot(t,ef(1,:),t,ef(2,:))


% simulacija na CLS
%F=minreal(ss(Aa,Ba,Ca,Da)*F); % filtracija so Besel
Y=minreal(dss([A,Bu; Cz,Dzu; bC, bD], [zeros(n,pz); -eye(pz); V], [zeros(mu,n), eye(mu)], zeros(mu,pz), [E,zeros(n,mu); zeros(mu,n+mu)] ));
%proverka
minreal(Gzu*Y-eye(pz))
% %%%%%%%%%%%%%%%%%%%%%
w=[f;d];
Bw=[Bf,Bd];
Dw=[Df,Dd];
Dzw=[Dzf,Dzd];

Tzf=minreal(Gzf+Gzu*R*Gf);
RF=minreal(R-Y*Tzf*F);

pom=(eye(mu)+R*Gu)^(-1)*[R,Y];
pom=balreal(minreal(pom));

[Axi,Byr,Cxi,Dyr]=ssdata(pom);
By=Byr(:,1:p);
Br=Byr(:,p+1:p+pz);
tDy=Dyr(:,1:p);
tDr=Dyr(:,p+1:p+pz);

L=(eye(mu)-tDy*Du)^(-1);
L1=(eye(p)-Du*tDy)^(-1);
pom=[Axi+By*Du*L*Cxi, By*L1*C; Bu*L*Cxi, A+Bu*L*tDy*C];

[NN,NN1]=size(pom);
A11=pom(1:NN-1,1:NN-1);
A12=pom(1:NN-1,NN);
A21=pom(NN,1:NN-1);
A22=pom(NN,NN);

pom=[By*L1*Dw; Bw+Bu*L*tDy*Dw];
Bw1=pom(1:NN-1,:);
Bw2=pom(NN,:);
pom=[Br+By*Du*L*tDr; Bu*L*tDr];
Br1=pom(1:NN-1,:);
Br2=pom(NN,:);

bA=A11-A12*A22^(-1)*A21;
bBw=Bw1-A12*A22^(-1)*Bw2;
bBr=Br1-A12*A22^(-1)*Br2;

pom=[Dzu*L*Cxi, Cz+Dzu*L*tDy*C];
Cz1=pom(:,1:NN-1);
Cz2=pom(:,NN);
bCz=Cz1-Cz2*A22^(-1)*A21;
bDw=Dzw+Dzu*L*tDy*Dw-Cz2*A22^(-1)*Bw2;
bDr=Dzu*L*tDr-Cz2*A22^(-1)*Br2;

eA=expm(bA*dt);
hxxx(:,1)=[zeros(1,NN-5),-3,2,-2,1]';
for i=1:N/2
    hxxx(:,i+1)=eA*hxxx(:,i)+bA^(-1)*(eA-eye(NN-1))*(bBw*w(:,i)+bBr*r(i));
end    
for i=1:N/2
    hz(:,i)=bCz*hxxx(:,i)+bDw*w(:,i)+bDr*r(i);
end

pom=(eye(mu)+RF*Gu)^(-1)*[RF,Y];
pom=balreal(minreal(pom));

[Axi,Byr,Cxi,Dyr]=ssdata(pom);
By=Byr(:,1:p);
Br=Byr(:,p+1:p+pz);
tDy=Dyr(:,1:p);
tDr=Dyr(:,p+1:p+pz);

L=(eye(mu)-tDy*Du)^(-1);
L1=(eye(p)-Du*tDy)^(-1);
pom=[Axi+By*Du*L*Cxi, By*L1*C; Bu*L*Cxi, A+Bu*L*tDy*C];

[NN,NN1]=size(pom);
A11=pom(1:NN-1,1:NN-1);
A12=pom(1:NN-1,NN);
A21=pom(NN,1:NN-1);
A22=pom(NN,NN);

pom=[By*L1*Dw; Bw+Bu*L*tDy*Dw];
Bw1=pom(1:NN-1,:);
Bw2=pom(NN,:);
pom=[Br+By*Du*L*tDr; Bu*L*tDr];
Br1=pom(1:NN-1,:);
Br2=pom(NN,:);

bA=A11-A12*A22^(-1)*A21;
bBw=Bw1-A12*A22^(-1)*Bw2;
bBr=Br1-A12*A22^(-1)*Br2;

pom=[Dzu*L*Cxi, Cz+Dzu*L*tDy*C];
Cz1=pom(:,1:NN-1);
Cz2=pom(:,NN);
bCz=Cz1-Cz2*A22^(-1)*A21;
bDw=Dzw+Dzu*L*tDy*Dw-Cz2*A22^(-1)*Bw2;
bDr=Dzu*L*tDr-Cz2*A22^(-1)*Br2;

eA=expm(bA*dt);
%hxxxx(:,N/2+1)=[zeros(19,1); hxxx(9:12,N/2+1)];
hxxxx(:,N/2+1)=zeros(19,1); % za 'ric', 19, za 'lmi' 23, za Besel 24 (So Besel golem overshoot za t=10 
for i=N/2+1:N
    hxxxx(:,i+1)=eA*hxxxx(:,i)+bA^(-1)*(eA-eye(NN-1))*(bBw*w(:,i)+bBr*r(i));
end    
for i=N/2+1:N
    hz(:,i)=bCz*hxxxx(:,i)+bDw*w(:,i)+bDr*r(i);
end

plot(t,hz(1,:),t,r)

Tzd=minreal(Gzd+Gzu*R*Gd);
bet=norm(Tzd,'inf')
Tzf=minreal(Gzf+Gzu*R*Gf);
al=norm(Tzf,'inf')