clear all;
close all;
clc;
format long;
global K;
global alpha;
global L;
global B;
global thre;
global c;
global dd0;
global v;
global f1total;
global newk1;
global CL;
global f2;
%% parameter setting
K=5;%the number of users
N=128;%number of RIS elements
B=2*10^6;%bandwidth 10MHz
delta=10^(-3)*B*10^(-174/10);%noise power in watt
c=700;%cycles/bit
% threshold=0.440:0.006:0.5;
threshold=0.464;
thres=1000*threshold; %the threshold of delay
for k=1:K
    pt(k)=rand+1;%transmit power 
    v(k)=pt(k)/delta;%ratio between transmit power and noise power
    L(k)=5*10^4*(rand+10^2);%computation task of user k 10^6-2*10^6 bits
    f2(k)=2*10^8*(rand+5);%CPU frequncy of user k 10^8-2*10^8 cycle/s
end
f1total=K*10^10;%total CPU frequency of MEC node
%% location of nodes and large scale fading
locamec=[-120,0];locaris=[-120,50];localuser=circle(0,0,10,K);%location of nodes
for k=1:K
   dd(k)=sqrt((localuser(1,k)+120)^2+(localuser(2,k))^2);
   dg(k)=sqrt((localuser(1,k)+120)^2+(localuser(2,k)-50)^2);
   kg(k)=10^(1.3-0.003*dg(k));
end
dr=50;%distances between nodes
for k=1:K
    bd(k)=10^((-30-35*log10(dd(k)))/10);
    bg(k)=10^((-30-22*log10(dg(k)))/10);
end
br=10^((-30-22*log10(dr))/10);%large scale fading coefficients of 3 paths
kr=10^(1.3-0.003*dr);%Rician factors 
for k=1:K
K1(k)=(kg(k)+1)*(kr+1);K2(k)=kg(k)+kr+1;%auxiliary expressions
end
%% generating LOS components
for k=1:K
dd0(k)=bd(k)+N*bg(k)*br*K2(k)/K1(k); %influence from small scale fading
end
hg_ba=zeros(N,K);hr_ba=zeros(N,1);
for k=1:K
for r=1:N
    hg_ba(r,k)=sqrt(kg(k)*bg(k)/(kg(k)+1))*exp(j*2*pi*sqrt(0.2)*rand(1));%LOS components of incident channels
    hr_ba(r)=sqrt(kr*br/(kr+1))*exp(j*2*pi*sqrt(0.2)*rand(1));%LOS components of reflecting channels
end
end
%% adding NLOS components
for k=1:K
    hg(:,k)=hg_ba(:,k)+sqrt(bg(k)/(kg(k)+1))*sqrt(1/2)*(randn(N,1)+j*randn(N,1)); %the incident channels
end
    hr=hr_ba+sqrt(br/(kr+1))*sqrt(1/2)*(randn(N,1)+j*randn(N,1)); %the reflecting cahnnels
for k=1:K
a=diag(hr_ba')*hg_ba(:,k);
A(:,:,k)=a*a';
end
y1=zeros(4,1);
y2=zeros(4,1);
y3=zeros(4,1);
%% improved LMM optimization
phi1L=zeros(N,N);
for m=1:N
phi1L(m,m)=exp(j*2*pi*rand(1));%generate RIS reflecting coefficients randomly
end   
for k=1:K
    f1L(k)=f1total/K;%MEC CPU frequenyc assigned to user k
end

%% initial solve offloading ratios
for k=1:K
eta1L(k)=abs(hr_ba'*phi1L*hg_ba(:,k))^2;
CL(k)=B*log2(1+v(k)*eta1L(k));
alpha1(k)=CL(k)*c*f1L(k)/(f1L(k)*f2(k)+CL(k)*c*(f1L(k)+f2(k)));%optimum offloading ratio of user k
end

uu(1)=5;
uu(2)=6;
coo=2;
while(abs(uu(coo)-uu(coo-1))>=0.00001)
%% modified edge computing resource;
options = optimset('MaxFunEvals',1e3);
 x0=[f1L*0.99 0.6];
[resultmecL,fval,exitflag]=fsolve('yexL',x0);
f1L=resultmecL(1:K);

%% solve offloading ratios
for k=1:K
alpha1(k)=CL(k)*c*f1L(k)/(f1L(k)*f2(k)+CL(k)*c*(f1L(k)+f2(k)));%optimum offloading ratio of user k
end

%% RIS reflecting matrix
for k=1:K
    ttt(k)=alpha1(k)*c/f1L(k);
end
cvx_precision best
cvx_begin sdp
% cvx_solver mosek
variable VL(N,N) complex semidefinite
variable tau nonnegative
minimize(tau);
subject to
for k=1:K
   pow_p(2,inv_pos((B*(tau-ttt(k)))/(alpha1(k)*L(k))))-1<=v(k)*trace(A(:,:,k)*VL);
end
for ii=1:N
    VL(ii,ii)==1;
end
cvx_end
%% SDR
term=1000000;
xxL=zeros(N,1);
for ite=1:5000
    [U,D,~] = svd(VL);
    r = (randn(N,1)+1j*randn(N,1))./sqrt(2);
    v_hat = U*sqrt(D)*r;
    v_hat=exp(1j*angle(v_hat));
    for k=1:K
   eta(k)=v(k)*v_hat'*A(:,:,k)*v_hat;    
   CL(k)=B*log2(1+eta(k));
   to(k)=alpha1(k)*L(k)/CL(k)+alpha1(k)*L(k)*c/f1L(k);
   tl(k)=(1-alpha1(k))*L(k)*c/f2(k);
   t(k)=max(to(k),tl(k));
    end
    if max(t)<term
        term=max(t);
        xxL=v_hat;
    end
end
phi1L=diag(xxL');

for k=1:K
   eta(k)=v(k)*xxL'*A(:,:,k)*xxL;    
   CL(k)=B*log2(1+abs(eta(k)));
   to(k)=alpha1(k)*L(k)/CL(k)+alpha1(k)*L(k)*c/f1L(k);
   tl(k)=(1-alpha1(k))*L(k)*c/f2(k);
   t(k)=max(to(k),tl(k));
end 
coo=coo+1;
uu(coo)=max(t);
end
f1Li=f1L;
phi1Li=phi1L;
alpha1i=alpha1;

%% LMM 
phi1L=zeros(N,N);
for m=1:N
phi1L(m,m)=exp(j*2*pi*rand(1));%generate RIS reflecting coefficients randomly
end   
for k=1:K
    f1L(k)=f1total/K;%MEC CPU frequenyc assigned to user k
end
%% initial solve offloading ratios
for k=1:K
eta1L(k)=abs(hr_ba'*phi1L*hg_ba(:,k))^2;
CL(k)=B*log2(1+v(k)*eta1L(k));
alpha1(k)=CL(k)*c*f1L(k)/(f1L(k)*f2(k)+CL(k)*c*(f1L(k)+f2(k)));%optimum offloading ratio of user k
end
%% initial solve offloading ratios
for k=1:K
eta1L(k)=abs(hr_ba'*phi1L*hg_ba(:,k))^2;
CL(k)=B*log2(1+v(k)*eta1L(k));
alpha1(k)=CL(k)*c*f1L(k)/(f1L(k)*f2(k)+CL(k)*c*(f1L(k)+f2(k)));%optimum offloading ratio of user k
end

uu(1)=5;
uu(2)=6;
coo=2;
while(abs(uu(coo)-uu(coo-1))>=0.00001)
%% Update MEC computation resources LMM
tau=10;
taut=9;
while(abs(tau-taut)>=0.1)
    taut=tau;
 for k=1:K
 late(k)=(L(k)*c^2*CL(k)+L(k)*c*f1L(k))/(f1L(k)*f2(k)+c*CL(k)*(f1L(k)+f2(k)));
 end
 tau=max(late);
cvx_precision best
cvx_begin 
variable f(K)
find(f);
subject to
for k=1:K
(L(k)*c^2*CL(k)+L(k)*c*f(k))-tau*(f(k)*f2(k)+c*CL(k)*(f(k)+f2(k)))<=0;
f(k)>=0;
end
sum(f)<=f1total;
cvx_end
end
f1L=f;

%% solve offloading ratios
for k=1:K
alpha1(k)=CL(k)*c*f1L(k)/(f1L(k)*f2(k)+CL(k)*c*(f1L(k)+f2(k)));%optimum offloading ratio of user k
end

%% RIS reflecting matrix
for k=1:K
    ttt(k)=alpha1(k)*c/f1L(k);
end
cvx_precision best
cvx_begin sdp
% cvx_solver mosek
variable VL(N,N) complex semidefinite
variable tau nonnegative
minimize(tau);
subject to
for k=1:K
   pow_p(2,inv_pos((B*(tau-ttt(k)))/(alpha1(k)*L(k))))-1<=v(k)*trace(A(:,:,k)*VL);
end
for ii=1:N
    VL(ii,ii)==1;
end
cvx_end
%% SDR
term=1000000;
xxL=zeros(N,1);
for ite=1:5000
    [U,D,~] = svd(VL);
    r = (randn(N,1)+1j*randn(N,1))./sqrt(2);
    v_hat = U*sqrt(D)*r;
    v_hat=exp(1j*angle(v_hat));
    for k=1:K
   eta(k)=v(k)*v_hat'*A(:,:,k)*v_hat;    
   CL(k)=B*log2(1+eta(k));
   to(k)=alpha1(k)*L(k)/CL(k)+alpha1(k)*L(k)*c/f1L(k);
   tl(k)=(1-alpha1(k))*L(k)*c/f2(k);
   t(k)=max(to(k),tl(k));
    end
    if max(t)<term
        term=max(t);
        xxL=v_hat;
    end
end
phi1L=diag(xxL');

for k=1:K
   eta(k)=v(k)*xxL'*A(:,:,k)*xxL;    
   CL(k)=B*log2(1+abs(eta(k)));
   to(k)=alpha1(k)*L(k)/CL(k)+alpha1(k)*L(k)*c/f1L(k);
   tl(k)=(1-alpha1(k))*L(k)*c/f2(k);
   t(k)=max(to(k),tl(k));
end 
coo=coo+1;
uu(coo)=max(t);
end

%% AOCMMF optimization
% record maximum outage probabilty as delay threshold changes
con=1; %index of delay threshold
for thre=threshold;
for k=1:K
    f1(k)=f1total/K;%MEC CPU frequenyc assigned to user k
end
%% solve offloading ratios
for k=1:K
alpha(k)=1-thre*f2(k)/(L(k)*c);%optimum offloading ratio of user k
z(k)=(2^(f1(k)*alpha(k)*L(k)/(B*(thre*f1(k)-alpha(k)*L(k)*c))))-1;
omega(k)=sqrt(2*z(k)/(dd0(k)*v(k)));
end

%% record maximum outage probability when nonminmax
phi1=zeros(N,N);
for m=1:N
phi1(m,m)=exp(j*2*pi*rand(1));%generate RIS reflecting coefficients randomly
end
for k=1:K
     if ((thre*f1(k)-alpha(k)*L(k)*c<=0))
         nonp(k)=1;
         continue
     end
    para1(k)=(abs(hr'*phi1*hg(:,k))^2)/dd0(k);
    nonp(k)=1-marcumq(sqrt(2*para1(k)),omega(k));%outage probability of all users
end
result2(con)=max(nonp);%the maximum outage probability

if (thre==0.464)
    y1(2,1)=max(nonp);
    y2(2,1)=min(nonp);
    y3(2,1)=mean(nonp);
end
%% record maximum outage probability when minmax
for k=1:K
neweta1(k)=abs(hr_ba'*phi1*hg_ba(:,k))^2/dd0(k); %initialization with random RIS phase shift
newk1(k)=sqrt(2*neweta1(k));
end
count=2;
newtau1(1)=1;newtau1(2)=0.95;
 while(abs(newtau1(count)-newtau1(count-1))>0.001&&count<=100)
%% Update MEC computation resources
 options = optimset('MaxFunEvals',1e3);
 x0=[f1*0.99 0.5];
[resultmec,fval,exitflag]=fsolve('yex',x0);
f1=resultmec(1:K);
for k=1:K
    qwe(k)=(2^(f1(k)*alpha(k)*L(k)/(B*(thre*f1(k)-alpha(k)*L(k)*c))))-1;
    pa2(k)=sqrt(2*qwe(k)/(dd0(k)*v(k)));
    newp2(k)=1-marcumq(newk1(k),pa2(k));
end
newtau1(count)=max(newp2);
%% RIS reflecting matrix
i=2;
p(1)=newtau1(count);
p(2)=newtau1(count)-0.0011;
for k=1:K
x(k,1)=0.01;%initial of eta(k)
end
%% CCCP procedure
while(abs(p(i)-p(i-1))>0.001&&i<=100)
cvx_precision best
cvx_begin sdp
% cvx_solver mosek
variable V(N,N) complex semidefinite
variable tau nonnegative
variable eta(K) nonnegative
minimize(tau);
subject to
for k=1:K
sump=0;
for opo=0:10
    sump=sump+eta(k)^opo*(gamma(1+opo)-igamma(1+opo,qwe(k)/(v(k)*dd0(k))))/(factorial(opo)*gamma(1+opo));%series presentetion of outage probability
end
 2*sump+tau^2+exp(2*eta(k))-((p(i)+exp(x(k,i-1)))^2+2*(p(i)+exp(x(k,i-1)))*exp(x(k,i-1))*(eta(k)-x(k,i-1))+2*(p(i)+exp(x(k,i-1)))*(tau-p(i)))<=0;%outage probability constraints
 eta(k)==trace(A(:,:,k)*V)/dd0(k);
end
for ii=1:N
    V(ii,ii)==1;
end
cvx_end
x(:,i)=eta;
i=i+1;
p(i)=tau;
end
for k=1:K
para3(k)=trace(A(:,:,k)*V)/dd0(k);
iop(k)=1-marcumq(abs(sqrt(para3(k)*2)),pa2(k));%outage probability of all users
end
%% SDR
term=1.1;
xx=zeros(N,1);
for ite=1:5000
    [U,D,~] = svd(V);
    r = (randn(N,1)+1j*randn(N,1))./sqrt(2);
    v_hat = U*sqrt(D)*r;
    v_hat=exp(1j*angle(v_hat));
    for k=1:K
    neweta(k)=v_hat'*A(:,:,k)*v_hat/dd0(k);
    newk(k)=abs(sqrt(2*neweta(k)));
    newtau(k)=1-marcumq(newk(k),pa2(k));
    end
    if max(newtau)<term
        term=max(newtau);
        xx=v_hat;
    end
end
phi=diag(xx');

for k=1:K
neweta1(k)=abs(hr_ba'*phi*hg_ba(:,k))^2/dd0(k);
newk1(k)=sqrt(2*neweta1(k));
newp1(k)=1-marcumq(newk1(k),pa2(k));%computing outage probabilities after updating RIS matrix
end
count=count+1;
newtau1(count)=max(newp1);
 end
if(thre==0.464)
y1(1,1)=max(newp1);
y2(1,1)=min(newp1);
y3(1,1)=mean(newp1);
end
 result(con)=newtau1(count);

%% compute modifiedLMM
for k=1:K
z1(k)=(2^(f1Li(k)*alpha(k)*L(k)/(B*(thre*f1Li(k)-alpha(k)*L(k)*c))))-1;
omega1(k)=sqrt(2*z1(k)/(dd0(k)*v(k)));
para11(k)=(abs(hr'*phi1Li*hg(:,k))^2)/dd0(k);
if(z1(k)>0)
poutLm(k)=1-marcumq(sqrt(2*para11(k)),omega1(k));%outage probability of all users
else
    poutLm(k)=1;
end  
end
poutlmm(con)=max(poutLm);
if (thre==0.464)
    y1(3,1)=max(poutLm);
    y2(3,1)=min(poutLm);
    y3(3,1)=mean(poutLm);
end

for k=1:K
z1(k)=(2^(f1L(k)*alpha1(k)*L(k)/(B*(thre*f1L(k)-alpha1(k)*L(k)*c))))-1;
omega1(k)=sqrt(2*z1(k)/(dd0(k)*v(k)));
para11(k)=(abs(hr'*phi1L*hg(:,k))^2)/dd0(k);
if(z1(k)>0&&to(k)<thre)
poutLma(k)=1-marcumq(sqrt(2*para11(k)),omega1(k));%outage probability of all users
else
    poutLma(k)=1;
end  
end
poutlmma(con)=max(poutLma);
 con=con+1;
if (thre==0.464)
    y1(4,1)=max(poutLma);
    y2(4,1)=min(poutLma);
    y3(4,1)=mean(poutLma);
end
end
y=[y1,y2,y3];

%% plot
figure
plot(thres,result,'-r*',thres,result2,'-bo',thres,poutlmma,'-gx',thres,poutlmm,'-mp','LineWidth',1.5);
grid on;
legend('AOCMMF','EECR','LMM','Improved-LMM');
xlabel('Latency Threshold /ms');
ylabel('Maximum COOP');
hold on

figure
bar(y);
set(gca,'XTickLabel',{'AOCMMF','EECR','Improved-LMMF','LMMF in [31]'});
grid on;
legend('maximum','minimum','average');
ylabel('COOP');
hold on


