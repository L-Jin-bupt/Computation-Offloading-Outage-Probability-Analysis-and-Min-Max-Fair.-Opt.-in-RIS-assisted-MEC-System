close all;
clear all;
clc;

%% setting parameters
M=128;  %number of RIS elements
B=10^7; %bandwidth 10MHz
delta=10^-3*B*10^(-174/10); %noise power -174dBm/Hz
threshold=0.34:0.003:0.4;
thres=1000*threshold;
mont_times=5000;  %Monte Carlo number
pt=1;   %transmit power w
v=pt/delta;
% % location of nodes
localdestination=[0,0,0];localRIS=[0,30,40];
localsource=[80,90,0];

% % distance
dd=sqrt((localsource(1)-localdestination(1))^2+(localsource(2)-localdestination(2))^2+(localsource(3)-localdestination(3))^2); %distance between source and destination
dg=sqrt((localsource(1)-localRIS(1))^2+(localsource(2)-localRIS(2))^2+(localsource(3)-localRIS(3))^2);  %distance between source and RIS
dr=50; %distance between RIS and destination
    
% % large-scale fading coefficients
bd=10^((-30-35*log10(dd))/10);    %large-scale fading coefficients of direct channel
bg=10^((-30-22*log10(dg))/10);    %large-scale fading coefficients of incident channels
br=10^((-30-22*log10(dr))/10);    %large-scale fading coefficients of reflecting channels

% % Rician factors
kg=10^(1.3-0.003*dg); %Rician factors of incident channel
kr=10^(1.3-0.003*dr); %Rician factors of reflecting channel
K1=(kg+1)*(kr+1);   %Auxiliary expression
K2=kg+kr+1; %Auxiliary expression

%% computation setting
f1total=10^10;
f2=10^9;
L=5*10^6;
c=700;


%% LOS components setting
hg_ba=zeros(M,1);hr_ba=zeros(M,1);
for r=1:M
    hg_ba(r)=sqrt(kg*bg/(kg+1))*exp(j*2*pi*rand(1));%The LOS components of incident channels
end

for r=1:M
    hr_ba(r)=sqrt(kr*br/(kr+1))*exp(j*2*pi*rand(1));%The LOS components of reflecting channels
end

%% setting RIS phase shift based on LOS components
phi=zeros(M,M);
for r=1:M
    theta(r)=exp(j*(angle(hr_ba(r))-(angle(hg_ba(r)))));%design RIS reflecting coefficients based on LOS components
end
phi=diag(theta);
phi1=zeros(M,M);
for m=1:M
phi1(m,m)=exp(j*2*pi*rand(1));%generate RIS reflecting coefficients randomly
end    
%% computing coverage probability with given target rate 
mm=1; %index of target rate variation
for thre=threshold
     alpha=1-f2*thre/(L*c);
%%%%%%%%% when latency of edge computing is above threshold       
     if ((thre*f1total-alpha*L*c<=0))
         pcov_Monte(mm)=1;
         pcov_Monte1(mm)=1;
         pcov_expression(mm)=1; 
         pcov_expression1(mm)=1; 
         pcov_series(mm)=1;
         continue
     end
     
%%%%%%% Monte Carlo simulation %%%%%%%
    z=2^(f1total*alpha*L/(B*(thre*f1total-alpha*L*c)))-1;
    SNR=zeros(mont_times,1);
    SNR1=zeros(mont_times,1);
    for con=1:mont_times

        % % Combined channels of LOS and NLOS components
        hd=sqrt(bd/2)*(randn(1)+j*randn(1));
        hg=hg_ba+sqrt(bg/(kg+1))*sqrt(1/2)*(randn(M,1)+j*randn(M,1)); %the incident channels
        hr=hr_ba+sqrt(br/(kr+1))*sqrt(1/2)*(randn(M,1)+j*randn(M,1)); %the reflecting cahnnels
        
        SNR(con)=abs(hd+hr'*phi*hg); %Received SNR 
        delayt(con)=alpha*L/(B*log2(1+v*SNR(con)));
        SNR1(con)=abs(hd+hr'*phi1*hg);
    end
    delayc=alpha*L*c/f1total;
    pcov_Monte(mm)=mean(SNR<=sqrt(z/v));
    pcov_Monte1(mm)=mean(SNR1<=sqrt(z/v));
%%%%%%% Closed-form expression %%%%%%%
eta=abs(hr_ba'*phi*hg_ba)^2;
eta1=abs(hr_ba'*phi1*hg_ba)^2;
dd0=bd+M*bg*br*K2/K1;
    % % computing coverage probability
%      pcov_expression(mm)=1-marcumq(sqrt(2/dd0)*eta,sqrt(2*epsi/(rou*dd0))); 
 pcov_expression(mm)=1-marcumq(sqrt(2*eta/dd0),sqrt(2*z/(dd0*v))); 
  pcov_expression1(mm)=1-marcumq(sqrt(2*eta1/dd0),sqrt(2*z/(dd0*v))); 
 pcov_series(mm)=0;
 for k=0:10
% pcov_series(mm)=pcov_series(mm)+igamma(1+k,z/(dd0*v))*(eta/dd0)^k/(factorial(k)*gamma(1+k));
 pcov_series(mm)=pcov_series(mm)+(gamma(1+k)-igamma(1+k,z/(dd0*v)))*(eta/dd0)^k/(factorial(k)*gamma(1+k));

 end
 pcov_series(mm)=exp(-eta/dd0)*pcov_series(mm);
 
 pcov_series1(mm)=0;
 for k=0:10
% pcov_series(mm)=pcov_series(mm)+igamma(1+k,z/(dd0*v))*(eta/dd0)^k/(factorial(k)*gamma(1+k));
 pcov_series1(mm)=pcov_series1(mm)+(gamma(1+k)-igamma(1+k,z/(dd0*v)))*(eta1/dd0)^k/(factorial(k)*gamma(1+k));

 end
%  pcov_series(mm)=1-exp(-eta/dd0)*pcov_series(mm);
 pcov_series1(mm)=exp(-eta1/dd0)*pcov_series1(mm);
    mm=mm+1;
 end
%% plot
figure
plot(thres,pcov_expression,'-b','LineWidth',1.5);
hold on
plot(thres,pcov_Monte,'ro',thres,pcov_series,'gp',thres,pcov_expression1,'-c',thres,pcov_Monte1,'md',thres,pcov_series1,'kx','LineWidth',1.5);
grid on;
legend('Expression-opt','Monte-opt','Series-opt','Expression-ran','Monte-ran','Series-ran');
xlabel('Latency Threshold /ms');
ylabel('Outage Probability');





