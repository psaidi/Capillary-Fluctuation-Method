clear
savedir='/Users/peymansaidi/Desktop/StepFE-Comp/crystallayers/6-plot of stiffness/';

Amplitude=importdata([savedir 'ToTAmplitude_60_211.txt']);
%Thou=[881,483,151,136,59, 37,34,25, 20, 17]; %90 percent
%Thou=[811,424,267,246,225, 186,19,18, 16, 16]; %60 percent
%Thou=[711,534,397,257,289, 234,19,36,25,25]; %30

Thou=[681,383,151,136,59, 37,34,25, 20, 17,11,10,10]; %90 percent
%Thou=[811,424,267,246,225, 186,19,18, 16, 16,11,10,10]; %60 percent
%Thou=[711,534,397,257,289, 234,19,36,25,25,11,11, 10]; %30


%corr=6480;%110-90
%corr=30000;%112-90
%corr=1920;%110-60
%corr=9330;%112-60
%corr=920;%110-30
%corr=930;%112-30
corr=10;


%NN=21;  L=92.65; %110
NN=27;   L=133.74; %211


Deltat=50;
[NSnapshots,col] = size(Amplitude);
Average=zeros(1,col-2);

for i=3:col
Average(1,i-2)=L*mean(Amplitude(1:end,i))/(NN^2);
end

k=linspace(1,(col-2),col-2);
k(1,:)=2*pi.*k(1,:)./L;
figure
plot(log(k(1,:)),log(Average(1,:)),'--*')


%fit of a line 
logk=log(k);
logAverage=log(Average);


fit1=3;
fitend=8;

p=polyfit(logk(1,fit1:fitend),logAverage(1,fit1:fitend),1)

hold on

plot(linspace(logk(1,1),logk(1,end),100),linspace(p(1,1)*logk(1,1)+p(1,2),p(1,1)*logk(1,end)+p(1,2),100),'r')
 

AmplitudeDiff2=zeros(NSnapshots,col-2);
for i=3:col
  AmplitudeDiff2(:,i-2)=abs(Amplitude(:,i)-Average(1,i-2)).^2; 
end


SumDiff2=zeros(1,col-2);

for i=1:col-2
SumDiff2(1,i)=2*sqrt(sum(AmplitudeDiff2(:,i)))/NSnapshots;
end


FinalError=zeros(1,col-2);
for i=1:col-2
FinalError(1,i)=Thou(1,i)/Deltat*SumDiff2(1,i)/(NN);
end



figure

loglog(k,Average,'o')

figure
errorbar(k,Average,FinalError,'o')


