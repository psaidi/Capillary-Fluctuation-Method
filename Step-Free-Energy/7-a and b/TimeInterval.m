clear
savedir='/Users/peymansaidi/Desktop/StepFE-Comp/crystallayers/7-a and b/';
a1=importdata([savedir 'a190_110.txt']);
b1=importdata([savedir 'b190_110.txt']);

[NAmplitudes,cola1] = size(a1);
%[90 110]=start from 3:  50,20,10, 7, 3,2 ,2,2,2,2
%[90 110]=Results of Thou:881,483,151,136,59, 37,34,25, 25, 25
%[60 110]=start from 3:  100,100,80, 35, 35,35 ,5,5,2,2
%[60 110]=Results of Thou:811,424,267,246,225, 186,19,18, 16, 16
%[30 110]=start from 3:  120,80,40, 40, 5,5 ,5,5,2,2
%[30 110]=Results of Thou:711,534,397,257,289, 234,19,36, 25, 25
K=4; %first K value is 3 since the first one is counter and the second is for k=0
tMax=5;
Timeinterval=10;

    for j=3:cola1
        a1Average(1,j-2)=mean(a1(:,j));
        b1Average(1,j-2)=mean(b1(:,j));
    end
    
    figure
    plot(a1Average(1,2:end),'*')
    figure
    plot(b1Average(1,2:end),'r*')
    

for j=0:NAmplitudes-1
for i=1:NAmplitudes-j
  Power(i,j+1)=a1(i,K)*a1(i+j,K)+b1(i,K)*b1(i+j,K);
end
end

for i=1:NAmplitudes
SumPower(1,i)=sum(Power(:,i));
end


average=zeros(1,NAmplitudes);

for i=0:NAmplitudes-1
average(1,i+1)=SumPower(1,i+1)/(NAmplitudes-i);
end
figure
plot (average,'g*')


Time=linspace(1,tMax*Timeinterval,tMax);
AverageMod0=average(1:tMax);
AverageMod1=AverageMod0(1,:)./AverageMod0(1,1);


AverageMod2=log(AverageMod1);
figure
plot (Time,AverageMod2,'s')
p1=(AverageMod2(1,end)-AverageMod2(1,1))/(Time(1,end))

liney=linspace(0,AverageMod2(1,end),tMax);
hold on
plot(Time,liney,'--r')

p = polyfit(Time,AverageMod2,1)
endPoint=p(1,1)*(Time(1,end))+p(1,2);

liney2=linspace(p(1,2),endPoint,tMax);
plot(Time,liney2)
AverageMode3=linspace(p(1,2),endPoint,tMax);
AverageMode4=exp(AverageMode3);
AverageMode4p=exp(liney);
figure
plot (Time,AverageMod1,'rs')
hold on 
plot (Time,AverageMode4p,'--r')
plot(Time,AverageMode4);

Thou=-2/(p1+p(1,1))



A(1,:)=1:30;
B(1,:)=average(1,1:30);


