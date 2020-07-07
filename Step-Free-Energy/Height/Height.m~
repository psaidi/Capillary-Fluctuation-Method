clear
comp=90;
temp=1570;
orientation=110;

savedir=['/Users/peymansaidi/Desktop/StepFE-Comp/crystallayers/RESAULTS/' int2str(comp) '-' int2str(temp) ...
     '-' int2str(orientation) '/'];

height=importdata([savedir 'Height' int2str(comp) '_'  int2str(orientation) '.txt']);

[total,col]=size(height);

for i=1:total
    average(i,1)=mean(height(i,2:end));
end

deviation=zeros(total,1);
for i=1:total
    for j=2:col
   deviation(i,1)=deviation(i,1)+abs(average(i,1)-height(i,j));
    end
    deviation(i,1)=deviation(i,1)/(col-1);
end
hold on

plot(deviation)

mean(deviation)