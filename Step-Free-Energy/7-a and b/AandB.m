
 
%this file extract the acurate position of interface

clear
a=5.46;

CriticalValue=0.1;

delX=0.8*a;
delZ=0.5*a;



savedir='/Users/peymansaidi/Desktop/StepFE-Comp/crystallayers/5-Boundary/';



height1=importdata([savedir 'Height30_110.txt']);
fid1 = fopen([savedir 'a130_110.txt'], 'w');
fid2 = fopen([savedir 'b130_110.txt'], 'w');


height=height1(:,2:end);
[TotDump,col] = size(height);


 a1=zeros(TotDump,col);
 b1=zeros(TotDump,col);







for i=1:TotDump
for k=0:col-1
    for n=0:col-1
        
    a1(i,k+1)=a1(i,k+1)+(height(i,n+1)*cos(2*n*pi*k/(col)));
    b1(i,k+1)=b1(i,k+1)+(height(i,n+1)*sin(2*n*pi*k/(col)));
    
    end
    
    AmpSQR1(i,k+1)=a1(i,k+1)^2+b1(i,k+1)^2;
    AAA=AmpSQR1(i,k+1);
    
    
end




fprintf(fid1,'%d %g %g %g %g %g %g %g %g %g %g %g \n', i...
    ,a1(i,1),a1(i,2),a1(i,3)...
    ,a1(i,4),a1(i,5),a1(i,6)...
    ,a1(i,7),a1(i,8),a1(i,9)...
    ,a1(i,10),a1(i,11));

fprintf(fid2,'%d %g %g %g %g %g %g %g %g %g %g %g \n', i...
    ,b1(i,1),b1(i,2),b1(i,3)...
    ,b1(i,4),b1(i,5),b1(i,6)...
    ,b1(i,7),b1(i,8),b1(i,9)...
    ,b1(i,10),b1(i,11));
end


%%

%figure
% plot(AmpSQR1(DumpNo,:))
 
% figure
 
% plot(AmpMatlab1(DumpNo,2:end), 'r')

fclose(fid1);
fclose(fid2);


for i=3:11
Average(1,i-2)=mean(AmpSQR1(1:end,i));
%err(1,i-2)=std(Amplitude(1:end,i));
end

k=linspace(1,(11-2),11-2);
figure
plot(log(k(1,:)),log(Average(1,:)),'--*')





 
