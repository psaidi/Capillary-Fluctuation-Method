clear
TotDump=77;
for DumpNo=1:TotDump
savedir='/Users/peymansaidi/Desktop/StepFE-Comp/crystallayers/1-dumptrack/';
savedir1='/Users/peymansaidi/Desktop/StepFE-Comp/crystallayers/2-crystallayers/';
positionfirst=importdata([savedir 'ModDump' int2str(DumpNo) '.txt']);
domain=importdata([savedir 'ModDomain' int2str(DumpNo) '.txt']);

[natoms,col] = size(positionfirst);
Xlo=domain(1,1);
Xhi=domain(1,2);
XSlope=domain(1,3);
Ylo=domain(2,1);
Yhi=domain(2,2);
YSlope=domain(2,3);
Zlo=domain(3,1);
Zhi=domain(3,2);
ZSlope=domain(3,3);


a=5.4627538;
limitinitial=0.7;
InitialPoint=Ylo+limitinitial;
Ref1=[1,1,1];
Ldis=a*sqrt(3)/3;








Final=zeros(natoms,5);
Final(:,1)=positionfirst(:,1);
Final(:,2)=positionfirst(:,3);
Final(:,3)=positionfirst(:,4);
Final(:,4)=positionfirst(:,5);
Final(:,5)=floor((positionfirst(:,4)-InitialPoint)/Ldis)+1;





%line1='ITEM: TIMESTEP';
%line2='0';
%line3='ITEM: NUMBER OF ATOMS';
%[N1,col2] = size(positionM);
%line4=num2str(N1);
%line5='ITEM: BOX BOUNDS pp pp pp';
%line6=[num2str(min(positionM(:,3))) ' ' num2str(min(positionM(:,3))+Lx)];
%line7=[num2str(Ylo) ' ' num2str(Yhi)];
%line8=[num2str(Zlo) ' ' num2str(Zhi)];
%line9='ITEM: ATOMS id type x y z ix iy iz';

fid = fopen([savedir1 'IDTypePositionLayer' int2str(DumpNo) '.txt'], 'w');
%fprintf(fid,'%s\n', line1);
%fprintf(fid,'%s\n', line2);
%fprintf(fid,'%s\n', line3);
%fprintf(fid,'%s\n', line4);
%fprintf(fid,'%s\n', line5);
%fprintf(fid,'%s\n', line6);
%fprintf(fid,'%s\n', line7);
%fprintf(fid,'%s\n', line8);
%fprintf(fid,'%s \n', line9);



for i=1:natoms			 
    fprintf(fid,'%d %d %f %f %f %d \n', positionfirst(i,1),positionfirst(i,2),Final(i,2),Final(i,3),Final(i,4),Final(i,5));
end	
fclose(fid);

fid = fopen([savedir1 'Domain' int2str(DumpNo) '.txt'], 'w');
line1=[num2str(Xlo) ' ' num2str(Xhi) ' ' num2str(XSlope)];
line2=[num2str(Ylo) ' ' num2str(Yhi) ' ' num2str(YSlope)];
line3=[num2str(Zlo) ' ' num2str(Zhi) ' ' num2str(ZSlope)];

fprintf(fid,'%s\n', line1);
fprintf(fid,'%s\n', line2);
fprintf(fid,'%s\n', line3);

end

%%


for j=1:natoms
    p=Final(j,5);
    if (rem(p,2)==0)
       plot3(Final(j,2),Final(j,3),Final(j,4),'o','MarkerEdgeColor',[1 1 0], 'MarkerFaceColor',[1 1 0] ,'MarkerSize',6);
       grid on
       box on
       hold on
    else
    plot3(Final(j,2),Final(j,3),Final(j,4),'o','MarkerEdgeColor',[1 0 1], 'MarkerFaceColor',[1 0 1] ,'MarkerSize',4);
       grid on
       box on
       hold on
    end
    
end
view (90,0)

