clear
TotDump=77;
StateCrit=0.15;
for DumpNo=1:TotDump
    DumpNo
%reading the position of the atoms 
savedir='/Users/peymansaidi/Desktop/StepFE-Comp/crystallayers/2-crystallayers/';
savedir1='/Users/peymansaidi/Desktop/StepFE-Comp/crystallayers/3-interfacefinder/';

AtomPos=importdata([savedir 'IDTypePositionLayer' int2str(DumpNo) '.txt']);
Domain=importdata([savedir 'Domain' int2str(DumpNo) '.txt']);

%Orientation of the crystal
OrgDir=[-1,1,1         %1 1 1
        1,1,1
        0,1,-2];


%OrgDir=[3,-1,-4         %321
%        2,2,2
%        1,-1,8];
    
%OrgDir=[5,-1,-10    %543
%        4,2,2
%        3,-1,14];
    

    
%OrgDir=[7,-1,-16    %765
%        6,2,2
%        5,-1,20];

    
%OrgDir=[9,-1,-22    %987
%        8,2,2
%        7,-1,26];
    
    
  

%Mapping the original orientation to the Global orientation
NewDir=[1,0,0
        0,1,0
        0,0,1];
   
%normalizing the original crystal orientation    
    OrgDir1=zeros(3,3);
    for i=1:3
        OrgDir1(:,i)=OrgDir(:,i)/norm(OrgDir(:,i));
    end
%finding the mapping matrix   
B=NewDir*inv(OrgDir1);
conv=inv(B);


[Natoms,col] = size(AtomPos);

%calculating the domain of the system
xlo=Domain(1,1);
xhi=Domain(1,2);
xSlope=Domain(1,3);
Lx=xhi-xlo;
ylo=Domain(2,1);
yhi=Domain(2,2);
ySlope=Domain(2,3);
Ly=yhi-ylo;
zlo=Domain(3,1);
zhi=Domain(3,2);
zSlope=Domain(3,3);
Lz=zhi-zlo;

%LXActual=max(AtomPos(:,3))-min(AtomPos(:,3));
%LYActual=max(AtomPos(:,4))-min(AtomPos(:,4));
%LZActual=max(AtomPos(:,5))-min(AtomPos(:,5));

%calculating the lattice parameter and the position of the first and second
%nearest neighbours
%V=LXActual*LYActual*LZActual;
%a=(8*V/Natoms)^(1/3);

a=5.46;


FirstNeighb=a*sqrt(3)/4;
SecNeighb=a*sqrt(2)/2;

%Difining the cut of distance, it should be increased if the dendity in
%liquid decreases
MaxCOD=SecNeighb+2.5;


%q vectores which are the standard vectors to second nearest neighbour
%atoms, the coefficient takes care of the size of the rj and it includes
%2*pi
QQ=2*pi*sqrt(2)*[1 1 0
    -1 1 0
    1 0 1
    -1 0 1
    0 1 1 
    0 -1 1];

%tilting the system to apply the priodicity
tanTET=zSlope/Lz;
AtomPos(:,4)=AtomPos(:,4)-AtomPos(:,5)*tanTET;



%periodicity of the system
AtomPostot = neighb(MaxCOD, xlo, xhi, ylo, yhi, zlo, zhi, zSlope, AtomPos);
[Natomstot,coltot] = size(AtomPostot);



%tilting the system back after applying the priodicity
AtomPos(:,4)=AtomPos(:,4)+AtomPos(:,5)*tanTET;
AtomPostot(:,4)=AtomPostot(:,4)+AtomPostot(:,5)*tanTET;




Distanceindex=zeros(1,1);
DistInd=zeros(Natoms,12);
DistInddis=zeros(Natoms,12);
DisToT=zeros(Natoms,1);
FNeighb=zeros(Natoms,1);
Rj=zeros(1,1,1);
Q=zeros(1,1,1);
Qfinal=zeros(Natoms,3,12);
ksi=zeros(Natoms,1);
KSI=zeros(Natoms,1);
counter=1;

%finding the nearest neighbours of all atoms
for i=1:Natoms
    
    
 x1=AtomPos(i,3);
 y1=AtomPos(i,4);
 z1=AtomPos(i,5);

 
 counter=1;
 e=2;
 Distance=zeros(1,1);
  

 for j=1:Natomstot
     
     x2=AtomPostot(j,3);
     y2=AtomPostot(j,4);
     z2=AtomPostot(j,5);
     

     
     delx1=x2-x1;
     dely1=y2-y1;
     delz1=z2-z1;
     
     if (delx1<=MaxCOD && dely1<=MaxCOD && delz1<=MaxCOD)
     
     DIST=sqrt(delx1^2+dely1^2+delz1^2);
    
         if (DIST<=MaxCOD && DIST~=0) %all neighbours
             Distance(1,counter)=DIST;
             Distance(2,counter)=j; 
             counter=counter+1;

             
             if (DIST<=3)  %first nearest neighbours (atoms in 3A CUD)
                 
                 FNeighb(i,1)=FNeighb(i,1)+1;
                 FNeighb(i,e)=j;
                e=e+1;              
             end

             
         end
     end
     
     
 end
 
 outsort = sortrows(Distance',1)';
 AtInGame=FNeighb(i,1)+12;
 
 DistInddis(i,1:AtInGame)=outsort(1,1:AtInGame);
 
 DistInd(i,1:AtInGame)=outsort(2,1:AtInGame);
 
 
 %finding rj(12 nearest aotms after first nearest neighbour atoms) 
  for q=FNeighb(i,1)+1:FNeighb(i,1)+12
     xsecond=AtomPostot(DistInd(i,q),3);
     ysecond=AtomPostot(DistInd(i,q),4);
     zsecond=AtomPostot(DistInd(i,q),5);
     
     delxq=xsecond-x1;
     delyq=ysecond-y1;
     delzq=zsecond-z1;
     
     Sizeq=sqrt(delxq^2+delyq^2+delzq^2);

     
     Q(i,1,q-FNeighb(i,1))=delxq/Sizeq;
     Q(i,2,q-FNeighb(i,1))=delyq/Sizeq;
     Q(i,3,q-FNeighb(i,1))=delzq/Sizeq;


    
 end
    

end

%calculating the direction of the rj in global coordination system
for i=1:Natoms
    for j=1:12
    Qfinal(i,:,j)=conv*Q(i,:,j)';

    end

end

%calculating the sum(exp(iq.r)) term
for i=1:Natoms
    for j=1:6
            Xq=QQ(j,1);
            Yq=QQ(j,2);
            Zq=QQ(j,3);
        for k=1:12
            Xr=Qfinal(i,1,k);
            Yr=Qfinal(i,2,k);
            Zr=Qfinal(i,3,k);

            
            QdotR=Xr*Xq+Yr*Yq+Zr*Zq;
            EXP=exp(complex(0,QdotR));
            ksi(i,1)=ksi(i,1)+EXP;

        end
    end
end

%calculation of the absolute value and square of that
Nq=12;
for i=1:Natoms
KSI(i,1)=abs(ksi(i,1)/Nq/6)^2;
end
%FiNneighbCriteria=zeros(Natoms,1);

%for i=1:Natoms
%    FNDeviation=abs(FNeighb(i,1)-4);
%    if (FNDeviation==0)
%    FiNneighbCriteria(i,1)=1;
%    else
%        FiNneighbCriteria(i,1)=1/FNDeviation;
%    end
        
%end

%for i=1:Natoms
%KSI(i,1)=KSI(i,1)*FiNneighbCriteria(i,1);
%end

KSIaverage=zeros(Natoms,1);

%calculation of 1/(Z+1)(sum(ksi(i)+ksi(i neighbours))

%for i=1:Natoms
%    sum=0;
%    ee=0;
%    for q=FNeighb(i,1)+1:FNeighb(i,1)+12
%        if (DistInd(i,q)<=Natoms)
%        sum=sum+KSI(DistInd(i,q),1);
%        ee=ee+1;
%        end
%    end
%    sum=sum+KSI(i,1);
%    KSIaverage(i,1)=sum/(ee+1);
    
%end

for i=1:Natoms
    sum=0;
    ee=0;
    for q=2:FNeighb(i,1)+1
        if (DistInd(i,q)<=Natoms)
        sum=sum+KSI(DistInd(i,q),1);
        ee=ee+1;
        end
    end
    sum=sum+KSI(i,1);
    KSIaverage(i,1)=sum/(ee+1);
    
end

fid = fopen([savedir1 'Final' int2str(DumpNo) '.txt'], 'w');
for i=1:Natoms	
    
    fprintf(fid,'%d %d %f %f %f %f \n', AtomPos(i,1),AtomPos(i,2),AtomPos(i,3),AtomPos(i,4),AtomPos(i,5),KSI(i,1));

end	

fclose(fid);

fid = fopen([savedir1 'DumpModFinal' int2str(DumpNo) '.txt'], 'w');
for i=1:Natoms	
    if (KSI(i,1)>StateCrit && AtomPos(i,2)~=1)
    
    fprintf(fid,'%d %d %f %f %f %d %d %d \n', AtomPos(i,1),3,AtomPos(i,3),AtomPos(i,4),AtomPos(i,5),0,0,0);
    else
        fprintf(fid,'%d %d %f %f %f %d %d %d \n', AtomPos(i,1),AtomPos(i,2),AtomPos(i,3),AtomPos(i,4),AtomPos(i,5),0,0,0);
    end
end	
fclose(fid);

end

%%

clear
TotDump=56;
for DumpNo=1:TotDump
savedir='/Users/peymansaidi/Desktop/untitled folder/Step Kinetic Coeficient/crystallayers/3-interfacefinder/';
AtomPos=importdata([savedir 'Final' int2str(DumpNo) '.txt']);
[Natoms,col] = size(AtomPos);
figure
for i=1:Natoms
  %  if (AtomPos(i,6)==20)
        if (AtomPos(i,7)>=0.2)
            plot3(AtomPos(i,3),AtomPos(i,4),AtomPos(i,5),'o','MarkerEdgeColor','r',...
                'MarkerFaceColor','r',...
               'MarkerSize',3 );
        else
            plot3(AtomPos(i,3),AtomPos(i,4),AtomPos(i,5),'o','MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
               'MarkerSize',2 );
        end
  %  end
    hold on

end
    set(gca,'DataAspectratio',[1 1 1])
    %pause(0.05);
    title(['Time= ' num2str(DumpNo)])
    view (0,0)
    print([savedir '/' num2str(DumpNo) '.png'],'-dpng','-r200',gcf)

end


%%

figure

plot(AtomPos(:,3),KSI(:,1));
figure
plot(AtomPos(:,3),KSIaverage(:,1));
%%

figure

for i=1:Natoms
    
  
      if (KSIaverage(i,1)>=0.2)
            plot3(AtomPos(i,3),AtomPos(i,4),AtomPos(i,5),'o','MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
               'MarkerSize',2 );
      end
       
      hold on
      view (0, -90)
    

end
        
            
