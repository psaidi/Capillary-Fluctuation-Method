%this file extract the acurate position of interface

clear
TotDump=1;
a=5.46;

CriticalValue=0.1;

delX=0.8*a;
delZ=0.5*a;

%37 samples from 525 file
 FileTitle='-525-37';
 Ymin1=11;
 Ymax1=14;
 Ymin2=58.7;
 Ymax2=61.7;

% %183 samples from 883 file (just the second part is valid
%  FileTitle='-243-91';
%  Ymin1=70;
%  Ymax1=73.5;
%  Ymin2=105;
%  Ymax2=108;

%179 samples from 705 file (just the second part is valid
% FileTitle='-218-81';
% Ymin1=100;
% Ymax1=103;
% Ymin2=26.2;
% Ymax2=29.3;

savedir='/Users/peymansaidi/Desktop/StepFE-Comp/crystallayers/2-crystallayers/';
savedir1='/Users/peymansaidi/Desktop/StepFE-Comp/crystallayers/3-interfacefinder/';
savedir2='/Users/peymansaidi/Desktop/StepFE-Comp/crystallayers/5-Boundary/';



Domain=importdata([savedir 'Domain1.txt']);
fid1 = fopen([savedir2 'BoundaryCoord1' FileTitle '.txt'], 'w');
fid2 = fopen([savedir2 'BoundaryCoord2' FileTitle '.txt'], 'w');
fid3 = fopen([savedir2 'Amplitude1' FileTitle '.txt'], 'w');
fid4 = fopen([savedir2 'Amplitude2' FileTitle '.txt'], 'w');

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

Totx=ceil(Lx/delX);
Totz=ceil(Lz/delZ);

OrderParaVal1=zeros(TotDump,Totx,Totz);
OrderParaCount1=zeros(TotDump,Totx,Totz);
OrderParaValAve1=zeros(TotDump,Totx*Totz);

OrderParaVal2=zeros(TotDump,Totx,Totz);
OrderParaCount2=zeros(TotDump,Totx,Totz);
OrderParaValAve2=zeros(TotDump,Totx*Totz);
AcceptedFinalRot1=zeros(TotDump,Totz,Totx);
AcceptedHeight1=zeros(TotDump,Totx-1);
AcceptedFinalRot2=zeros(TotDump,Totz,Totx);
AcceptedHeight2=zeros(TotDump,Totx-1);


 a1=zeros(TotDump,Totx-1);
 b1=zeros(TotDump,Totx-1);
 a2=zeros(TotDump,Totx-1);
 b2=zeros(TotDump,Totx-1);
 
 xFourier=linspace(delX/2,(Totx-2)*delX+delX/2,Totx-1);
 
 xFourier(1,:)=xFourier(1,:)-(Totx-1)*delX/2;
 xFourier(1,:)=linspace(-Lx/2,Lx/2,Totx-1);
 LxFourier=xFourier(1,end)-xFourier(1,1);

%Slicing the domain, the same discritization for both interface planes

IdSlice=0;
NeighbourDomainID=zeros(Totx*Totz,8);
for k=1:Totx
    for j=1:Totz
        NoOfNeigh=0;
        IdSlice=(j-1)*Totx+k;
        
            for J=-1:1
                for K=-1:1
                    CounterNeigh=(j+J-1)*(Totx)+(K+k);

                    if((J+j-1)>=0 && (J+j-1)<=Totz && (K+k)>0 && (K+k)<=Totx && CounterNeigh~=IdSlice && CounterNeigh<=Totx*Totz)
                        NoOfNeigh=NoOfNeigh+1;
                        NeighbourDomainID(IdSlice,1)=NoOfNeigh;
                        NeighbourDomainID(IdSlice,NoOfNeigh+1)=CounterNeigh;
                    end

                end
            end
   end
end



for DumpNo=1:TotDump
    DumpNo
%reading the position of the atoms 

AtomPos=importdata([savedir1 'Final' int2str(DumpNo) '.txt']);
Domain=importdata([savedir 'Domain' int2str(DumpNo) '.txt']);

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

Totx=ceil(Lx/delX);
Totz=ceil(Lz/delZ);


%finding the atoms in the target range
counter1=0;
counter2=0;

for i=1:Natoms
    
 x1=AtomPos(i,3);
 y1=AtomPos(i,4);
 z1=AtomPos(i,5);
 
 if (y1>Ymin1 && y1<Ymax1)
     Nx=ceil((x1-xlo)/delX);
     Nz=ceil((z1-zlo)/delZ);
 
         
     OrderParaVal1(DumpNo,Nx,Nz)=OrderParaVal1(DumpNo,Nx,Nz)+AtomPos(i,6);
     OrderParaCount1(DumpNo,Nx,Nz)=OrderParaCount1(DumpNo,Nx,Nz)+1;

 end

  if (y1>Ymin2 && y1<Ymax2)
     Nx=ceil((x1-xlo)/delX);
     Nz=ceil((z1-zlo)/delZ);
     OrderParaVal2(DumpNo,Nx,Nz)=OrderParaVal2(DumpNo,Nx,Nz)+AtomPos(i,6);
     OrderParaCount2(DumpNo,Nx,Nz)=OrderParaCount2(DumpNo,Nx,Nz)+1;


          
  end

end




for k=1:Totx
    for j=1:Totz
        IdSlice=(j-1)*Totx+k;
  OrderParaValAve1(DumpNo,IdSlice)=OrderParaVal1(DumpNo,k,j)/OrderParaCount1(DumpNo,k,j);
  OrderParaValAve2(DumpNo,IdSlice)=OrderParaVal2(DumpNo,k,j)/OrderParaCount2(DumpNo,k,j);
    end
end




for k=1:Totx
    for j=1:Totz
        IdSlice=(j-1)*Totx+k;
        if (OrderParaValAve1(DumpNo,IdSlice)>CriticalValue)
            State1(DumpNo,IdSlice)=1;
        else
            State1(DumpNo,IdSlice)=0;
        end
        if (OrderParaValAve2(DumpNo,IdSlice)>CriticalValue)
            State2(DumpNo,IdSlice)=1;
        else
            State2(DumpNo,IdSlice)=0;
        end
  
    end
end


AcceptedSolids1=zeros(Totx*Totz,1);
AcceptedSolids1(1:Totx,1)=1;
State1(DumpNo,1:Totx)=1;

for i=Totx+1:Totx*Totz
    if (State1(DumpNo,i)==1)
        
        
    for I=2:NeighbourDomainID(i,1)+1
        if (AcceptedSolids1(NeighbourDomainID(i,I),1)==1 && NeighbourDomainID(i,I)<i)
            AcceptedSolids1(i,1)=1;
        end
    end
            
    end
end


AcceptedLiquids2=zeros(Totx*Totz,1);
AcceptedLiquids2(:,1)=1;
AcceptedLiquids2(1:Totx,1)=0;
State2(DumpNo,1:Totx)=0;

for i=Totx+1:Totx*Totz
    if (State2(DumpNo,i)==0)
        
        
    for I=2:NeighbourDomainID(i,1)+1
        if (AcceptedLiquids2(NeighbourDomainID(i,I),1)==0 && NeighbourDomainID(i,I)<i)
            AcceptedLiquids2(i,1)=0;
        end
    end
            
    end
end





AcceptedLiquids1=zeros(Totx*Totz,1);
AcceptedLiquids1(:,1)=1;
AcceptedLiquids1((Totx-1)*Totz:Totx*Totz,1)=0;
State1(DumpNo,(Totx-1)*Totz:Totx*Totz,1)=0;

StartPoint=(Totx-1)*Totz-1;
for i=StartPoint:-1:1
    if (State1(DumpNo,i)==0)       
    for I=2:NeighbourDomainID(i,1)+1
        if (AcceptedLiquids1(NeighbourDomainID(i,I),1)==0 && NeighbourDomainID(i,I)>i)
            AcceptedLiquids1(i,1)=0;
        end
    end
            
    end
end


AcceptedSolids2=zeros(Totx*Totz,1);
AcceptedSolids2((Totx-1)*Totz:Totx*Totz,1)=1;
State2(DumpNo,(Totx-1)*Totz:Totx*Totz,1)=1;

StartPoint=(Totx-1)*Totz-1;
for i=StartPoint:-1:1
    if (State2(DumpNo,i)==1)       
    for I=2:NeighbourDomainID(i,1)+1
        if (AcceptedSolids2(NeighbourDomainID(i,I),1)==1 && NeighbourDomainID(i,I)>i)
            AcceptedSolids2(i,1)=1;
        end
    end
            
    end
end





%making the continus structure


AcceptedFinal1=zeros(TotDump,Totx*Totz);

for i=1:Totx*Totz
    if (State1(DumpNo,i)==0 && AcceptedLiquids1(i,1)==0 && AcceptedSolids1(i,1)==0)
        AcceptedFinal1(DumpNo,i)=0;
    elseif (State1(DumpNo,i)==1 && AcceptedLiquids1(i,1)==1 && AcceptedSolids1(i,1)==1)
        AcceptedFinal1(DumpNo,i)=1;
    elseif (State1(DumpNo,i)==1 && AcceptedLiquids1(i,1)==1 && AcceptedSolids1(i,1)==0)
        AcceptedFinal1(DumpNo,i)=0;
    elseif (State1(DumpNo,i)==0 && AcceptedLiquids1(i,1)==1 && AcceptedSolids1(i,1)==0)
        AcceptedFinal1(DumpNo,i)=1;    
    end
end


AcceptedFinal2=zeros(TotDump,Totx*Totz);

for i=1:Totx*Totz
    if (State2(DumpNo,i)==0 && AcceptedLiquids2(i,1)==0 && AcceptedSolids2(i,1)==0)
        AcceptedFinal2(DumpNo,i)=0;
    elseif (State2(DumpNo,i)==1 && AcceptedLiquids2(i,1)==1 && AcceptedSolids2(i,1)==1)
        AcceptedFinal2(DumpNo,i)=1;
    elseif (State2(DumpNo,i)==1 && AcceptedLiquids2(i,1)==1 && AcceptedSolids2(i,1)==0)
        AcceptedFinal2(DumpNo,i)=0;
    elseif (State2(DumpNo,i)==0 && AcceptedLiquids2(i,1)==1 && AcceptedSolids2(i,1)==0)
        AcceptedFinal2(DumpNo,i)=1;    
    end
end




%making a 2D picture of each snapshot

 for j=1:Totz
   AcceptedFinalRot1(DumpNo,j,:)=AcceptedFinal1(DumpNo,(j-1)*Totx+1:j*Totx);
 end
 
 for j=1:Totz
   AcceptedFinalRot2(DumpNo,j,:)=AcceptedFinal2(DumpNo,(j-1)*Totx+1:j*Totx);
 end
 
 
  
 %making it continus in z direction
 
 switchList=zeros(1,10);
 for i=1:Totx-1
     NumofSwitch=0;
  for j=1:Totz-1
  if (AcceptedFinalRot1(DumpNo,j,i)~=AcceptedFinalRot1(DumpNo,j+1,i))
      NumofSwitch=NumofSwitch+1;
      switchList(1,1)=NumofSwitch;
      switchList(NumofSwitch+1,1)=j;
  end
  end
    if (switchList(1,1)==1)
   AcceptedHeight1(DumpNo,i)=switchList(2,1)*delZ;
  elseif (switchList(1,1)==3)
      Switch1=switchList(3,1)-switchList(2,1);
      Switch2=switchList(4,1)-switchList(3,1);
      if (Switch2>Switch1)
      AcceptedFinalRot1(DumpNo,1:switchList(4,1),i)=1;
      AcceptedHeight1(DumpNo,i)=switchList(4,1)*delZ;
      else
      AcceptedFinalRot1(DumpNo,switchList(2,1)+1:end,i)=0;
      AcceptedHeight1(DumpNo,i)=(switchList(2,1))*delZ; 
      end
    else
      AcceptedFinalRot1(DumpNo,1:sum(AcceptedFinalRot1(DumpNo,:,i)),i)=1;
      AcceptedFinalRot1(DumpNo,sum(AcceptedFinalRot1(DumpNo,:,i)):end,i)=0;
      AcceptedHeight1(DumpNo,i)=sum(AcceptedFinalRot1(DumpNo,:,i))*delZ;                  
    end
 end
 
 
 switchList=zeros(1,10);
 for i=1:Totx-1
     NumofSwitch=0;
  for j=1:Totz-1
  if (AcceptedFinalRot2(DumpNo,j,i)~=AcceptedFinalRot2(DumpNo,j+1,i))
      NumofSwitch=NumofSwitch+1;
      switchList(1,1)=NumofSwitch;
      switchList(NumofSwitch+1,1)=j;
  end
  end
    if (switchList(1,1)==1)
   AcceptedHeight2(DumpNo,i)=switchList(2,1)*delZ;
  elseif (switchList(1,1)==3)
      Switch1=switchList(3,1)-switchList(2,1);
      Switch2=switchList(4,1)-switchList(3,1);
      if (Switch2>Switch1)
      AcceptedFinalRot2(DumpNo,1:switchList(4,1),i)=0;
      AcceptedHeight2(DumpNo,i)=switchList(4,1)*delZ;
      else
      AcceptedFinalRot2(DumpNo,switchList(2,1):end,i)=1;
      AcceptedHeight2(DumpNo,i)=(switchList(2,1))*delZ; 
      end
    else
      AcceptedFinalRot2(DumpNo,sum(AcceptedFinalRot2(DumpNo,:,i)):end,i)=1;
      AcceptedFinalRot2(DumpNo,1:sum(AcceptedFinalRot2(DumpNo,:,i))-1,i)=0;
      AcceptedHeight2(DumpNo,i)=sum(AcceptedFinalRot2(DumpNo,:,i))*delZ;                  
    end
 end
  

    figure('Position',[300 300 250 480]);
 
for i=1:Totz
    for k=1:Totx-1
        rtr=(i-1)*Totx+k;
        if (State1(DumpNo,rtr)==0)
        plot(k,i,'s','MarkerEdgeColor','g','MarkerFaceColor','g')
         
        else
        plot(k,i,'s','MarkerEdgeColor','r','MarkerFaceColor','r')
        end
        hold on
         

    end
end




    figure('Position',[300 300 250 480]);
 
for i=1:Totz
    for k=1:Totx-1
        rtr=(i-1)*Totx+k;
        if (AcceptedSolids1(rtr,1)==0)
        plot(k,i,'s','MarkerEdgeColor','g','MarkerFaceColor','g')
         
        else
        plot(k,i,'s','MarkerEdgeColor','r','MarkerFaceColor','r')
        end
        hold on
         

    end
end




    figure('Position',[300 300 250 480]);
 
for i=1:Totz
    for k=1:Totx-1
        rtr=(i-1)*Totx+k;
        if (AcceptedLiquids1(rtr,1)==0)
        plot(k,i,'s','MarkerEdgeColor','g','MarkerFaceColor','g')
         
        else
        plot(k,i,'s','MarkerEdgeColor','r','MarkerFaceColor','r')
        end
        hold on
         

    end
end


     figure('Position',[300 300 250 480]);
     
 
for i=1:Totz
    for k=1:Totx-1
        rtr=(i-1)*Totx+k;
        if (AcceptedFinal1(DumpNo,rtr)==0)
        plot(k,i,'s','MarkerEdgeColor','g','MarkerFaceColor','g')
         
        else
        plot(k,i,'s','MarkerEdgeColor','r','MarkerFaceColor','r')
         end
        hold on
         

     end
 end



 



     figure('Position',[300 300 250 480]);
     
  
 for i=1:Totz
     for k=1:Totx-1
         
         if (AcceptedFinalRot1(DumpNo,i,k)==0)
         plot(k,i,'s','MarkerEdgeColor','g','MarkerFaceColor','g')
          
         else
         plot(k,i,'s','MarkerEdgeColor','r','MarkerFaceColor','r')
          end
         hold on
          
 
      end
  end
 
 hold on
 
plot(AcceptedHeight1(DumpNo,:)/delZ, 'black')

title(['Time= ' num2str(DumpNo)])
set(gca,'DataAspectratio',[1 1 1])
print([savedir2 '/' num2str(1000+DumpNo) '.png'],'-dpng','-r200',gcf)





fprintf(fid1,'%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n', DumpNo...
    ,AcceptedHeight1(DumpNo,1),AcceptedHeight1(DumpNo,2),AcceptedHeight1(DumpNo,3)...
    ,AcceptedHeight1(DumpNo,4),AcceptedHeight1(DumpNo,5),AcceptedHeight1(DumpNo,6)...
    ,AcceptedHeight1(DumpNo,7),AcceptedHeight1(DumpNo,8),AcceptedHeight1(DumpNo,9)...
    ,AcceptedHeight1(DumpNo,10),AcceptedHeight1(DumpNo,11),AcceptedHeight1(DumpNo,12)...
    ,AcceptedHeight1(DumpNo,13),AcceptedHeight1(DumpNo,14),AcceptedHeight1(DumpNo,15)...
    ,AcceptedHeight1(DumpNo,16),AcceptedHeight1(DumpNo,17),AcceptedHeight1(DumpNo,18)...
    ,AcceptedHeight1(DumpNo,19),AcceptedHeight1(DumpNo,20),AcceptedHeight1(DumpNo,21)...
    );

fprintf(fid2,'%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n', DumpNo...
    ,AcceptedHeight2(DumpNo,1),AcceptedHeight2(DumpNo,2),AcceptedHeight2(DumpNo,3)...
    ,AcceptedHeight2(DumpNo,4),AcceptedHeight2(DumpNo,5),AcceptedHeight2(DumpNo,6)...
    ,AcceptedHeight2(DumpNo,7),AcceptedHeight2(DumpNo,8),AcceptedHeight2(DumpNo,9)...
    ,AcceptedHeight2(DumpNo,10),AcceptedHeight2(DumpNo,11),AcceptedHeight2(DumpNo,12)...
    ,AcceptedHeight2(DumpNo,13),AcceptedHeight2(DumpNo,14),AcceptedHeight2(DumpNo,15)...
    ,AcceptedHeight2(DumpNo,16),AcceptedHeight2(DumpNo,17),AcceptedHeight2(DumpNo,18)...
    ,AcceptedHeight2(DumpNo,19),AcceptedHeight2(DumpNo,20),AcceptedHeight2(DumpNo,21)...
    );

% manual method of calculation for the first interface
for k=0:Totx-2
    for n=0:Totx-2
    a1(DumpNo,k+1)=a1(DumpNo,k+1)+(AcceptedHeight1(DumpNo,n+1)*cos(2*n*pi*k/(Totx-1)));
    b1(DumpNo,k+1)=b1(DumpNo,k+1)+(AcceptedHeight1(DumpNo,n+1)*sin(2*n*pi*k/(Totx-1)));
    
    end
    AmpSQR1(DumpNo,k+1)=a1(DumpNo,k+1)^2+b1(DumpNo,k+1)^2
    
end


% manual method of calculation for the second interface
for k=0:Totx-2
    for n=0:Totx-2
    a2(DumpNo,k+1)=a2(DumpNo,k+1)+(AcceptedHeight1(DumpNo,n+1)*cos(2*n*pi*k/(Totx-1)));
    b2(DumpNo,k+1)=b2(DumpNo,k+1)+(AcceptedHeight1(DumpNo,n+1)*sin(2*n*pi*k/(Totx-1)));
    
    end

    AmpSQR2(DumpNo,k+1)=a2(DumpNo,k+1)^2+b2(DumpNo,k+1)^2;
    
end


% Matlab method of calculation for the first interface
AmpMatlab1(DumpNo,:)=abs(fft(AcceptedHeight1(DumpNo,:))).^2%we ignore the last column

% Matlab method of calculation for the second interface
AmpMatlab2(DumpNo,:)=abs(fft(AcceptedHeight2(DumpNo,:))).^2;



fprintf(fid3,'%d %g %g %g %g %g %g %g %g %g %g %g \n', DumpNo...
    ,AmpMatlab1(DumpNo,1),AmpMatlab1(DumpNo,2),AmpMatlab1(DumpNo,3)...
    ,AmpMatlab1(DumpNo,4),AmpMatlab1(DumpNo,5),AmpMatlab1(DumpNo,6)...
    ,AmpMatlab1(DumpNo,7),AmpMatlab1(DumpNo,8),AmpMatlab1(DumpNo,9)...
    ,AmpMatlab1(DumpNo,10),AmpMatlab1(DumpNo,11));

fprintf(fid4,'%d %g %g %g %g %g %g %g %g %g %g %g \n', DumpNo...
    ,AmpMatlab2(DumpNo,1),AmpMatlab2(DumpNo,2),AmpMatlab2(DumpNo,3)...
    ,AmpMatlab2(DumpNo,4),AmpMatlab2(DumpNo,5),AmpMatlab2(DumpNo,6)...
    ,AmpMatlab2(DumpNo,7),AmpMatlab2(DumpNo,8),AmpMatlab2(DumpNo,9)...
    ,AmpMatlab2(DumpNo,10),AmpMatlab2(DumpNo,11));





%figure
% plot(AmpSQR1(DumpNo,:))
 
% figure
 
% plot(AmpMatlab1(DumpNo,2:end), 'r')
end
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);





 %%
 
 
 
 
 
     figure('Position',[300 300 310 590]);
 
 for i=1:Totz
     for k=1:Totx
         rtr=(i-1)*Totx+k;
         if (State1(DumpNo,rtr)==0)
         plot(k,i,'s','MarkerEdgeColor','g','MarkerFaceColor','g')
         
         else
         plot(k,i,'s','MarkerEdgeColor','r','MarkerFaceColor','r')
         end
         hold on
         

     end
 end
 
 
 
 
 
 
 
 
 

 
     figure('Position',[300 300 310 590]);
 
 for i=1:Totz
     for k=1:Totx
         rtr=(i-1)*Totx+k;
         if (State2(DumpNo,rtr)==0)
         plot(k,i,'s','MarkerEdgeColor','g','MarkerFaceColor','g')
         
         else
         plot(k,i,'s','MarkerEdgeColor','r','MarkerFaceColor','r')
         end
         hold on
         

     end
 end
 

 
      figure('Position',[300 300 310 590]);
 
 for i=1:Totz
     for k=1:Totx
         
         if (AcceptedFinalRot2(DumpNo,i,k)==0)
         plot(k,i,'s','MarkerEdgeColor','g','MarkerFaceColor','g')
         
         else
         plot(k,i,'s','MarkerEdgeColor','r','MarkerFaceColor','r')
         end
         hold on
         

     end
 end
 
 hold on
 
plot(AcceptedHeight2(DumpNo,:)/delZ, '--black')