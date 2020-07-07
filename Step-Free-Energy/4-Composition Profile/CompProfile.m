clear
TotDump=119;
Composition=zeros(TotDump,1);

figure

dX=2;
savedir='/Users/peymansaidi/Desktop/untitled folder/Step Kinetic Coeficient/crystallayers/4-Composition Profile/';
savedir1='/Users/peymansaidi/Desktop/untitled folder/Step Kinetic Coeficient/crystallayers/3-interfacefinder/';
savedir2='/Users/peymansaidi/Desktop/untitled folder/Step Kinetic Coeficient/crystallayers/2-crystallayers/';

for DumpNo=1:TotDump
    DumpNo
    
positionfirst=importdata([savedir1 'Final' int2str(DumpNo) '.txt']);
domain=importdata([savedir2 'Domain' int2str(DumpNo) '.txt']);

[natoms,col] = size(positionfirst);
Xlo=domain(1,1);
Xhi=domain(1,2);
Ylo=domain(2,1);
Yhi=domain(2,2);
Zlo=domain(3,1);
Zhi=domain(3,2);

Lx=Xhi-Xlo;

Nlayers=floor(Lx/dX)+1;
%Concentration=zeros(TotDump,Nlayers,3);
Concentration=zeros(Nlayers,3);

criterion=0.2;

COMPAl=0;
COMPSiS=0;
COMPSiL=0;


for i=1:natoms
    
    ColumnNu=floor((positionfirst(i,3)-Xlo)/dX)+1;
    
    if (positionfirst(i,2)==1)
     %   Concentration(DumpNo,ColumnNu,1)=Concentration(DumpNo,ColumnNu,1)+1;
                Concentration(ColumnNu,1)=Concentration(ColumnNu,1)+1;
                COMPAl=COMPAl+1;

        
    else
        if (positionfirst(i,7)>criterion)
            %Concentration(DumpNo,ColumnNu,2)=Concentration(DumpNo,ColumnNu,2)+1;
            Concentration(ColumnNu,2)=Concentration(ColumnNu,2)+1;
            COMPSiS=COMPSiS+1;
        else
           % Concentration(DumpNo,ColumnNu,3)=Concentration(DumpNo,ColumnNu,3)+1;
                        Concentration(ColumnNu,3)=Concentration(ColumnNu,3)+1;
            COMPSiL=COMPSiL+1;
        end
    end
end

Composition(DumpNo,1)=COMPSiL/(COMPSiL+COMPAl);


x=linspace(0,Xhi-Xlo,Nlayers);


plot (x,Concentration(:,3),'--o','LineWidth',2)
hold on 

plot (x,Concentration(:,2),'--r*','LineWidth',2)

plot (x,Concentration(:,1),'--green+','LineWidth',2)
plot (x,Concentration(:,1)+Concentration(:,3),'--black','LineWidth',2)

xlim ([Xlo Xhi])
    pause(0.05);
    title(['Composition Of Liquid= ' num2str(Composition(DumpNo,1))])
    print([savedir '/' num2str(1000+DumpNo) '.png'],'-dpng','-r200')


hold off

end
%%

figure 
Composition(1,1)=0.9;
plot (Composition(2:end,1), 'r*')
%%

x=linspace(0,Xhi-Xlo,Nlayers);
plot (x,Concentration(:,3),'--o')
hold on 

plot (x,Concentration(:,2),'--r*')

plot (x,Concentration(:,1),'--green+')
plot (x,Concentration(:,1)+Concentration(:,3),'--black')



