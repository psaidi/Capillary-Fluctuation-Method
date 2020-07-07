clear
TotDump=56;
savedir='/Users/peymansaidi/Desktop/untitled folder/Step Kinetic Coeficient/crystallayers/3-interfacefinder/';
savedir1='/Users/peymansaidi/Desktop/untitled folder/Step Kinetic Coeficient/crystallayers/2-crystallayers/';
fid = fopen([savedir 'Finalonelayer.txt'], 'w');
for DumpNo=1:TotDump

positionfirst=importdata([savedir 'Final' int2str(DumpNo) '.txt']);
domain=importdata([savedir1 'Domain' int2str(DumpNo) '.txt']);

[natoms,col] = size(positionfirst);
Xlo=domain(1,1);
Xhi=domain(1,2);
Ylo=domain(2,1);
Yhi=domain(2,2);
Zlo=domain(3,1);
Zhi=domain(3,2);





line1='ITEM: TIMESTEP';
line2=num2str(DumpNo);
line3='ITEM: NUMBER OF ATOMS';
line4=num2str(natoms);
line5='ITEM: BOX BOUNDS pp pp pp';
line6=[num2str(Xlo) ' ' num2str(Xhi)];
line7=[num2str(Ylo) ' ' num2str(Yhi)];
line8=[num2str(Zlo) ' ' num2str(Zhi)];
line9='ITEM: ATOMS id type x y z ix iy iz';


fprintf(fid,'%s\n', line1);
fprintf(fid,'%s\n', line2);
fprintf(fid,'%s\n', line3);
fprintf(fid,'%s\n', line4);
fprintf(fid,'%s\n', line5);
fprintf(fid,'%s\n', line6);
fprintf(fid,'%s\n', line7);
fprintf(fid,'%s\n', line8);
fprintf(fid,'%s \n', line9);


criterion=0.2;
layernumber=20;
for i=1:natoms	
    if (positionfirst(i,7)>criterion && positionfirst(i,6)==layernumber)
    fprintf(fid,'%d %d %f %f %f %d %d %d \n', positionfirst(i,1),2,positionfirst(i,3),positionfirst(i,4),positionfirst(i,5),0,0,0);
    else
    fprintf(fid,'%d %d %f %f %f %d %d %d \n', positionfirst(i,1),3,positionfirst(i,3),positionfirst(i,4),positionfirst(i,5),0,0,0);
    end
end	


end
fclose(fid);

