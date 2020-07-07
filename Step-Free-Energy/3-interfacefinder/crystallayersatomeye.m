clear
TotDump=183;
criterion=0.08;
savedir='/Users/peymansaidi/Desktop/StepFE-Comp/crystallayers/3-interfacefinder/';
savedir1='/Users/peymansaidi/Desktop/StepFE-Comp/crystallayers/2-crystallayers/';



for DumpNo=180:TotDump
    
fid = fopen([savedir 'dumpfinal' int2str(DumpNo) '.txt'], 'w');

positionfirst=importdata([savedir 'Final' int2str(DumpNo) '.txt']);
domain=importdata([savedir1 'Domain' int2str(DumpNo) '.txt']);

[natoms,col] = size(positionfirst);
Xlo=domain(1,1);
Xhi=domain(1,2);
Xslope=domain(1,3);
Ylo=domain(2,1);
Yhi=domain(2,2);
Yslope=domain(2,3);
Zlo=domain(3,1);
Zhi=domain(3,2);
Zslope=domain(3,3);





line1='ITEM: TIMESTEP';
line2=num2str(DumpNo);
line3='ITEM: NUMBER OF ATOMS';
line4=num2str(natoms);
line5='ITEM: BOX BOUNDS xy xz yz pp pp pp';
line6=[num2str(Xlo) ' ' num2str(Xhi) ' ' num2str(Xslope)];
line7=[num2str(Ylo) ' ' num2str(Yhi) ' ' num2str(Yslope)];
line8=[num2str(Zlo) ' ' num2str(Zhi) ' ' num2str(Zslope)];
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



for i=1:natoms
    if (positionfirst(i,2)==1)
    fprintf(fid,'%d %d %f %f %f %d %d %d \n', positionfirst(i,1),1,positionfirst(i,3),positionfirst(i,4),positionfirst(i,5),0,0,0);
    else
        if (positionfirst(i,6)>criterion)
        fprintf(fid,'%d %d %f %f %f %d %d %d \n', positionfirst(i,1),2,positionfirst(i,3),positionfirst(i,4),positionfirst(i,5),0,0,0);
        else
        fprintf(fid,'%d %d %f %f %f %d %d %d \n', positionfirst(i,1),3,positionfirst(i,3),positionfirst(i,4),positionfirst(i,5),0,0,0);
        end
    end
end
fclose(fid);

end



%%


fid = fopen([savedir 'Finaltot.txt'], 'w');
for DumpNo=1:TotDump
    DumpNo

positionfirst=importdata([savedir 'Final' int2str(DumpNo) '.txt']);
domain=importdata([savedir1 'Domain' int2str(DumpNo) '.txt']);

[natoms,col] = size(positionfirst);
Xlo=domain(1,1);
Xhi=domain(1,2);
Xslope=domain(1,3);
Ylo=domain(2,1);
Yhi=domain(2,2);
Yslope=domain(2,3);
Zlo=domain(3,1);
Zhi=domain(3,2);
Zslope=domain(3,3);





line1='ITEM: TIMESTEP';
line2=num2str(DumpNo);
line3='ITEM: NUMBER OF ATOMS';
line4=num2str(natoms);
line5='ITEM: BOX BOUNDS xy xz yz pp pp pp';
line6=[num2str(Xlo) ' ' num2str(Xhi) ' ' num2str(Xslope)];
line7=[num2str(Ylo) ' ' num2str(Yhi) ' ' num2str(Yslope)];
line8=[num2str(Zlo) ' ' num2str(Zhi) ' ' num2str(Zslope)];
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



for i=1:natoms
    if (positionfirst(i,2)==1)
    fprintf(fid,'%d %d %f %f %f %d %d %d \n', positionfirst(i,1),1,positionfirst(i,3),positionfirst(i,4),positionfirst(i,5),0,0,0);
    else
        if (positionfirst(i,6)>criterion)
        fprintf(fid,'%d %d %f %f %f %d %d %d \n', positionfirst(i,1),2,positionfirst(i,3),positionfirst(i,4),positionfirst(i,5),0,0,0);
        else
        fprintf(fid,'%d %d %f %f %f %d %d %d \n', positionfirst(i,1),3,positionfirst(i,3),positionfirst(i,4),positionfirst(i,5),0,0,0);
        end
    end
end


end
fclose(fid);
