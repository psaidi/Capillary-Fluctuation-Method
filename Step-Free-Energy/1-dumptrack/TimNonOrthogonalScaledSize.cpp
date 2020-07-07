#include <stdio.h>
#include <iostream>
#include <stdlib.h> //or cstdlib.h page175
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>


using namespace std;
#define MAXLINE 1024

int nrhomax,nrmax;
int nfrho,nrhor,nphi;
int *nrrho, *nrphi, *type2frho;
int **type2rhor, **type2phi;
double *drrho, *drphi;
double **frho,**rhor,**phi, **extend;
double ***rhor_spline,***frho_spline,***phi_spline;
double potorg;
double *syspotential, *delU, *BultzmanFac;



struct Dump {
    int natoms, natomsT;
	int *id, *Atype, *idT, *AtypeT;
	double *X, *Y, *Z, *XT,*YT, *ZT;
	double xlo, ylo, zlo, xhi, yhi,zhi, slopex, slopey, slopez;

};







void readDumporg(char const *, int);
void readDump(char *, Dump &);

int main ()
{
	int i, c, p, j, r, count, count1, comp , T, limit, NofDumps,refid,Natoms;
	double A, B, C,xref,yref,zref,criteria, delx, dely, delz, Lx, Ly, Lz, YY0, ZZ0, XX, YY,ZZ, TET, tanTET;
	Dump Dumpf;
	const char *argDumporg={"dump_step_30_110_167.txt"};
	
	NofDumps=77;
	criteria=1.6;
	readDumporg(argDumporg,NofDumps);
	p=1;
	char modfile[20];
	char moddomain[20];
	char argDump[20];
	sprintf(argDump,"Dump%d.txt",p);
	readDump(argDump, Dumpf);				//making seprate dump files
	xref=(Dumpf.xhi-Dumpf.xlo)/2;			//reference position
//	yref=(Dumpf.yhi-Dumpf.ylo)/3;
	yref=13;
	zref=(Dumpf.zhi-Dumpf.zlo)/2;
	i=0;
	A=xref-Dumpf.X[i];
	B=yref-Dumpf.Y[i];
	C=zref-Dumpf.Z[i];


	

	while (A>criteria || B>criteria || C>criteria || A<-criteria || B<-criteria || C<-criteria || i==Dumpf.natoms) {
		i=i+1;
		A=xref-Dumpf.X[i];					//finding the coordination of the reference atom
		B=yref-Dumpf.Y[i];
		C=zref-Dumpf.Z[i];
	}
		
	xref=Dumpf.X[i];						//reference atom and id
	yref=Dumpf.Y[i];
	zref=Dumpf.Z[i];
	refid=Dumpf.id[i];
	Natoms=Dumpf.natoms;
//	cout << xref<<endl;
//	cout << yref<<endl;
//	cout << zref<<endl;
//	cout << A<<endl;
//	cout << B<<endl;
//	cout << C<<endl;
//	cout << Dumpf.X[i]<<endl;
//	cout << Dumpf.Y[i]<<endl;
//	cout << Dumpf.Z[i]<<endl;
//	cout << Dumpf.id[i]<<endl;
//	cout << i<<endl;

	for (p=1; p<=NofDumps; p++) {			//shifting all of the atoms back for all of the dump files
		
		sprintf(argDump,"Dump%d.txt",p);	//reading dump files one by one
		readDump(argDump, Dumpf);
		Lx=Dumpf.xhi-Dumpf.xlo;
		Ly=(Dumpf.yhi-Dumpf.ylo)-Dumpf.slopez;
		Lz=Dumpf.zhi-Dumpf.zlo;
	//	TET=atan(Dumpf.slopez/(Dumpf.zhi-Dumpf.zlo));
		tanTET=Dumpf.slopez/(Dumpf.zhi-Dumpf.zlo);
	//  Lz=(Dumpf.zhi-Dumpf.zlo)/cos(TET);
		i=0;
		
		
		while (refid!=Dumpf.id[i])			//finding the reference atom in the new sump file
			i=i+1;
		
		delx=xref-Dumpf.X[i];
		delz=zref-Dumpf.Z[i];
		dely=yref-Dumpf.Y[i]-tanTET*delz;

		
		
	//	delz=delz*cos(TET)+dely*sin(TET);
	//	dely=-delz*sin(TET)+dely*cos(TET);
		
		sprintf(modfile,"ModDump%d.txt",p);//dumpfile with the shifted back position of the atoms
		ofstream myfile (modfile);
		
		sprintf(moddomain,"ModDomain%d.txt",p);//modified domain of tne system, the same as the original one
		ofstream myfileD (moddomain);
		if (myfileD.is_open()) myfileD <<setprecision(10)<<Dumpf.xlo<<" "<<Dumpf.xhi<<" "<<Dumpf.slopex<<endl;
		if (myfileD.is_open()) myfileD <<setprecision(10)<<Dumpf.ylo<<" "<<Dumpf.yhi<<" "<<Dumpf.slopey<<endl;
		if (myfileD.is_open()) myfileD <<setprecision(10)<<Dumpf.zlo<<" "<<Dumpf.zhi<<" "<<Dumpf.slopez<<endl;
		
		
				 
		
		for (j=0; j<Dumpf.natoms; j++) {		//shifting back the atoms and taking care of the periodicity of the system
			
		//	YY0=Dumpf.Y[j];
		//	ZZ0=Dumpf.Z[j];
			Dumpf.Y[j]=Dumpf.Y[j]-Dumpf.Z[j]*tanTET;
		//	Dumpf.Z[j]=ZZ0*cos(TET)+YY0*sin(TET);				//rotate with (-tet) degree
		//	Dumpf.Y[j]=-ZZ0*sin(TET)+YY0*cos(TET);
			XX=Dumpf.X[j]+delx;
			YY=Dumpf.Y[j]+dely;
			ZZ=Dumpf.Z[j]+delz;
			if (XX>Dumpf.xhi) XX=XX-Lx;
			if (XX<Dumpf.xlo) XX=XX+Lx;
			if (YY>Dumpf.yhi-Dumpf.slopez) YY=YY-Ly;
			if (YY<Dumpf.ylo) YY=YY+Ly;
			if (ZZ>Dumpf.zhi) ZZ=ZZ-Lz;
			if (ZZ<Dumpf.zlo) ZZ=ZZ+Lz;	
			Dumpf.X[j]=XX;
			Dumpf.Y[j]=YY;
			Dumpf.Z[j]=ZZ;
		//	Dumpf.Z[j]=ZZ*cos(TET)-YY*sin(TET);			//rotateback with (tet) degree
		//	Dumpf.Y[j]=ZZ*sin(TET)+YY*cos(TET);
			Dumpf.Y[j]=Dumpf.Y[j]+Dumpf.Z[j]*tanTET;
			if (myfile.is_open()) myfile <<setprecision(10)<<Dumpf.id[j]<<" "<<Dumpf.Atype[j]<<" "<<Dumpf.X[j]<<" "<<Dumpf.Y[j]<<" "<<Dumpf.Z[j]<<" "<<endl;
			
			
			
		}
		myfile.close();
		myfileD.close();

		
		
		
	}
	
	int nheader3=9;							//making the modified dump file total.
	char line1[MAXLINE];
	ofstream myfileToT;
	myfileToT.open ("DumpModTot.txt");
	
	FILE * pFile;
	FILE * ppFile;
	for (r=1; r<=NofDumps; r++) {
		sprintf(argDump,"Dump%d.txt",r);
		pFile = fopen (argDump , "r");
		
		for (i=0; i<nheader3; i++){
			fgets (line1 , MAXLINE , pFile);
			
			if (myfileToT.is_open()) 
				myfileToT <<line1;
		}
		
		sprintf(modfile,"ModDump%d.txt",r);
		ppFile = fopen (modfile , "r");
		
		for (j=0; j<Natoms; j++){
			fgets (line1 , MAXLINE , ppFile);
			
			if (myfileToT.is_open()) 
				myfileToT <<line1;
		}
	}

	
	myfileToT.close();
	


	
	return 0;
}




//////////////////////Functions//////////////////////////////

void readDumporg(char const *Dumporgfilename, int Noffiles)
{
	// reading the dump file and extracting atom number, domain and coordination of atoms
	int i, j,q, natoms, IDTEMP,ATYPETEMP;
	double XTEMP, YTEMP, ZTEMP, tantet, Zs, ymin, ymax;
	double HiLo[9];
	int nheader1=4;
	int nheader2=1;
	int nheader3=3;
	int nheader4=1;
	char line[MAXLINE];
	FILE * dumporgFile;
	char filename[20];
	
	dumporgFile = fopen (Dumporgfilename,"r");
	
	if (dumporgFile!=NULL)
	{
		for (q=1; q<=Noffiles; q++) {
			sprintf(filename,"Dump%d.txt",q);
			ofstream myfile1 (filename);
			
			for (i=0; i<nheader1; i++){
				fgets (line , MAXLINE , dumporgFile);
				if (myfile1.is_open()) 
					myfile1 <<line;
			}
			
			sscanf(line,"%d",&natoms);
			
			
			for (i=0; i<nheader2; i++){
				fgets (line , MAXLINE , dumporgFile);
				if (myfile1.is_open()) 
					myfile1 <<line;
			}
			
			for (i=0; i<nheader3; i++){
				fgets (line , MAXLINE , dumporgFile);
				sscanf(line,"%lg %lg %lg",&HiLo[i*3], &HiLo[i*3+1] , &HiLo[i*3+2]);
				if (myfile1.is_open()) 
					myfile1 <<line;
			}
			
			tantet=HiLo[8]/(HiLo[7]-HiLo[6]);
			ymin=HiLo[3];
			ymax=HiLo[4];
			
			
			
			for (i=0; i<nheader4; i++){
				fgets (line , MAXLINE , dumporgFile);
				if (myfile1.is_open()) 
					myfile1 <<"ITEM: ATOMS id type x y z "<<endl;
			}
			
			
			for (i=0; i<natoms; i++) {
				fgets (line , MAXLINE , dumporgFile);
				sscanf(line,"%d %d %lg %lg %lg",&IDTEMP ,&ATYPETEMP,&XTEMP, &YTEMP, &ZTEMP);
				Zs=tantet*(HiLo[6]*(1-ZTEMP)+ZTEMP*HiLo[7]);
				
				HiLo[3]=ymin+Zs;
				HiLo[4]=ymax-HiLo[8]+Zs;
				
				if (myfile1.is_open()) 
					myfile1<<IDTEMP<<" "<<ATYPETEMP<<" "<<HiLo[0]*(1-XTEMP)+XTEMP*HiLo[1]<<" "<<HiLo[3]*(1-YTEMP)+YTEMP*HiLo[4]<<" "<<HiLo[6]*(1-ZTEMP)+ZTEMP*HiLo[7]<<endl;
				
			}
			
			
			
			myfile1.close();
			
			
			
		}
		fclose (dumporgFile);
			}
}



void readDump(char *Dumpfilename, Dump &file)
{
// reading the dump file and extracting atom number, domain and coordination of atoms
int i, j;
int nheader1=4;
char line[MAXLINE];
FILE * dumpFile;

dumpFile = fopen (Dumpfilename,"r");

if (dumpFile!=NULL)
{

	for (i=0; i<nheader1; i++) fgets (line , MAXLINE , dumpFile);

		sscanf(line,"%d",&file.natoms);
		file.id=new int[file.natoms];
		file.Atype=new int[file.natoms];
		file.X=new double[file.natoms];
		file.Y=new double[file.natoms];
		file.Z=new double[file.natoms];
		
		
		fgets (line , MAXLINE , dumpFile);
		fgets (line , MAXLINE , dumpFile);
		sscanf(line,"%lg  %lg %lg",&file.xlo ,&file.xhi, &file.slopex);
		
		fgets (line , MAXLINE , dumpFile);
		sscanf(line,"%lg  %lg %lg",&file.ylo ,&file.yhi, &file.slopey);
		
		fgets (line , MAXLINE , dumpFile);
		sscanf(line,"%lg  %lg %lg",&file.zlo ,&file.zhi, &file.slopez);
		
		fgets (line , MAXLINE , dumpFile);
		for (i=0; i<file.natoms; i++) {
			fgets (line , MAXLINE , dumpFile);
			sscanf(line,"%d %d %lg %lg %lg",&file.id[i] ,&file.Atype[i], &file.X[i], &file.Y[i], &file.Z[i]);
			
		}
	fclose (dumpFile);
	
}
}

