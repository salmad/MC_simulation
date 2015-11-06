# include "molecule.h"
# include "global.h"
# include "include/polymer.h"
# include "include/integrate.h"

# include <vector>
# include <cstdlib>
# include <iomanip>
# include <ctime>
# include <cmath>
# include <string>
# include <sys/stat.h>
# include <iostream>
# include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

inline double Distance(molecule m1, molecule m2){
// minimum image convention implemented
//	if (geometry==0){

//			double dx=(m1.x-m2.x);
//			double dy=(m1.y-m2.y);
//			double dz=(m1.z-m2.z);
//			// periodic conditions
//			if (dx>L/2.){  dx-=L;} else if (dx<-L/2.){  dx+=L;}
//			if (dy>L/2.){  dy-=L;} else if (dy<-L/2.){  dy+=L;}
//			if (dz>L/2.){  dz-=L;} else if (dz<-L/2.){  dz+=L;}
//			return sqrt(dx*dx+dy*dy+dz*dz);

//	}
//	else if (geometry==1){
			return sqrt((m1.x-m2.x)*(m1.x-m2.x)+(m1.y-m2.y)*(m1.y-m2.y)+(m1.z-m2.z)*(m1.z-m2.z));
//	}

}

inline double electrostatic_energy(molecule m1, molecule m2){

	if(m1.type==0 || m2.type==0){
		return 0.0;}
	else {
		return m1.q*m2.q*bjerrum/Distance(m1,m2);}


}

 inline double lj_energy(molecule m1, molecule m2){
	double r=Distance(m1,m2);
	if( ( m1.type == 0) || ( m2.type == 0) || ( r > lj_cut )   ){
		return 0.0;}
	else {

		double e= pow(lj_sigma/r,12)-pow(lj_sigma/r,6) + lj_shift;

		return 4.0*lj_eps*e;}
}


double constraint_energy(molecule m1){

	if (geometry==1){
	double r=L-m1.radius;
	if( ( r > m1.lj_cut_constraint ) ){
		return 0.0;}
	else if( ( r < -time_step) ){
		cout << "\n\n!!! Warning: atom out of box! Check the constraint." << endl;
		cout << "Atom # R = " << m1.radius <<"; xyz = " << m1.x <<" " << m1.y<<" " << m1.z << "; Type = " << m1.type << endl;
//		m1.x=0.1; m1.y=0.1;m1.z=0.1;
        m1.print();
		return 10000.;
	}
	else {
		if (m1.type==1){
//				set offset for anions
				r=r+r_off_anions;
                double e= pow(lj_sigma/r,12)-pow(lj_sigma/r,6) + lj_shift_anions ;
                return 4.0*m1.lj_eps_constraint*e;
                }
            else if (m1.type==2){
				r=r+0*r_off_anions;
                double e= pow(lj_sigma/r,12)-pow(lj_sigma/r,6) + lj_shift_cations ;
                return 4.0*m1.lj_eps_constraint*e;}
            else if (m1.type==0){
//                double e= pow(0.5*lj_sigma/r,12)-pow(0.5*lj_sigma/r,6) + lj_shift ;
                return 0.0;}
            else{
            	cout << "No such molecule type" << m1.type << endl;
            	return 10000.;
            }

		 }
	}
	else if (geometry==0){
		return 0.0;
//	if( ( m1.x > m1.lj_cut_constraint )&&( m1.y > m1.lj_cut_constraint )&&( m1.z > m1.lj_cut_constraint ) && ( L-m1.x > m1.lj_cut_constraint )&&(L- m1.y > m1.lj_cut_constraint )&&( L-m1.z > m1.lj_cut_constraint ) ) {
//		return 0.0;}
//	else {
//		double e=0.0;double r = L/2.;
////		r=m1.x ; e+= pow(lj_sigma/r,12)-pow(lj_sigma/r,6) + lj_shift_anions ;
////		r=m1.y ; e+= pow(lj_sigma/r,12)-pow(lj_sigma/r,6) + lj_shift_anions ;
////		r=m1.z ; e+= pow(lj_sigma/r,12)-pow(lj_sigma/r,6) + lj_shift_anions ;
////		r=L-m1.x ; e+= pow(lj_sigma/r,12)-pow(lj_sigma/r,6) + lj_shift_anions ;
////		r=L-m1.y ; e+= pow(lj_sigma/r,12)-pow(lj_sigma/r,6) + lj_shift_anions ;
////		r=L-m1.z ; e+= pow(lj_sigma/r,12)-pow(lj_sigma/r,6) + lj_shift_anions ;
//        return 4.0*m1.lj_eps_constraint*e;
//	}
	}
	else{
		cout << "unknown geometry! check the constraint"<< endl;
		return 1000.;
	}

}


double calc_energy(molecule M[]){
    double e=0.0;
    int size = N_molecules;
    for(int i=0;i<size-1;i++){
    molecule m1=M[i];
    e+=constraint_energy(m1);
    for (int j=i+1; j<size; j++){
        molecule m2=M[j];
        //coulomb interaction
        e+=electrostatic_energy(m1,m2);
        //Lennard-jones
        e+=lj_energy(m1,m2);
    }
}
    return e+constraint_energy(M[size-1]);
}

double recalc_energy(int i, molecule M[]){
    double e=0.0;
    int size = N_molecules;
    molecule m1=M[i];

    e+=constraint_energy(m1);
//    cout<< "recalc_energy: constraint energy = "<<e<< endl;
    for (int j=0; j<size; j++){
    	if(j!=i){
        molecule m2=M[j];
        //coulomb interaction
        e+=electrostatic_energy(m1,m2);
        //Lennard-jones
        e+=lj_energy(m1,m2);
    	}
    }
    return e;
}


int mc_step(molecule M[], int mol_ind){
//    double eold=calc_energy(M,size);
    double eold=recalc_energy(mol_ind,M);
    double x1=M[mol_ind].x;
    double y1=M[mol_ind].y;
    double z1=M[mol_ind].z;

    M[mol_ind].advance(1);
//    cout << "\n\nmc step index = " << mol_ind <<endl;
//    M[mol_ind].print();
//    double enew=calc_energy(M,size);
    double enew=recalc_energy(mol_ind,M);



    double prob = exp((-enew+eold)/kT);
    if(drand48() < prob){
    	return 1;}
    else {
    	M[mol_ind].move_to_position(x1,y1,z1);
//        cout << "probability " << drand48() << " and "<< exp(-enew+eold) << "; E new = "<< enew << " ; E old = "<< eold << endl;
    	return 0;	 }

}


double mc_steps(molecule M[],  int N){
	int success=0;
	int size = N_molecules;
//	for (int i=0;i<size;i++){
//    cout <<  "\n\n mc_steps"  << endl;
//    M[i].print();
//	}

    for (int i=0;i<N;i++){
        int index = (int)(drand48()*(double)(size));
//        cout <<  "index = "<<index << endl;
        success += mc_step(M, index);
//        if (success==1){
////            M[index].print();
//            }
           // cout << "index = " <<index;
    }

    return (double)success/N;

}


int gnuplot(molecule M[],int j){
	// Write a file
		int size = N_molecules;
	    ofstream positions_file;
        char str[80];
        sprintf(str, "/pos%d.dat",j);
//        printf(str);
	    positions_file.open( (foldername+str).c_str(), ofstream::out | ofstream::trunc );
	    for(int i=0;i<size;i++){
	    positions_file << M[i].x << "  " << M[i].y << "  " << M[i].z << "  " << M[i].radius << " " << M[i].q << endl;
	    }
	    positions_file.close();

//	    char * commandsForGnuplot[] = {"set title 'Graph'; set pointsize 0.5;set view equal xyz", " "};


	    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");

	fprintf(gnuplotPipe, "set terminal postscript eps enhanced color font 'Helvetica,20'\n");
	fprintf(gnuplotPipe, "set output '%s/3D_droplet%d.eps'\n",(foldername).c_str(),j );
	fprintf(gnuplotPipe, "set xlabel 'X '\n" );
	fprintf(gnuplotPipe, "set ylabel 'Y '\n" );
	fprintf(gnuplotPipe, "set zlabel 'Z '\n" );
	fprintf(gnuplotPipe, "set pointsize 1; \n" );
	fprintf(gnuplotPipe, "set ticslevel 0 \n" );
	fprintf(gnuplotPipe, "set title 'Graph #%d'; set pointsize 0.5; set view equal xyz \n",j);
	fprintf(gnuplotPipe, "splot '%s' u 1:2:3:5 w p palette pointtype 7 \n", (foldername+str).c_str());
    pclose(gnuplotPipe);
    return 0;
}

int gnuplot(polymer poly[],int j){
	// Write a file
		int size = N_polymers;
	    ofstream positions_file;
        char str[80];
        sprintf(str, "/polymer_pos%d.dat",j);
//        printf(str);
	    positions_file.open( (foldername+str).c_str(), ofstream::out | ofstream::trunc );
	    for(int i=0;i<size;i++)
	    {
            for(int j=0;j<poly[i].N;j++)
            {
                positions_file << poly[i].M[j].x << "  " << poly[i].M[j].y << "  " << poly[i].M[j].z << "  " << poly[i].M[j].radius  << endl;
            }
	    }
	    positions_file.close();

//	    char * commandsForGnuplot[] = {"set title 'Graph'; set pointsize 0.5;set view equal xyz", " "};


	    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");

	fprintf(gnuplotPipe, "set terminal postscript eps enhanced color font 'Helvetica,20'\n");
	fprintf(gnuplotPipe, "set output '%s/3D_droplet%d.eps'\n",(foldername).c_str(),j );
	fprintf(gnuplotPipe, "set xlabel 'X '\n" );
	fprintf(gnuplotPipe, "set ylabel 'Y '\n" );
	fprintf(gnuplotPipe, "set zlabel 'Z '\n" );
	fprintf(gnuplotPipe, "set pointsize 1; \n" );
	fprintf(gnuplotPipe, "set ticslevel 0 \n" );
	fprintf(gnuplotPipe, "set title 'Graph #%d'; set pointsize 0.5; set view equal xyz \n",j);
	fprintf(gnuplotPipe, "splot '%s' u 1:2:3 w lp palette pointtype 7 \n", (foldername+str).c_str());
    pclose(gnuplotPipe);
    return 0;
}

//int molecule::id = 0;


int radial_distribution(molecule M[], int counts[], int type){
//	int counts[nbins];
//	fill(counts, counts+nbins, 0);
	int size = N_molecules;
	for (int i=0;i<size;i++){
	if (type==M[i].type){
		if (M[i].radius>L){
			cout << "Error: Particle out of box" << endl;
			cout << "scaled radii = " << M[i].radius << " to 100 = "<<  (int) (M[i].radius/L*nbins) << endl;}

		counts[(int) (nbins*M[i].radius/(R)) ]++;}
	}
//	for (int i = nbins - 1; i >= 0; i--){
//	    cout << counts[i] << " " << i << endl;}
	return 0;
}


int axis_distribution(molecule M[], int counts[], int type, int axis){
		int size = N_molecules;
		for (int i=0;i<size;i++){

			if (type==M[i].type){
				if(axis==0){counts[(int) (nbins*M[i].x/L) ]++;}
				if(axis==1){counts[(int) (nbins*M[i].y/L) ]++;}
				if(axis==2){counts[(int) (nbins*M[i].z/L) ]++;}
			}
		}
		return 0;
}

void plot_radius(int  counts[], int j){
	// Write a file
	ofstream radius_file;
    char str[80];
    sprintf(str,"/radius_prob%d.dat",j);
//    printf(str);
	radius_file.open( (foldername+str).c_str() , ofstream::out | ofstream::trunc);
	double dr = R/nbins;
	for(int i=0;i<nbins;i++){
	   	radius_file << (double)counts[i]/4.0/M_PI/((1.0+i)*dr)/((1.0+i)*dr)/(dr)/(j+1) << "  " << (double)(1.0+i)*dr  << endl;
	    }
    radius_file.close();

	//	    char * commandsForGnuplot[] = {"set title 'Graph'; set pointsize 0.5;set view equal xyz", " "};


		    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");

		fprintf(gnuplotPipe, "set terminal postscript eps enhanced color font 'Helvetica,20'\n");
		fprintf(gnuplotPipe, "set output 'radial%d.eps'\n",j );
		fprintf(gnuplotPipe, "set xlabel 'X '\n" );
		fprintf(gnuplotPipe, "set ylabel 'Y '\n" );
		fprintf(gnuplotPipe, "set ticslevel 0 \n" );
		fprintf(gnuplotPipe, "set title 'Radial Distribution'; set pointsize 0.70; \n");
		fprintf(gnuplotPipe, "plot '%s' u 2:1 w l lt 1 lc -1,'%s' u 2:1 w p pt 7 lc 1  \n",(foldername+str).c_str(),(foldername+str).c_str() );

}
// make some return values for analysis functions using bulk concentration or charge?
double plot_radius(int  anion[], int cation[], int j){
	// Write a file
	ofstream radius_file;
    char str[80];
    sprintf(str,"/radius_prob%d.dat",N_molecules/2);
//    printf(str);
	radius_file.open( (foldername+str).c_str() , ofstream::out | ofstream::trunc);
	double dr = R/nbins; int charge1=0; int charge2=0;
	for(int i=0;i<nbins;i++){
	   	radius_file <<  (double)(1.0+i)*dr << "  " << (double)anion[i]/4.0/M_PI/((1.0+i)*dr)/((1.0+i)*dr)/dr/(j+1) << "  " << (double)cation[i]/4.0/M_PI/((1.0+i)*dr)/((1.0+i)*dr)/dr/(j+1)  << endl;
	   	if (R-((1.0+i)*dr) < 0.5){
	   		charge1+=anion[i]-cation[i];
	   	} else {
	   		charge2+=anion[i]-cation[i];
	   	}
	 }
    radius_file.close();

//    cout << "Distribution, charges. Q1 = " << (double)charge1/(j+1) << "; Q2 = " << (double)charge2/(j+1) << endl;
	//	    char * commandsForGnuplot[] = {"set title 'Graph'; set pointsize 0.5;set view equal xyz", " "};


		    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");

		fprintf(gnuplotPipe, "set terminal postscript eps enhanced color font 'Helvetica,20'\n");
		fprintf(gnuplotPipe, "set output '%s/radial%d.eps'\n",(foldername).c_str(),N_molecules/2 );
		fprintf(gnuplotPipe, "set xlabel 'r/nm '\n" );
		fprintf(gnuplotPipe, "set ylabel 'c nm^3 '\n" );
		fprintf(gnuplotPipe, "set ticslevel 0 \n" );
		fprintf(gnuplotPipe, "set title 'Radial Distribution'; set pointsize 0.40; \n");
		fprintf(gnuplotPipe, "plot[][0:3*%f]'%s' u 1:2 w lp pt 7 lt 1 lc -1 t 'anion density','%s' u 1:3 w lp lt 1 pt 7 lc 1 t 'cation density', '' u 1:(sqrt($2*$3)) w l lt 2 lc 3 t 'bulk conc' \n",(double)N_molecules/V/2.,(foldername+str).c_str(),(foldername+str).c_str() );
		fprintf(gnuplotPipe, "plot[:%f][] '%s' u 1:2 w lp pt 7 lt 1 lc -1 t 'anion density','%s' u 1:3 w lp lt 1 pt 7 lc 1 t 'cation density', '' u 1:(sqrt($2*$3)) w l lt 2 lc 3 t 'bulk conc' \n",L-lj_cut,(foldername+str).c_str(),(foldername+str).c_str() );
		pclose(gnuplotPipe);

	return -(double)charge1/(j+1);
}

//void plot_file(string filename){
//	// Write a file
////        char str[200];
////        sprintf(str," grep chemical %s | sed 's/;/ /g' | awk '{print $6,$9,$13,$17,$20}' > chem_pot.dat", (foldername+"/"+"output.dat").c_str());
//        system((" grep chemical " + foldername+"/"+ "output.dat" + " | sed 's/;/ /g' | awk '{print $6,$9,$13,$17,$20}' > chem_pot.dat").c_str());
//	    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
//
//		fprintf(gnuplotPipe, "set terminal postscript eps enhanced color font 'Helvetica,20'\n");
//		fprintf(gnuplotPipe, "set output '%s/%s.eps'\n",(foldername).c_str(),(filename).c_str() );
//		fprintf(gnuplotPipe, "set xlabel 'N'\n" );
//		fprintf(gnuplotPipe, "set ylabel  'Free energy' offset 2 \n" );
//        fprintf(gnuplotPipe, "set y2label 'Chemical potential'  offset -1 \n" );
//		fprintf(gnuplotPipe, "set ytics nomirror; set y2tics nomirror;y=0;x1=0;x2=0 \n" );
//		fprintf(gnuplotPipe, "set title 'Chemical potential'; set pointsize 0.40;V=%f**3.0 \n", L);
//		fprintf(gnuplotPipe, "plot[][] 'chem_pot.data' u 1:3 w p pt 7 t '-k(N-N0)+{/Symbol m}' axes x1y2, '' u 1:(x2=$1,y=y+(x2-x1)*$3,x1=$1,y  ) w lp pt 7 lc -1  t 'Free energy F(N)',  log(x/V/2./0.001)  t 'ideal gas'  axes x1y2 \n" );
//		pclose(gnuplotPipe);
//}

void plot_axis(int  anion[], int cation[], int j){

	// Write a file
	ofstream axis_file;
    char str[80];
    sprintf(str,"/axis_prob%d.dat",j);

    axis_file.open( (foldername+str).c_str() , ofstream::out | ofstream::trunc);
	double dx = L/nbins;double Lz=L; double Ly=L;
	for(int i=0;i<nbins;i++){
		axis_file <<  (double)(1.0+i)*dx << "  " << (double)anion[i]/Lz/Ly/dx/(j+1)<< "  " << (double)cation[i]/Lz/Ly/dx/(j+1)  << endl;
	    }
	axis_file.close();

	    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");

		fprintf(gnuplotPipe, "set terminal postscript eps enhanced color font 'Helvetica,20'\n");
		fprintf(gnuplotPipe, "set output '%s/axis%d.eps'\n",(foldername).c_str(),j );
		fprintf(gnuplotPipe, "set xlabel 'x/nm '\n" );
		fprintf(gnuplotPipe, "set ylabel 'c nm^3 '\n" );
		fprintf(gnuplotPipe, "set ticslevel 0 \n" );
		fprintf(gnuplotPipe, "set title 'Ionic profiles'; set pointsize 0.40; \n");
		fprintf(gnuplotPipe, "plot[][0:3*%f] '%s' u 1:2 w lp pt 7 lt 1 lc -1 t 'anion density','%s' u 1:3 w lp lt 1 pt 7 lc 1 t 'cation density', '' u 1:(sqrt($2*$3)) w l lt 2 lc 3 t 'bulk conc' \n",(double)N_molecules/V/2.,(foldername+str).c_str(),(foldername+str).c_str() );
		fprintf(gnuplotPipe, "plot[:%f][] '%s' u 1:2 w lp pt 7 lt 1 lc -1 t 'anion density','%s' u 1:3 w lp lt 1 pt 7 lc 1 t 'cation density', '' u 1:(sqrt($2*$3)) w l lt 2 lc 3 t 'bulk conc' \n",L-lj_cut,(foldername+str).c_str(),(foldername+str).c_str() );
		pclose(gnuplotPipe);
}


void plot_axis(int  counts[], int j){
	// Write a file
	ofstream axis_file;
    char str[80];
    sprintf(str,"/axis_prob%d.dat",j);
//    printf(str);
    axis_file.open( (foldername+str).c_str() , ofstream::out | ofstream::trunc);
    double dx = L/nbins;double Lz=L; double Ly=L;
	for(int i=0;i<nbins;i++){
		axis_file << (double)(1.0+i)*dx << "  " << (double)counts[i]/Lz/Ly/dx/(j+1)<< "  " << endl;
	    }
	axis_file.close();

	//	    char * commandsForGnuplot[] = {"set title 'Graph'; set pointsize 0.5;set view equal xyz", " "};


		    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");

		fprintf(gnuplotPipe, "set terminal postscript eps enhanced color font 'Helvetica,20'\n");
		fprintf(gnuplotPipe, "set output 'axis%d.eps'\n",j );
		fprintf(gnuplotPipe, "set xlabel 'x/nm '\n" );
		fprintf(gnuplotPipe, "set ylabel 'c nm^3 '\n" );
		fprintf(gnuplotPipe, "set ticslevel 0 \n" );
		fprintf(gnuplotPipe, "set title 'Radial Distribution'; set pointsize 0.70; \n");
		fprintf(gnuplotPipe, "plot[][0:3*%f] '%s' u 2:1 w l lt 1 lc -1,'%s' u 2:1 w p pt 7 lc 1  \n",(double)N_molecules/V/2.,(foldername+str).c_str(),(foldername+str).c_str() );

}




void return_molecules(molecule M[]){
	int size = N_molecules;
    cout << "\t\tReturning molecule ... "  << endl;
    int j=0;
    for (int i=0; i<size;i++){

    	if (geometry==1){
	   		double r=M[i].radius;
    		if (r > L){
    			j++;
    			r=(L-2.)*pow(drand48(),1./3.);
    			double theta=acos(2*(drand48()-0.5));
            	double phi=2.0*M_PI*drand48();
            	double x1=r*sin(theta)*cos(phi);
            	double y1=r*sin(theta)*sin(phi);
            	double z1=r*cos(theta);
            	M[i].move_to_position(x1,y1,z1);
            	cout << "\t\tSpherical geometry! Particle fixed! Rold = " << r << "; Rnew = " <<  M[i].radius << endl;
    		}
    	}
    	else if (geometry == 0){
    		if((M[i].x >L)||(M[i].y >L)||(M[i].z >L)||(M[i].x <-0.01)||(M[i].y <-0.01)||(M[i].z <-0.01)){
    			cout << "\t\t Cubic geometry! Particle fixed! xyz old = " <<  M[i].x<< " ;"<< M[i].y << " ;" << M[i].z   << endl;
    			j++;
    	    	double x1=drand48()*(L);
    	    	double y1=drand48()*(L);
    	    	double z1=drand48()*(L);
    			M[i].move_to_position(x1,y1,z1);
    			cout << "\t\t Cubic geometry! Particle fixed! xyz new = " <<  M[i].x<< " ;"<< M[i].y << " ;" << M[i].z   << endl;
    		}
    	}
    }
    cout << "\t\t number of returned molecules = " << j<<endl;
}

molecule delete_molecule(molecule M[],int type){

	if (N_molecules<2){
		return M[0];
	}

    int index = (int)(drand48()*(double)(N_molecules));
    while(M[index].type!=type){
        index = (int)(drand48()*(double)(N_molecules));
    }
    molecule deleted_mol = M[index];
    if (index==N_molecules-1){
        M[index].type=0;M[index].q=0.0;}
    else{
        M[index]=M[N_molecules-1];
        M[N_molecules-1] = deleted_mol;
        M[N_molecules-1].type=0;M[N_molecules-1].q=0.0;}
    N_molecules--;
    return deleted_mol;
}

molecule delete_molecule(molecule M[],int type,int id){

    molecule deleted_mol = M[id];

    if(M[id].type!=type){
        cout<< "Deleting molecule. Warning! inconsistent types."<<endl;
        cout<< "actual type = " << M[id].type << " assumed type = " << type<<endl;
    }

    if (id==N_molecules-1){
        M[id].type=0;M[id].q=0.0;}
    else{
        M[id]=M[N_molecules-1];
        M[N_molecules-1] = deleted_mol;
        M[N_molecules-1].type=0;M[N_molecules-1].q=0.0;}
    N_molecules--;
    return deleted_mol;
}

void add_molecule(molecule M[], molecule molecule_to_add){
    M[N_molecules] = molecule_to_add;
    N_molecules++;
}

molecule add_molecule(molecule M[], int type){
    molecule molecule_to_add;
    if(type==1){
        molecule_to_add.q=-1.0;
        molecule_to_add.type=1;
        molecule_to_add.lj_cut_constraint=lj_cut_constraint_anions;
        molecule_to_add.lj_eps_constraint=lj_eps_constraint_anions;
    }
    else if (type==2){
        molecule_to_add.q=1.0;
        molecule_to_add.type=2;
    } else{
     cout << "Adding a molecule. Warning! trying to add unknown or zero type molecule"<<endl;
    }

    M[N_molecules] = molecule_to_add;
    N_molecules++;

    return molecule_to_add;
}

int mc_delete_ions_step(molecule M[], double mu /* = 0.0 */ , double N0 /* = N_molecules/2 */ ){

	if (N_molecules<2){
		return 0;
	} else{

    int index1 = (int)(drand48()*(double)(N_molecules));
    while(M[index1].type!=1){
        index1 = (int)(drand48()*(double)(N_molecules));
    }
    int index2 = (int)(drand48()*(double)(N_molecules));
    while(M[index2].type!=2){
        index2 = (int)(drand48()*(double)(N_molecules));
    }
//    cout << "print indices : " << index1 << " ; " << index2 << " ;Nmol " << N_molecules <<endl;
    double eold=recalc_energy(index1,M)+recalc_energy(index2,M)+0.5*k*(N0-N_molecules/2)*(N0-N_molecules/2);
    double excess_energy=electrostatic_energy(M[index1],M[index2])+lj_energy(M[index1],M[index2]);
    eold-=excess_energy;
//    delete pair of molecules
    molecule cation  = delete_molecule(M,2,index2);
    molecule anion = delete_molecule(M,1,index1);
//    cout << "print indices : " << index1 << " ; " << index2 << " ;Nmol " << N_molecules <<endl;

    double enew=0.5*k*(N0-N_molecules/2)*(N0-N_molecules/2);

    double prob = exp(-2*mu+2*log((double)(N_molecules+2.0)/2.0*pow(lambda,3)/V)+(-enew+eold)/kT);
    if(drand48() < prob){
    	return 1;}
    else {
    	add_molecule(M, anion);
    	add_molecule(M,cation);
    	return 0;	 }
	}
}

double mc_delete_ions_steps(molecule M[],double mu /* = 0.0 */ , double N0 /* =N_molecules/2 */,  int N /* =1 */ ){
	int success=0;

    for (int i=0;i<N;i++){
        success += mc_delete_ions_step(M, mu,N0);
    }

    return (double)success/N;

}


int mc_add_ions_step(molecule M[], double mu /* = 0.0 */ ,double N0 /* = N_molecules/2 */){

    double eold=0.5*k*(N0-N_molecules/2)*(N0-N_molecules/2);
//    add pair of molecules
    add_molecule(M,1);
    add_molecule(M,2);

    double enew=recalc_energy(N_molecules-2,M)+recalc_energy(N_molecules-1,M)+0.5*k*(N0-N_molecules/2)*(N0-N_molecules/2);
    double excess_energy=electrostatic_energy(M[N_molecules-2],M[N_molecules-1])+lj_energy(M[N_molecules-2],M[N_molecules-1]);
    enew-=excess_energy;



    double prob = exp(2*mu-2*log((double)(N_molecules)/2.0*pow(lambda,3)/V)+(-enew+eold)/kT);
    if(drand48() < prob){
    	return 1;}
    else {
    	delete_molecule(M,2,N_molecules-1);
    	delete_molecule(M,1, N_molecules-1);

    	return 0;	 }
}

double mc_add_ions_steps(molecule M[],double mu /* = 0.0 */, double N0 /* = N_molecules/2*/,  int N /*=1 */){
	int success=0;

    for (int i=0;i<N;i++){
        success += mc_add_ions_step(M, mu,N0);
    }

    return (double)success/N;

}


double bulk_conc(double mu){

//	 other dependencies may be introduced here

	 //concentration of ideal gas with zero chem. potential
	 double c0=0.001;

//	 spherical bulk conc

	 return c0*exp(mu)*(0.0121086*mu+1.0) -0.000190477*mu-0.000206683;

//			 planar geometry bulk
//return c0*exp(mu/kT)*(1.0+0.0449513*mu);
}

double Gamma(double mu,double N_av /* = (double)N_molecules */ ){

	 // other dependencies may be introduced here

	 double cb = bulk_conc(mu);

	 return N_av/(4.*M_PI*R*R)-cb*R/3.0;

}


 double mc_any_steps(molecule M[], double mu /* = 0.0 */ , double N0 /* = N_molecules*1.0/2.0 */ , int N /* = 1 */ ){
	double acceptance_advance=0.0;
	double acceptance_delete=0.0;
	double acceptance_add=0.0;
	int trials_advance=0;int trials_delete=0;int trials_add=0;
	long int Nmol=0;
	long int NmolSq=0;
	for (int i=0;i<N;i++){
		if(drand48()<0.5){
			acceptance_advance+=mc_steps(M,10);
			trials_advance++;
		} else{

		if(drand48()<0.5){
			acceptance_delete+=mc_delete_ions_steps(M,mu,N0,2);
			trials_delete++;
			Nmol+=N_molecules;
			NmolSq+=N_molecules*N_molecules;

		}
		else{
			acceptance_add+=mc_add_ions_steps(M,mu,N0,2);
			trials_add++;
			Nmol+=N_molecules;
			NmolSq+=N_molecules*N_molecules;
		}
		}

	}
	double N_av=(double)Nmol/(trials_add+trials_delete)/2.0;double NSq_av=(double)NmolSq/(trials_add+trials_delete)/4.0;
	double dN=sqrt(NSq_av-N_av*N_av)*1.96/sqrt((trials_add+trials_delete));
	double mu_mes = mu-0.5*k*( N_av - N0 );
	double cb = bulk_conc(mu_mes) ;
	cout << "Advance trials = " << trials_advance<< " Delete trials = " << trials_delete  << " Addition trials = " << trials_add << "; dN = " << dN << "; dmu = " << 0.5*k*dN <<endl;
	cout << "Advance acceptance = " << acceptance_advance*2./N << " Delete acceptance = " << acceptance_delete*4./N  << " Addition acceptance = " << acceptance_add*4.0/N <<endl;
	cout << "Average number of molecules = " << N_av  << "; N0 = " << N0 << "; chemical potential = " << mu_mes << " ; Ideal gas = "<< log(N_av*pow(lambda,3)/V) << ";Conc = "<< N_av/V << "; lmabda="<< 1.0/sqrt(8.0*M_PI*bjerrum*N_av/V) << " Cb = " << cb << " ; Ideal gas correct = "<< log(cb*pow(lambda,3)) << " Gamma = " << Gamma(mu_mes, N_av) << endl;
	return N_av ;
}



 double bulk_conc2(double mu){

 //	 other dependencies may be introduced here

 	 //concentration of ideal gas with zero chem. potential
 	 double c0=0.001;

 //	 spherical bulk conc

 	 return c0*exp(mu)*(0.0121086*mu+1.0) -0.000190477*mu-0.000206683;

 //			 planar geometry bulk
 //return c0*exp(mu/kT)*(1.0+0.0449513*mu);
 }

 double Gamma2(double mu, double N_av /* = (double)N_molecules */ ){

 	 // other dependencies may be introduced here

 	 double cb = bulk_conc(mu);

 	 return N_av/(4.*M_PI*R*R)-cb*R/3.0;

 }
