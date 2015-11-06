#include "../include/sim_system.h"
#include "../global.h"
#include "../molecule.h"
#include "../include/polymer.h"
#include "../include/integrate.h"

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

#define NUM_COMMANDS 2
#define MAX_PART 5000
// #define MAX_MOL  1000
// #define MAX_POL  100

//max_part is the maximum number of simple single particles


using namespace std;

sim_system::sim_system()
{
    M = NULL;
    poly = NULL;
    other = NULL;
    // change list into list of references
    mol_list = new molecule * [MAX_PART];  //molecule[MAX_PART];
    temperature = kT;
    box         = L;
    volume      = V;
    max_part_id = 0;

    //bond_list;
    //verlet_lsit;

}


// make list of pointer to refer to main particle pointer M and poly

int sim_system::create_particle_list()
{
    max_part_id = 0;
    // for non-bonded interaction energy
    cout << " ## Creating molecules list ... " << endl;

    for (int i = 0; i < N_molecules; ++i)
    {
        //mol_list[i] = new molecule;
        mol_list[i] = &M[i];
        max_part_id ++;
    }

    cout << "Molecule list created. number of simple molecules = " << max_part_id << " ; Nmol = " << N_molecules << endl;
    cout << " ## Creating polymers list ... " << endl;
    for (int i = 0; i < N_polymers; i++)
    {
        cout << " \t polymer number  #" << i << " with " << poly[i].N << " molecules"<< endl;
        for (int j = 0; j < poly[i].N; ++j)
        {
            cout << "mol_list = " << max_part_id << endl;
            //mol_list[max_part_id] = new molecule;
//            molecule m_tmp = poly[i].M[j];
            mol_list[max_part_id] = &(poly[i].M[j]);
            max_part_id ++;
        }

    }

    cout << " List created. max_part_id = " << max_part_id << endl;

    return 0;
}


    //  takes a molecule from a list and calculates its energy
//  energy for a simple mol move, no donded part.

double sim_system::recalc_energy_mol(int i)
{
    double e=0.0;
    molecule m1 = *mol_list[i];

    e+=constraint_energy(m1);
//    cout<< "recalc_energy: constraint energy = "<<e<< endl;

    // calculate non-bonded part
    for (int j=0; j<max_part_id; j++){
    	if(j!=i){
        molecule m2 = *mol_list[j];
        //coulomb interaction
        e+=electrostatic_energy(m1,m2);
        //Lennard-jones
        e+=lj_energy(m1,m2);
    	}
    }
    return e;
}

// here one can include different energies to penalise some system moves
// e.g. harmonic potentials on some molecules
double sim_system::extra_energy_pol()
{
    double e = 0.0;
    // an example of constraint
    //springs on polymers' centers of masses.
    poly[0].update_COM();
    e += 100.0*(poly[0].xc)*(poly[0].xc)+(poly[0].yc)*(poly[0].yc)+(poly[0].zc+7.0)*(poly[0].zc+7.0);
    poly[1].update_COM();
    e += 100.0*( poly[1].xc)*( poly[1].xc)+( poly[1].yc)*( poly[1].yc)+( poly[1].zc-7.0)*( poly[1].zc-7.0);
    poly[2].update_COM();
    e += 100.0*(poly[2].xc-7.)*(poly[2].xc-7.)+(poly[2].yc)*(poly[2].yc)+(poly[2].zc)*(poly[2].zc);
    poly[3].update_COM();
    e += 100.0*(poly[3].xc+7.)*(poly[3].xc+7.)+(poly[3].yc)*(poly[3].yc)+(poly[3].zc)*(poly[3].zc);
//    cout << "error here? " << e << endl;
    return e;
}

// energy for a polymer bead move
double sim_system::recalc_energy_pol(int pol_i, int mol_i)
{

    double e=0.0;
    molecule m1 = (poly[pol_i].M[mol_i]);

    e+=constraint_energy(m1);
//    cout<< "recalc_energy: constraint energy = "<<e<< endl;

    // calculate non-bonded part
    for (int j = 0; j < max_part_id; j++){
        molecule m2 = *mol_list[j];
        if(m1.id != m2.id){

        //coulomb interaction
        e += electrostatic_energy(m1,m2);
        //Lennard-jones
        e += lj_energy(m1,m2);
        }
    }

    // go through bonded ... to be done

    e += poly[pol_i].bond_energy(mol_i);

    // we can add some extra energies if necessary, e.g.
    e += sim_system::extra_energy_pol();


//    for (int i = 0; i < N_polymers; i++){
//        molecule m = poly[i].M[0];
//        e += 100.0*(m.x)*(m.x)+(m.y)*(m.y)+(m.z)*(m.z);
//    }
    return e;
}


double sim_system::calc_total_energy()
{
    double e = 0.0;
    // nonbonded part
    for(int i=0;i<max_part_id-1;i++){
        molecule m1=*mol_list[i];
        e+=constraint_energy(m1);

        for (int j = i+1; j < max_part_id; j++){
            molecule m2 = *mol_list[j];
            //coulomb interaction
            e += electrostatic_energy(m1,m2);
            //Lennard-jones
            e += lj_energy(m1,m2);
        }
    }

    e += constraint_energy(*mol_list[max_part_id-1]);
    // bobnded part of energy
    for (int i = 0; i < N_polymers; i++)
    {
        e += poly[i].bond_energy();
    }
    return e;
}


int sim_system::mc_step_mol( int mol_ind){
//    double eold=calc_energy(M,size);
    double eold= recalc_energy_mol(mol_ind);
    double x1 = M[mol_ind].x;
    double y1 = M[mol_ind].y;
    double z1 = M[mol_ind].z;

    M[mol_ind].advance(1);
//    cout << "\n\nmc step index = " << mol_ind <<endl;
//    M[mol_ind].print();
//    double enew=calc_energy(M,size);
    double enew = recalc_energy_mol(mol_ind);



    double prob = exp((-enew+eold)/kT);
    if(drand48() < prob){
        return 1;}
    else {
        M[mol_ind].move_to_position(x1,y1,z1);
//        cout << "probability " << drand48() << " and "<< exp(-enew+eold) << "; E new = "<< enew << " ; E old = "<< eold << endl;
        return 0;    }

}


double sim_system::mc_steps_mol( int Nsteps){
    int success=0;

//  for (int i=0;i<size;i++){
//    cout <<  "\n\n mc_steps"  << endl;
//    M[i].print();
//  }

    for (int i = 0; i < Nsteps; i++){
        int index = (int)(drand48()*(double)(N_molecules));
//        cout <<  "index = "<<index << endl;
        success += mc_step_mol( index);
//        if (success==1){
////            M[index].print();
//            }
           // cout << "index = " <<index;
    }

    return (double)success/Nsteps;

}



int sim_system::mc_step_pol( int pol_ind){
//    double eold=calc_energy(M,size);
    int mol_ind = (int)(drand48()*(double)(poly[pol_ind].N));
    double eold = recalc_energy_pol(pol_ind,mol_ind);

    double x1   = poly[pol_ind].M[mol_ind].x;
    double y1   = poly[pol_ind].M[mol_ind].y;
    double z1   = poly[pol_ind].M[mol_ind].z;

    poly[pol_ind].M[mol_ind].advance(1);
//    cout << "\n\nmc step index = " << mol_ind <<endl;
//    M[mol_ind].print();
//    double enew=calc_energy(M,size);
    double enew = recalc_energy_pol(pol_ind, mol_ind);



    double prob = exp((-enew+eold)/kT);
    if(drand48() < prob){
        return 1;}
    else {
        poly[pol_ind].M[mol_ind].move_to_position(x1,y1,z1);
//        cout << "probability " << drand48() << " and "<< exp(-enew+eold) << "; E new = "<< enew << " ; E old = "<< eold << endl;
        return 0;    }

}


double sim_system::mc_steps_pol( int Nsteps){
    int success=0;

//  for (int i=0;i<size;i++){
//    cout <<  "\n\n mc_steps"  << endl;
//    M[i].print();
//  }

    for (int i = 0; i < Nsteps; i++){
        int index = (int)(drand48()*(double)(N_polymers));
//        cout <<  "index = "<<index << endl;
        success += mc_step_pol( index);
//        if (success==1){
////            M[index].print();
//            }
           // cout << "index = " <<index;
    }

    return (double)success/Nsteps;

}






int sim_system::gnuplot(int j /*, char * file_format*/){
	// Write a file

	    ofstream positions_file1;
        char str1[80];
        sprintf(str1, "/pos%d.dat",j);
//        printf(str);
	    positions_file1.open( (foldername+str1).c_str(), ofstream::out | ofstream::trunc );
	    for(int i=0; i < N_molecules;i++){
	    positions_file1 << M[i].x << "  " << M[i].y << "  " << M[i].z << "  " << (M[i].q+1)/2*N_polymers << endl;
	    }
	    positions_file1.close();


	    //ofstream positions_file2;
        char str2[80];
        sprintf(str2, "/polymer_pos%d.pdb",j);
//        printf(str);
	    //positions_file2.open( (foldername+str2).c_str(), ofstream::out | ofstream::trunc );
	    FILE * positions_file2;
	    positions_file2 = fopen ((foldername+str2).c_str(),"w");
	    int cnt = 1;
	    fprintf(positions_file2,"REMARK polymer configuration in MC code\n");

	    for(int i=0; i < 40;i++)
	    {
            for(int j=0; j < 30;j++)
                {
            //double theta=acos(2*(drand48()-0.5));
            double theta=acos(-1.0+i*2.0/40.);
            //double phi=M_PI*drand48();
            double phi=M_PI*j*1.0/30.;
            double x = 10*sin(theta)*cos(phi);
            double y = 10*sin(theta)*sin(phi);
            double z = 10*cos(theta);
                //positions_file2 << poly[i].M[j].x << "  " << poly[i].M[j].y << "  " << poly[i].M[j].z << "  " << poly[i].id  << endl;
            fprintf(positions_file2,"ATOM %6d%4s  UNX F%4d    %8.3f%8.3f%8.3f  0.00  0.00      T%03d \n",cnt, "FE",cnt-1,x,y,z , 10 );
            cnt++;

	    }}

	    	    for(int i=0; i < N_molecules;i++)
	    {
                //positions_file2 << poly[i].M[j].x << "  " << poly[i].M[j].y << "  " << poly[i].M[j].z << "  " << poly[i].id  << endl;
            fprintf(positions_file2,"ATOM %6d%4s  UNX F%4d    %8.3f%8.3f%8.3f  0.00  0.00      T%03d \n",cnt, "FE",cnt-1,M[i].x,M[i].y,M[i].z , int(M[i].q/2.+1.0) );
            cnt++;

	    }

	    for(int i=0; i < N_polymers; i++)
	    {
            for(int j=0;j<poly[i].N;j++)
            {
                //positions_file2 << poly[i].M[j].x << "  " << poly[i].M[j].y << "  " << poly[i].M[j].z << "  " << poly[i].id  << endl;
                fprintf(positions_file2,"ATOM %6d%4s  UNX F%4d    %8.3f%8.3f%8.3f  0.00  0.00      T%03d \n",cnt, "FE",cnt-1,poly[i].M[j].x,poly[i].M[j].y,poly[i].M[j].z ,  1+poly[i].id);
                cnt++;
            }
	    }
	    //positions_file2.close();
	    fclose(positions_file2);
//	    char * commandsForGnuplot[] = {"set title 'Graph'; set pointsize 0.5;set view equal xyz", " "};


	    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");

	//fprintf(gnuplotPipe, "set terminal postscript eps enhanced color font 'Helvetica,20'\n");
	fprintf(gnuplotPipe, "set term pngcairo size 1024,1024\n");
	fprintf(gnuplotPipe, "set output '%s/3D_droplet%d.png'\n",(foldername).c_str(),j );
	fprintf(gnuplotPipe, "set xlabel 'X '; \n" );
	fprintf(gnuplotPipe, "set ylabel 'Y '\n" );
	fprintf(gnuplotPipe, "set zlabel 'Z '\n" );
	fprintf(gnuplotPipe, "set pointsize 1; \n" );
	fprintf(gnuplotPipe, "set ticslevel 0 \n" );
	fprintf(gnuplotPipe, "set title 'Graph #%d'; set pointsize 0.5; set view equal xyz; \n set palette rgb 33,13,10; \n",j);
	fprintf(gnuplotPipe, "splot '%s' u 1:2:3:4 w p palette pointtype 7  ps 1.,\\\n", (foldername+str1).c_str());
	fprintf(gnuplotPipe, "'%s' u 1:2:3:4 w p palette pt 7 ps 4.0, '' u 1:2:3 w p  pt 6 lc rgb 'gray' lt 2 ps 4.0  \n", (foldername+str2).c_str());
    pclose(gnuplotPipe);
    return 0;
}



sim_system::~sim_system()
{
    //dtor
//    for (int i = 0; i < MAX_PART ; i++)
//    {
//        delete mol_list[i];
//    }
    delete[] mol_list;
}
