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

#include "molecule.h"
#include "global.h"
#include "include/integrate.h"
#include "include/polymer.h"
#include "include/sim_system.h"

using namespace std;

//time_step
extern string foldername="polymer_";
const double time_step = 0.5;
const int nbins=200;
//Droplet size
const int geometry=1;  // 0 if cubic , 1 if spherical
const double L = 10.15; // Nanometers
const double R = L-0.15;
const double V=4./3.*M_PI*pow(R,3); // for sphere 4./3.*M_PI*pow(L-2.0,3)

//Electrostatic parameter
const double bjerrum = 0.7;
//Lennard-jones parameters
const double lj_eps=1.4;
const double lj_sigma=0.5;            // Sodium atomic radius 2.5 A
const double lj_shift=0.25;
const double lj_cut=1.12246204831*lj_sigma;
//const double lj_cut=2.5;
//Parameters for anions
const int zero_adsorbption=1;
const double lj_eps_constraint_anions=4.0;
const double r_off_anions=0.5*lj_sigma;
//const double lj_shift_anions=0.25;
//const double lj_cut_constraint_anions=1.12246204831*lj_sigma;
const double lj_shift_anions=0.004079223;
const double lj_cut_constraint_anions=2.*lj_sigma;
const double lj_shift_cations=0.25;

const double kT=1.0;
const double lambda=10.0;
double k=2.0;


# Number of polymers and simple (single) molecules
int N_molecules = 200;
int N_polymers  = 4;
int mol_id      = 0;
int pol_id      = 0;
int main(int argc, char* argv[]) {

// initial and final number of molecues, number of steps
	int n0=102;int n_steps=500; int n_end=602;
	if (argc>1){n0=atoi(argv[1]);N_molecules=2*n0+2;cout<< " Start with N0 = "<<n0<<endl;}
	if (argc>2){n_end=atoi(argv[2]);cout<< " End with N_end = "<<n_end<<endl;}
	if (argc>3){n_steps=atoi(argv[3]);cout<< " with  "<<n_steps<< " steps"<<endl;}

	cout << "argc = " << argc << endl;
	for(int i = 0; i < argc; i++)
		cout << "argv[" << i << "] = " << argv[i] << endl;

//  ofstream output;
//	output.open( "output.dat" , ofstream::out | ofstream::trunc);
    char str0[80];
    sprintf(str0,"N0_%d_End_%d",n0,n_end);
    foldername=foldername+str0;
    mkdir((foldername).c_str(),0777);
	freopen((foldername+"/output.dat").c_str(),"w",stdout);
//  random number initialisator

    srand (time(0));
    srand48(time(0));
//    srand (1);
//    srand48(1);

    cout << "########################################" << endl;
    cout << "# Surface tension of a charged droplet #" << endl;
    cout << "# Bulk properties #" << endl;
    cout << "########################################" << endl;
    cout << "\nComputer simulation" << endl;
    cout << "\tSimulation parameters: " << endl;
    cout << "\t\tBjerrum length = " << bjerrum << " nm" << endl;
    cout << "\t\tBox size = " << L << " nm; Droplet volume =  " << V << " nm^3 ;"<< endl;
    cout << "\t\tNumber of particles = " << N_molecules << "; Concentration = " << (double)N_molecules/(V) << " nm^-3 ; " <<  (double)N_molecules/(V)/0.6 << "M" << endl;
    cout << "\t\tTime step = " << time_step << " s" << endl;
    cout << "\t\tScreening length = " << 1.0/(4*M_PI*bjerrum*(double)N_molecules/(V)) << " nm" << endl;
    cout << "\t\t kR = " << R*(4*M_PI*bjerrum*(double)N_molecules/V) << endl;
    cout << "\t LJ Interaction parameters: " << endl;
    cout << " \t\t anion-surface eps = " << lj_eps_constraint_anions << "; sigma = " << lj_sigma << "; rcut = " << lj_cut_constraint_anions << "; roff = " << r_off_anions << "; shift = " << lj_shift_anions << endl;
    cout << " \t\t cation-surface eps = " << lj_eps << "; sigma = " << lj_sigma << "; rcut = " << lj_cut << "; roff = " << 0 << "; shift = " << lj_shift_cations << endl;
    cout << " \t\t cation-anion eps = " << lj_eps << "; sigma = " << lj_sigma << "; rcut = " << lj_cut << "; roff = " << 0 << "; shift = " << lj_shift << endl;


    //create array of molecules
    molecule M[N_molecules];

    //setting anions
    cout << "\tSetting ionic properties..." << endl;
    for (int i=0;i<N_molecules;i++){
        M[i].type=2; M[i].q=1.0;
        if (i%2==0){
            M[i].q=-1.0;
            M[i].type=1;
            M[i].lj_cut_constraint=lj_cut_constraint_anions;
            M[i].lj_eps_constraint=lj_eps_constraint_anions;
        }
        cout << i << endl;
    M[i].print();
     }

     polymer poly[N_polymers];

//     system initialization;
     sim_system sys;
     sys.M = M;
     sys.poly = poly;

    cout << "\tSetting ionic properties..." << endl;
    for (int i=0;i<N_polymers;i++){
        cout << i << endl;
       poly[i].polymer_RW(2.0);
        cout << "from system print" << endl;
        sys.poly[i].print();
     }

     sys.create_particle_list();
     molecule ** part_list = sys.mol_list;
    //densities initialised to zero


     sys.gnuplot(2);

    int start_s=clock();

    return_molecules(M);

    cout << "\tWarming steps:" << endl;
    double acceptance = sys.mc_steps_mol(10000);
    acceptance = sys.mc_steps_pol(10000);

    cout << "\Main simulation steps:" << endl;
    for (int i=0;i<100;i++){
       	double acceptance = sys.mc_steps_mol(10000);
    	cout << "\t\t Step #" << i << "; Energy = " << sys.calc_total_energy() << "; Acceptance = " << acceptance << endl;
        acceptance = sys.mc_steps_pol(10000);
    	cout << "\t\t Step #" << i << "; Energy = " << sys.calc_total_energy() << "; Acceptance = " << acceptance << endl;
        //return_molecules(M);
        sys.gnuplot(i);


    }

    printf("\n Finish!\n");

    int stop_s=clock();
    cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " sec"<< endl;

    return 0;

}


