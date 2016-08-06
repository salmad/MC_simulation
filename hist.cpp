# include <vector>
# include <cstdlib>
# include <iomanip>
# include <ctime>
# include <cmath>
# include <string>
# include <sys/stat.h>
# include <iostream>
# include <fstream>
# include <stdio.h>
# include <stdlib.h>
# include <mpi.h>
#define NUM_COMMANDS 2

#include "molecule.h"
#include "global.h"
#include "include/integrate.h"
#include "include/polymer.h"
#include "include/sim_system.h"

using namespace std;

//time_step
extern string foldername="stars_";
const double time_step = 0.25;
const int nbins=80;
//Droplet size
const int geometry=0;  // 0 if cubic , 1 if spherical
const double L = 20.; // Nanometers
const double R = 8.0;
const double V=4./3.*M_PI*pow(R,3); // for sphere 4./3.*M_PI*pow(L-2.0,3)
//const double V=pow(L,3);
//Electrostatic parameter
const double bjerrum = 0.7;
//Lennard-jones parameters
const double lj_eps=0.05;
const double lj_sigma=0.1;            // Sodium atomic radius 2.5 A
const double lj_shift=0.0;
//const double lj_cut=1.12246204831*lj_sigma;
const double lj_cut=0.3;
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
double k=200.0;
double disntance_between_COM = 3.50;
double n0=0.0;
double n_end=2.5;


int N_arms1 = 2;
int N_arms2 = 2;
int N_monomers1[  ] = {17,17};
int N_monomers2[  ] = {17,17};
int N_monomers = 10;
double A_gauss = 0.65*pow(1.0*N_arms1*N_arms2,0.7);





//# Number of polymers and simple (single) molecules
int N_molecules = 2;
int N_polymers  = 0;
int N_stars     = 2;
int mol_id      = 0;
int pol_id      = 0;
int star_id     = 0;
int main(int argc, char* argv[]) {

    int np, pid ;
    double eTime,sTime,pTime;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    cout << "Starting simulation with " << np << " processors" << "; process id = "<< pid <<endl;
// initial and final number of molecues, number of steps
	double n_steps=50.0;
	if (argc>1){n0=atof(argv[1]);cout<< " Start with N0 = "<<n0<<endl;}
	if (argc>2){n_end=atof(argv[2]);cout<< " End with N_end = "<<n_end<<endl;}
	if (argc>3){n_steps=atof(argv[3]);cout<< " with  "<<n_steps<< " steps"<<endl;}

	cout << "argc = " << argc << endl;
	for(int i = 0; i < argc; i++)
		cout << "argv[" << i << "] = " << argv[i] << endl;


	srand (time(0)*(pid+1));
    srand48(time(0)*(pid+1));
    cout << "random seed = " << time(0)*(pid+1) << "pid = " << pid <<endl;

//  ofstream output;
//	output.open( "output.dat" , ofstream::out | ofstream::trunc);
    char str0[80];
    sprintf(str0,"Narm%d_%d_N%d_%d_%f_End_%f",N_arms1,N_arms2,N_monomers1[0],N_monomers2[0],n0,n_end);
    foldername=foldername+str0;
    mkdir((foldername).c_str(),0777);
	// freopen((foldername+"/output.dat").c_str(),"w",stdout);
//  random number initialisator
	MPI_Barrier(MPI_COMM_WORLD);

//    srand (1);
//    srand48(1);
    if (pid==0)
    {
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
    }

    //create array of molecules
    molecule M[N_molecules];

    //setting anions
    if (pid==0){cout << "\tSetting ionic properties..." << endl;}

    for (int i=0;i<N_molecules;i++){
        M[i].type=0; M[i].q=0;
        if (i%2==0){
            M[i].q=0;
            M[i].type=0;
            M[i].lj_cut_constraint=lj_cut_constraint_anions;
            M[i].lj_eps_constraint=lj_eps_constraint_anions;
        }
        cout << i << endl;
    M[i].print();
     }

     polymer poly[N_polymers];


     star    stars[ ]={star(N_arms1,L/2.,L/2.,L/2.),star(N_arms2,L/2.,L/2.,L/2.+disntance_between_COM)};
     stars[0].set_arms(N_arms1,N_monomers1,L/2.,L/2.,L/2.);
     cout << "ERROR" << endl;
     stars[1].set_arms(N_arms2,N_monomers2,L/2.,L/2.,L/2.+disntance_between_COM);

//     system initialization;
     sim_system sys;
     sys.M = M;
     sys.poly = NULL;
     sys.stars = stars;

    cout << "\tSetting ionic properties..." << endl;
    for (int i=0;i<N_polymers;i++){
        cout << i << endl;
      // poly[i].N = 3;
       poly[i].polymer_SAW(0.1,0.1,0.1);
       poly[i].k = 0;
       for (int j=0;j<poly[i].N;j++){
            poly[i].M[j].type = 0;
       }
        cout << "from system print" << endl;
        sys.poly[i].print();
     }


     for (int i=0;i<N_stars;i++){
        cout << i << endl;
        cout << "from system print" << endl;
        sys.stars[i].print();
     }





    sys.create_particle_list();
    molecule ** part_list = sys.mol_list;

    //densities initialised to zero


    int start_s=clock();


    return_molecules(M);
    cout << "\tWarming steps:" << endl;

    sys.create1D_linked_list();

    cout << "\tMain simulation steps:" << endl;

    MPI_Barrier(MPI_COMM_WORLD);


////////////////////////////////////////////////////////////
   int nsteps = 10000; int ntimes=40000;


   // Equilibrate configuration
	double acceptance = sys.mc_steps_star(100000);
	cout << "total energy " <<sys.calc_total_energy();
	while (sys.calc_total_energy()>1000.)
		{sys.mc_steps_star(10000);
			cout << "\ntotal energy " <<sys.calc_total_energy();}
	int hist_p[nbins]={};
	int hist[nbins] = {};
    double Rg = 0.0; double Rg_err = 0.;
    double Rg_p = 0.0; double Rg_err_p = 0.;
    double pos1 = 0.,pos2 = 0., pos_av=0.0;double pols_energy=0.;
    double r1 = 0.,r2 = 0.;


   for (int j=0;j<ntimes;j++){

           // Equilibrate configuration
       acceptance = sys.mc_steps_star(nsteps);
       stars[0].update_COM();
       stars[1].update_COM();
       // Data after equilibration
       pos2     = sqrt((stars[1].xc-stars[0].xc)*(stars[1].xc-stars[0].xc)+(stars[1].yc-stars[0].yc)*(stars[1].yc-stars[0].yc)+(stars[1].zc-stars[0].zc)*(stars[1].zc-stars[0].zc));
       hist_p[int(pos2/abs(n_end-n0)*nbins)]++;
       pols_energy += sys.calc_total_energy();
       Rg_p         += 0.5*(stars[1].gyration_radius()+stars[0].gyration_radius());
       Rg_err_p       = 0.5*(stars[1].gyration_radius()-stars[0].gyration_radius());




		if (j%(ntimes/100)==0){

			for(int ps = 0; ps < np; ps++) {
		    	MPI_Barrier(MPI_COMM_WORLD);
		    	if (ps == pid) {
					int stop_s=clock();
					cout << "pid " << pid<< " time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " sec"<< endl;
					cout<<"\t\t completion = " << j*100/ntimes << "% ; " << "; dist = " <<  pos2 << " ;Energy = " << pols_energy/(j+1) << " ;<z> = "<< stars[1].zc<<endl;
   					fflush(stdout);
   					plot_radius(hist_p,j+pid);
   					// sys.gnuplot(j+pid);
		    	}
			}

            MPI_Allreduce(&Rg_p, &Rg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&Rg_err_p, &Rg_err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            cout << "\t\t " << "; Rg = " << Rg/(j+1)/np << "; error +-   " << Rg_err/np  << endl;


			MPI_Allreduce(&hist_p, &hist, nbins, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			plot_radius(hist,j+np);
			// sys.gnuplot(j+np);
		}

    }


    cout << "\t\t Finish #" << "; Etot = " << pols_energy/ntimes << "; Accept. = " << acceptance  << endl;
   	plot_radius(hist,ntimes+1);

    return_molecules(M);
    sys.gnuplot(ntimes+1);

    int stop_s=clock();
    cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " sec"<< endl;
//


    MPI_Finalize();
    return 0;

}

