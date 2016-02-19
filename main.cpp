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
const double time_step = 0.25;
const int nbins=150;
//Droplet size
const int geometry=1;  // 0 if cubic , 1 if spherical
const double L = 20.; // Nanometers
const double R = 19.0;
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
double k=2.0;
double disntance_between_COM = 3.0;
double n0=2.5;
double n_end=0.0;

//# Number of polymers and simple (single) molecules
int N_molecules = 2;
int N_polymers  = 2;
int mol_id      = 0;
int pol_id      = 0;
int main(int argc, char* argv[]) {

// initial and final number of molecues, number of steps
	double n_steps=30.0;
	if (argc>1){n0=atof(argv[1]);cout<< " Start with N0 = "<<n0<<endl;}
	if (argc>2){n_end=atof(argv[2]);cout<< " End with N_end = "<<n_end<<endl;}
	if (argc>3){n_steps=atof(argv[3]);cout<< " with  "<<n_steps<< " steps"<<endl;}

	cout << "argc = " << argc << endl;
	for(int i = 0; i < argc; i++)
		cout << "argv[" << i << "] = " << argv[i] << endl;

//  ofstream output;
//	output.open( "output.dat" , ofstream::out | ofstream::trunc);
    char str0[80];
    sprintf(str0,"N0_%f_End_%f",n0,n_end);
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

//     system initialization;
     sim_system sys;
     sys.M = M;
     sys.poly = poly;

    cout << "\tSetting ionic properties..." << endl;
    for (int i=0;i<N_polymers;i++){
        cout << i << endl;
       poly[i].polymer_RW(0.1,0.1,0.1);
        cout << "from system print" << endl;
        sys.poly[i].print();
     }

     sys.create_particle_list();
     molecule ** part_list = sys.mol_list;

    //densities initialised to zero


     //sys.gnuplot(2);

    int start_s=clock();


    return_molecules(M);
    cout << "\tWarming steps:" << endl;

    poly[0].update_COM();
    poly[1].update_COM();
        cout << "\tHello" << endl;
    sys.create1D_linked_list();


    //double acceptance = sys.mc_steps_mol(1000);

    double acceptance = sys.mc_steps_pol(10000);

    return_molecules(M);
    cout << "\Main simulation steps:" << endl;

    sys.create1D_linked_list();
//    poly[0].polymer_RW_WI(0.,0.,0.);
//    poly[1].polymer_RW_WI(0.,0.,0.+disntance_between_COM);
//     acceptance = sys.mc_steps_pol(100000);
//    for (int i =0 ; i<200 ; i++)
//    {
//        sys.move_pivot_pol(1);
//        sys.gnuplot(i);
//        cout<< "i = " << i << endl;
//    }


////////////////////////////////////////////////////////////
//    int nsteps = 10000; int ntimes=100000;
//    // Generate configurations :
//    poly[0].polymer_RW_WI(0.,0.,0.);
//    poly[1].polymer_RW_WI(0.,0.,0.+disntance_between_COM);
//    // Equilibrate configuration
//    acceptance = sys.mc_steps_pol(500000);
//    int hist[nbins]={};
//
//
//    double pos2 = 0.;double pols_energy=0.;
//
//
//    for (int j=0;j<ntimes;j++){
//
//            // Equilibrate configuration
//            acceptance = sys.mc_steps_pol(nsteps);
//            poly[0].update_COM();
//            poly[1].update_COM();
//            // Data after equilibration
//            pos2     = sqrt(poly[1].zc* poly[1].zc+ poly[1].yc* poly[1].yc+poly[1].xc* poly[1].xc);
//            hist[int(pos2/(R)*nbins)]++;
//            pols_energy += sys.calc_total_energy();
//
//            if (j%(ntimes/100)==0)
//            {
//                int stop_s=clock();
//                cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " sec"<< endl;
//                cout<<"\t\t completion = " << j*100/ntimes << "% ; " << "; dist = " <<  pos2 << " ;Energy = " << pols_energy/j << " ;<z> = "<< poly[1].zc<<endl;
//
//                plot_radius(hist,j);
//            }
//
//        }
//
//    	cout << "\t\t Finish #" << "; Etot = " << pols_energy/ntimes << "; Accept. = " << acceptance  << endl;
//    	plot_radius(hist,10001);
//
//        return_molecules(M);
//        sys.gnuplot(1);
//
//        int stop_s=clock();
//        cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " sec"<< endl;
//




///////////////////////////////////////////////////////////////////////////////////
//
    int nsteps = 10000; int ntimes=3000;
    cout << "nsteps " << nsteps << " ntimes= " << ntimes << endl;
    for (int i=0;i<n_steps;i++){

        disntance_between_COM = n0-i*(n0-n_end)/(n_steps-1);
        cout<<"i = " << i<<endl;

        double pos1 = 0.,pos2 = 0., pos_av=0.0;double pols_energy=0.;
        double r1 = 0.,r2 = 0.;

        // Generate configurations :
        poly[0].polymer_SAW(0.,0.,0.);
        poly[1].polymer_SAW(0.,0.,0.+disntance_between_COM);
        // Equilibrate configuration
        acceptance = sys.mc_steps_pol(3000);
        for (int j=0;j<ntimes;j++){

//            // Generate configurations :
//            poly[0].polymer_SAW(0.,0.,0.);
//            poly[1].polymer_SAW(0.,0.,0.+disntance_between_COM);
            // Initial radius
            r1 += (pow(Distance(poly[0].M[0], poly[0].M[poly[0].N-1]),2)+pow(Distance(poly[1].M[0], poly[1].M[poly[1].N-1]),2));

            cout << " ;energy = "<<sys.calc_total_energy()  << " ;acceptance= " << acceptance << endl;


            // Equilibrate configuration
            acceptance = sys.mc_steps_pol(nsteps);
            poly[0].mindist();
            poly[1].mindist();
            // Data after equilibration
            pos1    += poly[0].zc;
            pos2    += poly[1].zc;
            pos_av  += 0.5*(poly[1].zc-poly[0].zc-1.*disntance_between_COM);
            pols_energy += sys.calc_total_energy();
            r2 += (pow(Distance(poly[0].M[0], poly[0].M[poly[0].N-1]),2)+pow(Distance(poly[1].M[0], poly[1].M[poly[1].N-1]),2));

            if (j%(ntimes/100)==0)
            {
                int stop_s=clock();
                cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " sec"<< endl;
                cout<<"\t\t completion = " << j*100/ntimes << "% ; " << "; Rad = " << 0.5*r1/(j+1) << " ; " << 0.5*r2/(j+1)<< "; Etot = " << pols_energy/(j+1) <<  " ;z = "<<disntance_between_COM << " ;<z> = "<<-poly[0].zc<< " ; " << poly[1].zc-disntance_between_COM <<endl;

            }

        }

    	cout << "\t\t Step #" << i << "; Etot = " << pols_energy/ntimes << "; Accept. = " << acceptance << "; <dx> = " << pos_av/ntimes << " ;z0 = " << disntance_between_COM << " ;<z1> = "<< pos1/ntimes << " ; <z2> = " << pos2/ntimes << endl;
    	cout << " Radii= " << 0.5*r1/ntimes << " ; " << 0.5*r2/ntimes<< endl;

        return_molecules(M);
        sys.gnuplot(i);

        int stop_s=clock();
        cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " sec"<< endl;

    }

//
//
//    for (int i=0;i<n_steps;i++){
//        disntance_between_COM = n_end+(i)*(n0-n_end)/(n_steps-1);
//
//        double pos1 = 0.,pos2 = 0., pos_av=0.0;double pols_energy=0.;
//        double r1 = 0.,r2 = 0.;
//
//        for (int j=0;j<ntimes;j++){
//
//            // Generate configurations :
//            poly[0].polymer_RW_WI(0.,0.,0.);
//            poly[1].polymer_RW_WI(0.,0.,0.+disntance_between_COM);
//            // Initial radius
//            r1 += (pow(Distance(poly[0].M[0], poly[0].M[poly[0].N-1]),2)+pow(Distance(poly[0].M[0], poly[0].M[poly[0].N-1]),2));
//
//            cout << " ;energy = "<<sys.calc_total_energy()  << " ;acceptance= " << acceptance << endl;
//
//            // Equilibrate configuration
//            acceptance = sys.mc_steps_pol(nsteps);
//            // Data after equilibration
//            pos1    += poly[0].zc;
//            pos2    += poly[1].zc;
//            pos_av  += 0.5*(poly[1].zc-poly[0].zc-1.*disntance_between_COM);
//            pols_energy += sys.calc_total_energy();
//            r2 += (pow(Distance(poly[0].M[0], poly[0].M[poly[0].N-1]),2)+pow(Distance(poly[0].M[0], poly[0].M[poly[0].N-1]),2));
//
//            if (j%(ntimes/100)==0){cout<<"\t\t completion = " << j*100/ntimes << "% ; " << "; Rad = " << 0.5*r1/(j+1) << " ; " << 0.5*r2/(j+1)<< " ;Energy = " << sys.calc_total_energy()<<  " ;z = "<<disntance_between_COM << " ;<z> = "<<-poly[0].zc<< " ; " << poly[1].zc-disntance_between_COM <<endl;}
//
//        }
//    	cout << "\t\t Step #" << i << "; Etot = " << pols_energy/ntimes << "; Accept. = " << acceptance << "; <dx> = " << pos_av/ntimes << " ;z0 = " << disntance_between_COM << " ;<z1> = "<< pos1/ntimes << " ; <z2> = " << pos2/ntimes << endl;
//    	cout << " Radii= " << 0.5*r1/ntimes << " ; " << 0.5*r2/ntimes<< endl;
//        return_molecules(M);
//        sys.gnuplot(i);
//
//
//    }
//
//
//    printf("\n Finish!\n");
//
//    int stop_s=clock();
//    cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " sec"<< endl;
//
    return 0;

}

//        poly[0].polymer_RW(0.,0.,0.);
//        poly[1].polymer_RW(0.,0.,0.+disntance_between_COM);
//        int counts = 0;
//        disntance_between_COM = 0.1;
//        for (int j=0;j<1000;j++)
//        {
//
//            acceptance = sys.mc_steps_pol(20);
//            cout<< "check energy Etot = " << sys.calc_total_energy() << " ; Nsteps = "<< 20*j<<endl;
//            if(j>200){pols_energy += sys.calc_total_energy();counts++;}
//        }
//
//        cout << "\t\t Finish checking... " << "; <E> " << pols_energy/(1.*counts) << endl;
