#include "../include/star.h"
#include "../include/polymer.h"
#include "../molecule.h"
#include "../include/integrate.h"
#include "../global.h"
# include <iostream>

# include <vector>
# include <cstdlib>
# include <iomanip>
# include <ctime>
# include <cmath>
# include <fstream>
using namespace std;



star::star()
{
    N_arms = 2;
    poly = new polymer [N_arms];


    star_id++;
    id = star_id;

    cout << "\ninit star ..." << endl;

    //  set star coordinates to zero
    xc=yc=zc=0;
    for (int i = 0; i < N_arms; i++){
        xc+=poly[i].xc; yc+=poly[i].yc ; zc+=poly[i].zc;
        }

    xc/=(1.0*N_arms);yc/=(1.0*N_arms);zc/=(1.0*N_arms);

}

star::star(int N, double x, double y, double z)
{
    N_arms = N;
    poly = new polymer [N_arms];

    star_id++;
    id = star_id;

    cout << "\ninit star ..." << endl;

    //  set star coordinates to zero
    xc=yc=zc=0;
    for (int i = 0; i < N_arms; i++){
    	poly[i].polymer_SAW(x,y,z);
        xc+=poly[i].xc; yc+=poly[i].yc ; zc+=poly[i].zc;
        }

    xc/=(1.0*N_arms);yc/=(1.0*N_arms);zc/=(1.0*N_arms);

}

void star::displace( double dx, double dy, double dz)
{
    for (int i = 0; i < N_arms; ++i)
    {
    	poly[i].displace(dx,dy,dz);
    }
    xc += dx; yc += dy ; zc += dz;
}


void star::update_COM()
{
    //  set star coordinates to zero
    xc=yc=zc=0;
    for (int i = 0; i < N_arms; i++){
    	poly[i].update_COM();   // for safety to ensure that COM of polymers is up to date
        xc+=poly[i].xc; yc+=poly[i].yc ; zc+=poly[i].zc;
        }

    xc/=(1.0*N_arms);yc/=(1.0*N_arms);zc/=(1.0*N_arms);
}


void star::set_COM( double x ,double y ,double z ){
	star::update_COM();
	star::displace(x-xc,y-yc,z-zc);
}


// connected to a central polymer
// double star::bond_energy()
// {
//     double en = 0.0;
//     for (int i = 0; i< N_arms;i++){
//     	en 		 		   += poly[i].bond_energy();
//     }

//     /*the center of the star is the first polymer*/
//     molecule mc 			= poly[0].M[(poly[0].N+1)/2] ; 
//     double k 				= poly[0].k;
//     double delta 			= poly[0].delta;


//     for (int i = 1; i< N_arms;i++){
//         molecule m_center   = poly[i].M[(poly[i].N+1)/2];  /* The central molecule to each polymer*/
//         double d     	    = Distance(mc,m_center);
//         en          	   += k*(d-delta)*(d-delta)*0.5;
//     }
//     //dtor
//     return en;
// }

//connected by to the end of the first polymer
double star::bond_energy()
{
    double en = 0.0;
    for (int i = 0; i< N_arms;i++){
    	en 		 		   += poly[i].bond_energy();
    }

    /*the center of the star is the first polymer*/
    molecule mc 			= poly[0].M[0] ; 
    double k 				= poly[0].k;
    double delta 			= poly[0].delta;


    for (int i = 1; i< N_arms;i++){
        molecule m_center   = poly[i].M[0];  /* The central molecule to each polymer*/
        double d     	    = Distance(mc,m_center);
        en          	   += k*(d-delta)*(d-delta)*0.5;
    }
    //dtor
    return en;
}


//connected to the center of the polymer
// double star::bond_energy(int i /*polymer index*/,int j /* monomer index*/)
// {
//     double en = 0.0;
//     en 		 += poly[i].bond_energy(j);

//         /*the center of the star is the first polymer*/
//     molecule mc 			= poly[0].M[(poly[0].N+1)/2] ; 
//     double k 				= poly[0].k;
//     double delta 			= poly[0].delta;


//     for (int i = 1; i< N_arms;i++){
//         molecule m_center   = poly[i].M[(poly[i].N+1)/2];  /* The central molecule to each polymer*/
//         double d     	    = Distance(mc,m_center);
//         en          	   += k*(d-delta)*(d-delta)*0.5;
//     }
//     //dtor
//     return en;

// }



//connected by ends
double star::bond_energy(int i /*polymer index*/,int j /* monomer index*/)
{
    double en = 0.0;
    en 		 += poly[i].bond_energy(j);

        /*the center of the star is the first polymer*/
    molecule mc 			= poly[0].M[0] ; 
    double k 				= poly[0].k;
    double delta 			= poly[0].delta;


    for (int i = 1; i< N_arms;i++){
        molecule m_center   = poly[i].M[0];  /* The central molecule to each polymer*/
        double d     	    = Distance(mc,m_center);
        en          	   += k*(d-delta)*(d-delta)*0.5;
    }
    //dtor
    return en;

}

void star::print(){

    cout << "\n Print star ID " << id << " ; COM : x = " << xc << " y = " << yc << " z = " << zc << endl;
    cout << " N arms = " << N_arms << " ;  "<< endl;
    for (int i = 0; i < N_arms; ++i)
    {
        /* print all molecules */
        poly[i].print();
    }
    cout << " Star printed ..." << endl;
}


star::~star()
{
    delete [] poly;
}
