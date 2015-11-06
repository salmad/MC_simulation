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

polymer::polymer()
{
    N = 20;
    M = new molecule [N];
//    molecule m2[N];
//    M = m1 ;//ctor
    delta = 1.0;
    k = 20.;
    pol_id++;
    id = pol_id;
    cout << "\ninit polymer ..." << endl;

//  set polymer coordinates to zero
    xc=yc=zc=0;
    for (int i = 0; i < N; i++){
//        M[i] = molecule();
        M[i].q=-1.0;
        M[i].type=1;
        //center of mass initialization
        xc+=M[i].x; yc+=M[i].y ; zc+=M[i].z;
        }

    xc/=(1.0*N);yc/=(1.0*N);zc/=(1.0*N);
    //polymer id? we need it
}

polymer::polymer(int Nset /*double q , other parameters can be included to set up polymer*/)
{
    N = Nset;
    M = new molecule [N];
//    molecule m2[N];
//    M = m1 ;//ctor
    delta = 1.0;
    k = 20.;
    pol_id++;
    id = pol_id;
    cout << "\ninit polymer ..." << endl;

//  set polymer coordinates to zero
    xc=yc=zc=0;
    for (int i = 0; i < N; i++){
//        M[i] = molecule();
        M[i].q=-1.0;
        M[i].type=1;
        //center of mass initialization
        xc+=M[i].x; yc+=M[i].y ; zc+=M[i].z;
        }

    xc/=(1.0*N);yc/=(1.0*N);zc/=(1.0*N);

}


// function to create polymer by random walk... to be done. DONE!
void polymer::polymer_RW( /* int Nset  = 10 , */ double r /* = 0.0 */)
{

    // delta is the spring length
    delta = 1.0;
    k = 20.;


    double theta=acos(2*(drand48()-0.5));
    double phi=2.0*M_PI*drand48();
    double x = r*sin(theta)*cos(phi);
    double y = r*sin(theta)*sin(phi);
    double z = r*cos(theta);
    M[0].move_to_position(x,y,z);


    for (int i = 1; i < N; ++i)
    {

        M[i].move_to_position(M[i-1]);
        M[i].advance(1,delta);
    }

}

void polymer::update_COM()
{
    xc = yc = zc = 0.;
    for (int i = 0; i < N; i++){
        xc+=M[i].x; yc+=M[i].y ; zc+=M[i].z;
        }

    xc/=(1.0*N);yc/=(1.0*N);zc/=(1.0*N);
}

//double polymer::advance(int n, double step)
//{
//
//}



// implement here polymer move functions...

double polymer::bond_energy()
{
    double en = 0.0;
    for (int i = 0; i< N-1;i++){
        molecule mi   = M[i];
        molecule next = M[i+1];
        double d      = Distance(mi,next);
        en           += k*(d-delta)*(d-delta)*0.5;
    }
    //dtor
    return en;
}

// We can accelerate energy calculations by calculating
// energy of two springs that are connected with molecule #j which belongs to 0 ... N-1

double polymer::bond_energy(int j)
{
    double en = 0.0;
    if (j==0)
    {
        molecule mj         = M[j];
        molecule next       = M[j+1];
        double d            = Distance(mj,next);
        en                 += k*(d-delta)*(d-delta)*0.5;
        return en;
    }
    else if (j==N-1)
    {
        molecule previous   = M[j-1];
        molecule mj         = M[j];
        double d            = Distance(mj,previous);
        en                 += k*(d-delta)*(d-delta)*0.5;
        return en;
    }
    else
    {
        molecule previous   = M[j-1];
        molecule mj         = M[j];
        molecule next       = M[j+1];
        double d1           = Distance(mj,next);
        double d2           = Distance(mj,previous);
        en                 += (k*(d1-delta)*(d1-delta)*0.5+k*(d2-delta)*(d2-delta)*0.5);
        return en;
    }

}

void polymer::print(){

    cout << "\n Print polymer ID " << id << " ; COM : x = " << xc << " y = " << yc << " z = " << zc << endl;
    cout << " N = " << N << " ;  "<< endl;
    for (int i = 0; i < N; ++i)
    {
        /* print all molecules */
        M[i].print();
    }
    cout << " Polymer printed ..." << endl;
}

polymer::~polymer()
{
    delete [] M;
    //dtor
}
