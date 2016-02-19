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
    delta = 0.1;
    k = 200.;
    pol_id++;
    id = pol_id;
    cout << "\ninit polymer ..." << endl;

//  set polymer coordinates to zero
    xc=yc=zc=0;
    for (int i = 0; i < N; i++){
//        M[i] = molecule();
        M[i].q = 0.;
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
    delta = 0.1;
    k = 200.;
    pol_id++;
    id = pol_id;
    cout << "\ninit polymer ..." << endl;

//  set polymer coordinates to zero
    xc=yc=zc=0;
    for (int i = 0; i < N; i++){
//        M[i] = molecule();
        M[i].q=0.0;
        M[i].type=1;
        //center of mass initialization
        xc+=M[i].x; yc+=M[i].y ; zc+=M[i].z;
        }

    xc/=(1.0*N);yc/=(1.0*N);zc/=(1.0*N);

}


double polymer::mindist( )
{

    double min_dist=lj_sigma*2.;
    for (int i = 0; i < N-1; ++i)
    {
        for (int j = i+1; j < N; ++j)
        {
            double r_ij = Distance(M[i],M[j]);
            if (r_ij<min_dist){min_dist = r_ij;}
        }
    }
    //cout << "Minimal distance between monomers = " <<min_dist << endl;
    return min_dist;
}

// function to create polymer by random walk... to be done. DONE!
void polymer::polymer_RW( /* int Nset  = 10 , */ double r /* = 0.0 */)
{

    // delta is the spring length
    delta = 0.1;
    k = 200.;


    double theta=acos(2*(drand48()-0.5));
    double phi=2.0*M_PI*drand48();
    double x = r*sin(theta)*cos(phi);
    double y = r*sin(theta)*sin(phi);
    double z = r*cos(theta);
    if (geometry == 0){
    if (x>L){x-=L;} else if (x<0.0){x+=L;}
    if (y>L){y-=L;} else if (y<0.0){y+=L;}
    if (z>L){z-=L;} else if (z<0.0){z+=L;}}
    M[0].move_to_position(x,y,z);


    for (int i = 1; i < N; ++i)
    {

        M[i].move_to_position(M[i-1]);
        M[i].advance(1,delta);
    }

}

void polymer::polymer_RW_WI( double x, double y, double z)
{

    // delta is the spring length
    delta = 0.1;
    k = 200.;

    M[0].move_to_position(x,y,z);
    double min_dist = 0.0;

    while (min_dist<0.045)
    {
        for (int i = 1; i < N; ++i)
        {

            M[i].move_to_position(M[i-1]);
            M[i].advance(1,delta);
        }

        polymer::update_COM();

        polymer::displace_polymer(x-xc,y-yc,z-zc);

        min_dist = polymer::mindist( );
    }
}

void polymer::polymer_SAW( double x, double y, double z)
{

    // delta is the spring length
    delta = 0.1;
    k = 200.;

    M[0].move_to_position(x,y,z);
    double min_dist = 0.0;
    M[1].move_to_position(M[0]);
    M[1].advance(1,delta*1.12);

    for (int i = 2; i < N; ++i)
    {
        min_dist = 0.1;
        while (min_dist<0.13)
        {
            min_dist = 0.14;
            M[i].move_to_position(M[i-1]);
            M[i].advance(1,delta*1.12);
            //check self-avoidance with previous molecules
            for (int j=0; j<i-1; j++)
            {
                double r_ij = Distance(M[i],M[j]);
                if (r_ij<min_dist){min_dist = r_ij;}
            }
        }
    }

        polymer::update_COM();

        polymer::displace_polymer(x-xc,y-yc,z-zc);

        min_dist = polymer::mindist( );
        cout << "Mindist = "<<min_dist << endl;


}

void polymer::displace_polymer( double x, double y, double z)
{
    for (int i = 0; i < N; ++i)
    {

        M[i].x+=x;
        M[i].y+=y;
        M[i].z+=z;
            if (geometry == 0){
    if (x>L){x-=L;} else if (x<0.0){x+=L;}
    if (y>L){y-=L;} else if (y<0.0){y+=L;}
    if (z>L){z-=L;} else if (z<0.0){z+=L;}}
    }
    xc+=x; yc+=y ; zc+=z;
}


//void polymer::inertia( double x, double y, double z)
//{
//    for (int i = 0; i < N; ++i)
//    {
//
//        M[i].x+=x;
//        M[i].y+=y;
//        M[i].z+=z;
//    }
//    xc+=x; yc+=y ; zc+=z;
//}

void polymer::polymer_RW( double x, double y, double z)
{

    // delta is the spring length
    delta = 0.1;
    k = 200.;

    M[0].move_to_position(x,y,z);


    for (int i = 1; i < N; ++i)
    {

        M[i].move_to_position(M[i-1]);
        M[i].advance(1,delta);
    }

    polymer::update_COM();

    polymer::displace_polymer(x-xc,y-yc,z-zc);


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
        //M[i].print();
    }
    cout << " Polymer printed ..." << endl;
}

polymer::~polymer()
{
    delete [] M;
    //dtor
}
