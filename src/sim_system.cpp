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
# include <cstring>
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


    //cells initialization

    dd.ncells = 1;
//    dd.hoc   ;
    dd.cell_neighbors=NULL;
    dd.linked_list = new int  [MAX_PART];
    dd.cell_size = 1.0;
    memset(dd.linked_list, -1, sizeof(dd.linked_list));

    //bond_list;
    //verlet_lsit;

}

// linked list should be created after particle lists (which refer to pointers)
void sim_system::create1D_linked_list()
{
    double z_range = 2*(R+2);
    dd.ncells = int(z_range/lj_cut);
    dd.cell_size = z_range/int(z_range/lj_cut);

    //set head of chains to zero

    memset( dd.hoc, -1, 600*sizeof(int) );
//    for (int i=0; i<dd.ncells; i++)
//    {
//        dd.hoc[i]=-1;// -1 means nothing is there
//    }

    // fill arrays with particles
    for (int i=0; i<max_part_id+1; i++)
    {
        int icell = int( ( (*mol_list[i]).z+z_range/2.)/dd.cell_size);

    //cout << "icell "<<icell << " i "<< i << " z = "<< (*mol_list[i]).z << " id" <<(*mol_list[i]).id  << endl;
        if (icell<dd.ncells )
        {
            dd.linked_list[i] = dd.hoc[icell];
            dd.hoc[icell] = (*mol_list[i]).id; // id and i should coincide
        }

    }



    //create neighbor list : to be done
//    dd.cell_neighbors = new int * ;
//    int neighbors[dd.ncells][3];
//    memset( neighbors, -1, dd.ncells*3*sizeof(int) );
//    cout << "neighbor of cell " << 0 << " is cell <<" << neighbors <<endl;
//    for (int i = 0; i<dd.ncells ; i++)
//    {
//        if (i==0) {neighbors[i][0] = i;neighbors[i][1] = i+1; }
//        else if (i==dd.ncells-1) {neighbors[i][0] = i-1;neighbors[i][1] = i;}
//        else {neighbors[i][0] = i-1;neighbors[i][1] = i;neighbors[i][2] = i+1;}
//        cout << "neighbor of cell " << i << " is cell " << neighbors[i][0] <<endl;
//    }
//    dd.cell_neighbors = neighbors;
}

//energy using cell list
double sim_system::recalc_energy_pol_ll(int pol_i, int mol_i)
{
    double e        = 0.0;
    molecule m1     = poly[pol_i].M[mol_i];
    int id1         = poly[pol_i].M[mol_i].id;
    double z_range  = dd.ncells*dd.cell_size;
    int icell       = int((m1.z+z_range/2.)/dd.cell_size);


    // calculate non-bonded part
    create1D_linked_list();
    //run through neigbouring cells
    if (icell>0)
    {
        int jcell = icell-1; //neighboring cell
        int id2 = dd.hoc[jcell]; //take the hoc of cell

        while (id2 != -1)
        {
            molecule m2 = *mol_list[id2];
            e += lj_energy(m1,m2);
//cout << "1 icell "<<icell << " jcell "<< jcell << " e = "<< e << " id " << id1 << " "<< id2  << " D "<< Distance(m1,m2)<< endl;
            id2 = dd.linked_list[id2];
        }

    }

    if (icell<dd.ncells-1 )
    {
        int jcell = icell+1; //neighboring cell
        int id2 = dd.hoc[jcell]; //take the hoc of cell
        while (id2 != -1)
        {
            molecule m2 = *mol_list[id2];
            e += lj_energy(m1,*mol_list[id2]);
//cout << "2 icell "<<icell << " jcell "<< jcell << " e = "<< e << " id " << id1 << " "<< id2  << " D "<< Distance(m1,m2)<< endl;

            id2 = dd.linked_list[id2];

        }
    }

    int jcell = icell; //neighboring cell
    int id2 = dd.hoc[jcell]; //take the hoc of cell

    while (id2 != -1)
    {
        if(id1!=id2)
        {
            e += lj_energy(m1,*mol_list[id2]);
//cout << "3 icell "<<icell << " jcell "<< jcell << " e = "<< e << " id " << id1 << " "<< id2  << " D "<< Distance(m1,*mol_list[id2])<< endl;

        }
        id2 = dd.linked_list[id2];

    }

    // go through bonded ... to be done

    e += poly[pol_i].bond_energy(mol_i);
    // we can add some extra energies if necessary, e.g.
    e += sim_system::extra_energy_pol();

    return e;
}



// make list of pointer to refer to main particle pointer M and poly

int sim_system::create_particle_list()
{
    max_part_id = -1;
    // for non-bonded interaction energy
    cout << " ## Creating molecules list ... " << endl;

    for (int i = 0; i < N_molecules; ++i)
    {
        max_part_id ++;
        //mol_list[i] = new molecule;
        mol_list[i] = &M[i];

    }

    cout << "Molecule list created. max molecule id= " << max_part_id << " ; Nmol = " << N_molecules << endl;
    cout << " ## Creating polymers list ... " << endl;
    for (int i = 0; i < N_polymers; i++)
    {
        cout << " \t polymer number  #" << i << " with " << poly[i].N << " molecules"<< endl;
        for (int j = 0; j < poly[i].N; ++j)
        {
            max_part_id ++;
            cout << "mol_list = " << max_part_id << endl;
            //mol_list[max_part_id] = new molecule;
//            molecule m_tmp = poly[i].M[j];
            mol_list[max_part_id] = &(poly[i].M[j]);

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
    for (int j=0; j<max_part_id+1; j++){
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
//    poly[0].update_COM();
    e += 50.0*((poly[0].xc)*(poly[0].xc)+(poly[0].yc)*(poly[0].yc)+(poly[0].zc)*(poly[0].zc));
//    poly[1].update_COM();
   e += 50.0*((poly[1].xc)*(poly[1].xc)+(poly[1].yc)*(poly[1].yc)+(poly[1].zc-disntance_between_COM)*(poly[1].zc-disntance_between_COM ));
//double rc = sqrt((poly[1].xc)*(poly[1].xc)+(poly[1].yc)*(poly[1].yc)+(poly[1].zc)*(poly[1].zc));
//    e += -2.*exp(-rc*rc/0.7/0.7/2.);
//    e += exp(500.*(rc-n_end));
//    e += exp(500.*(n0-rc));
//    e += log(4*M_PI*rc*rc);
//    poly[2].update_COM();
//    e += 100.0*(poly[2].xc-7.)*(poly[2].xc-7.)+(poly[2].yc)*(poly[2].yc)+(poly[2].zc)*(poly[2].zc);
//    poly[3].update_COM();
//    e += 100.0*(poly[3].xc+7.)*(poly[3].xc+7.)+(poly[3].yc)*(poly[3].yc)+(poly[3].zc)*(poly[3].zc);
//    cout << "error here? " << e << endl;
    return e;
}



// energy for a polymer bead move
double sim_system::recalc_energy_pol(int pol_i, int mol_i)
{

    double e=0.0;
    molecule m1 = (poly[pol_i].M[mol_i]);

 //   e+=constraint_energy(m1);
//    cout<< "recalc_energy: constraint energy = "<<e<< endl;

    // calculate non-bonded part
    for (int j = 0; j < max_part_id+1; j++){
        molecule m2 = *mol_list[j];
        if(m1.id != m2.id){

        //coulomb interaction
        //e += electrostatic_energy(m1,m2);
        //Lennard-jones
        e += lj_energy(m1,m2);
        }
    }

    // go through bonded ... to be done

    e += poly[pol_i].bond_energy(mol_i);

    // we can add some extra energies if necessary, e.g.
    e += sim_system::extra_energy_pol();

    return e;
}


double sim_system::calc_total_energy()
{
    double e = 0.0;
    // nonbonded part
    for(int i=0;i<max_part_id;i++){
        molecule m1=*mol_list[i];
        //e+=constraint_energy(m1);

        for (int j = i+1; j < max_part_id+1; j++){
            molecule m2 = *mol_list[j];
            //coulomb interaction
            //e += electrostatic_energy(m1,m2);
            //Lennard-jones
            e += lj_energy(m1,m2);
        }
    }
//    cout << "lj = " << e<<endl;

    //e += constraint_energy(*mol_list[max_part_id-1]);
    // bobnded part of energy
    for (int i = 0; i < N_polymers; i++)
    {
        e += poly[i].bond_energy();
//        cout << "bonded = " << poly[i].bond_energy()<<endl;
        poly[i].update_COM();
    }

    e += sim_system::extra_energy_pol();

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




int sim_system::move_COM_pol( int pol_ind){
//    double eold=calc_energy(M,size);
    double eold = calc_total_energy();

    double theta=acos(2*(drand48()-0.5));
    double phi=2.0*M_PI*drand48();
    double step = 0.5*drand48();
    double dx   = step*sin(theta)*cos(phi);
    double dy   = step*sin(theta)*sin(phi);
    double dz   = step*cos(theta);

    poly[pol_ind].displace_polymer(dx,dy,dz);

    double enew = calc_total_energy();


    double prob = exp((-enew+eold)/kT);
    if(drand48() < prob){
        return 1;}
    else {
        poly[pol_ind].displace_polymer(-dx,-dy,-dz);

//        cout << "probability " << drand48() << " and "<< exp(-enew+eold) << "; E new = "<< enew << " ; E old = "<< eold << endl;
        return 0;    }

}

int sim_system::move_pivot_pol( int pol_ind){
//    double eold=calc_energy(M,size);
    int mol_ind     = 1+ (int)(drand48()*(double)(poly[pol_ind].N-1));
    double eold     = calc_total_energy();

    double eold_s   = 0.;
    double enew_s   = 0.;

    double theta=0*acos(2*(drand48()-0.5));
    double phi=2.0*M_PI*drand48();
    molecule pivot = poly[pol_ind].M[mol_ind];
    for (int i = mol_ind+1 ; i<poly[pol_ind].N; i++)
    {
//        eold_s          += recalc_energy_pol(pol_ind,i);
        molecule  turned = poly[pol_ind].M[i];
//        double rij       = Distance(pivot,turned);
//        double phi_p     = atan((turned.y-pivot.y)/(turned.x-pivot.x));
//        double theta_p   = acos((turned.z-pivot.z)/rij);
        double x         = pivot.x+(turned.x-pivot.x)*cos(phi)-(turned.y-pivot.y)*sin(phi);
        double y         = pivot.y+(turned.x-pivot.x)*sin(phi)+(turned.y-pivot.y)*cos(phi);
        double z         = turned.z;
        poly[pol_ind].M[i].move_to_position(x,y,z);


    }
//    eold_s      +=  sim_system::extra_energy_pol();
    //keep COM constant; be careful to put it back when the move is rejected
    double xc=poly[pol_ind].xc;
    double yc=poly[pol_ind].yc;
    double zc=poly[pol_ind].zc;


    poly[pol_ind].update_COM();
    double Dxc=xc-poly[pol_ind].xc;
    double Dyc=yc-poly[pol_ind].yc;
    double Dzc=zc-poly[pol_ind].zc;
    poly[pol_ind].displace_polymer(Dxc,Dyc,Dzc); // put the polymer into the old COM

//    for (int i = mol_ind+1 ; i<poly[pol_ind].N; i++)
//    {
//        enew_s          += recalc_energy_pol(pol_ind,i);
//    }
//    enew_s      +=  sim_system::extra_energy_pol();
    double enew = calc_total_energy();

//    cout << "single " << enew_s-eold_s << " total " << enew-eold<<endl;

    double prob = exp((-enew+eold)/kT);
    if(drand48() < prob){
    //cout << "probability " << drand48() << " and "<< exp(-enew+eold) << "; E new = "<< enew << " ; E old = "<< eold << endl;
        return 1;}
    else {

        poly[pol_ind].displace_polymer(-Dxc,-Dyc,-Dzc);

        for (int i = mol_ind+1 ; i<poly[pol_ind].N; i++)
        {
            molecule turned  = poly[pol_ind].M[i];
            double x         = pivot.x+(turned.x-pivot.x)*cos(-phi)-(turned.y-pivot.y)*sin(-phi);
            double y         = pivot.y+(turned.x-pivot.x)*sin(-phi)+(turned.y-pivot.y)*cos(-phi);
            double z         = turned.z;
             poly[pol_ind].M[i].move_to_position(x,y,z);
        }
        poly[pol_ind].update_COM();

        //cout << "probability " << drand48() << " and "<< exp(-enew+eold) << "; E new = "<< enew << " ; E old = "<< eold << endl;
        return 0;    }

}

//Switch linked cell list here
int sim_system::mc_step_pol( int pol_ind){
//    double eold=calc_energy(M,size);
    int mol_ind = (int)(drand48()*(double)(poly[pol_ind].N));
    double eold = recalc_energy_pol_ll(pol_ind,mol_ind); //change into recalc_energy_pol() to use without linked lists
//    cout << " ll E = "<<eold <<endl;
//    eold = recalc_energy_pol(pol_ind,mol_ind); //change into recalc_energy_pol() to use without linked lists
//    cout << " " << eold <<endl;
   // eold       +=  sim_system::extra_energy_pol();

    double x1   = poly[pol_ind].M[mol_ind].x;
    double y1   = poly[pol_ind].M[mol_ind].y;
    double z1   = poly[pol_ind].M[mol_ind].z;

    poly[pol_ind].M[mol_ind].advance(1);

    double Dxc  = (poly[pol_ind].M[mol_ind].x-x1)/poly[pol_ind].N;
    double Dyc  = (poly[pol_ind].M[mol_ind].y-y1)/poly[pol_ind].N;
    double Dzc  = (poly[pol_ind].M[mol_ind].z-z1)/poly[pol_ind].N;

    poly[pol_ind].xc += Dxc;
    poly[pol_ind].yc += Dyc;
    poly[pol_ind].zc += Dzc;

    double enew = recalc_energy_pol_ll(pol_ind, mol_ind);
   // enew       +=  sim_system::extra_energy_pol();


    double prob = exp((-enew+eold)/kT);
    if(drand48() < prob){
        //cout << enew-eold<<endl;
        return 1;
        }
    else {
        poly[pol_ind].M[mol_ind].move_to_position(x1,y1,z1);

        poly[pol_ind].xc -= Dxc;
        poly[pol_ind].yc -= Dyc;
        poly[pol_ind].zc -= Dzc;
//        cout << "probability " << drand48() << " and "<< exp(-enew+eold) << "; E new = "<< enew << " ; E old = "<< eold << endl;
        return 0;    }

}


double sim_system::mc_steps_pol( int Nsteps){
    int success1=0;
    int success2=0;
    int success3=0;

//  for (int i=0;i<size;i++){
//    cout <<  "\n\n mc_steps"  << endl;
//    M[i].print();
//  }

    for (int i = 0; i < Nsteps; i++)
        {
            int index = (int)(drand48()*(double)(N_polymers));
//        cout <<  "index = "<<index << endl;
            if (drand48()<0.87){success1 += mc_step_pol( index);}
            else if (drand48()>0.93) {success2 += move_pivot_pol( index);}
            else {success3+=move_COM_pol( index);}

        }

    cout << "acceptances of moves: pivot = " << (double)success2/Nsteps/0.07 <<  " ; local = " << (double)success1/Nsteps/0.87 << " ; COM = " << (double)success3/Nsteps/0.05 << endl;
//        if (success==1){
////            M[index].print();
//            }
           // cout << "index = " <<index;


    return (double)(success1+success2+success3)/Nsteps;

}



double sim_system::mc_steps_pol_ll( int Nsteps){
    int success = 0;
    for (int i = 0; i < Nsteps; i++)
    {
        create1D_linked_list();
        int index = (int)(drand48()*(double)(N_polymers));
        success += mc_step_pol( index);
     }

    return (double)(success)/Nsteps;

}





int sim_system::gnuplot(int j /*, char * file_format*/){
	// Write a file

	    ofstream positions_file1;
        char str1[80];
        sprintf(str1, "/pos%3.3d.dat",j);
//        printf(str);
	    positions_file1.open( (foldername+str1).c_str(), ofstream::out | ofstream::trunc );
	    for(int i=0; i < N_molecules;i++){
	    positions_file1 << M[i].x << "  " << M[i].y << "  " << M[i].z << "  " << (M[i].q+1)/2*N_polymers << endl;
	    }
	    positions_file1.close();


	    //ofstream positions_file2;
        char str2[80];
        sprintf(str2, "/polymer_pos%3.3d.pdb",j);
//        printf(str);
	    //positions_file2.open( (foldername+str2).c_str(), ofstream::out | ofstream::trunc );
	    FILE * positions_file2;
	    positions_file2 = fopen ((foldername+str2).c_str(),"w");
	    int cnt = 1;
	    fprintf(positions_file2,"REMARK polymer configuration in MC code\n");

//	    for(int i=0; i < 20;i++)
//	    {
//            for(int kk=0; kk < 30;kk++)
//                {
//            //double theta=acos(2*(drand48()-0.5));
//            double theta=acos(-1.0+i*2.0/20.);
//            //double phi=M_PI*drand48();
//            double phi=2.0*M_PI*kk*1.0/20.;
//            double x = R*sin(theta)*cos(phi);
//            double y = R*sin(theta)*sin(phi);
//            double z = R*cos(theta);
//                //positions_file2 << poly[i].M[j].x << "  " << poly[i].M[j].y << "  " << poly[i].M[j].z << "  " << poly[i].id  << endl;
//            fprintf(positions_file2,"ATOM %6d%4s  UNX F%4d    %8.3f%8.3f%8.3f  0.00  0.00      T%03d \n",cnt, "FE",cnt-1,x,y,z , 10 );
//            cnt++;
//
//	    }}

//	    	    for(int i=0; i < N_molecules;i++)
//	    {
//                //positions_file2 << poly[i].M[j].x << "  " << poly[i].M[j].y << "  " << poly[i].M[j].z << "  " << poly[i].id  << endl;
//            fprintf(positions_file2,"ATOM %6d%4s  UNX F%4d    %8.3f%8.3f%8.3f  0.00  0.00      T%03d \n",cnt, "FE",cnt-1,M[i].x,M[i].y,M[i].z , int(M[i].q+1.0) );
//            cnt++;
//
//	    }

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
//	fprintf(gnuplotPipe, "splot '%s' u 1:2:3:4 w p palette pointtype 7  ps 1.,\\\n", (foldername+str1).c_str());
	fprintf(gnuplotPipe, "splot '%s' u 7:8:9:6 w lp palette pt 7 ps 2.0, '' u 7:8:9 w p  pt 6 lc rgb 'gray' lt 2 ps 4.0  \n", (foldername+str2).c_str());
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
    delete[] dd.linked_list;
}
