#ifndef SIM_SYSTEM_H
#define SIM_SYSTEM_H
#include "../molecule.h"
#include "polymer.h"
#include <vector>

using namespace std;
#define MAX_PART 5000

class sim_system
{
    public:
        molecule * M;
        polymer  * poly;
        molecule *other;
        // how to change mol_list into list of references , e.g. (molecule *) mol_list[]
        molecule ** mol_list;

        int max_part_id;
        //capsule cap;
        double temperature;
        double box;
        double volume;

        //Linked cells, decomposition

        struct dom_decomposition
        {
            int ncells;
            int ** cell_neighbors;
            int  hoc[600];
            int * linked_list;
            double cell_size;
        } ;
        dom_decomposition dd;
        void create1D_linked_list();


        sim_system();

        int create_particle_list();         //creates list of all particles
        double recalc_energy_mol(int );     // calculates energy of a single molecule
        double recalc_energy_pol(int , int); // calculatels energy when a molecule
                                            //  within a polymer is specified

        double recalc_energy_pol_ll(int pol_i, int mol_i); //calculates energy of a molecule using linked lists
        int move_COM_pol( int pol_ind);                    // attempt of moving polymer as a whole
        double extra_energy_pol();                         // general extra energy terms that can be added
                                                            // to apply constraints to the system
        double calc_total_energy();                         //calculates total energy of the system
        int move_pivot_pol( int );                          //attempt of a povit move

        int mc_step_mol( int );
        double mc_steps_mol( int );

        int mc_step_pol( int );
        double mc_steps_pol( int );
        double mc_steps_pol_ll( int );

        int gnuplot(int j);

        virtual ~sim_system();
    protected:
    private:

    //bond_list;
    //verlet_lsit;

};

#endif // SIM_SYSTEM_H
