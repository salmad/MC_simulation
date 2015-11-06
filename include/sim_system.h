#ifndef SIM_SYSTEM_H
#define SIM_SYSTEM_H
#include "../molecule.h"
#include "polymer.h"

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


        sim_system();
        int create_particle_list();
        double recalc_energy_mol(int );
        double recalc_energy_pol(int , int);
        double extra_energy_pol();
        double calc_total_energy();

        int mc_step_mol( int );
        double mc_steps_mol( int );
        int mc_step_pol( int );
        double mc_steps_pol( int );

        int gnuplot(int j);

        virtual ~sim_system();
    protected:
    private:

    //bond_list;
    //verlet_lsit;

};

#endif // SIM_SYSTEM_H
