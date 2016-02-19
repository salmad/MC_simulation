#ifndef POLYMER_H
#define POLYMER_H
#include "../molecule.h"
# include "../global.h"
using namespace std;

class polymer
{
    public:
        int N;

        // polymer center of mass
        double xc,yc,zc;

        molecule * M;

        // most likely we need polymer id

        //bonds
        double k;   // spring constant
        double delta; //inter-monomer distance
        int id ;

        polymer();
        polymer(int Nset);
        void polymer_RW( /* int Nset = 10, */ double r = 0.0);
        void polymer_RW( double x ,double y ,double z );
        double mindist( );
        void polymer_RW_WI( double x ,double y ,double z );
        void polymer_SAW( double x, double y, double z);
        void displace_polymer( double x ,double y ,double z );
        virtual ~polymer();
        void update_COM();

        // function to create polymer by random walk... to be done..done
        // functions should be private in future
        // polymer moves implemented here
        double bond_energy(int j);
        double bond_energy();
        void print();

    protected:
    private:

};

#endif // POLYMER_H
