# ifndef STAR_H
# define STAR_H
# include "../molecule.h"
# include "../global.h"
# include "polymer.h"
using namespace std;

class star
{
    public:

    	int N_arms;
    	polymer * poly;
    	// star center of mass
		double xc,yc,zc;

		int id ;

		star();
		star( int N_arms );

		void displace( double dx ,double dy ,double dz );
		void set_COM( double x ,double y ,double z );
		void update_COM();

		// bond energies derived from polymer bond energy
		double bond_energy(int i, int j /*jth atom from ith polymer*/);
		double bond_energy();
		void print(); /*prints properties of the star*/
        virtual ~star();



    protected:
    private:
};

#endif // STAR_H
