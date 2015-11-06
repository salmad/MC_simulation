# ifndef MOLECULE_H
# define MOLECULE_H
# include "global.h"
using namespace std;



// Class for keeping track of the properties for a particle/molecule
class molecule{
    public:
        double x;		// position in x axis
        double y;		// position in y axis
        double z;
        double q;
        double radius;
        double lj_cut_constraint;
        double lj_eps_constraint;
        int type ;
        int id ;
  //static int id;
        molecule();
//double radius(){
//    return  sqrt(x*x+y*y+z*z);
//}
        double move_to_position(molecule m);
        double advance(int );
        double move_to_position(double , double , double );
        double advance(int n, double step);

        int print();

};


#endif // MOLECULE_H
