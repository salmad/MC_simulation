# include "../molecule.h"
# include "../global.h"
# include <vector>
# include <cstdlib>
//# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>
# include <iostream>
# include <fstream>


using namespace std;


molecule::molecule(){
      q=0.0;lj_cut_constraint=lj_cut;lj_eps_constraint=lj_eps;type=0;

      mol_id++;
      id = mol_id;
      // type 0 switches off all interactions with the particle.
      // it is not propagated

      if(geometry==1){
//	  id++;:
	  radius=(L)*pow(drand48(),1./3.);
	  double theta=acos(2*(drand48()-0.5));
	  double phi=2.0*M_PI*drand48();
	  x=radius*sin(theta)*cos(phi);
	  y=radius*sin(theta)*sin(phi);
	  z=radius*cos(theta);
      }
      else if (geometry == 0){
    	  x=drand48()*(L);
    	  y=drand48()*(L);
    	  z=drand48()*(L);
    	  radius=sqrt(x*x+y*y+z*z);
      }
      else {
    	  cout << "unknown geometry factor! check molecule setup!" << endl;
      }
  }

double molecule::advance(int n){
    // We are using a 1.0 fs timestep, this is converted
    double dt = time_step;
    radius =  sqrt(x*x+y*y+z*z);
    double r=0.0,dr;
     for(int i=0;i<n;i++){
    double dx=dt*(drand48()-0.5);
    double dy=dt*(drand48()-0.5);
    double dz=dt*(drand48()-0.5);
    x += dx;
    y += dy;
    z += dz;
//    periodic boundary condtions
//    if (x>L){x-=L;} else if (x<0.0){x+=L;}
//    if (y>L){y-=L;} else if (y<0.0){y+=L;}
//    if (z>L){z-=L;} else if (z<0.0){z+=L;}
    dr=sqrt(dx*dx+dy*dy+dz*dz);
    r+=dr;
    }
     radius = sqrt(x*x+y*y+z*z);
//cout << "\n distance moved = " << r;
    return r;
}

double molecule::advance(int n, double step){
    // We are using a 1.0 fs timestep, this is converted

    radius =  sqrt(x*x+y*y+z*z);
    double r=0.0,dr;
     for(int i=0;i<n;i++){

// unit vector*step move
    double theta=acos(2*(drand48()-0.5));
    double phi=2.0*M_PI*drand48();
    double dx = step*sin(theta)*cos(phi);
    double dy = step*sin(theta)*sin(phi);
    double dz = step*cos(theta);

    x += dx;
    y += dy;
    z += dz;
//    periodic boundary condtions
//    if (x>L){x-=L;} else if (x<0.0){x+=L;}
//    if (y>L){y-=L;} else if (y<0.0){y+=L;}
//    if (z>L){z-=L;} else if (z<0.0){z+=L;}
    dr=sqrt(dx*dx+dy*dy+dz*dz);
    r+=dr;
    }
     radius = sqrt(x*x+y*y+z*z);
//cout << "\n distance moved = " << r;
    return r;
}

double molecule::move_to_position(double x1, double y1, double z1){
    x = x1;
    y = y1;
    z = z1;
    radius = sqrt(x*x+y*y+z*z);
    return radius;
}

double molecule::move_to_position(molecule m){
    x = m.x;
    y = m.y;
    z = m.z;
    radius = sqrt(x*x+y*y+z*z);
    return radius;
}

int molecule::print(){
  cout <<"\n ID = " << id << ";  Molecule type = " << type <<" ; q = "<< q << "; Radius = " << radius << endl;
  cout <<" x = "<< x << "; y = " << y << "; z =  " << z << endl;
  cout <<" LJ parameters: sigma = " << lj_cut_constraint << "; eps = "<< lj_eps_constraint<<endl;
  return 0;
  }
