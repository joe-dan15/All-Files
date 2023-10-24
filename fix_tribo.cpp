/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_tribo.h"
#include "atom.h"
#include "pair_gran.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "modify.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "fix_multisphere.h"  
#include "error.h"
#include "force.h"
#include "neigh_list.h"
#include "properties.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{CHUTE,SPHERICAL,VECTOR,LUNAR};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixTribo::FixTribo(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix tribo command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;

  mstr = vstr = pstr = tstr = xstr = ystr = zstr = NULL;
  mstyle = vstyle = pstyle = tstyle = xstyle = ystyle = zstyle = CONSTANT;

  int iarg = 3;

  if (strcmp(arg[3],"lunar") == 0) { // check what environment is wanted
    if (strstr(arg[iarg+1], "v_") == arg[3]){ // check if the value that follows is a variable 
      int n = strlen(&arg[iarg+1][2]) + 1; // creates an integer the same as the length of the variable name
      costr = new char[n]; // creates a char array of length n
      strcpy(costr,&arg[iarg+1][2]); // copies values of variable name into char array
    }else{
      cut_off = force->numeric(FLERR, arg[iarg+1]); // if value is not a name, put value straight into the "cut_off" variable
      style = LUNAR; // set style to lunar
    }
  } else{
    error->all(FLERR, "Please use valid environment... Currently Lunar is only supported"); // if requirements not met then throw error
  }
  /*
  if (strcmp(arg[4],"chute") == 0) {
    if (narg != 6) error->all(FLERR,"Illegal fix gravity command");
    style = CHUTE;
    if (strstr(arg[5],"v_") == arg[5]) {
      int n = strlen(&arg[5][2]) + 1;
      vstr = new char[n];
      strcpy(vstr,&arg[5][2]);
      vstyle = EQUAL;
    } else {
      vert = force->numeric(FLERR,arg[5]);
      vstyle = CONSTANT;
    }

  } else if (strcmp(arg[4],"spherical") == 0) {
    if (narg != 7) error->all(FLERR,"Illegal fix gravity command");
    style = SPHERICAL;
    if (strstr(arg[5],"v_") == arg[5]) {
      int n = strlen(&arg[5][2]) + 1;
      pstr = new char[n];
      strcpy(pstr,&arg[5][2]);
      pstyle = EQUAL;
    } else {
      phi = force->numeric(FLERR,arg[5]);
      pstyle = CONSTANT;
    }
    if (strstr(arg[6],"v_") == arg[6]) {
      int n = strlen(&arg[6][2]) + 1;
      tstr = new char[n];
      strcpy(tstr,&arg[6][2]);
      tstyle = EQUAL;
    } else {
      theta = force->numeric(FLERR,arg[6]);
      tstyle = CONSTANT;
    }

  } else if (strcmp(arg[4],"vector") == 0) {
    if (narg != 8) error->all(FLERR,"Illegal fix gravity command");
    style = VECTOR;
    if (strstr(arg[5],"v_") == arg[5]) {
      int n = strlen(&arg[5][2]) + 1;
      xstr = new char[n];
      strcpy(xstr,&arg[5][2]);
      xstyle = EQUAL;
    } else {
      xdir = force->numeric(FLERR,arg[5]);
      xstyle = CONSTANT;
    }
    if (strstr(arg[6],"v_") == arg[6]) {
      int n = strlen(&arg[6][2]) + 1;
      ystr = new char[n];
      strcpy(ystr,&arg[6][2]);
      ystyle = EQUAL;
    } else {
      ydir = force->numeric(FLERR,arg[6]);
      ystyle = CONSTANT;
    }
    if (strstr(arg[7],"v_") == arg[7]) {
      int n = strlen(&arg[7][2]) + 1;
      zstr = new char[n];
      strcpy(zstr,&arg[7][2]);
      zstyle = EQUAL;
    } else {
      zdir = force->numeric(FLERR,arg[7]);
      zstyle = CONSTANT;
    }

  } else error->all(FLERR,"Illegal fix gravity command");*/

  degree2rad = MY_PI/180.0;
  time_origin = update->ntimestep;

  eflag = 0;
  egrav = 0.0;

  fm = NULL; 
}

/* ---------------------------------------------------------------------- */

FixTribo::~FixTribo()
{
  delete [] mstr;
  delete [] vstr;
  delete [] pstr;
  delete [] tstr;
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] costr;
}

/* ---------------------------------------------------------------------- */

int FixTribo::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTribo::init()
{
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));

  // check variables
/*
  if (mstr) {
    mvar = input->variable->find(mstr);
    if (mvar < 0)
      error->all(FLERR,"Variable name for fix tribo does not exist");
    if (!input->variable->equalstyle(mvar))
      error->all(FLERR,"Variable for fix tribo is invalid style");
  }
  if (vstr) {
    vvar = input->variable->find(vstr);
    if (vvar < 0)
      error->all(FLERR,"Variable name for fix tribo does not exist");
    if (!input->variable->equalstyle(vvar))
      error->all(FLERR,"Variable for fix tribo is invalid style");
  }
  if (pstr) {
    pvar = input->variable->find(pstr);
    if (pvar < 0)
      error->all(FLERR,"Variable name for fix tribo does not exist");
    if (!input->variable->equalstyle(pvar))
      error->all(FLERR,"Variable for fix tribo is invalid style");
  }
  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Variable name for fix tribo does not exist");
    if (!input->variable->equalstyle(tvar))
      error->all(FLERR,"Variable for fix tribo is invalid style");
  }*/
  if (costr) {
    covar = input->variable->find(costr);
    if (!input->variable->equalstyle(covar))
      error->all(FLERR,"Variable for fix tribo is invalid style");
  }/*
  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix tribo does not exist");
    if (!input->variable->equalstyle(xvar))
      error->all(FLERR,"Variable for fix tribo is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix tribo does not exist");
    if (!input->variable->equalstyle(yvar))
      error->all(FLERR,"Variable for fix tribo is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix tribo does not exist");
    if (!input->variable->equalstyle(zvar))
      error->all(FLERR,"Variable for fix tribo is invalid style");
  }*/

  varflag = CONSTANT;
  if (mstyle != CONSTANT || vstyle != CONSTANT || pstyle != CONSTANT ||
      tstyle != CONSTANT || xstyle != CONSTANT || ystyle != CONSTANT ||
      zstyle != CONSTANT) varflag = EQUAL;

  // set gravity components once and for all

  if (varflag == CONSTANT) set_acceleration();

  fm = NULL;
  int nms = modify->n_fixes_style("multisphere");
  if(nms > 1)
    error->fix_error(FLERR,this,"support for more than one fix multisphere not implemented");
  if(nms)
    fm = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));

  int n_relax = modify->n_fixes_style("relax");
  if(n_relax > 1)
        error->fix_error(FLERR,this,"does not work with more than 1 fix relax");
  else if(1 == n_relax)
        fix_relax = static_cast<FixRelaxContacts*>(modify->find_fix_style("relax",0));
  else
        fix_relax = 0;
}

/* ---------------------------------------------------------------------- */

void FixTribo::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixTribo::post_force(int vflag)
{
  // update gravity due to variables

  cut_off_dist = input->variable->compute_equal(covar);
  /*
  if (varflag != CONSTANT) {
    modify->clearstep_compute();
    if (mstyle == EQUAL) magnitude = input->variable->compute_equal(mvar);
    if (vstyle == EQUAL) vert = input->variable->compute_equal(vvar);
    if (pstyle == EQUAL) phi = input->variable->compute_equal(pvar);
    if (tstyle == EQUAL) theta = input->variable->compute_equal(tvar);
    if (xstyle == EQUAL) xdir = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) ydir = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) zdir = input->variable->compute_equal(zvar);
    modify->addstep_compute(update->ntimestep + 1);

    set_acceleration();
  }*/


  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r,rinv;
  double meff,ccel;
  int *ilist,*jlist,*numneigh,**firstneigh;


  double **x = atom->x;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double massone;
  double *radius = atom->radius;

  // contact area variables
  double overlap;
  double a;
  double sqrtb;
  double b;
  double contactArea;

  // charge transfer variables
  double *q = atom->q;
  double cut_off_norm = cut_off; //* 1e-6; // turn vaue in input script to micrometers
  double **previous_icontacts = atom->previous_icontacts;
  double current_icontact;
  double r_star;

  eflag = 0;
  egrav = 0.0;

  // Start of cell division algorithm
  int natoms = atom->natoms;
  int ai;
  int **N = atom->N;

  /* iterating over all particles and assigning a cell for each - what is the performace decrease for doing this every timestep
  split_box();

  for (ai = 0; ai < natoms; ai++){
    N[ai][0] = ceil(x[ai][0]/cell_width[0]);
    N[ai][1] = ceil(x[ai][1]/cell_width[1]);
    N[ai][2] = ceil(x[ai][2]/cell_width[2]);
  } */

  


  // iterate over particle neighbours

  inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
  firstneigh = pair_gran->list->firstneigh;

  for (ii = 0; ii < inum; ii++)
  {
    i = ilist[ii];
    xtmp = x[i][0]; // position x value
    ytmp = x[i][1]; // position y value
    ztmp = x[i][2]; // position z value
    radi = radius[i]; // radius of i

    jlist = firstneigh[i]; // first neighbour of i atom
    jnum = numneigh[i]; // number of neighbours

    for (jj = 0; jj < jnum; jj++)
    {
      bool finished_c = false;
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0]; //delta vector pointing from j to i
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz; 
      r = sqrt(rsq); // distance from i to j

      radj = radius[j];
      radsum = radi + radj;

      double delta_contactArea;
      // Calculate contact area:
      

      //std::cout << "Checking particle: " << i << " and particle: " << j << "\n";
      //std::cout << "Distance between i and j is: " << r <<"\n"; 

      // Calculate charge transfer:
      if (r < cut_off_norm){

        if (radsum > r){ // if the particles are in contact
        std::cout << "Distance between " << i << " and " << j << " is: " << r << "\n";
        std::cout << "Overlap of particles is: " << radsum - r << "\n"; 
        

        std::cout << "Previous before calc: " << *previous_icontacts[j] << "\n";
        
        current_icontact; // variable with current contact area

        //calcuate contact area
        // Version 1

        //overlap = radsum - r;
        //a = radi - (overlap * 0.5);
        //sqrtb = (r*r) - (a*a); // square root of radius of contact circle
        //b = sqrt(sqrtb); // radius of contact circle (pythagorian theorum)
        //contactArea = M_PI * (b*b);
        

        // Version 2 - rasera et al 2022
        
        r_star = (radius[i] * radius[j]) / (radius[i] + radius[j]);
        overlap = r - radsum;
        contactArea = M_PI * r_star * overlap;

        // std::cout << "Current: " << contactArea << "\n";
        

        // assigns contact area to variable 
        current_icontact = contactArea;

        
        // get deta contact area
        delta_contactArea = current_icontact - *previous_icontacts[j];

        std::cout << "Delta: " << delta_contactArea << "\n";
        
        // assign previous to current contact area
        *previous_icontacts[j] = current_icontact;

        std::cout << "previous after calc: " << *previous_icontacts[j] << "\n";
        std:: cout << "_______________________________________________________\n";
        
        }
        if (delta_contactArea > 0){
          // std::cout << "atoms: " << i << " and " << j << " made contact... \n";
          // std::cout << "Distance between them is: " << r << "\n" << "Distance for cut-off is: " << cut_off_norm << "\n";
          if (q[i] > q[j] && !finished_c){
            // std::cout << "[q]i > [q]j\n";
            // std::cout << "charge on i before: " << q[i] << "\n";
            // std::cout << "charge on j before: " << q[j] << "\n";
            
            // Testing
            q[j] += (q[i] * 0.2);
            q[i] *= 0.8;

            // std::cout << "charge on i after: " << q[i] << "\n";
            // std::cout << "charge on j after: " << q[j] << "\n";
            // std::cout << "________________________________________" << "\n";
            finished_c = true;
          }else if (q[i] < q[j] && !finished_c){
            // std::cout << "q[i] < q[j]\n";
            // std::cout << "charge on i before: " << q[i] << "\n";
            // std::cout << "charge on j before: " << q[j] << "\n";
            q[i] += (q[j] * 0.2);
            q[j] *= 0.8;
            // std::cout << "charge on i after: " << q[i] << "\n";
            // std::cout << "charge on j after: " << q[j] << "\n";
            // std::cout << "________________________________________" << "\n";
            finished_c = true;
          }else if (q[i] == q[j] && !finished_c){
            // std::cout << "q[i] == q[j]\n";
            // std::cout << "charge on i before: " << q[i] << "\n";
            // std::cout << "charge on j before: " << q[j] << "\n";
            q[i] += (q[j] / 4);
            q[j] -= (1.25 * q[i]);
            // std::cout << "charge on i after: " << q[i] << "\n";
            // std::cout << "charge on j after: " << q[j] << "\n";
            // std::cout << "________________________________________" << "\n";
            finished_c = true;
          }
        }
      }
    }
  }

  // f is force
  double coul_cut = neighbor->skin + cut_off_norm
  ;
  // loop over all atoms
  int ic, jc;
  double xtmpc, ytmpc, ztmpc, radic, radjc, delxc, delyc, delzc, rsqc, rc;

  double Fx, Fy, Fz;
  double Ke = 9e9;
  

  // iterate over particles on that processor
  for (ic = 0; ic < nlocal; ic++) 
  {
    if (mask[ic] & groupbit){
      xtmpc = x[ic][0]; // position x value
      ytmpc = x[ic][1]; // position y value
      ztmpc = x[ic][2]; // position z value
      radic = radius[ic]; // radius of i
    }

    for (jc = 0; jc < nlocal; jc++)
    {
      if (mask[jc] & groupbit){

      //std::cout << "checking particle: " << ic << " and particle: " << jc << "\n";

      delxc = xtmpc - x[jc][0]; //delta vector pointing from j to i
      delyc = ytmpc - x[jc][1];
      delzc = ztmpc - x[jc][2];
      rsqc = delxc*delxc + delyc*delyc + delzc*delzc; 
      rc = sqrt(rsqc); // distance from i to j

      //std::cout << "distance from " << ic << " to " << jc << " is: " << rc << "\n";

      radjc = radius[jc];
      
      if (q[jc] != 0 && rc < coul_cut && ic != jc){
        // coulombs law
        
        //std::cout << "coulomb interaction initiated for particles " << ic << " and " << jc << "\n"; 

        Fx = Ke*(q[ic] * q[jc])*(1/(delxc*delxc));
        Fy = Ke*(q[ic] * q[jc])*(1/(delyc*delyc));
        Fz = Ke*(q[ic] * q[jc])*(1/(delzc*delzc));

        //std::cout << "force experienced: " << Fx << " " << Fy << " " << Fz << "\n";

        // only apply force on 'i' particle because the loop will pass over 'j' particle and compute forces in opposing direction
        f[ic][0] += Fx;
        f[ic][1] += Fy;
        f[ic][2] += Fz;
        }
      }
    }
  }
}

// Split simulation box
/* Optimisation to be added later
void FixTribo::split_box(){

  // (x, y, z) of lower bound of domain
  xboxlo = domain->boxlo[0];
  yboxlo = domain->boxlo[1];
  zboxlo = domain->boxlo[2];

  // (x, y, z) of upper bound of domain
  xboxhi = domain->boxhi[0];
  yboxhi = domain->boxhi[1];
  zboxhi = domain->boxhi[2];

  // size of simulation box
  delta_box_x = xboxhi - xboxlo;
  delta_box_y = yboxhi - yboxlo;
  delta_box_z = zboxhi - zboxlo;

  

  // change these to change divisions -> update later in order to make the divisions dynamic
  int ncells_x = 4;
  int ncells_y = 2;
  int ncells_z = 4;

  // array storing values of cell width
    
  cell_width[0] = delta_box_x / ncells_x;
  cell_width[1] = delta_box_y / ncells_y;
  cell_width[2] = delta_box_z / ncells_z;
}*/

/* ---------------------------------------------------------------------- */

void FixTribo::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTribo::set_acceleration()
{
  if (style == CHUTE || style == SPHERICAL) {
    if (style == CHUTE) {
      phi = 0.0;
      theta = 180.0 - vert;
    }
    if (domain->dimension == 3) {
      xgrav = sin(degree2rad * theta) * cos(degree2rad * phi);
      ygrav = sin(degree2rad * theta) * sin(degree2rad * phi);
      zgrav = cos(degree2rad * theta);
    } else {
      xgrav = sin(degree2rad * theta);
      ygrav = cos(degree2rad * theta);
      zgrav = 0.0;
    }
  } else if (style == VECTOR) {
    if (domain->dimension == 3) {
      double length = sqrt(xdir*xdir + ydir*ydir + zdir*zdir);
      if(length == 0.)
        error->one(FLERR,"Gravity direction vector = 0");
      xgrav = xdir/length;
      ygrav = ydir/length;
      zgrav = zdir/length;
    } else {
      double length = sqrt(xdir*xdir + ydir*ydir);
      if(length == 0.)
        error->one(FLERR,"Gravity direction vector = 0");
      xgrav = xdir/length;
      ygrav = ydir/length;
      zgrav = 0.0;
    }
  }

  xacc = magnitude*xgrav;
  yacc = magnitude*ygrav;
  zacc = magnitude*zgrav;
}

/* ----------------------------------------------------------------------
   potential energy in gravity field
------------------------------------------------------------------------- */

double FixTribo::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(&egrav,&egrav_all,1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return egrav_all;
}

/* ---------------------------------------------------------------------- */

void FixTribo::get_tribo(double *grav)
{
    grav[0] = xgrav * magnitude;
    grav[1] = ygrav * magnitude;
    grav[2] = zgrav * magnitude;
}
/*
void FixTribo::varparse(char *name, int m)
{
  varflag = 1;

  name = &name[2];
  int n = strlen(name) + 1;
  char *str = new char[n];
  strcpy(str,name);

  int ivar = input->variable->find(str); // find variable with that name
  delete [] str; // save memory

  if (ivar < 0) // cant find variable with that name
    error->all(FLERR,"Variable name for set command does not exist");
  if (!input->variable->atomstyle(ivar))
    error->all(FLERR,"Variable for set command is invalid style");

  if (m == 1) {
    varflag1 = 1; ivar1 = ivar;
  } else if (m == 2) {
    varflag2 = 1; ivar2 = ivar;
  } else if (m == 3) {
    varflag3 = 1; ivar3 = ivar;
  } else if (m == 4) {
    varflag4 = 1; ivar4 = ivar;
  }
}*/
