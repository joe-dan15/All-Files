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

#ifdef FIX_CLASS

FixStyle(tribo,FixTribo)

#else

#ifndef LMP_FIX_TRIBO_H
#define LMP_FIX_TRIBO_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTribo : public Fix {
  friend class FixPour;

 public:
  FixTribo(class LAMMPS *, int, char **);
  ~FixTribo();
  int setmask();
  void init();
  void setup(int);
  virtual void post_force(int);
  virtual void post_force_respa(int, int, int);
  double compute_scalar();
  void split_box();
  void build_nei(int);

  void get_tribo(double*); 
  void varparse(char *, int); // parse a var name to a value

private:
  class PairGran *pair_gran;

 protected:

  // new variables
  char material;
  double cut_off;
  char *costr;
  double covar;
  double cut_off_dist;

  bool found_max_rad = false;

  double xboxlo, yboxlo, zboxlo, xboxhi, yboxhi, zboxhi, delta_box_x, delta_box_y, delta_box_z;
  double *cell_width;

  double **prev_contact;
  double **current_contact;
  


  // previous variables
  int style;
  double magnitude;
  double vert,phi,theta;
  double xdir,ydir,zdir;
  double xgrav,ygrav,zgrav,xacc,yacc,zacc;
  double degree2rad;
  int nlevels_respa;
  int time_origin;
  int eflag;
  double egrav,egrav_all;

  int varflag;
  int mstyle,vstyle,pstyle,tstyle,xstyle,ystyle,zstyle;
  int mvar,vvar,pvar,tvar,xvar,yvar,zvar;
  char *mstr,*vstr,*pstr,*tstr,*xstr,*ystr,*zstr;

  void set_acceleration();
  class FixMultisphere *fm;
  class FixRelaxContacts *fix_relax;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for fix gravity does not exist

Self-explanatory.

E: Variable for fix gravity is invalid style

Only equal-style variables can be used.

*/
