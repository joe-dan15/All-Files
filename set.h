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

#ifdef COMMAND_CLASS

CommandStyle(set,Set)

#else

#ifndef LMP_SET_H
#define LMP_SET_H

#include "pointers.h"

namespace LAMMPS_NS {

class Set : protected Pointers {
 public:
  Set(class LAMMPS *lmp) : Pointers(lmp) {};
  void command(int, char **);

 private:
  char *id;
  int *select;
  int style,ivalue,newtype,count,index_custom;
  int ximage,yimage,zimage,ximageflag,yimageflag,zimageflag;
  double wfvalue, dvalue,xvalue,yvalue,zvalue,wvalue,fraction;

  int varflag,varflag1,varflag2,varflag3,varflag4;
  int ivar1,ivar2,ivar3,ivar4;
  double *vec1,*vec2,*vec3,*vec4;

  class FixPropertyAtom* updFix; 
  int nUpdValues; 
  double *updValues; 
  int add;
  bigint until; 
  bigint currentTimestep; 

  void selection(int);
  void set(int);
  void setrandom(int);
  void topology(int);
  void varparse(char *, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Set command before simulation box is defined

The set command cannot be used before a read_data, read_restart,
or create_box command.

E: Set command with no atoms existing

No atoms are yet defined so the set command cannot be used.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid value in set command

The value specified for the setting is invalid, likely because it is
too small or too large.

E: Invalid random number seed in set command

Random number seed must be > 0.

E: Cannot set this attribute for this atom style

The attribute being set does not exist for the defined atom style.

E: Invalid mass in set command

Self-explanatory.

E: Invalid shape in set command

Self-explanatory.

E: Invalid length in set command

Self-explanatory.

E: Invalid dipole length in set command

Self-explanatory.

E: Invalid diameter in set command

Self-explanatory.

E: Cannot set non-zero image flag for non-periodic dimension

Self-explanatory.

E: Cannot set meso_rho for this atom style

Self-explanatory.

E: Cannot use set atom with no atom IDs defined

Atom IDs are not defined, so they cannot be used to identify an atom.

E: Cannot use set mol with no molecule IDs defined

Self-explanatory.

E: Could not find set group ID

Group ID specified in set command does not exist.

E: Set region ID does not exist

Region ID specified in set command does not exist.

E: Cannot set quaternion for atom that has none

Self-explanatory.

E: Cannot set theta for atom that is not a line

Self-explanatory.

E: Bond atom missing in set command

The set command cannot find one or more atoms in a particular bond on
a particular processor.  The pairwise cutoff is too short or the atoms
are too far apart to make a valid bond.

E: Angle atom missing in set command

The set command cannot find one or more atoms in a particular angle on
a particular processor.  The pairwise cutoff is too short or the atoms
are too far apart to make a valid angle.

E: Dihedral atom missing in set command

The set command cannot find one or more atoms in a particular dihedral
on a particular processor.  The pairwise cutoff is too short or the
atoms are too far apart to make a valid dihedral.

E: Improper atom missing in set command

The set command cannot find one or more atoms in a particular improper
on a particular processor.  The pairwise cutoff is too short or the
atoms are too far apart to make a valid improper.

*/
