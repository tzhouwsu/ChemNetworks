// -*- c++ -*-

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "fix_chemnetwork.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "universe.h"
#include "update.h"
#include "citeme.h"

#define NUM_CN_ARGS 4

static const char colvars_pub[] =
  "fix chemnetwork command:\n\n"
  "@Article{chemnetwork,\n"
  " author =  {authors},\n"
  " title =   {title},\n"
  " journal = {journal},\n"
  " year =    year,\n"
  " note =    {doi: }\n"
  "}\n\n";

/***************************************************************/

using namespace LAMMPS_NS;
using namespace FixConst;

FixChemNetwork::FixChemNetwork(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // if(comm->nprocs > 1)
  //   error->all(FLERR,"Multiple MPI ranks not yet supported");
  
  if (narg < 5)
    error->all(FLERR,"Illegal fix chemnetwork command: too few arguments");

  int iarg = 3;
  
  filename_cn_input = new char;
  strcpy(filename_cn_input, arg[iarg++]);

  nevery = force->inumeric(FLERR,arg[iarg++]);

  int ntypes = atom->ntypes;
  int found_ntypes = narg - iarg;

  if(ntypes != found_ntypes) error->all(FLERR,"Incorrect ntypes in fix command");

  memory->create(map, ntypes+1, CN_ELEMENT_SIZE, "fix_chemnetwork::map");
  for(int i=1; i<ntypes+1; ++i) strcpy(map[i], arg[iarg++]);

  me = comm->me;

  // create instance of ChemNetworks on master rank

  cn = NULL;
  if(me == 0) cn = new ChemNetworkNew();
  
  cn_argv = NULL;
  
  if (lmp->citeme) lmp->citeme->add(colvars_pub);

  num_atoms = 0;
  type_all = NULL;
  x_all = NULL;

  memory->create(t_time,  NUM_TIMERS, "fix_chemnetwork::t_time");
  memory->create(t_count, NUM_TIMERS, "fix_chemnetwork::t_count");

  for(int i=0; i<NUM_TIMERS; ++i) {
    t_time[i]  = 0.0;
    t_count[i] = 0;
  }
}

/*********************************
 * Clean up on deleting the fix. *
 *********************************/
FixChemNetwork::~FixChemNetwork()
{
  if(cn) delete cn;
  delete [] cn_argv;

  if(x_all) memory->destroy(x_all);
  if(type_all) memory->destroy(type_all);

  write_timing_summary();
  
  memory->destroy(t_time);
  memory->destroy(t_count);
}

/* ---------------------------------------------------------------------- */

int FixChemNetwork::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

// initial checks for colvars run.

void FixChemNetwork::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use fix colvars without atom IDs");

  if (atom->map_style == 0)
    error->all(FLERR,"Fix colvars requires an atom map, see atom_modify");
}

/* ---------------------------------------------------------------------- */

void FixChemNetwork::setup(int vflag)
{
  const double t0 = MPI_Wtime();
  
  // read ChemNetworks input file
  
  cn_argv = new char * [NUM_CN_ARGS];
  cn_argv[0] = (char *) "cn.exe";
  cn_argv[1] = (char *) filename_cn_input;
  cn_argv[2] = (char *) "lmp_fix";
  cn_argv[3] = (char *) "-new";
  
  if(cn) cn->input_read(NUM_CN_ARGS, cn_argv);

  // allocate buffer for all particles in system (yeah, I know not scalable...)

  num_atoms = (int) atom->natoms;
  
  if(me == 0) memory->create(type_all, num_atoms, CN_ELEMENT_SIZE, "fix_chemnetwork::type_all");
  memory->create(x_all, 3 * num_atoms, "fix_chemnetwork::x_all");

  const tagint * const tag = atom->tag;
  const int * const type = atom->type;

  // Since we're not expecting total number of particles to change nor atoms
  // to change type, we'll populate this type array once at start

  int * buf;
  memory->create(buf, num_atoms+1, "fix_chemnetwork::buf");

  for(int i=0; i<num_atoms+1; ++i) buf[i] = 0;
  for(int i=0; i<atom->nlocal; ++i) buf[ tag[i] ] = type[i];

  if(me == 0) MPI_Reduce(MPI_IN_PLACE, buf, num_atoms+1, MPI_INT, MPI_SUM, 0, world);
  else MPI_Reduce(buf, buf, num_atoms+1, MPI_INT, MPI_SUM, 0, world);
  
  if(me == 0) for(int i=1; i<num_atoms+1; ++i) strcpy(type_all[i-1], map[ buf[i] ]);
  
  memory->destroy(buf);

  const double t1 = MPI_Wtime();

  t_time[SETUP] += t1 - t0;
  t_count[SETUP]++;
  
  // process initial configuration in ChemNetworks

  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixChemNetwork::end_of_step()
{
  const double t0 = MPI_Wtime();
  
  // map coordinates into ChemNetwork format
  
  const tagint * const tag = atom->tag;
  const double * const * const x = atom->x;

  // pack coordinates into global buffer

  for(int i=0; i<num_atoms*3; ++i) x_all[i] = 0.0;
  
  for(int i=0; i<atom->nlocal; ++i) {
    const int id = (int) tag[i] - 1;
    x_all[id*3  ] = x[i][0];
    x_all[id*3+1] = x[i][1];
    x_all[id*3+2] = x[i][2];
  }

  if(me == 0) MPI_Reduce(MPI_IN_PLACE, x_all, 3*num_atoms, MPI_DOUBLE, MPI_SUM, 0, world);
  else MPI_Reduce(x_all, x_all, 3*num_atoms, MPI_DOUBLE, MPI_SUM, 0, world);

  const double t1 = MPI_Wtime();

  t_time[PACK] += t1 - t0;
  t_count[PACK]++;
  
  if(cn) {
    
    // load configuration into ChemNetworks
    
    cn->load_xyz_solvent1(num_atoms, type_all, x_all);

    // execute ChemNetworks on this config
    
    cn->process_config(cn_argv);
    
    const double t2 = MPI_Wtime();

    t_time[COMPUTE]+= t2 - t1;
    t_count[COMPUTE]++;
  }

}

/* ---------------------------------------------------------------------- */
/* local memory usage. approximately. */
double FixChemNetwork::memory_usage(void)
{
  double bytes = 3 * num_atoms * sizeof(double);
  bytes += num_atoms * CN_ELEMENT_SIZE * sizeof(char);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixChemNetwork::write_timing_summary()
{
  if(me > 0) return;

  if(screen) {
    fprintf(screen,"\nChemNetwork Timings on Rank 0\n");
    fprintf(screen,"SETUP::   time= %f  count= %i  per_call= %f \n",
	    t_time[SETUP],t_count[SETUP],t_time[SETUP]/t_count[SETUP]);
    
    fprintf(screen,"PACK::    time= %f  count= %i  per_call= %f \n",
	    t_time[PACK],t_count[PACK],t_time[PACK]/t_count[PACK]);
    
    fprintf(screen,"COMPUTE:: time= %f  count= %i  per_call= %f \n",
	    t_time[COMPUTE],t_count[COMPUTE],t_time[COMPUTE]/t_count[COMPUTE]);
  }

  if(logfile) {
    fprintf(logfile,"\nChemNetwork Timings on Rank 0\n");
    fprintf(logfile,"SETUP::   time= %f  count= %i  per_call= %f \n",
	    t_time[SETUP],t_count[SETUP],t_time[SETUP]/t_count[SETUP]);
    
    fprintf(logfile,"PACK::    time= %f  count= %i  per_call= %f \n",
	    t_time[PACK],t_count[PACK],t_time[PACK]/t_count[PACK]);
    
    fprintf(logfile,"COMPUTE:: time= %f  count= %i  per_call= %f \n",
	    t_time[COMPUTE],t_count[COMPUTE],t_time[COMPUTE]/t_count[COMPUTE]);
  }
}
