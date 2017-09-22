/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef FIX_CLASS

FixStyle(chemnetwork,FixChemNetwork)

#else

#ifndef LMP_FIX_CHEMNETWORK_H
#define LMP_FIX_CHEMNETWORK_H

#include "fix.h"
#include "chemnetworks_new.h"

#define CN_ELEMENT_SIZE 32

using namespace CN_NS;

namespace LAMMPS_NS {

class FixChemNetwork : public Fix {

 public:
  FixChemNetwork(class LAMMPS *, int, char **);
  virtual ~FixChemNetwork();

  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void end_of_step();
  virtual double memory_usage();

 protected:
  enum ttype {SETUP, PACK, COMPUTE, NUM_TIMERS};

  void write_timing_summary();
  
  int me;
  
  class ChemNetwork * cn;

  char ** cn_argv;
  char * filename_cn_input;
  char ** map;

  int num_atoms;
  char ** type_all;
  double * x_all;

  double * t_time;
  int * t_count;
};

}

#endif
#endif

/* ERROR/WARNING messages: */
