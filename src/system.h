// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// system.h contains a list of properties needed to define a system of protein
// hard spheres

#ifndef HYBRIDMC_SYSTEM_H
#define HYBRIDMC_SYSTEM_H

#include "vec3.h"
#include "molecule.h"
#include <vector>

class System {
  public:
      //
      int numAtoms;
      // positions of all beads at time t
      std::vector<Vec3> pos;
      // velocities of all beads at time t
      std::vector<Vec3> vel;
      // bead clocks
      std::vector<double> times;
      // collision counters of beads
      std::vector<uint64_t> counter;
      // entropy of beads
      std::vector<double> s_bias;

      // flag to indicate whether to write distance data
      bool distanceWrite;

      bool useEnsemble;
      bool recordEnsemble;
      int ensembleSize;
      std::vector<Molecule> ensemble;
      std::vector<Molecule> nextEnsemble;
      std::vector<Molecule> trappedEnsemble;
      std::vector<bool> trapped; // flag to classify each member of "ensemble" as trapped or not

      System(unsigned int nbeads)
          : numAtoms(nbeads), pos(nbeads), vel(nbeads), times(nbeads), counter(nbeads), useEnsemble(false) {}

};

#endif
