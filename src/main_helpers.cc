// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// main.cc simulates the dynamics of protein folding using a coarse-grained
// model; each amino acid in the protein is represented by a bead, and the
// beads are connected by local and nonlocal bonds; the dynamics is
// event-driven, and the program converges when each state of the protein has
// been visited roughly equally, giving the final value of the entropy of each
// state;

#include "main_helpers.h"
#include "internalState.h"


double max_time = std::numeric_limits<double>::max();

std::string extractBaseName(const std::string& filename) {
    // Find the position of the last period character
    size_t periodPos = filename.find_last_of('.');
    // Check if the period character is found
    if (periodPos != std::string::npos) {
        // Find the position of the last underscore character before the period
        size_t underscorePos = filename.find_last_of('_', periodPos);
        // Check if the underscore character is found
        if (underscorePos != std::string::npos) {
            // Extract the substring from index 0 to underscorePos
            return filename.substr(0, underscorePos);
        }
    }
    // Return an empty string if extraction fails
    return "";
}


void writeFrustratedConfigurations(const System &sys, std::string h5_name)
{

    std::string base_name = extractBaseName(h5_name);

    std::string untrapped_name = base_name + ".untrapped.csv";
    std::string trapped_name = base_name + ".trapped.csv";
    std::ofstream untrappedFile(untrapped_name, std::ios_base::app);
    std::ofstream trappedFile(trapped_name, std::ios_base::app);

    std::cout << " Writing out frustrated and trapped structure to files " << untrapped_name << " and " << trapped_name << std::endl;
    for (size_t i = 0; i < sys.trappedEnsemble.size(); ++i){
        DihedralState dihedrals(sys.numAtoms, sys.trappedEnsemble[i].getPositions() ); // form internal angles from Cartesian positions of configuration x
        dihedrals.outputInternalAngles(trappedFile, true);
    }
    //
    //  After replacement, none of the initial ensemble of states is trapped
    //
    for (size_t i = 0; i < sys.ensemble.size(); ++i){

            DihedralState dihedrals(sys.numAtoms, sys.ensemble[i].getPositions() ); // form internal angles from Cartesian positions of configuration x
            dihedrals.outputInternalAngles(untrappedFile, false);
    }
}

double checkFrustration(System &sys, Random &mt, const Param &p, const Box &box)
{
    //  Check to see if a long trajectory starting from each ensemble state can reach
    //  the target state in the next layer
    //
    int numFails = 0;
    UpdateConfig update_config;
    Molecule current(p.nbeads);
    current.setPositions(sys.pos);
    std::cout << std::endl << "  Checking ensemble states for frustration: num states is "
        << sys.ensembleSize << " = " << sys.ensemble.size() << std::endl;

    //  find maximum distance of transient bond in ensemble
    double max_init_distance = 0.0;
    for (int i=0;i<sys.ensembleSize;i++){
        std::vector<Vec3> pos = sys.ensemble[i].getVec3Positions();
        double d_i = compute_transient_dist(pos, p, box);
        if (d_i > max_init_distance) max_init_distance = d_i;
        //std::cout << " initial d_i = " << d_i << " and max = " << max_init_distance << std::endl;
    }


    std::vector<bool> trapped(sys.ensembleSize);
    sys.trappedEnsemble.clear();
    sys.trapped.clear();

    std::vector<int> goodIndices;
    for (int e=0;e<sys.ensembleSize;e++){

        sys.pos = sys.ensemble[e].getVec3Positions();
        init_update_config(sys.pos, update_config, box, p.transient_bonds);

        bool trappedIndicator = run_anneal(sys, mt, p, box, max_init_distance);

        if (trappedIndicator){
            numFails++;
            sys.trappedEnsemble.push_back(sys.ensemble[e]);
            sys.trapped.push_back(true);
            std::cout << " Ensemble member " << e << " does not reach target." << std::endl;
        } else {
            std::cout << " Ensemble member " << e << " can flip state" << std::endl;
            sys.trapped.push_back(false);
            goodIndices.push_back(e); // save index if not trapped
        }

    }

    double frustration_factor = double(numFails)/sys.ensembleSize;
    std::cout << "  Frustration fraction is " <<  frustration_factor << std::endl;

    //
    // replace trapped states with untrapped states in ensemble
    //
    if (frustration_factor == 1.0){
	    throw std::runtime_error("No states in ensemble reach the target.");
    }

    if (frustration_factor > 0.0){

        int s_index = 0;
        for (int e=0;e<sys.ensembleSize;e++){
            if (sys.trapped[e]){
                std::cout << " Replacing trapped state " << e << " with untrapped state " << goodIndices[s_index]
                    << std::endl;

                sys.ensemble[e] = sys.ensemble[ goodIndices[s_index++] ];
                sys.trapped[e] = false;
                if (s_index >= (int )goodIndices.size() ) s_index = 0;

            }
        }
    }

    sys.pos = current.getVec3Positions(); // restore position
    init_update_config(sys.pos, update_config, box, p.transient_bonds);

    for (unsigned int i = 0; i < p.nbeads; i++)
    {
        sys.times[i] = 0.0;
        sys.counter[i] = 0.0;
    }
    return frustration_factor;

}
void checkEnsemble(System &sys, const Box &box, const NonlocalBonds &transient_bonds)
{
   //std::cout << " In check Ensemble, update_config is " << update_config.config << std::endl;

   UpdateConfig test_config;
   double count_bonded = 0.0;
   for (int e=0;e<sys.ensembleSize;e++)
   {
        std::vector<Vec3> pos_e = sys.ensemble[e].getVec3Positions();
        init_update_config(pos_e, test_config, box, transient_bonds) ;
        if (test_config.config == 1) count_bonded += 1.0;
        //std::cout << " Ensemble member " << i << " has config = " << test_config.config << std::endl;
   }
   std::cout << "  Fraction of bonded states in ensemble was " << count_bonded/sys.ensembleSize << std::endl;
}

void generateEnsemble(System &sys, Random &mt, const Param &p, const Box &box)
{
    std::cout << " generateEnsemble called with ensemble size of " << p.ensembleSize
        << " for molecule of length " << p.nbeads << std::endl;
    std::vector<Molecule> ensemble(p.ensembleSize);
    std::vector<Vec3> pos( p.nbeads );

    UpdateConfig test_config;

    for (int i=0;i<p.ensembleSize;i++){
        // set flag for random initialization success
        bool found = false;
        // set the maximum tries for a random initialization to 10
        int max_init_count = 100;
        // while random initialization not found and 10 attempts not hit, try random init
        while (found == false && max_init_count){
            found = init_pos(pos, box, mt, p);
            max_init_count--;
            test_config = config_int(pos, box, p.transient_bonds);
            if (test_config.config != 0) {
                found = false;
                std::cout << " Generated initially bound configuration.  Skipping." << std::endl;
            }
        }
        //std::cout << " Found configuration starting at " << pos[0] << " ending at " << pos[p.nbeads-1]
        //    << std::endl;
        Molecule mol(p.nbeads);
        mol.setPositions(pos);

        ensemble[i] = mol;
        //mol.printPositions(i);
    }

    sys.ensemble = ensemble;
    checkEnsemble(sys, box, p.transient_bonds);
    //printEnsemble(sys.ensemble);

}


void initialize_pos(System &sys, Random &mt, const Param &p, const Box &box,
                    UpdateConfig &update_config,
                    std::optional<std::string> input_name,
                    std::optional<std::string> snapshot_name,
                    const unsigned int t_bonds) {

  if (snapshot_name && std::filesystem::exists(*snapshot_name)) {
    // overwrite existing entries in pos and s_bias vectors with read-in values
    // from hdf5 file
    std::cout << " Reading in snapshot " << *snapshot_name << std::endl;
    read_snapshot(*snapshot_name, sys.pos, sys.s_bias, mt, update_config);
  } else if (input_name && std::filesystem::exists(*input_name)) {
    std::cout << "Reading in input file " << *input_name << std::endl;
    read_input(*input_name, sys.pos);

    init_update_config(sys.pos, update_config, box, p.transient_bonds);
    init_s(sys.s_bias, t_bonds);
    double f_factor = 0.0;
    if (sys.useEnsemble)
    {
        sys.ensemble = readMoleculesFromHDF5ByName(input_name->c_str());
        checkEnsemble(sys, box, p.transient_bonds);
        f_factor = checkFrustration(sys, mt, p, box); // only do if a problem appears?
        init_update_config(sys.pos, update_config, box, p.transient_bonds);
    }

    if (f_factor > 0.0) std::cerr << "  Note:  frustration exists." << std::endl;


  } else {

    // set flag for random initialization success
    bool found = false;
    // set the maximum tries for a random initialization to 10
    int max_init_count = 10;
    // while random initialization not found and 10 attempts not hit, try random init
    while (found == false && max_init_count)
    {
        found = init_pos(sys.pos, box, mt, p);
        max_init_count--;
    };

    if (sys.useEnsemble)
    {
        generateEnsemble(sys, mt, p, box); // generate the initial ensemble of configurations
        checkEnsemble(sys, box, p.transient_bonds);
    }

    // do linear draw if above procedure failed

    if (found == false) {

        std::cout << " Calling linear chain draw." << std::endl;
        draw_linear_chain(sys.pos,p);

    }

    init_s(sys.s_bias, t_bonds);
  }
}

void initialize_system(System &sys, Random &mt, const Param &p, const Box &box,
                       UpdateConfig &update_config, Cells &cells,
                       EventQueue &event_queue) {

  LOG_DEBUG("New step initialization");
  if (!check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                        p.nnear_max2)) {
    throw std::runtime_error("local beads overlap in initialize_system.");
  }

  if (!check_nonlocal_dist(sys.pos, box, p.rh2, p.stair2,
                           p.transient_bonds, p.permanent_bonds)) {
    throw std::runtime_error("nonlocal beads overlap in initialize_system");
  }

  if (cells.ncell < 4) {
    throw std::invalid_argument("bead ncell must be at least 4");
  }

  double check_rc = 0;
  for (unsigned int i=0; i < p.nonlocal_bonds.get_nbonds(); i++) {

      if (p.nonlocal_bonds.getrc(i) > check_rc) {
          check_rc = p.nonlocal_bonds.getrc(i);
      }
  }

  if (p.stair) {
    check_rc = *p.stair;
  }

  if (cells.lcell < check_rc) {
    std::cout << cells.lcell << " is cells.lcell and the check_rc is " << check_rc << std::endl;
    float req_length = cells.ncell * check_rc + 0.0001; // NOLINT(cppcoreguidelines-narrowing-conversions)
    std::cout << "Length has to be at least " << req_length << std::endl;
    throw std::invalid_argument("bead lcell must be at least rc. Increase length to at least above value");
  }

  // split the box into smaller cells, and store the beads in each cell
  init_cells(sys.pos, box, cells);

  init_vel(sys.vel, mt, p.temp, p.m);

  // fill priority queue with cell crossings of all particles
  init_cell_events(sys.pos, sys.vel, p.nbeads, box, sys.counter, event_queue,
                   sys.times, cells);

  // fill priority queue with nearest bond events between all particle pairs
  init_nearest_bond_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                           event_queue, sys.times, p.near_min2, p.near_max2);

  // fill priority queue with next-nearest bond events between all particle
  // pairs
  init_nnearest_bond_events(sys.pos, sys.vel, p.nbeads, box, sys.counter,
                            event_queue, sys.times, p.nnear_min2, p.nnear_max2);

  // fill priority queue with collisions of all particle pairs (executed once,
  // with an initial prediction that fills the entire priority queue)
  add_events_for_all_beads(sys.pos, sys.vel, p.nbeads, p.rh2, p.stair2, box, sys.counter, event_queue, sys.times,
                           cells, p.transient_bonds, p.permanent_bonds,
                           update_config, p.max_nbonds);

}

void run_step(System &sys, const Param &p, const Box &box,
              UpdateConfig &update_config, CountBond &count_bond,
              Cells &cells, EventQueue &event_queue,
              const unsigned int step, double del_t) {
  LOG_DEBUG("step = " << step);

  double zero_time = 0.0;
  // the current time interval the events are occurring in
  double step_time = step * del_t;

//  std::cout << " At start of step " << step << " local clock time is " << sys.times[0] << std::endl;

  // while events are occurring in step_time,
  while (!event_queue.empty()) {
    // access the minimum time to the next collision, and the indices and
    // collision counters of the beads associated with this collision
    const Event event = event_queue.top();

    if (std::visit(
            [=](auto &&ev) {
              // check for monotonically increasing event times
              assert(ev.t >= zero_time);
              return ev.t > step_time;
            },
            event))
      break;

    event_queue.pop();

    // process collision or cell crossing event
    std::visit(
        [&](auto &&ev) {
          zero_time = ev.t;
          LOG_DEBUG("wall time " << zero_time << " Queue Size is " << event_queue.size());
          process_event(ev, sys, p, box, event_queue, cells, update_config,
                        count_bond);
        },
        event);
  }

  // update positions at the moment the last collision in del_t occurred
  for (unsigned int i = 0; i < p.nbeads; i++) {
    update_pos(sys.pos[i], sys.vel[i], sys.times[i], step_time);
    assert(check_overlap(i, sys.pos, sys.vel, sys.times, p.rh2, box));
  }

  // check bond distances of nearest and next-nearest beads
  assert(check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                          p.nnear_max2));
  // check bond distances of nonlocal beads
  assert(check_nonlocal_dist(sys.pos, box, p.rh2, p.stair2,
                             p.transient_bonds, p.permanent_bonds));

}

void run_trajectory_eq(System &sys, Random &mt, const Param &p, const Box &box,
                       UpdateConfig &update_config, CountBond &count_bond,
                       unsigned int iter,
                       DistWriter &dist_writer, std::vector<double> &dist) {

  LOG_DEBUG("run_trajectory_eq");
  for (unsigned int step = iter * p.nsteps; step < (iter + 1) * p.nsteps;
       step++) {
    EventQueue event_queue;
    Cells cells{p.ncell, p.length / p.ncell};

    //set max time
    if (step != 0) {max_time = (step * p.del_t) + 0.001;}

    LOG_DEBUG("max_time = " << max_time);

    initialize_system(sys, mt, p, box, update_config, cells, event_queue);
       // to check energy conservation
    const double tot_E_before =
        compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);

    run_step(sys, p, box, update_config, count_bond, cells,
             event_queue, step, p.del_t);

    const double tot_E_during =
        compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);
    const double E_diff = std::abs(1 - (tot_E_during / tot_E_before));
    if (E_diff >= 1e-6) {
      std::cout << E_diff << " energy difference in equilibration iter " << iter << " step = " << step << std::endl;
      throw std::runtime_error("energy is not conserved in equilibration: steps = ");
    }

  }

  dist_between_nonlocal_beads(sys.pos, box, p.nonlocal_bonds, dist);
  if (sys.distanceWrite) dist_writer.append(dist);

  assert(check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                          p.nnear_max2));
  assert(check_nonlocal_dist(sys.pos, box, p.rh2, p.stair2,
                             p.transient_bonds, p.permanent_bonds));
}

void run_trajectory(System &sys, Random &mt, const Param &p, const Box &box,
                    std::vector<double> &dist, UpdateConfig &update_config,
                    UpdateConfigWriter &update_config_writer,
                    PosWriter &pos_writer, VelWriter &vel_writer,
                    ConfigWriter &config_writer, DistWriter &dist_writer,
                    std::set<Config> &store_config, ConfigInt &store_config_int,
                    CountBond &count_bond,
                    unsigned int iter,bool storeTrajectory) {

  LOG_DEBUG("run_trajectory");
   //  Do swap MC move if active
  if (sys.useEnsemble and (update_config.config == 0) )
  {
      //std::cout << "  Doing swap MC move. " << std::endl;
      swapMC(sys, mt, box, p);
      UpdateConfig trial_config = config_int(sys.pos, box, p.transient_bonds);
      if (trial_config.config != update_config.config)
      {

        std::ostringstream errorMessage;
        errorMessage << "Possible error with swap move. State is " << trial_config.config
                     << " and was " << update_config.config << std::endl;

          throw std::runtime_error(errorMessage.str());
      }
  }



  // Do Molecular Dynamics moves

  // assume that the entire time during which the beads are undergoing events
  // can be divided into intervals of length p.del_t; the total number of such
  // intervals is p.nsteps (thus, the variable called step marks the intervals
  // until the p.nsteps-th interval is reached); each iteration of the
  // BIIIIIIIG loop will run until update_time and then dump the output to the
  // hdf5 file



  for (unsigned int step = iter * p.nsteps; step < (iter + 1) * p.nsteps;
       step++) {
    EventQueue event_queue;
    Cells cells{p.ncell, p.length / p.ncell};

    //set max time
    if (step != 0) {max_time = (step * p.del_t) + 0.001;}

/*     std::cout << " step = " << step
        << " max_time = " << max_time << " with nsteps = " << p.nsteps
        << std::endl;*/

    initialize_system(sys, mt, p, box, update_config, cells, event_queue);

    // to check energy conservation
    const double tot_E_before =
        compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);

    run_step(sys, p, box, update_config, count_bond, cells,
             event_queue, step, p.del_t);
    if (sys.distanceWrite)
    {
        dist_between_nonlocal_beads(sys.pos, box, p.nonlocal_bonds, dist);
        dist_writer.append(dist);
    }

    if (step % p.write_step == 0) {
      // store the integer of the configuration and the time of the event

      //store_config_int.emplace_back(update_config.config);
      //update_config_writer.config_int.emplace_back(update_config.config);

      // if a transient bond forms, check if configuration has been previously
      // visited by comparing configuration to set of saved configurations; if
      // not visited, save configuration to the set and write positions of
      // beads to file
      if (update_config.config == 1 &&
          store_config.insert(update_config.config).second) {
        //
        //  The insertion into the set of Config above fails if store_config already contains update_config.config (i.e 1 here)
        //   The insert(obj).second is a boolean that is true if insertion is successful.
        //
        pos_writer.append(sys.pos);
        vel_writer.append(sys.vel);
        config_writer.append(update_config.config);
      }

    }

    const double tot_E_during =
        compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);
    const double E_diff = std::abs(1 - (tot_E_during / tot_E_before));
    if (E_diff >= 1e-6) {
      std::cout << E_diff << " energy difference" << std::endl;
      throw std::runtime_error("energy is not conserved");
    }
  }
  //
  //  Now add to ensemble if needed

  if ( ((int)sys.nextEnsemble.size() < sys.ensembleSize) and sys.recordEnsemble )
  {
      if ( (update_config.config == 1) and (iter % p.ensemble_write_step == 0) ){
            //std::cout << " @@@@@ Adding configuration to stored ensemble at iter " << iter
            //    << " size is " << sys.nextEnsemble.size() << " of " << sys.ensembleSize << std::endl;
            Molecule mol(p.nbeads);
            mol.setPositions( sys.pos );
            sys.nextEnsemble.push_back(mol);
      }
  }

  // Do Monte Carlo moves (mc moves)
  for (unsigned int i = 0; i < p.mc_moves; i++) {
      crankshaft(sys.pos, update_config, box, p.near_min2, p.near_max2,
                 p.nnear_min2, p.nnear_max2, p.rh2, p.stair2,
                 p.transient_bonds, p.permanent_bonds, mt, sys.s_bias);
      /* for (unsigned int j = 1; j < p.mc_write; j++) {
          if (i == ((j * p.mc_moves) / p.mc_write)) {
              update_config_writer.config_int.emplace_back(update_config.config);
          }
      }*/
  }

  store_config_int.emplace_back(update_config.config); // record the current state in vector list

  if (storeTrajectory)
  {

      update_config_writer.config_int.emplace_back(update_config.config);

      assert(check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                              p.nnear_max2));
      assert(check_nonlocal_dist(sys.pos, box, p.rh2, p.stair2,
                                 p.transient_bonds, p.permanent_bonds));

      // store the configurations in the hdf5 file
      update_config_writer.append();
      update_config_writer.clear();
  }
}

//Config getState_from_pos()
//{
//}

Config run_trajectory_wl(System &sys, Random &mt, const Param &p,
                         const Box &box, UpdateConfig &update_config,
                         CountBond &count_bond,
                         unsigned int iter_wl,
                         bool record_dists,
                         std::vector<double>* dist,
                         DistWriter* dist_writer) {

  LOG_DEBUG("run_trajectory_wl");

  if (sys.useEnsemble and (update_config.config == 0) )
  {
      //std::cout << "  Doing swap MC move. " << std::endl;
      swapMC(sys, mt, box, p);
      UpdateConfig trial_config = config_int(sys.pos, box, p.transient_bonds);
      if (trial_config.config != update_config.config)
      {

        std::ostringstream errorMessage;
        errorMessage << "Possible error with swap move. State is " << trial_config.config
                     << " and was " << update_config.config << std::endl;

          throw std::runtime_error(errorMessage.str());
      }
  }


  for (unsigned int step = iter_wl * p.nsteps_wl; step < (iter_wl + 1) * p.nsteps_wl; step++) {

    EventQueue event_queue;
    Cells cells{p.ncell, p.length / p.ncell};

    //set max time
    if (step != 0) {max_time = (step * p.del_t_wl) + 0.001;}

    initialize_system(sys, mt, p, box, update_config, cells, event_queue);

    const double tot_E_before =
        compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);

    run_step(sys, p, box, update_config, count_bond, cells,
             event_queue, step, p.del_t_wl);


    const double tot_E_during =
        compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);

    const double E_diff = std::abs(1 - (tot_E_during / tot_E_before));

    if (E_diff >= 1e-6) {
      std::cout << E_diff << " energy difference" << std::endl;
      throw std::runtime_error("energy is not conserved");
    }

    if (step % p.write_step == 0 and record_dists) {
      dist_between_nonlocal_beads(sys.pos, box, p.nonlocal_bonds, *dist);
      dist_writer->append(*dist);
    }

  }

  return update_config.config;
}

std::vector<double> createUnevenIntervalPoints(int n, double x_min, double x_max) {
    // Generate a set of n points, starting at x_min and ending at x_max, where they are separated
    // by 1/sqrt(x)
    //
    std::vector<double> points;

    // Step 1: Compute the total integral of 1/sqrt(x) over [x_min, x_max]
    double total_integral = 2.0 * (sqrt(x_max) - sqrt(x_min));

    // Step 2: Compute the boundaries of each segment
    for (int i = 1; i < n; ++i) {
        double ratio = static_cast<double>(i) / n;
        double boundary = pow((ratio * total_integral / 2.0 + sqrt(x_min)), 2);
        points.push_back(boundary);
    }

    // Step 3: Concatenate x_min and x_max with the segment boundaries
    points.insert(points.begin(), x_min);
    points.push_back(x_max);
    //  Now reverse the order
    std::reverse(points.begin(), points.end());

    return points;
}


bool run_anneal(System &sys, Random &mt, const Param &p, const Box &box, double max_d)
{
    bool trapped = true;
    for (unsigned int i = 0; i < p.nbeads; i++){
            sys.times[i] = 0.0;
            sys.counter[i] = 0.0;
    }
    //
    // initialize both members of CountBond struct to 0
    CountBond count_bond = {};
    unsigned int numIter = p.nsteps_max/p.nsteps;
    unsigned int nsteps = p.nsteps;


    // make a local copy of p to modify
    Param p_local = p;
    // set the transient bond to be the target rc for this bond squared
    p_local.transient_bonds.setrc(0, p.rc_target);

    //  Annealing schedule: follow 1/sqrt(x) scheme for interval
    std::vector<double> current_outer_rc = createUnevenIntervalPoints(numIter, p.rc_target + 0.1, max_d + 0.1);

    //double bond_distance_i = compute_transient_dist(sys, p, box);
    //std::cout << " starting ensemble with outer distance = " << current_outer_rc[rc_index] << "  initial d = " << bond_distance_i << std::endl;

    unsigned int rc_index = 0;
    p_local.stair = current_outer_rc[rc_index];
    p_local.stair2 = current_outer_rc[rc_index] * current_outer_rc[rc_index];


    UpdateConfig update_config;
    init_update_config(sys.pos, update_config, box, p_local.transient_bonds);
    if (update_config.config == 1){
        std::cout << " Initial configuration is already in bound state." << std::endl;
	    return false;
    }

    for (unsigned int iter = 0; iter < numIter; iter++){

      double bond_distance = compute_transient_dist(sys, p_local, box);
      // check to verify that the bond distance is not greater than next value of rc in annealing list
      if (bond_distance < current_outer_rc[rc_index+1]){
        // Move outer wall in
        rc_index++;
        if (rc_index >= numIter -1 ) rc_index = numIter-2;

        p_local.stair = current_outer_rc[rc_index];
        p_local.stair2 = current_outer_rc[rc_index] * current_outer_rc[rc_index];
        //std::cout << " Current outer rc is " << current_outer_rc[rc_index] << " d = " << bond_distance << std::endl;
      }


      // run trajectory for nsteps with timestep dt
      run_trajectory_anneal(sys, mt, p_local, box, update_config,
                                         count_bond, iter, nsteps, p.del_t);

	//std::cout << " After iter " << iter << " of " << numIter << " bonds formed is " << count_bond.formed << std::endl;

      auto flips = double(count_bond.formed + count_bond.broken);
      if (flips > 0) {
            std::cout << " Found the target in iteration " << iter << std::endl;
            trapped = false;
            return trapped;
      }
    }
    return trapped;
}


double compute_transient_dist(System &sys, const Param &p, const Box &box) {

    std::tuple<unsigned int, unsigned int, double> t_bond_tuple = p.transient_bonds.getBond(0);

    int bead_i = std::get<0>(t_bond_tuple);
    int bead_j = std::get<1>(t_bond_tuple);

    double dx = sys.pos[bead_i].x - sys.pos[bead_j].x;
    double dy = sys.pos[bead_i].y - sys.pos[bead_j].y;
    double dz = sys.pos[bead_i].z - sys.pos[bead_j].z;
    box.mindist(dx, dy, dz);

    const double dist2 = dx * dx + dy * dy + dz * dz;

    return std::sqrt(dist2);

}

double compute_transient_dist(std::vector<Vec3> &pos, const Param &p, const Box &box) {

    std::tuple<unsigned int, unsigned int, double> t_bond_tuple = p.transient_bonds.getBond(0);

    int bead_i = std::get<0>(t_bond_tuple);
    int bead_j = std::get<1>(t_bond_tuple);

    double dx = pos[bead_i].x - pos[bead_j].x;
    double dy = pos[bead_i].y - pos[bead_j].y;
    double dz = pos[bead_i].z - pos[bead_j].z;
    box.mindist(dx, dy, dz);

    const double dist2 = dx * dx + dy * dy + dz * dz;

    return std::sqrt(dist2);

}


void run_trajectory_anneal(System &sys, Random &mt, Param &p,
                         const Box &box, UpdateConfig &update_config,
                         CountBond &count_bond,
                         unsigned int iter, unsigned int nsteps, double dt) {


    init_update_config(sys.pos, update_config, box, p.transient_bonds);

    for (unsigned int step = iter * nsteps; step < (iter + 1) * nsteps; step++) {

        EventQueue event_queue;
        Cells cells{p.ncell, p.length / p.ncell};

        //set max time
        if (step != 0) {max_time = (step * dt) + 0.001;}

        initialize_system(sys, mt, p, box, update_config, cells, event_queue);

        const double tot_E_before =
            compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);

        run_step(sys, p, box, update_config, count_bond, cells,
                 event_queue, step, dt);


        const double tot_E_during =
            compute_hamiltonian(sys.vel, sys.s_bias, update_config.config, p.m);

        const double E_diff = std::abs(1 - (tot_E_during / tot_E_before));

        if (E_diff >= 1e-6) {
          std::cout << E_diff << " energy difference" << std::endl;
          throw std::runtime_error("energy is not conserved");
    }
  }

}

// Wang-Landau algorithm for estimating entropy
void wang_landau_process(System &sys, Random &mt, const Param &p, const Box &box,
                 UpdateConfig &update_config, CountBond &count_bond,
                 const unsigned int nstates, std::vector<double> &s_bias,
                 DistWriter &dist_writer, std::vector<double> &dist) {

  // amount by which entropy is adjusted
  double gamma = p.gamma;
  // initialize iteration counter at 0
  unsigned int iter_wl = 0;
  unsigned int native_ind = nstates - 1;
  bool record_dists = true;

  while (gamma > p.gamma_f) {
    // iterate over the 2 states
    for (unsigned int i = 0; i < nstates; i++) {
      // run trajectory to get final state
      Config state = run_trajectory_wl(sys, mt, p, box, update_config,
                                       count_bond, iter_wl, record_dists, &dist, &dist_writer);
      if (state == 0) {
        s_bias[native_ind] -= gamma;
      } else {
        s_bias[native_ind] += gamma;
      }

      // update the iteration counter
      iter_wl += 1;
    }
    // normalized gamma updated
    gamma = double(nstates) / double(iter_wl);
  }
}

void from_json(const nlohmann::json &json, Param &p) {
    p.m = json["m"];
    p.sigma = json["sigma_bb"];
    p.sigma2 = p.sigma * p.sigma;
    p.near_min = json["near_min"];
    p.near_max = json["near_max"];
    p.near_min2 = p.near_min * p.near_min;
    p.near_max2 = p.near_max * p.near_max;
    p.nnear_min = json["nnear_min"];
    p.nnear_max = json["nnear_max"];
    p.nnear_min2 = p.nnear_min * p.nnear_min;
    p.nnear_max2 = p.nnear_max * p.nnear_max;
    p.rh = json["rh"];
    p.rh2 = p.rh * p.rh;
    p.rc_target = json["rc_target"];
    p.rc_target2 = p.rc_target * p.rc_target;
    // if stair var is not false aka 0
    if (json.count("stair") != 0) {
        p.stair = json["stair"];
        p.stair2 = *p.stair * *p.stair;
    }

    p.nonlocal_bonds = json["nonlocal_bonds"];
    p.transient_bonds = json["transient_bonds"];
    p.permanent_bonds = json["permanent_bonds"];

    // define the stair bonds
    /*if (json.count("stair_bonds") != 0) {
        p.stair_bonds = json["stair_bonds"];
    }*/

    p.tries = json["tries"];
    p.nbeads = json["nbeads"];
    p.length = json["length"];
    p.ncell = json["ncell"];
    p.nsteps = json["nsteps"];
    p.nsteps_eq = json["nsteps_eq"];
    p.del_t = json["del_t"];
    p.nsteps_wl = json["nsteps_wl"];
    p.del_t_wl = json["del_t_wl"];
    p.gamma = json["gamma"];
    p.gamma_f_screening = json["gamma_f_screening"];
    p.gamma_f = json["gamma_f"];
    p.write_step = json["write_step"];
    p.seeds = json["seeds"].get<std::vector<unsigned int>>();
    p.temp = json["temp"];
    p.mc_moves = json["mc_moves"];
    p.total_iter_initial = json["total_iter_initial"];
    p.total_iter = json["total_iter"];
    p.total_iter_eq = json["total_iter_eq"];
    p.pos_scale = json["pos_scale"];
    p.neg_scale = json["neg_scale"];
    p.sig_level = json["sig_level"];
    p.max_nbonds = json["max_nbonds"];
    p.max_g_test_count = json["max_g_test_count"];
    p.flip_req = json["flip_req"];
    p.fail_max = json["fail_max"];
    p.req_dists = json["req_dists"];
    p.nsteps_max = json["nsteps_max"];
    p.useEnsemble = json["useEnsemble"];
    p.ensembleSize = json["ensembleSize"];
    p.ensemble_write_step = json["ensemble_write_step"];
}
