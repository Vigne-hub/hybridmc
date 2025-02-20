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
//#include "molecule.h"
namespace po = boost::program_options;

//#define LOCAL_DEBUG

void mean_stdev(const std::vector<double> v,double& mean, double& stdev)
{
    int nterms = v.size();
    if (nterms == 0)
    {
        mean = 0.0;
        stdev = 0.0;
    }
    else if (nterms == 1)
    {
        mean = v[0];
        stdev = 0.0;
    }
    else
    {
        double sum = std::accumulate(v.begin(), v.end(), 0.0);
        mean = sum / nterms;
        std::vector<double> diff(v.size());
        std::transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return x - mean; });
        double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        stdev = std::sqrt(sq_sum / (nterms - 1.0) );
    }
}

int main(int argc, char *argv[]) {
  // restart program from point of interruption
  // see Prus, 2002, Boost.Program_options doc
  //
  // declare the supported options
  po::options_description desc("allowed options");
  desc.add_options()("help,h", "produce help message")(
      "json-file", po::value<std::string>()->required(), "json")(
      "output-file", po::value<std::string>()->required(),
      "hdf5 output")("input-file", po::value<std::string>(), "hdf5 input")(
      "snapshot-file", po::value<std::string>(), "hdf5 snapshot");

  po::positional_options_description pod;
  pod.add("json-file", 1);
  pod.add("output-file", 1);

  po::variables_map vm;
  try {
    po::store(
        po::command_line_parser(argc, argv).options(desc).positional(pod).run(),
        vm);

    if (vm.count("help")) {
      std::cout << "usage: hybridmc [options]... json-file output-file"
                << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 0;
    }

    po::notify(vm);
  } catch (const po::error &e) {
    std::cerr << "hybridmc: " << e.what() << std::endl;
    std::cerr << "try hybridmc --help" << std::endl;
    return 1;
  }

  std::optional<std::string> input_name;
  if (vm.count("input-file")) {
    input_name = vm["input-file"].as<std::string>();
  }

  std::optional<std::string> snapshot_name;
  if (vm.count("snapshot-file")) {
    snapshot_name = vm["snapshot-file"].as<std::string>();
  }

  const std::string json_name = vm["json-file"].as<std::string>();
  const std::string output_name = vm["output-file"].as<std::string>();

  std::cout << "git commit " << VERSION << std::endl;

  std::ifstream input(json_name);
  nlohmann::json json;
  input >> json;

  Param p = json;
  Param p_eq = p;
  const Box box{p.length};
  UpdateConfig update_config;

  const unsigned int t_bonds = p.transient_bonds.get_nbonds();
  const unsigned int nbonds = p.nonlocal_bonds.get_nbonds();
  const unsigned int nstates = std::pow(2, t_bonds);
  ConfigInt store_config_int;
  std::vector<uint64_t> config_count(nstates);
  std::vector<double> dist(nbonds);

  // initialize both members of CountBond struct to 0
  CountBond count_bond = {};

  System sys(p.nbeads);
  sys.distanceWrite = true;

  if (p.useEnsemble)
  {
    std::cout << " Will use ensemble of size " << p.ensembleSize << std::endl;
    sys.useEnsemble = true;
    sys.ensemble.reserve(p.ensembleSize);
  }
  sys.ensembleSize = p.ensembleSize;
  sys.nextEnsemble.reserve(p.ensembleSize);

  std::seed_seq seq(p.seeds.begin(), p.seeds.end());
  Random mt(seq);

  // open output file
  H5::H5File file(output_name + ".tmp", H5F_ACC_TRUNC);

  UpdateConfigWriter update_config_writer(file);

  H5::Group group{file.createGroup("unique_config")};
  H5::Group group_ensemble{file.createGroup("ensemble")};

  PosWriter pos_writer{group, "pos", p.nbeads};
  VelWriter vel_writer{group, "vel", p.nbeads};
  ConfigWriter config_writer{group, "config"};
  ConfigCountWriter config_count_writer{file, "config_count", nstates};

  // Initialize distance writer with the H5 file, given nbonds and dataset name "dist"
  DistWriter dist_writer{file, "dist", nbonds};

  config_count_writer.append(config_count);
  std::set<Config> store_config;

  H5::Attribute attr_sigma =
      file.createAttribute("sigma", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
  attr_sigma.write(H5::PredType::NATIVE_DOUBLE, &p.sigma);
  H5::Attribute attr_length =
      file.createAttribute("length", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
  attr_length.write(H5::PredType::NATIVE_DOUBLE, &box.l);
  H5::Attribute attr_nsteps =
      file.createAttribute("nsteps", H5::PredType::NATIVE_UINT, H5S_SCALAR);
  attr_nsteps.write(H5::PredType::NATIVE_UINT, &p.nsteps);

  p.nonlocal_bonds.write_hdf5(file, "nonlocal_bonds");
  p.transient_bonds.write_hdf5(file, "transient_bonds");
  p.permanent_bonds.write_hdf5(file, "permanent_bonds");

  // equilibrate initial structure
  //
  // reset bead clocks, counters, and wall time
  for (unsigned int i = 0; i < p.nbeads; i++) {
    sys.times[i] = 0.0;
    sys.counter[i] = 0.0;
  }


  initialize_pos(sys, mt, p_eq, box, update_config, input_name,
                 snapshot_name, t_bonds);

  init_update_config(sys.pos, update_config, box, p.transient_bonds);

  sys.distanceWrite = false; // don't write distance data during equilibration step since states aren't equally likely
  sys.recordEnsemble = false; // don't record configurations for next ensemble yet

  for (unsigned int iter = 0; iter < p_eq.total_iter_eq; iter++) {
     run_trajectory_eq(sys, mt, p, box, update_config, count_bond,
                      iter, dist_writer, dist);
  }

  init_update_config(sys.pos, update_config, box, p.transient_bonds);

  for (unsigned int i = 0; i < p.nbeads; i++) {
    sys.times[i] = 0.0;
    sys.counter[i] = 0.0;
  }
  // Wang-Landau
  wang_landau_process(sys, mt, p, box, update_config, count_bond, nstates, sys.s_bias, dist_writer, dist);
  double dS_start = sys.s_bias[0] - sys.s_bias[1];
  if (dS_start > 4.5)
  {
    std::cout << "  Initial dS = " << dS_start << " is suspiciously large after WL routine.  Running WL longer." << std::endl;
    init_update_config(sys.pos, update_config, box, p.transient_bonds);
    for (unsigned int i = 0; i < p.nbeads; i++)
    {
        sys.times[i] = 0.0;
        sys.counter[i] = 0.0;
    }
    wang_landau_process(sys, mt, p, box, update_config, count_bond, nstates, sys.s_bias, dist_writer, dist);
    std::cout << " After repeat, dS is now " << sys.s_bias[0] - sys.s_bias[1] << std::endl;
  }
  //
  // Initialize output binary files
  //
  const hsize_t mem_dims[1] = {sys.s_bias.size()};
  H5::DataSpace mem_space(1, mem_dims);
  H5::DataSet dataset_s_bias =
      file.createDataSet("s_bias", H5::PredType::NATIVE_FLOAT, mem_space);

  // parameters for checking if g test is validated or not
  unsigned int g_test_count = 0;
  double flipping_rate = 0.0;
  bool done_flip = false;
  bool done_g_test = false;
  bool done_distances = false;
  bool done_ensemble = false;
  sys.recordEnsemble = true; // turn on recording of next ensemble of configurations

  int fail_counter = 0;
  int current_distances = dist_writer.get_size();
  std::cout << " Starting convergence loop with Sbias values = "
        << sys.s_bias[0] << " " << sys.s_bias[1] << " recorded distances = "
        << current_distances << std::endl;
  //
  //  Carry out quick convergence with small number of total counts to get approximately uniform
  // probability of states to run for accurate results with large number of counts.
  //

  std::vector<double> Sbias_array; // store all converged values according to g-test
  while (!done_g_test or !done_flip or !done_distances or !done_ensemble){
        // reset bead clocks, counters and wall time
        for (unsigned int i = 0; i < p.nbeads; i++) {
          sys.times[i] = 0.0;
          sys.counter[i] = 0.0;
        }

        // reset configuration count
        store_config_int.clear();
        std::fill(config_count.begin(), config_count.end(), 0);

        // add this to see to only count bond events in each iteration
        count_bond.formed = 0;
        count_bond.broken = 0;
        sys.distanceWrite = false; // default to no distance data storage
        current_distances = dist_writer.get_size();
        if (current_distances < p.req_dists)
        {
            sys.distanceWrite = true; // only record distance data if needed
            std::cout << " Recording distances since size  = " << current_distances
                        << "  req_dists = " << p.req_dists << std::endl;
        }

        //
        //  Generate p.total_iter samples of state
        //
        for (unsigned int iter = 0; iter < p.total_iter_initial; iter++)
        {
            run_trajectory(sys, mt, p, box, dist, update_config, update_config_writer,
                         pos_writer, vel_writer, config_writer, dist_writer,
                         store_config, store_config_int, count_bond,
                         iter);
        }

        // count the number of times each configuration is visited
        for (unsigned long i: store_config_int) {
          config_count[i]++;
        }

        // Now recompute the s_bias
        compute_entropy(config_count, sys.s_bias, nstates, p.pos_scale, p.neg_scale);


        auto flips = double(count_bond.formed + count_bond.broken);
        auto  stateCount = std::reduce(config_count.begin(), config_count.end());
        flipping_rate = flips / stateCount;

        done_g_test = g_test(config_count,nstates, p.sig_level);
        if (done_g_test)  Sbias_array.push_back(sys.s_bias[0]);

        if (flipping_rate > p.flip_req) {
          done_flip = true;
        } else {
          done_flip = false;
          fail_counter++;
          // check if too many fails have occurred
          if (fail_counter >= p.fail_max) {
            //std::cout << " Too many fails have occurred." << std::endl;
            // increase the number of steps in trajectory
            if (flipping_rate > 0)
            {
                p.nsteps = int(p.flip_req * p.nsteps/flipping_rate) + 1;
                if (p.nsteps > p.nsteps_max) p.nsteps = p.nsteps_max;
            }
            else
            {
                // flip rate is zero
                throw std::runtime_error("Flip rate was zero with maximum trajectory! Increase del_t to get longer trajectories.");
            }

            std::cout << " Will increase the number of steps in trajectory to " << p.nsteps << " time of propagation ="
                    << p.nsteps*p.del_t << " for next iteration." << std::endl;
            // reset fail counter
            fail_counter = 0;
          } // end of fail_counter MD adjustment
        } // end of flipping rate check

        int numD = dist_writer.get_size();
        if (numD < p.req_dists)
        {
            std::cout << "  Too few distances recorded for convergence: number is now " << numD << std::endl;
        }
        else
        {
            done_distances = true;
            sys.distanceWrite = false;
        }

        int numRecordedEnsemble = sys.nextEnsemble.size();
        if (numRecordedEnsemble < p.ensembleSize)
        {
            done_ensemble = false;
           // std::cout << " Have recorded " << numRecordedEnsemble << " configurations of next ensemble size "
            //    << p.ensembleSize << std::endl;
            sys.recordEnsemble = true;
        }
        else
        {
            done_ensemble = true;
            sys.recordEnsemble = false;
        }

        //std::cout << "Working on output_file: " << output_name << std::endl;
        std::cout << " In iteration " << g_test_count << " flip rate = " << flipping_rate
                  << " done_flip = " << done_flip
                  << " done_g = " << done_g_test << " done_distances = " << done_distances
                  << " done_ensemble = " << done_ensemble
                  << std::endl << "      ensemble size = " << numRecordedEnsemble
                  << "  distances recorded = " << numD << std::endl << std::endl;
        g_test_count++;
        if (g_test_count >= p.max_g_test_count)
        {
            std::cout << "  Done the maximum number of iterations: " << g_test_count
                    <<  " without convergence." << std::endl;
            break;
        }
  } // end of while loop

  double mean_s = 0.0;
  double stdev_s = 0.0;
  mean_stdev(Sbias_array, mean_s, stdev_s);

  std::cout << " Average s_bias is " << mean_s << "  standard deviation = " << stdev_s << " for "
            << Sbias_array.size() << " values." << std::endl;

  sys.s_bias[0] = mean_s;
  done_g_test = false;

  std::cout << std::endl << " Starting refined sampling with " << p.total_iter << " states to get good statistics with Sbias[0] = "
        << sys.s_bias[0] << " time steps per trajectory = " << p.nsteps << std::endl;
  sys.distanceWrite = false;
  sys.recordEnsemble = false;

  while (!done_g_test){
      //  Store value of s_bias used to generate dataset
      dataset_s_bias.write(&sys.s_bias[0], H5::PredType::NATIVE_DOUBLE, mem_space);

      //  Now sample a long trajectory of states
      for (unsigned int i = 0; i < p.nbeads; i++)
      {
              sys.times[i] = 0.0;
              sys.counter[i] = 0.0;
      }
      // reset configuration count
      store_config_int.clear();
      std::fill(config_count.begin(), config_count.end(), 0);
      //


      for (unsigned int iter = 0; iter < p.total_iter; iter++)
      {
          run_trajectory(sys, mt, p, box, dist, update_config, update_config_writer,
                         pos_writer, vel_writer, config_writer, dist_writer,
                         store_config, store_config_int, count_bond,
                         iter, true);
      }

      // count the number of times each configuration is visited
      for (unsigned long i: store_config_int) {
          config_count[i]++;
      }
      config_count_writer.append(config_count); // store trajectory of sampled states in binary .h5 file
      //
      //std::cout << std::endl << "  Next ensemble is:" << std::endl;
      //printEnsemble(sys.nextEnsemble);
      //

      //
      // Output data for analysis:  Note that sys.s_bias recorded in .h5 file is the base sbias used to generate data
      //  The diff_s_bias.py routine will later recompute the difference in entropy based on the trajectory
      //  and bootstrap the errors.
      //
      // final unbiased entropy estimate from long run with pos_scale = 1.0, neg_scale = 1.0
      //  This value is not recorded in h5 file, but it computed with diff_s_bias.py with bootstraps
      //

      done_g_test = g_test(config_count,nstates, p.sig_level);
       if (!done_g_test)
       {
            // minor adjustment to s_bias for next iteration
            compute_entropy(config_count, sys.s_bias, nstates, p.pos_scale, p.neg_scale);
       }

  }
  compute_entropy(config_count, sys.s_bias, nstates, 1.0, 1.0);  // this will be true estimate of s_bias
  write_ensemble(group_ensemble, sys);
  file.close();
  std::filesystem::rename(output_name + ".tmp", output_name);

  return 0;
}

