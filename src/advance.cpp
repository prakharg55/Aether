// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

// -----------------------------------------------------------------------------
// main function to increment model states by one iteration. It needs
// so many inputs because it alters all of the states in the model.
// -----------------------------------------------------------------------------

int advance(Planets &planet,
            Grid &gGrid,
            Times &time,
            Euv &euv,
            Neutrals &neutrals,
            Ions &ions,
            Chemistry &chemistry,
            Electrodynamics &electrodynamics,
            Indices &indices,
            Logfile &logfile) {

  int iErr = 0;

  std::string function = "advance";
  static int iFunction = -1;
  enter(function, iFunction);

  if (time.check_time_gate(get_dt_report()) &&
      test_verbose(0))
    time.display();

  if (get_is_student())
    print(-1, "(1) What function is this " +
		  get_student_name() + "?");

  gGrid.calc_sza(planet, time);
  neutrals.calc_mass_density();
  neutrals.calc_specific_heat();
  time.calc_dt();

  iErr = calc_euv(planet,
                  gGrid,
                  time,
                  euv,
                  neutrals,
                  ions,
                  indices);

  iErr = electrodynamics.update(planet,
                                gGrid,
                                time,
                                ions);
  calc_ion_neutral_coll_freq(neutrals, ions);
  ions.calc_ion_drift(neutrals, gGrid, time.get_dt());

  calc_aurora(gGrid, neutrals, ions);

  // Calculate some neutral source terms:
  neutrals.calc_conduction(gGrid, time);
  chemistry.calc_chemistry(neutrals, ions, time, gGrid);
  if (get_O_cooling())
    neutrals.calc_O_cool();
  if (get_NO_cooling())
    neutrals.calc_NO_cool();
  neutrals.add_sources(time);

  // Calculate Ion and Electron Temperatures:
  ions.calc_ion_temperature(neutrals, gGrid, time);
  ions.calc_electron_temperature(neutrals, gGrid);

  neutrals.set_bcs(gGrid, time, indices);
  neutrals.fill_with_hydrostatic(gGrid);

  neutrals.exchange(gGrid);

  time.increment_time();

  if (time.check_time_gate(get_dt_write_restarts())) {
    print(3, "Writing restart files");
    neutrals.restart_file(get_restartout_dir(), DoWrite);
    ions.restart_file(get_restartout_dir(), DoWrite);
    time.restart_file(get_restartout_dir(), DoWrite);
  }

  iErr = output(neutrals, ions, gGrid, time, planet);

  exit(function);

  logfile.write_logfile(indices, neutrals, ions, gGrid, time);

  return iErr;
}
