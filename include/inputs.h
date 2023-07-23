// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_INPUTS_H_
#define INCLUDE_INPUTS_H_

#include <vector>
#include <string>

void Inputs(Times &time);
int read(Times &time);
bool read_inputs_json(Times &time);
int get_verbose();
int get_verbose_proc();
precision_t get_dt_euv();
bool get_include_photoelectrons();
precision_t get_dt_report();
precision_t get_n_outputs();
precision_t get_dt_output(int iOutput);
std::string get_type_output(int iOutput);
precision_t get_euv_heating_eff_neutrals();
std::string get_euv_model();
std::string get_euv_file();
std::string get_aurora_file();
std::string get_chemistry_file();
std::string get_indices_lookup_file();
std::vector<std::string> get_omniweb_files();
int get_number_of_omniweb_files();
std::string get_f107_file();
std::string get_planet();
std::string get_planetary_file();
std::string get_planet_species_file();
std::string get_collision_file();
bool get_do_calc_bulk_ion_temp();
std::string get_bfield_type();
std::string get_electrodynamics_file();
bool get_do_restart();
std::string get_restartout_dir();
std::string get_restartin_dir();
precision_t get_dt_write_restarts();
int get_original_seed();
int get_updated_seed();
void set_seed(int seed);
bool write_restart();
json get_perturb_values(); 
bool get_do_lat_dependent_radius();
bool get_do_J2();

bool get_is_cubesphere();

bool get_NO_cooling();
bool get_O_cooling();

bool get_cent_acc();

std::string get_student_name();
bool get_is_student();

json get_initial_condition_types();
json get_boundary_condition_types();

// ------------------------------
// Grid inputs:

struct grid_input_struct {
  std::string alt_file;
  bool IsUniformAlt;
  precision_t alt_min;
  precision_t dalt;
  precision_t lat_min;
  precision_t lat_max;
  precision_t lon_min;
  precision_t lon_max;
};

grid_input_struct get_grid_inputs();

int get_nLonsGeo();
int get_nLatsGeo();
int get_nAltsGeo();

int get_nBlocksLonGeo();
int get_nBlocksLatGeo();

int get_nMembers();

extern int iVerbose;
extern int iVerboseProc;
extern int iTimingDepth;

std::string get_logfile();
std::vector<std::string> get_species_vector();
bool get_logfile_append();
precision_t get_logfile_dt();

// Satellites
std::vector<std::string> get_satellite_files();
std::vector<std::string> get_satellite_names();
std::vector<precision_t> get_satellite_dts();

std::string get_settings_str(std::string key1);
std::string get_settings_str(std::string key1, std::string key2);
std::vector<int> get_settings_timearr(std::string key1);
std::vector<int> get_settings_intarr(std::string key1);

/**********************************************************************
    \brief Check to see if internal state of class is ok
  **/

bool is_ok();


extern json settings;

extern std::string euv_file;
extern std::string aurora_file;
extern std::string chemistry_file;
extern std::string collision_file;
extern std::string input_file;
extern std::string euv_model;
extern std::string planetary_file;
extern std::vector<std::string> omniweb_files;
extern std::string planet;
extern std::string f107_file;
extern std::string planet_species_file;
extern std::string electrodynamics_file;

extern std::string bfield;

extern grid_input_struct geo_grid_input;

extern precision_t euv_heating_eff_neutrals;
extern precision_t euv_heating_eff_electrons;

extern std::vector<float> dt_output;
extern std::vector<std::string> type_output;
extern std::string output_directory;
extern std::string restart_out_directory;
extern std::string restart_in_directory;

extern bool DoRestart;

extern precision_t dt_euv;
extern precision_t dt_report;

extern int nLonsGeo;
extern int nLatsGeo;
extern int nAltsGeo;

extern int updated_seed;

/// An internal variable to hold the state of the class
extern bool IsOk;

#endif  // INCLUDE_INPUTS_H_
