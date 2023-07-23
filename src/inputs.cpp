// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "aether.h"

int iVerbose;
int iVerboseProc;
int iTimingDepth;

json settings;

std::string euv_file = "UA/inputs/euv.csv";
std::string aurora_file = "UA/inputs/aurora_earth.csv";
std::string chemistry_file = "UA/inputs/chemistry_earth.csv";
std::string collision_file =
  "UA/inputs/ion_neutral_collision_frequencies.csv";
std::string input_file = "aether.in";
std::string euv_model = "euvac";
std::string planetary_file = "UA/inputs/orbits.csv";
std::vector<std::string> omniweb_files;
std::string planet = "earth";
std::string f107_file = "";
std::string planet_species_file = "";
std::string electrodynamics_file = "";

std::string bfield = "none";

grid_input_struct geo_grid_input;

precision_t euv_heating_eff_neutrals;
precision_t euv_heating_eff_electrons;

std::vector<float> dt_output;
std::vector<std::string> type_output;
std::string output_directory = "UA/output";
std::string restart_out_directory = "UA/restartOut";
std::string restart_in_directory = "UA/restartIn";

bool DoRestart;

precision_t dt_euv;
precision_t dt_report;

int nLonsGeo;
int nLatsGeo;
int nAltsGeo;

int updated_seed;

bool IsOk;

// -----------------------------------------------------------------------
// Initialize the Inputs class.  This also sets some initial values.
// The setting of initial values should probably be moved.
// -----------------------------------------------------------------------

void Inputs(Times &time) {

  // ------------------------------------------------
  // Set some defaults:

  iVerbose = 0;
  iTimingDepth = 3;
  euv_model = "euvac";
  planet = "Earth";

  // ------------------------------------------------
  // Grid Defaults:
  geo_grid_input.alt_file = "";
  geo_grid_input.IsUniformAlt = true;
  geo_grid_input.alt_min = 100.0 * 1000.0;
  geo_grid_input.dalt = 5.0 * 1000.0;

  nLonsGeo = 12;
  nLatsGeo = 20;
  nAltsGeo = 40;

  if (nLonsGeo == 1) {
    geo_grid_input.lon_min = 0.0;
    geo_grid_input.lon_max = 0.0;
  } else {
    geo_grid_input.lon_min = 0.0;
    geo_grid_input.lon_max = 2.0 * cPI;
  }

  if (nLatsGeo == 1) {
    geo_grid_input.lat_min = 0.0;
    geo_grid_input.lat_max = 0.0;
  } else {
    geo_grid_input.lat_min = -cPI / 2;
    geo_grid_input.lat_max = cPI / 2;
  }

  euv_heating_eff_neutrals = 0.40;
  euv_heating_eff_electrons = 0.05;

  dt_output.push_back(300.0);
  type_output.push_back("states");
  dt_euv = 60.0;
  dt_report = 60.0;

  // ------------------------------------------------
  // Now read the input file:
  IsOk = read_inputs_json(time);

  if (!IsOk && iProc == 0)
    std::cout << "Error in reading input file!\n";
}

// -----------------------------------------------------------------------
// output the settings json to a file (for restart)
// -----------------------------------------------------------------------

bool write_restart() {
  bool DidWork = true;

  if (iProc == 0) {
    std::string filename = settings["Restart"]["OutDir"];
    filename = filename + "/settings.json";
    DidWork = write_json(filename, settings);
  }

  DidWork = sync_across_all_procs(DidWork);
  return DidWork;
}

// -----------------------------------------------------------------------
// Return log file name
// -----------------------------------------------------------------------

std::string get_logfile() {
  return settings["Logfile"]["name"];
}

// -----------------------------------------------------------------------
// Return how oftern to write log file
// -----------------------------------------------------------------------

precision_t get_logfile_dt() {
  return settings["Logfile"]["dt"];
}

// -----------------------------------------------------------------------
// Return whether to append or rewrite
// -----------------------------------------------------------------------

bool get_logfile_append() {
  return settings["Logfile"]["append"];
}

// -----------------------------------------------------------------------
// Return the name of specified variables as a vector
// -----------------------------------------------------------------------

std::vector<std::string> get_species_vector() {
  std::vector<std::string> species;
  const json &json_species = settings["Logfile"]["species"];

  for (size_t iOutput = 0; iOutput < json_species.size(); iOutput++) {
    species.push_back(json_species.at(iOutput));
  }

  return species;
}

// -----------------------------------------------------------------------
// Return the name of satellite files as a vector
// -----------------------------------------------------------------------

std::vector<std::string> get_satellite_files() {
    std::vector<std::string> files;
    const json &json_files = settings["Satellites"]["files"];

    for (size_t i = 0; i < json_files.size(); ++i) {
        files.push_back(json_files.at(i));
    }

    return files;
}

// -----------------------------------------------------------------------
// Return the output file names of satellites as a vector
// -----------------------------------------------------------------------

std::vector<std::string> get_satellite_names() {
    std::vector<std::string> names;
    const json &json_names = settings["Satellites"]["names"];

    for (size_t i = 0; i < json_names.size(); ++i) {
        names.push_back(json_names.at(i));
    }

    return names;
}

// -----------------------------------------------------------------------
// Return how oftern to write log file for satellites as a vector
// -----------------------------------------------------------------------

std::vector<precision_t> get_satellite_dts() {
    std::vector<precision_t> dts;
    const json &json_dts = settings["Satellites"]["dts"];

    for (size_t i = 0; i < json_dts.size(); ++i) {
        dts.push_back(json_dts.at(i));
    }

    return dts;
}


// -----------------------------------------------------------------------
// Return value of a key in the json formatted inputs
// -----------------------------------------------------------------------

std::string get_settings_str(std::string key1) {
  std::string value = "unknown";

  if (settings.contains(key1))
    value = settings.at(key1);

  return value;
}

std::string get_settings_str(std::string key1,
                                     std::string key2) {
  std::string value = "unknown";

  if (settings.contains(key1))
    if (settings.at(key1).contains(key2))
      value = settings.at(key1).at(key2);

  return value;
}

std::vector<int> get_settings_intarr(std::string key1) {
  std::vector<int> value;

  if (settings.contains(key1)) {
    int nPts = settings.at(key1).size();

    for (int i = 0; i < nPts; i++)
      value.push_back(settings.at(key1).at(i));
  }

  return value;
}

std::vector<int> get_settings_timearr(std::string key1) {
  int nPtsTime = 7;
  std::vector<int> outarr(nPtsTime, 0);
  std::vector<int> timearr = get_settings_intarr(key1);

  int nPts = timearr.size();

  if (nPts > nPtsTime)
    nPts = nPtsTime;

  for (int i = 0; i < nPts; i++)
    outarr[i] = timearr[i];

  return outarr;
}

// -----------------------------------------------------------------------
// Return characteristics of the grid that are entered by the user
// -----------------------------------------------------------------------

grid_input_struct get_grid_inputs() {
  // First Get Values:
  geo_grid_input.alt_file = settings["GeoGrid"]["AltFile"];
  geo_grid_input.IsUniformAlt = settings["GeoGrid"]["IsUniformAlt"];
  geo_grid_input.alt_min = settings["GeoGrid"]["MinAlt"];
  geo_grid_input.dalt = settings["GeoGrid"]["dAlt"];
  geo_grid_input.lat_min = settings["GeoGrid"]["MinLat"];
  geo_grid_input.lat_max = settings["GeoGrid"]["MaxLat"];
  geo_grid_input.lon_min = settings["GeoGrid"]["MinLon"];
  geo_grid_input.lon_max = settings["GeoGrid"]["MaxLon"];

  // Second Change Units
  geo_grid_input.alt_min = geo_grid_input.alt_min * cKMtoM;
  geo_grid_input.lat_min = geo_grid_input.lat_min * cDtoR;
  geo_grid_input.lat_max = geo_grid_input.lat_max * cDtoR;
  geo_grid_input.lon_min = geo_grid_input.lon_min * cDtoR;
  geo_grid_input.lon_max = geo_grid_input.lon_max * cDtoR;

  // If the grid is uniform, dalt is in km, else it is in fractions of
  // scale height:
  if (geo_grid_input.IsUniformAlt)
    geo_grid_input.dalt = geo_grid_input.dalt * cKMtoM;

  return geo_grid_input;
}

// -----------------------------------------------------------------------
// Return whether user is student
// -----------------------------------------------------------------------

bool get_is_student() {
  return settings["Student"]["is"];
}

// -----------------------------------------------------------------------
// Return student name
// -----------------------------------------------------------------------

std::string get_student_name() {
  return settings["Student"]["name"];
}

// -----------------------------------------------------------------------
// Return whether grid is cubesphere or spherical
// -----------------------------------------------------------------------

bool get_is_cubesphere() {
  return settings["CubeSphere"]["is"];
}

// -----------------------------------------------------------------------
// Return whether to restart or not
// -----------------------------------------------------------------------

bool get_do_restart() {
  return settings["Restart"]["do"];
}

// -----------------------------------------------------------------------
// Return NO cooling
// -----------------------------------------------------------------------

bool get_NO_cooling() {
  return settings["Sources"]["Neutrals"]["NO_cool"];
}

// -----------------------------------------------------------------------
// Return O cooling
// -----------------------------------------------------------------------

bool get_O_cooling() {
  return settings["Sources"]["Neutrals"]["O_cool"];
}

// -----------------------------------------------------------------------
// Return centripetal acceleration
// -----------------------------------------------------------------------

bool get_cent_acc() {
  return settings["Sources"]["Grid"]["Cent_acc"];
}

// -----------------------------------------------------------------------
// Return restart OUT directory
// -----------------------------------------------------------------------

std::string get_restartout_dir() {
  return settings["Restart"]["OutDir"];
}

// -----------------------------------------------------------------------
// Return restart OUT directory
// -----------------------------------------------------------------------

std::string get_restartin_dir() {
  return settings["Restart"]["InDir"];
}

// -----------------------------------------------------------------------
// dt for writing restart files
// -----------------------------------------------------------------------

precision_t get_dt_write_restarts() {
  return settings["Restart"]["dt"];
}

// -----------------------------------------------------------------------
// Return magnetic field type (dipole and none defined now.)
// -----------------------------------------------------------------------

std::string get_bfield_type() {
  return settings["BField"];
}

// -----------------------------------------------------------------------
// Return the EUV model used (EUVAC only option now)
// -----------------------------------------------------------------------

std::string get_euv_model() {
  return settings["Euv"]["Model"];
}

// -----------------------------------------------------------------------
// Return the heating efficiency of the neutrals for EUV
// -----------------------------------------------------------------------

precision_t get_euv_heating_eff_neutrals() {
  return settings["Euv"]["HeatingEfficiency"];
}

// -----------------------------------------------------------------------
// Return whether to include the photoelectron ionization
// -----------------------------------------------------------------------

bool get_include_photoelectrons() {
  return settings["Euv"]["IncludePhotoElectrons"];
}

// -----------------------------------------------------------------------
// Return how often to calculate EUV energy deposition
// -----------------------------------------------------------------------

precision_t get_dt_euv() {
  return settings["Euv"]["dt"];
}

// -----------------------------------------------------------------------
// Return how often to report progress of simulation
// -----------------------------------------------------------------------

precision_t get_dt_report() {
  return settings["Debug"]["dt"];
}

// -----------------------------------------------------------------------
// Return number of output types
// -----------------------------------------------------------------------

precision_t get_n_outputs() {
  return settings["Outputs"]["type"].size();
}

// -----------------------------------------------------------------------
// Return original random number seed
// -----------------------------------------------------------------------

int get_original_seed() {
  return settings["Seed"];
}

// -----------------------------------------------------------------------
// Set random number seed
// -----------------------------------------------------------------------

void set_seed(int seed) {
  settings["Seed"] = seed;
  updated_seed = seed;
}

// -----------------------------------------------------------------------
// Return random number seed that has been updated
// -----------------------------------------------------------------------

int get_updated_seed() {
  std::default_random_engine get_random(updated_seed);
  updated_seed = get_random();
  return updated_seed;
}

// -----------------------------------------------------------------------
// Return number of longitudes, latitudes, and altitudes in Geo grid
// -----------------------------------------------------------------------

int get_nLonsGeo() {
  return settings["GeoBlockSize"]["nLons"];
}

int get_nLatsGeo() {
  return settings["GeoBlockSize"]["nLats"];;
}

int get_nAltsGeo() {
  return settings["GeoBlockSize"]["nAlts"];
}

// -----------------------------------------------------------------------
// Return number of Blocks of longitudes and latitudes in Geo grid
// -----------------------------------------------------------------------

int get_nBlocksLonGeo() {
  return settings["GeoBlockSize"]["nBlocksLon"];
}

int get_nBlocksLatGeo() {
  return settings["GeoBlockSize"]["nBlocksLat"];;
}

// -----------------------------------------------------------------------
// Return number of ensemble members
// -----------------------------------------------------------------------

int get_nMembers() {
  return settings["Ensembles"]["nMembers"];
}

// -----------------------------------------------------------------------
// Return verbose variables
// -----------------------------------------------------------------------

int get_verbose() {
  return settings["Debug"]["iVerbose"];
}

int get_verbose_proc() {
  return settings["Debug"]["iProc"];
}

// -----------------------------------------------------------------------
// Return how often to output a given output type
// -----------------------------------------------------------------------

precision_t get_dt_output(int iOutput) {
  precision_t value = 0.0;
  int nOutputs = settings.at("Outputs").at("type").size();

  if (iOutput < nOutputs)
    value = settings.at("Outputs").at("dt").at(iOutput);

  return value;
}

// -----------------------------------------------------------------------
// Return the output type
// -----------------------------------------------------------------------

std::string get_type_output(int iOutput) {
  std::string value = "";
  int nOutputs = settings.at("Outputs").at("type").size();

  if (iOutput < nOutputs)
    value = settings.at("Outputs").at("type").at(iOutput);

  return value;
}

// -----------------------------------------------------------------------
// Return EUV file name
// -----------------------------------------------------------------------

std::string get_euv_file() {
  return settings["Euv"]["File"];
}

// -----------------------------------------------------------------------
// Return aurora file name
// -----------------------------------------------------------------------

std::string get_aurora_file() {
  return settings["AuroraFile"];
}

// -----------------------------------------------------------------------
// Return Chemistry file name
// -----------------------------------------------------------------------

std::string get_chemistry_file() {
  return settings["ChemistryFile"];
}

// -----------------------------------------------------------------------
// Return Collision file name
// -----------------------------------------------------------------------

std::string get_collision_file() {
  return settings["CollisionsFile"];
}

// -----------------------------------------------------------------------
// Return Indices Lookup Filename
// -----------------------------------------------------------------------

std::string get_indices_lookup_file() {
  return settings["IndicesLookupFile"];
}

// -----------------------------------------------------------------------
// Return total number of OMNIWeb files to read
// -----------------------------------------------------------------------

int get_number_of_omniweb_files() {
  return settings["OmniwebFiles"].size();
}

// -----------------------------------------------------------------------
// Return OMNIWeb file names as a vector
// -----------------------------------------------------------------------

std::vector<std::string> get_omniweb_files() {
  std::vector<std::string> omniweb_files;
  int nFiles = settings["OmniwebFiles"].size();

  for (int i = 0; i < nFiles; i++)
    omniweb_files.push_back(settings.at("OmniwebFiles").at(i));

  return omniweb_files;
}

// -----------------------------------------------------------------------
// Return F107 file to read
// -----------------------------------------------------------------------

std::string get_f107_file() {
  return settings["F107File"];
}

// -----------------------------------------------------------------------
// Return planet name
// -----------------------------------------------------------------------

std::string get_planet() {
  return settings["Planet"];
}

// -----------------------------------------------------------------------
// Return file that contains (all) planetary characteristics
// -----------------------------------------------------------------------

std::string get_planetary_file() {
  return settings["PlanetCharacteristicsFile"];
}

// -----------------------------------------------------------------------
// Return planetary file name that describes the species and such for
// a given planet
// -----------------------------------------------------------------------

std::string get_planet_species_file() {
  return settings["PlanetSpeciesFile"];
}

std::string get_electrodynamics_file() {
  return settings["ElectrodynamicsFile"];
}

// -----------------------------------------------------------------------
// Flag to do the bulk ion temperature calculation instead
// of individual ion specie temperature calculations
// -----------------------------------------------------------------------

bool get_do_calc_bulk_ion_temp() {
  return settings["DoCalcBulkIonTemp"];
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

json get_perturb_values() {
  json values;

  if (settings.contains("Perturb"))
    values = settings["Perturb"];

  return values;
}

// -----------------------------------------------------------------------
// Flag to have a latitude dependent radius, and by extension gravity
// -----------------------------------------------------------------------

bool get_do_lat_dependent_radius() {
  return settings["Oblate"]["isOblate"];
}

// -----------------------------------------------------------------------
// Flag to include J2 term in the gravity calculation
// -----------------------------------------------------------------------

bool get_do_J2() {
  return settings["Oblate"]["isJ2"];
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

json get_initial_condition_types() {
  json values;

  if (settings.contains("InitialConditions"))
    values = settings["InitialConditions"];

  return values;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

json get_boundary_condition_types() {
  json values;

  if (settings.contains("BoundaryConditions"))
    values = settings["BoundaryConditions"];

  return values;
}

// --------------------------------------------------------------------------
// check to see if class is ok
// --------------------------------------------------------------------------

bool is_ok() {
  return IsOk;
}
