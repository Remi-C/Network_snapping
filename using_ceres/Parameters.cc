//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////
/**
  This file holds the global paramters.
  When mature code is reached, should put every parameters in a parameter class.

  */

/*
//#include "Parameters.h"

#include <string>
#include <string.h>

#include "ceres/ceres.h"

#include "Data.h"


//! @TODO @WARNING very ugly : should be put in a "parameter class"
///////////////////////PARAMETERS///////////////////////
const std::string input_file_path("../data/data_in_reduced_export_area/reduced_area.csv");
const std::string  output_file_path("../data/data_in_reduced_export_area/snapping_output.csv");

const double K_origin = 1  ; //! this parameter scale the distance to origin for a node
const double K_obs= 1 ; //! this parameter scale the measure of distance between observation and line (n_i,n_j)
const double K_spacing=  1 ; //! this parameter scale the measure of similarity between [n_i,n_j] original and current


const bool use_initial_position_constraint = false;
const bool use_initial_spacing_constraint = false;
const bool use_distance_to_proj_constraint = true;

const bool useLoss = true;//! shall we use a loss function to reduce outliers weight
const double lossScale = 3.0; //! what shall be the loss function scale (after this scale, outliers mode)
////////////////////////////////////////////////////////


//////////////// gloable variable. Should use a singleton instead
DataStorage * data_pointer ;
int addConstraintsOnInitialPosition(DataStorage *, ceres::Problem * );
int addConstraintsOnInitialspacing( DataStorage *, ceres::Problem * );
int addConstraintsOnOrthDistToObservation(DataStorage * , ceres::Problem * ) ;
////////////////
*/
