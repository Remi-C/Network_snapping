#ifndef PARAMETERS_H
#define PARAMETERS_H
//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////
/**
  This file holds the global paramters.
  When mature code is reached, should put every parameters in a parameter class.

  */

#include <string>

#include "ceres/ceres.h"

#include "Data.h"

struct Parameter{
   public :
    Parameter(){
        input_file_path = "../data/data_in_reduced_export_area/reduced_area.csv";
        output_file_path =  "../data/data_in_reduced_export_area/snapping_output.csv" ;
        K_origin = 100;
        K_obs= 1 ;
        K_spacing=  1 ;
        use_initial_position_constraint = false;
        use_initial_spacing_constraint = true;
        use_distance_to_proj_constraint = false;
        use_manual_distance_to_proj_constraint = true;

         useLoss = true;
         lossScale = 3.0;
    }

    std::string input_file_path ;
    std::string  output_file_path;

    double K_origin ; //! this parameter scale the distance to origin for a node
    double K_obs ; //! this parameter scale the measure of distance between observation and line (n_i,n_j)
    double K_spacing ; //! this parameter scale the measure of similarity between [n_i,n_j] original and current


    bool use_initial_position_constraint;
    bool use_initial_spacing_constraint;
    bool use_distance_to_proj_constraint;
    bool use_manual_distance_to_proj_constraint;

    bool useLoss  ;//! shall we use a loss function to reduce outliers weight
    double lossScale ; //! what shall be the loss function scale (after this scale, outliers mode)
};


#endif //PARAMETERS_H

