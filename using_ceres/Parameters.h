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
#include <cstdlib>

#include "ceres/ceres.h"

#include "Data.h"

struct Parameter{
   public :
    Parameter(){
        input_file_path ="";// "../data/data_in_reduced_export_area/reduced_area.csv";
        output_file_path ="";// "../data/data_in_reduced_export_area/snapping_output.csv" ;
        class_definition_path = "";
        parameters_file_path = "../parameters.txt";
        objects_path = "";
        K_origin =0;// 100;
        K_obs= 0; //1 ;
        K_spacing= 0 ;// 1 ;
        K_angle= 0 ;// 1 ;
        K_obj = 1;
        use_initial_position_constraint =true; // false;
        use_initial_spacing_constraint = true; //false;
        use_distance_to_proj_constraint = true; //false;
        use_manual_distance_to_proj_constraint = true;// false;
        use_manual_distance_to_original_angle = true ;
        use_manual_initial_spacing_constraint = true ;
        use_manual_Surf_Dist_To_Objects_constraint = true ;

         useLoss = false;//true;
         lossScale = 0;//3.0;
    }

    void readParameters();
    void setParameters(std::string , std::string );
    std::string printParameters();


    std::string input_file_path ;
    std::string output_file_path;
    std::string class_definition_path;
    std::string parameters_file_path;
    std::string objects_path;

    double K_origin ; //! this parameter scale the distance to origin for a node
    double K_obs ; //! this parameter scale the measure of distance between observation and line (n_i,n_j)
    double K_spacing ; //! this parameter scale the measure of similarity between [n_i,n_j] original and current
    double K_angle ; //! this parameter scale the measure of distance between angle at node and original angle
    double K_obj ; /// this param'eter scale the cost of surfacique distance between an object and the edge

    bool use_initial_position_constraint;
    bool use_initial_spacing_constraint;
    bool use_distance_to_proj_constraint;
    bool use_manual_distance_to_proj_constraint;
    bool use_manual_distance_to_original_angle;
    bool use_manual_initial_spacing_constraint;
    bool use_manual_Surf_Dist_To_Objects_constraint;

    bool useLoss  ;//! shall we use a loss function to reduce outliers weight
    double lossScale ; //! what shall be the loss function scale (after this scale, outliers mode)
};


#endif //PARAMETERS_H

