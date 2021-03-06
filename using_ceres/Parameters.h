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
#include "enum_functions.h"


struct Parameter{

private :
    static Parameter *s_instance;
public :
    Parameter(){
        input_file_path ="";// "../data/data_in_reduced_export_area/reduced_area.csv";
        output_file_path ="";// "../data/data_in_reduced_export_area/snapping_output.csv" ;
        class_definition_path = "";
        parameters_file_path = "./parameters.txt";
        objects_path = "";
        K_origin =0;// 100;
        K_obs= 0; //1 ;
        K_spacing= 0 ;// 1 ;

        K_original_width = 0;

        K_angle= 0 ;// 1 ;
        K_obj = 0;
        K_obs_width= 0; //1 ;
        K_obj_width= 0; //1 ;
		K_slope = 0 ;
		
        use_manual_initial_position_constraint =false; // false;
        use_manual_distance_to_original_angle = false ;
        use_manual_initial_spacing_constraint = false ;

        use_manual_initial_width_constraint = false ;

        use_manual_distance_to_proj_constraint = false;// false;
        use_manual_Surf_Dist_To_Objects_constraint = false ;
		use_manual_target_slope = false ; 
        use_manual_distance_to_proj_constraint_width = false;// false;
        use_manual_Surf_Dist_To_Objects_constraint_width = false ;
        useLoss = false;//true;
        lossScale = 0;//3.0;
        optimisation_type = SnapEnums::POSITION ;
        optimisation_method = ceres::LINE_SEARCH;

        geom_bound = 5;
        width_bound_minimal = 2 ;
        width_bound_maximal = 20 ;
        width_bound_range = 5 ;
    }

    static Parameter *instance()
    {
        if (!s_instance)
            s_instance = new Parameter;
        return s_instance;
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
    double K_spacing ; //! this parameter scale the measure of similarity between [n_i,n_j] original and current
    double K_angle ; //! this parameter scale the measure of distance between angle at node and original angle

    double K_original_width ; //! this parameter scale the constraint to original width of an edge

    double K_obs ; //! this parameter scale the measure of distance between observation and line (n_i,n_j)
    double K_obj ; /// this param'eter scale the cost of surfacique distance between an object and the edge
    double K_obs_width ; //! this parameter scale the chzange of road width based on sidewalk observation
    double K_obj_width ; /// this param'eter scale the cost of surfacique distance between an object and the edge to change width
	double K_slope; 
	
    bool use_manual_initial_position_constraint;
    bool use_manual_distance_to_original_angle;
    bool use_manual_initial_spacing_constraint;

    bool use_manual_initial_width_constraint;

    bool use_manual_distance_to_proj_constraint;
    bool use_manual_Surf_Dist_To_Objects_constraint;
	bool use_manual_target_slope ; 
    bool use_manual_distance_to_proj_constraint_width;
    bool use_manual_Surf_Dist_To_Objects_constraint_width;

    bool useLoss  ;//! shall we use a loss function to reduce outliers weight
    double lossScale ; //! what shall be the loss function scale (after this scale, outliers mode)

    SnapEnums::optimisation_target optimisation_type;
    ceres::MinimizerType optimisation_method;

    double geom_bound ;
    double width_bound_minimal ;
    double  width_bound_maximal ;
    double width_bound_range ;

};


#endif //PARAMETERS_H

