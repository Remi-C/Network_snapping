//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////



//! reading the parameters in the parameter file
#include "Parameters.h"

Parameter *Parameter::s_instance = 0;

void Parameter::readParameters(){
    //open the parameter file
    //do untill reaching en dof file :
    //read a line
    //parse the line
    //update the parameter if necessary

    //opening files
    FILE* p_fptr = fopen(parameters_file_path.c_str(), "r");

    //std::cout << "reading config file" <<std::endl ;
    char line[1000];//input line buffer
    char comment[1000];//buffer to hold comment line
    char key[1000];//name of the parameter to set
    char value[1000];//value of the parameter to set

    if (p_fptr == NULL) {
        std::cerr << "Error: unable to open parameter file " << parameters_file_path;
        return;
    }

    while ( fgets(line, sizeof line, p_fptr) != NULL) {

        if (	(sscanf(line, "#%s",comment) == 1) ) {
            //std::cerr<< "reading a commentary" << std::endl ;

        }
        else{

            if  (sscanf(line, "%s = %s %s",key, value,comment) >= 2) {
                // std::cout<< "reading a parameter : " << std::endl ;
                // std::cout <<"key :"""<< key << """, value """ << value <<""" "<< std::endl;
                setParameters(key,value);
            }
        }
    }
    fclose(p_fptr);

    return;
}


void Parameter::setParameters(std::string key, std::string value){

    key.compare("input_file_path")==0?input_file_path=value:"NULL";
    key.compare("output_file_path")==0?output_file_path=value:"NULL";
    key.compare("class_definition_path")==0?class_definition_path=value:"NULL";
    key.compare("objects_path")==0?objects_path=value:"NULL";

    key.compare("K_origin")==0?K_origin=atof(value.c_str()):0;
    key.compare("K_spacing")==0?K_spacing=atof(value.c_str()):0;
    key.compare("K_angle")==0?K_angle=atof(value.c_str()):0;

    key.compare("K_original_width")==0?K_original_width=atof(value.c_str()):0;

    key.compare("K_obs")==0?K_obs=atof(value.c_str()):0;
    key.compare("K_obj")==0?K_obj=atof(value.c_str()):0;
    key.compare("K_obs_width")==0?K_obs_width=atof(value.c_str()):0;
    key.compare("K_obj_width")==0?K_obj_width=atof(value.c_str()):0;

    key.compare("use_manual_initial_position_constraint")==0?use_manual_initial_position_constraint=bool(value.compare("false")):false;
    key.compare("use_manual_distance_to_original_angle")==0?use_manual_distance_to_original_angle=bool(value.compare("false")):false;
    key.compare("use_manual_initial_spacing_constraint")==0?use_manual_initial_spacing_constraint=bool(value.compare("false")):false;

    key.compare("use_manual_initial_width_constraint")==0?use_manual_initial_width_constraint=bool(value.compare("false")):false;

    key.compare("use_manual_distance_to_proj_constraint")==0?use_manual_distance_to_proj_constraint=bool(value.compare("false")):false;
    key.compare("use_manual_Surf_Dist_To_Objects_constraint")==0?use_manual_Surf_Dist_To_Objects_constraint=bool(value.compare("false")):false;

    key.compare("use_manual_distance_to_proj_constraint_width")==0?use_manual_distance_to_proj_constraint_width=bool(value.compare("false")):false;
    key.compare("use_manual_Surf_Dist_To_Objects_constraint_width")==0?use_manual_Surf_Dist_To_Objects_constraint_width=bool(value.compare("false")):false;

    key.compare("useLoss")==0?useLoss=bool(value.compare("false")):false;
    key.compare("lossScale")==0?lossScale=atof(value.c_str()):0;

    key.compare("optimisation_type")==0?optimisation_type=SnapEnums::String_toot(value):0;

    key.compare("geom_bound")==0?geom_bound=atof(value.c_str()):0;
    key.compare("width_bound_minimal")==0?width_bound_minimal=atof(value.c_str()):0;
    key.compare("width_bound_maximal")==0?width_bound_maximal=atof(value.c_str()):0;
    key.compare("width_bound_range")==0?width_bound_range=atof(value.c_str()):0;

    return;

}


std::string Parameter::printParameters(){

    std::ostringstream nstring;
    //nstring.precision(10);
    nstring << " input_file_path : " << input_file_path  << std::endl
            << " output_file_path : " << output_file_path  << std::endl
            << " class_definition_path : " << class_definition_path  << std::endl
            << " objects_path : " << objects_path  << std::endl
            << " parameters_file_path : " << parameters_file_path  << std::endl
            << " K_origin : " << K_origin  << std::endl
            << " K_spacing : " << K_spacing  << std::endl
            << " K_angle : " << K_angle  << std::endl
            << " K_original_width : " << K_original_width  << std::endl

            << " K_obs : " << K_obs  << std::endl
            << " K_obj : " << K_obj  << std::endl
            << " K_obs_width : " << K_obs_width  << std::endl
            << " K_obj_width : " << K_obj_width  << std::endl

            << " use_manual_initial_position_constraint : " << use_manual_initial_position_constraint  << std::endl
            << " use_manual_distance_to_original_angle : "
                << use_manual_distance_to_original_angle  << std::endl
            << " use_manual_initial_spacing_constraint : "
                << use_manual_initial_spacing_constraint << std::endl
            << " use_manual_initial_width_constraint : "
                << use_manual_initial_width_constraint << std::endl

            << " use_manual_distance_to_proj_constraint : "
                << use_manual_distance_to_proj_constraint  << std::endl
            << " use_manual_Surf_Dist_To_Objects_constraint : "
                << use_manual_Surf_Dist_To_Objects_constraint_width << std::endl
            << " use_manual_distance_to_proj_constraint_width : "
                << use_manual_distance_to_proj_constraint_width << std::endl
            << " use_manual_Surf_Dist_To_Objects_constraint_width : "
                << use_manual_Surf_Dist_To_Objects_constraint_width << std::endl
            << " useLoss : " << useLoss  << std::endl
            << " lossScale : " << lossScale  << std::endl
            << " optimisation_type : " << SnapEnums::ot_toString(optimisation_type)  << std::endl

            << " geom_bound : " << geom_bound  << std::endl
            << " width_bound_minimal : " << width_bound_minimal  << std::endl
            << " width_bound_maximal : " << width_bound_maximal  << std::endl
            << " width_bound_range : " << width_bound_range  << std::endl  ;

    return nstring.str() ;


}
