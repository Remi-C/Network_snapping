//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////



//! reading the parameters in the parameter file
#include "Parameters.h"


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
    key.compare("K_origin")==0?K_origin=atof(value.c_str()):0;
    key.compare("K_spacing")==0?K_spacing=atof(value.c_str()):0;
    key.compare("K_obs")==0?K_obs=atof(value.c_str()):0;
    key.compare("K_angle")==0?K_angle=atof(value.c_str()):0;
    key.compare("use_initial_position_constraint")==0?use_initial_position_constraint=bool(value.compare("false")):false;
    key.compare("use_initial_spacing_constraint")==0?use_initial_spacing_constraint=bool(value.compare("false")):false;
    key.compare("use_distance_to_proj_constraint")==0?use_distance_to_proj_constraint=bool(value.compare("false")):false;
    key.compare("use_manual_distance_to_proj_constraint")==0?use_manual_distance_to_proj_constraint=bool(value.compare("false")):false;
    key.compare("use_manual_distance_to_original_angle")==0?use_manual_distance_to_original_angle=bool(value.compare("false")):false;
    key.compare("use_manual_initial_spacing_constraint")==0?use_manual_initial_spacing_constraint=bool(value.compare("false")):false;
    key.compare("useLoss")==0?useLoss=bool(value.compare("false")):false;
    key.compare("lossScale")==0?lossScale=atof(value.c_str()):0;
    return;

}


 std::string Parameter::printParameters(){

     std::ostringstream nstring;
     //nstring.precision(10);
     nstring << " input_file_path : " << input_file_path  << std::endl
        << " output_file_path : " << output_file_path  << std::endl
        << " class_definition_path : " << class_definition_path  << std::endl
        << " parameters_file_path : " << parameters_file_path  << std::endl
        << " K_origin : " << K_origin  << std::endl
        << " K_obs : " << K_obs  << std::endl
        << " K_spacing : " << K_spacing  << std::endl
        << " K_angle : " << K_angle  << std::endl
        << " use_initial_position_constraint : " << use_initial_position_constraint  << std::endl
        << " use_initial_spacing_constraint : " << use_initial_spacing_constraint  << std::endl
        << " use_distance_to_proj_constraint : " << use_distance_to_proj_constraint  << std::endl
        << " use_manual_distance_to_proj_constraint : " << use_manual_distance_to_proj_constraint  << std::endl
        << " use_manual_distance_to_original_angle : " << use_manual_distance_to_original_angle  << std::endl
        << " use_manual_initial_spacing_constraint : " << use_manual_initial_spacing_constraint << std::endl
        << " useLoss : " << useLoss  << std::endl
        << " lossScale : " << lossScale  << std::endl  ;
     return nstring.str() ;


 }

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
