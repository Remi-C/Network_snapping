//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////
/**
  @TODO
    _ use Eigen to perfom all geometrical computation ( in the autodiff functor, seelibmv_homeography))
    DONE _ real using of id of element (not just index in array)
    DONE serious parameter reading
    NOT POSSIBLE use LocalParameterization instead of what I use now : should greatly improve quality, robustness and speed

    _ use a 2D cost function output for orthogonal distance, put a term penalyzing being inside an observation width zone
    _ add a constraint on distance to orginal angles between segements.

    _use boost for config files instead of lib found on internet

  @DEBUG
        remvoe this line from manualdistorth :" const double* w_i_j_;//! the width of the edge, in meter"


  */

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "glog/logging.h"
#include "Eigen/Core"

#include <string>
#include "Parameters.h"
#include "Data.h"                       //for data input/storage
#include "Constraints.h"                //define all the cost functor to be used in otpimization
#include "utils_function.h"             //for easier vect operation @TODO : we should work on templated eighen like in libmv_homography.cc
#include "WritingTempResultCallback.h"  //functor to write temp result in file at each iteration.
#include "geometry_function.h"

using ceres::NumericDiffCostFunction;
using ceres::CENTRAL;
//using ceres::FORWARD;
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
typedef Eigen::Vector3d Vec3;
typedef const Eigen::Vector3d Vec3const;
typedef Eigen::Map<const Eigen::Vector3d> mcVec3 ;
//typedef Eigen::Map<Eigen::Vector3f> VectorRef;
//typedef Eigen::Map<const Eigen::Vector3f> ConstVectorRef;


DataStorage * g_data_pointer ;
Parameter * g_param;

int main(int argc, char** argv) {


    //creating the set of parameters (could be read from file)
    g_param = new Parameter();
    g_param->readParameters();
    //std::cout << g_param->printParameters();





    std::cout << "  \E[34;1mReading data\E[m \n" ;
    DataStorage * data = new DataStorage(  g_param->input_file_path,g_param->output_file_path) ;
    g_data_pointer = data;

    //reading the classification file :
    data->readClassifications() ;

    //reading the network data ;
    data->readData();

    //setting the mapping beetween node_id and node*
    std::cout << "mapping between id of objects and objects" <<"\n";
    data->setMap();

    //testing the map on classifications :

    //string s = "LINESTRING(0 0 0 , 1 1 1, 2 2 2 )" ;
    string s = "POLYGON((1 3 50 , 3 3 50, 3 4 50, 1 4 50 ,1 3 50) )" ;
    read_WKT(s);

    return 1; /// @temp @debug


    //creating the problem to be solved
    Problem problem;


    //setting constraint on initial position for each node.
    if (g_param->use_initial_position_constraint == true) {
        addConstraintsOnInitialPosition( data, &problem) ;
    }

    //setting constraint on initial spacing between nodes for each pair of node.
    if(g_param->use_initial_spacing_constraint==true){
        addConstraintsOnInitialspacing( data, &problem) ;
    }

    // constraints based on observation : oth distance from observation to segment
    if(g_param->use_distance_to_proj_constraint == true){
        addConstraintsOnOrthDistToObservation(data , &problem) ;
    }

    // constraints based on observation : oth distance from observation to segment
    if(g_param->use_manual_distance_to_proj_constraint == true){
        addManualConstraintsOnOrthDistToObservation(data, &problem);
    }

    // constraints based on initial angle between edges
    if(g_param->use_manual_distance_to_original_angle == true){
        addManualConstraintsOnDistanceToOriginalAngle(data, &problem);
    }

    //manual constraints regularisation on distance between nodes
    if(g_param->use_manual_initial_spacing_constraint == true){
        addManualConstraintsOnInitialspacing(data, &problem);
    }

    Solver::Options options;
    options.max_num_iterations = 50;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;

    options.minimizer_type = ceres::LINE_SEARCH ; //can also be : TRUST_REGION or LINE_SEARCH
    //options.num_threads = 2; /// @todo : handy for speed, but makes it hard to understand cout

    //options.trust_region_strategy_type = ceres::DOGLEG ;
    //options.dogleg_type = ceres::SUBSPACE_DOGLEG ;
    //options.use_inner_iterations =true ;
    //options.use_approximate_eigenvalue_bfgs_scaling = true;

    //when stop the solver :

    //options.function_tolerance = 0.1;
    //options.gradient_tolerance= 0.001*0.001 ;
    //options.parameter_tolerance = 0.005;//std::pow(10,-10) ; // stop when the improvment is less than a millimeter

    ////output writing
    options.update_state_every_iteration= true ;
    WritingTempResultCallback callback(data->output_file_path(),0);
    options.callbacks.push_back(&callback);

    Solver::Summary summary;
    Solve(options, &problem, &summary);

    // std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";
    std::cout << summary.IsSolutionUsable() << "\n";

    return 0;
}

