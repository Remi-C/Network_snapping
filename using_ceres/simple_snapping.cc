//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////
/**
  @todo
    _ use Eigen to perfom all geometrical computation, in the autodiff functor. (seelibmv_homeography)
    _ real using of id of element (not just index in array)
    _ serious parameter reading

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

  //getting the data ;
  std::cout << "  \E[34;1mReading data\E[m \n" ;
  DataStorage * data = new DataStorage(  g_param->input_file_path,g_param->output_file_path) ;
  g_data_pointer = data;
  data->readData();
  //setting the mapping beetween node_id and node*
  std::cout << "mapping between node_id and node *" <<"\n";
  data->setMap();
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

  Solver::Options options;
  options.max_num_iterations = 500;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = true;

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

