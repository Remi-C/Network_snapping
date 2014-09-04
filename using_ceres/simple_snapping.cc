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


//! @TODO @WARNING very ugly : should be set in a "parameter class"
///////////////////////PARAMETERS///////////////////////
const string input_file_path("../data/data_in_reduced_export_area/reduced_area.csv");
const string  output_file_path("../data/data_in_reduced_export_area/snapping_output.csv");

const double K_origin = 1  ; //! this parameter scale the distance to origin for a node
const double K_obs= 1 ; //! this parameter scale the measure of distance between observation and line (n_i,n_j)
const double K_spacing=  1 ; //! this parameter scale the measure of similarity between [n_i,n_j] original and current

const bool use_initial_position_constraint = true;
const bool use_initial_spacing_constraint = false;
const bool use_distance_to_proj_constraint = true;
////////////////////////////////////////////////////////

//////////////// gloable variable. Should use a singleton instead
DataStorage * data_pointer ;
int addConstraintsOnInitialPosition(DataStorage *, Problem * );
int addConstraintsOnInitialspacing( DataStorage *, Problem * );
int addConstraintsOnOrthDistToObservation(DataStorage * , Problem * ) ;
////////////////



int main(int argc, char** argv) {

  std::cout << "  \E[34;1mReading data\E[m \n" ;


  //getting the data ;
  DataStorage * data = new DataStorage(  input_file_path,output_file_path) ;
  data_pointer = data;
  data->readData();
  //setting the mapping beetween node_id and node*
  std::cout << "mapping between node_id and node *" <<"\n";
  data->setMap();
  //creating the problem to be solved
  Problem problem;


    //setting constraint on initial position for each node.
  if (use_initial_position_constraint == true) {
    addConstraintsOnInitialPosition( data, &problem) ;
  }

//setting constraint on initial spacing between nodes for each pair of node.
  if(use_initial_spacing_constraint==true){
      addConstraintsOnInitialspacing( data, &problem) ;
  }

  //! constraint based on observation
  if(use_distance_to_proj_constraint == true){
  for (int i = 0; i < data->num_observations(); ++i) {
       //finding the 2 nodes concerned by this observations
      edge * relativ_edge = data->ebe(data->observations(i)->edge_id) ;
      node * start_node = data->nbn(relativ_edge->start_node)  ;
      node * end_node = data->nbn(relativ_edge->end_node)  ;

      DistanceToProjectionResidual* distance_functor =
              new DistanceToProjectionResidual( data->observations(i)->position, &relativ_edge->width, data->observations(i)  ) ;
      CostFunction* distance_cost_function
          = new AutoDiffCostFunction<DistanceToProjectionResidual, 1, 3, 3>(
              distance_functor);
        problem.AddResidualBlock(
            distance_cost_function
            ,NULL
            ,start_node->position
            ,end_node->position
            ); //note : both observations are referring to these nodes.


  }
  }

  Solver::Options options;
  options.max_num_iterations = 500;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = true;

  ////output writing option
  options.update_state_every_iteration= true ;
  WritingTempResultCallback callback(output_file_path,0);
  options.callbacks.push_back(&callback);

  Solver::Summary summary;
  Solve(options, &problem, &summary);

  std::cout << summary.BriefReport() << "\n";
  std::cout << summary.FullReport() << "\n";
  std::cout << summary.IsSolutionUsable() << "\n";

  return 0;
}

