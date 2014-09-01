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

const int kNumObservations = 2;
const int jNumNodes = 2;
const int uNumPairs =1 ;
const double observation_position[] = { //! the initial position of observations . No change
  4, 2, 0,
  2, -1, 0
};
double node_position[] = { //! the initial node position. This parameters are going to be optiized
  0, 0, 0,
  6, 0, 0
};
const int node_pair[1][2]={ //! the coupling between nodes to form segments
    {0,1}
};

//! @TODO @WARNING very ugly : should be set in a "parameter class"
///////////////////////PARAMETERS///////////////////////
const string output_file_path("/media/sf_E_RemiCura/PROJETS/snapping/visu/visu_inqgis_timemanager/simple_test.csv") ;
const string input_file_path("/media/sf_E_RemiCura/PROJETS/snapping/using_ceres/data/simple_example.csv") ;

const double K_origin = 1 ; //! this parameter scale the distance to origin for a node
const double K_obs= 1 ; //! this parameter scale the measure of distance between observation and line (n_i,n_j)
const double K_spacing= 1 ; //! this parameter scale the measure of similarity between [n_i,n_j] original and current

////////////////////////////////////////////////////////




int main(int argc, char** argv) {

  std::cout << "  \E[34;1mReading data\E[m \n" ;


  //getting the data ;
  DataStorage * data = new DataStorage(  input_file_path,output_file_path) ;
  data->readData();
  //setting the mapping beetween node_id and node*
  std::cout << "mapping between node_id and node *" <<"\n";
  data->setMap();

  std::cout << data->nbn(2)->nodeToString() <<endl;

  //creating the problem to be solved
  Problem problem;

  //filling the problem with constraints to be optimized

  //setting constraint on initial position for each node.
  for (int k = 0; k <jNumNodes; ++k ){

      DistanceToInitialPosition* self_distance_functor =
                new DistanceToInitialPosition( data->nodes(k)->position) ;
      CostFunction* distance_cost_function
          = new AutoDiffCostFunction<DistanceToInitialPosition, 1,3>(
              self_distance_functor);
        problem.AddResidualBlock(
            distance_cost_function
            ,NULL
            ,data->nodes(k)->position
            );
  }

  //setting constraint on initial spacing between nodes for each pair of node.
  for (int u = 0; u <uNumPairs; ++u ){
      node * start_node =  data->nbn(data->edges(u)->start_node);
      node * end_node =  data->nbn(data->edges(u)->end_node);
      double o_s[3] = {
          start_node->position[0]  - end_node->position[0]
          ,start_node->position[1]  - end_node->position[1]
          ,start_node->position[2]  - end_node->position[2]};
      DistanceToInitialSpacing* original_spacing_distance_functor = new DistanceToInitialSpacing( o_s) ;

      CostFunction* original_spacing_distance_cost_function
          = new AutoDiffCostFunction<DistanceToInitialSpacing, 3,3,3>(
              original_spacing_distance_functor);
        problem.AddResidualBlock(
            original_spacing_distance_cost_function
            ,NULL
            , start_node->position
            , end_node->position
            );
  }

  for (int i = 0; i < kNumObservations; ++i) {

      DistanceToProjectionResidual* distance_functor =
                new DistanceToProjectionResidual( data->observations(i)->position  ) ;
      CostFunction* distance_cost_function
          = new AutoDiffCostFunction<DistanceToProjectionResidual, 1, 3, 3>(
              distance_functor);
        problem.AddResidualBlock(
            distance_cost_function
            ,NULL
            ,data->nodes(0)->position
            ,data->nodes(1)->position
            ); //note : both observations are referring to these nodes.
  }

  Solver::Options options;
  options.max_num_iterations = 50;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = true;

  //output writing option
  options.update_state_every_iteration= true ;
  WritingTempResultCallback callback(output_file_path,0);
  options.callbacks.push_back(&callback);
  //
  Solver::Summary summary;
  Solve(options, &problem, &summary);

  std::cout << summary.BriefReport() << "\n";
  std::cout << summary.FullReport() << "\n";
  std::cout << summary.IsSolutionUsable() << "\n";

    for (int k = 0; k <jNumNodes; ++k ){
        std::cout << "result :"<<endl
          << data->nodes(k)->nodeToString() <<endl;

   }



  return 0;
}

