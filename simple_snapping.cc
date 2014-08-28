//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////

/**
  Abstract : snapping a street axis to observations

  _ Original Data :
    street axis, along with width information
    observations of sidewalk, along with confidence and weight information

  _ Input data :
    street axis node (n_i)
    connectivity information, along with width w_i_j of the part of the axis (e_i_j)

    observations of sidewalk, related to 2 nodes ; points o_k(n_i, n_j )

  _ optimisation model data
    n_i(X,Y,Z)
    w_i_j(n_i, n_j)

  _ optimisation constraints (K are normalised weight,Sum(K)=1, so that we can parametrize importance of constraints)
    _ K_n * intial position of n_i
    _ K_n_n * initial direction between n_i, n_j
    _ K_w * initial value of w_i_j
    _ K_d * distance of o_k to line (n_i,n_j)

  _ optimisation cost function
    _ a straightforward euclidian distance function (we can use the squared euclidian distance for faster computation) d_eucl(R_3,R_3)->R
    _ a

_ optimisation results :

  */
    /*
  simple file to put together all the idea about how to use ceres

  */

/*


  */

//////// creating the variables to be optimized : setting ParametterBlock

    /*
      //from bundle_adjuster.cc

      ceres::ParameterBlockOrdering* ordering =
          new ceres::ParameterBlockOrdering;

          //add the position of p_i
            for (int i = 0; i < num_points; ++i) {
              ordering->AddElementToGroup(points + point_block_size * i, 0);

          //add the width of the axis
            for (int i = 0; i < num_edges; ++i) {
              ordering->AddElementToGroup( 1 * i, 0);


       options->linear_solver_ordering.reset(ordering);
    }
    */


/*

  ///////////////PSEUDO CODE  :

  _ Structure of the model
    _ define Loss function
          // Configure the loss function.
          LossFunction* loss = NULL;
          if (FLAGS_robust_threshold) {
            loss = new CauchyLoss(FLAGS_robust_threshold);
          }


    _ Functor : class performing the error computation in a function and old some values
        create functor

    _ Cost function
        create cost function using functor



    add residual block (function on parameters)

  */


///////// creating a fist constraint on n_i_j : it shouldn't go too far away from it's original position
    /*


    */

 
// Author: Rémi Cura

#include "ceres/ceres.h"
#include "glog/logging.h"
#include "Eigen/Core"


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
const double observation_position[] = {
  4, 2, 0,
  2, -1, 0
};
double node_initialisation[] = {
  1, 1, 1,
  5, -1, -1
};
double node_position[] = {
  0, 0, 0,
  6, 0, 0
};

///////////////////////PARAMETERS///////////////////////
const double K_origin = 0.1 ; //! this parameter scale the distance to origin for a node

////////////////////////////////////////////////////////


/** functor to compute cost between one observation and the segment between 2 nodes
  */
struct DistanceToProjectionResidual {

   //! this is the constructor, it expects an array of at least 3 doubles.
    DistanceToProjectionResidual(const double* input_vect)
        :position_(input_vect) {}

    //@TODO : can't make it work , grrr
   // ~DistanceToProjectionResidual(){
   //    //destructor
   //    delete position_ ;
   //}


    //! this is the operator computing the cost, ie the distance projeted on normal of (n_i,n_j)
    bool operator()(const double* const n_i, /**< the first node of the network*/
                      const double* const n_j,/**< the second node of the network*/
                      double* distance_to_axis) const /**< this is the cost to be optimised*/
        {

        //converting input double array to EIgen 3D vector, to be able to use poweerfull eigen functions
        mcVec3 observation(position_);
        mcVec3 n_i_vect(n_i);
        mcVec3 n_j_vect(n_j);

        //here is the distance to axis (n_i, n_j) following (n_i,n_j normal)
            //Note that there is an offset that is the road width.
        distance_to_axis[0] = (observation-n_i_vect).cross(observation-n_j_vect).squaredNorm()
               /  (n_i_vect-n_j_vect).squaredNorm() -1;

        // @TEST : test distance : sum of distance to both nodes. Simpler distance for test purpose.
        //distance_to_axis[0] = (n_i_vect-observation).squaredNorm() + (observation-n_j_vect).squaredNorm() ;
        return true;
      }



 private:
    const double* position_; /**< store the 3D position of the observation point */
};


/** functor to compute cost between one node position and this node original position
  */

struct DistanceToInitialPosition {
   //! this is the constructor, it expects an array of 3 doubles to fill the original position
    DistanceToInitialPosition(const double* input_vect)
        :initial_position_(input_vect) {
    }

    //this is the operator computing the cost
//    bool operator()(const double* const n_i, /**< the current position of this node*/
//                      double* distance_to_origin) const /**< this is the cost to be optimised*/
//        {
//        //converting input double array to EIgen 3D vector, to be able to use poweerfull eigen functions
//        mcVec3 observation(initial_position_);
//        mcVec3 n_i_vect(n_i);
//        // note : eucl distance to the origin of this node.
//        distance_to_origin[0] = (observation-n_i_vect).squaredNorm()+1;
//        std::cout <<"\n : distance to origin : "<<distance_to_origin[0] <<"\n";
//      }

    //! distance computing function (eucl distance), this is templated.
    template <typename T> bool operator()(const T* const n_i,
                                        T* distance_to_origin) const {
    T s[3] ;
    s[0] =n_i[0]-T(initial_position_[0]);
    s[1] =n_i[1]-T(initial_position_[1]);
    s[2] =n_i[2]-T(initial_position_[2]);
    distance_to_origin[0] = K_origin * (s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);

    return true;
  }
 private:
    const double* initial_position_;
};


int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);


  // the node of the network : it should be an nput, this is just a test case

  //creating the problem to be solved
  Problem problem;

  //filling the problem with constraints to be optimized

  //setting constraint on initial position for each node.
  for (int k = 0; k <jNumNodes; ++k ){

      //        DistanceToInitialPosition* self_distance_functor = new DistanceToInitialPosition( &node_position[3 * k]) ;
      //        CostFunction* distance_cost_function
      //            = new NumericDiffCostFunction<DistanceToInitialPosition, FORWARD, 1, 3>(
      //                self_distance_functor);

      //          problem.AddResidualBlock(
      //              distance_cost_function
      //              ,NULL
      //              ,&node_position[3 * k]
      //              ); //note : both observations are referring to these nodes.

      DistanceToInitialPosition* self_distance_functor = new DistanceToInitialPosition( &node_position[3 * k]) ;
      CostFunction* distance_cost_function
          = new AutoDiffCostFunction<DistanceToInitialPosition, 1,3>(
              self_distance_functor);
        problem.AddResidualBlock(
            distance_cost_function
            ,NULL
            ,&node_position[3 * k]
            );
  }

  for (int i = 0; i < kNumObservations; ++i) {


      DistanceToProjectionResidual* distance_functor = new DistanceToProjectionResidual( &observation_position[3 * i]) ;

      CostFunction* distance_cost_function
          = new NumericDiffCostFunction<DistanceToProjectionResidual, CENTRAL, 1, 3, 3>(
              distance_functor);

        problem.AddResidualBlock(
            distance_cost_function
            ,NULL
            ,&node_position[0]
            ,&node_position[3]
            ); //note : both observations are referring to these nodes.
  }

  Solver::Options options;
  options.max_num_iterations = 50;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = true;

  Solver::Summary summary;
  Solve(options, &problem, &summary);

  std::cout << summary.BriefReport() << "\n";
  std::cout << summary.FullReport() << "\n";
  std::cout << summary.IsSolutionUsable() << "\n";

    for (int k = 0; k <jNumNodes; ++k ){
     std::cout << "NEw n_"<<k<<": \n"
                 << node_position[3 * k] <<": \n"
                << node_position[3 * k+1]  <<": \n"
                << node_position[3 * k+2]
               << "\n ";

   }
  return 0;
}

