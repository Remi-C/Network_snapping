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
//using ceres::AutoDiffCostFunction;
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
const double data[] = {
  4, 2, 0,
  2, -1, 0,
};

struct DistanceToProjectionResidual {

   //this is the constructor, it expects an array of at least 3 doubles.
    DistanceToProjectionResidual(const double* input_vect)
        :position_(input_vect) {}

    //@TODO : can't make it work , grrr
   // ~DistanceToProjectionResidual(){
   //    //destructor
   //    delete position_ ;
   //}


    //this is the operator computing the cost
    bool operator()(const double* const n_i,
                      const double* const n_j,
                      double* distance_to_axis) const {

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
    const double* position_;
};


int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);


  // the node of the network
  Vec3 n_1 ; n_1 << 0,0,0 ;
  Vec3 n_2 ; n_2 << 6,0,0 ;
  Vec3 initial_n_1 = n_1;
  Vec3 initial_n_2 = n_2;


  Problem problem;
  for (int i = 0; i < kNumObservations; ++i) {


  DistanceToProjectionResidual* distance_functor = new DistanceToProjectionResidual( &data[3 * i]) ;

  CostFunction* distance_cost_function
      = new NumericDiffCostFunction<DistanceToProjectionResidual, CENTRAL, 1, 3, 3>(
          distance_functor);

    problem.AddResidualBlock(
        distance_cost_function
        ,NULL
        ,( double * ) n_1.data()
        ,( double * ) n_2.data()
        ); //note : both observations are referring to these nodes.
  }

  Solver::Options options;
  options.max_num_iterations = 25;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = true;

  Solver::Summary summary;
  Solve(options, &problem, &summary);
  std::cout << summary.BriefReport() << "\n";
  std::cout << summary.FullReport() << "\n";
  std::cout << summary.IsSolutionUsable() << "\n";
  std::cout << "Initial n1: \n" << initial_n_1 << "\n n2: \n" << initial_n_2 << "\n";
  std::cout << "Final   n1: \n" << n_1 << "\n _n n2: \n" << n_2 << "\n";
  return 0;
}

