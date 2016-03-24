//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////
/**

  @STATE 2014 11 05
    Some serious issues
        DONE for sidewalk : should exclude the observations that are close to the intersections
        for objects, there is a design flaw in surf dist: the jacobian should be geometrical and return the change to do to go to 0 cost.
        for the moment it returns a somehting proportional to the shared surface !
        For objects : started a new pure distance to object function, but it seems to be not working
    It is impossible to properly add new measures and balances cost without a visualization of everything at each step.
    There is no easy way to visualize.
    At each iteration, we have to manually evaluate each cost function with proper arguments, then output the processed result into a file to see what are the constraints
    We woul need : for each cost function : to represent the jacobian
    for cost function on objects & suidewalk
        a linestring (oriented) from border of axis (includiong width) with jac direction and norm(jac) size
        We can create such line  : obs+jac should be a point exactly on axis dilated border. We use this point as first point, of the line, then offset it by -jac to create the line representing jac.
        we do not represent regularization constraints.




  @TODO
    _ use Eigen to perfom all geometrical computation ( in the autodiff functor, seelibmv_homeography))
    DONE _ real using of id of element (not just index in array)
    DONE serious parameter reading
    NOT POSSIBLE use LocalParameterization instead of what I use now : should greatly improve quality, robustness and speed

    _ use a 2D cost function output for orthogonal distance, put a term penalyzing being inside an observation width zone

    DONE add a constraint on distance to orginal angles between segements.
     - allow same observation id (allow an observation to be affected to multiple edges)
     - Generate visualisation data for every cost and jacobian
        DONE DistToObs
        DONE OriginalSpacing
        DONE OriginalAngle
        TODO OjectDistSurf

     - Allow fixed nodes for border of the network
        - Allow different behaviour based on the node being in intersection or not.

  @DEBUG
    _cosntraint on angle --> jacobians are way way too big.
     probably a math error somewhere

    _Border /in/out/ non working. Sign problem(again)?
    _ surface dist : big design mistake : the jacobian and residual should be proportionnal to the distance between
       object and edge ofsseted by width, and not proportionnal to area, because
       The power of the induced moves are not realted to geometric reality

    _OrthDist : big problem with all edges becoming verticals. A mistake somewhere in the jac evaluation?
       Cost and jac remains small even if the network is radically modified, what's wrong?

    _sharedArea :
        for OUT behavior, seems to be buggy when fully inside

    _objectDist
        some Z coordinates of the Jac are put to Nan

    _

  */

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "glog/logging.h"
#include "Eigen/Core"

#include <ctime>
#include <string>
#include "Parameters.h"
#include "Data.h"                       //for data input/storage
#include "Constraints.h"                //define all the cost functor to be used in otpimization
#include "utils_function.h"             //for easier vect operation @TODO : we should work on templated eighen like in libmv_homography.cc
#include "WritingTempResultCallback.h"  //functor to write temp result in file at each iteration.
#include "geometry_function.h"
#include "enum_functions.h"

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
    initialize_geom_computation();
     std::time_t time_begining = std::time(nullptr);

    //creating the set of parameters (could be read from file)
    std::cout << "  \E[34;1mReading parameters\E[m \n" ;
    g_param =Parameter::instance() ;
    g_param->readParameters();
    //std::cout << g_param->printParameters();


    std::cout << "  \E[34;1mReading data\E[m \n" ;
    DataStorage * data = new DataStorage(  g_param->input_file_path,g_param->output_file_path) ;
    g_data_pointer = data;

    //reading the classification file :
    std::cout << "  \E[34;1m \t Reading classification\E[m \n" ;
    data->readClassifications() ;

    //reading the network data ;
    std::cout << "  \E[34;1m \t Reading network\E[m \n" ;
    data->readData();

    //setting the mapping beetween node_id and node*
    std::cout << "mapping between id of objects and objects" <<"\n";
    data->setMap();

    //reading the objects for snapping
    std::cout << "  \E[34;1m \t Reading Objects\E[m \n" ;
    data->readObjects();

    //constructing problem
    std::cout << "  \E[34;1m \tConstructing Problem\E[m \n" ;
    /// clean version to not allow ceres to destroy the memory itself
    Problem::Options pb_options;
    pb_options.local_parameterization_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
    //creating the problem to be solved
    Problem problem(pb_options);

    addAllParameterBlocks(data,&problem,g_param) ;
    addAllConstraints(data,&problem,g_param);

    if(g_param->optimisation_method==ceres::TRUST_REGION){ // only possible to bound with trust_region
        boundConstraints(data,&problem,g_param) ;
    }

    Solver::Options options;
    Solver::Summary summary;


    options.max_num_iterations = 500;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    options.minimizer_progress_to_stdout = true;

    options.minimizer_type = g_param->optimisation_method ; //can also be : TRUST_REGION or LINE_SEARCH
    options.num_threads = 6; /// @todo : handy for speed, but makes it hard to understand cout

    options.line_search_direction_type = ceres::BFGS ;//   BFGS and LBFGS
    options.trust_region_strategy_type = ceres::DOGLEG ;
    options.dogleg_type = ceres::SUBSPACE_DOGLEG ;
    options.use_inner_iterations = true ;
    //options.use_approximate_eigenvalue_bfgs_scaling = true;

    //options.use_nonmonotonic_steps = true ;


    //when stop the solver :

    //options.function_tolerance = 0.1;
    //options.gradient_tolerance= 0.001*0.001 ;
    //options.parameter_tolerance = std::pow(10,-10) ; // stop when the improvment is less than a millimeter
    //options.max_num_consecutive_invalid_steps = 3;
    //options.max_num_iterations = 8;

    ////output writing
    options.update_state_every_iteration= true ;
    WritingTempResultCallback callback(data->output_file_path(),0);
    options.callbacks.push_back(&callback);

    //solving
    int n_iter = 0 ;
    std::cout << "  \E[34;1m \tSolving Problem\E[m \n" ;
    //Solve(options, &problem, &summary);
    //n_iter +=summary.iterations.size()-1 ;
    //depending on optimizing on width or position, it will be efficient to "freeze" the parameters that are not used (position, or width)
    //void Problem::SetParameterBlockConstant(double* values)
    //void Problem::SetParameterBlockVariable(double* values)
    g_data_pointer->writeData(1);
    g_data_pointer->writeConstraints(1);

    bool mixed_optim = g_param->optimisation_type==SnapEnums::MIXED ;
    if(mixed_optim==true)
         g_param->optimisation_type = SnapEnums::WIDTH;

    int stop_optim = 0;
    int n_loop = 0 ;
    do{
        std::cout << "\n optimizing on " << g_param->optimisation_type << "\n" ;
        activate_desactivate_ParameterBlocks( data,&problem,g_param);
        Solve(options, &problem, &summary);
        n_iter +=summary.iterations.size();
        if(summary.final_cost == summary.initial_cost) // this solving hasn't solved anything!
            stop_optim ++ ;

        if(mixed_optim == true)
            g_param->optimisation_type = g_param->optimisation_type==SnapEnums::WIDTH?SnapEnums::POSITION:SnapEnums::WIDTH ;
        std::cout << " summary.iterations.size() " << summary.iterations.size() << " n_iter " <<n_iter << " stop_optim " << stop_optim << " nlopp "  << n_loop << "\n" ;
        n_loop ++ ;
    }while((summary.iterations.size()-1)<=250
           && n_iter < 600
           && (stop_optim <2 || g_param->use_manual_distance_to_proj_constraint_width !=true )
           && n_loop < 1
          );
    // std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";

    std::cout <<n_iter<< std::endl ;

    g_data_pointer->writeData(2);
    g_param->optimisation_type = SnapEnums::POSITION ; //putting the position, so as to have vector in visu
    g_data_pointer->writeConstraints(2);
    g_param->optimisation_type = SnapEnums::WIDTH ; //putting the position, so as to have vector in visu
    g_data_pointer->writeConstraints(2);



    finish_geom_computation();
    std::time_t time_end  = std::time(nullptr);
    const std::time_t time_duration = time_end-time_begining ;
    printf("elapsed time : %ld hour %ld minutes %ld seconds \n ",time_duration/3600%24,time_duration/60%60,time_duration%60);
    return 0;
}

