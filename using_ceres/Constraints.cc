//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////
/** File to create constraint block from cost function

  */
#include "ceres/ceres.h"
#include "ceres/rotation.h"

#include "Constraints.h" //the header for constraint inclusion
#include "Data.h"
#include "utils_function.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;


//setting constraint on initial position for each node.
int addConstraintsOnInitialPosition(DataStorage * data, Problem * problem){
    for(const auto& element : data->nodes_by_node_id()){
        //std::cout << element.second->end_node << std::endl;
        node * n = element.second;

        double n_p[3] = { //! @todo : use eighen to hide this ugliness!
                          //! @note : we slighty pertubate the original position to try to improve initial solution
                          n->position[0] +0.001
                          ,n->position[1]+0.001
                          ,n->position[2]+0.001};
        DistanceToInitialPosition* self_distance_functor =
                new DistanceToInitialPosition( n->position) ;
        CostFunction* distance_cost_function
                = new AutoDiffCostFunction<DistanceToInitialPosition, 1,3>(
                    self_distance_functor);
        problem->AddResidualBlock(
                    distance_cost_function
                    ,NULL
                    ,n->position
                    );
    }
}



//setting constraint on initial spacing between nodes for each pair of node.
int addConstraintsOnInitialspacing(DataStorage * data, Problem * problem){
    for(const auto& element : data->edges_by_edge_id()){
        //std::cout << element.second->end_node << std::endl;
        edge * edge_to_output = element.second;
        node * start_node = data->nbn(edge_to_output->start_node) ;
        node * end_node = data->nbn(edge_to_output->end_node) ;


        double o_s[3] = { //! @todo : use eighen to hide this ugliness!
                          //! @note : we slighty pertubate the original position to try to improve initial solution
                          start_node->position[0]   - end_node->position[0] +0.001
                          ,start_node->position[1]  - end_node->position[1]+0.001
                          ,start_node->position[2]  - end_node->position[2]+0.001};
        DistanceToInitialSpacing* original_spacing_distance_functor = new DistanceToInitialSpacing( o_s) ;

        CostFunction* original_spacing_distance_cost_function
                = new AutoDiffCostFunction<DistanceToInitialSpacing, 1,3,3>(
                    original_spacing_distance_functor);
        problem->AddResidualBlock(
                    original_spacing_distance_cost_function
                    ,NULL
                    , start_node->position
                    , end_node->position
                    );
    }
}

//! constraint based on observation
int addConstraintsOnOrthDistToObservation(DataStorage * data, Problem * problem){
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
        problem->AddResidualBlock(
                    distance_cost_function
                    ,NULL
                    ,start_node->position
                    ,end_node->position
                    ); //note : both observations are referring to these nodes.


    }
}

