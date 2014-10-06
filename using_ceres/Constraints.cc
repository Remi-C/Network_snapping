//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////
/** File to create constraint block from cost function

  */
#include "Constraints.h" //the header for constraint inclusion


extern Parameter* g_param;



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

        //untill 2.0 meters of distance, normal behavior. after that outliers behavior (not square)
        LossFunction* loss = NULL;
        loss = new ceres::ScaledLoss(
                    //g_param->useLoss?(new ceres::SoftLOneLoss(g_param->lossScale)):NULL
                    /// @DEBUG : temporary
                    NULL
                    ,g_param->K_origin,ceres::DO_NOT_TAKE_OWNERSHIP) ;

        problem->AddResidualBlock(
                    distance_cost_function
                    ,loss
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

        //untill 2.0 meters of distance, normal behavior. after that outliers behavior (not square)
        LossFunction* loss = NULL;
        loss = new ceres::ScaledLoss( g_param->useLoss?(new ceres::SoftLOneLoss(g_param->lossScale)):NULL
                                        ,g_param->K_spacing,ceres::DO_NOT_TAKE_OWNERSHIP) ;


        problem->AddResidualBlock(
                    original_spacing_distance_cost_function
                    ,loss
                    , start_node->position
                    , end_node->position
                    );
    }
}

//! manual constraint based on regularisation of distance between nodes
int addManualConstraintsOnInitialspacing(DataStorage * data, Problem * problem){
    for(const auto& element : data->edges_by_edge_id()){
        //std::cout << element.second->end_node << std::endl;

        edge * edge_to_output = element.second;
        node * start_node = data->nbn(edge_to_output->start_node) ;
        node * end_node = data->nbn(edge_to_output->end_node) ;

        VectorRef sNode(start_node->position,3);
        VectorRef eNode(end_node->position,3);
        //std::cout << "sNode : " <<sNode.transpose() << std::endl ;
        //std::cout << "eNode : " <<eNode.transpose() << std::endl ;

        Eigen::Vector3d o_s = sNode-eNode  ;//+ VectorRef::Random(); /// @DEBUG : warning : this should be maybe 0.001
        CostFunction* original_spacing_distance_functor = new ManualOriginalSpacing( o_s) ;


        LossFunction* loss = NULL;
        loss = new ceres::ScaledLoss( g_param->useLoss?(new ceres::SoftLOneLoss(g_param->lossScale)):NULL
                                        ,g_param->K_spacing,ceres::DO_NOT_TAKE_OWNERSHIP) ;

        problem->AddResidualBlock(
                     original_spacing_distance_functor
                    ,loss
                    ,start_node->position
                    ,end_node->position
                    ); //note : both observations are referring to these nodes.
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

        //untill 2.0 meters of distance, normal behavior. after that outliers behavior (not square)
        LossFunction* loss = NULL;
        loss = new ceres::ScaledLoss( g_param->useLoss?(new ceres::SoftLOneLoss(g_param->lossScale)):NULL
                                        ,g_param->K_obs,ceres::DO_NOT_TAKE_OWNERSHIP) ;

        problem->AddResidualBlock(
                    distance_cost_function
                    ,loss
                    ,start_node->position
                    ,end_node->position
                    ); //note : both observations are referring to these nodes.
    }
}



//! manual constraint based original angles in network
int addManualConstraintsOnDistanceToOriginalAngle(DataStorage * data, Problem * problem){

    for (int i = 0 ; i< data->num_nodes(); ++i){//for every node,
        int current_node_id = data->nodes(i)->node_id ;
        //std::cout << "curr node " << current_node_id <<std::endl;
        auto range = data->edges_by_node_id()->equal_range(current_node_id);

        for (auto it = range.first; it != range.second; ++it) {//for every node, loop on all edges concerned
            edge * first_edge = data->ebe(it->second->edge_id);
            node * first_node = first_edge->start_node!=current_node_id?
                        data->nbn(first_edge->start_node)
                      :data->nbn(first_edge->end_node);


            for (auto it2 = it; it2 != range.second; it2++) {//generating all unique pair of edges for a node
                if(it->second->edge_id!=it2->second->edge_id) {
                    //std::cout << "curr edge " << it->second->edge_id <<", sec edge pair : " << it2->second->edge_id <<std::endl;

                    edge * sec_edge = data->ebe(it2->second->edge_id);
                    node * sec_node = sec_edge->start_node!=current_node_id?
                                data->nbn(sec_edge->start_node)
                              :data->nbn(sec_edge->end_node);

                    //note : the node to input are then : Nc : data->nbn(current_node_id) , Ni : first_node , Nj : sec_node

                    //computing the original angle  :
                    ConstVectorRef Nc( data->nbn(current_node_id)->position,3 );
                    ConstVectorRef Ni( first_node->position ,3 );
                    ConstVectorRef Nj( sec_node->position ,3 );
                    Eigen::Vector3d pe;//initial perturbation
                    pe << 0.001,0.001,0.001 ;
                    double scalar_a = (Nc-Ni+pe).dot(Nc-Nj+pe)/((Nc-Ni+pe).norm() * (Nc-Nj+pe).norm());
                    double cross_a = ((Nc-Ni+pe).cross(Nc-Nj+pe)/((Nc-Ni+pe).norm() * (Nc-Nj+pe).norm())).norm();

                    std::cout << "input for block : "<< scalar_a << "," << cross_a << std::endl;
                    std::cout << "nodes : " << current_node_id <<","<<first_node->node_id <<"," << sec_node->node_id << std::endl;
                    CostFunction* distance_cost_function=
                            new  ManualDistanceToOriginalAngle(scalar_a,cross_a ) ;
                    //untill 2.0 meters of distance, normal behavior. after that outliers behavior (not square)
                    LossFunction* loss = NULL;
                    loss = new ceres::ScaledLoss( g_param->useLoss?(new ceres::SoftLOneLoss(g_param->lossScale)):NULL
                                                                   ,g_param->K_angle,ceres::DO_NOT_TAKE_OWNERSHIP) ;

                    problem->AddResidualBlock(
                                distance_cost_function
                                ,loss
                                ,data->nbn(current_node_id)->position
                                ,first_node->position
                                ,sec_node->position
                                );
                }
            }
        }
    }
}





//! manual constraint based on regularisation of angles between edges
int addManualConstraintsOnOrthDistToObservation(DataStorage * data, Problem * problem){
    for (int i = 0; i < data->num_observations(); ++i) {

        //finding the 2 nodes concerned by this observations
        edge * relativ_edge = data->ebe(data->observations(i)->edge_id) ;
        node * start_node = data->nbn(relativ_edge->start_node)  ;
        node * end_node = data->nbn(relativ_edge->end_node)  ;

        CostFunction* distance_cost_function=
            new  ManualOrthDistanceToObservation( data->observations(i)->position, &relativ_edge->width, data->observations(i)  ) ;
        //untill 2.0 meters of distance, normal behavior. after that outliers behavior (not square)
        LossFunction* loss = NULL;
        loss = new ceres::ScaledLoss( g_param->useLoss?(new ceres::SoftLOneLoss(g_param->lossScale)):NULL
                                        ,g_param->K_obs,ceres::DO_NOT_TAKE_OWNERSHIP) ;

        problem->AddResidualBlock(
                    distance_cost_function
                    ,loss
                    ,start_node->position
                    ,end_node->position
                    ); //note : both observations are referring to these nodes.
    }
}

