//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////
/** File to create constraint block from cost function

  */
#include "Constraints.h" //the header for constraint inclusion


extern Parameter* g_param;


//creating all constraints
int addAllConstraints(DataStorage * data, ceres::Problem * problem, Parameter* param){
    //setting constraint on initial position for each node.
    //    if (param->use_initial_position_constraint == true) {
    //        addConstraintsOnInitialPosition( data, problem) ;
    //    }

    if (param->use_manual_initial_position_constraint == true) {
        addManualConstraintsOnInitialPosition(  data, problem) ;
    }

    if(param->use_manual_initial_spacing_constraint == true){
        addManualConstraintsOnInitialspacing(data, problem);
    }

    // constraints based on initial angle between edges
    if(param->use_manual_distance_to_original_angle == true){
        addManualConstraintsOnDistanceToOriginalAngle(data, problem);
    }


    if (param->use_manual_initial_width_constraint == true) {
        addManualConstraintsOnInitialWidth(data, problem) ;
    }

    // constraints based on observation : oth distance from observation to segment
    if(param->use_manual_distance_to_proj_constraint == true){
        addManualConstraintsOnOrthDistToObservation(data, problem);
    }

    //manual constraints regularisation on shared surface (plus maybe distance to border)
    if(param->use_manual_Surf_Dist_To_Objects_constraint == true){
        addManualConstraintsOnSurfDistToObjects(data, problem);
    }
	
	if(param->use_manual_target_slope == true){
		addManualTargetSlope(data, problem);
	}

    // constraints based on observation for width : oth distance from observation to segment
    if(param->use_manual_distance_to_proj_constraint_width == true){
        addManualConstraintsOnOrthDistToObservation_width(data, problem);
    }

    //manual constraints regularisation on shared surface (plus maybe distance to border)
    if(param->use_manual_Surf_Dist_To_Objects_constraint_width == true){
        addManualConstraintsOnSurfDistToObjects_width(data, problem);
    }
    return 0;
}

/** create bounds for optimisation parameter
  the coordinates cannot vary more than +- specified bound (geom_bound)
  the width cannot vary more than the +6 specified bound (width_bound)
  */
int boundConstraints(DataStorage * data, ceres::Problem * problem, Parameter* param ){

    //loop on all nodes and all dimensions
    for(int i=0;i<data->num_nodes();++i){
        for(int j=0;j<3;++j){
            problem->SetParameterLowerBound(data->nodes(i)->position, j
                                            , data->nodes(i)->position[j]-param->geom_bound
                                            ) ;
            problem->SetParameterUpperBound(data->nodes(i)->position, j
                                            , data->nodes(i)->position[j]+param->geom_bound
                                            ) ;
        }
    }

    for(int i=0;i<data->num_edges();++i){
        double t_lb = std::min(std::max(data->edges(i)->width[0]-param->width_bound_range,param->width_bound_minimal),param->width_bound_maximal );
        double t_ub = std::max(std::min(data->edges(i)->width[0]+param->width_bound_range,param->width_bound_maximal),param->width_bound_minimal);
        //        printf(" edge_id : %d , width : %f , lower bound : %f, upper bound : %f \n"
        //               ,data->edges(i)->edge_id,data->edges(i)->width[0]
        //               ,t_lb
        //               ,t_ub) ;
        problem->SetParameterLowerBound( data->edges(i)->width , 0
                                         , t_lb
                                         ) ;
        problem->SetParameterUpperBound(  data->edges(i)->width, 0
                                          , t_ub
                                          ) ;
    }
    return 0;
}


/** create all the variable for the optimisation. It ensure that variable exists even when not used in any constraints
  */
int addAllParameterBlocks(DataStorage * data, ceres::Problem * problem, Parameter* param ){

    //loop on all nodes and all dimensions
    for(int i=0;i<data->num_nodes();++i){
        problem->AddParameterBlock(data->nodes(i)->position,3);
    }

    for(int i=0;i<data->num_edges();++i){
        problem->AddParameterBlock(data->edges(i)->width,1);
    }
    return 0;
}

/** create all the variable for the optimisation. It ensure that variable exists even when not used in any constraints
  */
int activate_desactivate_ParameterBlocks(DataStorage * data, ceres::Problem * problem, Parameter* param ){

    //loop on all nodes and all dimensions
    for(int i=0;i<data->num_nodes();++i){
        //for(int j=0;j<3;++j){
//            param->optimisation_type==SnapEnums::POSITION?
//                        problem->SetParameterBlockVariable(&data->nodes(i)->position[j])
//                      :problem->SetParameterBlockConstant(&data->nodes(i)->position[j]) ;
        //}
        if(param->optimisation_type == SnapEnums::POSITION || param->optimisation_type == SnapEnums::MIXED){
                    problem->SetParameterBlockVariable(data->nodes(i)->position );
        }else{
        problem->SetParameterBlockConstant(data->nodes(i)->position ) ;
        }
    }
    for(int i=0;i<data->num_edges();++i){
        if(param->optimisation_type==SnapEnums::WIDTH || param->use_manual_distance_to_proj_constraint_width == true || param->use_manual_Surf_Dist_To_Objects_constraint_width == true ){
           problem->SetParameterBlockVariable(data->edges(i)->width);
        }else{
           problem->SetParameterBlockConstant(data->edges(i)->width) ;
        }
    }
    return 0;
}

//setting constraint on initial position for each node.
int addManualConstraintsOnInitialPosition(DataStorage * data, Problem * problem){
    double* neg = new double(-1);
    for(const auto& element : data->nodes_by_node_id()){
        //std::cout << element.second->end_node << std::endl;
        node * n = element.second;
        VectorRef sNode(n->position,3);
        sNode += Eigen::Vector3d::Constant(TOLERANCE);

        CostFunction* origin_distance_functor = new ManualDistanceToInitialPosition( sNode) ;

        LossFunction* loss = NULL;
        loss = new ceres::ScaledLoss( g_param->useLoss?(new ceres::SoftLOneLoss(g_param->lossScale)):NULL
                                                       ,g_param->K_origin,ceres::DO_NOT_TAKE_OWNERSHIP) ;

        problem->AddResidualBlock(
                    origin_distance_functor
                    ,loss
                    ,n->position
                    );

        Constraint * n_constraint = new Constraint_origin(
                    n->position
                    ,n->position//useless
                    ,n->position//useless
                    ,neg
                    ,origin_distance_functor
                    ,n->position
                    ) ;
        data->constraints()->push_back(n_constraint);
    }
    return 0;
}



//! manual constraint based on regularisation of distance between nodes
int addManualConstraintsOnInitialspacing(DataStorage * data, Problem * problem){
    double* neg = new double(-1);
    double* pos = new double(+1);
    int i = 0 ;
    for(const auto& element : data->edges_by_edge_id()){
        i++;
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

        //saving the constraint for reuse to output at each step
        //creating the constraint 2 times, so visualisation is easier
        Constraint * n_constraint = new Constraint_spacing(
                    start_node->position
                    ,end_node->position
                    ,start_node->position//useless
                    ,neg
                    ,original_spacing_distance_functor
                    ,start_node->position
                    ) ;

        Constraint * n_constraint2 = new Constraint_spacing(
                    start_node->position
                    ,end_node->position
                    ,start_node->position//useless
                    ,pos
                    ,original_spacing_distance_functor
                    ,end_node->position
                    ) ;
        //adding it to list of constraints
        data->constraints()->push_back(n_constraint);
        data->constraints()->push_back(n_constraint2);
    }
    return i;
}

//! manual constraint based original angles in network
int addManualConstraintsOnDistanceToOriginalAngle(DataStorage * data, Problem * problem){
    double * zer = new double(0);
    int count = 0 ;
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
                    pe << 0,0,0 ; //0.001,0.001,0.001 ;
                    double scalar_a = (Nc-Ni+pe).dot(Nc-Nj+pe)/((Nc-Ni+pe).norm() * (Nc-Nj+pe).norm());
                    double cross_a = ((Nc-Ni+pe).cross(Nc-Nj+pe)/((Nc-Ni+pe).norm() * (Nc-Nj+pe).norm())).norm();

                    //                    std::cout << "input for block : "<< scalar_a << "," << cross_a << std::endl;
                    //                    std::cout << "nodes : " << current_node_id <<","<<first_node->node_id <<"," << sec_node->node_id << std::endl;
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



                    //saving the constraint for reuse to output at each step
                    Constraint * n_constraint = new Constraint_angle(
                                data->nbn(current_node_id)->position
                                ,first_node->position
                                ,sec_node->position
                                ,zer //useless
                                ,distance_cost_function
                                ,data->nbn(current_node_id)->position
                                ) ;

                    //adding it to list of constraints
                    data->constraints()->push_back(n_constraint);
                    count++;
                }
            }
        }
    }
    return count ;
}

//setting constraint on initial position for each node.
int addManualConstraintsOnInitialWidth(DataStorage * data, Problem * problem){
    int count = 0;
    for(const auto& element : data->edges_by_edge_id()){
        //std::cout << element.second->end_node << std::endl;
        edge * e = element.second;

        CostFunction* origin_distance_functor = new ManualDistanceToInitialWidth(*e->width) ;

        LossFunction* loss = NULL;
        loss = new ceres::ScaledLoss( g_param->useLoss?(new ceres::SoftLOneLoss(g_param->lossScale)):NULL
                                                       ,g_param->K_original_width,ceres::DO_NOT_TAKE_OWNERSHIP) ;

        problem->AddResidualBlock(
                    origin_distance_functor
                    ,loss
                    ,e->width
                    );

        //computing the middle of the edge as the application point
        Constraint * n_constraint = new Constraint_original_width(
                    data->nbn(e->start_node)->position//used
                    ,data->nbn(e->end_node)->position//used
                    ,data->nbn(e->start_node)->position//useless
                    ,e->width
                    ,origin_distance_functor
                    ,data->nbn(e->start_node)->position
                    ) ;
        data->constraints()->push_back(n_constraint);
        count++;
    }
    return count;
}


//! manual constraint based on distance between observation and edges
int addManualConstraintsOnOrthDistToObservation(DataStorage * data, Problem * problem){
    int count = 0;
    for (int i = 0; i < data->num_observations(); ++i) {

        //finding the 2 nodes concerned by this observations
        edge * relativ_edge = data->ebe(data->observations(i)->edge_id) ;
        node * start_node = data->nbn(relativ_edge->start_node)  ;
        node * end_node = data->nbn(relativ_edge->end_node)  ;

        CostFunction* distance_cost_function=
                new  ManualOrthDistanceToObservation( data->observations(i)->position, relativ_edge->width, data->observations(i)  ) ;
        //untill 2.0 meters of distance, normal behavior. after that outliers behavior (not square)
        LossFunction* loss = NULL;
        loss = new ceres::ScaledLoss( g_param->useLoss?(new ceres::CauchyLoss(g_param->lossScale)):NULL
                                                       ,g_param->K_obs
                                      *data->observations(i)->confidence * data->observations(i)->weight /// @FIXEME @TEMP @TODO :
                                      ,ceres::DO_NOT_TAKE_OWNERSHIP) ;

        problem->AddResidualBlock(
                    distance_cost_function
                    ,loss
                    ,start_node->position
                    ,end_node->position
                    ,relativ_edge->width
                    ); //note : both observations are referring to these nodes.

        //saving the constraint for reuse to output at each step
        //creating the constraint
        Constraint * n_constraint = new Constraint_sidewalk(
                    start_node->position
                    ,end_node->position
                    ,start_node->position
                    ,relativ_edge->width
                    ,distance_cost_function
                    ,data->observations(i)->position) ;
        //adding it to list of constraints
        data->constraints()->push_back(n_constraint);
        count++;
    }
    return count ;
}




//! manual constraint based on surf distance between object and and edge
int addManualConstraintsOnSurfDistToObjects(DataStorage * data, Problem * problem){
    int count = 0;
    for (int i = 0; i < data->num_street_objects(); ++i) {

        //finding the street_object
        //finding the 2 nodes concerned by this observations

        street_object * obj = data->street_objects(i) ;

        edge * relativ_edge = data->ebe(obj->edge_id) ;
        node * start_node = data->nbn(relativ_edge->start_node)  ;
        node * end_node = data->nbn(relativ_edge->end_node)  ;

        CostFunction* distance_cost_function=
                new  ManualAttr_Rep_Object(i, data ) ;
        //untill 2.0 meters of distance, normal behavior. after that outliers behavior (not square)
        LossFunction* loss = NULL;
        loss = new ceres::ScaledLoss( g_param->useLoss?(new ceres::CauchyLoss(g_param->lossScale)):NULL
                                                       ,g_param->K_obj,ceres::DO_NOT_TAKE_OWNERSHIP) ;

        problem->AddResidualBlock(
                    distance_cost_function
                    ,loss
                    ,start_node->position
                    ,end_node->position
                    ,relativ_edge->width
                    ); //note : obj is referring to these nodes.

        //saving the constraint for reuse to output at each step
        double* obj_centroid_double = new double[3] ;
        obj_centroid_double[0]=0; obj_centroid_double[1]=0; obj_centroid_double[2]=0;
        Geometry::geomPoint2Double(obj->geom_centroid,obj_centroid_double) ;



        //creating the constraint
        Constraint * n_constraint = new Constraint_objects(
                    start_node->position
                    ,end_node->position
                    ,start_node->position//useless
                    ,relativ_edge->width
                    ,distance_cost_function
                    ,obj_centroid_double) ;
        //adding it to list of constraints
        data->constraints()->push_back(n_constraint);
        count++;

    }
    return count ;
}

/* adding constraint to meet the target slope*/
int addManualTargetSlope(DataStorage * data, Problem * problem){
	for (int i = 0; i < data->num_slopes(); ++i) {
		//for each target slope, adding a constraint
		slope*  slo = data->slopes(i) ;
		edge * relativ_edge = data->ebe(slo->edge_id) ;
        node * start_node = data->nbn(relativ_edge->start_node)  ;
        node * end_node = data->nbn(relativ_edge->end_node)  ;
		//K_slope
		CostFunction* distance_cost_function=
                new  Distance_target_slope(i, data) ;
				
		LossFunction* loss = NULL;
        loss = new ceres::ScaledLoss( g_param->useLoss?(new ceres::CauchyLoss(g_param->lossScale)):NULL
                                                       ,g_param->K_slope,ceres::DO_NOT_TAKE_OWNERSHIP) ;
		
		problem->AddResidualBlock(
                    distance_cost_function
                    ,loss
                    ,start_node->position
                    ,end_node->position 
                    );  
	}
};




//! manual constraint based on distance between observation and edges
int addManualConstraintsOnOrthDistToObservation_width(DataStorage * data, Problem * problem){
    int count = 0;
    for (int i = 0; i < data->num_observations(); ++i) {

        //finding the 2 nodes concerned by this observations
        edge * relativ_edge = data->ebe(data->observations(i)->edge_id) ;
        node * start_node = data->nbn(relativ_edge->start_node)  ;
        node * end_node = data->nbn(relativ_edge->end_node)  ;

        CostFunction* distance_cost_function=
                new  ManualOrthDistanceToObservation_width( data->observations(i)->position, relativ_edge->width, data->observations(i)  ) ;
        //untill 2.0 meters of distance, normal behavior. after that outliers behavior (not square)
        LossFunction* loss = NULL;
        loss = new ceres::ScaledLoss( g_param->useLoss?(new ceres::SoftLOneLoss(g_param->lossScale)):NULL
                                                       ,g_param->K_obs_width ,ceres::DO_NOT_TAKE_OWNERSHIP) ;

        problem->AddResidualBlock(
                    distance_cost_function
                    ,loss
                    ,start_node->position
                    ,end_node->position
                    ,relativ_edge->width
                    ); //note : both observations are referring to these nodes.

        //saving the constraint for reuse to output at each step
        //creating the constraint
        Constraint * n_constraint = new Constraint_sidewalk_width(
                    start_node->position
                    ,end_node->position
                    ,start_node->position //useless
                    ,relativ_edge->width
                    ,distance_cost_function
                    ,data->observations(i)->position) ;
        //adding it to list of constraints
        data->constraints()->push_back(n_constraint);
        count ++ ;
    }
    return count;
}




//! manual constraint based on surf distance between object and and edge
int addManualConstraintsOnSurfDistToObjects_width(DataStorage * data, Problem * problem){
    int count = 0;
    for (int i = 0; i < data->num_street_objects(); ++i) {

        //finding the street_object
        //finding the 2 nodes concerned by this observations
        street_object * obj = data->street_objects(i) ;

        edge * relativ_edge = data->ebe(obj->edge_id) ;
        node * start_node = data->nbn(relativ_edge->start_node)  ;
        node * end_node = data->nbn(relativ_edge->end_node)  ;

        CostFunction* distance_cost_function=
                new  ManualAttr_Rep_Object_width(i, data ) ;
        //untill 2.0 meters of distance, normal behavior. after that outliers behavior (not square)
        LossFunction* loss = NULL;
        loss = new ceres::ScaledLoss( g_param->useLoss?(new ceres::SoftLOneLoss(g_param->lossScale/5)):NULL
                                                       ,g_param->K_obj_width,ceres::DO_NOT_TAKE_OWNERSHIP) ;

        problem->AddResidualBlock(
                    distance_cost_function
                    ,loss
                    ,start_node->position
                    ,end_node->position
                    ,relativ_edge->width
                    ); //note : obj is referring to these nodes.

        //saving the constraint for reuse to output at each step
        double* obj_centroid_double = new double[3] ;
        obj_centroid_double[0]=0; obj_centroid_double[1]=0; obj_centroid_double[2]=0;
        Geometry::geomPoint2Double(obj->geom_centroid,obj_centroid_double) ;



        //creating the constraint
        Constraint * n_constraint = new Constraint_objects_width(
                    start_node->position
                    ,end_node->position
                    ,start_node->position//useless
                    ,relativ_edge->width
                    ,distance_cost_function
                    ,obj_centroid_double) ;
        //adding it to list of constraints
        data->constraints()->push_back(n_constraint);
        count ++;
    }
    return count ;
}

