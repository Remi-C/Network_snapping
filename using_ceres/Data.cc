//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////
/**
  This file holds the Data.h class implementation

  */


#include "Data.h"
using std::endl;
using std::cout;
extern Parameter* g_param;

/** Constructor for data input/storage
  @param name of the file containing the input data : the format is very specific
  @param name of the file where to write the temporary result for each iteration
  */
DataStorage::DataStorage(const string i_filename, const string o_filename )
    : input_file_path_(i_filename), output_file_path_(o_filename) {
    //this->readData(i_filename);
    nodes_ = 0x0;
    edges_ = 0x0;
    observations_ = 0x0;
    classifications_ = 0x0;
	street_objects_ = 0x0;
	slopes_ = 0x0 ; 

    numConstraintsWritten_ = 0;

}
DataStorage::~DataStorage(){
    delete[] nodes_;
    delete[] edges_;
    delete[] observations_;
    delete[] classifications_;
    delete[] street_objects_;
	delete[] slopes_;
    /// @TODO loop on constraint vector, for surface_objects, free applciation point
}

//! unit test for node to string
bool TEST_nodeToString() {
    node * n = new node();
    n->node_id = 1;
    n->is_in_intersection = 1;
    double d[3] = {0.0,0.0,0.0};
    std::copy(n->position, n->position + 3, d);
    std::cout << n->nodeToString().c_str() ;
    return true;
}

//! unit test for edge to string
bool TEST_edgeToString() {
    edge * n = new edge();
    n->edge_id  =1;
    n->start_node = 1 ;
    n->end_node = 2;
    n->width[0] = 4.0  ;
    std::cout << n->edgeToString().c_str() ;
    return true;
}

//! unit test for observation to string
bool TEST_observationToString() {
    observation * n = new observation();
    n->obs_id = 1;
    n->confidence = 1.0;
    n->weight = 1.0;
    double d[3] = {0.0,0.0,0.0};
    std::copy(n->position, n->position + 3, d);
    std::cout << n->observationToString().c_str() ;
    return true;
}





void DataStorage::readClassifications(){
    //opening files //cout << "here is the file path read in parameters.txt : "  << g_param->class_definition_path <<endl;
    FILE* i_fptr = fopen( g_param->class_definition_path.c_str() , "r");
    char line[1000];//input line buffer

    if (i_fptr == NULL) {
        std::cerr << "Error: unable to open file " << g_param->class_definition_path <<endl;
        return;
    }

    int class_id = 0;
    char class_name[1000] = "#######"; char geom_type[1000] ="########"; char road_surface_relation[1000] = "#####";
    double precision = 0;  double importance = 0; double dist_to_border = 0;


    int class_number = 0;
    //finding the number of class : parsing the file a first time :

    // char line[1000] = "34;wire;LINESTRING;OUT;0.2;1;0.1\n" ;

    while (!feof(i_fptr)) {
        fgets(line, sizeof line, i_fptr);
        if (
                sscanf( line , "%d;%[^;];%[^;];%[^;];%lG;%lG;%lG"
                        ,&class_id
                        ,class_name ,geom_type ,road_surface_relation
                        ,&precision ,&importance ,&dist_to_border
                        )
                != 7) {
            //reading a comment line
            //std::cout << "reading comment line : " << line ;
        } else {
            class_number++ ;
        }
    }
    //allocating memory for a new classification
    classifications_ = new classification[class_number];  //! @TODO : do alligned memory allocation it's better

    //new parsing
    class_number = 0 ;
    rewind(i_fptr);

    while (!feof(i_fptr)) {
        fgets(line, sizeof line, i_fptr);
        if (
                sscanf( line , "%d;%[^;];%[^;];%[^;];%lG;%lG;%lG"
                        ,&class_id
                        ,class_name ,geom_type ,road_surface_relation
                        ,&precision ,&importance ,&dist_to_border
                        )
                != 7) {
            //reading a comment line
            //std::cout << "reading comment line : " << line ;
        } else {
            //cout << "filling classification  " << class_number <<endl;
            //filling the classification :
            classifications_[class_number].setClassification(
                        class_id
                        ,string(class_name)
                        ,string(geom_type)
                        ,string(road_surface_relation)
                        ,precision
                        ,importance
                        ,dist_to_border
                        );
            //             cout <<class_id << endl
            //                << "\t"<<  string(class_name) << endl
            //                << "\t"<< string(geom_type) << " " <<   endl
            //                  << "\t"<< string(road_surface_relation)  <<endl
            //                << "\t"<< precision << endl
            //                 << "\t"<< importance <<endl
            //                << "\t"<< dist_to_border <<endl ;

            //             cout << classifications_[class_number].classificationToString()<<endl ;
            class_number++ ;

        }
    }

    //write class_number into parameters
    this->num_classifications_ = class_number ;
    fclose(i_fptr);
    return;
}

void DataStorage::readObjects(){
    //opening files //cout << "here is the file path read in parameters.txt : "  << g_param->class_definition_path <<endl;
    FILE* i_fptr = fopen( g_param->objects_path.c_str() , "r");
    char line[100000];//input line buffer

    if (i_fptr == NULL) {
        std::cerr << "Error: unable to open file " << g_param->objects_path <<endl;
        return;
    }
    //object_id;class_id;class_name;edge_id;geom;confidence;
    int obj_number = 0;
    int obj_id= 0;
    int class_id = 0;
    char class_name[1000] = "#########";
    int edge_id = 0;
    char geom_wkt[100000] ="#########";
    double confidence = 0;

    //char line[1000] = "1;75;car;1;POLYGON((0 0, 1 0, 1 1 , 0 1 , 0 0 ));0.89";

    while (fgets(line, sizeof line, i_fptr) != NULL ) {

        //cout << "readed line: " << line ;
        if (
                sscanf( line , "%d;%d;%[^;];%d;%[^;];%lG"
                        ,&obj_id
                        ,&class_id
                        ,class_name
                        ,&edge_id
                        ,geom_wkt
                        ,&confidence
                        )
                != 6) { //reading a comment line
            //std::cout << "reading comment line : " << line ;
        } else {
            obj_number++ ;
        }
    }
    //allocating memory for a new classification
    street_objects_ = new street_object[obj_number];  //! @TODO : do alligned memory allocation it's better
    obj_number = 0 ;
    rewind(i_fptr);
    while (fgets(line, sizeof line, i_fptr) != NULL ) {

        //        cout << "readed line: " << line ;
        if (
                sscanf( line , "%d;%d;%[^;];%d;%[^;];%lG"
                        ,&obj_id
                        ,&class_id
                        ,class_name
                        ,&edge_id
                        ,geom_wkt
                        ,&confidence
                        )
                != 6) {
            //reading a comment line
            //std::cout << "reading comment line : " << line ;
        } else {

            //            cout <<obj_id << endl
            //                << "\t"<<  class_id << endl
            //                << "\t"<< string(class_name) << " " <<   endl
            //                << "\t"<< edge_id  <<endl
            //                << "\t"<< geom_wkt << endl
            //                << "\t"<< confidence <<endl ;

            classification * object_classification = this->cbn(class_name);

            //cout << "width :" << object_classification->dist_to_border <<endl ;
            street_objects_[obj_number].setStreet_Object(
                        obj_id
                        ,class_id
                        ,class_name
                        ,edge_id
                        ,geom_wkt
                        , object_classification->dist_to_border
                        ,confidence
                        );

            obj_number++ ;
        }
    }


    //write class_number into parameters
    this->num_street_objects_= obj_number ;
    fclose(i_fptr);
    return;
}

void DataStorage::readData(){
    //opening files
    FILE* i_fptr = fopen(input_file_path_.c_str(), "r");
    char line[1000];//input line buffer

    if (i_fptr == NULL) {
        std::cerr << "Error: unable to open file " << input_file_path_;
        return;
    }

    //skipping the first line of comment
    fgets(line, sizeof line, i_fptr);
    //std::cout << "header1 : "<< line ;
    //reading data parameter : how much we have to read after this.
    fgets(line, sizeof line, i_fptr);
    if (	sscanf(line, "%d;%d;%d;%d",&num_nodes_,&num_edges_,&num_observations_, &num_slopes_) != 4) {
        std::cerr<< "error when trying to read the numbero of nodes, number of edges, number of observations, number of target_slopes wrong format" ;
    }

    //some debug output
    //std::cout << "num nodes : " << num_nodes() << "," << "num edges : " << num_edges() << ","<< "num observations : " << num_observations() << "\n";

    //reading node data:
    //reading the nodes commentary :
    fgets(line, sizeof line, i_fptr);
    //std::cout << "header2 : "<< line ;
    //allocating nodes memory
    nodes_ = new node[num_nodes_];  //! @TODO : do alligned memory allocation it's better

    //we try to write the following line format :
    //node_id::int;X::double;Y::double;Z::double;is_in_intersection::int
    for (int i = 0; i < num_nodes_; ++i) {
        fgets(line, sizeof line, i_fptr);
        if (sscanf(line, "%d;%lG;%lG;%lG;%hd",&nodes_[i].node_id,&nodes_[i].position[0],&nodes_[i].position[1],&nodes_[i].position[2],&nodes_[i].is_in_intersection) != 5) {
            std::cerr<< "error when trying ot read the nodes, wrong format\n" ;
        } else {

            //std::copy(t_node.position, t_node.position + 3, coor);
            //nodes_[i] = &t_node;
            //std::cout << "node readed : " << nodes_[i].nodeToString().c_str() << " \n"  ;

        }
    }
    //reading the edge data :
    //reading the edge commentary :
    fgets(line, sizeof line, i_fptr);
    //std::cout << "header3 : "<< line ;

    //allocating the memory for edges%hd
    edges_ = new edge[num_edges_];  //! @TODO : do alligned memory allocation it's better
    //#edge_id::int;start_node::int;end_node::int;width::double
    for (int i = 0; i < num_edges_; ++i) {
        fgets(line, sizeof line, i_fptr);
        if (sscanf(line, "%d;%d;%d;%lG",&edges_[i].edge_id,&edges_[i].start_node,&edges_[i].end_node ,edges_[i].width) != 4) {
            std::cerr<< "error when trying ot read the edges, wrong format\n" ;
        } else {
            //std::cout << "edge readed : " << edges_[i].edgeToString().c_str() << " \n" ;
        }
    }

    //reading the commentary for observation
    fgets(line, sizeof line, i_fptr);
    //std::cout << "header4 : " <<line ;

    //alocating hte memory for observations
    observations_ = new observation[num_observations()];//! @TOOO

    //#obs_id::int;edge_id::int;X::double;Y::double;Z::double;confidence::double;weight::double
    for (int i = 0; i < num_observations_ ; ++i) {
        fgets(line, sizeof line, i_fptr);
        //std::cout << line  << "\n";
        if (sscanf(line, "%d;%d;%lG;%lG;%lG;%lG;%lG"
                   ,&observations_[i].obs_id
                   ,&observations_[i].edge_id
                   ,&observations_[i].position[0],&observations_[i].position[1],&observations_[i].position[2]
                   ,&observations_[i].confidence, &observations_[i].weight) != 7) {
            std::cerr<< "error when trying to read the observations, wrong format\n" ;
            std::cerr<< "was expecting : obs_id::int;edge_id::int;X::double;Y::double;Z::double;confidence::double;weight::double, \n recevied : " << line ;
        } else {
            //std::copy(t_observation->position, t_observation->position + 3, coor);
            //std::cout << "observation readed : " << observations_[i].observationToString().c_str() << " \n"  ;

        }
    }
	
    //reading the slope commentary :
    fgets(line, sizeof line, i_fptr);

	slopes_ = new slope[num_slopes()];//! @TOOO
	//#edge_id::int;slope::double;confidence::double;weight::double
    for (int i = 0; i < num_slopes_ ; ++i) {
        fgets(line, sizeof line, i_fptr);
        //std::cout << line  << "\n";
        if (sscanf(line, "%d;%lG;%lG;%lG"
                   ,&slopes_[i].edge_id
                   ,&slopes_[i].slope[0]
                   ,&slopes_[i].confidence, &slopes_[i].weight) != 4) {
            std::cerr<< "error when trying to read the target slopes, wrong format\n" ;
            std::cerr<< "was expecting : edge_id::int;slope::double;confidence::double;weight::double \n recevied : " << line ;
        } else {
            //std::copy(t_observation->position, t_observation->position + 3, coor);
            //std::cout << "observation readed : " << observations_[i].observationToString().c_str() << " \n"  ;

        }
    }
	
	

    std::cout << num_nodes() <<" nodes, " << num_edges() << " edges, " << num_observations() << " observations "
    << num_slopes()  << " target slopes read, " <<"  from file " << input_file_path_ << " \n" ;
    fclose(i_fptr);
    return;
}



void DataStorage::writeData(int iteration){
    //opening files
    FILE* o_fptr ;

    if(iteration == 1){ //we clean the file at first writting, after we append
        o_fptr = fopen(output_file_path_.c_str(), "w");
    }
    else {
        o_fptr = fopen(output_file_path_.c_str(), "a+");
    }
    if (o_fptr == NULL) {
        std::cerr << "Error: unable to open file " << output_file_path_;
        return;
    }

    //writing the header if needed :
    if(iteration == 1){ // no need to write the header each time
        fprintf(o_fptr, "#geom;width;cost;start_time;end_time;iteration\n") ;
    }

    //writing data . Example of what we want :
    //LINESTRINGZ(X1 Y1 Z1, X2 Y2 Z2);12.98;YYYY-MM-DD HH:MM:SS.ssssss;YYYY-MM-DD HH:MM:SS.ssssss

    //this allows to represent 60*60*60 iterations
    //we compute the H:M:S based on number of iteration.
    //Note : we rely on cast to int : int(A/B) is going to return rest of euclidian div.
    int hours_c = int( iteration/(60*60) ) ;
    int minutes_c = int( iteration/(60) ) ;
    int seconds_c = iteration%60 ;
    //same based on iteration+1, for next time.
    int hours_n = int( (iteration+1)/(60*60) ) ;
    int minutes_n = int( (iteration+1)/(60) ) ;
    int seconds_n = (iteration+1)%60 ;

    for(const auto& element : this->edges_by_edge_id_){
        //std::cout << element.second->end_node << std::endl;
        edge * edge_to_output = element.second;
        node * start_node = nbn(edge_to_output->start_node) ;
        node * end_node = nbn(edge_to_output->end_node) ;
        double cost = 10.0 ;
        fprintf(o_fptr,"%d;LINESTRINGZ(%.4lf %.4lf %.4lf, %.4lf %.4lf %.4lf);%lG;%lG;2014-08-30 %02d:%02d:%02d;2014-08-30 %02d:%02d:%02d;%d\n"
                , edge_to_output->edge_id
                , start_node->position[0]
                , start_node->position[1]
                , start_node->position[2]
                , end_node->position[0]
                , end_node->position[1]
                , end_node->position[2]
                , edge_to_output->width[0]
                , cost
                , hours_c
                , minutes_c
                , seconds_c
                , hours_n
                , minutes_n
                , seconds_n
                , iteration
                );
    }
    std::cout << num_edges() <<" edges written to file : " << output_file_path_ << "\n";
    fclose(o_fptr);
    return;
}

/** fill the hash_table that allow to reference node by node_id rather than by
  index in nodes_
  */
void DataStorage::setMap(){

    //setting the map to find node by node_id
    for (int i = 0; i < num_nodes_; ++i) {
        nodes_by_node_id_[nodes_[i].node_id] = &nodes_[i] ;
    }//loop on all node from nodes_

    //seting the map to find edge by edge_id
    for (int i = 0; i < num_edges(); ++i) {
        edges_by_edge_id_[edges_[i].edge_id] = &edges_[i] ;
    }//loop on all edge from edges_

    //setting the map to find edege by node_id (each edge has 2 node id)
    for (int i = 0; i < num_edges(); ++i) {
        edges_by_node_id_.insert(std::pair<int,edge *>(edges_[i].start_node,&edges_[i] ));
        edges_by_node_id_.insert(std::pair<int,edge *>(edges_[i].end_node,&edges_[i] ));
    }//loop on all edge from edges_

    //setting the map to find classification by class id or by class_name
    for (int i = 0; i < num_classifications(); ++i) {
        classification_by_id_.insert(std::pair<int,classification*>(classifications_[i].class_id,&classifications_[i] ));
        classification_by_name_.insert(std::pair<string,classification*>(classifications_[i].class_name,&classifications_[i] ));
    }//loop on all classifications

    return;
}



void DataStorage::writeConstraints(int iteration){
    //opening files
    FILE* o_fptr ;
    string output_file_path_ = "../data/data_in_reduced_export_area/snapping_constraints.csv" ; /// @todo : put a parameter here

    if(iteration == 1){ //we clean the file at first writting, after we append
        o_fptr = fopen(output_file_path_.c_str(), "w");
    }
    else {
        o_fptr = fopen(output_file_path_.c_str(), "a+");
    }
    if (o_fptr == NULL) {
        std::cerr << "Error: unable to open file " << output_file_path_;
        return;
    }

    //writing the header if needed :
    if(iteration == 1){ // no need to write the header each time
        fprintf(o_fptr, "#gid;edge_id;type;cost;geom;start_time;end_time;iteration\n") ;
    }

    //writing data . Example of what we want :
    //1;1;SURF_DIST;3.76;LINESTRINGZ(X1 Y1 Z1, X2 Y2 Z2);YYYY-MM-DD HH:MM:SS.ssssss;YYYY-MM-DD HH:MM:SS.ssssss;1
    int hours_c = int( iteration/(60*60) ) ;
    int minutes_c = int( iteration/(60) ) ;
    int seconds_c = iteration%60 ;
    //same based on iteration+1, for next time.
    int hours_n = int( (iteration+1)/(60*60) ) ;
    int minutes_n = int( (iteration+1)/(60) ) ;
    int seconds_n = (iteration+1)%60 ;

    //loop on all constraints to be outputted
    int old_gid = numConstraintsWritten();
    int gid = numConstraintsWritten()+1;

    double t_c[2] = {0,0} ;
    double * t_g = new double[6];

    for(auto& constraint : constraints_){

        //allocate memory for cost and geom

        //get the cost and geom associated to this constraint for this step
        constraint->get_graphical_constraint(t_c, t_g );
        //getting the name of the class of constraint :
        int status = 0;
        char* constraint_name = abi::__cxa_demangle(
                    typeid( *constraint ).name()
                                ,NULL
                                ,0
                                ,&status) ;

        //this allows to represent 60*60*60 iterations
        //we compute the H:M:S based on number of iteration.
        //Note : we rely on cast to int : int(A/B) is going to return rest of euclidian div.

        fprintf(o_fptr,"%d;edge_id;%s;%lG;LINESTRINGZ(%.4lf %.4lf %.4lf, %.4lf %.4lf %.4lf);2014-08-30 %02d:%02d:%02d;2014-08-30 %02d:%02d:%02d;%d\n"
                , gid
                , constraint_name
                , t_c[0]
                , t_g[0], t_g[1], t_g[2], t_g[3], t_g[4], t_g[5]
                , hours_c
                , minutes_c
                , seconds_c
                , hours_n
                , minutes_n
                , seconds_n
                , iteration
                );
        gid ++;
        numConstraintsWritten_++;
        // for
        //    for(const auto& element : this->edges_by_edge_id_){
        //        //std::cout << element.second->end_node << std::endl;
        //        edge * edge_to_output = element.second;
        //        node * start_node = nbn(edge_to_output->start_node) ;
        //        node * end_node = nbn(edge_to_output->end_node) ;
        //        double cost = 10.0 ;
        //        fprintf(o_fptr,"LINESTRINGZ(%lG %lG %lG, %lG %lG %lG);%lG;%lG;2014-08-30 %02d:%02d:%02d;2014-08-30 %02d:%02d:%02d;%d\n"
        //                , start_node->position[0]
        //                , start_node->position[1]
        //                , start_node->position[2]
        //                , end_node->position[0]
        //                , end_node->position[1]
        //                , end_node->position[2]
        //                , edge_to_output->width
        //                , cost
        //                , hours_c
        //                , minutes_c
        //                , seconds_c
        //                , hours_n
        //                , minutes_n
        //                , seconds_n
        //                , iteration
        //                );
        //    }
    }
    std::cout << (gid - old_gid)  <<" constraints written to file : " << output_file_path_ << "\n";
    fclose(o_fptr);
    return;
}
