//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////
/**
  This file holds the Data.h class implementation

  */

#include "Data.h"

///////////////// 	  PARAMETER     ////////////////////////////
const string input_file_path("../../data/simple_example.csv") ;
const string output_file_path("../../data/optimisation_output.csv") ;
////////////////////////////////////////////////////////////////

/** Constructor for data input/storage
  @param name of the file containing the input data : the format is very specific
  @param name of the file where to write the temporary result for each iteration
  */
DataStorage::DataStorage(const string i_filename, const string o_filename ) : input_file_path_(i_filename), output_file_path_(o_filename) {
    //this->readData(i_filename);
    nodes_ = 0x0;
    edges_ = 0x0;
    observations_ = 0x0;
    }
DataStorage::~DataStorage(){
	delete[] nodes_;
	delete[] edges_;
	delete[] observations_;
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
	n->width = 4.0  ; 
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
    
    //reading data parameter : how much we have to read after this.
    fgets(line, sizeof line, i_fptr);
    if (	sscanf(line, "%d;%d;%d",&num_nodes_,&num_edges_,&num_observations_) != 3) {
            std::cerr<< "error when trying to read the numbero of nodes, number of edges, number of observations, wrong format" ;
        }
    
    //some debug output
    std::cout << "num nodes : " << num_nodes() << ","
		<< "num edges : " << num_edges() << ","
        << "num observations : " << num_observations() << "\n";

	//reading node data: 
    //reading the nodes commentary :
    fgets(line, sizeof line, i_fptr);
	
	//allocating nodes memory
    nodes_ = new node[num_nodes_];  //! @TODO : do alligned memory allocation it's better
	
	//we try to write the following line format : 
    //node_id::int;X::double;Yi_filename::double;Z::double;is_in_intersection::int
    for (int i = 0; i < num_nodes_; ++i) {  
        if (fscanf(i_fptr, "%d;%lG;%lG;%lG;%hd",&nodes_[i].node_id,&nodes_[i].position[0],&nodes_[i].position[1],&nodes_[i].position[2],&nodes_[i].is_in_intersection) != 5) {
            std::cerr<< "error when trying ot read the nodes, wrong format\n" ;
        } else {
			 
			//std::copy(t_node.position, t_node.position + 3, coor); 
            //nodes_[i] = &t_node; 
            std::cout << "node readed : " << nodes_[i].nodeToString().c_str() << " \n"  ;
            
        }
    }
    //reading the edge data :
    //reading the edge commentary :
    fgets(line, sizeof line, i_fptr);
    fgets(line, sizeof line, i_fptr);
    std::cout << line ;

    //allocating the memory for edges
    edges_ = new edge[num_edges_];  //! @TODO : do alligned memory allocation it's better
    //#edge_id::int;start_node::int;end_node::int;width::double 
    for (int i = 0; i < num_edges_; ++i) {

        if (fscanf(i_fptr, "%d;%d;%d;%lG",&edges_[i].edge_id,&edges_[i].start_node,&edges_[i].end_node ,&edges_[i].width) != 4) {
            std::cerr<< "error when trying ot read the edges, wrong format\n" ;
        } else {
            std::cout << "edge readed : " << edges_[i].edgeToString().c_str() << " \n" ;
        }
    }
    
    //reading the commentary for observation
    fgets(line, sizeof line, i_fptr);
    fgets(line, sizeof line, i_fptr);
    std::cout << line ;

    //alocating hte memory for observations
    observations_ = new observation[num_observations()];//! @TOOO

    //#obs_id::int;X::double;Y::double;Z::double;confidence::double;weight::double
    for (int i = 0; i < num_observations_ ; ++i) {
        
        if (fscanf(i_fptr, "%d;%lG;%lG;%lG;%lG;%lG",&observations_[i].obs_id
                ,&observations_[i].position[0],&observations_[i].position[1],&observations_[i].position[2]
                ,&observations_[i].confidence, &observations_[i].weight) != 6) {
            std::cerr<< "error when trying to read the observations, wrong format\n" ;
        } else {
            //std::copy(t_observation->position, t_observation->position + 3, coor);
            std::cout << "observation readed : " << observations_[i].observationToString().c_str() << " \n"  ;
            
        }
    } 
    

    
    std::cout << fclose(i_fptr); 
    return;
 }



void DataStorage::writeData(int iteration){
    //opening files
    FILE* o_fptr = fopen(output_file_path_.c_str(), "w");

    if (o_fptr == NULL) {
      std::cerr << "Error: unable to open file " << output_file_path_;
      return;
    }
	
	//writing the header :
    fprintf(o_fptr, "#geom;cost;start_time;end_time\n") ;
    
    std::cout << "#geom;cost;start_time;end_time\n" ;
    std::cout << "\nbeggingn writting the rest \n";

    //writing data . Example of what we want : 
    //LINESTRINGZ(X1 Y1 Z1, X2 Y2 Z2);12.98;YYYY-MM-DD HH:MM:SS.ssssss;YYYY-MM-DD HH:MM:SS.ssssss
    
    
    for (int i = 0; i < num_edges_; ++i) { 
		double cost = 10.0 ; 
        fprintf(o_fptr,"LINESTRINGZ(%lG %lG %lG, %lG %lG %lG);%lG;2014-08-30 00:00:00.%06d;2014-08-30 00:00:00.%06d"
            , nodes_[edges_[i].start_node-1].position[0]
            , nodes_[edges_[i].start_node-1].position[1]
            , nodes_[edges_[i].start_node-1].position[2]
            , nodes_[edges_[i].end_node-1].position[0]
            , nodes_[edges_[i].end_node-1].position[1]
            , nodes_[edges_[i].end_node-1].position[2]
            , cost
            , iteration
            , (iteration+1)
            );
    }
    
    std::cout << fclose(o_fptr); 
    return;
 }



int main(int argc, char** argv) { 
	
	//getting the data ;
	DataStorage * data = new DataStorage(input_file_path,output_file_path) ;
	data->readData();
	data->writeData(1);
  return 0;
}



