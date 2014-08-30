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
    }

//! function to get an idea of what is in the node
string nodeToString(node * n){ 
	//#node_id::int;X::double;Yi_filename::double;Z::double;is_in_intersection::int
	 
	std::ostringstream nstring;
	nstring << "(node_id : " << n->node_id  
		<< " , position : (" << n->position[0] << "," << n->position[1] <<"," << n->position[2] 
		<< "), is_in_intersection : "<< n->is_in_intersection << ")";  
	//return nstring ;
	return nstring.str() ;
}

//! function to get an idea of what is in the edge
string edgeToString(edge * n){ 
	//#edge_id::int;start_node::int;end_node::int;width::double
	 
	std::ostringstream nstring;
	nstring << "(edge_id : " << n->edge_id  
		<< " , start_node : " << n->start_node << ", end_node : " << n->end_node <<", width : " << n->width <<")"; 
	return nstring.str() ;
}

//! function to get an idea of what is in the observation
string observationToString(observation * n){ 
	//#obs_id::int;X::double;Y::double;Z::double;confidence::double;weight::double
	 
	std::ostringstream nstring;
	nstring << "(obs_id : " << n->obs_id  
		<< " , position : (" << n->position[0] << "," << n->position[1] <<"," << n->position[2] 
		<< "), confidence : " << n->confidence << ", weight : " << n->weight <<")"; 
	return nstring.str() ;
}

//! unit test for node to string
bool TEST_nodeToString() { 
	node * n = new node();
	n->node_id = 1;
	n->is_in_intersection = 1;
	double d[3] = {0.0,0.0,0.0};
	n->position = d;
	std::cout << nodeToString(n).c_str() ; 
	return true;
}

//! unit test for edge to string
bool TEST_edgeToString() { 
	edge * n = new edge();
	n->edge_id  =1;
	n->start_node = 1 ;
	n->end_node = 2;
	n->width = 4.0  ; 
	std::cout << edgeToString(n).c_str() ; 
	return true;
}

//! unit test for observation to string
bool TEST_observationToString() { 
	observation * n = new observation();
	n->obs_id = 1;
	n->confidence = 1.0;
	n->weight = 1.0;
	double d[3] = {0.0,0.0,0.0};
	n->position = d;
	std::cout << observationToString(n).c_str() ; 
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
    std::cout << "num nodes : " << num_nodes() << "\n"
		<< "num edges : " << num_edges() << "\n"
        << "num observations : " << num_observations() << "\n";

    //reading the nodes commentary :
    fgets(line, sizeof line, i_fptr);

    //node_id::int;X::double;Yi_filename::double;Z::double;is_in_intersection::int
    for (int i = 0; i < num_nodes_; ++i) {
        fgets(line, sizeof line, i_fptr);
        node* t_node = new node();
        double * coor = new double[3];
        if (sscanf(line, "%d;%lG;%lG;%lG;%hd",&t_node->node_id,&coor[0],&coor[1],&coor[2],&t_node->is_in_intersection) != 5) {
            std::cerr<< "error when trying ot read the nodes, wrong format" ;
        } else {
			t_node->position = coor;
            nodes_[i] = t_node; 
            std::cout << "node readed : " << nodeToString(nodes_[i]).c_str() << " \n"  ; 
            
        }
    }
    //reading the edge data :
    //reading the edge commentary :
    fgets(line, sizeof line, i_fptr);

    //#edge_id::int;start_node::int;end_node::int;width::double 
    for (int i = 0; i < num_edges_; ++i) {
        fgets(line, sizeof line, i_fptr);
        edge* t_edge = new edge();
        if (sscanf(line, "%d;%d;%d;%lG",&t_edge->edge_id,&t_edge->start_node,&t_edge->end_node ,&t_edge->width) != 4) {
            std::cerr<< "error when trying ot read the edges, wrong format" ;
        } else {
            edges_[i] = t_edge; 
            std::cout << "edge readed : " << edgeToString(edges_[i]).c_str() << " \n" ; 
        }
    }
    
    
    //reading the observation commentary :
    fgets(line, sizeof line, i_fptr);
    
    //#obs_id::int;X::double;Y::double;Z::double;confidence::double;weight::double
    for (int i = 0; i < num_observations_ ; ++i) {
        fgets(line, sizeof line, i_fptr);
        observation* t_observation = new observation();
        double * coor = new double[3];
        if (sscanf(line, "%d;%lG;%lG;%lG;%lG;%lG",&t_observation->obs_id,&coor[0],&coor[1],&coor[2],&t_observation->confidence, &t_observation->weight) != 6) {
            std::cerr<< "error when trying ot read the observations, wrong format" ;
        } else {
			t_observation->position = coor;
            observations_[i] = t_observation; 
            std::cout << "observation readed : " << observationToString(observations_[i]).c_str() << " \n"  ; 
            
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
    
    std::cout << "\nbeggingn writting the rest \n";
    std::cout << edges_[0]->start_node <<"\n";
    std::cout << nodeToString(nodes_[0]) <<"\n";
    std::cout << edgeToString(edges_[0]) <<"\n";
    std::cout << observationToString(observations_[0]);
    std::cout << "\n" ;
    //writing data . Example of what we want : 
    //LINESTRINGZ(X1 Y1 Z1, X2 Y2 Z2);12.98;YYYY-MM-DD HH:MM:SS.ssssss;YYYY-MM-DD HH:MM:SS.ssssss
    
    
    for (int i = 0; i < num_edges_; ++i) { 
		double cost = 10.0 ; 
		//fprintf(o_fptr,"LINESTRINGZ(%lG %lG %lG, %lG %lG %lG);%lG;2014-08-30 00:00:00.%6d;2014-08-30 00:00:00.%6d"
			std::cout << nodes_[edges_[i]->start_node]->position[0] << "\n";
			std::cout << nodes_[edges_[i]->start_node]->position[1] << "\n";
			std::cout << nodes_[edges_[i]->start_node]->position[2] << "\n";
			std::cout << nodes_[edges_[i]->end_node]->position[0] << "\n";
			std::cout << nodes_[edges_[i]->end_node]->position[1] << "\n";
			std::cout << nodes_[edges_[i]->end_node]->position[2] << "\n";
			std::cout << cost << "\n";
			std::cout << iteration << "\n";
			std::cout << (iteration+1) << "\n";
			//);
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



