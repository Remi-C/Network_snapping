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
const string output_file_path("/media/sf_E_RemiCura/PROJETS/snapping/visu/visu_inqgis_timemanager/simple_test.csv") ;
////////////////////////////////////////////////////////////////

/** Constructor for data input/storage
  @param name of the file containing the input data : the format is very specific
  @param name of the file where to write the temporary result for each iteration
  */
DataStorage::DataStorage(const string i_filename, const string o_filename ) {
    this->readData(i_filename);
    }

//! function to get an idea of what is in the node
string nodeToString(node * n){
	string nstring;
	//#node_id::int;X::double;Yi_filename::double;Z::double;is_in_intersection::int
	nstring = '(node_id : ' + n->node_id 
		+ ' , (' + n->position[0] + ',' + n->position[1] + ',' + n->position[2] 
		+ '), is_in_intersection : ' + n->is_in_intersection + ' )';  
	return nstring ;
}

void DataStorage::readData(const  string i_filename){
    //opening files
    FILE* i_fptr = fopen(i_filename.c_str(), "r");
 
	char line[1000];//input line buffer 
	
    if (i_fptr == NULL) {
      std::cerr << "Error: unable to open file " << i_filename;
      return;
    };

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
            std::cout << " toto \n";
            std::cout << "data readed : " << nodeToString(nodes_[i]).c_str() << " \n"  ; 
            
        }
    }
    //reading the edge data :
    fclose(i_fptr);
 }




int main(int argc, char** argv) { 
	
	//test of node output
	
	//creating a node
	node * n = new node();
	n->node_id = 1;
	n->is_in_intersection = 1;
//	double * d = new double[3] ;
	double d[3] = {0.0,0.0,0.0};
	n->position = d;
	
	std::cout << nodeToString(n);
	
	//getting the data ;
	//DataStorage * data = new DataStorage(input_file_path,output_file_path) ;
 
  return 0;
}


