#ifndef DATA_H
#define DATA_H
//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////
/**
  This file holds the class that deals with data input and storage during optimisation


  format of the file to be read
""""""""""""""""""""""""""""""""""""" 
  #header
  num_nodes
  num_edges
  num_observations
  #header node
  node_id::int;X::double;Yi_filename::double;Z::double;is_in_intersection::int
  ...
  #header edge
  edge_id::int;start_node::int;end_node::int;width::double
  ...
  #header observations
  obs_id::int;X::double;Y::double;Z::double;confidence::double;weight::double
"""""""""""""""""""""""""""""""""""""  

  format of the file to be written
""""""""""""""""""""""""""""""""""""" 
	#geom;cost;start_time;end_time
	LINESTRINGZ(X1 Y1 Z1, X2 Y2 Z2);12.98;YYYY-MM-DD HH:MM:SS.ssssss;YYYY-MM-DD HH:MM:SS.ssssss
	LINESTRINGZ(X1 Y1 Z1, X2 Y2 Z2);12.98;YYYY-MM-DD HH:MM:SS.ssssss;YYYY-MM-DD HH:MM:SS.ssssss
	...
"""""""""""""""""""""""""""""""""""""  
  
  */
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map> //the hash_table


using std::string;
using std::basic_string;


//! a simple structure to hold observation data : that is a point with some attributes
struct node{
  int node_id;            //! unique id per node
  double position[3];      //! 3 coordinates X,Y,Z. Will be changed upon optimization
  short is_in_intersection;//! is this node part of an intersection, or is this an intermediary node?
  
  
  node(int id,double X,double Y,double Z,bool in_inter) {
		node_id = id;
        position[0] = X ;
		position[1] = Y ;
		position[2] = Z;
		is_in_intersection = in_inter ;
  }

  node() {
        node_id = 0;
        position[0] = 0.0 ;
        position[1] = 0.0 ;
        position[2] = 0.0 ;
        is_in_intersection = false ;
  }

  string nodeToString(){ 
	//#node_id::int;X::double;Yi_filename::double;Z::double;is_in_intersection::int
	 
	std::ostringstream nstring;
	nstring << "(node_id : " << node_id  
		<< " , position : (" << position[0] << "," << position[1] <<"," << position[2] 
		<< "), is_in_intersection : "<< is_in_intersection << ")";  
	//return nstring ;
	return nstring.str() ;
}

};



struct edge{
  int edge_id;      //! unique id per edge
  int start_node;   //! link to the id of the node that starts this edge
  int end_node;     //! link to the id of the node that ends this edge
  double width;     //! width of the edge, in meters.
  
  //! function to get an idea of what is in the edge
string edgeToString(){ 
	//#edge_id::int;start_node::int;end_node::int;width::double
	 
	std::ostringstream nstring;
	nstring << "(edge_id : " << edge_id  
		<< " , start_node : " << start_node << ", end_node : " << end_node <<", width : " << width <<")"; 
	return nstring.str() ;
}

};

struct observation{
  int obs_id;           //! unique id per observations
  double position[3];    //! 3 coordinates X,Y,Z
  //int edge_id;          //! this edge should be snapped to this observation
  double confidence;    //! confidence between 0 and 1. 0-> very unlikely ; 1-> certain
  double weight;        //! statistical weight of this observation (observation are points extracted from lines, weight = length of the 2 lines/2)
  
  //! function to get an idea of what is in the observation
string observationToString(){ 
	//#obs_id::int;X::double;Y::double;Z::double;confidence::double;weight::double
	 
	std::ostringstream nstring;
	nstring << "(obs_id : " << obs_id  
		<< " , position : (" << position[0] << "," << position[1] <<"," << position[2] 
		<< "), confidence : " << confidence << ", weight : " << weight <<")"; 
	return nstring.str() ;
}

};
 

class DataStorage {
public:
    DataStorage(const string, const string );

    ~DataStorage();
	void readData();
	void writeData(int);
    void setMap();

    int num_nodes()             { return num_nodes_;  }
    int num_edges()             { return num_edges_;    }
    int num_observations()       { return num_observations_;  }
    node* nodes()                { return nodes_;}
    node* nodes(int i)           { return &nodes_[i];}
    edge* edges()                { return edges_;}
    edge* edges(int i)                { return &edges_[i];}
    observation* observations()  { return observations_;}
    observation* observations(int i)  { return &observations_[i];}
    node* nbn(int i ) { return nodes_by_node_id_.at(i); }
    std::unordered_map <int /*node_id*/, node *> nodes_by_node_id() {return nodes_by_node_id_;}

private:
    int num_nodes_; //! total num of nodes we are going to read
    int num_edges_; //! total num of edges we are going to read
    int num_observations_;//! total num of observations we are going to read

    node* nodes_;//! an array of node
    edge* edges_;//! an array og edges
    observation* observations_;//! an array of observations
    
    const string input_file_path_;//! name of the file containing the input data
    const string output_file_path_;//! name of the file where to write results

    //! a hads table
    std::unordered_map <int /*node_id*/, node *> nodes_by_node_id_;

};

#endif // DATA_H
