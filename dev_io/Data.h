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
  num_nodes
  num_edges
  num_observations
  node_id::int;X::double;Yi_filename::double;Z::double;is_in_intersection::int
  ...
  edge_id::int;start_node::int;end_node::int;width::double
  ...
  obs_id::int;X::double;Y::double;Z::double;confidence::double;weight::double
  */
#include <string>
#include <iostream>
#include <fstream>

using std::string;
using std::basic_string;


//! a simple structure to hold observation data : that is a point with some attributes
struct node{
  int node_id;            //! unique id per node
  double * position;      //! 3 coordinates X,Y,Z. Will be changed upon optimization
  short is_in_intersection;//! is this node part of an intersection, or is this an intermediary node?
};



struct edge{
  int edge_id;      //! unique id per edge
  int start_node;   //! link to the id of the node that starts this edge
  int end_node;     //! link to the id of the node that ends this edge
  double width;     //! width of the edge, in meters.
};

struct observation{
  int obs_id;           //! unique id per observations
  double * position;    //! 3 coordinates X,Y,Z
  double confidence;    //! confidence between 0 and 1. 0-> very unlikely ; 1-> certain
  double weight;        //! statistical weight of this observation (observation are points extracted from lines, weight = length of the 2 lines/2)
};
 

class DataStorage {
public:
    DataStorage(const  string, const  string );

    ~DataStorage();
	void readData(const string);
    void WriteToFile(const string filename) const;

    int num_nodes()             { return num_nodes_;  }
    int num_edges()             { return num_edges_;    }
    int num_observations()       { return num_observations_;  }
    node** nodes()                { return nodes_;}
    edge** edges()                { return edges_;}
    observation** observations()  { return observations_;}


private:
    int num_nodes_;
    int num_edges_;
    int num_observations_;

    node* nodes_[];
    edge* edges_[];
    observation* observations_[];
};

#endif // DATA_H
