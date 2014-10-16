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
  num_nodes;num_edges;num_observations
  #header node
  node_id::int;X::double;Yi_filename::double;Z::double;is_in_intersection::int
  ...
  #header edge
  edge_id::int;start_node::int;end_node::int;width::double
  ...
  #header observations
  obs_id::int;edge_id::int;X::double;Y::double;Z::double;confidence::double;weight::double
"""""""""""""""""""""""""""""""""""""  

  format of the file to be written
""""""""""""""""""""""""""""""""""""" 
    #geom;width;cost;start_time;end_time;iteration
    LINESTRINGZ(X1 Y1 Z1, X2 Y2 Z2);2.25;12.98;YYYY-MM-DD HH:MM:SS.ssssss;YYYY-MM-DD HH:MM:SS.ssssss;1
    LINESTRINGZ(X1 Y1 Z1, X2 Y2 Z2);2.25;12.98;YYYY-MM-DD HH:MM:SS.ssssss;YYYY-MM-DD HH:MM:SS.ssssss;2
    ...
"""""""""""""""""""""""""""""""""""""  

  */
#include "Parameters.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map> //the hash_table
#include <math.h> //for

#include "Parameters.h"
#include "enum_functions.h"
#include "geometry_function.h"

using std::string;
using std::basic_string;
//using SnapEnums;



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
    int edge_id;          //! id of the edge that shoud be snapped to this observation
    double position[3];    //! 3 coordinates X,Y,Z
    //int edge_id;          //! this edge should be snapped to this observation
    double confidence;    //! confidence between 0 and 1. 0-> very unlikely ; 1-> certain
    double weight;        //! statistical weight of this observation (observation are points extracted from lines, weight = length of the 2 lines/2)

    //! function to get an idea of what is in the observation
    string observationToString(){
        //#obs_id::int;X::double;Y::double;Z::double;confidence::double;weight::double

        std::ostringstream nstring;
        //nstring.precision(10);
        nstring << "(obs_id : " << obs_id  << " , edge_id : (" << edge_id
                << " , position : (" << position[0] << "," << position[1] <<"," << position[2]
                << "), confidence : " << confidence << ", weight : " << weight <<")";
        return nstring.str() ;
    }


};









struct classification{

    unsigned char class_id;           //! unique id per class, positive, 8 bits long, contains kind of information on object type grouping
    std::string class_name;          //! name of the class, in english
    SnapEnums::geom_type_enum geom_type;   //! geometry type, folowwing geos  : 1 = point, 2 = line, 3 = polygon
    SnapEnums::road_relation_enum road_surface_relation;//! how does the object relate to road surface : allowed : IN/OUT/BORDER/UNDEF
    double precision ;
    double importance ;
    double dist_to_border ;

    //! function to get an idea of what is in the observation
    string classificationToString(){
        //#obs_id::int;X::double;Y::double;Z::double;confidence::double;weight::double
        std::ostringstream nstring;
        //nstring.precision(10);
        nstring << "class_id : " << int(class_id)  << std::endl
                <<"\t class_name : " << class_name << std::endl
                << "\t geom_type : " << geom_type<< ":"<<SnapEnums::gt_toString(geom_type) << std::endl
                << "\t road_surface_relation :"<<road_surface_relation<< ":"<<SnapEnums::rre_toString(road_surface_relation) << std::endl
                << "\t precision : " << precision << std::endl
                << "\t importance : " << importance << std::endl
                << "\t dist_to_border : " << dist_to_border <<std::endl;
        return nstring.str() ;
    }

    friend std::ostream &operator<<(std::ostream &os, classification const &c) {
        return os << c.class_id  << "\t"
                  << c.class_name << "\t"
                  << c.geom_type    << "\t"
                  << c.road_surface_relation   << "\t"
                  << c.precision   << "\t"
                  << c.importance   << "\t"
                  << c.dist_to_border;
    }

    void setClassification(int class_id_
                   ,std::string class_name_
                   ,std::string geom_type_
                   ,std::string road_surface_relation_
                   ,double precision_
                   ,double importance_
                   ,double dist_to_border_ ){

        this->class_id = int(class_id_);
        this->class_name = class_name_;
        this->geom_type = SnapEnums::String_togt(geom_type_);
        this->road_surface_relation = SnapEnums::String_torre(road_surface_relation_);
        this->precision = abs(precision_) ;
        this->importance = abs(importance_) ;
        this->dist_to_border = abs(dist_to_border_);
    }
};


typedef std::unordered_multimap <int /*node_id*/, edge *> ummap_e;

class DataStorage {
public:
    DataStorage(const string, const string );

    ~DataStorage();
    void readData();
    void readClassifications();
    void writeData(int);
    void setMap();


    int num_nodes()             { return num_nodes_;  }
    int num_edges()             { return num_edges_;    }
    int num_observations()       { return num_observations_;  }
    int num_classifications()   {return num_classifications_; }
    const std::string output_file_path()     { return output_file_path_;  }
    node* nodes()                { return nodes_;}
    node* nodes(int i)           { return &nodes_[i];}
    edge* edges()                { return edges_;}
    edge* edges(int i)                { return &edges_[i];}
    observation* observations()  { return observations_;}
    observation* observations(int i)  { return &observations_[i];}
    classification* classifications()  { return classifications_;}
    classification* classifications(int i)  { return &classifications_[i];}
    node* nbn(int i ) { return nodes_by_node_id_.at(i); }
    std::unordered_map <int /*node_id*/, node *> nodes_by_node_id() {return nodes_by_node_id_;}
    edge* ebe(int i ) { return edges_by_edge_id_.at(i); }
    std::unordered_map <int /*edge_id*/, edge *> edges_by_edge_id()
    {return edges_by_edge_id_;}
    std::pair <ummap_e::iterator, ummap_e::iterator> ebn(int i ) { return edges_by_node_id_.equal_range(i); }
    ummap_e * edges_by_node_id()
    {return &edges_by_node_id_;}
    std::unordered_map <int /*class_id*/, classification *>*  classification_by_id()
    {return &classification_by_id_;}
     classification *  cbi(int i){
        return classification_by_id_.at(i);
    }

    std::unordered_map <std::string /*class_name*/, classification *>* classification_by_name()
    {return &classification_by_name_;}
    classification *  cbn(string cname){
        return classification_by_name_.at(cname);
    }

private:
    int num_nodes_; //! total num of nodes we are going to read
    int num_edges_; //! total num of edges we are going to read
    int num_observations_;//! total num of observations we are going to read
    int num_classifications_;//! total num of classifications we have read

    node* nodes_;//! an array of node
    edge* edges_;//! an array og edges
    observation* observations_;//! an array of observations
    classification* classifications_; //! an array of classification class.
    
    const string input_file_path_;//! name of the file containing the input data
    const string output_file_path_;//! name of the file where to write results

    //! a hash table
    std::unordered_map <int /*node_id*/, node *> nodes_by_node_id_;
    std::unordered_map <int /*edge_id*/, edge *> edges_by_edge_id_;
    ummap_e edges_by_node_id_;
    std::unordered_map <int /*class_id*/, classification *> classification_by_id_;
    std::unordered_map <std::string /*class_name*/, classification *> classification_by_name_;
};

#endif // DATA_H
