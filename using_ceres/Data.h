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
#include <vector>
#include <typeinfo> //to get the name of a class
#include <cxxabi.h> // to translate name of class into human language


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
    double width[1];     //! width of the edge, in meters.

    //! function to get an idea of what is in the edge
    string edgeToString(){
        //#edge_id::int;start_node::int;end_node::int;width::double

        std::ostringstream nstring;
        nstring << "(edge_id : " << edge_id
                << " , start_node : " << start_node << ", end_node : " << end_node <<", width : " << width[0] <<")";
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
    string classificationToString() const {
        //#obs_id::int;X::double;Y::double;Z::double;confidence::double;weight::double
        std::ostringstream nstring;
        //nstring.precision(10);
        nstring  << "class_id : " << int(class_id)  << std::endl
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
        this->precision = std::abs(precision_) ;
        this->importance = std::abs(importance_) ;
        this->dist_to_border = std::abs(dist_to_border_);
    }
};



struct street_object{

    unsigned char object_id;           //! unique id per object, positive, 8 bits long, contains kind of information on object type grouping
    unsigned char class_id;
    std::string class_name;          //! name of the class, in english
    int edge_id ;
    geometry geom;
    geometry geom_centroid;
    geometry geom_border_surface;
    double geom_border_area;
    double confidence;


    //! function to get an idea of what is in the observation
    string street_objectsToString() const {

        std::ostringstream nstring;
        //nstring.precision(10);
        nstring << "object_id : " << int(object_id)  << std::endl;
        nstring <<"\t class_id : " << int(class_id) << std::endl;
        nstring <<"\t class_name : " << class_name << std::endl;
        nstring << "\t edge_id : " << edge_id<< ":"<<  std::endl;
        nstring << "\t geom :"<< write_WKT(geom,3) << std::endl;
        nstring << "\t geom_centroid :"<< write_WKT(geom_centroid,3) << std::endl;
        nstring << "\t geom_border_surface :"<< write_WKT(geom_border_surface,3) << std::endl;
        nstring << "\t geom_border_area :"<< geom_border_area << std::endl;
        nstring << "\t confidence : " << confidence   << std::endl;
        return nstring.str() ;
        //return string("toto");
    }


    void setStreet_Object(
            unsigned char object_id_,
            unsigned char class_id_,
            std::string class_name_,
            int edge_id_,
            std::string geom_wkt,
            double radius,
            double confidence_){

        this->object_id = int(object_id_);
        this->class_id = int(class_id_);
        this->class_name =  string(class_name_) ;

        //finding the classification associated with this object



        this->edge_id =  int(edge_id_);
        this->geom = read_WKT(geom_wkt) ;
        this->geom_centroid = Geometry::centroid(this->geom) ;
        this->geom_border_surface =
                Geometry::BufferWithStyle(
                    this->geom
                    ,radius//width
                    ,0 //quadseg
                    ,CAP_FLAT //endcapstyle
                    ,JOIN_ROUND //joinStyle
                    ,0 //mitreLimit
                    );
        this->geom_border_area = Geometry::area(this->geom_border_surface);
        this->confidence =  confidence_ ;
    }

};






/** A class to keep info about cost function so to be able to output constraints at each step
*/
class Constraint {

public:
    //put destructor to virtual
    Constraint(){};
    //constructor , do almost nothing
    Constraint(double node_1[3],double node_2[3],double node_3[3], double * width, ceres::CostFunction * cost_function , double application_point[3] ){
        node_1_ = node_1;
        node_2_ = node_2;
        node_3_ = node_3;
        width_= width;
        cost_function_ = cost_function;
        application_point_ = application_point;
    }
    virtual ~Constraint(){
    };

    // this function should be overcharged by heriting classes to be adapted to each constraints
   virtual int get_graphical_constraint(double * cost, double * geom ) = 0;


    double* node_1_;//nodes used by the cost function
    double* node_2_;
    double* node_3_;
    double * width_; // width of the edge, used by some constraints
    ceres::CostFunction * cost_function_; //cost function of the constaint, contains the "evaluate" method
    double * application_point_; //graphical : where should visually the vector start?
private :
};

/// heriting from constraints to be adapted to sidewlak distance cost function
class Constraint_sidewalk : public Constraint{
public :
    Constraint_sidewalk(double node_1[3],double node_2[3],double node_3[3], double * width, ceres::CostFunction * cost_function , double application_point[3] )
        :
    Constraint( node_1 ,node_2, node_3,  width,  cost_function ,application_point )
    {
    }

    int get_graphical_constraint(double * cost, double * geom ){

        double temp_residuals ;
        double* * temp_jacobian = new double*[3]  ;
        double* * temp_parameters= new double*[3] ;
        temp_jacobian[0]=new double[3];temp_jacobian[1]=new double[3];temp_jacobian[2]=new double[3];
       // temp_parameters[0] = new double[3];temp_parameters[1] = new double[3];temp_parameters[2] = new double[3];
        temp_parameters[0] =node_1_ ;
        temp_parameters[1] =node_2_  ;
        temp_parameters[2]=width_ ;

        cost_function_->Evaluate(temp_parameters,
                                 &temp_residuals,
                                 temp_jacobian) ;
        *cost = temp_residuals;
        geom[0] = application_point_[0]+temp_jacobian[0][0];
        geom[1] = application_point_[1]+temp_jacobian[0][1];
        geom[2] = application_point_[2]+temp_jacobian[0][2];

        geom[3] =  application_point_[0] ;
        geom[4] =  application_point_[1] ;
        geom[5] =  application_point_[2] ;

        delete[] temp_jacobian[0];delete[] temp_jacobian[1];delete[] temp_jacobian[2];
        delete[] temp_jacobian;
        //delete[] temp_parameters[0];delete[] temp_parameters[1];delete[] temp_parameters[2];
        delete[] temp_parameters;
        return -1;
    };
};


/// heriting from constraints to be adapted to sidewlak distance cost function
class Constraint_objects : public Constraint{
public :
    Constraint_objects(double node_1[3],double node_2[3],double node_3[3], double * width, ceres::CostFunction * cost_function , double application_point[3] )
        :
    Constraint( node_1 ,node_2, node_3,  width,  cost_function ,application_point )
    {
    }

    int get_graphical_constraint(double * cost, double * geom ){

        double temp_residuals ;
        double* * temp_jacobian = new double*[3]  ;
        double* * temp_parameters= new double*[3] ;
        temp_jacobian[0]=new double[3];temp_jacobian[1]=new double[3];temp_jacobian[2]=new double[3];
       // temp_parameters[0] = new double[3];temp_parameters[1] = new double[3];temp_parameters[2] = new double[3];
        temp_parameters[0] = node_1_ ;
        temp_parameters[1] = node_2_  ;
        temp_parameters[2] = width_ ;

        cost_function_->Evaluate(temp_parameters,
                                 &temp_residuals,
                                 temp_jacobian) ;

        //projecting the centroid on axis :
        //creating a geom for axis
        geometry axis ;
        axis_to_rectangle(temp_parameters[0], temp_parameters[1], 1, &axis ) ;
        //projecting the centroid on axis
        double* closest_on_axis = new double[3] ;
        for(int i = 0; i<3;++i){closest_on_axis[i]=0;}
        Geometry::closestPoint(axis, Geometry::double2geomPoint(application_point_),closest_on_axis);
        closest_on_axis[2]=0;
        //translating the closest point on axis by width/ 2 toward object
        //rceating a normalized eigen vector out of centroid and closest point
        VectorRef caxis( closest_on_axis ,3 ) ;
        VectorRef centroid( application_point_ ,3 ) ;
        Eigen::Vector3d n_axis =  centroid-caxis   ;
        n_axis = n_axis.normalized() ;
        //translate by
        double w = *width_ ;
        caxis = caxis + n_axis * w / 2.0 ;


        *cost = temp_residuals;
        geom[0] =  closest_on_axis[0];
        geom[1] =  closest_on_axis[1];
        geom[2] =  closest_on_axis[2];

        geom[3] =  closest_on_axis[0]-temp_jacobian[0][0];
        geom[4] =  closest_on_axis[1]-temp_jacobian[0][1];
        geom[5] =  closest_on_axis[2]-temp_jacobian[0][2];

        delete[] temp_jacobian[0];delete[] temp_jacobian[1];delete[] temp_jacobian[2];
        delete[] temp_jacobian;
        //delete[] temp_parameters[0];delete[] temp_parameters[1];delete[] temp_parameters[2];
        delete[] temp_parameters;
        return -1;
    };
};


/// heriting from constraints to be adapted to initial position cost function
class Constraint_origin : public Constraint{
public :
    Constraint_origin(double node_1[3],double node_2[3],double node_3[3], double * width, ceres::CostFunction * cost_function , double application_point[3] )
        :
    Constraint( node_1 ,node_2, node_3,  width,  cost_function ,application_point )
    {
    }

    int get_graphical_constraint(double * cost, double * geom ){

        double temp_residuals ;
        double* * temp_jacobian = new double*[3]  ;
        double* * temp_parameters= new double*[3] ;
        temp_jacobian[0]=new double[3];temp_jacobian[1]=new double[3];temp_jacobian[2]=new double[3];
       // temp_parameters[0] = new double[3];temp_parameters[1] = new double[3];temp_parameters[2] = new double[3];
        temp_parameters[0] =node_1_ ;

        cost_function_->Evaluate(temp_parameters,
                                 &temp_residuals,
                                 temp_jacobian) ;
        *cost = temp_residuals;
        geom[0] = application_point_[0] ;
        geom[1] = application_point_[1] ;
        geom[2] = application_point_[2] ;

        geom[3] =  application_point_[0]+temp_jacobian[0][0]  ;
        geom[4] =  application_point_[1]+temp_jacobian[0][1]  ;
        geom[5] =  application_point_[2]+temp_jacobian[0][2]  ;

        delete[] temp_jacobian[0];delete[] temp_jacobian[1];delete[] temp_jacobian[2];
        delete[] temp_jacobian;
        //delete[] temp_parameters[0];delete[] temp_parameters[1];delete[] temp_parameters[2];
        delete[] temp_parameters;
        return -1;
    };
};




/// heriting from constraints to be adapted to initial position cost function
class Constraint_spacing : public Constraint{
public :
    Constraint_spacing(double node_1[3],double node_2[3],double node_3[3], double * width, ceres::CostFunction * cost_function , double application_point[3] )
        :
    Constraint( node_1 ,node_2, node_3,  width,  cost_function ,application_point )
    {
    }

    int get_graphical_constraint(double * cost, double * geom ){

        double temp_residuals ;
        double* * temp_jacobian = new double*[3]  ;
        double* * temp_parameters= new double*[3] ;
        temp_jacobian[0]=new double[3];temp_jacobian[1]=new double[3];temp_jacobian[2]=new double[3];
       // temp_parameters[0] = new double[3];temp_parameters[1] = new double[3];temp_parameters[2] = new double[3];
        temp_parameters[0] =node_1_ ;
        temp_parameters[1] =node_2_  ;
        temp_parameters[2]=width_ ;

        cost_function_->Evaluate(temp_parameters,
                                 &temp_residuals,
                                 temp_jacobian) ;
        *cost = temp_residuals;
        geom[0] = application_point_[0] ;
        geom[1] = application_point_[1] ;
        geom[2] = application_point_[2] ;

        geom[3] =  application_point_[0]+temp_jacobian[0][0]  ;
        geom[4] =  application_point_[1]+temp_jacobian[0][1]  ;
        geom[5] =  application_point_[2]+temp_jacobian[0][2]  ;

        delete[] temp_jacobian[0];delete[] temp_jacobian[1];delete[] temp_jacobian[2];
        delete[] temp_jacobian;
        //delete[] temp_parameters[0];delete[] temp_parameters[1];delete[] temp_parameters[2];
        delete[] temp_parameters;
        return -1;
    };
};


/// heriting from constraints to be adapted to initial angle cost function
class Constraint_angle : public Constraint{
public :
    Constraint_angle(double node_1[3],double node_2[3],double node_3[3], double * width, ceres::CostFunction * cost_function , double application_point[3] )
        :
    Constraint( node_1 ,node_2, node_3,  width,  cost_function ,application_point )
    {
    }

    int get_graphical_constraint(double * cost, double * geom ){

        double  temp_residuals[2] = {0,0};
        double* * temp_jacobian = new double*[3]  ;
        double* * temp_parameters= new double*[3] ;
        temp_jacobian[0]=new double[6];temp_jacobian[1]=new double[6];temp_jacobian[2]=new double[6];

        temp_parameters[0] =node_1_ ; //center node
        temp_parameters[1] =node_2_  ; //left node
        temp_parameters[2] =node_3_  ; //right node

        cost_function_->Evaluate(temp_parameters,
                                 temp_residuals,
                                 temp_jacobian) ;

        *cost =  temp_residuals[0]+temp_residuals[1];
        geom[0] = application_point_[0] ;
        geom[1] =  application_point_[1] ;
        geom[2] =  application_point_[2] ;

        geom[3] =  application_point_[0]+temp_jacobian[0][0] ;
        geom[4] =  application_point_[1]+temp_jacobian[0][1] ;
        geom[5] =  application_point_[2]+temp_jacobian[0][2] ;

        delete[] temp_jacobian[0];delete[] temp_jacobian[1];//delete[] temp_jacobian[2];
        //delete[] temp_residuals;
        delete[] temp_jacobian;
        //delete[] temp_parameters[0];delete[] temp_parameters[1];delete[] temp_parameters[2];
        delete[] temp_parameters;
        return -1;
    };
};


/// heriting from constraints to be adapted to sidewalk distance cost function
class Constraint_sidewalk_width : public Constraint{
public :
    Constraint_sidewalk_width(double node_1[3],double node_2[3],double node_3[3], double * width, ceres::CostFunction * cost_function , double application_point[3] )
        :
    Constraint( node_1 ,node_2, node_3,  width,  cost_function ,application_point )
    {
    }

    int get_graphical_constraint(double * cost, double * geom ){

        double temp_residuals ;
        double* * temp_jacobian = new double*[3]  ;
        double* * temp_parameters= new double*[3] ;
        temp_jacobian[0]=new double[3];temp_jacobian[1]=new double[3];temp_jacobian[2]=new double[3];
       // temp_parameters[0] = new double[3];temp_parameters[1] = new double[3];temp_parameters[2] = new double[3];
        temp_parameters[0] =node_1_ ;
        temp_parameters[1] =node_2_  ;
        temp_parameters[2]=width_ ;

        cost_function_->Evaluate(temp_parameters,
                                 &temp_residuals,
                                 temp_jacobian) ;

        //this force starts from observation position and goes in the normal direction
        //compute  normal for NiNj
        ConstVectorRef Ni( temp_parameters[0],3 );
        ConstVectorRef Nj( temp_parameters[1],3 );
        Eigen::Vector3d Np = Eigen::Vector3d::UnitZ() ;
        //compute director vector of (NiNj) ie : NiNj/norm(NiNj) = u
        Eigen::Vector3d U = (Nj-Ni)/(Nj-Ni).norm();
        //the direction of change should be :
        Eigen::Vector3d Vja = (U.cross(Np)).normalized();
        double sign = Geometry::orientationIndex(node_1_,node_2_,application_point_) ;
        *cost = temp_residuals;
        geom[0] = application_point_[0];
        geom[1] = application_point_[1];
        geom[2] = application_point_[2];

        geom[3] = application_point_[0] - sign * Vja[0] * temp_residuals ;
        geom[4] = application_point_[1] - sign * Vja[1] * temp_residuals ;
        geom[5] = application_point_[2] - sign * Vja[2] * temp_residuals ;

        delete[] temp_jacobian[0];delete[] temp_jacobian[1];delete[] temp_jacobian[2];
        delete[] temp_jacobian;
        //delete[] temp_parameters[0];delete[] temp_parameters[1];delete[] temp_parameters[2];
        delete[] temp_parameters;
        return -1;
    };
};


/// heriting from constraints to be adapted to sidewlak distance cost function
class Constraint_objects_width : public Constraint{
public :
    Constraint_objects_width(double node_1[3],double node_2[3],double node_3[3], double * width, ceres::CostFunction * cost_function , double application_point[3] )
        :
    Constraint( node_1 ,node_2, node_3,  width,  cost_function ,application_point )
    {
    }

    int get_graphical_constraint(double * cost, double * geom ){

        double temp_residuals ;
        double* * temp_jacobian = new double*[3]  ;
        double* * temp_parameters= new double*[3] ;
        temp_jacobian[0]=new double[3];temp_jacobian[1]=new double[3];temp_jacobian[2]=new double[3];
       // temp_parameters[0] = new double[3];temp_parameters[1] = new double[3];temp_parameters[2] = new double[3];
        temp_parameters[0] = node_1_ ;
        temp_parameters[1] = node_2_  ;
        temp_parameters[2] = width_ ;

        cost_function_->Evaluate(temp_parameters,
                                 &temp_residuals,
                                 temp_jacobian) ;

        //projecting the centroid on axis :
        //creating a geom for axis
        geometry axis ;
        axis_to_rectangle(temp_parameters[0], temp_parameters[1], 1, &axis ) ;
        //projecting the centroid on axis
        double* closest_on_axis = new double[3] ;
        for(int i = 0; i<3;++i){closest_on_axis[i]=0;}
        Geometry::closestPoint(axis, Geometry::double2geomPoint(application_point_),closest_on_axis);
        closest_on_axis[2]=0;
        //translating the closest point on axis by width/ 2 toward object
        //rceating a normalized eigen vector out of centroid and closest point
        VectorRef caxis( closest_on_axis ,3 ) ;
        VectorRef centroid( application_point_ ,3 ) ;
        Eigen::Vector3d n_axis =  centroid-caxis   ;
        n_axis = n_axis.normalized() ;
        //translate by
        double w = *width_ ;
        caxis = caxis + n_axis * w / 2.0 ;


        *cost = temp_residuals;
        geom[0] =  centroid[0];
        geom[1] =  centroid[1];
        geom[2] =  centroid[2];

        geom[3] =  centroid[0];
        geom[4] =  centroid[1];
        geom[5] =  centroid[2];

        delete[] temp_jacobian[0];delete[] temp_jacobian[1];delete[] temp_jacobian[2];
        delete[] temp_jacobian;
        //delete[] temp_parameters[0];delete[] temp_parameters[1];delete[] temp_parameters[2];
        delete[] temp_parameters;
        return -1;
    };
};

typedef std::unordered_multimap <int /*node_id*/, edge *> ummap_e;

class DataStorage {
public:
    DataStorage(const string, const string );

    ~DataStorage();
    void readData();
    void readClassifications();
    void readObjects();
    void writeData(int);
    void writeConstraints(int);
    void setMap();


    int num_nodes()             { return num_nodes_;  }
    int num_edges()             { return num_edges_;    }
    int num_observations()       { return num_observations_;  }
    int num_classifications()   {return num_classifications_; }
    int num_street_objects() {return num_street_objects_;}

    const std::string output_file_path()     { return output_file_path_;  }
    node* nodes()                { return nodes_;}
    node* nodes(int i)           { return &nodes_[i];}
    edge* edges()                { return edges_;}
    edge* edges(int i)                { return &edges_[i];}
    observation* observations()  { return observations_;}
    observation* observations(int i)  { return &observations_[i];}
    classification* classifications()  { return classifications_;}
    classification* classifications(int i)  { return &classifications_[i];}
    street_object* street_objects()  { return street_objects_;}
    street_object* street_objects(int i)  { return &street_objects_[i];}
    std::vector<Constraint*>*  constraints() {return &constraints_;} ;
    int numConstraintsWritten(){return numConstraintsWritten_ ; };
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
    // Copy and assignment declared private and not defined to prevent from copying
    // instances of this class
    DataStorage(const DataStorage& nocopy);
    DataStorage& operator=(const DataStorage& nocopy);

    int num_nodes_; //! total num of nodes we are going to read
    int num_edges_; //! total num of edges we are going to read
    int num_observations_;//! total num of observations we are going to read
    int num_classifications_;//! total num of classifications we have read
    int num_street_objects_;

    int numConstraintsWritten_;

    node* nodes_;//! an array of node
    edge* edges_;//! an array og edges
    observation* observations_;//! an array of observations
    classification* classifications_; //! an array of classification class.
    street_object * street_objects_;

    const string input_file_path_;//! name of the file containing the input data
    const string output_file_path_;//! name of the file where to write results

    //! a hash table
    std::unordered_map <int /*node_id*/, node *> nodes_by_node_id_;
    std::unordered_map <int /*edge_id*/, edge *> edges_by_edge_id_;
    ummap_e edges_by_node_id_;
    std::unordered_map <int /*class_id*/, classification *> classification_by_id_;
    std::unordered_map <std::string /*class_name*/, classification *> classification_by_name_;

    std::vector<Constraint*> constraints_ ;  //used to store constraints for output at each step
};

#endif // DATA_H
