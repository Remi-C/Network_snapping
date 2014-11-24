#include "geometry_function.h"


using std::endl ;
using std::cout ;
using std::string;

Geometry::Geometry()
{
}


geometry read_WKT(std::string s ){
    //    cout << "reading wkt \n" ;

    //    GEOSMessageHandler notice_function;
    //    GEOSMessageHandler error_function;
    //    initGEOS(notice_function,error_function);
    GEOSWKTReader GEOS_DLL* reader ;
    reader = GEOSWKTReader_create();
    geometry readed_geom ;
    readed_geom = GEOSWKTReader_read(reader, s.c_str());

    GEOSWKTReader_destroy(reader);
    //    finishGEOS();

    //    printf("end of reading wkt \n");
    return readed_geom;
}

/// write a geom into wkt
char *write_WKT(const geometry input_geom,const int dim){

    //    GEOSMessageHandler notice_function;
    //    GEOSMessageHandler error_function;
    //    initGEOS(notice_function,error_function);


    GEOSWKTWriter GEOS_DLL * writer ;
    writer = GEOSWKTWriter_create();

    GEOSWKTWriter_setOutputDimension( writer, dim);
    GEOSWKTWriter_setTrim(writer,1);
    GEOSWKTWriter_setRoundingPrecision(writer,10);

    char* wkt_geom;
    wkt_geom = GEOSWKTWriter_write(writer, input_geom) ;

    GEOSWKTWriter_destroy( writer);
    //    finishGEOS();

    return wkt_geom;
}

void EigenToCoordinate_Seq(Eigen::Vector3d tmp_point, GEOSCoordSequence GEOS_DLL * coor, int position){
    GEOSCoordSeq_setX(coor,position, tmp_point.x() ) ;
    GEOSCoordSeq_setY(coor,position, tmp_point.y() ) ;
    //GEOSCoordSeq_setZ(coor,position, tmp_point.z() ) ;
}

/** compute the width of the geometry orthogonaly to the axis
  */
double geometry_width_regarding_axis(const double* pt1,const double* pt2, geometry object, geometry centroid){
    //create proxy lines with big offset
    //measure distance between lines offsted and object
    //differnce of distance  is width
    double max_width = 100;
    ConstVectorRef Ni( pt1 ,3 );
    ConstVectorRef Nj( pt2 ,3 );
    Eigen::Vector3d NormalVect = ((Ni-Nj).cross(Eigen::Vector3d::UnitZ())).normalized();

    Eigen::Vector3d Nil = Ni + max_width * NormalVect ;
    Eigen::Vector3d Njl = Nj + max_width * NormalVect ;

    Eigen::Vector3d Nir = Ni - max_width * NormalVect ;
    Eigen::Vector3d Njr = Nj - max_width * NormalVect ;

    GEOSCoordSequence GEOS_DLL * axis_point = GEOSCoordSeq_create(2,2);
    EigenToCoordinate_Seq(Ni,axis_point,0);
    EigenToCoordinate_Seq(Nj,axis_point,1);
    geometry axis  = GEOSGeom_createLineString(axis_point);

    GEOSCoordSequence GEOS_DLL * axisl_points= GEOSCoordSeq_create(2,2);
    EigenToCoordinate_Seq(Nil,axisl_points,0);
    EigenToCoordinate_Seq(Njl,axisl_points,1);

    GEOSCoordSequence GEOS_DLL * axisr_points= GEOSCoordSeq_create(2,2);
    EigenToCoordinate_Seq(Nir,axisr_points,0);
    EigenToCoordinate_Seq(Njr,axisr_points,1);

    geometry axisl = GEOSGeom_createLineString(axisl_points);
    geometry axisr = GEOSGeom_createLineString(axisr_points);

    double distl;double distr;double dista ;
    GEOSDistance(object, axisl, &distl);
    GEOSDistance(object, axisr, &distr);
    GEOSDistance(centroid,axis, &dista);

    return  std::abs( 2.0*max_width - distl-distr);
}


geometry axis_to_rectangle(const double * pt1, const double * pt2, double axis_width, geometry * axis_to_be_filled){

    //cout << "computing the rectangle" << endl;
    //    GEOSMessageHandler notice_function;
    //    GEOSMessageHandler error_function;
    //    initGEOS(notice_function,error_function);


    ConstVectorRef Ni( pt1 ,3 );
    ConstVectorRef Nj( pt2 ,3 );
    GEOSCoordSequence GEOS_DLL * axis_points;
    GEOSCoordSequence GEOS_DLL * rectangle_points;
    Eigen::Vector3d tmp_point;
    geometry rectangle;

    axis_points = GEOSCoordSeq_create(2,2);

    EigenToCoordinate_Seq(Ni,axis_points,0);
    EigenToCoordinate_Seq(Nj,axis_points,1);

    * axis_to_be_filled = GEOSGeom_createLineString(axis_points);

    Eigen::Vector3d u = (Ni-Nj).normalized(); //director vector of segment, normalized
    Eigen::Vector3d normal = u.cross(Eigen::Vector3d::UnitZ()); //normal of segment

    //cout << "u : " << u.transpose() << " normal : " << normal.transpose()  << endl;
    //computing the 4 points of the rectangle :
    //pts1 + n*width, pts2 + n*width, p2-n*width, p1-n*width

    rectangle_points = GEOSCoordSeq_create(5,2);

    //(!rectangle_points)?printf("coordinate seq is null, wtf?\n"):printf("coordinate seq is not null\n");

    tmp_point = Ni + axis_width*normal;
    EigenToCoordinate_Seq(tmp_point,rectangle_points,0);
    tmp_point = Nj + axis_width*normal;
    EigenToCoordinate_Seq(tmp_point,rectangle_points,1);
    tmp_point = Nj - axis_width*normal;
    EigenToCoordinate_Seq(tmp_point,rectangle_points,2);
    tmp_point = Ni - axis_width*normal;
    EigenToCoordinate_Seq(tmp_point,rectangle_points,3);
    tmp_point = Ni + axis_width*normal;
    EigenToCoordinate_Seq(tmp_point,rectangle_points,4);

    rectangle = GEOSGeom_createPolygon(GEOSGeom_createLinearRing(rectangle_points)//shell
                                       ,NULL//holes
                                       ,0// nholes
                                       );
    //cout << "rectangle created : " <<  write_WKT(rectangle,3) << endl ;
    //finishGEOS();
    return rectangle;

}



double shared_area_cost(SnapEnums::road_relation_enum road_relation, const double* pt1, const double* pt2, double axis_width, geometry object_snapping_surface, double object_snapping_surface_area, geometry centroid  ){
    /**
      @param road_relation : what is the behaviour of the object toward road surface
      @param pt1 : first node
      @param pt2 : second node
      @param axis_to_rectangle(const double * pt1, const double * pt2, double axis_width){ pt2 : second node
      @param axis_width : width of the segment [pt1,pt2]
      @param object_snapping_surface : the surface of the object dilated by distance to border
      @param object_snapping_surface_area : area in square meter of the object_surface dilated
    */
    //    GEOSMessageHandler notice_function;
    //    GEOSMessageHandler error_function;
    //    initGEOS(notice_function,error_function);

    geometry street_rectangle;
    geometry axis;
    int intersects ; // 1 = true
    double distance_to_shell =  100 ;
    SnapEnums::attractive_repulsive attractive ;
    if(road_relation==SnapEnums::IN || road_relation==SnapEnums::BORDER_IN ){attractive = SnapEnums::ATTRACTIVE;}
    if(road_relation==SnapEnums::OUT || road_relation==SnapEnums::BORDER_OUT ){attractive = SnapEnums::REPULSIVE;}
    if(road_relation==SnapEnums::BORDER){attractive = SnapEnums::ATTR_AND_REP;};
    //{IN=1 ,OUT=-1 ,BORDER=0, BORDER_IN = 10, BORDER_OUT= -10, UNDEF=-110 } ;
    double object_width = 0;

    object_width = geometry_width_regarding_axis(pt1,pt2,object_snapping_surface,centroid) ;

    //compute the rectangle from pts
    street_rectangle = axis_to_rectangle(pt1,pt2, axis_width/2.0,&axis) ;
    double dist_to_axis = 0;
    double cost_surface  = 0 ;
    double cost_distance = 0;
    double shared_area = 0 ;

    //cout << "axis : " <<  write_WKT(axis,2) << endl;
    //cout << "working on geom : rectangle, obj : " << write_WKT(street_rectangle,3)
    //     << " , " << write_WKT(object_snapping_surface,3) << endl;

    //possibilities
    if(road_relation == SnapEnums::UNDEF){return 0 ; } //nothing to do


    //road relation must be IN=1 ,OUT=-1 ,BORDER=0, BORDER_IN = 10, BORDER_OUT

    //need to compute dist_to_axis : NOTE : could be done only when half within
    GEOSDistance(axis, street_rectangle, &dist_to_axis);

    intersects = GEOSIntersects(street_rectangle , object_snapping_surface) ;

    //this is a shortcut trick to avoid computing intersection and/or distance when not necessary
    if(intersects==1){  //fully or half within

        GEOSArea(GEOSIntersection(street_rectangle, object_snapping_surface)
                 , &shared_area) ;

        if( TOLERANCY_EQUAL(object_snapping_surface_area,shared_area) ){//fully within, need to compute distance

            GEOSDistance( GEOSGetExteriorRing(street_rectangle) , object_snapping_surface, &distance_to_shell);

        }else{//half within, no need to compute distance, it's 0
            distance_to_shell = 0;

        }
    }else{ //fully outside
        //, no need to compute intersection, but need to compute distance
        shared_area = 0 ;
        GEOSDistance( GEOSGetExteriorRing(street_rectangle) , object_snapping_surface, &distance_to_shell);
    }
    int sign=0;


    //putting a cost based on shared area.
    //we take some shortcuts to avoid computing intersection if it's not necessary

    if(distance_to_shell!=0){//the object is either fully inside or fully outside
        if(intersects==1){//the object is fully inside
            if(road_relation==SnapEnums::BORDER){
                cost_surface = object_snapping_surface_area;//object_snapping_surface_area;
                sign= +1;
            }
            if(attractive==SnapEnums::ATTRACTIVE){
                cost_surface = 0;
                if(road_relation == SnapEnums::BORDER_IN){
                    //adding a cost proportionnal to distance
                    cost_distance = distance_to_shell * object_snapping_surface_area ;
                }
                sign = 1;//is no 0 because of BORDER_IN case
            }
            if(attractive==SnapEnums::REPULSIVE){
                cost_surface = object_snapping_surface_area  ;//object_snapping_surface_area;
                cost_distance = (distance_to_shell+object_width) * object_snapping_surface_area ; //object_snapping_surface_area ;
                sign = +1 ;
            }
        }else{//the object is fully outside
            if(road_relation==SnapEnums::BORDER){
                cost_surface = 1 * object_snapping_surface_area;
                sign = -1;
            }
            if(attractive==SnapEnums::ATTRACTIVE){
                cost_surface = 1 * object_snapping_surface_area;
                cost_distance = (distance_to_shell+object_width) * object_snapping_surface_area;//object_snapping_surface_area ;
                sign = -1;
            }
            if(attractive==SnapEnums::REPULSIVE){
                cost_surface = 0;
                if(road_relation == SnapEnums::BORDER_OUT){
                    //adding a cost proportionnal to distance
                    cost_distance = distance_to_shell * object_snapping_surface_area ;
                }
                sign = -1;//is not 0 because of border_out
            }

        }
    }else{//the object intersects the border
        //we must compute the shared surface

        if(road_relation==SnapEnums::BORDER){//cost is 0 when object is centered on border
            cost_surface = (object_snapping_surface_area-2.0*shared_area)* object_width/2.0  ;
            sign = -1*SIGN(cost_surface);//this is the sign function
            cost_surface = std::abs(cost_surface);
        }
        if(attractive==SnapEnums::ATTRACTIVE){//Cost is high when object is outside
            cost_surface = (object_snapping_surface_area-shared_area) * object_width ;
            sign = -1;
        }
        if(attractive==SnapEnums::REPULSIVE){//cost is high when object is inside
            cost_surface = shared_area * object_width;
            sign = +1;
        }

    }

    if(road_relation == SnapEnums::BORDER ){
        //we add a term to cost that is proprortionnal to the distance to border

        if(distance_to_shell==0){
            //the distance to border is null, no modification of the cost
        }else{//the object is either fully inside or fully outside
            //the cost is proportionnal to its distance to border
            cost_distance = (distance_to_shell+object_width/2.0) * object_snapping_surface_area ;
        }
    }

    /** we compute a sign that gives the direction into witch make the modification

        CASE(object left of axis) : A : object fully outside. B : object hlaf inside, C: object fully inside
        CASE                A   B   C
        INTERSECTION==TRUE  0   1   1
        DISTANCE==0         0   1   0

        ATTRACTIVE==TRUE    1   1   0
        REPULSIVE==TRUE     1   0   0

        Thus, ATTRACTIVE = ((INTERSECTION==TRUE) != (DISTANCE==0))
        REPULSIVE = (!(INTERSECTION)) && (!(DISTANCE==0))
    */
    //bool att = (attractive==SnapEnums::ATTRACTIVE) && ((intersects==1)!=(distance_to_shell==0));//0 or 1

    //int sign = !att;//only one of the 2 may contribute at the same time
    //sign= (sign*2)-1;//putting sign to value -1 or 1

//    cout << "\t  total_cost : " << sign *( cost_surface + cost_distance) <<endl;
//    cout << "\t  cost_surface : " << cost_surface
//         <<" , cost_distance : " << cost_distance <<endl ;
    //    cout << "\t \t road_relation type :" << road_relation  <<endl ;
    //    cout << "\t \t attractivity ; " <<  attractive <<endl;
//    cout << "\t \t distance_to_shell : " << distance_to_shell
            //         << " , intersects? "<< intersects
//         <<endl ;

    return sign* (cost_surface+cost_distance) ;

}


void initialize_geom_computation(){
    GEOSMessageHandler notice_function;
    GEOSMessageHandler error_function;
    initGEOS(notice_function,error_function);
}
void finish_geom_computation(){
    finishGEOS();
}



double signed_dist_to_border(SnapEnums::road_relation_enum road_relation, const double* pt1, const double* pt2, double axis_width, geometry object_snapping_surface , double object_snapping_surface_area   ){
    /**
      @param road_relation : what is the behaviour of the object toward road surface
      @param pt1 : first node
      @param pt2 : second node
      @param axis_to_rectangle(const double * pt1, const double * pt2, double axis_width){ pt2 : second node
      @param axis_width : width of the segment [pt1,pt2]
      @param object_snapping_surface : the surface of the object dilated by distance to border
    */
    //    GEOSMessageHandler notice_function;
    //    GEOSMessageHandler error_function;
    //    initGEOS(notice_function,error_function);

    geometry street_rectangle;
    geometry axis;
    int intersects ; // 1 = true
    double distance_to_shell =  100 ;
    SnapEnums::attractive_repulsive attractive ;
    if(road_relation==SnapEnums::IN || road_relation==SnapEnums::BORDER_IN ){attractive = SnapEnums::ATTRACTIVE;}
    if(road_relation==SnapEnums::OUT || road_relation==SnapEnums::BORDER_OUT ){attractive = SnapEnums::REPULSIVE;}
    if(road_relation==SnapEnums::BORDER){attractive = SnapEnums::ATTR_AND_REP;};
    //{IN=1 ,OUT=-1 ,BORDER=0, BORDER_IN = 10, BORDER_OUT= -10, UNDEF=-110 } ;

    //compute the rectangle from pts
    street_rectangle = axis_to_rectangle(pt1,pt2, axis_width/2.0,&axis) ;
    double dist_to_axis = 0;
    double cost_surface  = 0 ;
    double cost_distance = 0;
    double shared_area = 0 ;

    //cout << "axis : " <<  write_WKT(axis,2) << endl;
    //cout << "working on geom : rectangle, obj : " << write_WKT(street_rectangle,3)
    //     << " , " << write_WKT(object_snapping_surface,3) << endl;

    //possibilities
    if(road_relation == SnapEnums::UNDEF){return 0 ; } //nothing to do


    //road relation must be IN=1 ,OUT=-1 ,BORDER=0, BORDER_IN = 10, BORDER_OUT

    //need to compute dist_to_axis : NOTE : could be done only when half within
    GEOSDistance(axis, object_snapping_surface, &dist_to_axis);
    intersects = GEOSIntersects(street_rectangle , object_snapping_surface) ;
    GEOSDistance( GEOSGetExteriorRing(street_rectangle) , object_snapping_surface, &distance_to_shell);


    double intensity = 0;
    int activating =0 ;
    int sign;

    /// Only dealing with REPULSIVE  : OUT and BORDER_OUT for the moment
    if(attractive == SnapEnums::REPULSIVE){
       //out or border_out
        intensity =std::max( std::abs(dist_to_axis-axis_width/2), distance_to_shell)  ;
        activating  = 1 ;
        if(road_relation==SnapEnums::OUT){//we need to say that being fully outside is 0
            if( (dist_to_axis-axis_width/2) > 0 ){//we are fully outside
                activating  =0 ;
            }
        }
    }

    if(intersects==1){
        sign = -1;
    }else{
        sign=+1;
    }

    /*



    //this is a shortcut trick to avoid computing intersection and/or distance when not necessary
    if(intersects==1){  //fully or half within

        GEOSArea(GEOSIntersection(street_rectangle, object_snapping_surface)
                 , &shared_area) ;

        if( TOLERANCY_EQUAL(object_snapping_surface_area,shared_area) ){//fully within, need to compute distance

            GEOSDistance( GEOSGetExteriorRing(street_rectangle) , object_snapping_surface, &distance_to_shell);

        }else{//half within, no need to compute distance, it's 0
            distance_to_shell = 0;

        }
    }else{ //fully outside
        //, no need to compute intersection, but need to compute distance
        shared_area = 0 ;
        GEOSDistance( GEOSGetExteriorRing(street_rectangle) , object_snapping_surface, &distance_to_shell);
    }
    int sign=0;


    //putting a cost based on shared area.
    //we take some shortcuts to avoid computing intersection if it's not necessary

    if(distance_to_shell!=0){//the object is either fully inside or fully outside
        if(intersects==1){//the object is fully inside
            if(road_relation==SnapEnums::BORDER){
                cost_surface = object_snapping_surface_area;
                sign= +1;
            }
            if(attractive==SnapEnums::ATTRACTIVE){
                cost_surface = 0;
                sign = 1;//is no 0 because of BORDER_IN case
            }
            if(attractive==SnapEnums::REPULSIVE){
                cost_surface = object_snapping_surface_area;
                cost_distance = (distance_to_shell) * object_snapping_surface_area ;
                sign = +1 ;
            }
        }else{//the object is fully outside
            if(road_relation==SnapEnums::BORDER){
                cost_surface = object_snapping_surface_area;
                sign = -1;
            }
            if(attractive==SnapEnums::ATTRACTIVE){
                cost_surface = object_snapping_surface_area;
                cost_distance = (distance_to_shell) * object_snapping_surface_area ;
                sign = -1;
            }
            if(attractive==SnapEnums::REPULSIVE){
                cost_surface = 0;
                sign = -1;//is not 0 because of border_out
            }

        }
    }else{//the object intersects the border
        //we must compute the shared surface

        if(road_relation==SnapEnums::BORDER){//cost is 0 when object is centered on border
            cost_surface = object_snapping_surface_area-2*shared_area ;
            sign = -1*SIGN(cost_surface);//this is the sign function
            cost_surface = std::abs(cost_surface);
        }
        if(attractive==SnapEnums::ATTRACTIVE){//Cost is high when object is outside
            cost_surface = object_snapping_surface_area-shared_area;
            sign = -1;
        }
        if(attractive==SnapEnums::REPULSIVE){//cost is high when object is inside
            cost_surface = shared_area;
            sign = +1;
        }

    }

    if(road_relation == SnapEnums::BORDER || road_relation == SnapEnums::BORDER_IN || road_relation == SnapEnums::BORDER_OUT){
        //we add a term to cost that is proprortionnal to the distance to border

        if(distance_to_shell==0){
            //the distance to border is null, no modification of the cost
        }else{//the object is either fully inside or fully outside
            //the cost is proportionnal to its distance to border
            cost_distance = (distance_to_shell) * object_snapping_surface_area ;
        }
    }
    */
    /** we compute a sign that gives the direction into witch make the modification

        CASE(object left of axis) : A : object fully outside. B : object hlaf inside, C: object fully inside
        CASE                A   B   C
        INTERSECTION==TRUE  0   1   1
        DISTANCE==0         0   1   0

        ATTRACTIVE==TRUE    1   1   0
        REPULSIVE==TRUE     1   0   0

        Thus, ATTRACTIVE = ((INTERSECTION==TRUE) != (DISTANCE==0))
        REPULSIVE = (!(INTERSECTION)) && (!(DISTANCE==0))
    */
    //bool att = (attractive==SnapEnums::ATTRACTIVE) && ((intersects==1)!=(distance_to_shell==0));//0 or 1

    //int sign = !att;//only one of the 2 may contribute at the same time
    //sign= (sign*2)-1;//putting sign to value -1 or 1

//    cout << "\t  total_cost : " << sign *( cost_surface + cost_distance) <<endl;
//    cout << "\t  cost_surface : " << cost_surface
//         <<" , cost_distance : " << cost_distance <<endl ;
    //    cout << "\t \t road_relation type :" << road_relation  <<endl ;
    //    cout << "\t \t attractivity ; " <<  attractive <<endl;
//    cout << "\t \t distance_to_shell : " << distance_to_shell
            //         << " , intersects? "<< intersects
//         <<endl ;

    //return sign* (cost_surface+cost_distance) ;

    return sign*activating * intensity;

}


