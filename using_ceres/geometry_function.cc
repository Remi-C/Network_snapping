#include "geometry_function.h"


using std::endl ;
using std::cout ;
using std::string;

geometry_function::geometry_function()
{
}


geom read_WKT(std::string s ){
//    cout << "reading wkt \n" ;

//    GEOSMessageHandler notice_function;
//    GEOSMessageHandler error_function;
//    initGEOS(notice_function,error_function);
    GEOSWKTReader GEOS_DLL* reader ;
    reader = GEOSWKTReader_create();
    geom readed_geom;
    readed_geom = GEOSWKTReader_read(reader, s.c_str());

    GEOSWKTReader_destroy(reader);
//    finishGEOS();

//    printf("end of reading wkt \n");
    return readed_geom;
}

/// write a geom into wkt
char *write_WKT(geom input_geom, int dim){

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


geom axis_to_rectangle(const double * pt1, const double * pt2, double axis_width){

    //cout << "computing the rectangle" << endl;
//    GEOSMessageHandler notice_function;
//    GEOSMessageHandler error_function;
//    initGEOS(notice_function,error_function);


    ConstVectorRef Ni( pt1 ,3 );
    ConstVectorRef Nj( pt2 ,3 );
    GEOSCoordSequence GEOS_DLL * rectangle_points;
    Eigen::Vector3d tmp_point;
    geom rectangle;

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



double shared_area_cost(attractive_repulsive attractive, const double* pt1, const double* pt2, double axis_width, geom object_snapping_surface, double object_snapping_surface_area  ){
    /**
      @param attractive : 1 if the object is attractive, 0 if the object is repulsive
      @param pt1 : first node
      @param axis_to_rectangle(const double * pt1, const double * pt2, double axis_width){ pt2 : second node
      @param axis_width : width of the segment [pt1,pt2]
      @param object_snapping_surface : the surface of the object dilated by distance to border
      @param object_snapping_surface_area : area in square meter of the object_surface dilated

      compute a polygon from pt1,pt2,axis_width aka rectangle
      if rectangle and obj_surface don't intersect, cost = area(obj_surface)(1+dist(obj_surface, rectangle))
      if rectangle and obj_surface intersect, cost = shared surface
    */
//    GEOSMessageHandler notice_function;
//    GEOSMessageHandler error_function;
//    initGEOS(notice_function,error_function);

    geom street_rectangle;
    int intersects ; // 1 = true
    double distance_to_shell =  100 ;

    //compute the rectangle from pts
    street_rectangle = axis_to_rectangle(pt1,pt2, axis_width) ;
    double cost  = 0 ;
    double shared_area = 0 ;

    GEOSSetSRID(street_rectangle, 0);
    GEOSSetSRID(object_snapping_surface, 0);

    cout << "working on geom : rectangle, obj : " << write_WKT(street_rectangle,3)
         << " , " << write_WKT(object_snapping_surface,3) << endl;
    cout << "validity : " << int(GEOSisValid(street_rectangle))
         << " " << int(GEOSisValid(object_snapping_surface)) << endl ;

    //check if rectangle and other geom intersects
    intersects = GEOSIntersects(street_rectangle , object_snapping_surface) ;
    cout <<"intersects : " << intersects <<endl ;
    //check distance to rectangle shell :
    GEOSDistance( GEOSGetExteriorRing(street_rectangle) , object_snapping_surface, &distance_to_shell);
    cout <<"distance_to_shell : " << distance_to_shell <<endl ;



        if(distance_to_shell== 0 ){
            //the object surface is halw inside, half outside of the road surface
            //the cost is shared surface or tot_surf - shared surf, depending on attractiv or repulsiv

            GEOSArea(GEOSIntersection(street_rectangle, object_snapping_surface)
                     , &shared_area) ;
            cout << "\t on the border" << endl;
            if(attractive==ATTRACTIVE){
                //cost is high when shared surface is low
                cost = object_snapping_surface_area - shared_area;
                cout << "\t \t attractive : tot-shared" << endl;
            }else{//cost is high when surface shared is high
                cost = shared_area;
                cout << "\t \t repulsive : shared" << endl;
            }

        }else{ // the object is either totally inside or totally outside
            cout << "\t strictly inside or outside" << endl;
            if(attractive==ATTRACTIVE){//attractive case : we want the surface to be inside the rectangle
                if(intersects==1){//the obj is inside the rectangle, perfect
                    cout << "\t \t attractive, inside, good" << endl;
                    cost = 0;
                }else{//the object is outside the rectangle, cost is high
                    cout << "\t \t attractive, outside, bad" << endl;
                    cost = (1+distance_to_shell) * object_snapping_surface_area;
                }
            }else{//repulsive case
                if(intersects==1){//the obj is inside the rectangle, cost is high
                    cout << "\t \t repuslive, inside, bad" << endl;
                    cost = (1+distance_to_shell) * object_snapping_surface_area;
                }else{//the obj is outside the rectangle , perfect
                    cout << "\t \t repuslive, outside, good" << endl;
                    cost = 0 ;
                }
            }
        }
        cout << "\t \t  cost : " << cost <<endl ;

//    finishGEOS();
    return cost ;

}





double test_geos(){


    std::string rec = "POLYGON((0 0.5, 5 0.5, 5 -0.5, 0 -0.5, 0 0.5))";
    std::string obj  = "POLYGON((1 3, 3 3, 3 4, 1 4, 1 3))";


    cout << "input string: " << rec << " ,  " << obj << endl ;
    geom rectangle;
    geom object ;

    GEOSWKTReader GEOS_DLL* reader ;
    reader = GEOSWKTReader_create();
    rectangle = GEOSWKTReader_read(reader, rec.c_str());
    object = GEOSWKTReader_read(reader, obj.c_str());
    GEOSWKTReader_destroy(reader);

    rectangle = read_WKT(rec);
    object = read_WKT(obj);

    cout << "having read the geoms " << endl ;
    cout << "rectangle :" <<  write_WKT(rectangle,2) <<endl ;
    cout << "object :" <<  write_WKT(object,2) <<endl ;

    double distance = 20 ;
    GEOSDistance( rectangle , object, &distance);

    cout << "distance : " << distance << endl ;


    return 0 ;
}

void initialize_geom_computation(){
    GEOSMessageHandler notice_function;
    GEOSMessageHandler error_function;
    initGEOS(notice_function,error_function);
}
void finish_geom_computation(){
    finishGEOS();
}


//bool read_WKT(std::string s){

//    GEOSContextHandle_t handle;
//    GEOSMessageHandler notice_function;
//    GEOSMessageHandler error_function;
//    GEOSWKTReader GEOS_DLL* reader ;
//    GEOSWKTWriter GEOS_DLL * writer ;

//    GEOSGeometry GEOS_DLL* parsed_geom;
//    GEOSGeometry GEOS_DLL* obj_buff;
//    GEOSGeometry GEOS_DLL* street_rectangle;
//    GEOSGeometry GEOS_DLL* intersection;

//    initGEOS(notice_function,error_function);
//    printf("begginnign of reading WKT \n");
//    std::cout << "wkt_geom : " << s << endl;


//    /// simulating the reading of an object (polygon) : we want to buffer it using the distance to border of 0.2

//    //parsing the geom
//    reader = GEOSWKTReader_create();
//    writer = GEOSWKTWriter_create_r(handle);
//    GEOSWKTWriter_setOutputDimension_r( handle, writer, 3);
//    GEOSWKTWriter_setTrim_r(handle,writer,1);
//    GEOSWKTWriter_setRoundingPrecision_r(handle,writer,10);


//    parsed_geom = GEOSWKTReader_read(reader, s.c_str());


//    //processing the geom
//    //buffer:

//    obj_buff = GEOSBufferWithStyle(parsed_geom
//                                   ,0.2 //width
//                                   ,1 //quadsegs
//                                   ,CAP_FLAT//endCapStyle
//                                   ,0 // int joinStyle
//                                   ,0.2//double mitreLimit
//                                   );


//    //parsed_geom = GEOSBuffer(parsed_geom,0.2 , 2);



//    /// A  : now we get the rectangle of the street axis dilated with street width
//    std::string street = "POLYGON((0 -4, 7 -4, 7 4, 0 4 , 0 -4))";
//    street_rectangle = GEOSWKTReader_read(reader, street.c_str());

//    /// C : computing area of intersection between object geom and street_rectangle geom

//    //only computing if the geom intersects :
//    if( GEOSIntersects(street_rectangle , obj_buff)==1){
//        cout << "the obejct and street_ractangle intersects, computing the shared surface" <<endl;
//        intersection = GEOSIntersection(street_rectangle, obj_buff);
//        cout << " \t intersection " << GEOSWKTWriter_write_r( handle, writer, intersection)<< endl;
//        double* area = new double();
//        int area_return  =  GEOSArea_r(handle,intersection, area) ;
//        cout << " \t area return : " <<  area_return << ", area computed : " << *area << endl;
//    }
//    else{//the two geom are not intersecting.
//        //computing the distance between the 2 geoms, this is the cost modulo area of object buff
//        cout << "the two geoms are not intersecting" << endl ;
//        cout << " \t cost = area(object) * (x 1+dist(object,erctangle) obj_buff: " ;

//    }
//    //writting the geom
//    //cleaning the reader



//    //cout <<  GEOSWKTWriter_write_r( handle, writer, obj_buff)<< endl;
//    GEOSWKTReader_destroy(reader);
//    GEOSWKTWriter_destroy_r(handle, writer);
//    //mandatory
//    finishGEOS();
//    cout << "end of reading_wkt"  << endl;
//    return false;
//}
