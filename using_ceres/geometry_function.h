#ifndef GEOMETRY_FUNCTION_H
#define GEOMETRY_FUNCTION_H

#define TOLERANCE 0.0004
#define TOLERANCY_EQUAL(A,B) ( std::abs((A)-(B))<(TOLERANCE))
#define SIGN(X) (((X) > 0) - ((X) < 0))

#include "geos_c.h"

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
//#include <algorithm> //for max


#include "eigen3/Eigen/Eigen"
typedef Eigen::Map<Eigen::Vector3d> VectorRef;
typedef Eigen::Map<const Eigen::Vector3d> ConstVectorRef;


#include "time_measurement.h"
#include "enum_functions.h"
//#include "Constraints.h"
/// type to abstract geos lib for the rest of the application
typedef GEOSGeometry* GEOS_DLL geometry;



//enum In_Out_Border{ATTRACTIVE=1 ,REPULSIVE=-1 } ;



/// read a string containing wkt geom. Returns a geom
geometry read_WKT(std::string );

/// write a geom into wkt . Don't forget to free the return after !
char *write_WKT(geometry, int dim);

/// free memory of char buffer allocated by writting
//void free_writter_buffer(void *buffer){
//    GEOSFree( buffer);
//}

/** function that compute a rectangle given 2 points and a width.
    The 2 points are a segment, the rectangle is obtained by dilating segment by width
*/
void axis_to_rectangle(const double* pt1,const double* pt2, double axis_width, geometry * axis_to_be_filled, geometry * rectangle);

double geometry_width_regarding_axis(const double* pt1,const double* pt2 , geometry object, geometry centroid);


/** Main function to compute cost that is shared area if overlaps, area*(1+dist) if no overlaps
  this geometric function takes a segment (2 points), this segment width? It computes a rectangle
  then it computes the shared area between this rectangle and the object_surface
    If this shared area is 0, it returns the distance between object_surface and rectangle.
    Else, it returns the shared area
  */
double shared_area_cost( SnapEnums::road_relation_enum road_relation , const double* pt1, const double* pt2, double axis_width, geometry object_snapping_surface, double object_snapping_surface_area, geometry centroid);
double signed_dist_to_border( SnapEnums::road_relation_enum road_relation , const double* pt1, const double* pt2, double axis_width, geometry object_snapping_surface, double object_snapping_surface_area);


void initialize_geom_computation();
void finish_geom_computation() ;


class Geometry
{
public:
    Geometry();
    //bool read_WKT(string);

    static geometry BufferWithStyle(const geometry g,
                                    double width, int quadsegs, int endCapStyle, int joinStyle,
                                    double mitreLimit){
        return   GEOSBufferWithStyle( g,
                                      width, quadsegs, endCapStyle,  joinStyle,
                                      mitreLimit);
    }
    static double area(const geometry g){
        double a;
        GEOSArea(g, &a) ;
        return a;
    }

    static geometry centroid(geometry g){
        return GEOSGetCentroid(g);
    }
    static const double * return_double(int dim, const double * coor){
        return coor;
    }
    static int  coorPoint2Double(const GEOSCoordSequence GEOS_DLL * s, double * coordinates){
        //input must be a point
        //convert point ot coordinate sequence

        unsigned int dims=4;
        //filling dims
        GEOSCoordSeq_getDimensions(s,&dims) ;

        double * temp= new double[dims];
        // Initialize all elements to zero.
        for (unsigned int i=0; i<dims; ++i) {
            temp[i] = 0;
        }
        for(unsigned int i=0;i<dims;i++){//writing x, y, z
            GEOSCoordSeq_getOrdinate(s,0,i,&temp[i]);
        }
        //std::cout << temp[0] << " " << temp[1] << std::endl;
        //std::cout  << dims << std::endl ;

        for (unsigned int i=0; i<dims; ++i) {
            coordinates[i] = temp[i];
        }
        return dims;
    }
    static int  geomPoint2Double(const geometry g, double * coordinates){
        const GEOSCoordSequence GEOS_DLL * s = GEOSGeom_getCoordSeq( g);
        return coorPoint2Double(  s,  coordinates) ;
    }
    static geometry double2geomPoint(const double * coordinates ){

        GEOSCoordSequence * s=  GEOSCoordSeq_create( 1,3);
        GEOSCoordSeq_setX(s,0,coordinates[0]);
        GEOSCoordSeq_setY(s,0,coordinates[1]);
        GEOSCoordSeq_setZ(s,0,coordinates[2]);
        return  GEOSGeom_createPoint( s);
    }

    static int orientationIndex(const double * A,const  double * B,const double * P)  {
        return GEOSOrientationIndex(
                    A[0],A[1]
                    , B[0],B[1]
                    , P[0],P[1]);
    }

    static int closestPoint( geometry g1,  geometry g2,double * closest_on_g1){
        GEOSCoordSequence * s= GEOSNearestPoints( g1,  g2);
        coorPoint2Double( s,  closest_on_g1) ;
        return 1;
    }



};








enum EndCapStyle {

    /// Specifies a round line buffer end cap style.
    CAP_ROUND=1,

    /// Specifies a flat line buffer end cap style.
    CAP_FLAT=2,

    /// Specifies a square line buffer end cap style.
    CAP_SQUARE=3
};

/// Join styles
enum JoinStyle {

    /// Specifies a round join style.
    JOIN_ROUND=1,

    /// Specifies a mitre join style.
    JOIN_MITRE=2,

    /// Specifies a bevel join style.
    JOIN_BEVEL=3
};


#endif // GEOMETRY_FUNCTION_H
