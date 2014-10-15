#ifndef GEOMETRY_FUNCTION_H
#define GEOMETRY_FUNCTION_H

#include "geos_c.h"

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "time_measurement.h"
#include "Constraints.h"




//enum In_Out_Border{ATTRACTIVE=1 ,REPULSIVE=-1 } ;

/// type to abstract geos lib for the rest of the application
typedef GEOSGeometry GEOS_DLL* geom;

/// read a string containing wkt geom. Returns a geom
geom read_WKT(std::string );

/// write a geom into wkt . Don't forget to free the return after !
char *write_WKT(geom, int dim);

/// free memory of char buffer allocated by writting
//void free_writter_buffer(void *buffer){
//    GEOSFree( buffer);
//}

/** function that compute a rectangle given 2 points and a width.
    The 2 points are a segment, the rectangle is obtained by dilating segment by width
*/
void axis_to_rectangle(const double* pt1,const double* pt2, double axis_width, geom rectangle);

/** Main function to compute cost that is shared area if overlaps, area*(1+dist) if no overlaps
  this geometric function takes a segment (2 points), this segment width? It computes a rectangle
  then it computes the shared area between this rectangle and the object_surface
    If this shared area is 0, it returns the distance between object_surface and rectangle.
    Else, it returns the shared area
  */
double shared_area_cost( road_relation_enum road_relation , const double* pt1, const double* pt2, double axis_width, geom object_snapping_surface, double object_snapping_surface_area);



double test_geos();

void initialize_geom_computation();
void finish_geom_computation() ;


class geometry_function
{
public:
    geometry_function();
    //bool read_WKT(string);
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
