#include "geometry_function.h"


using std::endl ;
using std::cout ;
using std::string;

geometry_function::geometry_function()
{
}


bool read_WKT(std::string s){

    GEOSContextHandle_t handle;
    GEOSMessageHandler notice_function;
    GEOSMessageHandler error_function;
    GEOSWKTReader GEOS_DLL* reader ;
    GEOSWKTWriter GEOS_DLL * writer ;

    GEOSGeometry GEOS_DLL* parsed_geom;
    GEOSGeometry GEOS_DLL* buff_geom;

    initGEOS(notice_function,error_function);
    printf("begginnign of reading WKT \n");
    std::cout << "wkt_geom : " << s << endl;


    /// simulating the reading of an object (polygon) : we want to buffer it using the distance to border of 0.2

    //parsing the geom
    reader = GEOSWKTReader_create();
    parsed_geom = GEOSWKTReader_read(reader, s.c_str());
    GEOSWKTReader_destroy(reader);

    //processing the geom
    //buffer:

    buff_geom = GEOSBufferWithStyle(parsed_geom
                                  ,0.2 //width
                                  ,1 //quadsegs
                                  ,CAP_FLAT//endCapStyle
                                  ,0 // int joinStyle
                                  ,0.2//double mitreLimit
                                  );


    //parsed_geom = GEOSBuffer(parsed_geom,0.2 , 2);






   //writting the geom
   writer = GEOSWKTWriter_create_r(handle);
   GEOSWKTWriter_setOutputDimension_r( handle, writer, 3);
   GEOSWKTWriter_setTrim_r(handle,writer,1);
   GEOSWKTWriter_setRoundingPrecision_r(handle,writer,10);

   cout <<  GEOSWKTWriter_write_r( handle, writer, buff_geom)<< endl;
    GEOSWKTWriter_destroy_r(handle, writer);
    //mandatory
    finishGEOS();
    cout << "end of reading_wkt"  << endl;
    return false;
}
