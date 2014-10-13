#ifndef GEOMETRY_FUNCTION_H
#define GEOMETRY_FUNCTION_H

#include "geos_c.h"

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstdlib>

#include "time_measurement.h"

bool read_WKT(std::string s);


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
