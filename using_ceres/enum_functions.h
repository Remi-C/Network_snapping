#ifndef ENUM_FUNCTIONS_H
#define ENUM_FUNCTIONS_H

#include <string>
#include <iostream>

template <typename T>
struct EnumString {

    typedef T enum_type;

    EnumString(T val, const std::string& str)
        : _val(val), _str(str) {
    }

    T _val;
    const std::string _str;
};

template <typename T>
EnumString<T> make_enum_string(T val, const std::string& str)
{
    return EnumString<T>(val, str);
}


template< class TEnumMapping >
typename TEnumMapping::enum_type ToEnum(TEnumMapping enumStrings[], const std::string name)
{
    TEnumMapping* ptr = enumStrings;
    unsigned int i = 0 ;

    while (!ptr->_str.empty())
    {
        if( ! ptr->_str.compare(name ) )
            break ;
        ++i ;
        ++ptr ;
    }
    return enumStrings[i]._val;
}

template< typename TEnumMapping >
std::string ToString(TEnumMapping enumStrings[], const typename TEnumMapping::enum_type enu)
{
    TEnumMapping* ptr = enumStrings;
    unsigned int i = 0 ;

    while(i<100)// (!ptr->_str.empty())
    {
        if(  ptr->_val == enu )
            break ;
        ++i ;
        ++ptr ;
    }
    return enumStrings[i]._str;
}





class SnapEnums
{
public:
    enum road_relation_enum{IN=1 ,OUT=-1 ,BORDER=0, BORDER_IN = 10, BORDER_OUT= -10, UNDEF=-110 } ;
    enum geom_type_enum{POINT=1,LINESTRING=2,POLYGON=3,COLLECTION=4} ;
    enum attractive_repulsive{ATTRACTIVE=1 ,REPULSIVE=-1,ATTR_AND_REP =-11,  NEUTRAL=0 } ;

    static EnumString<road_relation_enum> road_relation_string[];
    static EnumString<geom_type_enum> geom_type_string[];
    static EnumString<attractive_repulsive> attractive_repulsive_string[];



//    static road_relation_enum  string_to_road_relation_enum(std::string s);
//    static std::string road_relation_enum_to_string(road_relation_enum e);
//    static geom_type_enum string_to_geom_type_enum(std::string s);
//    static std::string geom_type_enum_to_string(geom_type_enum e);

    static void test_enum();

    static std::string rre_toString(road_relation_enum);
    static std::string gt_toString(geom_type_enum);
    static std::string ar_toString(attractive_repulsive);

    static road_relation_enum String_torre(std::string s);
    static geom_type_enum String_togt(std::string s);
    static attractive_repulsive String_toar(std::string s);
};




#endif // ENUM_FUNCTIONS_H
