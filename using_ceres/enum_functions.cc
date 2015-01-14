#include "enum_functions.h"

using std::string;
using std::cout;
using std::endl;

EnumString<SnapEnums::road_relation_enum>   SnapEnums::road_relation_string[]  =
{
    make_enum_string(IN, "IN"),
    make_enum_string(OUT, "OUT"),
    make_enum_string(BORDER, "BORDER"),
    make_enum_string(BORDER_IN, "BORDER_IN"),
    make_enum_string(BORDER_OUT, "BORDER_OUT"),
    make_enum_string(UNDEF, "UNDEF"),
    make_enum_string(UNDEF, "") //the last string must be empty to stop the search
};

EnumString<SnapEnums::geom_type_enum>   SnapEnums::geom_type_string[]  =
{
    //POINT=1,LINESTRING=2,POLYGON=3,COLLECTION=4
    make_enum_string(POINT, "POINT"),
    make_enum_string(LINESTRING, "LINESTRING"),
    make_enum_string(POLYGON, "POLYGON"),
    make_enum_string(COLLECTION, "COLLECTION"),
    make_enum_string(COLLECTION, "") //the last string must be empty to stop the search
};


EnumString<SnapEnums::attractive_repulsive>   SnapEnums::attractive_repulsive_string[]  =
{
    //ATTRACTIVE=1 ,REPULSIVE=-1,ATTR_AND_REP =-11,  NEUTRAL=0
    make_enum_string(ATTRACTIVE, "ATTRACTIVE"),
    make_enum_string(REPULSIVE, "REPULSIVE"),
    make_enum_string(ATTR_AND_REP, "ATTR_AND_REP"),
    make_enum_string(NEUTRAL, "NEUTRAL"),
    make_enum_string(NEUTRAL, "") //the last string must be empty to stop the search
};

EnumString<SnapEnums::optimisation_target>   SnapEnums::optimisation_target_string[]  =
{
    //ATTRACTIVE=1 ,REPULSIVE=-1,ATTR_AND_REP =-11,  NEUTRAL=0
    make_enum_string(POSITION, "POSITION"),
    make_enum_string(WIDTH, "WIDTH"),
    make_enum_string(MIXED, "MIXED"),
    make_enum_string(OTHER, "") //the last string must be empty to stop the search
};

//SnapEnums::road_relation_enum SnapEnums::string_to_road_relation_enum(string s){
//    if(s.compare("IN")==0)
//        return IN;
//    if(s.compare("OUT")==0)
//        return OUT;
//    if(s.compare("BORDER")==0)
//        return BORDER;
//    if(s.compare("BORDER_IN")==0)
//        return BORDER_IN;
//    if(s.compare("BORDER_OUT")==0)
//        return BORDER_OUT;
//    if(s.compare("UNDEF")==0)
//        return UNDEF;
//}
//string SnapEnums::road_relation_enum_to_string(road_relation_enum e){
//    if(e==IN)
//        return "IN";
//    if(e==OUT)
//        return "OUT";
//    if(e==BORDER)
//        return "BORDER";
//    if(e==BORDER_IN)
//        return "BORDER_IN";
//    if(e==BORDER_OUT)
//        return "BORDER_OUT";
//    if(e==UNDEF)
//        return "UNDEF";
//    return std::to_string( e );
//}


//SnapEnums::geom_type_enum SnapEnums::string_to_geom_type_enum(string s){
//    geom_type_enum e;
//    if(s.compare("POINT")==0){
//        e= POINT;}
//    if(s.compare("LINESTRING")==0){
//        e= LINESTRING;}
//    if(s.compare("POLYGON")==0){
//        e= POLYGON;}
//    if(s.compare("COLLECTION")==0){
//        e= COLLECTION;}
//    return e;
//}
//string SnapEnums::geom_type_enum_to_string(geom_type_enum e){
//    if(e==POINT)
//        return "POINT";
//    if(e==LINESTRING)
//        return "LINESTRING";
//    if(e==POLYGON)
//        return "POLYGON";
//    if(e==COLLECTION)
//        return "COLLECTION";
//    return std::to_string(e);
//}

std::string SnapEnums::rre_toString(SnapEnums::road_relation_enum e){
    return ToString(SnapEnums::road_relation_string , e);
}
std::string SnapEnums::gt_toString(SnapEnums::geom_type_enum e){
    return ToString(SnapEnums::geom_type_string , e);
}
std::string SnapEnums::ar_toString(SnapEnums::attractive_repulsive e){
    return ToString(SnapEnums::attractive_repulsive_string , e);
}
std::string SnapEnums::ot_toString(SnapEnums::optimisation_target e){
    return ToString(SnapEnums::optimisation_target_string , e);
}

SnapEnums::road_relation_enum SnapEnums::String_torre(std::string s){
    return ToEnum(SnapEnums::road_relation_string , s);
}
SnapEnums::geom_type_enum SnapEnums::String_togt(std::string s){
    return ToEnum(SnapEnums::geom_type_string , s);
}
SnapEnums::attractive_repulsive SnapEnums::String_toar(std::string s){
    return ToEnum(SnapEnums::attractive_repulsive_string , s);
}
SnapEnums::optimisation_target SnapEnums::String_toot(std::string s){
    return ToEnum(SnapEnums::optimisation_target_string , s);
}



void SnapEnums::test_enum(){
    //testing the enum
    SnapEnums::road_relation_enum rre = SnapEnums::BORDER_OUT ;
    string readed_param = "BORDER_OUT";
    std::cout << ToEnum( SnapEnums::road_relation_string , readed_param) <<std::endl;
    std::cout << ToString(SnapEnums::road_relation_string , rre) << endl;
    string readed_param_2 = "REPULSIVE";

    cout << ToString(attractive_repulsive_string,ToEnum(attractive_repulsive_string,readed_param_2)) <<endl;
    cout << ToString(geom_type_string,ToEnum(geom_type_string,readed_param_2)) <<endl;


    cout << SnapEnums::rre_toString(SnapEnums::BORDER_OUT)<<endl ;
    cout << SnapEnums::gt_toString(SnapEnums::POINT)<<endl ;
    cout << SnapEnums::ar_toString(SnapEnums::ATTRACTIVE)<<endl ;
    cout << SnapEnums::ot_toString(SnapEnums::POSITION)<<endl  ;

    cout << SnapEnums::rre_toString(SnapEnums::String_torre("BORDER_IN")) << endl;
    cout << SnapEnums::gt_toString(SnapEnums::String_togt("LINESTRING")) << endl;
    cout << SnapEnums::ar_toString(SnapEnums::String_toar("REPULSIVE")) << endl;
    cout << SnapEnums::ot_toString(SnapEnums::String_toot("POSITION")) << endl;
}


