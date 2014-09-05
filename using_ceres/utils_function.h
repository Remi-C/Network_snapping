#ifndef UTILS_FUNCTION_H
#define UTILS_FUNCTION_H
//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////


//! fucntion to perform soustraction on double[3] vector
template<typename T> inline
void soustraction(const T x[3], const T y[3], T result[3]) {
  result[0] = x[0]- y[0];
  result[1] = x[1] -y[1];
  result[2] = x[2] -y[2];
}

//! fucntion to perform squared norm evaluation
template<typename T> inline
T squaredNorm(const T x[3] ) {
  return (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}



#endif //UTILS_FUNCTION_H
