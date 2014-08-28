/*



template<typename T> inline
void CrossProduct(const T x[3], const T y[3], T x_cross_y[3]) {
  x_cross_y[0] = x[1] * y[2] - x[2] * y[1];
  x_cross_y[1] = x[2] * y[0] - x[0] * y[2];
  x_cross_y[2] = x[0] * y[1] - x[1] * y[0];
}
template<typename T> inline
void vect_addition(const T x[3], const T y[3], T x_plus_y[3]) {
  x_plus_y[0] = x[0] + y[0]  ;
  x_plus_y[1] = x[1] + y[1]  ;
  x_plus_y[2] = x[2] + y[2]  ;
}
template<typename T> inline
void vect_soustraction(const T x[3], const T y[3], T x_minus_y[3]) {
  x_minus_y[0] = x[0] - y[0]  ;
  x_minus_y[1] = x[1] - y[1]  ;
  x_minus_y[2] = x[2] - y[2]  ;
}
template<typename T> inline
void vect_squaredNorm(const T x[3], T x_squared) {
  x_squared  = x[0] * x[0] + x[1] * x[1] + x[2] * x[2]  ;
}



T vect_product[3];
T i_j[3] ;
vect_soustraction(n_i,n_j,i_j);
 T i_o[3];
 T obs[3]; obs[0] = T(position_[0]) ;  obs[1] = T(position_[1]) ;  obs[2] = T(position_[2]) ;

vect_soustraction(n_i,obs,i_o);

CrossProduct(i_j,i_o,vect_product);

T c_squared ;
vect_squaredNorm(vect_product, c_squared );

T ij_squared ;
vect_squaredNorm(i_j, ij_squared );


distance_to_axis[0] = c_squared/ij_squared ;


//distance computing function, this is templated.
template <typename T> bool operator()(const T* const n_i,
                                    const T* const n_j,
                                    T* distance_to_axis) const {

*/
/**
 * Here the type are expected to be :
 *  n_i,n_j : node
 *  distance_to_axis : a double/float
 *
 *   distance_to_axis : one of the quantities the solver will try to minimize
 *      this is the distance from the observation position
 *      to the observation projection on the line defined by n_i, n_j
 *      the distance is then norm( crossproduct(ON1,ON2))/norm(N1N2)
 **/
  /*

    //distance_to_axis[0] = (observation-n_i_vect).cross(observation-n_j_vect).squaredNorm()
    //        /  (n_i_vect-n_j_vect).squaredNorm();
return true;
}



*/
