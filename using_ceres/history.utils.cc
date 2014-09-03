
/**
  Abstract : snapping a street axis to observations

  _ Original Data :
    street axis, along with width information
    observations of sidewalk, along with confidence and weight information

  _ Input data :
    street axis node (n_i)
    connectivity information, along with width w_i_j of the part of the axis (e_i_j)

    observations of sidewalk, related to 2 nodes ; points o_k(n_i, n_j )

  _ optimisation model data
    n_i(X,Y,Z)
    w_i_j(n_i, n_j)

  _ optimisation constraints (K are normalised weight,Sum(K)=1, so that we can parametrize importance of constraints)
    _ K_n * intial position of n_i
    _ K_n_n * initial direction between n_i, n_j
    _ K_w * initial value of w_i_j
    _ K_d * distance of o_k to line (n_i,n_j)

  _ optimisation cost function
    _ a straightforward euclidian distance function (we can use the squared euclidian distance for faster computation) d_eucl(R_3,R_3)->R
    _ a

_ optimisation results :

  */
    /*
  simple file to put together all the idea about how to use ceres

  */

/*


  */

//////// creating the variables to be optimized : setting ParametterBlock

    /*
      //from bundle_adjuster.cc

      ceres::ParameterBlockOrdering* ordering =
          new ceres::ParameterBlockOrdering;

          //add the position of p_i
            for (int i = 0; i < num_points; ++i) {
              ordering->AddElementToGroup(points + point_block_size * i, 0);

          //add the width of the axis
            for (int i = 0; i < num_edges; ++i) {
              ordering->AddElementToGroup( 1 * i, 0);


       options->linear_solver_ordering.reset(ordering);
    }
    */


/*

  ///////////////PSEUDO CODE  :

  _ Structure of the model
    _ define Loss function
          // Configure the loss function.
          LossFunction* loss = NULL;
          if (FLAGS_robust_threshold) {
            loss = new CauchyLoss(FLAGS_robust_threshold);
          }


    _ Functor : class performing the error computation in a function and old some values
        create functor

    _ Cost function
        create cost function using functor



    add residual block (function on parameters)

  */


///////// creating a fist constraint on n_i_j : it shouldn't go too far away from it's original position
    /*


    */

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



        T soustraction[3] ;
        T i_pos[3] ;
        i_pos[0] = T(initial_position_[0]) ;
        i_pos[1] = T(initial_position_[1]) ;
        i_pos[2] = T(initial_position_[2]) ;
        //i_pos =  (T) initial_position_ ;
        vect_soustraction(n_i,i_pos, soustraction) ;
        distance_to_axis[0] = vect_squaredNorm(soustraction) ;




    template<typename T> inline
    void vect_soustraction(const T x[3], T y[3], T x_minus_y[3]) {
      x_minus_y[0] = x[0] - y[0]  ;
      x_minus_y[1] = x[1] - y[1]  ;
      x_minus_y[2] = x[2] - y[2]  ;
    }
    template<typename T> inline
    T vect_squaredNorm(const T x[3] ) {
      return x[0] * x[0] + x[1] * x[1] + x[2] * x[2]  ;
    }


*/
