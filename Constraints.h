#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H
//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////
/**
  This file defines the constaint (functor for cost)
  relevant to the problem
  */

extern const double K_origin;
extern const double K_obs;
extern const double K_spacing;

/** functor to compute cost between one observation and the segment between 2 nodes
  */
struct DistanceToProjectionResidual {
   //! this is the constructor, it expects an array of at least 3 doubles.
    DistanceToProjectionResidual(const double* input_vect)
        :position_(input_vect) {}

    //! this is the operator computing the cost, ie the distance projeted on normal of (n_i,n_j)
        template <typename T> bool operator()(const T* const n_i,/**< the first node */
                                              const T* const n_j,/**< the second node */
                                            T* distance_to_original_spacing) const {/**< this is the cost to be optimised*/

        //computting soustraction
        T obs[3] = {T(position_[0]),T(position_[1]),T(position_[2])};
        T n_i_minus_obs[3] ;
        T n_i_minus_n_j[3] ;
        soustraction(n_i,obs,n_i_minus_obs);
        soustraction(n_i,n_j,n_i_minus_n_j);
        T cross[3];
        ceres::CrossProduct(n_i_minus_obs,n_i_minus_n_j,cross);
        distance_to_original_spacing[0] =   squaredNorm(cross)/squaredNorm(n_i_minus_n_j) -T(1.0)/*width of this edge*/ ;

        return true;
      }
 private:
    const double* position_; /**< store the 3D position of the observation point */
};


/** functor to compute cost between one node position and this node original position
  */

struct DistanceToInitialPosition {
   //! this is the constructor, it expects an array of 3 doubles to fill the original position
    DistanceToInitialPosition(const double* input_vect)
        :initial_position_(input_vect) {
    }

    //! distance computing function (eucl distance), this is templated.
    template <typename T> bool operator()(const T* const n_i,
                                        T* distance_to_origin) const {
    T s[3] ;
    s[0] =n_i[0]-T(initial_position_[0]);
    s[1] =n_i[1]-T(initial_position_[1]);
    s[2] =n_i[2]-T(initial_position_[2]);
    distance_to_origin[0] = K_origin *  ceres::sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2] + 0.000001);

    return true;
  }
 private:
    const double* initial_position_;
};



/** functor to compute cost between 2 nodes ; they should keep about the same direction/distance
  */
struct DistanceToInitialSpacing{
   //! this is the constructor, it expects an array of at least 3 doubles.
    DistanceToInitialSpacing(const double* input_vect)
        :initial_spacing_(input_vect) {}

    //! this is the operator computing the cost, that is the difference to the original spacing
    template <typename T> bool operator()(const T* const n_i,/**< the first node */
                                          const T* const n_j,/**< the second node */
                                        T* distance_to_original_spacing) const {/**< this is the cost to be optimised*/

        //compute the difference with original spacing:
        distance_to_original_spacing[0] = K_spacing * T(initial_spacing_[0]) - (n_i[0] - n_j[0])  ;
        distance_to_original_spacing[1] = K_spacing * T(initial_spacing_[1]) - (n_i[1] - n_j[1])  ;
        distance_to_original_spacing[2] = K_spacing * T(initial_spacing_[2]) - (n_i[2] - n_j[2])  ;

        return true;
      }
 private:
    const double* initial_spacing_; /**< store the original 3D vector between the 2 nodes*/
};


#endif // CONSTRAINTS_H
