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

#include "Data.h"

//extern const double K_origin;
//extern const double K_obs;
//extern const double K_spacing;


#include "ceres/ceres.h"
#include "ceres/rotation.h"

typedef Eigen::Map<Eigen::Vector3d> VectorRef;
typedef Eigen::Map<const Eigen::Vector3d> ConstVectorRef;

#include "Parameters.h"
#include "utils_function.h"

using std::cout;
using std::endl;

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::LossFunction;



int addConstraintsOnInitialPosition(DataStorage *, ceres::Problem * );
int addConstraintsOnInitialspacing( DataStorage *, ceres::Problem * );
int addConstraintsOnOrthDistToObservation(DataStorage * , ceres::Problem * ) ;
int addManualConstraintsOnOrthDistToObservation(DataStorage * , Problem * ) ;
int addManualConstraintsOnDistanceToOriginalAngle(DataStorage * , Problem * ) ;
int addManualConstraintsOnInitialspacing(DataStorage * , Problem * ) ;

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
    s[0] = (n_i[0]-T(initial_position_[0]) );
    s[1] = (n_i[1]-T(initial_position_[1]) );
    s[2] = (n_i[2]-T(initial_position_[2]) );
    //distance_to_origin[0] =  squaredNorm(s) ;
    distance_to_origin[0] =  ceres::sqrt(squaredNorm(s)+0.0001) ;
    //distance_to_origin[0] = T(K_origin) *  ceres::sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2] + 0.000001);

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
        T spac[3] = {T(initial_spacing_[0]),T(initial_spacing_[1]),T(initial_spacing_[2])};
        T n_i_minus_n_j[3] ;
        soustraction(n_i,n_j,n_i_minus_n_j);

         distance_to_original_spacing[0]=  ceres::pow(squaredNorm(n_i_minus_n_j) - squaredNorm(spac),2)   ;
//        //compute the difference with original spacing:
//        distance_to_original_spacing[0] = T(K_spacing) * ( T(initial_spacing_[0]) - (n_i[0] - n_j[0]) ) ;
//        distance_to_original_spacing[1] = T(K_spacing) * ( T(initial_spacing_[1]) - (n_i[1] - n_j[1]) );
//        distance_to_original_spacing[2] = T(K_spacing) * ( T(initial_spacing_[2]) - (n_i[2] - n_j[2]) );

        return true;
      }
 private:
    const double* initial_spacing_; /**< store the original 3D vector between the 2 nodes*/
};


/** functor to compute cost between one observation and the segment between 2 nodes
  */
struct DistanceToProjectionResidual {
   //! this is the constructor, it expects an array of at least 3 doubles.
    DistanceToProjectionResidual(const double* input_vect,const double * input_w, const observation * input_obs)
        :position_(input_vect), w_i_j_(input_w), obs_(input_obs) {}

    //! this is the operator computing the cost, ie the distance projeted on normal of (n_i,n_j)
        template <typename T> bool operator()(const T* const n_i,/**< the first node */
                                              const T* const n_j,/**< the second node */
                                            T* distance_to_proj) const {/**< this is the cost to be optimised*/

        //computting soustraction
        T obs[3] = {T(position_[0]),T(position_[1]),T(position_[2])};
        T obs_minus_n_i[3] ;
        T obs_minus_n_j[3] ;
        T n_i_minus_n_j[3] ;
        soustraction(obs,n_i,obs_minus_n_i);
        soustraction(obs,n_j,obs_minus_n_j);
        soustraction(n_i,n_j,n_i_minus_n_j);
        T cross[3];
        ceres::CrossProduct(obs_minus_n_i,obs_minus_n_j,cross);


//        distance_to_proj[0] =   T(K_obs) * ceres::sqrt(squaredNorm(n_i_minus_obs)) - T(w_i_j_[0])/2.0 ;
//        distance_to_proj[1] =   T(K_obs) * ceres::sqrt(squaredNorm(n_j_minus_obs)) -T(w_i_j_[0])/2.0 ;

//        distance_to_proj[0] =   1.0 / T(2.0) * ( ceres::sqrt(squaredNorm(n_i_minus_obs)) - T(w_i_j_[0])/2.0  )  ;
//        distance_to_proj[1] =  1.0 /T(2.0) * ( ceres::sqrt(squaredNorm(n_j_minus_obs)) -T(w_i_j_[0])/2.0)  ;


//        distance_to_proj[0] =
//                ceres::sqrt( squaredNorm(n_i_minus_obs) + squaredNorm(n_j_minus_obs) +0.00001)  ;
//        distance_to_proj[0] =
//                    T(K_obs) * T(obs_->confidence) * T(obs_->weight) *
//                (  squaredNorm(cross)/squaredNorm(n_i_minus_n_j)
//                 -T(1)/T(2.0)) ;


        distance_to_proj[0] =
                     T(obs_->confidence) * T(obs_->weight) *
                (
                 ceres::sqrt(squaredNorm(cross)/squaredNorm(n_i_minus_n_j) )
                 -T(w_i_j_[0])/2.0
                ) ;

//        distance_to_proj[0] =
//                ceres::sqrt(
//                    T(K_obs) * T(obs_->confidence) * T(obs_->weight) *
//                    (cross[0] * cross[0]/ squaredNorm(n_i_minus_n_j)-T(w_i_j_[0])/T(2.0) )
//                 +T(0.00001)) ;
//        distance_to_proj[1] =
//                ceres::sqrt(T(K_obs) * T(obs_->confidence) * T(obs_->weight)* (cross[1]* cross[1] / squaredNorm(n_i_minus_n_j)-T(w_i_j_[0])/T(2.0) )+T(0.00001)) ;
//        distance_to_proj[2] =
//                ceres::sqrt( T(K_obs) * T(obs_->confidence) * T(obs_->weight) * ( cross[2]* cross[2] / squaredNorm(n_i_minus_n_j)-T(w_i_j_[0])/T(2.0) )+ T(0.000001)) ;

        return true;
      }
 private:
    const double* position_; /**< store the 3D position of the observation point */
    const double* w_i_j_;//! the width of the edge, in meter
    const observation * obs_;//link to the observation, for easy use of confidence and weight.
};


//! trying to replace the auto computation of Jacobian by a custom one
class ManualOrthDistanceToObservation  : public ceres::SizedCostFunction<1,3,3> {

public :
   //! this is the constructor, it expects an array of at least 3 doubles.
    ManualOrthDistanceToObservation(const double* input_vect,const double * input_w, const observation * input_obs)
        :position_(input_vect), w_i_j_(input_w), obs_(input_obs) {}

    //! this is the operator computing the cost, ie the distance projeted on normal of (n_i,n_j)
    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const {
        //the parameters is as follow : parameter[0-2] = n_i = first node;parameter[3-5] = n_j = second node;

        //map the input array into 2 eigen vectors, plus map observation position into Eigen
        //cout << "\nbeginning of evaluate" << endl;
        ConstVectorRef Ni( parameters[0],3 );
        ConstVectorRef Nj( parameters[1],3 );
        ConstVectorRef Ob(position_,3);

        //compute the normal of the plan defined by vect(NiO, NiNj)/norm(...) = n
        Eigen::Vector3d Np = (Ob-Ni).cross(Nj-Ni) ;
        //Np = Np.norm();
        //compute director vector of (NiNj) ie : NiNj/norm(NiNj) = u
        Eigen::Vector3d U = (Nj-Ni)/(Nj-Ni).norm();
        //compute residual = distance from O to NiNj : norm(vect(NiO,NiNj))/norm(NiNj)
        double d = Np.norm()/(Nj-Ni).norm()-  w_i_j_[0]/2.0;
        residuals[0]= pow(d,2) *   obs_->confidence * obs_->weight ;
        //compute Jacobian director vector : vect(u,n)
        Eigen::Vector3d Vja = U.cross(Np/Np.norm());
        //if(obs_->obs_id ==1 ) {Vja << -0.7,-0.7,0;}
        //else{Vja << +0.7,+0.7,0;}

        //compute the direction of movement : - = toward the obs, + = away from point
        int sign = ((  residuals[0]  >0) - (residuals[0] <0));
        //compute Jacobian norm for Ni : for test simply take d
        Eigen::Vector3d Ji = -1 * sign* Vja * d ;
        //compute Jacobian norm for Nj : for test simply take d
        Eigen::Vector3d Jj = -1 * sign * Vja *  d ;/// @warning : remove this 10 factor aspa !

//        cout << "  Observation_id : " <<  obs_->obs_id <<std::endl;
//         cout << "  Ni : " << Ni.transpose() <<std::endl;
//        cout << " Nj : " << Nj.transpose() <<std::endl;
//        cout << " Vja : " << Vja.transpose() <<std::endl;
//        cout << "  distance : " << residuals[0] <<std::endl;
//        cout << "   Ji :" << Ji.transpose() <<endl;
//        cout << "   Jj :" << Jj.transpose() <<endl;
        // std::cout << "\njac (eigen): \n" << jac << std::endl;

        if (jacobians == NULL) {
            //cout << "JACOBIAN NULL" <<endl;
          return 1;
        }

         if (jacobians != NULL && jacobians[0] != NULL) {
             //note: null jacobian means end of computation?
            jacobians[0][0] =  Ji(0);
            jacobians[0][1] =  Ji(1);
            jacobians[0][2]=   Ji(2);

        }
         if (jacobians != NULL && jacobians[1] != NULL) {
             //note: null jacobian means end of computation?
            jacobians[1][0] =  Jj(0);
            jacobians[1][1] =  Jj(1);
            jacobians[1][2]=   Jj(2);

        }

       //  cout<<"end of evaluate()\n";
        return true;
      }
 private:
    const double* position_; /**< store the 3D position of the observation point */
    const double* w_i_j_;//! the width of the edge, in meter
    const observation * obs_;//link to the observation, for easy use of confidence and weight.
};



/** this class compute the angle formed by 2 edges using a node, and measure the difference with the original angle
This is a regularisation class.
We return a 2D cost : sin(angle) and cos(angle), so has to avoid pi/2 roation invariance.

We expect the parameters to be in a fixed order: node with angle constraint, a node, another node
*/
class ManualDistanceToOriginalAngle  : public ceres::SizedCostFunction<2,3,3,3> {
public :
   //! this is the constructor, it expects 2 values : vect_1.vect_2/(norm(vect_1)*norm(vect_2)) and vect_1xvect_2/(norm(vect_1)*norm(vect_2))
    ManualDistanceToOriginalAngle(const double i_scalar_angle,const double i_cross_angle)
        : scalar_angle(i_scalar_angle) , cross_angle(i_cross_angle) {}



    //! this is the operator computing the costs, ie the distance of scalar and cross angle with the original values
    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const {
        //! @param parameters[0] : node on which to compute the constraint (center node)
        //! @param parameters[1] : one of the node forming the angle
        //! @param parameters[2] : other node forming the angle


//        cout << "begining of evaluate" <<endl ;
        //map the input array into 3 eigen vectors
        ConstVectorRef Nc( parameters[0],3 );
        ConstVectorRef Ni( parameters[1],3 );
        ConstVectorRef Nj( parameters[2],3 );

        double scalar_a = (Nc-Ni).dot(Nc-Nj)/((Nc-Ni).norm() * (Nc-Nj).norm());
        double cross_a = ((Nc-Ni).cross(Nc-Nj)/((Nc-Ni).norm() * (Nc-Nj).norm())).norm();


        //compute residuals (= cost)
        residuals[0] = pow(scalar_angle-scalar_a,2) ;
        residuals[1] = pow(cross_angle-cross_a,2) ;

        //compute jacobian :
        //only the center node (Nc) should be moved
        //the direction is along the bissect of angle (Ni,Nc,Nj)
        Eigen::Vector3d Vjc = ((Nc-Ni).normalized() + (Nc-Nj).normalized() ).normalized();

        /** value of displacement :
          lets consider a triangle formed on angle alpha (Ni,Nc,Nj) (upper summit is Nc)
          but with Ni,Nc and Nc,Nj normalized to 1.
          we want to find Nc_p, that is the new position of Nc after a move of d along the bissect.
          the final expected angle (Ni,Nc',Nj) is the original angle beta (caracterized by sin and cos)
          the distance d = sin(alpha)*tan(beta)-cos(alpha)
          Note that sign determin which direction it goes
          Note : in theory every angle shoud be divided by 2,thus we may have ot use sin2x=f(sinx^2, sinx) to have correct result
          */
       double d = cross_a * cross_angle/scalar_angle  -scalar_a ;
       double init[9] = {0,0,0,0,0,0,0,0,0};

       for(int i=0; i<3;++i){
            for(int j=0; j<6 ; ++j){
                jacobians[i][j] =0 ;
//                cout << i<<"," << j << endl;
            }
       }

        if (jacobians == NULL) {
            //cout << "JACOBIAN NULL" <<endl;
          return 0;
        }

         if (jacobians != NULL && jacobians[0] != NULL) {
             //cout << "filled first jac" <<endl;
             //note: null jacobian means end of computation?
            jacobians[0][0] =  Vjc(0) *d; /// @debug : put a d factor here
            jacobians[0][1] =  Vjc(1) * d;
            jacobians[0][2]=   Vjc(2) * d;

        }
         if (jacobians != NULL && jacobians[1] != NULL) {
             //note: null jacobian means end of computation?
              //cout << "filled second jac" <<endl;
            jacobians[1][0] =  init[0];
            jacobians[1][1] =  init[0];
            jacobians[1][2]=   init[0];
        }
         if (jacobians != NULL && jacobians[2] != NULL) {
             //note: null jacobian means end of computation?
             // cout << "filled third jac" <<endl;
            jacobians[2][0] =  init[0];
            jacobians[2][1] =  init[0];
            jacobians[2][2]=   init[0];
        }

//         cout << "end of evaluate" <<endl ;
        return true;
      }
 private:
    const double scalar_angle; //! vect_1.vect_2/(norm(vect_1)*norm(vect_2)) (original position)
    const double cross_angle;  //! vect_1xvect_2/(norm(vect_1)*norm(vect_2)) (original position)
};





class ManualOriginalSpacing  : public ceres::SizedCostFunction<1,3,3> {
public :
   //! this is the constructor, it expects an array of at least 3 doubles.
    ManualOriginalSpacing(Eigen::Vector3d input_vect)
        :initial_spacing_(input_vect) {}

    //! this is the operator computing the cost
    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const {
        //the parameters are as follow : parameter[0-2] = n_i = first node;parameter[3-5] = n_j = second node;

        //map the input array into 2 eigen vectors into Eigen
        //cout << "\nbeginning of evaluate" << endl;
        ConstVectorRef Ni( parameters[0],3 );
        ConstVectorRef Nj( parameters[1],3 );
        Eigen::Vector3d Is = initial_spacing_;

        //compute the cost, that is the change in distance between orignial spacing and new spacing.
        double cost = (Ni-Nj).norm() - Is.norm() ;

        //compute the sign of the jacobian : if NiNj are too close, negativ, if too far, positiv
        int sign = ((  cost > 0) - (cost < 0));

        //compute director vector of (NiNj) ie : NiNj/norm(NiNj) = u
        Eigen::Vector3d U = (Nj-Ni).normalized();

        //write residual
        residuals[0] = pow(cost,2);

        //compute Jacobian
        Eigen::Vector3d Ji = -1 * sign* U * cost/2 ;
        Eigen::Vector3d Jj = +1 * sign* U * cost/2 ;

//        cout << "initial spacing : " << Is.transpose() << endl;
//        cout << "  Ni : " << Ni.transpose()     << endl;
//        cout << " Nj : " << Nj.transpose()      << endl;
//        cout << "  Cost :" << cost     << endl;
//        cout << "   Ji :" << Ji.transpose()     <<endl;
//        cout << "   Jj :" << Jj.transpose()     << endl;


            if (jacobians == NULL) {
                //cout << "JACOBIAN NULL" <<endl;
              return 1;
            }

             if (jacobians != NULL && jacobians[0] != NULL) {
                jacobians[0][0] =  Ji(0);
                jacobians[0][1] =  Ji(1);
                jacobians[0][2]=   Ji(2);

            }
             if (jacobians != NULL && jacobians[1] != NULL) {
                 //note: null jacobian means end of computation?
                jacobians[1][0] =  Jj(0);
                jacobians[1][1] =  Jj(1);
                jacobians[1][2]=   Jj(2);

            }

//            cout<<"end of evaluate()\n";
            return true;
      }
 private:
     Eigen::Vector3d initial_spacing_; /**< store the original 3D vector between the 2 nodes*/
};





//! trying to replace the auto computation of Jacobian by a custom one
class ManualRepulsiveObject  : public ceres::SizedCostFunction<1,3,3> {

public :
   //! this is the constructor, it expects an array of at least 3 doubles.
    ManualRepulsiveObject(const double* input_vect,const double * input_w, const observation * input_obs)
        :position_(input_vect), w_i_j_(input_w), obs_(input_obs) {}

    //! this is the operator computing the cost, ie the distance projeted on normal of (n_i,n_j)
    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const {
        //the parameters is as follow : parameter[0-2] = n_i = first node;parameter[3-5] = n_j = second node;

        //map the input array into 2 eigen vectors, plus map observation position into Eigen
        //cout << "\nbeginning of evaluate" << endl;
        ConstVectorRef Ni( parameters[0],3 );
        ConstVectorRef Nj( parameters[1],3 );
        ConstVectorRef Ob(position_,3);

        //compute the normal of the plan defined by vect(NiO, NiNj)/norm(...) = n
        Eigen::Vector3d Np = (Ob-Ni).cross(Nj-Ni) ;
        //Np = Np.norm();
        //compute director vector of (NiNj) ie : NiNj/norm(NiNj) = u
        Eigen::Vector3d U = (Nj-Ni)/(Nj-Ni).norm();
        //compute residual = distance from O to NiNj : norm(vect(NiO,NiNj))/norm(NiNj)
        double d = Np.norm()/(Nj-Ni).norm()-  w_i_j_[0]/2.0;
        residuals[0]= pow(d,2) *   obs_->confidence * obs_->weight ;
        //compute Jacobian director vector : vect(u,n)
        Eigen::Vector3d Vja = U.cross(Np/Np.norm());
        //if(obs_->obs_id ==1 ) {Vja << -0.7,-0.7,0;}
        //else{Vja << +0.7,+0.7,0;}

        //compute the direction of movement : - = toward the obs, + = away from point
        int sign = ((  residuals[0]  >0) - (residuals[0] <0));
        //compute Jacobian norm for Ni : for test simply take d
        Eigen::Vector3d Ji = -1 * sign* Vja * d ;
        //compute Jacobian norm for Nj : for test simply take d
        Eigen::Vector3d Jj = -1 * sign * Vja *  d ;/// @warning : remove this 10 factor aspa !

//        cout << "  Observation_id : " <<  obs_->obs_id <<std::endl;
//         cout << "  Ni : " << Ni.transpose() <<std::endl;
//        cout << " Nj : " << Nj.transpose() <<std::endl;
//        cout << " Vja : " << Vja.transpose() <<std::endl;
//        cout << "  distance : " << residuals[0] <<std::endl;
//        cout << "   Ji :" << Ji.transpose() <<endl;
//        cout << "   Jj :" << Jj.transpose() <<endl;
        // std::cout << "\njac (eigen): \n" << jac << std::endl;

        if (jacobians == NULL) {
            //cout << "JACOBIAN NULL" <<endl;
          return 1;
        }

         if (jacobians != NULL && jacobians[0] != NULL) {
             //note: null jacobian means end of computation?
            jacobians[0][0] =  Ji(0);
            jacobians[0][1] =  Ji(1);
            jacobians[0][2]=   Ji(2);

        }
         if (jacobians != NULL && jacobians[1] != NULL) {
             //note: null jacobian means end of computation?
            jacobians[1][0] =  Jj(0);
            jacobians[1][1] =  Jj(1);
            jacobians[1][2]=   Jj(2);

        }

       //  cout<<"end of evaluate()\n";
        return true;
      }
 private:
    const double* position_; /**< store the 3D position of the observation point */
    const double* w_i_j_;//! the width of the edge, in meter
    const observation * obs_;//link to the observation, for easy use of confidence and weight.
};



#endif // CONSTRAINTS_H
