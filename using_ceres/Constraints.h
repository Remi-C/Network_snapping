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

#include "Parameters.h"
#include "utils_function.h"

#include "ceres/ceres.h"
#include "ceres/rotation.h"

#include <cmath>


typedef Eigen::Map<Eigen::Vector3d> VectorRef;
typedef Eigen::Map<const Eigen::Vector3d> ConstVectorRef;


using std::cout;
using std::endl;

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::LossFunction;


int addAllConstraints(DataStorage * data, ceres::Problem * problem, Parameter* param);
int addManualConstraintsOnInitialPosition(DataStorage *, ceres::Problem * );

int addManualConstraintsOnDistanceToOriginalAngle(DataStorage * , Problem * ) ;
int addManualConstraintsOnInitialspacing(DataStorage * , Problem * ) ;

int addManualConstraintsOnOrthDistToObservation(DataStorage * , Problem * ) ;
int addManualConstraintsOnSurfDistToObjects(DataStorage * , Problem * ) ;
int addManualConstraintsOnOrthDistToObservation_width(DataStorage * , Problem * ) ;
int addManualConstraintsOnSurfDistToObjects_width(DataStorage * , Problem * ) ;



/** functor to compute cost between one node position and this node original position
  */
class ManualDistanceToInitialPosition  : public ceres::SizedCostFunction<1,3> {
public :
    //! this is the constructor, it expects an array of at least 3 doubles.
    ManualDistanceToInitialPosition(Eigen::Vector3d input_vect)
        :initial_position_(input_vect) {}

    //! this is the operator computing the cost
    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const {
        //the parameters are as follow : parameter[0-2] = n_i = first node
        //we also get the global parameter of this programm
        int is_width = 0;
        Parameter * param = Parameter::instance() ;
        if(param->optimisation_type==SnapEnums::WIDTH){
            residuals[0] = 0 ;
            if ((jacobians != NULL) && (jacobians[0]!=NULL)){
                jacobians[0][0] =  0 ;
                jacobians[0][1] =  0 ;
                jacobians[0][2]=   0 ;

            }
            is_width = 1;
            return 1;
        }


        //map the input array into 2 eigen vectors into Eigen
        ConstVectorRef Ni( parameters[0],3 ); //the node position
        Eigen::Vector3d Is = initial_position_;

        //compute the cost, that is the eucl dist to original position
        double cost =  (Ni - Is).norm()  ;

        //compute director vector of cost : Ni-Is. If null, can't normalize!
        Eigen::Vector3d U ;
        if(std::abs(cost)<TOLERANCE ){
            U = Eigen::Vector3d::Constant(1) ;
        }else{
            U = (Ni - Is).normalized();
        }
        //write residual
        residuals[0] = cost;

        //compute Jacobian
        Eigen::Vector3d Ji = - U * cost * (1-is_width);


        if (jacobians == NULL) {
            //cout << "JACOBIAN NULL" <<endl;
            return 1;
        }

        if (jacobians != NULL && jacobians[0] != NULL) {
            jacobians[0][0] =  Ji(0);
            jacobians[0][1] =  Ji(1);
            jacobians[0][2]=   0 ;//Ji(2);

        }
        //            cout<<"end of evaluate()\n";
        return true;
    }
private:
    Eigen::Vector3d initial_position_; /**< store the original 3D position of the node*/
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
        int is_width = 0;
        Parameter * param = Parameter::instance() ;
        if(param->optimisation_type==SnapEnums::WIDTH){
            residuals[0] = 0 ;
            residuals[1] = 0 ;
            if(jacobians != NULL){
                for(int i = 0; i<3; ++i){
                    if(jacobians[i] != NULL) {
                        for(int j=0;j<6;++j){
                            jacobians[i][j] = 0 ;
                        }
                    }
                }
            }
            is_width = 1 ;
            return 1;
        }

        //cout << "begining of evaluate" <<endl ;
        //map the input array into 3 eigen vectors
        ConstVectorRef Nc( parameters[0],3 );
        ConstVectorRef Ni( parameters[1],3 );
        ConstVectorRef Nj( parameters[2],3 );

        double scalar_2a = (Nc-Ni).dot(Nc-Nj)/((Nc-Ni).norm() * (Nc-Nj).norm());
        double cross_2a = ((Nc-Ni).cross(Nc-Nj)/((Nc-Ni).norm() * (Nc-Nj).norm())).norm();

        //compute jacobian :
        //only the center node (Nc) should be moved
        //the direction is along the bissect of angle (Ni,Nc,Nj).
        //in the case where NcNi and NcNj are colinear, we need an alternate method
        Eigen::Vector3d Vjc;
        if(std::abs(scalar_2a)>0.5){//regular case
            Vjc = ((Nc-Ni).normalized() + (Nc-Nj).normalized() ).normalized();
        }else{//computing the bisect of NcNi, NjNc , then 90° rotation regarding Z axis
            Vjc = ((Nc-Ni).normalized() + (Nj-Nc).normalized() ).normalized();
            Vjc = Vjc.cross(Eigen::Vector3d::UnitZ());
        }
        for (int i =0; i<3;++i){
            if(std::isnan(Vjc[i])==true){
                Vjc << 0,0,0 ;
            }
        }

        /** value of displacement :
                  lets consider a triangle formed on angle alpha (Ni,Nc,Nj) (upper summit is Nc)
                  but with Ni,Nc and Nc,Nj normalized to 1.
                  we want to find Nc_p, that is the new position of Nc after a move of d along the bissect.
                  the final expected angle (Ni,Nc',Nj) is the original angle beta (caracterized by sin and cos)
                  the distance d = sin(alpha)*tan(beta)-cos(alpha)
                  Note that sign determin which direction it goes
                  Note : in theory every angle shoud be divided by 2,thus we may have ot use sin2x=f(sinx^2, sinx) to have correct result
                  */
        //double d = cross_a * cross_angle/scalar_angle  -scalar_a ;
        //double d = cos_a - cos_o / sin_a ;
        //double d = SIGN(scalar_angle-scalar_2a + cross_angle -cross_2a ) *( std::abs(scalar_angle-scalar_2a) + std::abs(cross_angle -cross_2a) ) / 2.0 * std::max((Nc-Ni).norm()/(Nc-Nj).norm(),(Nc-Nj).norm()/(Nc-Ni).norm()) ;
        //compute residuals (= cost)
        //std::cout << "scalar_2a:  " << scalar_2a << " scalar_angle : " << scalar_angle
        //             << "cross_2a:  " << cross_2a << " cross_angle : " << cross_angle << std::endl;
        residuals[0] = scalar_2a - scalar_angle  ;// scalar_angle-scalar_a  ;
        residuals[1] = cross_2a - cross_angle ;// cross_angle-cross_a  ;

        double d  = std::abs(residuals[0]) + std::abs(residuals[1]) *(1-is_width);

        Eigen::Vector3d Vj1 = - Vjc * d *10;
        Eigen::Vector3d Vj2 = - Vjc * d *10 ;


        //Eigen::Vector3d Vj1 = Vjc * (residuals[0]+residuals[1])*10;
        //Eigen::Vector3d Vj2 =  Vjc *(residuals[0]+residuals[1])*10 ;

        if (jacobians == NULL) {
            //    cout << "JACOBIAN NULL" <<endl;
            return 1;
        }

        if (jacobians != NULL && jacobians[0] != NULL) {
            //    cout << "filled first jac" <<endl;
            //note: null jacobian means end of computation?
            jacobians[0][0] = Vj1(0) ;
            jacobians[0][1] = Vj1(1) ;
            jacobians[0][2] = 0 ; //Vj1(2) ;
            jacobians[0][3] = Vj2(0) ;
            jacobians[0][4] = Vj2(1) ;
            jacobians[0][5] = 0 ; //Vj2(2) ;

        }
        if (jacobians != NULL && jacobians[1] != NULL) {
            //note: null jacobian means end of computation?
            //   cout << "filled second jac" <<endl;
            jacobians[1][0] =  0;
            jacobians[1][1] =  0;
            jacobians[1][2]=   0;
            jacobians[1][3] =  0;
            jacobians[1][4]=   0;
            jacobians[1][5]=   0;
        }
        if (jacobians != NULL && jacobians[2] != NULL) {
            //note: null jacobian means end of computation?
            //    cout << "filled third jac" <<endl;
            jacobians[2][0] =  0;
            jacobians[2][1] =  0;
            jacobians[2][2]=   0;
            jacobians[2][3] =  0;
            jacobians[2][4] =  0;
            jacobians[2][5]=   0;
        }

        //cout << "end of evaluate" <<endl ;
        return 1;
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
        //this is a way to bypass all computation in the case where we don't optimise on position
        int is_width = 0;
        Parameter * param = Parameter::instance() ;
        if(param->optimisation_type==SnapEnums::WIDTH){
            residuals[0] = 0 ;
            if(jacobians != NULL){
                for(int i = 0; i<2; ++i){
                    if(jacobians[i] != NULL) {
                        for(int j=0;j<3;++j){
                            jacobians[i][j] = 0 ;
                        }
                    }
                }
            }
            is_width = 1 ;
            return 1;
        }

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
        Eigen::Vector3d Ji = -1 * sign* U * cost /2.0 *(1-is_width);// /2
        Eigen::Vector3d Jj = +1 * sign* U * cost /2.0 *(1-is_width);// /2


        if (jacobians == NULL) {
            //cout << "JACOBIAN NULL" <<endl;
            return 1;
        }

        if (jacobians != NULL && jacobians[0] != NULL) {
            jacobians[0][0] =  Ji(0);
            jacobians[0][1] =  Ji(1);
            jacobians[0][2]=   0 ;// Ji(2);

        }
        if (jacobians != NULL && jacobians[1] != NULL) {
            //note: null jacobian means end of computation?
            jacobians[1][0] =  Jj(0);
            jacobians[1][1] =  Jj(1);
            jacobians[1][2]= 0 ;//   Jj(2);

        }

        //            cout<<"end of evaluate()\n";
        return true;
    }
private:
    Eigen::Vector3d initial_spacing_; /**< store the original 3D vector between the 2 nodes*/
};






//! trying to replace the auto computation of Jacobian by a custom one
class ManualOrthDistanceToObservation  : public ceres::SizedCostFunction<1,3,3,1> {

public :
    //! this is the constructor, it expects an array of at least 3 doubles.
    ManualOrthDistanceToObservation(const double* input_vect,const double * input_w, const observation * input_obs)
        :position_(input_vect), w_i_j_(input_w), obs_(input_obs) {}

    //! this is the operator computing the cost, ie the distance projeted on normal of (n_i,n_j)
    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const {
        //the parameters is as follow : parameter[0-2] = n_i = first node;parameter[3-5] = n_j = second node;
        //parameter[6] : width of the edge

        //map the input array into 2 eigen vectors, plus map observation position into Eigen
        //cout << "\nbeginning of evaluate" << endl;
        ConstVectorRef Ni( parameters[0],3 );
        ConstVectorRef Nj( parameters[1],3 );
        ConstVectorRef Ob(position_,3);

        int is_width = 0;
        Parameter * param = Parameter::instance() ;
        if(param->optimisation_type==SnapEnums::WIDTH){
            is_width = 1;
        }
        //compute the normal of the plan defined by vect(NiO, NiNj)/norm(...) = n
        Eigen::Vector3d Np = (Ob-Ni).cross(Nj-Ni) ;
        //Np = Np.norm();
        //compute director vector of (NiNj) ie : NiNj/norm(NiNj) = u
        Eigen::Vector3d U = (Nj-Ni)/(Nj-Ni).norm();
        //compute residual = distance from O to NiNj : norm(vect(NiO,NiNj))/norm(NiNj)
        double d = (Np.norm()/(Nj-Ni).norm()-  parameters[2][0]/2.0) ;


        residuals[0]=  ceres::pow(d ,2)* obs_->confidence * obs_->weight;
        //compute Jacobian director vector : vect(u,n)
        Eigen::Vector3d Vja = U.cross(Np/Np.norm());
        //if(obs_->obs_id ==1 ) {Vja << -0.7,-0.7,0;}
        //else{Vja << +0.7,+0.7,0;}

        double t = std::sqrt( std::pow((Ob-Ni).norm(),2) - std::pow(Np.norm()/(Nj-Ni).norm(),2) )/(Nj-Ni).norm() ;
        t  = (10+t)/11.0 ;
        double coeff_i = 1 ;
        double coeff_j = 1 ;
        //        if(t > 0.5 ){
        //            //obs is closer to Nj
        //            coeff_i = 0.5+t ;
        //            coeff_j = 1-t ;
        //        } else{
        //            //obs is closer to Ni
        //            coeff_i = 1-t ;
        //            coeff_j = 0.5 + t ;
        //        }
        //compute the direction of movement : - = toward the obs, + = away from point
        //compute Jacobian norm for Ni : for test simply take d
        Eigen::Vector3d Ji = - 1 *  Vja *  d * (coeff_i) * (1-is_width);
        //compute Jacobian norm for Nj : for test simply take d
        Eigen::Vector3d Jj = - 1 *  Vja *  d *(coeff_j) * (1-is_width);


        //        cout << "  Observation_id : " <<  obs_->obs_id <<std::endl;
        //         cout << "  Ni : " << Ni.transpose() <<std::endl;
        //        cout << " Nj : " << Nj.transpose() <<std::endl;
        //        cout << " Vja : " << Vja.transpose() <<std::endl;
        //        cout << "  distance : " << residuals[0] <<std::endl;
        //        cout << "   Ji :" << Ji.transpose() <<endl;
        //        cout << "   Jj :" << Jj.transpose() <<endl;
        // std::cout << "\njac (eigen): \n" << jac << std::endl;

        if (jacobians == NULL) {
            return 1;
        }

        if (jacobians != NULL && jacobians[0] != NULL) {
            jacobians[0][0] =  Ji(0);
            jacobians[0][1] =  Ji(1);
            jacobians[0][2]= 0 ;
        }
        if (jacobians != NULL && jacobians[1] != NULL) {
            jacobians[1][0] = Jj(0);
            jacobians[1][1] = Jj(1);
            jacobians[1][2]= 0 ;
        }
        if (jacobians != NULL && jacobians[2] != NULL) {
            jacobians[2][0] = -1 * SIGN(d) * std::abs(d) * is_width;
        }
        //  cout<<"end of evaluate()\n";
        return true;
    }
private:
    const double* position_; /**< store the 3D position of the observation point */
    const double* w_i_j_;//! the width of the edge, in meter. point to the edge->width memory
    const observation * obs_;//link to the observation, for easy use of confidence and weight.
};




class ManualOrthDistanceToObservation_width  : public ceres::SizedCostFunction<1,3,3,1> {

public :
    //! this is the constructor, it expects an array of at least 3 doubles.
    ManualOrthDistanceToObservation_width(const double* input_vect,const double * input_w, const observation * input_obs)
        :position_(input_vect), w_i_j_(input_w), obs_(input_obs) {}

    //! the cost is proportionnal to the shortest distance between observation and segment minus segment width
    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const {
        //the parameters is as follow : parameter[0-2] = n_i = first node;parameter[3-5] = n_j = second node;
        //parameter[6] : width of the edge

        //map the input array into 2 eigen vectors, plus map observation position into Eigen
        //cout << "\nbeginning of evaluate" << endl;
        ConstVectorRef Ni( parameters[0],3 );
        ConstVectorRef Nj( parameters[1],3 );
        ConstVectorRef Ob(position_,3);

        //compute the normal of the plan defined by vect(NiO, NiNj)/norm(...) = n
        Eigen::Vector3d Np = (Ob-Ni).cross(Nj-Ni) ;
        //compute director vector of (NiNj) ie : NiNj/norm(NiNj) = u
        Eigen::Vector3d U = (Nj-Ni)/(Nj-Ni).norm();
        //compute residual = distance from O to NiNj : norm(vect(NiO,NiNj))/norm(NiNj)

        double d = (Np.norm()/(Nj-Ni).norm()-  parameters[2][0]/2.0) ;
        residuals[0]= std::pow(d,2)* obs_->confidence * obs_->weight ;


        //compute Jacobian director vector : vect(u,n)
        Eigen::Vector3d Vja = U.cross(Np/Np.norm());
        //if(obs_->obs_id ==1 ) {Vja << -0.7,-0.7,0;}
        //else{Vja << +0.7,+0.7,0;}

        //compute Jacobian norm for Ni : for test simply take d
        Eigen::Vector3d Ji = -1 * Vja * ceres::sqrt(residuals[0])*SIGN(d); ;
        //compute Jacobian norm for Nj : for test simply take d
        Eigen::Vector3d Jj = Ji ;

        double t = std::sqrt( std::pow((Ob-Ni).norm(),2) - std::pow(Np.norm()/(Nj-Ni).norm(),2) )/(Nj-Ni).norm() ;


        if (jacobians == NULL) {
            //cout << "JACOBIAN NULL" <<endl;
            return 1;
        }

        if (jacobians != NULL && jacobians[0] != NULL) {
            //note: null jacobian means end of computation?
            jacobians[0][0] = Ji(0)/2 * t * 0;
            jacobians[0][1] = -Ji(1)/2 * t * 0 ;
            jacobians[0][2] = 0;// -Ji(2);
        }
        if (jacobians != NULL && jacobians[1] != NULL) {
            //note: null jacobian means end of computation?
            jacobians[1][0] =  Jj(0)/2 * (1-t) * 0 ;
            jacobians[1][1] =  Jj(1)/2 * (1-t) * 0;
            jacobians[1][2] = 0;// -Jj(2);
        }
        if (jacobians != NULL && jacobians[2] != NULL) {
            //note: null jacobian means end of computation?
            jacobians[2][0] = -1 * SIGN(d) * std::abs(d) ;//0; //-1 *  residuals[0]/2  ;

        }

        //  cout<<"end of evaluate()\n";
        return true;
    }
private:
    const double* position_; /**< store the 3D position of the observation point */
    const double* w_i_j_;//! the width of the edge, in meter. point to the edge->width memory
    const observation * obs_;//link to the observation, for easy use of confidence and weight.
};




//! trying to replace the auto computation of Jacobian by a custom one
class ManualAttr_Rep_Object_width  : public ceres::SizedCostFunction<1,3,3,1> {

public :
    //! this is the constructor, it expects the index of the street_object

    ManualAttr_Rep_Object_width(const unsigned int index, DataStorage * data)
        : object_index_(index)
        , obj_(data->street_objects(index))
        ,classification_(data->cbn(obj_->class_name))
        ,centroid2D_{0,0}
    ,axis_width_(&data->ebe(obj_->edge_id)->width)
    {
        Geometry::geomPoint2Double( obj_->geom_centroid
                                    , this->centroid2D_  );
    }

    /** this is the operator computing the cost,
        Here the cost depends on the type of object, and is proportionnal to the shared surface (+ maybe the distance to border)
        Refer to geometry_function for better description of this cost.
    */
    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const {
        //the parameters is as follow : parameter[0-2] = n_i = first node;parameter[3-5] = n_j = second node;
        //parameter[6] = edge width

        //map the input array into 2 eigen vectors, plus map observation position into Eigen
        //        cout << "\nbeginning of evaluate" << endl;




        ConstVectorRef Ni( parameters[0],3 );
        ConstVectorRef Nj( parameters[1],3 );

        //cout << centroid2D_[0] << " " << centroid2D_[1] << endl;
        //center of object : ConstVectorRef Ob(position_,3);
        Eigen::Vector3d Obj_center;
        Obj_center[0] = centroid2D_[0] ;
        Obj_center[1] = centroid2D_[1] ;
        Obj_center[2] = 0 ;


        //the jacobian should be in the plan normal to Z vector
        // we can't correct Z because of the variety of possible objects. Some may not be at road height
        Eigen::Vector3d Np = Eigen::Vector3d::UnitZ() ;

        //compute director vector of (NiNj) ie : NiNj/norm(NiNj) = u
        Eigen::Vector3d U = (Nj-Ni)/(Nj-Ni).norm();

        //the direction of change should be :
        Eigen::Vector3d Vja = (U.cross(Np)).normalized();

        //compute the cost using the geometric distance
        double cost = shared_area_cost(classification_->road_surface_relation
                                       ,parameters[0]
                                       ,parameters[1]
                                       ,parameters[2][0]
                                       ,obj_->geom_border_surface
                                       ,obj_->geom_border_area
                                       ,obj_->geom_centroid);
        double d =   1.0 * cost /obj_->geom_border_area ;
        residuals[0] =pow(d,2)  ;
        //residuals[0] = pow(std::abs(cost/obj_->geom_border_area),2) ; /// @FIXME @TODO @DEBUG warning : should put the confidence here

        int sign =-1* Geometry::orientationIndex(parameters[0],parameters[1],centroid2D_);//depends on left or right !

        //compute Jacobian norm for Ni
        Eigen::Vector3d Ji =  sign * Vja * SIGN(cost)*  ceres::sqrt(residuals[0]) ;
        //Eigen::Vector3d Ji =  sign * Vja * cost / obj_->geom_border_area;
        //compute Jacobian norm for Nj
        Eigen::Vector3d Jj =  Ji ;
        if (jacobians == NULL) {
            //cout << "JACOBIAN NULL" <<endl;
            return 1;
        }

        if (jacobians != NULL && jacobians[0] != NULL) {
            //note: null jacobian means end of computation?
            jacobians[0][0] = 0 ;// Ji(0);
            jacobians[0][1] = 0 ; //Ji(1);
            jacobians[0][2]=  0 ;// Ji(2);

        }
        if (jacobians != NULL && jacobians[1] != NULL) {
            //note: null jacobian means end of computation?
            jacobians[1][0] = 0 ;//Jj(0);
            jacobians[1][1] = 0 ; // Jj(1);
            jacobians[1][2]=  0 ;// Jj(2);

        }
        if (jacobians != NULL && jacobians[2] != NULL) {
            //note: null jacobian means end of computation?
            jacobians[2][0] = d ; //-1 * SIGN(cost)*  ceres::sqrt(residuals[0])/10 ;

        }
        //        cout << "end of evaluate()\n";
        return true;
    }

private:
    const int object_index_ ; // store the index of this object into the object array
    const street_object * obj_;//link to the observation, for easy use of confidence and weight.
    const classification * classification_; // link to the class of the object
    double centroid2D_[2]; // X and Y of the centroid of the geom.
    const double* axis_width_; // points to the edge-> width memory
    ///  @TODO : should be const
};





//! trying to replace the auto computation of Jacobian by a custom one
class ManualAttr_Rep_Object  : public ceres::SizedCostFunction<1,3,3,1> {

public :
    //! this is the constructor, it expects the index of the street_object

    ManualAttr_Rep_Object(const unsigned int index, DataStorage * data)
        : object_index_(index)
        , obj_(data->street_objects(index))
        ,classification_(data->cbn(obj_->class_name))
        ,centroid2D_{0,0}
    ,axis_width_(&data->ebe(obj_->edge_id)->width)
    {
        Geometry::geomPoint2Double( obj_->geom_centroid
                                    , this->centroid2D_  );
    }

    /** this is the operator computing the cost,
        Here the cost depends on the type of object, and is proportionnal to the shared surface (+ maybe the distance to border)
        Refer to geometry_function for better description of this cost.
    */
    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const {
        //the parameters is as follow : parameter[0-2] = n_i = first node;parameter[3-5] = n_j = second node;
        int is_width = 0;
        Parameter * param = Parameter::instance() ;
        if(param->optimisation_type==SnapEnums::WIDTH){
            is_width = 1;
        }
        //map the input array into 2 eigen vectors, plus map observation position into Eigen
        //        cout << "\nbeginning of evaluate" << endl;
        ConstVectorRef Ni( parameters[0],3 );
        ConstVectorRef Nj( parameters[1],3 );

        //cout << centroid2D_[0] << " " << centroid2D_[1] << endl;
        //center of object : ConstVectorRef Ob(position_,3);
        Eigen::Vector3d Obj_center;
        Obj_center[0] = centroid2D_[0] ;
        Obj_center[1] = centroid2D_[1] ;
        Obj_center[2] = 0 ;


        //the jacobian should be in the plan normal to Z vector
        // we can't correct Z because of the variety of possible objects. Some may not be at road height
        Eigen::Vector3d Np = Eigen::Vector3d::UnitZ() ;

        //compute director vector of (NiNj) ie : NiNj/norm(NiNj) = u
        Eigen::Vector3d U = (Nj-Ni)/(Nj-Ni).norm();

        //the direction of change should be :
        Eigen::Vector3d Vja = (U.cross(Np)).normalized();

        //compute the cost using the geometric distance
        double cost = shared_area_cost(classification_->road_surface_relation
                                       ,parameters[0]
                                       ,parameters[1]
                                       ,parameters[2][0]
                                       ,obj_->geom_border_surface
                                       ,obj_->geom_border_area
                                       ,obj_->geom_centroid);
        double d =   1.0 * cost /obj_->geom_border_area ;
        residuals[0] =pow(d,2)  ; /// @FIXME @TODO @DEBUG warning : should put the confidence here

        int sign =-1* Geometry::orientationIndex(parameters[0],parameters[1],centroid2D_);//depends on left or right !
        //compute Jacobian norm for Ni
        Eigen::Vector3d Ji =  sign * Vja * d * (1-is_width);
        //Eigen::Vector3d Ji =  sign * Vja * cost / obj_->geom_border_area;
        //compute Jacobian norm for Nj
        Eigen::Vector3d Jj =  Ji ;
        if (jacobians == NULL) {
            //cout << "JACOBIAN NULL" <<endl;
            return 1;
        }

        if (jacobians != NULL && jacobians[0] != NULL) {
            //note: null jacobian means end of computation?
            jacobians[0][0] =  Ji(0);
            jacobians[0][1] =  Ji(1);
            jacobians[0][2]=  0 ;// Ji(2);

        }
        if (jacobians != NULL && jacobians[1] != NULL) {
            //note: null jacobian means end of computation?
            jacobians[1][0] =  Jj(0);
            jacobians[1][1] =  Jj(1);
            jacobians[1][2]=  0 ;// Jj(2);

        }
        if (jacobians != NULL && jacobians[2] != NULL) {
            //note: null jacobian means end of computation?
            jacobians[2][0] = d * is_width; //-1 * SIGN(cost)*  ceres::sqrt(residuals[0])/10 ;

        }


        //output :
        //        if (std::isnan(Ji(2)) == true ) {
        //            cout << "  Ni : " << Ni.transpose() <<std::endl;
        //            cout << " Nj : " << Nj.transpose() <<std::endl;
        //            cout << " Vja : " << Vja.transpose() <<std::endl;
        //            cout << "  residual : " << residuals[0] <<std::endl;
        //            cout << "  cost : " << cost <<std::endl;
        //            cout << "   Ji :" << Ji.transpose() <<endl;
        //            cout << "   Jj :" << Jj.transpose() <<endl;
        //        }

        //        cout << "end of evaluate()\n";
        return true;
    }

private:
    const int object_index_ ; // store the index of this object into the object array
    const street_object * obj_;//link to the observation, for easy use of confidence and weight.
    const classification * classification_; // li-2.0034nk to the class of the object
    double centroid2D_[2]; // X and Y of the centroid of the geom.
    const double* axis_width_; // points to the edge-> width memory
    ///  @TODO : should be const
};




#endif // CONSTRAINTS_H
