#ifndef _BCQL_H_
#define _BCQL_H_

#include <cblas.h>
#include <lapacke/lapacke.h>


//TODO Only one order is currently supported...
enum BCQL_Q_ORDER { BcqlWXYZ = 130 };
enum BCQL_V_ORDER { BcqlXYZ  = 140 };
enum BCQL_M_ORDER { BcqlColMajor = 102 }; //TODO Add Row Major


//copies first "len" entries of src to corresponding entries in dst
int dDub ( double *dst , double *src , unsigned int len ) ;
// Sets first "len" entries of array to 0;
int dZeros ( double *array , unsigned int len ) ;
// puts in dst the element-by-element sum of src1 and src2
int dSum ( double *dst , double alpha , double *src1 , double beta , double *src2 , unsigned int len ) ;

// Compute cross product between 2 3D vecs
int dCross ( double *c , double *a , double *b ) ;


// Set the square matrix mat to the identity matrix
int dEye ( double *mat , unsigned int order ) ;




// Stores in quaternion c the hamilton product between quaternions a and b
int dHamilton (const enum BCQL_Q_ORDER Order , double *c , double *a , double *b ) ;

// Store in b the inverse of quaternion a
int dQInv( const enum BCQL_Q_ORDER Order , double *b , double *a) ;

//TODO vec arma_quat_between_vecs(const vec& a,const vec& b);

/***********************************************
* Compute angular velocity given the orientation quaternion
* and its derivative.
***********************************************/
//int dQuatInv( const enum LIBQUAT_ORDER BOrder , double *v , const enum LIBQUAT_ORDER AOrder , double *a) ;
int dQD2Vel ( const enum BCQL_V_ORDER v_order    ,
              const enum BCQL_Q_ORDER q_order    ,
              double *v                          ,
              double *q                          ,
              double *dqdt                       );


/*TODO**********************************************
* Compute angular accelleration given the orientation quaternion
* and its derivatives.
***********************************************/
/*vec arma_quat_d_to_acc(const vec& q0,const vec& q1,const vec& q2){
	vec w1(4);
	w1 = 2 * arma_quat_hamilton( q2-arma_quat_hamilton(arma_quat_hamilton(q1,arma_quat_inv(q0)),q1) , arma_quat_inv(q0) );
	return(w1.subvec(1,3));
}
*/

/***********************************************
* Rotate vector v by rotation q.
* This formula should be faster than:
* v_new = q * v * q^-1
***********************************************/
//TODO vec arma_quat_rot ( const vec& q,const vec& v ) ;

/***********************************************
* Return matrix M from vector v s.t. for any x
* M*x = cross(v,x)
***********************************************/
int dCrossMat ( const enum BCQL_M_ORDER MOrder , double *M , const double *v ) ;

/***********************************************
* Applyes Rodriguez formula to quaternion q.
* Store result in matrix R
***********************************************/
int dQ2R( const enum BCQL_M_ORDER Order , const enum BCQL_Q_ORDER , double *R , double *q) ;


/***********************************************
* Conversion from rotation matrix to quaternion
***********************************************/
//TODO vec arma_mat_to_q ( const mat& R ) ;

#endif

