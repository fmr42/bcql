#include "libquat.h"


// Store in c the hamilton product between a and b
int dHamiltonProd (const enum LIBQUAT_ORDER COrder , double *c ,
                  const enum LIBQUAT_ORDER AOrder , double *a , 
                  const enum LIBQUAT_ORDER BOrder , double *b ) {

  if ( AOrder == LibQuatWXYZ && BOrder == LibQuatWXYZ ) {
    c[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3] ;
    c[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2] ;
    c[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1] ;
    c[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0] ;
  } else {
    // Not implemented
    return (1);
  }

  return(0);
}

// Store in b the inverse of quaternion a
int dQuatInv( const enum LIBQUAT_ORDER BOrder , double *b , const enum LIBQUAT_ORDER AOrder , double *a) {
  if ( AOrder == LibQuatWXYZ && BOrder == LibQuatWXYZ ) {
    double tmp_mod = cblas_dnrm2 ( 4, a , 1);
    b[0] =   a[0] / tmp_mod	;
    b[1] = - a[1] / tmp_mod	;
    b[2] = - a[2] / tmp_mod	;
    b[3] = - a[3] / tmp_mod	;
  } else {
    return (1) ;
  }
  return(0);
}




/*TODO
// Compute and store in q the rotatation quaternion from a to b
// TODO rotation of 180 degrees is a special case and needs to be implemented
vec arma_quat_between_vecs(const vec& a,const vec& b){
	vec q(4);
	// Compute axis/angle representation 
	double cos_teta = dot( (a),(b) ) / ( norm((a),2) * norm((b),2) )	;
	vec axis = cross((a),(b))			;
	axis = normalise(axis)				;
	// Compute quaternion 
	double half_teta = acos(cos_teta) / 2	;
	q(0) = cos(half_teta)		;
	q(1) = axis(0)	* sin(half_teta);
	q(2) = axis(1)	* sin(half_teta);
	q(3) = axis(2)	* sin(half_teta);
	q=normalise(q);
	return(q);
}
*/




/***********************************************
* Compute angular velocity given the orientation quaternion
* and its derivative.
***********************************************/
/*TODO
vec arma_quat_d_to_vel(const vec& q0,const vec& q1){
	vec w(4);
	w = arma_quat_hamilton( 2 * arma_quat_inv(q0) , q1 );
	return(w.subvec(1,3));
}
*/


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
/*TODO
vec arma_quat_rot ( const vec& q,const vec& v ) {
	return ( v + 2*cross( q.subvec(1,3) , cross ( q.subvec(1,3),v)+q(0)*v) );
}
*/



/***********************************************
* Return matrix M from vector v s.t. for any x
* M*x = cross(v,x)
***********************************************/
/*
mat arma_cross_to_mat ( const vec& v ) {
	mat M = zeros<mat>(3,3);
	M(0,1) = -v(2);
	M(0,2) =  v(1);
	M(1,0) =  v(2);
	M(1,2) = -v(0);
	M(2,0) = -v(1);
	M(2,1) =  v(0);
	return ( M );
}
*/

/***********************************************
* Applyes Rodriguez formula to quaternion q.
* Returns rotation matrix
***********************************************/
/*
mat arma_rodriguez ( const vec& q ) {
	mat tmp = arma_cross_to_mat( q.subvec(1,3) );
	return ( eye<mat>(3,3) + 2*q(0)*tmp +2*tmp*tmp );
}
*/






/***********************************************
* Conversion from rotation matrix to quaternion
***********************************************/
/*
vec arma_mat_to_q ( const mat& R ) {
	vec q = zeros<vec>(4) ;
	q(0) = 0.5 * sqrt ( 1 + R(0,0) + R(1,1) + R(2,2) ) ;

	double t =  R(0,0) + R(1,1) + R(2,2) ;
	double r = sqrt(1+t) ;
	q(0) = 0.5*r ;
	q(1) = copysign(0.5*sqrt(1+R(0,0)-R(1,1)-R(2,2)), R(2,1)-R(1,2)) ;
	q(2) = copysign(0.5*sqrt(1-R(0,0)+R(1,1)-R(2,2)), R(0,2)-R(2,0)) ;
	q(3) = copysign(0.5*sqrt(1-R(0,0)-R(1,1)+R(2,2)), R(1,0)-R(0,1)) ;
	if ( isnan(q(0)) )
		q(0)=0;
	if ( isnan(q(1)) )
		q(1)=0;
	if ( isnan(q(2)) )
		q(2)=0;
	if ( isnan(q(3)) )
		q(3)=0;

	return ( q ) ;
}
*/


