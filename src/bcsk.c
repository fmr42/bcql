#include "bcsk.h"


int dDub ( double *dst , double *src , unsigned int len ) {
  int i ;
  for ( i=0 ; i<len ; i++ ) {
    dst[i] = src[i] ;
  }
  // TODO Error check
  return 0 ;

}


int dSum2 ( double *dst , double alpha , double *src1 , double beta , double *src2 , unsigned int len ) {
  int i ;
  for ( i=0 ; i<len ; i++ ) {
    dst[i] = alpha * src1[i] + beta * src2[i] ;
  }
  // TODO Error check
  return 0 ;

}

// sum up 5 arrays
int dSum5 ( double *dst ,
            double a1 , double *src1 ,
            double a2 , double *src2 ,
            double a3 , double *src3 ,
            double a4 , double *src4 ,
            double a5 , double *src5 ,
            unsigned int len ) {
  int i ;
  for ( i=0 ; i<len ; i++ ) {
    dst[i] = a1 * src1[i] + a2 * src2[i] + a3 * src3[i] + a4 * src4[i] + a5 * src5[i]  ;
  }
  
  return(0) ;
}

// sum up 3 arrays
int dSum3 ( double *dst ,
            double a1 , double *src1 ,
            double a2 , double *src2 ,
            double a3 , double *src3 ,
            unsigned int len ) {
  int i ; 
  for ( i=0 ; i<len ; i++ ) {
    dst[i] = a1 * src1[i] + a2 * src2[i] + a3 * src3[i] ;
  }

  return(0) ;
}



// sum up 4 arrays
int dSum4 ( double *dst ,
            double a1 , double *src1 ,
            double a2 , double *src2 ,
            double a3 , double *src3 ,
            double a4 , double *src4 ,
            unsigned int len ) {
  int i ;
  for ( i=0 ; i<len ; i++ ) {
    dst[i] = a1 * src1[i] + a2 * src2[i] + a3 * src3[i] + a4 * src4[i] ;
  }

  return(0) ;
}



// Transpose square matrix
int dSqTr ( int n , double *src , double *dst ) {
  // if src == dst => have to use a tmp array
  if ( src == dst ) {
    double *tmp;
    tmp = malloc ( sizeof(double) * n * n ) ;
    dDub ( tmp , src , n*n ) ;
    src = tmp ;
  }

  // transpose
  int i,j;
  for ( i=0 ; i<n ; i++ ) {
    for ( j=0 ; j<n ; j++ ) {
      dst[i*n+j] = src[j*n+i] ;
    }
  }
  
  return (0) ;
}




int dZeros ( double *array , unsigned int len ) {
  int i ;
  for ( i=0 ; i<len ; i++ ) {
    array[i] = 0 ;
  }
  // TODO check for errors

  return 0;
}



int dCross ( double *c , double *a , double *b ) {

  c[0] = a[1] * b[2] - a[2] * b[1] ;
  c[1] = a[2] * b[0] - a[0] * b[2] ;
  c[2] = a[0] * b[1] - a[1] * b[0] ;

  return (0);

} ;


// Store in c the hamilton product between a and b
int dHamilton (const enum BCQL_Q_ORDER Order , double *c , double *a , double *b ) {

  if ( Order == BcqlWXYZ ) {
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
int dQInv( const enum BCQL_Q_ORDER Order ,  double *b , double *a) {
  if ( Order == BcqlWXYZ ) {
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
int dQD2Vel ( const enum BCQL_V_ORDER v_order    ,
              const enum BCQL_Q_ORDER q_order    ,
              double *v                          ,
              double *q                          ,
              double *dqdt                       ) {

  if ( q_order == BcqlWXYZ && v_order == BcqlXYZ ) {
    double adj_res[4];
    double tmp[4] ;

    dQInv ( BcqlWXYZ , tmp , q ) ;
    cblas_dscal ( 4 , 2 , tmp , 1 );

    dHamilton ( BcqlWXYZ , adj_res , tmp , dqdt ) ;
    v[0] = adj_res[1] ;
    v[1] = adj_res[2] ;
    v[2] = adj_res[3] ;
  } else {
    return(1);
  }
  
  return (0);  

}



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
int dCrossMat ( const enum BCQL_M_ORDER MOrder , double *M , const double *v ) {
//TODO check MOrder
  
  M[0] = 0     ;      
  M[1] =  v[2] ;
  M[2] = -v[1] ;
  M[3] = -v[2] ;
  M[4] = 0     ;
  M[5] =  v[0] ;
  M[6] =  v[1] ;
  M[7] = -v[0] ;
  M[8] = 0     ;
  return ( 0 );
}



// Set the square matrix mat to the identity matrix
int dEye ( unsigned int order , double scale , double *mat ) {
  dZeros(mat,order*order);
  int i;
  for ( i=0 ; i < (order*order) ; i=i+order+1 ) {
    mat[i]=scale;
  }
  return (0);
}



// Computes rotation matrix R from quaternion q
// R(q) = I + 2 η S(ǫ) + 2 (S(ǫ)^2)
// For details see:
// M.D. Shuster. A survey of attitude representation. The Journal of the
// Astronautical Sciences, 41(4):439–517, December 1993.
int dQ2R( const enum BCQL_M_ORDER QOrder , const enum BCQL_Q_ORDER MOrder , double *R , double *q) {

  double eye [9] ;
  dEye(3,1,eye);

  double cross[9];
  dCrossMat ( BcqlColMajor , cross , q+1 ) ;

  double tmp1[9];
  dDub(tmp1,cross,9);

  dSum2 ( R , 1 , eye , 2*q[0] , tmp1 , 9 );

  cblas_dgemm ( CblasColMajor , CblasNoTrans , CblasNoTrans , 3 , 3 , 3 , 2 , cross , 3 , cross , 3 , 0 , tmp1 , 3); 

  dSum2 ( R , 1 , R , 1 , tmp1 , 9 );
  
  return (0);
}


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


