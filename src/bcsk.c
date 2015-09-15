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
int dHamilton (const enum BCSK_Q_ORDER Order , double *c , double *a , double *b ) {

    c[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3] ;
    c[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2] ;
    c[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1] ;
    c[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0] ;
    
  return(0);
}

// Store in b the inverse of quaternion a
int dQInv( const enum BCSK_Q_ORDER Order ,  double *b , double *a) {
    double tmp_mod = a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3] ;
    b[0] =   a[0] / tmp_mod	;
    b[1] = - a[1] / tmp_mod	;
    b[2] = - a[2] / tmp_mod	;
    b[3] = - a[3] / tmp_mod	;
  return(0);
}





// Compute and store in q the rotatation quaternion from a to b
// NB for rotation of 180 degrees z will be assumed as axes of rotation

int dQbetVs ( double *Q , double *v1 , double *v2){
  double v1_n[3] ;
  double v2_n[3] ;
  double cos_teta;

  // compute cosine of the angle between the 2 vectors
  dDub ( v1_n , v1 , 3 ) ;
  dNormalise ( 3 , v1_n ) ;
  dDub ( v2_n , v2 , 3 ) ;
  dNormalise ( 3 , v2_n ) ;

  cos_teta = cblas_ddot ( 3 , v1_n , 1 , v2_n , 1 );

  // find rotation axis
  double epsilon = 0.00001 ; // TODO!!!! Hard-coded
  if ( 1 - cos_teta < epsilon ) { // no rotation
    Q[0] = 1 ;
    Q[1] = 0 ;
    Q[2] = 0 ;
    Q[3] = 0 ;
  } else if ( 1 + cos_teta < epsilon ) { // 180 degrees rotation
    Q[0] = 0 ;
    Q[1] = 0 ;
    Q[2] = 0 ;
    Q[3] = 1 ;
  } else { // not a special case
    double teta = acos ( cos_teta ) ;
    Q[0] = cos ( teta/2 );
    dCross ( Q+1 , v1_n , v2_n ) ;
    cblas_dscal ( 3 , sin ( teta/2 ) / sin(teta) , Q+1 , 1 );
  }
  
  return (0);
}


/***********************************************
* Compute angular velocity given the orientation quaternion
* and its derivative.
* For details see:
* https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions
***********************************************/
int dQD2Vel ( const enum BCSK_V_ORDER v_order    ,
              const enum BCSK_Q_ORDER q_order    ,
              double *v                          ,
              double *q                          ,
              double *dqdt                       ) {

// Formula:  [ 0 w ] = hamilton ( 2*(dq/dt) , q^-1)
//                               ____tmp1____
//                     ________tmp2_______________
    double tmp1[4];
    double tmp2[4];
    double q_inv[4] ;

    dQInv ( BcskWXYZ , q_inv , q ) ;

    dDub( tmp1 , dqdt , 4);
    cblas_dscal ( 4 , 2 , tmp1 , 1 );

    dHamilton ( BcskWXYZ , tmp2 , tmp1 , q_inv ) ;
    v[0] = tmp2[1] ;
    v[1] = tmp2[2] ;
    v[2] = tmp2[3] ;
  
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

int dQV ( double *q , double *v_src , double *v_dst ) {
  double tmp1 [4] ;
  double tmp2 [4] ;
  double tmp3 [4] ;
  double q_inv[4] ;
  // Algorithm:
  // [0;v_dst] = hamilton ( q , hamilton ( [0;v_src] ,  q^-1    ) )
  //                                   __tmp1__    _q_inv_
  //                        _________________tmp2____________

  dQInv( BcskWXYZ ,  q_inv , q ) ;
  tmp1[0] = 0 ;
  tmp1[1] = v_src[0] ;
  tmp1[2] = v_src[1] ;
  tmp1[3] = v_src[2] ;
  dHamilton (BcskWXYZ , tmp2  , tmp1 , q_inv );
  dHamilton (BcskWXYZ , tmp3  , q    , tmp2  );

  v_dst[0] = tmp3[1] ;
  v_dst[1] = tmp3[2] ;
  v_dst[2] = tmp3[3] ;

  return(0);

}




/***********************************************
* Return matrix M from vector v s.t. for any x
* M*x = cross(v,x)
***********************************************/
int dCrossMat ( const enum BCSK_M_ORDER MOrder , double *M , const double *v ) {
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
int dQ2R( const enum BCSK_M_ORDER QOrder , const enum BCSK_Q_ORDER MOrder , double *R , double *q) {

  double eye [9] ;
  dEye(3,1,eye);

  double cross[9];
  dCrossMat ( BcskColMajor , cross , q+1 ) ;

  double tmp1[9];
  dDub(tmp1,cross,9);

  dSum2 ( R , 1 , eye , 2*q[0] , tmp1 , 9 );

  cblas_dgemm ( CblasColMajor , CblasNoTrans , CblasNoTrans , 3 , 3 , 3 , 2 , cross , 3 , cross , 3 , 0 , tmp1 , 3); 

  dSum2 ( R , 1 , R , 1 , tmp1 , 9 );
  
  return (0);
}

// for implementation details see:
// http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
int dR2Q ( const enum BCSK_M_ORDER ROrder , const enum BCSK_Q_ORDER QOrder , double *R , double *q) {
  // Assuming q = [w x y z]'
  // Assuming R in col major order
  double trace = R[0] + R[4] + R[8] ;
  if( trace > 0 ) {
    double s = 0.5 / sqrt( trace + 1.0 );
    q[0] = 0.25 / s;
    q[1] = ( R[5] - R[7] ) * s;
    q[2] = ( R[6] - R[6] ) * s;
    q[3] = ( R[1] - R[3] ) * s;
  } else {
    if ( R[0] > R[4] && R[0] > R[8] ) {
      double s = 2.0 * sqrtf( 1.0 + R[0] - R[4] - R[8]);
      q[0] = (R[5] - R[7] ) / s;
      q[1] = 0.25f * s;
      q[2] = (R[3] + R[1] ) / s;
      q[3] = (R[6] + R[2] ) / s;
    } else if (R[4] > R[8]) {
      double s = 2.0f * sqrtf( 1.0f + R[4] - R[0] - R[8]);
      q[0] = (R[6] - R[2] ) / s;
      q[1] = (R[3] + R[1] ) / s;
      q[2] = 0.25f * s;
      q[3] = (R[7] + R[5] ) / s;
    } else {
      double s = 2.0f * sqrtf( 1.0f + R[8] - R[0] - R[4] );
      q[0] = (R[1] - R[3] ) / s;
      q[1] = (R[6] + R[2] ) / s;
      q[2] = (R[7] + R[5] ) / s;
      q[3] = 0.25f * s;
    }
  }
  dNormalise ( 4 , q );
  return ( 0 ) ;
}


int dSat( double *src , double *dst , double sat_level , int len ) {
  int i;
  if ( sat_level < 0 )
    sat_level = - sat_level;
  for ( i=0 ; i<len ; i++ ){
    if (src[i] > sat_level )
      dst[i] = sat_level ;
    else if (src[i] < -sat_level )
      dst[i] = -sat_level ;
    else
      dst[i] = src[i] ;
  }
  return 0;
}


int dNormalise ( const int N , double *X  ) {
  double norm = cblas_dnrm2( N , X , 1);
  int i;
  for ( i=0 ; i<N ; i++ ) {
    X[i] = X[i] / norm ;
  }

  return 0;
}


