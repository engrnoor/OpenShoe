/*
 * nav_util.h
 *
 * Created: 2016-05-05 14:30:03
 *  Author: jnil02
 */ 


#ifndef NAV_UTIL_H_
#define NAV_UTIL_H_

#include "nav_types.h"

/*! \brief Function for calculating the square root.


	\details Function for calculating the square root using the hardware multipliers if the precision used is float,
	otherwise use the standard square root function.
 */
static inline precision sqrt_hf(precision arg){


	if(sizeof(precision)==4)
	{
		float tmp1, tmp2, tmp3;
		__asm__ __volatile__ ( "frsqrta.s %0, %1" : "=r" (tmp3) : "r" (arg));
		tmp1 = tmp3*tmp3;
		tmp2 = 3.0f - tmp1*arg;
		tmp3 = 0.5f * (tmp2 * tmp3);
		tmp1 = tmp3*tmp3;
		tmp2 = 3.0f - tmp1*arg;
		tmp3 = 0.5f * (tmp2 * tmp3);
		tmp1 = tmp3*tmp3;
		tmp2 = 3.0f - tmp1*arg;
		tmp3 = 0.5f * (tmp2 * tmp3);
		return tmp3*arg;
	}
	else
	{
		return sqrt(arg);
	}

}

static inline precision invsqrt(precision arg){


	if(sizeof(precision)==4)
	{
		float tmp1, tmp2, tmp3;
		__asm__ __volatile__ ( "frsqrta.s %0, %1" : "=r" (tmp3) : "r" (arg));
		tmp1 = tmp3*tmp3;
		tmp2 = 3.0f - tmp1*arg;
		tmp3 = 0.5f * (tmp2 * tmp3);
		tmp1 = tmp3*tmp3;
		tmp2 = 3.0f - tmp1*arg;
		tmp3 = 0.5f * (tmp2 * tmp3);
		tmp1 = tmp3*tmp3;
		tmp2 = 3.0f - tmp1*arg;
		tmp3 = 0.5f * (tmp2 * tmp3);
		return tmp3;
	}
	else
	{
		return 1.0/sqrt(arg);
	}

}

/*! \brief Function that calculates the squared Euclidean norm of a vector.

	 @param[out] norm2				The squared Euclidean norm of the input vector.
	 @param[in] arg_vec				The input vectors
	 @param[in] len					The length of the input vector.
 */
static inline precision vecnorm2(precision *arg_vec, uint8_t len)
{
	precision norm2=0;
	uint8_t ctr;

	for (ctr=0; ctr<len; ctr++){
		norm2=norm2+arg_vec[ctr]*arg_vec[ctr];
	}

	return norm2;
}

/*! \brief Function for converting a rotation matrix \f$R_b^t\f$ to quaternions.

	 @param[out] q					Vector of quaternions.
	 @param[in] rotmat				Vector of representation of the rotation matrix.
 */
inline static void rotation2quat(quat_vec q,const mat3 R){

	// For checking robustness of DCM, diagonal elements
	precision T = (precision)1.0 + R[0]+R[4]+R[8];  // 1+diag(R)

	// Calculate quaternion based on largest diagonal element
	if(T > (0.00000001f))
	{
		precision S = invsqrt(T);// (precision)0.5/ sqrt_hf(T);
		q[3] = S*T;//(precision)0.25 / S;
		q[0] =(R[7]-R[5]) * S;
		q[1] =(R[2]-R[6]) * S;
		q[2] =(R[3]-R[1]) * S;
	}
	else if( (R[0] > R[4] ) && (R[0] > R[8]) )
	{
		T = (precision)1.0 + R[0]-R[4]-R[8];
		precision S = invsqrt(T);//sqrt_hf(T)*2;
		q[3] =(R[7]-R[5])*S;///S;
		q[0] = S*T;//(precision)0.25 * S;
		q[1] = (R[1]+R[3])*S;///S;
		q[2] = (R[2]+R[6])*S;///S;
	}
	else if (R[4]>R[8])
	{
		T = (precision)1.0 + R[4]-R[0]-R[8];
		precision S = invsqrt(T);//sqrt_hf(T)*2;
		q[3] = (R[2]-R[6])*S;///S;
		q[0] =	(R[1]+R[3])*S;///S;
		q[1] = S*T;//(precision)0.25 * S;
		q[2] = (R[5]+R[7])*S;///S;
	}
	else
	{
		T = (precision)1.0 + R[8]-R[0]-R[4];
		precision S = invsqrt(T);//sqrt_hf(T)*2;
		q[3] = (R[3]-R[1])*S;///S;
		q[0] = (R[2]+R[6])*S;///S;
		q[1] = (R[5]+R[7])*S;///S;
		q[2] = S*T;//(precision)0.25 * S;
	}

	//Normalize the quaternion
	T=invsqrt(vecnorm2(q,4));//1.0f/sqrt_hf(vecnorm2(q,4));

	q[0]=q[0]*T;
	q[1]=q[1]*T;
	q[2]=q[2]*T;
	q[3]=q[3]*T;
}

/*! \brief Function for converting quaternions to a rotation matrix \f$R_b^t\f$.

	@param[out] rotmat				Vector of representation of the rotation matrix.
	 @param[in] q					Vector of quaternions.

 */
inline static void quat2rotation(mat3 R,const quat_vec q){

	precision p[6];

	p[0]=q[0]*q[0];
	p[1]=q[1]*q[1];
	p[2]=q[2]*q[2];
	p[3]=q[3]*q[3];

	p[4]=p[1]+p[2];

	if (p[0]+p[3]+p[4] != 0)
	p[5]=(precision)2.0 / (p[0]+p[3]+p[4]);
	else
	p[5]=0;

	R[0]=1-p[5]*p[4];
	R[4]=1-p[5]*(p[0]+p[2]);
	R[8]=1-p[5]*(p[0]+p[1]);

	p[0]=p[5]*q[0];
	p[1]=p[5]*q[1];
	p[4]=p[5]*q[2]*q[3];
	p[5]=p[0]*q[1];

	R[1]=p[5]-p[4];
	R[3]=p[5]+p[4];

	p[4]=p[1]*q[3];
	p[5]=p[0]*q[2];

	R[2]=p[5]+p[4];
	R[6]=p[5]-p[4];

	p[4]=p[0]*q[3];
	p[5]=p[1]*q[2];

	R[5]=p[5]-p[4];
	R[7]=p[5]+p[4];
}


#endif /* NAV_UTIL_H_ */