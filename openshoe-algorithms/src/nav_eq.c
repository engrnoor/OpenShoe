/*! \file nav_eq.c
	\brief The c-file for the OpenShoe navigation algorithm.

	\details This is the c-file for the OpenShoe navigation algorithm. It includes
	the signal processing functions needed to implement a zero-velocity aided inertial navigation system using a nine
	state Kalman filter working in a complimentary feedback configuration. It also includes the functions needed for
	implementing a accelerometer bias calibration framework.

	\authors John-Olof Nilsson, Isaac Skog
 	\copyright Copyright (c) 2011 OpenShoe, ISC License (open source)
 */

///\addtogroup nav_eq
//@{

#include "nav_eq.h"
#include "conf_clock.h"
#include <string.h>
#include "inertial_frontend.h"
#include "nav_util.h"

#if defined(MIMU3333)
#	include "conf_mimu3333.h"
#elif defined(MIMU22BT)
#   include "conf_MIMU22BT.h"
#elif defined(MIMU4444)
#   include "conf_MIMU4444.h"
#elif defined(MIMU4444BT)
#   include "conf_MIMU4444BT.h"
#else
#	error No platform specified!
#endif


/*!
\name General control parameters.

  Parameters controlling general settings of the navigation algorithm

*/
//@{

/// Rough latitude of the system [\f$degrees\f$]. (Used to calculate the magnitude of the gravity vector)
precision latitude=13;

/// Rough altitude of the system [\f$m\f$]. (Used to calculate the magnitude of the gravity vector)
precision altitude=920;

/// Magnitude of the local gravity acceleration [\f$m/s^2\f$]
precision g=9.81f;//9.782940329221166;

/// Error signaling vector. If zero no error has occurred.
uint8_t error_signal;

//@}



/*!
\name Initialization control parameters

  Parameters controlling the initialization of the navigation algorithm,
  i.e., the initial states of the inertial navigation system equations
  and the initial Kalman filter covariance matrix.


*/
//@{
/// A flag that should be set to true when initialization is started and that becomes false when the initialization is finished.
bool init_done=false;

/// Initial heading [\f$rad\f$]
#define INITIAL_HEADING 0.0f

/// Initial position (North, East, Down) [\f$m\f$]
#define INITIAL_POSX 0.0f
#define INITIAL_POSY 0.0f
#define INITIAL_POSZ 0.0f

/// Standard deviations in the initial position uncertainties [\f$m\f$].
#define SIGMA2_INIT_POSX (0.00001f*0.00001f)
#define SIGMA2_INIT_POSY (0.00001f*0.00001f)
#define SIGMA2_INIT_POSZ (0.00001f*0.00001f)

/// Standard deviations in the initial velocity uncertainties [\f$m/s\f$].
#define SIGMA2_INIT_VELX (0.01f*0.01f)
#define SIGMA2_INIT_VELY (0.01f*0.01f)
#define SIGMA2_INIT_VELZ (0.01f*0.01f)

/// Standard deviations in the initial attitude uncertainties [\f$rad\f$].
#define SIGMA2_INIT_ROLL  (0.00174f*0.00174f)
#define SIGMA2_INIT_PITCH (0.00174f*0.00174f)
#define SIGMA2_INIT_YAW   (0.00174f*0.00174f)
//@}


/*!
\name Kalman filter control parameters

  Parameters controlling the behavior of the Kalman filter.
  The parameters can be changed while the filter is running
  to adapt the filter to the current motion dynamics.

  \note The default noise standard deviation figures are not set to
  reflect the true noise figures of the IMU sensors, but
  rather to model the sum of all the errors (biases,
  scale factors, nonlinearities, etc.) in the system and
  the measurement model.

*/
//@{
/// Accelerometer process noise standard deviation [\f$m/s^2\f$]
//precision sigma_acceleration=KALMAN_SIGMA_2_ACC;

/// Gyroscope process noise standard deviation [\f$rad/s\f$]
//precision sigma_gyroscope=KALMAN_SIGMA_2_GYRO;

/// Pseudo zero-velocity measurement noise standard deviations (north, east, down) [\f$m/s\f$]
//vec3 sigma_velocity={KALMAN_SIGMA_2_ZUPT,KALMAN_SIGMA_2_ZUPT,KALMAN_SIGMA_2_ZUPT};
//@}




/*!
\name Navigation and filter state variables.

  Vectors that holds the current navigation state estimate and the covariance and gain of the Kalman filter.

*/
//@{
///  Position estimate (North,East,Down) [\f$m\f$].
vec3 pos;
vec3 dpos;

/// Velocity estimate (North,East,Down) [\f$m/s\f$]
vec3 vel;
vec3 dvel;

/// Attitude (quaternions) estimate
quat_vec quat;
vec3 deuler;

/// Rotation matrix used as an "aiding" variable in the filter algorithm. Holds the same information as the quaternions.
mat3 Rb2t;

/// Vector representation of the Kalman filter covariance matrix.
mat9sym P;

/// Vector representation of the Kalman filter gain matrix.
mat9by3 K;
//@}

//@}


/*!
\name Step-wise dead-reckoning parameters and variables

  Parameters and variables controlling the resets of the step-wise dead-reckoning mode.

*/
//@{

/// Covariance threshold for filter reset
precision reset_cov_threshold = RESET_COV_THRESHOLD;

/// Minimum period between filter resets (time)
float min_time_between_resets = 0.75f; //[s]

/// Flag signaling that the filter has been reset
bool filter_reset_flag = false;

/// Maximum stationary period before a pending filter reset is done
float max_time_reset_pending = 0.25f; //[s]

/// Displacement and heading change vector
vec4 dx;

/// Displacement and heading change covariance vector
mat4sym dP;

uint16_t step_counter = 0;

//@}




//******************* INLINE FUNCTIONS ******************************//

/**
	\defgroup aux_func Auxilliary functions
	\brief Auxilliary functions for calculation various matrix and scalar operations.

	@{
*/

/*! \brief Function that converts Euler angles ([roll,pitch,yaw]) into a rotation matrix \f$R_b^t\f$.

	 @param[out] rotmat				Vector representation of the rotation matrix.
	 @param[in] euler				Vector of euler angles.
 */
inline void euler2rotation(mat3 R, const vec3 euler){

	// Trigonometric value variables
	precision sin_phi = sin(euler[0]);
	precision cos_phi = cos(euler[0]);
	precision sin_theta =sin(euler[1]);
	precision cos_theta =cos(euler[1]);
	precision sin_psi = sin(euler[2]);
	precision cos_psi = cos(euler[2]);

	R[0]=cos_psi*cos_theta;
	R[3]=sin_psi*cos_theta;
	R[6]=-sin_theta;
	R[1]=(-sin_psi*cos_phi) + cos_psi*(sin_theta*sin_phi);
	R[4]=(cos_psi*cos_phi) + sin_psi*(sin_theta*sin_phi);
	R[7]=cos_theta*sin_phi;
	R[2]=(sin_psi*sin_phi) + cos_psi*(sin_theta*cos_phi);
	R[5]=(-cos_psi*sin_phi) + sin_psi*(sin_theta*cos_phi);
	R[8]=cos_theta*cos_phi;
}

/*! \brief Function for converting a rotation matrix \f$R_b^t\f$ to Euler angles ([roll,pitch,yaw]).

	@param[out] euler				Vector of euler angles.
	@param[in] rotmat				Vector of representation of the rotation matrix.
 */
inline void  rotation2euler(vec3 euler, const mat3 R){

	// Compute Euler angles [WARNING! check atan2]
	euler[0] = atan2(R[7],R[8]);	//atan2( R(3,2), R(3,3) );
	euler[1] = asin(-R[6]);			//asin( -R(3,1) );
	euler[2] = -atan2(R[3],R[0]);	//atan2( R(2,1), R(1,1));
}




/*! \brief Function for calculating the Kalman filter innovation covariance.

	@param[out] re		Vector representation of the innovation covariance matrix.
	@param[in] pvec		Vector representation of the Kalman filter covariance matrix.
	@param[in] sigma	Vector representation of pseudo zero-velocity measurement noise standard deviations.
 */
inline void innovation_cov(mat3sym re,mat9sym pvec){
	re[0]=pvec[24]+KALMAN_SIGMA_2_ZUPT; // x
	re[1]=pvec[25];
	re[2]=pvec[26];
	re[3]=pvec[30]+KALMAN_SIGMA_2_ZUPT; // y
	re[4]=pvec[31];
	re[5]=pvec[35]+KALMAN_SIGMA_2_ZUPT; // z
}




/*! \brief Function for inverting a 3 by 3 matrix hermitian matrix.

	@param[out] ainv		Vector representation of the inverted matrix.
	@param[out] error		Error signaling vector.
	@param[in]	a			Vector representation of the matrix to be inverted.
 */
inline void invmat3sys(mat3sym ainv, mat3sym a){

	// Calculate the determinant of matrix
	precision det=-a[2]*(a[2]*a[3]) + 2.0f*a[1]*(a[2]*a[4]) - a[0]*(a[4]*a[4]) - a[1]*(a[1]*a[5]) + a[0]*(a[3]*a[5]);

	// Check the size of the determinant and send a error message if it is to small.
	if(absf(det)<0.0000001f)
	{
		error_signal=MATRIX_INVERSION_ERROR;
	}

	precision invdet = (precision)1.0/det;
	ainv[0]=(a[3]*a[5]-a[4]*a[4])*invdet;
	ainv[1]=(a[2]*a[4]-a[1]*a[5])*invdet;
	ainv[2]=(a[1]*a[4]-a[2]*a[3])*invdet;
	ainv[3]=(a[0]*a[5]-a[2]*a[2])*invdet;
	ainv[4]=(a[1]*a[2]-a[0]*a[4])*invdet;
	ainv[5]=(a[0]*a[3]-a[1]*a[1])*invdet;
}




/*! \brief	Function that calculates the maximum value of a vector and returns
			the max value and the index of the vector element holding the maximum value.

	@param[out] max_v		Largest value of the input vector.
	@param[out] index		Index of the vector element holding the largest value.
	@param[in]	arg_vec		The input vector.
 */
inline void max_value(precision *max_v,uint8_t *index, precision *arg_vec){

	// Set the initial max value and vector element index
	*max_v=arg_vec[0];
	*index=0;

	//Iterate through the vector
	for(uint8_t ctr=1;ctr<(sizeof(arg_vec)/sizeof(arg_vec[0]));ctr++)
	{
		//If the current element of the vector is larger than the previous element, update the max value and the index.
		if(arg_vec[ctr]>arg_vec[ctr-1])
		{

			*max_v=arg_vec[ctr];
			*index=ctr;
		}
	}
}

//@}

/**
	\defgroup miscel Miscellaneous
	\brief Miscellaneous routines.

	@{
*/

/*! \brief Function that calculates the magnitude of the local gravity vector based upon the WGS84 gravity model.

	 @param[out] g			The magnitude of the local gravity vector.
	 @param[in] latitude	The latitude of the navigation system.
	 @param[in] altitude	The altitude of the navigation system.

 */
inline void gravity(void){

	precision lambda=M_PI/180.0*latitude;  //latitude [rad]

	g=9.780327*(1+0.0053024*(sin(lambda)*sin(lambda))-0.0000058*(sin(2*lambda)*sin(2*lambda)))-(0.0000030877-0.000000004*(sin(lambda)*sin(lambda)))*altitude+0.000000000000072*(altitude*altitude);

}

//@}

/**
	\defgroup main_func ZUPT aided INS functions
	\brief the main functions for running the ZUPT aided INS.

	@{
*/

//********************* ORDINARY FUNCTIONS *********************************//

void strapdown_mechanisation_equations(void){
	/*
		The inputs and outputs are the following:

		position - position vector in the n-frame [m]
		velocity - velocity vector in the n-frame [m]
		quaternions - quaternions describing the orientation of the navigation platform
		Rb2t - The rotation matrix Rb2t (Includes the same information as the quaternionsuaternions but we use it other functions, i.e., we saves computations)
		accelerations - Measured accelerations vector (acceleration) in the b-frame [m/s^2]
		angular_rates - Measured angular-rate vector in the b-frame [rad/s]

		The following global variables most be defined in the main function:

		g - The gravity constant [m/s^2]
		dt - The sampling period [s]

   Note that this is a rudimentary mechanization of the navigation equaternionsuations which only is suitable for
   use with systems that uses low-cost/low-performance inertial sensors and where only short periods of
   free inertial navigation is expected.
*/



	// Working variables
	vec3 an_hat;
    vec3 angular_rates_dt;
	quat_vec quat_tmp;
	precision cos_v;
	precision sin_v;


	// Convert quaternions to rotation matrix
	quat2rotation(Rb2t,quat);  //Rb2t

	// Compute acceleration in navigation coordinate frame and subtract the acceleration due to the earth gravity force.
	an_hat[0]=Rb2t[0]*u_k.f.x+Rb2t[1]*u_k.f.y+Rb2t[2]*u_k.f.z;
	an_hat[1]=Rb2t[3]*u_k.f.x+Rb2t[4]*u_k.f.y+Rb2t[5]*u_k.f.z;
	an_hat[2]=Rb2t[6]*u_k.f.x+Rb2t[7]*u_k.f.y+Rb2t[8]*u_k.f.z+g;
	
	//******** Update the position and velocity *******//

	// Integrate the velocity to get the position
	if(!zaru){
		pos[0]=pos[0]+vel[0]*dt_k;
		pos[1]=pos[1]+vel[1]*dt_k;
		pos[2]=pos[2]+vel[2]*dt_k;
	}

	// Integrate the acceleration to get the velocity
	vel[0]=vel[0]+an_hat[0]*dt_k;
	vel[1]=vel[1]+an_hat[1]*dt_k;
	vel[2]=vel[2]+an_hat[2]*dt_k;
	//******************************************************//

	//********** Integrate angular rates  **********//

	// If we are measuring some angular rates, update the quaternion vector.
    if((u_k.g.x!=0)|(u_k.g.y!=0)|(u_k.g.z!=0))
	{

	// Multiply the angular rates with the sampling period dt.
	angular_rates_dt[0]=u_k.g.x*dt_k;
	angular_rates_dt[1]=u_k.g.y*dt_k;
	angular_rates_dt[2]=u_k.g.z*dt_k;

	if(zaru){
//		quat2rotation(Rb2t,quat);
		precision d1 = Rb2t[6]*Rb2t[6]*angular_rates_dt[0] + Rb2t[6]*Rb2t[7]*angular_rates_dt[1] + Rb2t[6]*Rb2t[8]*angular_rates_dt[2];
		precision d2 = Rb2t[7]*Rb2t[6]*angular_rates_dt[0] + Rb2t[7]*Rb2t[7]*angular_rates_dt[1] + Rb2t[7]*Rb2t[8]*angular_rates_dt[2];
		precision d3 = Rb2t[8]*Rb2t[6]*angular_rates_dt[0] + Rb2t[8]*Rb2t[7]*angular_rates_dt[1] + Rb2t[8]*Rb2t[8]*angular_rates_dt[2];
		angular_rates_dt[0]-=d1;
		angular_rates_dt[1]-=d2;
		angular_rates_dt[2]-=d3;
	}

	// Calculate the norm of the vector angular_rates_dt
	precision v=( sqrt_hf( vecnorm2(angular_rates_dt, 3) ) );

	cos_v=cos(v/2);
	sin_v=(sin(v/2)/v);

	// Time update of the quaternions
	quat_tmp[0]=cos_v*quat[0]+sin_v*(angular_rates_dt[2]*quat[1]-angular_rates_dt[1]*quat[2]+angular_rates_dt[0]*quat[3]);	// w_tb(2)*quaternions(1)-w_tb(1)*quaternions(2)+w_tb(0)*quaternions(3)
	quat_tmp[1]=cos_v*quat[1]+sin_v*(-angular_rates_dt[2]*quat[0]+angular_rates_dt[0]*quat[2]+angular_rates_dt[1]*quat[3]);  //-w_tb(2)*quaternions(0)+w_tb(0)*quaternions(2)+w_tb(1)*quaternions(3)
	quat_tmp[2]=cos_v*quat[2]+sin_v*(angular_rates_dt[1]*quat[0]-angular_rates_dt[0]*quat[1]+angular_rates_dt[2]*quat[3]);   //w_tb(1)*quaternions(0)-w_tb(0)*quaternions(1)+w_tb(2)*quaternions(3)
	quat_tmp[3]=cos_v*quat[3]+sin_v*(-angular_rates_dt[0]*quat[0]-angular_rates_dt[1]*quat[1]-angular_rates_dt[2]*quat[2]);  //-w_tb(0)*quaternions(0)-w_tb(1)*quaternions(1)-w_tb(2)*quaternions(2)



	// Re-normalize the quaternions and update the global variable
	v=sqrt_hf(vecnorm2(quat_tmp, 4));
	quat[0]=quat_tmp[0]/v;
	quat[1]=quat_tmp[1]/v;
	quat[2]=quat_tmp[2]/v;
	quat[3]=quat_tmp[3]/v;
	}
	//*****************************************************//
}


void time_up_data(void){

	//Working variables
	precision dt2_sigma2_acc= (dt_k*dt_k)*KALMAN_SIGMA_2_ACC;
	precision dt2_sigma2_gyro=(dt_k*dt_k)*KALMAN_SIGMA_2_GYRO;
//	mat9sym Pp;		//Temporary vector holding the update covariances
	vec3 s;				//Specific accelerations vector in the n-frame.

	// Convert quaternions to rotation matrix
//	quat2rotation(Rb2t,quat);  //Rb2t
	
	s[0]=(Rb2t[0]*u_k.f.x+Rb2t[1]*u_k.f.y+Rb2t[2]*u_k.f.z)*dt_k;
	s[1]=(Rb2t[3]*u_k.f.x+Rb2t[4]*u_k.f.y+Rb2t[5]*u_k.f.z)*dt_k;
	s[2]=(Rb2t[6]*u_k.f.x+Rb2t[7]*u_k.f.y+Rb2t[8]*u_k.f.z)*dt_k;
	
	if(!zaru){
		P[0]+=P[24]*dt_k*dt_k + (2.0f*P[3])*dt_k;
		P[1]+=P[25]*dt_k*dt_k + (P[4] + P[11])*dt_k;
		P[2]+=P[26]*dt_k*dt_k + (P[5] + P[18])*dt_k;
		P[3]+=(P[24] - P[28]*s[2] + P[29]*s[1])*dt_k + (P[8]*s[1] - P[7]*s[2]);
		P[4]+=(P[25] + P[27]*s[2] - P[29]*s[0])*dt_k + (P[6]*s[2] - P[8]*s[0]);
		P[5]+=(P[26] - P[27]*s[1] + P[28]*s[0])*dt_k + (P[7]*s[0] - P[6]*s[1]);
		P[6]+=P[27]*dt_k;
		P[7]+=P[28]*dt_k;
		P[8]+=P[29]*dt_k;
		P[9]+=P[30]*dt_k*dt_k + (2.0f*P[12])*dt_k;
		P[10]+=P[31]*dt_k*dt_k + (P[13] + P[19])*dt_k;
		P[11]+=(P[25] - P[33]*s[2] + P[34]*s[1])*dt_k + (P[16]*s[1] - P[15]*s[2]);
		P[12]+=(P[30] + P[32]*s[2] - P[34]*s[0])*dt_k + (P[14]*s[2] - P[16]*s[0]);
		P[13]+=(P[31] - P[32]*s[1] + P[33]*s[0])*dt_k + (P[15]*s[0] - P[14]*s[1]);
		P[14]+=P[32]*dt_k;
		P[15]+=P[33]*dt_k;
		P[16]+=P[34]*dt_k;
		P[17]+=P[35]*dt_k*dt_k + (2.0f*P[20])*dt_k;
		P[18]+=(P[26] - P[37]*s[2] + P[38]*s[1])*dt_k + (P[23]*s[1] - P[22]*s[2]);
		P[19]+=(P[31] + P[36]*s[2] - P[38]*s[0])*dt_k + (P[21]*s[2] - P[23]*s[0]);
		P[20]+=(P[35] - P[36]*s[1] + P[37]*s[0])*dt_k + (P[22]*s[0] - P[21]*s[1]);
		P[21]+=P[36]*dt_k;
		P[22]+=P[37]*dt_k;
		P[23]+=P[38]*dt_k;
		P[24]+=P[44]*s[1]*s[1] - 2.0f*P[43]*s[1]*s[2] + 2.0f*P[29]*s[1] + P[42]*s[2]*s[2] - 2.0f*P[28]*s[2] + dt2_sigma2_acc;
		P[25]+=P[34]*s[1] - P[33]*s[2] + s[2]*(P[27] - P[40]*s[2] + P[41]*s[1]) - s[0]*(P[29] - P[43]*s[2] + P[44]*s[1]);
		P[26]+=P[38]*s[1] - P[37]*s[2] - s[1]*(P[27] - P[40]*s[2] + P[41]*s[1]) + s[0]*(P[28] - P[42]*s[2] + P[43]*s[1]);
		P[27]+=P[41]*s[1] - P[40]*s[2];
		P[28]+=P[43]*s[1] - P[42]*s[2];
		P[29]+=P[44]*s[1] - P[43]*s[2];
		P[30]+=P[44]*s[0]*s[0] - 2.0f*P[41]*s[0]*s[2] - 2.0f*P[34]*s[0] + P[39]*s[2]*s[2] + 2.0f*P[32]*s[2] + dt2_sigma2_acc;
		P[31]+=P[36]*s[2] - P[38]*s[0] - s[1]*(P[32] + P[39]*s[2] - P[41]*s[0]) + s[0]*(P[33] + P[40]*s[2] - P[43]*s[0]);
		P[32]+=P[39]*s[2] - P[41]*s[0];
		P[33]+=P[40]*s[2] - P[43]*s[0];
		P[34]+=P[41]*s[2] - P[44]*s[0];
		P[35]+=P[42]*s[0]*s[0] - 2.0f*P[40]*s[0]*s[1] + 2.0f*P[37]*s[0] + P[39]*s[1]*s[1] - 2.0f*P[36]*s[1] + dt2_sigma2_acc;
		P[36]+=P[40]*s[0] - P[39]*s[1];
		P[37]+=P[42]*s[0] - P[40]*s[1];
		P[38]+=P[43]*s[0] - P[41]*s[1];
		P[39]+=dt2_sigma2_gyro;
		P[40]+=0;
		P[41]+=0;
		P[42]+=dt2_sigma2_gyro;
		P[43]+=0;
		P[44]+=dt2_sigma2_gyro;
	} else {
		P[0]+=0;
		P[1]+=0;
		P[2]+=0;
		P[3]+=P[8]*s[1] - P[7]*s[2];
		P[4]+=P[6]*s[2] - P[8]*s[0];
		P[5]+=P[7]*s[0] - P[6]*s[1];
		P[6]+=0;
		P[7]+=0;
		P[8]+=0;
		P[9]+=0;
		P[10]+=0;
		P[11]+=P[16]*s[1] - P[15]*s[2];
		P[12]+=P[14]*s[2] - P[16]*s[0];
		P[13]+=P[15]*s[0] - P[14]*s[1];
		P[14]+=0;
		P[15]+=0;
		P[16]+=0;
		P[17]+=0;
		P[18]+=P[23]*s[1] - P[22]*s[2];
		P[19]+=P[21]*s[2] - P[23]*s[0];
		P[20]+=P[22]*s[0] - P[21]*s[1];
		P[21]+=0;
		P[22]+=0;
		P[23]+=0;
		P[24]+=dt2_sigma2_acc - 2.0f*P[28]*s[2] + 2.0f*P[29]*s[1] + P[42]*s[2]*s[2] + P[44]*s[1]*s[1] - 2.0f*P[43]*s[1]*s[2];
		P[25]+=P[34]*s[1] - P[33]*s[2] + s[2]*(P[27] - P[40]*s[2] + P[41]*s[1]) - s[0]*(P[29] - P[43]*s[2] + P[44]*s[1]);
		P[26]+=P[38]*s[1] - P[37]*s[2] - s[1]*(P[27] - P[40]*s[2] + P[41]*s[1]) + s[0]*(P[28] - P[42]*s[2] + P[43]*s[1]);
		P[27]+=P[41]*s[1] - P[40]*s[2];
		P[28]+=P[43]*s[1] - P[42]*s[2];
		P[29]+=P[44]*s[1] - P[43]*s[2];
		P[30]+=dt2_sigma2_acc + 2.0f*P[32]*s[2] - 2.0f*P[34]*s[0] + P[39]*s[2]*s[2] + P[44]*s[0]*s[0] - 2.0f*P[41]*s[0]*s[2];
		P[31]+=P[36]*s[2] - P[38]*s[0] - s[1]*(P[32] + P[39]*s[2] - P[41]*s[0]) + s[0]*(P[33] + P[40]*s[2] - P[43]*s[0]);
		P[32]+=P[39]*s[2] - P[41]*s[0];
		P[33]+=P[40]*s[2] - P[43]*s[0];
		P[34]+=P[41]*s[2] - P[44]*s[0];
		P[35]+=dt2_sigma2_acc - 2.0f*P[36]*s[1] + 2.0f*P[37]*s[0] + P[39]*s[1]*s[1] + P[42]*s[0]*s[0] - 2.0f*P[40]*s[0]*s[1];
		P[36]+=P[40]*s[0] - P[39]*s[1];
		P[37]+=P[42]*s[0] - P[40]*s[1];
		P[38]+=P[43]*s[0] - P[41]*s[1];
		P[39]+=dt2_sigma2_gyro;
		P[40]+=0;
		P[41]+=0;
		P[42]+=dt2_sigma2_gyro;
		P[43]+=0;
		P[44]+=0;
	}
}



inline void gain_matrix(void){

mat3sym Re;			//Innovation matrix
mat3sym invRe;		//Inverse of the innovation matrix

/************ Calculate the Kalman filter innovation matrix *******************/
innovation_cov(Re,P);

/************ Calculate the inverse of the innovation matrix *********/
invmat3sys(invRe,Re);

/******************* Calculate the Kalman filter gain **************************/
// First row of the gain matrix
K[0]=P[3]*invRe[0] + P[4]*invRe[1] + P[5]*invRe[2];
K[1]=P[3]*invRe[1] + P[4]*invRe[3] + P[5]*invRe[4];
K[2]=P[3]*invRe[2] + P[4]*invRe[4] + P[5]*invRe[5];
// Second row of the gain matrix
K[3]=P[11]*invRe[0] + P[12]*invRe[1] + P[13]*invRe[2];
K[4]=P[11]*invRe[1] + P[12]*invRe[3] + P[13]*invRe[4];
K[5]=P[11]*invRe[2] + P[12]*invRe[4] + P[13]*invRe[5];
// Third row of the gain matrix
K[6]=P[18]*invRe[0] + P[19]*invRe[1] + P[20]*invRe[2];
K[7]=P[18]*invRe[1] + P[19]*invRe[3] + P[20]*invRe[4];
K[8]=P[18]*invRe[2] + P[19]*invRe[4] + P[20]*invRe[5];
//Forth row of the gain matrix
K[9] =P[24]*invRe[0] + P[25]*invRe[1] + P[26]*invRe[2];
K[10]=P[24]*invRe[1] + P[25]*invRe[3] + P[26]*invRe[4];
K[11]=P[24]*invRe[2] + P[25]*invRe[4] + P[26]*invRe[5];
//Fifth row of the gain matrix
K[12]=P[25]*invRe[0] + P[30]*invRe[1] + P[31]*invRe[2];
K[13]=P[25]*invRe[1] + P[30]*invRe[3] + P[31]*invRe[4];
K[14]=P[25]*invRe[2] + P[30]*invRe[4] + P[31]*invRe[5];
//Sixth row of the gain matrix
K[15]=P[26]*invRe[0] + P[31]*invRe[1] + P[35]*invRe[2];
K[16]=P[26]*invRe[1] + P[31]*invRe[3] + P[35]*invRe[4];
K[17]=P[26]*invRe[2] + P[31]*invRe[4] + P[35]*invRe[5];
//Seventh row of the gain matrix
K[18]=P[27]*invRe[0] + P[32]*invRe[1] + P[36]*invRe[2];
K[19]=P[27]*invRe[1] + P[32]*invRe[3] + P[36]*invRe[4];
K[20]=P[27]*invRe[2] + P[32]*invRe[4] + P[36]*invRe[5];
//Eight row of the gain matrix
K[21]=P[28]*invRe[0] + P[33]*invRe[1] + P[37]*invRe[2];
K[22]=P[28]*invRe[1] + P[33]*invRe[3] + P[37]*invRe[4];
K[23]=P[28]*invRe[2] + P[33]*invRe[4] + P[37]*invRe[5];
//Ninth row of the gain matrix
K[24]=P[29]*invRe[0] + P[34]*invRe[1] + P[38]*invRe[2];
K[25]=P[29]*invRe[1] + P[34]*invRe[3] + P[38]*invRe[4];
K[26]=P[29]*invRe[2] + P[34]*invRe[4] + P[38]*invRe[5];
}


inline void measurement_update(void){

mat9sym Pp;		//Temporary vector holding the update covariances

// First row
Pp[0]=P[0] - K[0]*P[3]  - K[1]*P[4]  - K[2]*P[5];
Pp[1]=P[1] - K[0]*P[11] - K[1]*P[12] - K[2]*P[13];
Pp[2]=P[2] - K[0]*P[18] - K[1]*P[19] - K[2]*P[20];
Pp[3]=P[3] - K[0]*P[24] - K[1]*P[25] - K[2]*P[26];
Pp[4]=P[4] - K[0]*P[25] - K[1]*P[30] - K[2]*P[31];
Pp[5]=P[5] - K[0]*P[26] - K[1]*P[31] - K[2]*P[35];
Pp[6]=P[6] - K[0]*P[27] - K[1]*P[32] - K[2]*P[36];
Pp[7]=P[7] - K[0]*P[28] - K[1]*P[33] - K[2]*P[37];
Pp[8]=P[8] - K[0]*P[29] - K[1]*P[34] - K[2]*P[38];
// Second row
Pp[9] =P[9]  - K[3]*P[11] - K[4]*P[12] - K[5]*P[13];
Pp[10]=P[10] - K[3]*P[18] - K[4]*P[19] - K[5]*P[20];
Pp[11]=P[11] - K[3]*P[24] - K[4]*P[25] - K[5]*P[26];
Pp[12]=P[12] - K[3]*P[25] - K[4]*P[30] - K[5]*P[31];
Pp[13]=P[13] - K[3]*P[26] - K[4]*P[31] - K[5]*P[35];
Pp[14]=P[14] - K[3]*P[27] - K[4]*P[32] - K[5]*P[36];
Pp[15]=P[15] - K[3]*P[28] - K[4]*P[33] - K[5]*P[37];
Pp[16]=P[16] - K[3]*P[29] - K[4]*P[34] - K[5]*P[38];
// Third row
Pp[17]=P[17] - K[6]*P[18] - K[7]*P[19] - K[8]*P[20];
Pp[18]=P[18] - K[6]*P[24] - K[7]*P[25] - K[8]*P[26];
Pp[19]=P[19] - K[6]*P[25] - K[7]*P[30] - K[8]*P[31];
Pp[20]=P[20] - K[6]*P[26] - K[7]*P[31] - K[8]*P[35];
Pp[21]=P[21] - K[6]*P[27] - K[7]*P[32] - K[8]*P[36];
Pp[22]=P[22] - K[6]*P[28] - K[7]*P[33] - K[8]*P[37];
Pp[23]=P[23] - K[6]*P[29] - K[7]*P[34] - K[8]*P[38];
// Forth row
Pp[24]=P[24] - K[9]*P[24] - K[10]*P[25] - K[11]*P[26];
Pp[25]=P[25] - K[9]*P[25] - K[10]*P[30] - K[11]*P[31];
Pp[26]=P[26] - K[9]*P[26] - K[10]*P[31] - K[11]*P[35];
Pp[27]=P[27] - K[9]*P[27] - K[10]*P[32] - K[11]*P[36];
Pp[28]=P[28] - K[9]*P[28] - K[10]*P[33] - K[11]*P[37];
Pp[29]=P[29] - K[9]*P[29] - K[10]*P[34] - K[11]*P[38];
// Fifth row
Pp[30]=P[30] - K[12]*P[25] - K[13]*P[30] - K[14]*P[31];
Pp[31]=P[31] - K[12]*P[26] - K[13]*P[31] - K[14]*P[35];
Pp[32]=P[32] - K[12]*P[27] - K[13]*P[32] - K[14]*P[36];
Pp[33]=P[33] - K[12]*P[28] - K[13]*P[33] - K[14]*P[37];
Pp[34]=P[34] - K[12]*P[29] - K[13]*P[34] - K[14]*P[38];
// Sixth row
Pp[35]=P[35] - K[15]*P[26] - K[16]*P[31] - K[17]*P[35];
Pp[36]=P[36] - K[15]*P[27] - K[16]*P[32] - K[17]*P[36];
Pp[37]=P[37] - K[15]*P[28] - K[16]*P[33] - K[17]*P[37];
Pp[38]=P[38] - K[15]*P[29] - K[16]*P[34] - K[17]*P[38];
// Seventh row
Pp[39]=P[39] - K[18]*P[27] - K[19]*P[32] - K[20]*P[36];
Pp[40]=P[40] - K[18]*P[28] - K[19]*P[33] - K[20]*P[37];
Pp[41]=P[41] - K[18]*P[29] - K[19]*P[34] - K[20]*P[38];
// Eight row
Pp[42]=P[42] - K[21]*P[28] - K[22]*P[33] - K[23]*P[37];
Pp[43]=P[43] - K[21]*P[29] - K[22]*P[34] - K[23]*P[38];
// Ninth row
Pp[44]=P[44] - K[24]*P[29] - K[25]*P[34] - K[26]*P[38];

memcpy(P,Pp,sizeof(P));
}



inline void correct_navigation_states(void){

// Calculate corrections
dpos[0]=K[0]*vel[0]+K[1]*vel[1]+K[2]*vel[2];
dpos[1]=K[3]*vel[0]+K[4]*vel[1]+K[5]*vel[2];
dpos[2]=K[6]*vel[0]+K[7]*vel[1]+K[8]*vel[2];
dvel[0]=K[9]*vel[0]+K[10]*vel[1]+K[11]*vel[2];
dvel[1]=K[12]*vel[0]+K[13]*vel[1]+K[14]*vel[2];
dvel[2]=K[15]*vel[0]+K[16]*vel[1]+K[17]*vel[2];
deuler[0]=K[18]*vel[0]+K[19]*vel[1]+K[20]*vel[2]; //roll
deuler[1]=K[21]*vel[0]+K[22]*vel[1]+K[23]*vel[2]; //pitch
deuler[2]=K[24]*vel[0]+K[25]*vel[1]+K[26]*vel[2]; //yaw

// Correct the position and velocity
pos[0]-=dpos[0];
pos[1]-=dpos[1];
pos[2]-=dpos[2];
vel[0]-=dvel[0];
vel[1]-=dvel[1];
vel[2]-=dvel[2];

// Correct the rotation matrix
mat3 new_R;
mat3 Rb2t_updated;
quat2rotation(Rb2t_updated,quat);
new_R[0]=Rb2t_updated[0]-deuler[2]*Rb2t_updated[3] + deuler[1]*Rb2t_updated[6];
new_R[1]=Rb2t_updated[1]-deuler[2]*Rb2t_updated[4] + deuler[1]*Rb2t_updated[7];
new_R[2]=Rb2t_updated[2]-deuler[2]*Rb2t_updated[5] + deuler[1]*Rb2t_updated[8];
new_R[3]=deuler[2]*Rb2t_updated[0] + Rb2t_updated[3] - deuler[0]*Rb2t_updated[6];
new_R[4]=deuler[2]*Rb2t_updated[1] + Rb2t_updated[4] - deuler[0]*Rb2t_updated[7];
new_R[5]=deuler[2]*Rb2t_updated[2] + Rb2t_updated[5] - deuler[0]*Rb2t_updated[8];
new_R[6]=-deuler[1]*Rb2t_updated[0] + deuler[0]*Rb2t_updated[3] + Rb2t_updated[6];
new_R[7]=-deuler[1]*Rb2t_updated[1] + deuler[0]*Rb2t_updated[4] + Rb2t_updated[7];
new_R[8]=-deuler[1]*Rb2t_updated[2] + deuler[0]*Rb2t_updated[5] + Rb2t_updated[8];

// Calculate the corrected quaternions
rotation2quat(quat,new_R);
}

//@}

/**
	\defgroup init Initialization routines
	\brief Routines for initializing the system. Only coarse initial alignment is implemented

	@{
*/


void frontend_initial_alignment(void){

	// Counter that counts the number of initialization samples that have been processed.
	static uint32_t initialize_sample_ctr = 0;

	// Mean acceleration and angular rate vectors used in the initial alignment.
	static vec3 acceleration_mean;

	// Make sure we capture a reasonably stable start period AND THAT THE BUFFERS ARE FILLED UP
	if (!zupt)
		initialize_sample_ctr=0;

	// Reset the mean if we start a new initialization
	if (initialize_sample_ctr==0)
	{
		acceleration_mean[0]=0;
		acceleration_mean[1]=0;
		acceleration_mean[2]=0;
	}

	//Calculate the cumulative sum of the IMU readings until we have the number
	//of samples specified in the initial alignment settings
	if(initialize_sample_ctr<NR_OF_INITIAL_ALIGNMENT_SAMPLES && zupt){

		//Cumulative sum
		acceleration_mean[0]+=u_k.f.x;
		acceleration_mean[1]+=u_k.f.y;
		acceleration_mean[2]+=u_k.f.z;

		//Increase the counter
		initialize_sample_ctr++;
	}


	// If we are at the last iteration of the initial alignment, do this:
	if(initialize_sample_ctr==NR_OF_INITIAL_ALIGNMENT_SAMPLES){

		vec3 initial_attitude;
		/************* Initialize the navigation states *************/

		// Calculate the mean acceleration and angular rates from the cumulative sums.
		const precision tmp = 1.0f/NR_OF_INITIAL_ALIGNMENT_SAMPLES;
		acceleration_mean[0]*=tmp;
		acceleration_mean[1]*=tmp;
		acceleration_mean[2]*=tmp;

		//Calculate the roll and pitch
		initial_attitude[0]=atan2(-acceleration_mean[1],-acceleration_mean[2]);		//roll
		initial_attitude[1]=atan2(acceleration_mean[0],sqrt_hf((acceleration_mean[1]*acceleration_mean[1])+(acceleration_mean[2]*acceleration_mean[2]))); //pitch

		//Set the initial heading
		initial_attitude[2]=INITIAL_HEADING;

		// Calculate the initial rotation matrix, used as temporary variable in the calculation of the initial quaternions
		mat3 initial_rotmat;
		euler2rotation(initial_rotmat,initial_attitude);

		//Set the initial quaternions using the initial rotation matrix
		rotation2quat(quat,initial_rotmat);

		//Set the initial velocity (assumed zero since we have a zupt)
		vel[0]=0.0f;
		vel[1]=0.0f;
		vel[2]=0.0f;

		//Set the initial position
		pos[0]=INITIAL_POSX;
		pos[1]=INITIAL_POSY;
		pos[2]=INITIAL_POSZ;

		/*************************************************************/

		/************** Initialize the filter covariance *************/
		memset(P,0,sizeof(P));
		P[0]=SIGMA2_INIT_POSX;
		P[9]=SIGMA2_INIT_POSY;
		P[17]=SIGMA2_INIT_POSZ;

		P[24]=SIGMA2_INIT_VELX;
		P[30]=SIGMA2_INIT_VELY;
		P[35]=SIGMA2_INIT_VELZ;

		P[39]=SIGMA2_INIT_ROLL;
		P[42]=SIGMA2_INIT_PITCH;
		P[44]=SIGMA2_INIT_YAW;
		/*************************************************************/

		//Reset the initialization ctr
		initialize_sample_ctr=0;

		//Set the initialization flag
		init_done=true;
	}
}

//@}

/**
	\addtogroup miscel
	@{
*/

/// Routine collecting the functions which need to be run to make a ZUPT update.
void zupt_update(void){
	if(zupt)
	{
		//Calculate the Kalman filter gain
		gain_matrix();

		//Correct the navigation states
		correct_navigation_states();

		//Update the covariance matrix
		measurement_update();
	}
}

//@}

/**
	\defgroup stepwise Step-wise dead reckoning routines
	\brief Routines for running step-wise dead reckoning. These routines serve
	to reset the filter at every step and transmit the displacement and heading
	change estimates together with related covariances to the user. Thereby the
	user can run "step-wise" dead reckoning based on this data.

	@{
*/

void stepwise_system_reset(void) {
	// Local counters
	static precision time_since_last_reset = 0;
	static precision time_reset_pending = 0;

	// No reset per default
	filter_reset_flag=false;
	time_since_last_reset+=dt_k;

	// Check if necessary conditions for reset are met (low covariance AND sufficient time since last reset OR reset pending)
	if ((P[24]<reset_cov_threshold && time_since_last_reset>min_time_between_resets) || time_reset_pending){
		time_reset_pending+=dt_k;
		// Check if sufficient conditions for reset are met (end of OR sufficiently long stationary period)
		if ( !zupt || time_reset_pending>max_time_reset_pending) {
			
			quat2rotation(Rb2t,quat);

			// Save displacement, heading change, and covariance estimates before reset.
			dx[0]=pos[0]; dx[1]=pos[1]; dx[2]=pos[2]; dx[3]=atan2(Rb2t[3],Rb2t[0]);
			dP[0]=P[0]; dP[1]=P[1]; dP[2]= P[2]; dP[3]= P[8];
						dP[4]=P[9]; dP[5]=P[10]; dP[6]=P[16];
									dP[7]=P[17]; dP[8]=P[23];
												 dP[9]=P[44];

			// Reset position, heading, and related covariances
			pos[0]=0; pos[1]=0; pos[2]=0;
			// TODO: Do this directly with the rotation matrix
			vec3 euler;
			rotation2euler(euler,Rb2t);
			euler[2] = 0;
			euler2rotation(Rb2t,euler);
			rotation2quat(quat,Rb2t);

			// Reset covariance estimates
			P[0]=0; P[1]=0; P[2] =0; P[3] =0; P[4] =0; P[5] =0; P[6] =0; P[7] =0; P[8] =0;
			        P[9]=0; P[10]=0; P[11]=0; P[12]=0; P[13]=0; P[14]=0; P[15]=0; P[16]=0;
			                P[17]=0; P[18]=0; P[19]=0; P[20]=0; P[21]=0; P[22]=0; P[23]=0;
			                      /* P[24]    P[25]    P[26]    P[27]    P[28] */ P[29]=0;
			                               /* P[30]    P[31]    P[32]    P[33] */ P[34]=0;
			                                        /* P[35]    P[36]    P[37] */ P[38]=0;
			                                                 /* P[39]    P[40] */ P[41]=0;
			                                                          /* P[42] */ P[43]=0;
			                                                                      P[44]=0;

			// Increase reset counter
			step_counter++;
			// Set reset flag to trigger transmission of above.
			filter_reset_flag = true;
			// Reset local counters
			time_since_last_reset = 0;
			time_reset_pending = 0;
		}
	}

}

//@}

//@}
