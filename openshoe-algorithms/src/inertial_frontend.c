
/** \file
	\authors John-Olof Nilsson
	\copyright Copyright (c) 2014 OpenShoe, Cre�ative Com�mons Attri�bu�tion 4.0 License
*/

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "inertial_frontend.h"
#include "conf_clock.h"

#if defined(MIMU3333)
#  include "conf_mimu3333.h"
#elif defined(MIMU22BT)
#  include "conf_MIMU22BT.h"
#elif defined(MIMU4444)
#  include "conf_MIMU4444.h"
#elif defined(MIMU4444BT)
#  include "conf_MIMU4444BT.h"
#else
#  error No platform specified!
#endif

#ifdef USER_CALIBRATION
#  undef FLOG2_NR_IMU
#  undef ACC_BIAS
#  undef GYRO_BIAS
#  undef CALIBRATION_MATRIX
#  if defined(MIMU3333)
#    include "user_calibration_MIMU3333.h"
#  elif defined(MIMU22BT)
#    include "user_calibration_MIMU22BT.h"
#  elif defined(MIMU4444)
#    include "user_calibration_MIMU4444.h"
#  elif defined(MIMU4444BT)
#    include "user_calibration_MIMU4444BT.h"
#  else
#    error No platform specified!
#  endif
#endif

typedef struct norms {
	uint32_t f_norm, g_norm;
} norms;
typedef struct inert_data {
	inert_int16 u16;
	norms norms;
	inert_int32 u32;
	uint32_t ts;
} inert_data;
typedef struct sum_meas {
	xyz_int32 f, g;
} sum_meas;
typedef struct norms sum_norms;

// IMU (ARRAY) CALIBRATION
static const xyz_int32 f_bias = ACC_BIAS;
static const xyz_int32 g_bias = GYRO_BIAS;
static const int32_t imu_calib[NR_IMU_SLOTS][18] = CALIBRATION_MATRIX;

static inline void comp_calib(inert_int32* restrict u,const int32_t C[18],int16_t imu[6]){
  u->f.x = C[0] *imu[0] + C[1] *imu[1] + C[2] *imu[2];
  u->f.y = C[3] *imu[0] + C[4] *imu[1] + C[5] *imu[2];
  u->f.z = C[6] *imu[0] + C[7] *imu[1] + C[8] *imu[2];
  u->g.x = C[9] *imu[3] + C[10]*imu[4] + C[11]*imu[5];
  u->g.y = C[12]*imu[3] + C[13]*imu[4] + C[14]*imu[5];
  u->g.z = C[15]*imu[3] + C[16]*imu[4] + C[17]*imu[5];
  u->f.x >>= FLOG2_NR_IMU;
  u->f.y >>= FLOG2_NR_IMU;
  u->f.z >>= FLOG2_NR_IMU;
  u->g.x >>= FLOG2_NR_IMU;
  u->g.y >>= FLOG2_NR_IMU;
  u->g.z >>= FLOG2_NR_IMU;
}
static inline void add_meas(inert_int32* sum,inert_int32 meas){
	sum->f.x += meas.f.x;
	sum->f.y += meas.f.y;
	sum->f.z += meas.f.z;
	sum->g.x += meas.g.x;
	sum->g.y += meas.g.y;
	sum->g.z += meas.g.z;
}
static inline void comp_bias(inert_int32* u,xyz_int32 fb,xyz_int32 gb){
	u->f.x -= fb.x;
	u->f.y -= fb.y;
	u->f.z -= fb.z;
	u->g.x -= gb.x;
	u->g.y -= gb.y;
	u->g.z -= gb.z;
}
static inline void sub_sumAB(sum_meas* AB,inert_int16 fg){
	AB->f.x -= fg.f.x;
	AB->f.y -= fg.f.y;
	AB->f.z -= fg.f.z;
	AB->g.x -= fg.g.x;
	AB->g.y -= fg.g.y;
	AB->g.z -= fg.g.z;
}
static inline void sub_sumA_(sum_meas* AB,inert_int16 fg){
	AB->f.x -= fg.f.x;
	AB->f.y -= fg.f.y;
	AB->f.z -= fg.f.z;
}
static inline void sub_norms(sum_norms* CDp,sum_norms* CDpp,norms fg_norms){
	CDp->f_norm -= fg_norms.f_norm & 0xFFFF;
	CDp->g_norm -= fg_norms.g_norm & 0xFFFF;
	CDpp->f_norm -= fg_norms.f_norm >> 16;
	CDpp->g_norm -= fg_norms.g_norm >> 16;
}
static inline void add_sumAB(sum_meas* AB,inert_int16 fg){
	AB->f.x += fg.f.x;
	AB->f.y += fg.f.y;
	AB->f.z += fg.f.z;
	AB->g.x += fg.g.x;
	AB->g.y += fg.g.y;
	AB->g.z += fg.g.z;
}
static inline void add_sumA_(sum_meas* A,inert_int16 fg){
	A->f.x += fg.f.x;
	A->f.y += fg.f.y;
	A->f.z += fg.f.z;
}
static inline void add_norms(sum_norms* CDp,norms* CDpp,norms fg_norms){
	CDp->f_norm += fg_norms.f_norm & 0xFFFF;
	CDp->g_norm += fg_norms.g_norm & 0xFFFF;
	CDpp->f_norm += fg_norms.f_norm >> 16;
	CDpp->g_norm += fg_norms.g_norm >> 16;
}
static inline void upd_buf(inert_data* buf_elem,inert_int32 fg,uint32_t ts){
	buf_elem->u16.f.x = fg.f.x >> BIT_REDU_F;
	buf_elem->u16.f.y = fg.f.y >> BIT_REDU_F;
	buf_elem->u16.f.z = fg.f.z >> BIT_REDU_F;
	buf_elem->norms.f_norm = buf_elem->u16.f.x*buf_elem->u16.f.x + buf_elem->u16.f.y*buf_elem->u16.f.y + buf_elem->u16.f.z*buf_elem->u16.f.z;
	buf_elem->u32.f.x = fg.f.x;
	buf_elem->u32.f.y = fg.f.y;
	buf_elem->u32.f.z = fg.f.z;
	buf_elem->u16.g.x = fg.g.x >> BIT_REDU_G;
	buf_elem->u16.g.y = fg.g.y >> BIT_REDU_G;
	buf_elem->u16.g.z = fg.g.z >> BIT_REDU_G;
	buf_elem->norms.g_norm = buf_elem->u16.g.x*buf_elem->u16.g.x + buf_elem->u16.g.y*buf_elem->u16.g.y + buf_elem->u16.g.z*buf_elem->u16.g.z;
	buf_elem->u32.g.x = fg.g.x;
	buf_elem->u32.g.y = fg.g.y;
	buf_elem->u32.g.z = fg.g.z;
	buf_elem->ts = ts;
}
static inline void ext_buf(inert_int32* fg,uint32_t* t,inert_data buf_elem){
	fg->f.x = buf_elem.u32.f.x;
	fg->f.y = buf_elem.u32.f.y;
	fg->f.z = buf_elem.u32.f.z;
	fg->g.x = buf_elem.u32.g.x;
	fg->g.y = buf_elem.u32.g.y;
	fg->g.z = buf_elem.u32.g.z;
	*t = buf_elem.ts;
}
static inline uint32_t comb_CD(uint32_t CDp,uint32_t CDpp,uint8_t log2){
	return ((CDp & 0xFFFF)>>log2) | ( ( (CDp >> 16) + CDpp ) << (16-log2));
}
static inline uint32_t norm_AB(xyz_int32 fg,uint8_t log2){
	uint32_t Ax12 = abs(fg.x);
	uint32_t Ax34 = Ax12 >> 16;
	Ax12 &= 0xFFFF;
	uint32_t Ay12 = abs(fg.y);
	uint32_t Ay34 = Ay12 >> 16;
	Ay12 &= 0xFFFF;
	uint32_t Az12 = abs(fg.z);
	uint32_t Az34 = Az12 >> 16;
	Az12 &= 0xFFFF;

	uint32_t Ax12_2 = Ax12*Ax12;
	uint32_t Ay12_2 = Ay12*Ay12;
	uint32_t Az12_2 = Az12*Az12;
	uint32_t a = (Ax12_2 & 0xFFFF) + (Ay12_2 & 0xFFFF) + (Az12_2 & 0xFFFF);
	uint32_t b = (Ax12_2 >> 16)    + (Ay12_2 >> 16)    + (Az12_2 >> 16);
	uint32_t c = (Ax12*Ax34 + Ay12*Ay34 + Az12*Az34) << 1;
	uint32_t d =  Ax34*Ax34 + Ay34*Ay34 + Az34*Az34;

	uint32_t tmp = (a >> 16) + b + (c & 0xFFFF);
	uint32_t A2_N2 = ((a & 0xFFFF) >> (2*log2)) | ((tmp & 0xFFFF) << (16-2*log2)) | (((tmp >> 16) + (c >> 16) + d) << (32-2*log2));
	return A2_N2;
}
static inline uint32_t ceillog2(uint32_t x){
	static const uint32_t t[5] = {
			0xFFFF0000u,
			0x0000FF00u,
			0x000000F0u,
			0x0000000Cu,
			0x00000002u
	};

	int y = (((x & (x - 1)) == 0) ? 0 : 1);
	int j = 16;
	int i;

	for (i = 0; i < 5; i++) {
		int k = (((x & t[i]) == 0) ? 0 : j);
		y += k;
		x >>= k;
		j >>= 1;
	}

	return y;
}
static inline void ud_cov_subsample(uint32_t* P){
	static uint8_t Q_cnt = 0;
	Q_cnt++;
	*P += Q_cnt >> NLOG2Q;
	Q_cnt &= (1<<NLOG2Q)-1;
	*P = *P < PMAX ? *P : PMAX;
}
static inline void ud_bias_single(int32_t* b,int32_t u_int,uint32_t K,uint32_t clog2K){
// TODO: Verify under which conditions this <<LOG2BS could cause overflow
// TODO: Verify that this cannot cause a deadlock in the system (too large bias estimate error causes K to be shifted to zero all the time)
	int32_t diffx = (u_int<<LOG2BS)-*b;
	uint32_t clog2diffx = ceillog2(abs(diffx));
	int32_t shiftx = clog2K+clog2diffx-31;
	shiftx = shiftx & (~(shiftx>>31));
	*b += ((int32_t)( (K>>shiftx) * diffx ))>>(LOG2GS - shiftx);
}
static inline void comp_est_bias(inert_int32* u_int,xyz_int32 b){
	u_int->g.x -= b.x>>LOG2BS;
	u_int->g.y -= b.y>>LOG2BS;
	u_int->g.z -= b.z>>LOG2BS;
}
static inline void conv_to_float(inert_float* u_float,inert_int32 u_int){
	u_float->f.x = ((precision)u_int.f.x)*FSCALE;
	u_float->f.y = ((precision)u_int.f.y)*FSCALE;
	u_float->f.z = ((precision)u_int.f.z)*FSCALE;
	u_float->g.x = ((precision)u_int.g.x)*GSCALE;
	u_float->g.y = ((precision)u_int.g.y)*GSCALE;
	u_float->g.z = ((precision)u_int.g.z)*GSCALE;
}
static inline void float_to_int16(inert_int16* u_int16,inert_float u_float){
	u_int16->f.x = ((int32_t)(u_float.f.x*13002079.3645f))>>16;
	u_int16->f.y = ((int32_t)(u_float.f.y*13002079.3645f))>>16;
	u_int16->f.z = ((int32_t)(u_float.f.z*13002079.3645f))>>16;
	u_int16->g.x = ((int32_t)(u_float.g.x*58444831.0618f))>>16;
	u_int16->g.y = ((int32_t)(u_float.g.y*58444831.0618f))>>16;
	u_int16->g.z = ((int32_t)(u_float.g.z*58444831.0618f))>>16;
}

extern int16_t mimu_data[32][10];
extern uint32_t ts_u;
inert_int32 u_new;
inert_int32 u_int_k;
inert_int16 u_int16_k;
inert_int32 u_bias;
uint32_t t_int_k;
uint32_t t_bias;
inert_float u_k;
precision dt_k;

static inert_data buf[N];
static int inza = N-1, outza = 0, kza = N/2;
static int inzu = (N/2+M/2-1+ZAZU_OFFSET)&(N-1);
static int outzu = (N/2-M/2+ZAZU_OFFSET)&(N-1);
static int kzu = N/2+ZAZU_OFFSET;

uint32_t T1s2f;
uint32_t T2s2f;
uint32_t th_zupt = TH_ZUPT;
uint32_t th_zaru = TH_ZARU;
bool zupt;
bool zaru;

static xyz_float gb = {0,0,0};

// Inefficient logarithm algorithm, intended for compile-time constants.
#define FLOOR_LOG2_8BIT(v)  (8 - 90/(((v)/4+14)|1) - 2/((v)/2+1))
#define FLOOR_LOG2_16BIT(v) (8*((v)>255) + FLOOR_LOG2_8BIT((v) >>8*((v)>255)))
#define FLOOR_LOG2_32BIT(v) (16*((v)>65535L) + FLOOR_LOG2_16BIT((v)*1L >>16*((v)>65535L)))
uint32_t sqrt_trk_T1(uint32_t x) {
	uint32_t a, m;
	static uint32_t b = 0;
	static uint32_t y = 0;
	if (y < x){
		a = b + 1;
		b = b + SQRT_TRK_C;
	} else {
		a = b - SQRT_TRK_C;
		a = ~((int32_t)a) >> 31 & a;
		a++;
	}
	y = x;
	for(int i=0; i<FLOOR_LOG2_32BIT(SQRT_TRK_C+1)+1; i++) { // TODO: Should be replaced with MREPEAT macro to ensure unrolling
		m = (a + b) >> 1;
		(m*m > x) ? (b = m - 1) : (a = m + 1);
	}
	b = a - 1;
	return b;
}

void frontend_preproc(void){
	struct inert_int32 u_tmp = {{0}};

	u_new.f.x = u_new.f.y = u_new.f.z = 0;
	u_new.g.x = u_new.g.y = u_new.g.z = 0;

	uint8_t i;
	for(i = 0; i < NR_IMUS; i++){
		comp_calib(&u_tmp,imu_calib[i],mimu_data[i]);
		add_meas(&u_new,u_tmp);
	}
	comp_bias(&u_new,f_bias,g_bias);
	conv_to_float(&u_k,u_new);
}
void frontend_statdet(void){
	static sum_meas ABza, ABzu;
	static sum_norms CDpza , CDppza, CDpzu , CDppzu;

	sub_sumAB(&ABza,buf[outza].u16);
	sub_norms(&CDpza,&CDppza,buf[outza].norms);
	sub_sumA_(&ABzu,buf[outzu].u16);
	sub_norms(&CDpzu,&CDppzu,buf[outzu].norms);

	++inza; inza &= N-1; ++outza; outza &= N-1; ++kza; kza &= N-1;
	++inzu; inzu &= N-1; ++outzu; outzu &= N-1;	++kzu; kzu &= N-1;
	upd_buf(&buf[inza],u_new,ts_u);
	ext_buf(&u_bias,&t_bias,buf[kza]);
	ext_buf(&u_int_k,&t_int_k,buf[kzu]);

	add_sumAB(&ABza,buf[inza].u16);
	add_norms(&CDpza,&CDppza,buf[inza].norms);
	add_sumA_(&ABzu,buf[inzu].u16);
	add_norms(&CDpzu,&CDppzu,buf[inzu].norms);

	uint32_t Cza_N = comb_CD(CDpza.f_norm,CDppza.f_norm,LOG2N);
	uint32_t Dza_N = comb_CD(CDpza.g_norm,CDppza.g_norm,LOG2N);
	uint32_t Aza2_N2 = norm_AB(ABza.f,LOG2N);
	uint32_t Bza2_N2 = norm_AB(ABza.g,LOG2N);
	T2s2f = Cza_N - Aza2_N2 + SR*(Dza_N - Bza2_N2);

	uint32_t Czu_N = comb_CD(CDpzu.f_norm,CDppzu.f_norm,LOG2M);
	uint32_t Dzu_N = comb_CD(CDpzu.g_norm,CDppzu.g_norm,LOG2M);
	uint32_t Azu2_N2 = norm_AB(ABzu.f,LOG2M);
	T1s2f = Czu_N + G2 - TWOG*sqrt_trk_T1(Azu2_N2) + SR*Dzu_N;

	zupt = T1s2f<th_zupt;
	zaru = zupt && T2s2f<th_zaru;
}
void frontend_biasest(void){
//	static uint32_t t_int_km1 = 0;
/*
	// Integer bias estimation
	static uint32_t P = PMAX;
	static xyz_int32 b = {0,0,0};
	

	ud_cov_subsample(&P);

	if(zaru){
		uint32_t K = (P<<LOG2GS)/(P+R_BIAS);
		uint32_t clog2K = ceillog2(K);
		ud_bias_single(&b.x,u_int_k.g.x,K,clog2K);
		ud_bias_single(&b.y,u_int_k.g.y,K,clog2K);
		ud_bias_single(&b.z,u_int_k.g.z,K,clog2K);

		P = ((1<<LOG2GS) - K)*P >> LOG2GS;
	}
	comp_est_bias(&u_int_k,b);
	conv_to_float(&u_k,u_int_k);
*/
	// Floating point bias estimation
	static precision P = 0.02f;
	P += QF;
	if(zaru){
		inert_float u_fbias;
		conv_to_float(&u_fbias,u_bias);
		precision K = P/(P+RF);
		gb.x=gb.x+K*(u_fbias.g.x-gb.x);
		gb.y=gb.y+K*(u_fbias.g.y-gb.y);
		gb.z=gb.z+K*(u_fbias.g.z-gb.z);
		P*=(1.0f-K);
	}
}
void frontend_convcomp(void){
	static uint32_t t_int_km1 = 0;
	
	conv_to_float(&u_k,u_int_k);
	u_k.g.x-=gb.x;
	u_k.g.y-=gb.y;
	u_k.g.z-=gb.z;
	
	dt_k = ((precision)(t_int_k - t_int_km1))*TIME_SCALE;
	t_int_km1 = t_int_k;
}
void frontend_conv16(void){
	float_to_int16(&u_int16_k,u_k);
}

//inert_float frontend_lp(uint redu){
	// I would like to have argument in [s] but
	// window size = 2^redu
	// loop over window around k (split in 16+16)
	// Divide and recombine
	// Convert and bias compensate result
	// return result
//}
