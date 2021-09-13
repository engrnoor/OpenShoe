
/** \file
	\authors John-Olof Nilsson
	\copyright Copyright (c) 2014 OpenShoe, Cre­ative Com­mons Attri­bu­tion 4.0 License
*/

#ifndef INERTIAL_FRONTEND_H_
#define INERTIAL_FRONTEND_H_

#include "nav_types.h"
#include <stdint.h>
#include <stdbool.h>

typedef struct xyz_int32 {
	int32_t x, y, z;
} xyz_int32;
typedef struct xyz_int16 {
	int16_t x, y, z;
} xyz_int16;
typedef struct xyz_float {
	precision x, y, z;
} xyz_float;
typedef struct inert_int16 {
	xyz_int16 f, g;
} inert_int16;
typedef struct inert_int32 {
	xyz_int32 f, g;
} inert_int32;
typedef struct inert_float {
	xyz_float f, g;
} inert_float;

extern inert_int32 u_new;
extern inert_int32 u_int_k;
extern inert_int16 u_int16_k;
extern uint32_t t_int_k;
extern inert_float u_k;
extern precision dt_k;

extern uint32_t T1s2f;
extern uint32_t T2s2f;
extern uint32_t th_zupt;
extern uint32_t th_zaru;
extern bool zupt;
extern bool zaru;

void frontend_preproc(void);
void frontend_statdet(void);
void frontend_biasest(void);
void frontend_convcomp(void);
void frontend_conv16(void);

#endif /* INERTIAL_FRONTEND_H_ */
