#ifndef __c2_UKF_h__
#define __c2_UKF_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc2_UKFInstanceStruct
#define typedef_SFc2_UKFInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c2_sfEvent;
  boolean_T c2_isStable;
  boolean_T c2_doneDoubleBufferReInit;
  uint8_T c2_is_active_c2_UKF;
  real_T c2_P[25];
  boolean_T c2_P_not_empty;
  real_T c2_x[55];
  boolean_T c2_x_not_empty;
  real_T c2_W1[25];
  boolean_T c2_W1_not_empty;
  real_T c2_V1[4];
  boolean_T c2_V1_not_empty;
} SFc2_UKFInstanceStruct;

#endif                                 /*typedef_SFc2_UKFInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c2_UKF_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c2_UKF_get_check_sum(mxArray *plhs[]);
extern void c2_UKF_method_dispatcher(SimStruct *S, int_T method, void *data);

#endif
