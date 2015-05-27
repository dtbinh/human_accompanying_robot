#ifndef __c16_IMM_UKF_h__
#define __c16_IMM_UKF_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc16_IMM_UKFInstanceStruct
#define typedef_SFc16_IMM_UKFInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c16_sfEvent;
  boolean_T c16_isStable;
  boolean_T c16_doneDoubleBufferReInit;
  uint8_T c16_is_active_c16_IMM_UKF;
} SFc16_IMM_UKFInstanceStruct;

#endif                                 /*typedef_SFc16_IMM_UKFInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c16_IMM_UKF_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c16_IMM_UKF_get_check_sum(mxArray *plhs[]);
extern void c16_IMM_UKF_method_dispatcher(SimStruct *S, int_T method, void *data);

#endif
