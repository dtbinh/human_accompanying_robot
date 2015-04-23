#ifndef __c17_IMM_UKF_h__
#define __c17_IMM_UKF_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc17_IMM_UKFInstanceStruct
#define typedef_SFc17_IMM_UKFInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c17_sfEvent;
  boolean_T c17_isStable;
  boolean_T c17_doneDoubleBufferReInit;
  uint8_T c17_is_active_c17_IMM_UKF;
} SFc17_IMM_UKFInstanceStruct;

#endif                                 /*typedef_SFc17_IMM_UKFInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c17_IMM_UKF_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c17_IMM_UKF_get_check_sum(mxArray *plhs[]);
extern void c17_IMM_UKF_method_dispatcher(SimStruct *S, int_T method, void *data);

#endif
