/* Include files */

#include <stddef.h>
#include "blas.h"
#include "IMM_UKF_sfun.h"
#include "c13_IMM_UKF.h"
#include <math.h>
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "IMM_UKF_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c13_debug_family_names[44] = { "T", "L", "alpha", "kappa",
  "beta", "lambda", "c", "Wm", "Wc", "H", "meas_noise1", "a_process_noise",
  "w_process_noise", "Sv1", "Sw1", "Sw2", "Sw", "Bw", "Dw", "x", "sP", "chi_p",
  "chi_m", "i", "x_m", "P_m", "psi_m", "y_m", "Pyy", "Pxy", "K", "innov",
  "c_out", "Wm_out", "nargin", "nargout", "xhat20", "meas", "xhat2", "Lambda2",
  "P_out", "P2", "W2", "V2" };

/* Function Declarations */
static void initialize_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void initialize_params_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct
  *chartInstance);
static void enable_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void disable_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_update_debugger_state_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct
  *chartInstance);
static void set_sim_state_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_st);
static void finalize_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void sf_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_chartstep_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void initSimStructsc13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c13_machineNumber, uint32_T
  c13_chartNumber);
static const mxArray *c13_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static void c13_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_b_V2, const char_T *c13_identifier, real_T c13_y[4]);
static void c13_b_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[4]);
static void c13_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static const mxArray *c13_b_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static void c13_c_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_b_W2, const char_T *c13_identifier, real_T c13_y[25]);
static void c13_d_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[25]);
static void c13_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static const mxArray *c13_c_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static void c13_e_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_b_P2, const char_T *c13_identifier, real_T c13_y[25]);
static void c13_f_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[25]);
static void c13_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static const mxArray *c13_d_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static void c13_g_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_P_out, const char_T *c13_identifier, real_T c13_y[25]);
static void c13_h_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[25]);
static void c13_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static const mxArray *c13_e_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static real_T c13_i_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_Lambda2, const char_T *c13_identifier);
static real_T c13_j_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId);
static void c13_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static const mxArray *c13_f_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static void c13_k_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_xhat2, const char_T *c13_identifier, real_T c13_y[5]);
static void c13_l_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[5]);
static void c13_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static const mxArray *c13_g_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static const mxArray *c13_h_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static void c13_m_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[11]);
static void c13_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static void c13_n_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[2]);
static void c13_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static const mxArray *c13_i_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static void c13_o_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[10]);
static void c13_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static const mxArray *c13_j_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static void c13_p_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[4]);
static void c13_j_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static const mxArray *c13_k_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static void c13_q_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[22]);
static void c13_k_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static const mxArray *c13_l_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static void c13_r_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[55]);
static void c13_l_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static const mxArray *c13_m_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static const mxArray *c13_n_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static const mxArray *c13_o_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static void c13_info_helper(const mxArray **c13_info);
static const mxArray *c13_emlrt_marshallOut(char * c13_u);
static const mxArray *c13_b_emlrt_marshallOut(uint32_T c13_u);
static void c13_b_info_helper(const mxArray **c13_info);
static void c13_c_info_helper(const mxArray **c13_info);
static real_T c13_mpower(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T
  c13_a);
static void c13_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance);
static real_T c13_sqrt(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T c13_x);
static void c13_eml_error(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_b_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_c_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_d_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_eye(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T c13_I[25]);
static void c13_b_eml_error(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_eml_matlab_zpotrf(SFc13_IMM_UKFInstanceStruct *chartInstance,
  real_T c13_A[25], real_T c13_b_A[25], int32_T *c13_info);
static void c13_check_forloop_overflow_error(SFc13_IMM_UKFInstanceStruct
  *chartInstance, boolean_T c13_overflow);
static void c13_eml_xgemv(SFc13_IMM_UKFInstanceStruct *chartInstance, int32_T
  c13_m, int32_T c13_n, int32_T c13_ia0, int32_T c13_ix0, real_T c13_y[25],
  int32_T c13_iy0, real_T c13_b_y[25]);
static void c13_c_eml_error(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_e_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_f_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_g_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_inv(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T c13_x[4],
                    real_T c13_y[4]);
static void c13_inv2x2(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T c13_x
  [4], real_T c13_y[4]);
static real_T c13_norm(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T c13_x
  [4]);
static void c13_eml_warning(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_b_eml_warning(SFc13_IMM_UKFInstanceStruct *chartInstance, char_T
  c13_varargin_2[14]);
static void c13_h_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_i_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_j_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_k_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance);
static void c13_l_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance);
static real_T c13_expm(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T c13_A);
static real_T c13_PadeApproximantOfDegree(SFc13_IMM_UKFInstanceStruct
  *chartInstance, real_T c13_A, real_T c13_m);
static void c13_s_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_sprintf, const char_T *c13_identifier, char_T c13_y[14]);
static void c13_t_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, char_T c13_y[14]);
static const mxArray *c13_p_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static int32_T c13_u_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId);
static void c13_m_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static uint8_T c13_v_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_b_is_active_c13_IMM_UKF, const char_T *c13_identifier);
static uint8_T c13_w_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId);
static void c13_b_sqrt(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T *c13_x);
static int32_T c13_b_eml_matlab_zpotrf(SFc13_IMM_UKFInstanceStruct
  *chartInstance, real_T c13_A[25]);
static void c13_b_eml_xgemv(SFc13_IMM_UKFInstanceStruct *chartInstance, int32_T
  c13_m, int32_T c13_n, int32_T c13_ia0, int32_T c13_ix0, real_T c13_y[25],
  int32_T c13_iy0);
static void init_dsm_address_info(SFc13_IMM_UKFInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
  chartInstance->c13_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c13_P2_not_empty = FALSE;
  chartInstance->c13_W2_not_empty = FALSE;
  chartInstance->c13_V2_not_empty = FALSE;
  chartInstance->c13_is_active_c13_IMM_UKF = 0U;
}

static void initialize_params_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct
  *chartInstance)
{
}

static void enable_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c13_update_debugger_state_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct
  *chartInstance)
{
}

static const mxArray *get_sim_state_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct
  *chartInstance)
{
  const mxArray *c13_st;
  const mxArray *c13_y = NULL;
  real_T c13_hoistedGlobal;
  real_T c13_u;
  const mxArray *c13_b_y = NULL;
  int32_T c13_i0;
  real_T c13_b_u[25];
  const mxArray *c13_c_y = NULL;
  int32_T c13_i1;
  real_T c13_c_u[5];
  const mxArray *c13_d_y = NULL;
  int32_T c13_i2;
  real_T c13_d_u[25];
  const mxArray *c13_e_y = NULL;
  int32_T c13_i3;
  real_T c13_e_u[4];
  const mxArray *c13_f_y = NULL;
  int32_T c13_i4;
  real_T c13_f_u[25];
  const mxArray *c13_g_y = NULL;
  uint8_T c13_b_hoistedGlobal;
  uint8_T c13_g_u;
  const mxArray *c13_h_y = NULL;
  real_T *c13_Lambda2;
  real_T (*c13_xhat2)[5];
  real_T (*c13_P_out)[25];
  c13_P_out = (real_T (*)[25])ssGetOutputPortSignal(chartInstance->S, 3);
  c13_Lambda2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c13_xhat2 = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
  c13_st = NULL;
  c13_st = NULL;
  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_createcellarray(7), FALSE);
  c13_hoistedGlobal = *c13_Lambda2;
  c13_u = c13_hoistedGlobal;
  c13_b_y = NULL;
  sf_mex_assign(&c13_b_y, sf_mex_create("y", &c13_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c13_y, 0, c13_b_y);
  for (c13_i0 = 0; c13_i0 < 25; c13_i0++) {
    c13_b_u[c13_i0] = (*c13_P_out)[c13_i0];
  }

  c13_c_y = NULL;
  sf_mex_assign(&c13_c_y, sf_mex_create("y", c13_b_u, 0, 0U, 1U, 0U, 2, 5, 5),
                FALSE);
  sf_mex_setcell(c13_y, 1, c13_c_y);
  for (c13_i1 = 0; c13_i1 < 5; c13_i1++) {
    c13_c_u[c13_i1] = (*c13_xhat2)[c13_i1];
  }

  c13_d_y = NULL;
  sf_mex_assign(&c13_d_y, sf_mex_create("y", c13_c_u, 0, 0U, 1U, 0U, 1, 5),
                FALSE);
  sf_mex_setcell(c13_y, 2, c13_d_y);
  for (c13_i2 = 0; c13_i2 < 25; c13_i2++) {
    c13_d_u[c13_i2] = chartInstance->c13_P2[c13_i2];
  }

  c13_e_y = NULL;
  if (!chartInstance->c13_P2_not_empty) {
    sf_mex_assign(&c13_e_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c13_e_y, sf_mex_create("y", c13_d_u, 0, 0U, 1U, 0U, 2, 5, 5),
                  FALSE);
  }

  sf_mex_setcell(c13_y, 3, c13_e_y);
  for (c13_i3 = 0; c13_i3 < 4; c13_i3++) {
    c13_e_u[c13_i3] = chartInstance->c13_V2[c13_i3];
  }

  c13_f_y = NULL;
  if (!chartInstance->c13_V2_not_empty) {
    sf_mex_assign(&c13_f_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c13_f_y, sf_mex_create("y", c13_e_u, 0, 0U, 1U, 0U, 2, 2, 2),
                  FALSE);
  }

  sf_mex_setcell(c13_y, 4, c13_f_y);
  for (c13_i4 = 0; c13_i4 < 25; c13_i4++) {
    c13_f_u[c13_i4] = chartInstance->c13_W2[c13_i4];
  }

  c13_g_y = NULL;
  if (!chartInstance->c13_W2_not_empty) {
    sf_mex_assign(&c13_g_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c13_g_y, sf_mex_create("y", c13_f_u, 0, 0U, 1U, 0U, 2, 5, 5),
                  FALSE);
  }

  sf_mex_setcell(c13_y, 5, c13_g_y);
  c13_b_hoistedGlobal = chartInstance->c13_is_active_c13_IMM_UKF;
  c13_g_u = c13_b_hoistedGlobal;
  c13_h_y = NULL;
  sf_mex_assign(&c13_h_y, sf_mex_create("y", &c13_g_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c13_y, 6, c13_h_y);
  sf_mex_assign(&c13_st, c13_y, FALSE);
  return c13_st;
}

static void set_sim_state_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_st)
{
  const mxArray *c13_u;
  real_T c13_dv0[25];
  int32_T c13_i5;
  real_T c13_dv1[5];
  int32_T c13_i6;
  real_T c13_dv2[25];
  int32_T c13_i7;
  real_T c13_dv3[4];
  int32_T c13_i8;
  real_T c13_dv4[25];
  int32_T c13_i9;
  real_T *c13_Lambda2;
  real_T (*c13_P_out)[25];
  real_T (*c13_xhat2)[5];
  c13_P_out = (real_T (*)[25])ssGetOutputPortSignal(chartInstance->S, 3);
  c13_Lambda2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c13_xhat2 = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c13_doneDoubleBufferReInit = TRUE;
  c13_u = sf_mex_dup(c13_st);
  *c13_Lambda2 = c13_i_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c13_u, 0)), "Lambda2");
  c13_g_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c13_u, 1)),
    "P_out", c13_dv0);
  for (c13_i5 = 0; c13_i5 < 25; c13_i5++) {
    (*c13_P_out)[c13_i5] = c13_dv0[c13_i5];
  }

  c13_k_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c13_u, 2)),
    "xhat2", c13_dv1);
  for (c13_i6 = 0; c13_i6 < 5; c13_i6++) {
    (*c13_xhat2)[c13_i6] = c13_dv1[c13_i6];
  }

  c13_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c13_u, 3)),
    "P2", c13_dv2);
  for (c13_i7 = 0; c13_i7 < 25; c13_i7++) {
    chartInstance->c13_P2[c13_i7] = c13_dv2[c13_i7];
  }

  c13_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c13_u, 4)), "V2",
                       c13_dv3);
  for (c13_i8 = 0; c13_i8 < 4; c13_i8++) {
    chartInstance->c13_V2[c13_i8] = c13_dv3[c13_i8];
  }

  c13_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c13_u, 5)),
    "W2", c13_dv4);
  for (c13_i9 = 0; c13_i9 < 25; c13_i9++) {
    chartInstance->c13_W2[c13_i9] = c13_dv4[c13_i9];
  }

  chartInstance->c13_is_active_c13_IMM_UKF = c13_v_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c13_u, 6)),
     "is_active_c13_IMM_UKF");
  sf_mex_destroy(&c13_u);
  c13_update_debugger_state_c13_IMM_UKF(chartInstance);
  sf_mex_destroy(&c13_st);
}

static void finalize_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static void sf_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
  int32_T c13_i10;
  int32_T c13_i11;
  int32_T c13_i12;
  int32_T c13_i13;
  real_T *c13_Lambda2;
  real_T (*c13_P_out)[25];
  real_T (*c13_xhat2)[5];
  real_T (*c13_meas)[2];
  real_T (*c13_xhat20)[5];
  c13_P_out = (real_T (*)[25])ssGetOutputPortSignal(chartInstance->S, 3);
  c13_Lambda2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c13_xhat2 = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
  c13_meas = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 1);
  c13_xhat20 = (real_T (*)[5])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 6U, chartInstance->c13_sfEvent);
  for (c13_i10 = 0; c13_i10 < 5; c13_i10++) {
    _SFD_DATA_RANGE_CHECK((*c13_xhat20)[c13_i10], 0U);
  }

  for (c13_i11 = 0; c13_i11 < 2; c13_i11++) {
    _SFD_DATA_RANGE_CHECK((*c13_meas)[c13_i11], 1U);
  }

  for (c13_i12 = 0; c13_i12 < 5; c13_i12++) {
    _SFD_DATA_RANGE_CHECK((*c13_xhat2)[c13_i12], 2U);
  }

  _SFD_DATA_RANGE_CHECK(*c13_Lambda2, 3U);
  for (c13_i13 = 0; c13_i13 < 25; c13_i13++) {
    _SFD_DATA_RANGE_CHECK((*c13_P_out)[c13_i13], 4U);
  }

  chartInstance->c13_sfEvent = CALL_EVENT;
  c13_chartstep_c13_IMM_UKF(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_IMM_UKFMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c13_chartstep_c13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
  int32_T c13_i14;
  real_T c13_xhat20[5];
  int32_T c13_i15;
  real_T c13_meas[2];
  uint32_T c13_debug_family_var_map[44];
  real_T c13_T;
  real_T c13_L;
  real_T c13_alpha;
  real_T c13_kappa;
  real_T c13_beta;
  real_T c13_lambda;
  real_T c13_c;
  real_T c13_Wm[11];
  real_T c13_Wc[11];
  real_T c13_H[10];
  real_T c13_meas_noise1;
  real_T c13_a_process_noise;
  real_T c13_w_process_noise;
  real_T c13_Sv1;
  real_T c13_Sw1;
  real_T c13_Sw2;
  real_T c13_Sw[9];
  real_T c13_Bw[15];
  real_T c13_Dw[4];
  real_T c13_x[5];
  real_T c13_sP[25];
  real_T c13_chi_p[55];
  real_T c13_chi_m[55];
  real_T c13_i;
  real_T c13_x_m[5];
  real_T c13_P_m[25];
  real_T c13_psi_m[22];
  real_T c13_y_m[2];
  real_T c13_Pyy[4];
  real_T c13_Pxy[10];
  real_T c13_K[10];
  real_T c13_innov[2];
  real_T c13_c_out;
  real_T c13_Wm_out[11];
  real_T c13_nargin = 2.0;
  real_T c13_nargout = 3.0;
  real_T c13_xhat2[5];
  real_T c13_Lambda2;
  real_T c13_P_out[25];
  int32_T c13_i16;
  int32_T c13_i17;
  int32_T c13_i18;
  static real_T c13_a[10] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };

  int32_T c13_i19;
  static real_T c13_dv5[9] = { 0.000225, 0.0, 0.0, 0.0, 0.000225, 0.0, 0.0, 0.0,
    2.25E-6 };

  int32_T c13_i20;
  static real_T c13_dv6[15] = { 0.05, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  int32_T c13_i21;
  static real_T c13_b[4] = { 1.0, 0.0, 0.0, 1.0 };

  int32_T c13_i22;
  static real_T c13_y[25] = { 5.625E-7, 1.125E-5, 0.0, 0.0, 0.0, 1.125E-5,
    0.000225, 0.0, 0.0, 0.0, 0.0, 0.0, 5.625E-7, 1.125E-5, 0.0, 0.0, 0.0,
    1.125E-5, 0.000225, 0.0, 0.0, 0.0, 0.0, 0.0, 2.25E-6 };

  int32_T c13_i23;
  static real_T c13_b_y[4] = { 2.25, 0.0, 0.0, 2.25 };

  real_T c13_dv7[25];
  int32_T c13_i24;
  int32_T c13_i25;
  int32_T c13_i26;
  real_T c13_hoistedGlobal[25];
  int32_T c13_j;
  int32_T c13_b_j;
  int32_T c13_info;
  int32_T c13_b_info;
  int32_T c13_c_info;
  int32_T c13_d_info;
  int32_T c13_jmax;
  int32_T c13_b_a;
  int32_T c13_b_jmax;
  int32_T c13_b_b;
  int32_T c13_c_b;
  boolean_T c13_overflow;
  int32_T c13_c_j;
  int32_T c13_c_a;
  int32_T c13_i27;
  int32_T c13_d_b;
  int32_T c13_e_b;
  boolean_T c13_b_overflow;
  int32_T c13_b_i;
  int32_T c13_c_i;
  int32_T c13_i28;
  int32_T c13_i29;
  real_T c13_d_a[5];
  int32_T c13_i30;
  int32_T c13_i31;
  int32_T c13_i32;
  int32_T c13_i33;
  real_T c13_f_b[25];
  int32_T c13_i34;
  int32_T c13_i35;
  int32_T c13_i36;
  int32_T c13_i37;
  int32_T c13_i38;
  real_T c13_c_y[25];
  int32_T c13_i39;
  real_T c13_g_b[25];
  int32_T c13_i40;
  int32_T c13_i41;
  int32_T c13_i42;
  int32_T c13_i43;
  int32_T c13_i44;
  int32_T c13_i45;
  int32_T c13_i46;
  int32_T c13_i47;
  int32_T c13_i48;
  int32_T c13_d_i;
  real_T c13_e_a;
  real_T c13_d_y;
  real_T c13_b_x;
  real_T c13_c_x;
  real_T c13_f_a;
  real_T c13_h_b;
  real_T c13_e_y;
  real_T c13_A;
  real_T c13_B;
  real_T c13_d_x;
  real_T c13_f_y;
  real_T c13_e_x;
  real_T c13_g_y;
  real_T c13_h_y;
  real_T c13_g_a;
  real_T c13_i_y;
  real_T c13_f_x;
  real_T c13_g_x;
  real_T c13_h_a;
  real_T c13_i_b;
  real_T c13_j_y;
  real_T c13_b_A;
  real_T c13_b_B;
  real_T c13_h_x;
  real_T c13_k_y;
  real_T c13_i_x;
  real_T c13_l_y;
  real_T c13_m_y;
  real_T c13_i_a;
  real_T c13_n_y;
  real_T c13_j_x;
  real_T c13_k_x;
  real_T c13_j_a;
  real_T c13_j_b;
  real_T c13_o_y;
  real_T c13_k_a;
  real_T c13_p_y;
  real_T c13_l_x;
  real_T c13_m_x;
  real_T c13_l_a;
  real_T c13_k_b;
  real_T c13_q_y;
  real_T c13_m_a;
  real_T c13_r_y;
  real_T c13_n_x;
  real_T c13_o_x;
  real_T c13_n_a;
  real_T c13_l_b;
  real_T c13_s_y;
  real_T c13_c_A;
  real_T c13_c_B;
  real_T c13_p_x;
  real_T c13_t_y;
  real_T c13_q_x;
  real_T c13_u_y;
  real_T c13_v_y;
  real_T c13_o_a;
  real_T c13_w_y;
  real_T c13_r_x;
  real_T c13_s_x;
  real_T c13_p_a;
  real_T c13_m_b;
  real_T c13_x_y;
  real_T c13_d_A;
  real_T c13_d_B;
  real_T c13_t_x;
  real_T c13_y_y;
  real_T c13_u_x;
  real_T c13_ab_y;
  real_T c13_bb_y;
  real_T c13_q_a;
  real_T c13_cb_y;
  real_T c13_v_x;
  real_T c13_w_x;
  real_T c13_r_a;
  real_T c13_n_b;
  real_T c13_db_y;
  real_T c13_s_a;
  real_T c13_eb_y;
  real_T c13_x_x;
  real_T c13_y_x;
  real_T c13_t_a;
  real_T c13_o_b;
  real_T c13_fb_y;
  int32_T c13_i49;
  real_T c13_u_a[55];
  int32_T c13_i50;
  real_T c13_p_b[11];
  int32_T c13_i51;
  int32_T c13_i52;
  int32_T c13_i53;
  int32_T c13_i54;
  int32_T c13_i55;
  int32_T c13_i56;
  int32_T c13_i57;
  int32_T c13_i58;
  int32_T c13_i59;
  int32_T c13_i60;
  int32_T c13_e_i;
  real_T c13_v_a;
  int32_T c13_f_i;
  int32_T c13_i61;
  int32_T c13_i62;
  int32_T c13_g_i;
  int32_T c13_i63;
  real_T c13_q_b[5];
  int32_T c13_i64;
  int32_T c13_i65;
  int32_T c13_i66;
  int32_T c13_i67;
  int32_T c13_i68;
  int32_T c13_i69;
  int32_T c13_i70;
  int32_T c13_i71;
  real_T c13_C[22];
  int32_T c13_i72;
  int32_T c13_i73;
  int32_T c13_i74;
  int32_T c13_i75;
  int32_T c13_i76;
  int32_T c13_i77;
  int32_T c13_i78;
  int32_T c13_i79;
  int32_T c13_i80;
  int32_T c13_i81;
  int32_T c13_i82;
  int32_T c13_i83;
  int32_T c13_i84;
  int32_T c13_i85;
  real_T c13_b_C[2];
  int32_T c13_i86;
  int32_T c13_i87;
  int32_T c13_i88;
  int32_T c13_i89;
  int32_T c13_i90;
  int32_T c13_i91;
  int32_T c13_i92;
  int32_T c13_i93;
  int32_T c13_h_i;
  real_T c13_w_a;
  int32_T c13_i_i;
  int32_T c13_i94;
  int32_T c13_i95;
  int32_T c13_j_i;
  int32_T c13_i96;
  real_T c13_x_a[2];
  int32_T c13_i97;
  int32_T c13_i98;
  int32_T c13_i99;
  real_T c13_r_b[4];
  int32_T c13_i100;
  real_T c13_y_a;
  int32_T c13_k_i;
  int32_T c13_i101;
  int32_T c13_i102;
  int32_T c13_l_i;
  int32_T c13_i103;
  int32_T c13_i104;
  int32_T c13_i105;
  int32_T c13_i106;
  real_T c13_c_C[10];
  int32_T c13_i107;
  int32_T c13_i108;
  real_T c13_ab_a[10];
  int32_T c13_i109;
  real_T c13_b_Pyy[4];
  int32_T c13_i110;
  int32_T c13_i111;
  int32_T c13_i112;
  int32_T c13_i113;
  int32_T c13_i114;
  int32_T c13_i115;
  int32_T c13_i116;
  int32_T c13_i117;
  int32_T c13_i118;
  int32_T c13_i119;
  int32_T c13_i120;
  int32_T c13_i121;
  int32_T c13_i122;
  int32_T c13_i123;
  int32_T c13_i124;
  int32_T c13_i125;
  int32_T c13_i126;
  int32_T c13_i127;
  int32_T c13_i128;
  int32_T c13_i129;
  int32_T c13_i130;
  int32_T c13_i131;
  int32_T c13_i132;
  int32_T c13_i133;
  int32_T c13_i134;
  int32_T c13_i135;
  int32_T c13_i136;
  int32_T c13_i137;
  int32_T c13_i138;
  int32_T c13_i139;
  real_T c13_s_b[10];
  int32_T c13_i140;
  int32_T c13_i141;
  int32_T c13_i142;
  int32_T c13_i143;
  int32_T c13_i144;
  int32_T c13_i145;
  int32_T c13_i146;
  int32_T c13_i147;
  real_T c13_t_b;
  real_T c13_gb_y;
  real_T c13_bb_a;
  real_T c13_u_b;
  real_T c13_hb_y;
  real_T c13_e_B;
  real_T c13_ib_y;
  real_T c13_jb_y;
  real_T c13_kb_y;
  int32_T c13_i148;
  int32_T c13_i149;
  real_T c13_c_Pyy[4];
  int32_T c13_i150;
  int32_T c13_i151;
  real_T c13_lb_y[2];
  int32_T c13_i152;
  int32_T c13_i153;
  real_T c13_mb_y;
  int32_T c13_k;
  int32_T c13_b_k;
  real_T c13_v_b;
  real_T c13_nb_y;
  real_T c13_cb_a;
  real_T c13_w_b;
  int32_T c13_i154;
  int32_T c13_i155;
  int32_T c13_i156;
  int32_T c13_i157;
  real_T (*c13_b_P_out)[25];
  real_T (*c13_b_xhat2)[5];
  real_T *c13_b_Lambda2;
  real_T (*c13_b_meas)[2];
  real_T (*c13_b_xhat20)[5];
  c13_b_P_out = (real_T (*)[25])ssGetOutputPortSignal(chartInstance->S, 3);
  c13_b_Lambda2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c13_b_xhat2 = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
  c13_b_meas = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 1);
  c13_b_xhat20 = (real_T (*)[5])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 6U, chartInstance->c13_sfEvent);
  for (c13_i14 = 0; c13_i14 < 5; c13_i14++) {
    c13_xhat20[c13_i14] = (*c13_b_xhat20)[c13_i14];
  }

  for (c13_i15 = 0; c13_i15 < 2; c13_i15++) {
    c13_meas[c13_i15] = (*c13_b_meas)[c13_i15];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 44U, 44U, c13_debug_family_names,
    c13_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_T, 0U, c13_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_L, 1U, c13_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_alpha, 2U, c13_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_kappa, 3U, c13_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_beta, 4U, c13_e_sf_marshallOut,
    c13_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_lambda, 5U, c13_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_c, 6U, c13_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_Wm, 7U, c13_h_sf_marshallOut,
    c13_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_Wc, 8U, c13_h_sf_marshallOut,
    c13_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c13_H, 9U, c13_o_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_meas_noise1, 10U, c13_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_a_process_noise, 11U, c13_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_w_process_noise, 12U, c13_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_Sv1, 13U, c13_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_Sw1, 14U, c13_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_Sw2, 15U, c13_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c13_Sw, 16U, c13_n_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c13_Bw, 17U, c13_m_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c13_Dw, 18U, c13_j_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_x, 19U, c13_f_sf_marshallOut,
    c13_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_sP, 20U, c13_d_sf_marshallOut,
    c13_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_chi_p, 21U, c13_l_sf_marshallOut,
    c13_l_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_chi_m, 22U, c13_l_sf_marshallOut,
    c13_l_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_i, 23U, c13_e_sf_marshallOut,
    c13_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_x_m, 24U, c13_f_sf_marshallOut,
    c13_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_P_m, 25U, c13_d_sf_marshallOut,
    c13_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_psi_m, 26U, c13_k_sf_marshallOut,
    c13_k_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_y_m, 27U, c13_g_sf_marshallOut,
    c13_h_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_Pyy, 28U, c13_j_sf_marshallOut,
    c13_j_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_Pxy, 29U, c13_i_sf_marshallOut,
    c13_i_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_K, 30U, c13_i_sf_marshallOut,
    c13_i_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_innov, 31U, c13_g_sf_marshallOut,
    c13_h_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_c_out, 32U, c13_e_sf_marshallOut,
    c13_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_Wm_out, 33U, c13_h_sf_marshallOut,
    c13_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_nargin, 34U, c13_e_sf_marshallOut,
    c13_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_nargout, 35U, c13_e_sf_marshallOut,
    c13_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c13_xhat20, 36U, c13_f_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c13_meas, 37U, c13_g_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_xhat2, 38U, c13_f_sf_marshallOut,
    c13_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_Lambda2, 39U, c13_e_sf_marshallOut,
    c13_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_P_out, 40U, c13_d_sf_marshallOut,
    c13_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(chartInstance->c13_P2, 41U,
    c13_c_sf_marshallOut, c13_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(chartInstance->c13_W2, 42U,
    c13_b_sf_marshallOut, c13_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(chartInstance->c13_V2, 43U,
    c13_sf_marshallOut, c13_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 3);
  c13_T = 0.05;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 4);
  c13_L = 5.0;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 6);
  c13_alpha = 0.001;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 7);
  c13_kappa = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 8);
  c13_beta = 2.0;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 10);
  c13_lambda = -4.999995;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 11);
  c13_c = 4.9999999998107114E-6;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 16);
  for (c13_i16 = 0; c13_i16 < 11; c13_i16++) {
    c13_Wm[c13_i16] = 100000.00000378577;
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 17);
  for (c13_i17 = 0; c13_i17 < 11; c13_i17++) {
    c13_Wc[c13_i17] = c13_Wm[c13_i17];
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 18);
  c13_Wm[0] = -999999.00003785768;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 19);
  c13_Wc[0] = (-999998.00003785768 - c13_mpower(chartInstance, 0.001)) +
    c13_beta;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 20);
  c13_c = 0.0022360679774574635;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 22);
  for (c13_i18 = 0; c13_i18 < 10; c13_i18++) {
    c13_H[c13_i18] = c13_a[c13_i18];
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 26);
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 27);
  if (CV_EML_IF(0, 1, 0, !chartInstance->c13_P2_not_empty)) {
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 32);
    c13_meas_noise1 = 1.5;
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 33);
    c13_a_process_noise = 0.015;
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 34);
    c13_w_process_noise = 0.0015;
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 36);
    c13_Sv1 = 2.25;
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 37);
    c13_Sw1 = 0.000225;
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 38);
    c13_Sw2 = 2.25E-6;
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 39);
    for (c13_i19 = 0; c13_i19 < 9; c13_i19++) {
      c13_Sw[c13_i19] = c13_dv5[c13_i19];
    }

    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 42);
    for (c13_i20 = 0; c13_i20 < 15; c13_i20++) {
      c13_Bw[c13_i20] = c13_dv6[c13_i20];
    }

    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 44);
    for (c13_i21 = 0; c13_i21 < 4; c13_i21++) {
      c13_Dw[c13_i21] = c13_b[c13_i21];
    }

    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 45);
    c13_b_eml_scalar_eg(chartInstance);
    c13_b_eml_scalar_eg(chartInstance);
    for (c13_i22 = 0; c13_i22 < 25; c13_i22++) {
      chartInstance->c13_W2[c13_i22] = c13_y[c13_i22];
    }

    chartInstance->c13_W2_not_empty = TRUE;
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 46);
    c13_d_eml_scalar_eg(chartInstance);
    c13_d_eml_scalar_eg(chartInstance);
    for (c13_i23 = 0; c13_i23 < 4; c13_i23++) {
      chartInstance->c13_V2[c13_i23] = c13_b_y[c13_i23];
    }

    chartInstance->c13_V2_not_empty = TRUE;
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 48);
    c13_eye(chartInstance, c13_dv7);
    for (c13_i24 = 0; c13_i24 < 25; c13_i24++) {
      chartInstance->c13_P2[c13_i24] = c13_dv7[c13_i24];
    }

    chartInstance->c13_P2_not_empty = TRUE;
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 56);
  for (c13_i25 = 0; c13_i25 < 5; c13_i25++) {
    c13_x[c13_i25] = c13_xhat20[c13_i25];
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 57);
  for (c13_i26 = 0; c13_i26 < 25; c13_i26++) {
    c13_hoistedGlobal[c13_i26] = chartInstance->c13_P2[c13_i26];
  }

  for (c13_j = 1; c13_j < 6; c13_j++) {
    c13_b_j = c13_j;
  }

  c13_info = c13_b_eml_matlab_zpotrf(chartInstance, c13_hoistedGlobal);
  c13_b_info = c13_info;
  c13_c_info = c13_b_info;
  c13_d_info = c13_c_info;
  if (c13_d_info == 0) {
    c13_jmax = 6;
  } else {
    c13_c_eml_error(chartInstance);
    c13_b_a = c13_d_info;
    c13_jmax = c13_b_a;
  }

  c13_b_jmax = c13_jmax - 1;
  c13_b_b = c13_b_jmax;
  c13_c_b = c13_b_b;
  if (2 > c13_c_b) {
    c13_overflow = FALSE;
  } else {
    c13_overflow = (c13_c_b > 2147483646);
  }

  if (c13_overflow) {
    c13_check_forloop_overflow_error(chartInstance, c13_overflow);
  }

  for (c13_c_j = 2; c13_c_j <= c13_b_jmax; c13_c_j++) {
    c13_b_j = c13_c_j;
    c13_c_a = c13_b_j - 1;
    c13_i27 = c13_c_a;
    c13_d_b = c13_i27;
    c13_e_b = c13_d_b;
    if (1 > c13_e_b) {
      c13_b_overflow = FALSE;
    } else {
      c13_b_overflow = (c13_e_b > 2147483646);
    }

    if (c13_b_overflow) {
      c13_check_forloop_overflow_error(chartInstance, c13_b_overflow);
    }

    for (c13_b_i = 1; c13_b_i <= c13_i27; c13_b_i++) {
      c13_c_i = c13_b_i;
      c13_hoistedGlobal[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c13_c_i), 1, 5, 1, 0) + 5 *
                         (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c13_b_j), 1, 5, 2, 0) - 1)) - 1] = 0.0;
    }
  }

  for (c13_i28 = 0; c13_i28 < 25; c13_i28++) {
    c13_sP[c13_i28] = c13_hoistedGlobal[c13_i28];
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 58);
  for (c13_i29 = 0; c13_i29 < 5; c13_i29++) {
    c13_d_a[c13_i29] = c13_x[c13_i29];
  }

  for (c13_i30 = 0; c13_i30 < 5; c13_i30++) {
    c13_i31 = 0;
    for (c13_i32 = 0; c13_i32 < 5; c13_i32++) {
      c13_hoistedGlobal[c13_i31 + c13_i30] = c13_d_a[c13_i30];
      c13_i31 += 5;
    }
  }

  for (c13_i33 = 0; c13_i33 < 25; c13_i33++) {
    c13_f_b[c13_i33] = c13_sP[c13_i33];
  }

  for (c13_i34 = 0; c13_i34 < 25; c13_i34++) {
    c13_f_b[c13_i34] *= 0.0022360679774574635;
  }

  for (c13_i35 = 0; c13_i35 < 5; c13_i35++) {
    c13_d_a[c13_i35] = c13_x[c13_i35];
  }

  for (c13_i36 = 0; c13_i36 < 5; c13_i36++) {
    c13_i37 = 0;
    for (c13_i38 = 0; c13_i38 < 5; c13_i38++) {
      c13_c_y[c13_i37 + c13_i36] = c13_d_a[c13_i36];
      c13_i37 += 5;
    }
  }

  for (c13_i39 = 0; c13_i39 < 25; c13_i39++) {
    c13_g_b[c13_i39] = c13_sP[c13_i39];
  }

  for (c13_i40 = 0; c13_i40 < 25; c13_i40++) {
    c13_g_b[c13_i40] *= 0.0022360679774574635;
  }

  for (c13_i41 = 0; c13_i41 < 5; c13_i41++) {
    c13_chi_p[c13_i41] = c13_x[c13_i41];
  }

  c13_i42 = 0;
  for (c13_i43 = 0; c13_i43 < 5; c13_i43++) {
    for (c13_i44 = 0; c13_i44 < 5; c13_i44++) {
      c13_chi_p[(c13_i44 + c13_i42) + 5] = c13_hoistedGlobal[c13_i44 + c13_i42]
        + c13_f_b[c13_i44 + c13_i42];
    }

    c13_i42 += 5;
  }

  c13_i45 = 0;
  for (c13_i46 = 0; c13_i46 < 5; c13_i46++) {
    for (c13_i47 = 0; c13_i47 < 5; c13_i47++) {
      c13_chi_p[(c13_i47 + c13_i45) + 30] = c13_c_y[c13_i47 + c13_i45] -
        c13_g_b[c13_i47 + c13_i45];
    }

    c13_i45 += 5;
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 61);
  for (c13_i48 = 0; c13_i48 < 55; c13_i48++) {
    c13_chi_m[c13_i48] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 62);
  c13_i = 1.0;
  c13_d_i = 0;
  while (c13_d_i < 11) {
    c13_i = 1.0 + (real_T)c13_d_i;
    CV_EML_FOR(0, 1, 0, 1);
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 63);
    c13_e_a = c13_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_d_y = c13_e_a * 0.05;
    c13_b_x = c13_d_y;
    c13_c_x = c13_b_x;
    c13_c_x = muDoubleScalarSin(c13_c_x);
    c13_f_a = c13_chi_p[1 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_h_b = c13_c_x;
    c13_e_y = c13_f_a * c13_h_b;
    c13_A = c13_e_y;
    c13_B = c13_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_d_x = c13_A;
    c13_f_y = c13_B;
    c13_e_x = c13_d_x;
    c13_g_y = c13_f_y;
    c13_h_y = c13_e_x / c13_g_y;
    c13_g_a = c13_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_i_y = c13_g_a * 0.05;
    c13_f_x = c13_i_y;
    c13_g_x = c13_f_x;
    c13_g_x = muDoubleScalarCos(c13_g_x);
    c13_h_a = c13_chi_p[3 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_i_b = 1.0 - c13_g_x;
    c13_j_y = c13_h_a * c13_i_b;
    c13_b_A = c13_j_y;
    c13_b_B = c13_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_h_x = c13_b_A;
    c13_k_y = c13_b_B;
    c13_i_x = c13_h_x;
    c13_l_y = c13_k_y;
    c13_m_y = c13_i_x / c13_l_y;
    c13_chi_m[5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m",
      (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)] = (c13_chi_p[5
      * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)] + c13_h_y) - c13_m_y;
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 64);
    c13_i_a = c13_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_n_y = c13_i_a * 0.05;
    c13_j_x = c13_n_y;
    c13_k_x = c13_j_x;
    c13_k_x = muDoubleScalarCos(c13_k_x);
    c13_j_a = c13_k_x;
    c13_j_b = c13_chi_p[1 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_o_y = c13_j_a * c13_j_b;
    c13_k_a = c13_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_p_y = c13_k_a * 0.05;
    c13_l_x = c13_p_y;
    c13_m_x = c13_l_x;
    c13_m_x = muDoubleScalarSin(c13_m_x);
    c13_l_a = c13_m_x;
    c13_k_b = c13_chi_p[3 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_q_y = c13_l_a * c13_k_b;
    c13_chi_m[1 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m",
      (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)] = c13_o_y -
      c13_q_y;
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 65);
    c13_m_a = c13_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_r_y = c13_m_a * 0.05;
    c13_n_x = c13_r_y;
    c13_o_x = c13_n_x;
    c13_o_x = muDoubleScalarCos(c13_o_x);
    c13_n_a = c13_chi_p[1 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_l_b = 1.0 - c13_o_x;
    c13_s_y = c13_n_a * c13_l_b;
    c13_c_A = c13_s_y;
    c13_c_B = c13_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_p_x = c13_c_A;
    c13_t_y = c13_c_B;
    c13_q_x = c13_p_x;
    c13_u_y = c13_t_y;
    c13_v_y = c13_q_x / c13_u_y;
    c13_o_a = c13_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_w_y = c13_o_a * 0.05;
    c13_r_x = c13_w_y;
    c13_s_x = c13_r_x;
    c13_s_x = muDoubleScalarSin(c13_s_x);
    c13_p_a = c13_chi_p[3 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_m_b = c13_s_x;
    c13_x_y = c13_p_a * c13_m_b;
    c13_d_A = c13_x_y;
    c13_d_B = c13_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_t_x = c13_d_A;
    c13_y_y = c13_d_B;
    c13_u_x = c13_t_x;
    c13_ab_y = c13_y_y;
    c13_bb_y = c13_u_x / c13_ab_y;
    c13_chi_m[2 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m",
      (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)] = (c13_v_y +
      c13_chi_p[2 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p",
      (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)]) + c13_bb_y;
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 66);
    c13_q_a = c13_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_cb_y = c13_q_a * 0.05;
    c13_v_x = c13_cb_y;
    c13_w_x = c13_v_x;
    c13_w_x = muDoubleScalarSin(c13_w_x);
    c13_r_a = c13_chi_p[1 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_n_b = c13_w_x;
    c13_db_y = c13_r_a * c13_n_b;
    c13_s_a = c13_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_eb_y = c13_s_a * 0.05;
    c13_x_x = c13_eb_y;
    c13_y_x = c13_x_x;
    c13_y_x = muDoubleScalarCos(c13_y_x);
    c13_t_a = c13_chi_p[3 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_o_b = c13_y_x;
    c13_fb_y = c13_t_a * c13_o_b;
    c13_chi_m[3 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m",
      (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)] = c13_db_y +
      c13_fb_y;
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 67);
    c13_chi_m[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m",
      (int32_T)_SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)] = c13_chi_p[4
      + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
              _SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1)];
    c13_d_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_EML_FOR(0, 1, 0, 0);
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 69);
  for (c13_i49 = 0; c13_i49 < 55; c13_i49++) {
    c13_u_a[c13_i49] = c13_chi_m[c13_i49];
  }

  for (c13_i50 = 0; c13_i50 < 11; c13_i50++) {
    c13_p_b[c13_i50] = c13_Wm[c13_i50];
  }

  c13_e_eml_scalar_eg(chartInstance);
  c13_e_eml_scalar_eg(chartInstance);
  for (c13_i51 = 0; c13_i51 < 5; c13_i51++) {
    c13_x_m[c13_i51] = 0.0;
  }

  for (c13_i52 = 0; c13_i52 < 5; c13_i52++) {
    c13_x_m[c13_i52] = 0.0;
  }

  for (c13_i53 = 0; c13_i53 < 5; c13_i53++) {
    c13_d_a[c13_i53] = c13_x_m[c13_i53];
  }

  for (c13_i54 = 0; c13_i54 < 5; c13_i54++) {
    c13_x_m[c13_i54] = c13_d_a[c13_i54];
  }

  for (c13_i55 = 0; c13_i55 < 5; c13_i55++) {
    c13_d_a[c13_i55] = c13_x_m[c13_i55];
  }

  for (c13_i56 = 0; c13_i56 < 5; c13_i56++) {
    c13_x_m[c13_i56] = c13_d_a[c13_i56];
  }

  for (c13_i57 = 0; c13_i57 < 5; c13_i57++) {
    c13_x_m[c13_i57] = 0.0;
    c13_i58 = 0;
    for (c13_i59 = 0; c13_i59 < 11; c13_i59++) {
      c13_x_m[c13_i57] += c13_u_a[c13_i58 + c13_i57] * c13_p_b[c13_i59];
      c13_i58 += 5;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 70);
  for (c13_i60 = 0; c13_i60 < 25; c13_i60++) {
    c13_P_m[c13_i60] = chartInstance->c13_W2[c13_i60];
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 71);
  c13_i = 1.0;
  c13_e_i = 0;
  while (c13_e_i < 11) {
    c13_i = 1.0 + (real_T)c13_e_i;
    CV_EML_FOR(0, 1, 1, 1);
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 72);
    c13_v_a = c13_Wc[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("Wc", (int32_T)
      _SFD_INTEGER_CHECK("i", c13_i), 1, 11, 1, 0) - 1];
    c13_f_i = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m", (int32_T)
      _SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1;
    for (c13_i61 = 0; c13_i61 < 5; c13_i61++) {
      c13_d_a[c13_i61] = c13_chi_m[c13_i61 + 5 * c13_f_i] - c13_x_m[c13_i61];
    }

    for (c13_i62 = 0; c13_i62 < 5; c13_i62++) {
      c13_d_a[c13_i62] *= c13_v_a;
    }

    c13_g_i = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m", (int32_T)
      _SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1;
    for (c13_i63 = 0; c13_i63 < 5; c13_i63++) {
      c13_q_b[c13_i63] = c13_chi_m[c13_i63 + 5 * c13_g_i] - c13_x_m[c13_i63];
    }

    for (c13_i64 = 0; c13_i64 < 5; c13_i64++) {
      c13_i65 = 0;
      for (c13_i66 = 0; c13_i66 < 5; c13_i66++) {
        c13_hoistedGlobal[c13_i65 + c13_i64] = c13_d_a[c13_i64] *
          c13_q_b[c13_i66];
        c13_i65 += 5;
      }
    }

    for (c13_i67 = 0; c13_i67 < 25; c13_i67++) {
      c13_P_m[c13_i67] += c13_hoistedGlobal[c13_i67];
    }

    c13_e_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_EML_FOR(0, 1, 1, 0);
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 74);
  for (c13_i68 = 0; c13_i68 < 55; c13_i68++) {
    c13_u_a[c13_i68] = c13_chi_m[c13_i68];
  }

  c13_f_eml_scalar_eg(chartInstance);
  c13_f_eml_scalar_eg(chartInstance);
  for (c13_i69 = 0; c13_i69 < 22; c13_i69++) {
    c13_psi_m[c13_i69] = 0.0;
  }

  for (c13_i70 = 0; c13_i70 < 22; c13_i70++) {
    c13_psi_m[c13_i70] = 0.0;
  }

  for (c13_i71 = 0; c13_i71 < 22; c13_i71++) {
    c13_C[c13_i71] = c13_psi_m[c13_i71];
  }

  for (c13_i72 = 0; c13_i72 < 22; c13_i72++) {
    c13_psi_m[c13_i72] = c13_C[c13_i72];
  }

  for (c13_i73 = 0; c13_i73 < 22; c13_i73++) {
    c13_C[c13_i73] = c13_psi_m[c13_i73];
  }

  for (c13_i74 = 0; c13_i74 < 22; c13_i74++) {
    c13_psi_m[c13_i74] = c13_C[c13_i74];
  }

  for (c13_i75 = 0; c13_i75 < 2; c13_i75++) {
    c13_i76 = 0;
    c13_i77 = 0;
    for (c13_i78 = 0; c13_i78 < 11; c13_i78++) {
      c13_psi_m[c13_i76 + c13_i75] = 0.0;
      c13_i79 = 0;
      for (c13_i80 = 0; c13_i80 < 5; c13_i80++) {
        c13_psi_m[c13_i76 + c13_i75] += c13_a[c13_i79 + c13_i75] *
          c13_u_a[c13_i80 + c13_i77];
        c13_i79 += 2;
      }

      c13_i76 += 2;
      c13_i77 += 5;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 75);
  for (c13_i81 = 0; c13_i81 < 22; c13_i81++) {
    c13_C[c13_i81] = c13_psi_m[c13_i81];
  }

  for (c13_i82 = 0; c13_i82 < 11; c13_i82++) {
    c13_p_b[c13_i82] = c13_Wm[c13_i82];
  }

  c13_g_eml_scalar_eg(chartInstance);
  c13_g_eml_scalar_eg(chartInstance);
  for (c13_i83 = 0; c13_i83 < 2; c13_i83++) {
    c13_y_m[c13_i83] = 0.0;
  }

  for (c13_i84 = 0; c13_i84 < 2; c13_i84++) {
    c13_y_m[c13_i84] = 0.0;
  }

  for (c13_i85 = 0; c13_i85 < 2; c13_i85++) {
    c13_b_C[c13_i85] = c13_y_m[c13_i85];
  }

  for (c13_i86 = 0; c13_i86 < 2; c13_i86++) {
    c13_y_m[c13_i86] = c13_b_C[c13_i86];
  }

  for (c13_i87 = 0; c13_i87 < 2; c13_i87++) {
    c13_b_C[c13_i87] = c13_y_m[c13_i87];
  }

  for (c13_i88 = 0; c13_i88 < 2; c13_i88++) {
    c13_y_m[c13_i88] = c13_b_C[c13_i88];
  }

  for (c13_i89 = 0; c13_i89 < 2; c13_i89++) {
    c13_y_m[c13_i89] = 0.0;
    c13_i90 = 0;
    for (c13_i91 = 0; c13_i91 < 11; c13_i91++) {
      c13_y_m[c13_i89] += c13_C[c13_i90 + c13_i89] * c13_p_b[c13_i91];
      c13_i90 += 2;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 78);
  for (c13_i92 = 0; c13_i92 < 4; c13_i92++) {
    c13_Pyy[c13_i92] = chartInstance->c13_V2[c13_i92];
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 79);
  for (c13_i93 = 0; c13_i93 < 10; c13_i93++) {
    c13_Pxy[c13_i93] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 80);
  c13_i = 1.0;
  c13_h_i = 0;
  while (c13_h_i < 11) {
    c13_i = 1.0 + (real_T)c13_h_i;
    CV_EML_FOR(0, 1, 2, 1);
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 81);
    c13_w_a = c13_Wc[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("Wc", (int32_T)
      _SFD_INTEGER_CHECK("i", c13_i), 1, 11, 1, 0) - 1];
    c13_i_i = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("psi_m", (int32_T)
      _SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1;
    for (c13_i94 = 0; c13_i94 < 2; c13_i94++) {
      c13_b_C[c13_i94] = c13_psi_m[c13_i94 + (c13_i_i << 1)] - c13_y_m[c13_i94];
    }

    for (c13_i95 = 0; c13_i95 < 2; c13_i95++) {
      c13_b_C[c13_i95] *= c13_w_a;
    }

    c13_j_i = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("psi_m", (int32_T)
      _SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1;
    for (c13_i96 = 0; c13_i96 < 2; c13_i96++) {
      c13_x_a[c13_i96] = c13_psi_m[c13_i96 + (c13_j_i << 1)] - c13_y_m[c13_i96];
    }

    for (c13_i97 = 0; c13_i97 < 2; c13_i97++) {
      c13_i98 = 0;
      for (c13_i99 = 0; c13_i99 < 2; c13_i99++) {
        c13_r_b[c13_i98 + c13_i97] = c13_b_C[c13_i97] * c13_x_a[c13_i99];
        c13_i98 += 2;
      }
    }

    for (c13_i100 = 0; c13_i100 < 4; c13_i100++) {
      c13_Pyy[c13_i100] += c13_r_b[c13_i100];
    }

    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 82);
    c13_y_a = c13_Wc[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("Wc", (int32_T)
      _SFD_INTEGER_CHECK("i", c13_i), 1, 11, 1, 0) - 1];
    c13_k_i = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m", (int32_T)
      _SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1;
    for (c13_i101 = 0; c13_i101 < 5; c13_i101++) {
      c13_d_a[c13_i101] = c13_chi_m[c13_i101 + 5 * c13_k_i] - c13_x_m[c13_i101];
    }

    for (c13_i102 = 0; c13_i102 < 5; c13_i102++) {
      c13_d_a[c13_i102] *= c13_y_a;
    }

    c13_l_i = (int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("psi_m", (int32_T)
      _SFD_INTEGER_CHECK("i", c13_i), 1, 11, 2, 0) - 1;
    for (c13_i103 = 0; c13_i103 < 2; c13_i103++) {
      c13_x_a[c13_i103] = c13_psi_m[c13_i103 + (c13_l_i << 1)] -
        c13_y_m[c13_i103];
    }

    for (c13_i104 = 0; c13_i104 < 5; c13_i104++) {
      c13_i105 = 0;
      for (c13_i106 = 0; c13_i106 < 2; c13_i106++) {
        c13_c_C[c13_i105 + c13_i104] = c13_d_a[c13_i104] * c13_x_a[c13_i106];
        c13_i105 += 5;
      }
    }

    for (c13_i107 = 0; c13_i107 < 10; c13_i107++) {
      c13_Pxy[c13_i107] += c13_c_C[c13_i107];
    }

    c13_h_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_EML_FOR(0, 1, 2, 0);
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 85);
  for (c13_i108 = 0; c13_i108 < 10; c13_i108++) {
    c13_ab_a[c13_i108] = c13_Pxy[c13_i108];
  }

  for (c13_i109 = 0; c13_i109 < 4; c13_i109++) {
    c13_b_Pyy[c13_i109] = c13_Pyy[c13_i109];
  }

  c13_inv(chartInstance, c13_b_Pyy, c13_r_b);
  c13_h_eml_scalar_eg(chartInstance);
  c13_h_eml_scalar_eg(chartInstance);
  for (c13_i110 = 0; c13_i110 < 10; c13_i110++) {
    c13_K[c13_i110] = 0.0;
  }

  for (c13_i111 = 0; c13_i111 < 10; c13_i111++) {
    c13_K[c13_i111] = 0.0;
  }

  for (c13_i112 = 0; c13_i112 < 10; c13_i112++) {
    c13_c_C[c13_i112] = c13_K[c13_i112];
  }

  for (c13_i113 = 0; c13_i113 < 10; c13_i113++) {
    c13_K[c13_i113] = c13_c_C[c13_i113];
  }

  for (c13_i114 = 0; c13_i114 < 10; c13_i114++) {
    c13_c_C[c13_i114] = c13_K[c13_i114];
  }

  for (c13_i115 = 0; c13_i115 < 10; c13_i115++) {
    c13_K[c13_i115] = c13_c_C[c13_i115];
  }

  for (c13_i116 = 0; c13_i116 < 5; c13_i116++) {
    c13_i117 = 0;
    c13_i118 = 0;
    for (c13_i119 = 0; c13_i119 < 2; c13_i119++) {
      c13_K[c13_i117 + c13_i116] = 0.0;
      c13_i120 = 0;
      for (c13_i121 = 0; c13_i121 < 2; c13_i121++) {
        c13_K[c13_i117 + c13_i116] += c13_ab_a[c13_i120 + c13_i116] *
          c13_r_b[c13_i121 + c13_i118];
        c13_i120 += 5;
      }

      c13_i117 += 5;
      c13_i118 += 2;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 86);
  for (c13_i122 = 0; c13_i122 < 10; c13_i122++) {
    c13_ab_a[c13_i122] = c13_K[c13_i122];
  }

  for (c13_i123 = 0; c13_i123 < 2; c13_i123++) {
    c13_b_C[c13_i123] = c13_meas[c13_i123] - c13_y_m[c13_i123];
  }

  c13_i_eml_scalar_eg(chartInstance);
  c13_i_eml_scalar_eg(chartInstance);
  for (c13_i124 = 0; c13_i124 < 5; c13_i124++) {
    c13_d_a[c13_i124] = 0.0;
    c13_i125 = 0;
    for (c13_i126 = 0; c13_i126 < 2; c13_i126++) {
      c13_d_a[c13_i124] += c13_ab_a[c13_i125 + c13_i124] * c13_b_C[c13_i126];
      c13_i125 += 5;
    }
  }

  for (c13_i127 = 0; c13_i127 < 5; c13_i127++) {
    c13_xhat2[c13_i127] = c13_x_m[c13_i127] + c13_d_a[c13_i127];
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 87);
  for (c13_i128 = 0; c13_i128 < 10; c13_i128++) {
    c13_ab_a[c13_i128] = c13_K[c13_i128];
  }

  for (c13_i129 = 0; c13_i129 < 4; c13_i129++) {
    c13_r_b[c13_i129] = c13_Pyy[c13_i129];
  }

  c13_h_eml_scalar_eg(chartInstance);
  c13_h_eml_scalar_eg(chartInstance);
  for (c13_i130 = 0; c13_i130 < 5; c13_i130++) {
    c13_i131 = 0;
    c13_i132 = 0;
    for (c13_i133 = 0; c13_i133 < 2; c13_i133++) {
      c13_c_C[c13_i131 + c13_i130] = 0.0;
      c13_i134 = 0;
      for (c13_i135 = 0; c13_i135 < 2; c13_i135++) {
        c13_c_C[c13_i131 + c13_i130] += c13_ab_a[c13_i134 + c13_i130] *
          c13_r_b[c13_i135 + c13_i132];
        c13_i134 += 5;
      }

      c13_i131 += 5;
      c13_i132 += 2;
    }
  }

  c13_i136 = 0;
  for (c13_i137 = 0; c13_i137 < 5; c13_i137++) {
    c13_i138 = 0;
    for (c13_i139 = 0; c13_i139 < 2; c13_i139++) {
      c13_s_b[c13_i139 + c13_i136] = c13_K[c13_i138 + c13_i137];
      c13_i138 += 5;
    }

    c13_i136 += 2;
  }

  c13_j_eml_scalar_eg(chartInstance);
  c13_j_eml_scalar_eg(chartInstance);
  for (c13_i140 = 0; c13_i140 < 5; c13_i140++) {
    c13_i141 = 0;
    c13_i142 = 0;
    for (c13_i143 = 0; c13_i143 < 5; c13_i143++) {
      c13_hoistedGlobal[c13_i141 + c13_i140] = 0.0;
      c13_i144 = 0;
      for (c13_i145 = 0; c13_i145 < 2; c13_i145++) {
        c13_hoistedGlobal[c13_i141 + c13_i140] += c13_c_C[c13_i144 + c13_i140] *
          c13_s_b[c13_i145 + c13_i142];
        c13_i144 += 5;
      }

      c13_i141 += 5;
      c13_i142 += 2;
    }
  }

  for (c13_i146 = 0; c13_i146 < 25; c13_i146++) {
    chartInstance->c13_P2[c13_i146] = c13_P_m[c13_i146] -
      c13_hoistedGlobal[c13_i146];
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 89);
  for (c13_i147 = 0; c13_i147 < 2; c13_i147++) {
    c13_innov[c13_i147] = c13_meas[c13_i147] - c13_y_m[c13_i147];
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 90);
  c13_t_b = c13_Pyy[0];
  c13_gb_y = 39.478417604357432 * c13_t_b;
  c13_bb_a = c13_gb_y;
  c13_u_b = c13_Pyy[3];
  c13_hb_y = c13_bb_a * c13_u_b;
  c13_e_B = c13_hb_y;
  c13_b_sqrt(chartInstance, &c13_e_B);
  c13_ib_y = c13_e_B;
  c13_jb_y = c13_ib_y;
  c13_kb_y = 1.0 / c13_jb_y;
  for (c13_i148 = 0; c13_i148 < 2; c13_i148++) {
    c13_x_a[c13_i148] = c13_innov[c13_i148];
  }

  for (c13_i149 = 0; c13_i149 < 4; c13_i149++) {
    c13_c_Pyy[c13_i149] = c13_Pyy[c13_i149];
  }

  c13_inv(chartInstance, c13_c_Pyy, c13_r_b);
  c13_k_eml_scalar_eg(chartInstance);
  c13_k_eml_scalar_eg(chartInstance);
  c13_i150 = 0;
  for (c13_i151 = 0; c13_i151 < 2; c13_i151++) {
    c13_lb_y[c13_i151] = 0.0;
    for (c13_i152 = 0; c13_i152 < 2; c13_i152++) {
      c13_lb_y[c13_i151] += c13_x_a[c13_i152] * c13_r_b[c13_i152 + c13_i150];
    }

    c13_i150 += 2;
  }

  for (c13_i153 = 0; c13_i153 < 2; c13_i153++) {
    c13_b_C[c13_i153] = c13_innov[c13_i153];
  }

  c13_l_eml_scalar_eg(chartInstance);
  c13_l_eml_scalar_eg(chartInstance);
  c13_mb_y = 0.0;
  for (c13_k = 1; c13_k < 3; c13_k++) {
    c13_b_k = c13_k;
    c13_mb_y += c13_lb_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c13_b_k), 1, 2, 1, 0) - 1] *
      c13_b_C[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c13_b_k), 1, 2, 1, 0) - 1];
  }

  c13_v_b = c13_mb_y;
  c13_nb_y = -0.5 * c13_v_b;
  c13_cb_a = c13_kb_y;
  c13_w_b = c13_expm(chartInstance, c13_nb_y);
  c13_Lambda2 = c13_cb_a * c13_w_b;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 97);
  for (c13_i154 = 0; c13_i154 < 25; c13_i154++) {
    c13_P_out[c13_i154] = chartInstance->c13_P2[c13_i154];
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 98);
  c13_c_out = c13_c;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 99);
  for (c13_i155 = 0; c13_i155 < 11; c13_i155++) {
    c13_Wm_out[c13_i155] = c13_Wm[c13_i155];
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, -99);
  _SFD_SYMBOL_SCOPE_POP();
  for (c13_i156 = 0; c13_i156 < 5; c13_i156++) {
    (*c13_b_xhat2)[c13_i156] = c13_xhat2[c13_i156];
  }

  *c13_b_Lambda2 = c13_Lambda2;
  for (c13_i157 = 0; c13_i157 < 25; c13_i157++) {
    (*c13_b_P_out)[c13_i157] = c13_P_out[c13_i157];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 6U, chartInstance->c13_sfEvent);
}

static void initSimStructsc13_IMM_UKF(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static void init_script_number_translation(uint32_T c13_machineNumber, uint32_T
  c13_chartNumber)
{
}

static const mxArray *c13_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i158;
  int32_T c13_i159;
  int32_T c13_i160;
  real_T c13_b_inData[4];
  int32_T c13_i161;
  int32_T c13_i162;
  int32_T c13_i163;
  real_T c13_u[4];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_i158 = 0;
  for (c13_i159 = 0; c13_i159 < 2; c13_i159++) {
    for (c13_i160 = 0; c13_i160 < 2; c13_i160++) {
      c13_b_inData[c13_i160 + c13_i158] = (*(real_T (*)[4])c13_inData)[c13_i160
        + c13_i158];
    }

    c13_i158 += 2;
  }

  c13_i161 = 0;
  for (c13_i162 = 0; c13_i162 < 2; c13_i162++) {
    for (c13_i163 = 0; c13_i163 < 2; c13_i163++) {
      c13_u[c13_i163 + c13_i161] = c13_b_inData[c13_i163 + c13_i161];
    }

    c13_i161 += 2;
  }

  c13_y = NULL;
  if (!chartInstance->c13_V2_not_empty) {
    sf_mex_assign(&c13_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 2, 2, 2),
                  FALSE);
  }

  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static void c13_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_b_V2, const char_T *c13_identifier, real_T c13_y[4])
{
  emlrtMsgIdentifier c13_thisId;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_b_V2), &c13_thisId, c13_y);
  sf_mex_destroy(&c13_b_V2);
}

static void c13_b_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[4])
{
  real_T c13_dv8[4];
  int32_T c13_i164;
  if (mxIsEmpty(c13_u)) {
    chartInstance->c13_V2_not_empty = FALSE;
  } else {
    chartInstance->c13_V2_not_empty = TRUE;
    sf_mex_import(c13_parentId, sf_mex_dup(c13_u), c13_dv8, 1, 0, 0U, 1, 0U, 2,
                  2, 2);
    for (c13_i164 = 0; c13_i164 < 4; c13_i164++) {
      c13_y[c13_i164] = c13_dv8[c13_i164];
    }
  }

  sf_mex_destroy(&c13_u);
}

static void c13_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_b_V2;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y[4];
  int32_T c13_i165;
  int32_T c13_i166;
  int32_T c13_i167;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_b_V2 = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_b_V2), &c13_thisId, c13_y);
  sf_mex_destroy(&c13_b_V2);
  c13_i165 = 0;
  for (c13_i166 = 0; c13_i166 < 2; c13_i166++) {
    for (c13_i167 = 0; c13_i167 < 2; c13_i167++) {
      (*(real_T (*)[4])c13_outData)[c13_i167 + c13_i165] = c13_y[c13_i167 +
        c13_i165];
    }

    c13_i165 += 2;
  }

  sf_mex_destroy(&c13_mxArrayInData);
}

static const mxArray *c13_b_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i168;
  int32_T c13_i169;
  int32_T c13_i170;
  real_T c13_b_inData[25];
  int32_T c13_i171;
  int32_T c13_i172;
  int32_T c13_i173;
  real_T c13_u[25];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_i168 = 0;
  for (c13_i169 = 0; c13_i169 < 5; c13_i169++) {
    for (c13_i170 = 0; c13_i170 < 5; c13_i170++) {
      c13_b_inData[c13_i170 + c13_i168] = (*(real_T (*)[25])c13_inData)[c13_i170
        + c13_i168];
    }

    c13_i168 += 5;
  }

  c13_i171 = 0;
  for (c13_i172 = 0; c13_i172 < 5; c13_i172++) {
    for (c13_i173 = 0; c13_i173 < 5; c13_i173++) {
      c13_u[c13_i173 + c13_i171] = c13_b_inData[c13_i173 + c13_i171];
    }

    c13_i171 += 5;
  }

  c13_y = NULL;
  if (!chartInstance->c13_W2_not_empty) {
    sf_mex_assign(&c13_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 2, 5, 5),
                  FALSE);
  }

  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static void c13_c_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_b_W2, const char_T *c13_identifier, real_T c13_y[25])
{
  emlrtMsgIdentifier c13_thisId;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_b_W2), &c13_thisId, c13_y);
  sf_mex_destroy(&c13_b_W2);
}

static void c13_d_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[25])
{
  real_T c13_dv9[25];
  int32_T c13_i174;
  if (mxIsEmpty(c13_u)) {
    chartInstance->c13_W2_not_empty = FALSE;
  } else {
    chartInstance->c13_W2_not_empty = TRUE;
    sf_mex_import(c13_parentId, sf_mex_dup(c13_u), c13_dv9, 1, 0, 0U, 1, 0U, 2,
                  5, 5);
    for (c13_i174 = 0; c13_i174 < 25; c13_i174++) {
      c13_y[c13_i174] = c13_dv9[c13_i174];
    }
  }

  sf_mex_destroy(&c13_u);
}

static void c13_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_b_W2;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y[25];
  int32_T c13_i175;
  int32_T c13_i176;
  int32_T c13_i177;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_b_W2 = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_b_W2), &c13_thisId, c13_y);
  sf_mex_destroy(&c13_b_W2);
  c13_i175 = 0;
  for (c13_i176 = 0; c13_i176 < 5; c13_i176++) {
    for (c13_i177 = 0; c13_i177 < 5; c13_i177++) {
      (*(real_T (*)[25])c13_outData)[c13_i177 + c13_i175] = c13_y[c13_i177 +
        c13_i175];
    }

    c13_i175 += 5;
  }

  sf_mex_destroy(&c13_mxArrayInData);
}

static const mxArray *c13_c_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i178;
  int32_T c13_i179;
  int32_T c13_i180;
  real_T c13_b_inData[25];
  int32_T c13_i181;
  int32_T c13_i182;
  int32_T c13_i183;
  real_T c13_u[25];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_i178 = 0;
  for (c13_i179 = 0; c13_i179 < 5; c13_i179++) {
    for (c13_i180 = 0; c13_i180 < 5; c13_i180++) {
      c13_b_inData[c13_i180 + c13_i178] = (*(real_T (*)[25])c13_inData)[c13_i180
        + c13_i178];
    }

    c13_i178 += 5;
  }

  c13_i181 = 0;
  for (c13_i182 = 0; c13_i182 < 5; c13_i182++) {
    for (c13_i183 = 0; c13_i183 < 5; c13_i183++) {
      c13_u[c13_i183 + c13_i181] = c13_b_inData[c13_i183 + c13_i181];
    }

    c13_i181 += 5;
  }

  c13_y = NULL;
  if (!chartInstance->c13_P2_not_empty) {
    sf_mex_assign(&c13_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 2, 5, 5),
                  FALSE);
  }

  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static void c13_e_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_b_P2, const char_T *c13_identifier, real_T c13_y[25])
{
  emlrtMsgIdentifier c13_thisId;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_b_P2), &c13_thisId, c13_y);
  sf_mex_destroy(&c13_b_P2);
}

static void c13_f_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[25])
{
  real_T c13_dv10[25];
  int32_T c13_i184;
  if (mxIsEmpty(c13_u)) {
    chartInstance->c13_P2_not_empty = FALSE;
  } else {
    chartInstance->c13_P2_not_empty = TRUE;
    sf_mex_import(c13_parentId, sf_mex_dup(c13_u), c13_dv10, 1, 0, 0U, 1, 0U, 2,
                  5, 5);
    for (c13_i184 = 0; c13_i184 < 25; c13_i184++) {
      c13_y[c13_i184] = c13_dv10[c13_i184];
    }
  }

  sf_mex_destroy(&c13_u);
}

static void c13_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_b_P2;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y[25];
  int32_T c13_i185;
  int32_T c13_i186;
  int32_T c13_i187;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_b_P2 = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_b_P2), &c13_thisId, c13_y);
  sf_mex_destroy(&c13_b_P2);
  c13_i185 = 0;
  for (c13_i186 = 0; c13_i186 < 5; c13_i186++) {
    for (c13_i187 = 0; c13_i187 < 5; c13_i187++) {
      (*(real_T (*)[25])c13_outData)[c13_i187 + c13_i185] = c13_y[c13_i187 +
        c13_i185];
    }

    c13_i185 += 5;
  }

  sf_mex_destroy(&c13_mxArrayInData);
}

static const mxArray *c13_d_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i188;
  int32_T c13_i189;
  int32_T c13_i190;
  real_T c13_b_inData[25];
  int32_T c13_i191;
  int32_T c13_i192;
  int32_T c13_i193;
  real_T c13_u[25];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_i188 = 0;
  for (c13_i189 = 0; c13_i189 < 5; c13_i189++) {
    for (c13_i190 = 0; c13_i190 < 5; c13_i190++) {
      c13_b_inData[c13_i190 + c13_i188] = (*(real_T (*)[25])c13_inData)[c13_i190
        + c13_i188];
    }

    c13_i188 += 5;
  }

  c13_i191 = 0;
  for (c13_i192 = 0; c13_i192 < 5; c13_i192++) {
    for (c13_i193 = 0; c13_i193 < 5; c13_i193++) {
      c13_u[c13_i193 + c13_i191] = c13_b_inData[c13_i193 + c13_i191];
    }

    c13_i191 += 5;
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 2, 5, 5), FALSE);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static void c13_g_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_P_out, const char_T *c13_identifier, real_T c13_y[25])
{
  emlrtMsgIdentifier c13_thisId;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_P_out), &c13_thisId,
    c13_y);
  sf_mex_destroy(&c13_P_out);
}

static void c13_h_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[25])
{
  real_T c13_dv11[25];
  int32_T c13_i194;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), c13_dv11, 1, 0, 0U, 1, 0U, 2, 5,
                5);
  for (c13_i194 = 0; c13_i194 < 25; c13_i194++) {
    c13_y[c13_i194] = c13_dv11[c13_i194];
  }

  sf_mex_destroy(&c13_u);
}

static void c13_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_P_out;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y[25];
  int32_T c13_i195;
  int32_T c13_i196;
  int32_T c13_i197;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_P_out = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_P_out), &c13_thisId,
    c13_y);
  sf_mex_destroy(&c13_P_out);
  c13_i195 = 0;
  for (c13_i196 = 0; c13_i196 < 5; c13_i196++) {
    for (c13_i197 = 0; c13_i197 < 5; c13_i197++) {
      (*(real_T (*)[25])c13_outData)[c13_i197 + c13_i195] = c13_y[c13_i197 +
        c13_i195];
    }

    c13_i195 += 5;
  }

  sf_mex_destroy(&c13_mxArrayInData);
}

static const mxArray *c13_e_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  real_T c13_u;
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_u = *(real_T *)c13_inData;
  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", &c13_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static real_T c13_i_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_Lambda2, const char_T *c13_identifier)
{
  real_T c13_y;
  emlrtMsgIdentifier c13_thisId;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_y = c13_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_Lambda2),
    &c13_thisId);
  sf_mex_destroy(&c13_Lambda2);
  return c13_y;
}

static real_T c13_j_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId)
{
  real_T c13_y;
  real_T c13_d0;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), &c13_d0, 1, 0, 0U, 0, 0U, 0);
  c13_y = c13_d0;
  sf_mex_destroy(&c13_u);
  return c13_y;
}

static void c13_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_Lambda2;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_Lambda2 = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_y = c13_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_Lambda2),
    &c13_thisId);
  sf_mex_destroy(&c13_Lambda2);
  *(real_T *)c13_outData = c13_y;
  sf_mex_destroy(&c13_mxArrayInData);
}

static const mxArray *c13_f_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i198;
  real_T c13_b_inData[5];
  int32_T c13_i199;
  real_T c13_u[5];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  for (c13_i198 = 0; c13_i198 < 5; c13_i198++) {
    c13_b_inData[c13_i198] = (*(real_T (*)[5])c13_inData)[c13_i198];
  }

  for (c13_i199 = 0; c13_i199 < 5; c13_i199++) {
    c13_u[c13_i199] = c13_b_inData[c13_i199];
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 1, 5), FALSE);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static void c13_k_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_xhat2, const char_T *c13_identifier, real_T c13_y[5])
{
  emlrtMsgIdentifier c13_thisId;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_l_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_xhat2), &c13_thisId,
    c13_y);
  sf_mex_destroy(&c13_xhat2);
}

static void c13_l_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[5])
{
  real_T c13_dv12[5];
  int32_T c13_i200;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), c13_dv12, 1, 0, 0U, 1, 0U, 1, 5);
  for (c13_i200 = 0; c13_i200 < 5; c13_i200++) {
    c13_y[c13_i200] = c13_dv12[c13_i200];
  }

  sf_mex_destroy(&c13_u);
}

static void c13_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_xhat2;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y[5];
  int32_T c13_i201;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_xhat2 = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_l_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_xhat2), &c13_thisId,
    c13_y);
  sf_mex_destroy(&c13_xhat2);
  for (c13_i201 = 0; c13_i201 < 5; c13_i201++) {
    (*(real_T (*)[5])c13_outData)[c13_i201] = c13_y[c13_i201];
  }

  sf_mex_destroy(&c13_mxArrayInData);
}

static const mxArray *c13_g_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i202;
  real_T c13_b_inData[2];
  int32_T c13_i203;
  real_T c13_u[2];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  for (c13_i202 = 0; c13_i202 < 2; c13_i202++) {
    c13_b_inData[c13_i202] = (*(real_T (*)[2])c13_inData)[c13_i202];
  }

  for (c13_i203 = 0; c13_i203 < 2; c13_i203++) {
    c13_u[c13_i203] = c13_b_inData[c13_i203];
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static const mxArray *c13_h_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i204;
  real_T c13_b_inData[11];
  int32_T c13_i205;
  real_T c13_u[11];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  for (c13_i204 = 0; c13_i204 < 11; c13_i204++) {
    c13_b_inData[c13_i204] = (*(real_T (*)[11])c13_inData)[c13_i204];
  }

  for (c13_i205 = 0; c13_i205 < 11; c13_i205++) {
    c13_u[c13_i205] = c13_b_inData[c13_i205];
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 1, 11), FALSE);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static void c13_m_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[11])
{
  real_T c13_dv13[11];
  int32_T c13_i206;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), c13_dv13, 1, 0, 0U, 1, 0U, 1,
                11);
  for (c13_i206 = 0; c13_i206 < 11; c13_i206++) {
    c13_y[c13_i206] = c13_dv13[c13_i206];
  }

  sf_mex_destroy(&c13_u);
}

static void c13_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_Wm_out;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y[11];
  int32_T c13_i207;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_Wm_out = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_m_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_Wm_out), &c13_thisId,
    c13_y);
  sf_mex_destroy(&c13_Wm_out);
  for (c13_i207 = 0; c13_i207 < 11; c13_i207++) {
    (*(real_T (*)[11])c13_outData)[c13_i207] = c13_y[c13_i207];
  }

  sf_mex_destroy(&c13_mxArrayInData);
}

static void c13_n_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[2])
{
  real_T c13_dv14[2];
  int32_T c13_i208;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), c13_dv14, 1, 0, 0U, 1, 0U, 1, 2);
  for (c13_i208 = 0; c13_i208 < 2; c13_i208++) {
    c13_y[c13_i208] = c13_dv14[c13_i208];
  }

  sf_mex_destroy(&c13_u);
}

static void c13_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_innov;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y[2];
  int32_T c13_i209;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_innov = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_n_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_innov), &c13_thisId,
    c13_y);
  sf_mex_destroy(&c13_innov);
  for (c13_i209 = 0; c13_i209 < 2; c13_i209++) {
    (*(real_T (*)[2])c13_outData)[c13_i209] = c13_y[c13_i209];
  }

  sf_mex_destroy(&c13_mxArrayInData);
}

static const mxArray *c13_i_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i210;
  int32_T c13_i211;
  int32_T c13_i212;
  real_T c13_b_inData[10];
  int32_T c13_i213;
  int32_T c13_i214;
  int32_T c13_i215;
  real_T c13_u[10];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_i210 = 0;
  for (c13_i211 = 0; c13_i211 < 2; c13_i211++) {
    for (c13_i212 = 0; c13_i212 < 5; c13_i212++) {
      c13_b_inData[c13_i212 + c13_i210] = (*(real_T (*)[10])c13_inData)[c13_i212
        + c13_i210];
    }

    c13_i210 += 5;
  }

  c13_i213 = 0;
  for (c13_i214 = 0; c13_i214 < 2; c13_i214++) {
    for (c13_i215 = 0; c13_i215 < 5; c13_i215++) {
      c13_u[c13_i215 + c13_i213] = c13_b_inData[c13_i215 + c13_i213];
    }

    c13_i213 += 5;
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 2, 5, 2), FALSE);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static void c13_o_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[10])
{
  real_T c13_dv15[10];
  int32_T c13_i216;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), c13_dv15, 1, 0, 0U, 1, 0U, 2, 5,
                2);
  for (c13_i216 = 0; c13_i216 < 10; c13_i216++) {
    c13_y[c13_i216] = c13_dv15[c13_i216];
  }

  sf_mex_destroy(&c13_u);
}

static void c13_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_K;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y[10];
  int32_T c13_i217;
  int32_T c13_i218;
  int32_T c13_i219;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_K = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_o_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_K), &c13_thisId, c13_y);
  sf_mex_destroy(&c13_K);
  c13_i217 = 0;
  for (c13_i218 = 0; c13_i218 < 2; c13_i218++) {
    for (c13_i219 = 0; c13_i219 < 5; c13_i219++) {
      (*(real_T (*)[10])c13_outData)[c13_i219 + c13_i217] = c13_y[c13_i219 +
        c13_i217];
    }

    c13_i217 += 5;
  }

  sf_mex_destroy(&c13_mxArrayInData);
}

static const mxArray *c13_j_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i220;
  int32_T c13_i221;
  int32_T c13_i222;
  real_T c13_b_inData[4];
  int32_T c13_i223;
  int32_T c13_i224;
  int32_T c13_i225;
  real_T c13_u[4];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_i220 = 0;
  for (c13_i221 = 0; c13_i221 < 2; c13_i221++) {
    for (c13_i222 = 0; c13_i222 < 2; c13_i222++) {
      c13_b_inData[c13_i222 + c13_i220] = (*(real_T (*)[4])c13_inData)[c13_i222
        + c13_i220];
    }

    c13_i220 += 2;
  }

  c13_i223 = 0;
  for (c13_i224 = 0; c13_i224 < 2; c13_i224++) {
    for (c13_i225 = 0; c13_i225 < 2; c13_i225++) {
      c13_u[c13_i225 + c13_i223] = c13_b_inData[c13_i225 + c13_i223];
    }

    c13_i223 += 2;
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 2, 2, 2), FALSE);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static void c13_p_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[4])
{
  real_T c13_dv16[4];
  int32_T c13_i226;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), c13_dv16, 1, 0, 0U, 1, 0U, 2, 2,
                2);
  for (c13_i226 = 0; c13_i226 < 4; c13_i226++) {
    c13_y[c13_i226] = c13_dv16[c13_i226];
  }

  sf_mex_destroy(&c13_u);
}

static void c13_j_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_Pyy;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y[4];
  int32_T c13_i227;
  int32_T c13_i228;
  int32_T c13_i229;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_Pyy = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_p_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_Pyy), &c13_thisId, c13_y);
  sf_mex_destroy(&c13_Pyy);
  c13_i227 = 0;
  for (c13_i228 = 0; c13_i228 < 2; c13_i228++) {
    for (c13_i229 = 0; c13_i229 < 2; c13_i229++) {
      (*(real_T (*)[4])c13_outData)[c13_i229 + c13_i227] = c13_y[c13_i229 +
        c13_i227];
    }

    c13_i227 += 2;
  }

  sf_mex_destroy(&c13_mxArrayInData);
}

static const mxArray *c13_k_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i230;
  int32_T c13_i231;
  int32_T c13_i232;
  real_T c13_b_inData[22];
  int32_T c13_i233;
  int32_T c13_i234;
  int32_T c13_i235;
  real_T c13_u[22];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_i230 = 0;
  for (c13_i231 = 0; c13_i231 < 11; c13_i231++) {
    for (c13_i232 = 0; c13_i232 < 2; c13_i232++) {
      c13_b_inData[c13_i232 + c13_i230] = (*(real_T (*)[22])c13_inData)[c13_i232
        + c13_i230];
    }

    c13_i230 += 2;
  }

  c13_i233 = 0;
  for (c13_i234 = 0; c13_i234 < 11; c13_i234++) {
    for (c13_i235 = 0; c13_i235 < 2; c13_i235++) {
      c13_u[c13_i235 + c13_i233] = c13_b_inData[c13_i235 + c13_i233];
    }

    c13_i233 += 2;
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 2, 2, 11),
                FALSE);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static void c13_q_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[22])
{
  real_T c13_dv17[22];
  int32_T c13_i236;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), c13_dv17, 1, 0, 0U, 1, 0U, 2, 2,
                11);
  for (c13_i236 = 0; c13_i236 < 22; c13_i236++) {
    c13_y[c13_i236] = c13_dv17[c13_i236];
  }

  sf_mex_destroy(&c13_u);
}

static void c13_k_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_psi_m;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y[22];
  int32_T c13_i237;
  int32_T c13_i238;
  int32_T c13_i239;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_psi_m = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_q_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_psi_m), &c13_thisId,
    c13_y);
  sf_mex_destroy(&c13_psi_m);
  c13_i237 = 0;
  for (c13_i238 = 0; c13_i238 < 11; c13_i238++) {
    for (c13_i239 = 0; c13_i239 < 2; c13_i239++) {
      (*(real_T (*)[22])c13_outData)[c13_i239 + c13_i237] = c13_y[c13_i239 +
        c13_i237];
    }

    c13_i237 += 2;
  }

  sf_mex_destroy(&c13_mxArrayInData);
}

static const mxArray *c13_l_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i240;
  int32_T c13_i241;
  int32_T c13_i242;
  real_T c13_b_inData[55];
  int32_T c13_i243;
  int32_T c13_i244;
  int32_T c13_i245;
  real_T c13_u[55];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_i240 = 0;
  for (c13_i241 = 0; c13_i241 < 11; c13_i241++) {
    for (c13_i242 = 0; c13_i242 < 5; c13_i242++) {
      c13_b_inData[c13_i242 + c13_i240] = (*(real_T (*)[55])c13_inData)[c13_i242
        + c13_i240];
    }

    c13_i240 += 5;
  }

  c13_i243 = 0;
  for (c13_i244 = 0; c13_i244 < 11; c13_i244++) {
    for (c13_i245 = 0; c13_i245 < 5; c13_i245++) {
      c13_u[c13_i245 + c13_i243] = c13_b_inData[c13_i245 + c13_i243];
    }

    c13_i243 += 5;
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 2, 5, 11),
                FALSE);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static void c13_r_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, real_T c13_y[55])
{
  real_T c13_dv18[55];
  int32_T c13_i246;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), c13_dv18, 1, 0, 0U, 1, 0U, 2, 5,
                11);
  for (c13_i246 = 0; c13_i246 < 55; c13_i246++) {
    c13_y[c13_i246] = c13_dv18[c13_i246];
  }

  sf_mex_destroy(&c13_u);
}

static void c13_l_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_chi_m;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y[55];
  int32_T c13_i247;
  int32_T c13_i248;
  int32_T c13_i249;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_chi_m = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_r_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_chi_m), &c13_thisId,
    c13_y);
  sf_mex_destroy(&c13_chi_m);
  c13_i247 = 0;
  for (c13_i248 = 0; c13_i248 < 11; c13_i248++) {
    for (c13_i249 = 0; c13_i249 < 5; c13_i249++) {
      (*(real_T (*)[55])c13_outData)[c13_i249 + c13_i247] = c13_y[c13_i249 +
        c13_i247];
    }

    c13_i247 += 5;
  }

  sf_mex_destroy(&c13_mxArrayInData);
}

static const mxArray *c13_m_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i250;
  int32_T c13_i251;
  int32_T c13_i252;
  real_T c13_b_inData[15];
  int32_T c13_i253;
  int32_T c13_i254;
  int32_T c13_i255;
  real_T c13_u[15];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_i250 = 0;
  for (c13_i251 = 0; c13_i251 < 3; c13_i251++) {
    for (c13_i252 = 0; c13_i252 < 5; c13_i252++) {
      c13_b_inData[c13_i252 + c13_i250] = (*(real_T (*)[15])c13_inData)[c13_i252
        + c13_i250];
    }

    c13_i250 += 5;
  }

  c13_i253 = 0;
  for (c13_i254 = 0; c13_i254 < 3; c13_i254++) {
    for (c13_i255 = 0; c13_i255 < 5; c13_i255++) {
      c13_u[c13_i255 + c13_i253] = c13_b_inData[c13_i255 + c13_i253];
    }

    c13_i253 += 5;
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 2, 5, 3), FALSE);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static const mxArray *c13_n_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i256;
  int32_T c13_i257;
  int32_T c13_i258;
  real_T c13_b_inData[9];
  int32_T c13_i259;
  int32_T c13_i260;
  int32_T c13_i261;
  real_T c13_u[9];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_i256 = 0;
  for (c13_i257 = 0; c13_i257 < 3; c13_i257++) {
    for (c13_i258 = 0; c13_i258 < 3; c13_i258++) {
      c13_b_inData[c13_i258 + c13_i256] = (*(real_T (*)[9])c13_inData)[c13_i258
        + c13_i256];
    }

    c13_i256 += 3;
  }

  c13_i259 = 0;
  for (c13_i260 = 0; c13_i260 < 3; c13_i260++) {
    for (c13_i261 = 0; c13_i261 < 3; c13_i261++) {
      c13_u[c13_i261 + c13_i259] = c13_b_inData[c13_i261 + c13_i259];
    }

    c13_i259 += 3;
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 2, 3, 3), FALSE);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static const mxArray *c13_o_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i262;
  int32_T c13_i263;
  int32_T c13_i264;
  real_T c13_b_inData[10];
  int32_T c13_i265;
  int32_T c13_i266;
  int32_T c13_i267;
  real_T c13_u[10];
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_i262 = 0;
  for (c13_i263 = 0; c13_i263 < 5; c13_i263++) {
    for (c13_i264 = 0; c13_i264 < 2; c13_i264++) {
      c13_b_inData[c13_i264 + c13_i262] = (*(real_T (*)[10])c13_inData)[c13_i264
        + c13_i262];
    }

    c13_i262 += 2;
  }

  c13_i265 = 0;
  for (c13_i266 = 0; c13_i266 < 5; c13_i266++) {
    for (c13_i267 = 0; c13_i267 < 2; c13_i267++) {
      c13_u[c13_i267 + c13_i265] = c13_b_inData[c13_i267 + c13_i265];
    }

    c13_i265 += 2;
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 2, 2, 5), FALSE);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

const mxArray *sf_c13_IMM_UKF_get_eml_resolved_functions_info(void)
{
  const mxArray *c13_nameCaptureInfo = NULL;
  c13_nameCaptureInfo = NULL;
  sf_mex_assign(&c13_nameCaptureInfo, sf_mex_createstruct("structure", 2, 161, 1),
                FALSE);
  c13_info_helper(&c13_nameCaptureInfo);
  c13_b_info_helper(&c13_nameCaptureInfo);
  c13_c_info_helper(&c13_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c13_nameCaptureInfo);
  return c13_nameCaptureInfo;
}

static void c13_info_helper(const mxArray **c13_info)
{
  const mxArray *c13_rhs0 = NULL;
  const mxArray *c13_lhs0 = NULL;
  const mxArray *c13_rhs1 = NULL;
  const mxArray *c13_lhs1 = NULL;
  const mxArray *c13_rhs2 = NULL;
  const mxArray *c13_lhs2 = NULL;
  const mxArray *c13_rhs3 = NULL;
  const mxArray *c13_lhs3 = NULL;
  const mxArray *c13_rhs4 = NULL;
  const mxArray *c13_lhs4 = NULL;
  const mxArray *c13_rhs5 = NULL;
  const mxArray *c13_lhs5 = NULL;
  const mxArray *c13_rhs6 = NULL;
  const mxArray *c13_lhs6 = NULL;
  const mxArray *c13_rhs7 = NULL;
  const mxArray *c13_lhs7 = NULL;
  const mxArray *c13_rhs8 = NULL;
  const mxArray *c13_lhs8 = NULL;
  const mxArray *c13_rhs9 = NULL;
  const mxArray *c13_lhs9 = NULL;
  const mxArray *c13_rhs10 = NULL;
  const mxArray *c13_lhs10 = NULL;
  const mxArray *c13_rhs11 = NULL;
  const mxArray *c13_lhs11 = NULL;
  const mxArray *c13_rhs12 = NULL;
  const mxArray *c13_lhs12 = NULL;
  const mxArray *c13_rhs13 = NULL;
  const mxArray *c13_lhs13 = NULL;
  const mxArray *c13_rhs14 = NULL;
  const mxArray *c13_lhs14 = NULL;
  const mxArray *c13_rhs15 = NULL;
  const mxArray *c13_lhs15 = NULL;
  const mxArray *c13_rhs16 = NULL;
  const mxArray *c13_lhs16 = NULL;
  const mxArray *c13_rhs17 = NULL;
  const mxArray *c13_lhs17 = NULL;
  const mxArray *c13_rhs18 = NULL;
  const mxArray *c13_lhs18 = NULL;
  const mxArray *c13_rhs19 = NULL;
  const mxArray *c13_lhs19 = NULL;
  const mxArray *c13_rhs20 = NULL;
  const mxArray *c13_lhs20 = NULL;
  const mxArray *c13_rhs21 = NULL;
  const mxArray *c13_lhs21 = NULL;
  const mxArray *c13_rhs22 = NULL;
  const mxArray *c13_lhs22 = NULL;
  const mxArray *c13_rhs23 = NULL;
  const mxArray *c13_lhs23 = NULL;
  const mxArray *c13_rhs24 = NULL;
  const mxArray *c13_lhs24 = NULL;
  const mxArray *c13_rhs25 = NULL;
  const mxArray *c13_lhs25 = NULL;
  const mxArray *c13_rhs26 = NULL;
  const mxArray *c13_lhs26 = NULL;
  const mxArray *c13_rhs27 = NULL;
  const mxArray *c13_lhs27 = NULL;
  const mxArray *c13_rhs28 = NULL;
  const mxArray *c13_lhs28 = NULL;
  const mxArray *c13_rhs29 = NULL;
  const mxArray *c13_lhs29 = NULL;
  const mxArray *c13_rhs30 = NULL;
  const mxArray *c13_lhs30 = NULL;
  const mxArray *c13_rhs31 = NULL;
  const mxArray *c13_lhs31 = NULL;
  const mxArray *c13_rhs32 = NULL;
  const mxArray *c13_lhs32 = NULL;
  const mxArray *c13_rhs33 = NULL;
  const mxArray *c13_lhs33 = NULL;
  const mxArray *c13_rhs34 = NULL;
  const mxArray *c13_lhs34 = NULL;
  const mxArray *c13_rhs35 = NULL;
  const mxArray *c13_lhs35 = NULL;
  const mxArray *c13_rhs36 = NULL;
  const mxArray *c13_lhs36 = NULL;
  const mxArray *c13_rhs37 = NULL;
  const mxArray *c13_lhs37 = NULL;
  const mxArray *c13_rhs38 = NULL;
  const mxArray *c13_lhs38 = NULL;
  const mxArray *c13_rhs39 = NULL;
  const mxArray *c13_lhs39 = NULL;
  const mxArray *c13_rhs40 = NULL;
  const mxArray *c13_lhs40 = NULL;
  const mxArray *c13_rhs41 = NULL;
  const mxArray *c13_lhs41 = NULL;
  const mxArray *c13_rhs42 = NULL;
  const mxArray *c13_lhs42 = NULL;
  const mxArray *c13_rhs43 = NULL;
  const mxArray *c13_lhs43 = NULL;
  const mxArray *c13_rhs44 = NULL;
  const mxArray *c13_lhs44 = NULL;
  const mxArray *c13_rhs45 = NULL;
  const mxArray *c13_lhs45 = NULL;
  const mxArray *c13_rhs46 = NULL;
  const mxArray *c13_lhs46 = NULL;
  const mxArray *c13_rhs47 = NULL;
  const mxArray *c13_lhs47 = NULL;
  const mxArray *c13_rhs48 = NULL;
  const mxArray *c13_lhs48 = NULL;
  const mxArray *c13_rhs49 = NULL;
  const mxArray *c13_lhs49 = NULL;
  const mxArray *c13_rhs50 = NULL;
  const mxArray *c13_lhs50 = NULL;
  const mxArray *c13_rhs51 = NULL;
  const mxArray *c13_lhs51 = NULL;
  const mxArray *c13_rhs52 = NULL;
  const mxArray *c13_lhs52 = NULL;
  const mxArray *c13_rhs53 = NULL;
  const mxArray *c13_lhs53 = NULL;
  const mxArray *c13_rhs54 = NULL;
  const mxArray *c13_lhs54 = NULL;
  const mxArray *c13_rhs55 = NULL;
  const mxArray *c13_lhs55 = NULL;
  const mxArray *c13_rhs56 = NULL;
  const mxArray *c13_lhs56 = NULL;
  const mxArray *c13_rhs57 = NULL;
  const mxArray *c13_lhs57 = NULL;
  const mxArray *c13_rhs58 = NULL;
  const mxArray *c13_lhs58 = NULL;
  const mxArray *c13_rhs59 = NULL;
  const mxArray *c13_lhs59 = NULL;
  const mxArray *c13_rhs60 = NULL;
  const mxArray *c13_lhs60 = NULL;
  const mxArray *c13_rhs61 = NULL;
  const mxArray *c13_lhs61 = NULL;
  const mxArray *c13_rhs62 = NULL;
  const mxArray *c13_lhs62 = NULL;
  const mxArray *c13_rhs63 = NULL;
  const mxArray *c13_lhs63 = NULL;
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mpower"), "name", "name", 0);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c13_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c13_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("ismatrix"), "name", "name",
                  2);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1331337258U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c13_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("power"), "name", "name", 3);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742680U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c13_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c13_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 5);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 5);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c13_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 6);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 6);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1358218540U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c13_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 7);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("floor"), "name", "name", 7);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742654U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c13_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 8);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c13_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 9);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851126U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c13_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 10);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 10);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c13_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 11);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mtimes"), "name", "name", 11);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c13_rhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 12);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 12);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c13_rhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "context", "context", 13);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mtimes"), "name", "name", 13);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c13_rhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "context", "context", 14);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mrdivide"), "name", "name",
                  14);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1373338908U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1319762366U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c13_rhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 15);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("rdivide"), "name", "name",
                  15);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742680U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c13_rhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 16);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c13_rhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 17);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c13_rhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_div"), "name", "name",
                  18);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742666U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c13_rhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "context", "context", 19);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("sqrt"), "name", "name", 19);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1343862786U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c13_rhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_error"), "name", "name",
                  20);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 20);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1343862758U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c13_rhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 21);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 21);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851138U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c13_rhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "context", "context", 22);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eye"), "name", "name", 22);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 22);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1368215430U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c13_rhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 23);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 23);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1368215430U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c13_rhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 24);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 24);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c13_rhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 25);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("isinf"), "name", "name", 25);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742656U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c13_rhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 26);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 26);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c13_rhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 27);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_is_integer_class"),
                  "name", "name", 27);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851182U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c13_rhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 28);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("intmax"), "name", "name", 28);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1362294282U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c13_rhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 29);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("intmin"), "name", "name", 29);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1362294282U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c13_rhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 30);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.indexIntRelop"), "name", "name", 30);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1326760722U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c13_rhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 31);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 31);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1326760396U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c13_rhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 32);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("intmin"), "name", "name", 32);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 32);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1362294282U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c13_rhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!numel_for_size"),
                  "context", "context", 33);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mtimes"), "name", "name", 33);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 33);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c13_rhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 34);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 34);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c13_rhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 35);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("intmax"), "name", "name", 35);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 35);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1362294282U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c13_rhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 36);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 36);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1346542740U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c13_rhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 37);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("intmax"), "name", "name", 37);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 37);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1362294282U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c13_rhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 38);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 38);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c13_rhs38, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs38, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 39);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 39);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 39);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c13_rhs39, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs39, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 40);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  40);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742670U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c13_rhs40, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs40, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 41);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 41);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1299109168U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c13_rhs41, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs41, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold"),
                  "context", "context", 42);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mtimes"), "name", "name", 42);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 42);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c13_rhs42, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs42, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 43);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 43);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c13_rhs43, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs43, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 44);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 44);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 44);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c13_rhs44, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs44, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 45);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_refblas_xgemm"), "name",
                  "name", 45);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1360314750U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c13_rhs45, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs45, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 46);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("min"), "name", "name", 46);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 46);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1311287718U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c13_rhs46, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs46, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "context",
                  "context", 47);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 47);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 47);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742670U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c13_rhs47, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs47, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 48);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 48);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 48);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c13_rhs48, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs48, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 49);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 49);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1358218540U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c13_rhs49, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs49, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 50);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 50);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c13_rhs50, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs50, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 51);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 51);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 51);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c13_rhs51, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs51, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 52);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 52);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c13_rhs52, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs52, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "context", "context", 53);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("chol"), "name", "name", 53);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m"), "resolved",
                  "resolved", 53);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1344504434U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c13_rhs53, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs53, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m"), "context",
                  "context", 54);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_tolower"), "name",
                  "name", 54);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_tolower.m"), "resolved",
                  "resolved", 54);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742670U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c13_rhs54, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs54, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 55);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 55);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c13_rhs55, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs55, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 56);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("ismatrix"), "name", "name",
                  56);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 56);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1331337258U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c13_rhs56, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs56, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 57);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 57);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 57);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1346542740U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c13_rhs57, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs57, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 58);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_error"), "name", "name",
                  58);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 58);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1343862758U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c13_rhs58, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs58, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 59);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_xpotrf"), "name", "name",
                  59);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xpotrf.m"),
                  "resolved", "resolved", 59);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851208U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c13_rhs59, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs59, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xpotrf.m"),
                  "context", "context", 60);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_lapack_xpotrf"), "name",
                  "name", 60);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 60);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xpotrf.m"),
                  "resolved", "resolved", 60);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851212U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c13_rhs60, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs60, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xpotrf.m"),
                  "context", "context", 61);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_matlab_zpotrf"), "name",
                  "name", 61);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "resolved", "resolved", 61);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851224U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c13_rhs61, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs61, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 62);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 62);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 62);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c13_rhs62, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs62, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 63);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 63);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 63);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c13_rhs63, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs63, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c13_rhs0);
  sf_mex_destroy(&c13_lhs0);
  sf_mex_destroy(&c13_rhs1);
  sf_mex_destroy(&c13_lhs1);
  sf_mex_destroy(&c13_rhs2);
  sf_mex_destroy(&c13_lhs2);
  sf_mex_destroy(&c13_rhs3);
  sf_mex_destroy(&c13_lhs3);
  sf_mex_destroy(&c13_rhs4);
  sf_mex_destroy(&c13_lhs4);
  sf_mex_destroy(&c13_rhs5);
  sf_mex_destroy(&c13_lhs5);
  sf_mex_destroy(&c13_rhs6);
  sf_mex_destroy(&c13_lhs6);
  sf_mex_destroy(&c13_rhs7);
  sf_mex_destroy(&c13_lhs7);
  sf_mex_destroy(&c13_rhs8);
  sf_mex_destroy(&c13_lhs8);
  sf_mex_destroy(&c13_rhs9);
  sf_mex_destroy(&c13_lhs9);
  sf_mex_destroy(&c13_rhs10);
  sf_mex_destroy(&c13_lhs10);
  sf_mex_destroy(&c13_rhs11);
  sf_mex_destroy(&c13_lhs11);
  sf_mex_destroy(&c13_rhs12);
  sf_mex_destroy(&c13_lhs12);
  sf_mex_destroy(&c13_rhs13);
  sf_mex_destroy(&c13_lhs13);
  sf_mex_destroy(&c13_rhs14);
  sf_mex_destroy(&c13_lhs14);
  sf_mex_destroy(&c13_rhs15);
  sf_mex_destroy(&c13_lhs15);
  sf_mex_destroy(&c13_rhs16);
  sf_mex_destroy(&c13_lhs16);
  sf_mex_destroy(&c13_rhs17);
  sf_mex_destroy(&c13_lhs17);
  sf_mex_destroy(&c13_rhs18);
  sf_mex_destroy(&c13_lhs18);
  sf_mex_destroy(&c13_rhs19);
  sf_mex_destroy(&c13_lhs19);
  sf_mex_destroy(&c13_rhs20);
  sf_mex_destroy(&c13_lhs20);
  sf_mex_destroy(&c13_rhs21);
  sf_mex_destroy(&c13_lhs21);
  sf_mex_destroy(&c13_rhs22);
  sf_mex_destroy(&c13_lhs22);
  sf_mex_destroy(&c13_rhs23);
  sf_mex_destroy(&c13_lhs23);
  sf_mex_destroy(&c13_rhs24);
  sf_mex_destroy(&c13_lhs24);
  sf_mex_destroy(&c13_rhs25);
  sf_mex_destroy(&c13_lhs25);
  sf_mex_destroy(&c13_rhs26);
  sf_mex_destroy(&c13_lhs26);
  sf_mex_destroy(&c13_rhs27);
  sf_mex_destroy(&c13_lhs27);
  sf_mex_destroy(&c13_rhs28);
  sf_mex_destroy(&c13_lhs28);
  sf_mex_destroy(&c13_rhs29);
  sf_mex_destroy(&c13_lhs29);
  sf_mex_destroy(&c13_rhs30);
  sf_mex_destroy(&c13_lhs30);
  sf_mex_destroy(&c13_rhs31);
  sf_mex_destroy(&c13_lhs31);
  sf_mex_destroy(&c13_rhs32);
  sf_mex_destroy(&c13_lhs32);
  sf_mex_destroy(&c13_rhs33);
  sf_mex_destroy(&c13_lhs33);
  sf_mex_destroy(&c13_rhs34);
  sf_mex_destroy(&c13_lhs34);
  sf_mex_destroy(&c13_rhs35);
  sf_mex_destroy(&c13_lhs35);
  sf_mex_destroy(&c13_rhs36);
  sf_mex_destroy(&c13_lhs36);
  sf_mex_destroy(&c13_rhs37);
  sf_mex_destroy(&c13_lhs37);
  sf_mex_destroy(&c13_rhs38);
  sf_mex_destroy(&c13_lhs38);
  sf_mex_destroy(&c13_rhs39);
  sf_mex_destroy(&c13_lhs39);
  sf_mex_destroy(&c13_rhs40);
  sf_mex_destroy(&c13_lhs40);
  sf_mex_destroy(&c13_rhs41);
  sf_mex_destroy(&c13_lhs41);
  sf_mex_destroy(&c13_rhs42);
  sf_mex_destroy(&c13_lhs42);
  sf_mex_destroy(&c13_rhs43);
  sf_mex_destroy(&c13_lhs43);
  sf_mex_destroy(&c13_rhs44);
  sf_mex_destroy(&c13_lhs44);
  sf_mex_destroy(&c13_rhs45);
  sf_mex_destroy(&c13_lhs45);
  sf_mex_destroy(&c13_rhs46);
  sf_mex_destroy(&c13_lhs46);
  sf_mex_destroy(&c13_rhs47);
  sf_mex_destroy(&c13_lhs47);
  sf_mex_destroy(&c13_rhs48);
  sf_mex_destroy(&c13_lhs48);
  sf_mex_destroy(&c13_rhs49);
  sf_mex_destroy(&c13_lhs49);
  sf_mex_destroy(&c13_rhs50);
  sf_mex_destroy(&c13_lhs50);
  sf_mex_destroy(&c13_rhs51);
  sf_mex_destroy(&c13_lhs51);
  sf_mex_destroy(&c13_rhs52);
  sf_mex_destroy(&c13_lhs52);
  sf_mex_destroy(&c13_rhs53);
  sf_mex_destroy(&c13_lhs53);
  sf_mex_destroy(&c13_rhs54);
  sf_mex_destroy(&c13_lhs54);
  sf_mex_destroy(&c13_rhs55);
  sf_mex_destroy(&c13_lhs55);
  sf_mex_destroy(&c13_rhs56);
  sf_mex_destroy(&c13_lhs56);
  sf_mex_destroy(&c13_rhs57);
  sf_mex_destroy(&c13_lhs57);
  sf_mex_destroy(&c13_rhs58);
  sf_mex_destroy(&c13_lhs58);
  sf_mex_destroy(&c13_rhs59);
  sf_mex_destroy(&c13_lhs59);
  sf_mex_destroy(&c13_rhs60);
  sf_mex_destroy(&c13_lhs60);
  sf_mex_destroy(&c13_rhs61);
  sf_mex_destroy(&c13_lhs61);
  sf_mex_destroy(&c13_rhs62);
  sf_mex_destroy(&c13_lhs62);
  sf_mex_destroy(&c13_rhs63);
  sf_mex_destroy(&c13_lhs63);
}

static const mxArray *c13_emlrt_marshallOut(char * c13_u)
{
  const mxArray *c13_y = NULL;
  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c13_u)), FALSE);
  return c13_y;
}

static const mxArray *c13_b_emlrt_marshallOut(uint32_T c13_u)
{
  const mxArray *c13_y = NULL;
  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", &c13_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c13_y;
}

static void c13_b_info_helper(const mxArray **c13_info)
{
  const mxArray *c13_rhs64 = NULL;
  const mxArray *c13_lhs64 = NULL;
  const mxArray *c13_rhs65 = NULL;
  const mxArray *c13_lhs65 = NULL;
  const mxArray *c13_rhs66 = NULL;
  const mxArray *c13_lhs66 = NULL;
  const mxArray *c13_rhs67 = NULL;
  const mxArray *c13_lhs67 = NULL;
  const mxArray *c13_rhs68 = NULL;
  const mxArray *c13_lhs68 = NULL;
  const mxArray *c13_rhs69 = NULL;
  const mxArray *c13_lhs69 = NULL;
  const mxArray *c13_rhs70 = NULL;
  const mxArray *c13_lhs70 = NULL;
  const mxArray *c13_rhs71 = NULL;
  const mxArray *c13_lhs71 = NULL;
  const mxArray *c13_rhs72 = NULL;
  const mxArray *c13_lhs72 = NULL;
  const mxArray *c13_rhs73 = NULL;
  const mxArray *c13_lhs73 = NULL;
  const mxArray *c13_rhs74 = NULL;
  const mxArray *c13_lhs74 = NULL;
  const mxArray *c13_rhs75 = NULL;
  const mxArray *c13_lhs75 = NULL;
  const mxArray *c13_rhs76 = NULL;
  const mxArray *c13_lhs76 = NULL;
  const mxArray *c13_rhs77 = NULL;
  const mxArray *c13_lhs77 = NULL;
  const mxArray *c13_rhs78 = NULL;
  const mxArray *c13_lhs78 = NULL;
  const mxArray *c13_rhs79 = NULL;
  const mxArray *c13_lhs79 = NULL;
  const mxArray *c13_rhs80 = NULL;
  const mxArray *c13_lhs80 = NULL;
  const mxArray *c13_rhs81 = NULL;
  const mxArray *c13_lhs81 = NULL;
  const mxArray *c13_rhs82 = NULL;
  const mxArray *c13_lhs82 = NULL;
  const mxArray *c13_rhs83 = NULL;
  const mxArray *c13_lhs83 = NULL;
  const mxArray *c13_rhs84 = NULL;
  const mxArray *c13_lhs84 = NULL;
  const mxArray *c13_rhs85 = NULL;
  const mxArray *c13_lhs85 = NULL;
  const mxArray *c13_rhs86 = NULL;
  const mxArray *c13_lhs86 = NULL;
  const mxArray *c13_rhs87 = NULL;
  const mxArray *c13_lhs87 = NULL;
  const mxArray *c13_rhs88 = NULL;
  const mxArray *c13_lhs88 = NULL;
  const mxArray *c13_rhs89 = NULL;
  const mxArray *c13_lhs89 = NULL;
  const mxArray *c13_rhs90 = NULL;
  const mxArray *c13_lhs90 = NULL;
  const mxArray *c13_rhs91 = NULL;
  const mxArray *c13_lhs91 = NULL;
  const mxArray *c13_rhs92 = NULL;
  const mxArray *c13_lhs92 = NULL;
  const mxArray *c13_rhs93 = NULL;
  const mxArray *c13_lhs93 = NULL;
  const mxArray *c13_rhs94 = NULL;
  const mxArray *c13_lhs94 = NULL;
  const mxArray *c13_rhs95 = NULL;
  const mxArray *c13_lhs95 = NULL;
  const mxArray *c13_rhs96 = NULL;
  const mxArray *c13_lhs96 = NULL;
  const mxArray *c13_rhs97 = NULL;
  const mxArray *c13_lhs97 = NULL;
  const mxArray *c13_rhs98 = NULL;
  const mxArray *c13_lhs98 = NULL;
  const mxArray *c13_rhs99 = NULL;
  const mxArray *c13_lhs99 = NULL;
  const mxArray *c13_rhs100 = NULL;
  const mxArray *c13_lhs100 = NULL;
  const mxArray *c13_rhs101 = NULL;
  const mxArray *c13_lhs101 = NULL;
  const mxArray *c13_rhs102 = NULL;
  const mxArray *c13_lhs102 = NULL;
  const mxArray *c13_rhs103 = NULL;
  const mxArray *c13_lhs103 = NULL;
  const mxArray *c13_rhs104 = NULL;
  const mxArray *c13_lhs104 = NULL;
  const mxArray *c13_rhs105 = NULL;
  const mxArray *c13_lhs105 = NULL;
  const mxArray *c13_rhs106 = NULL;
  const mxArray *c13_lhs106 = NULL;
  const mxArray *c13_rhs107 = NULL;
  const mxArray *c13_lhs107 = NULL;
  const mxArray *c13_rhs108 = NULL;
  const mxArray *c13_lhs108 = NULL;
  const mxArray *c13_rhs109 = NULL;
  const mxArray *c13_lhs109 = NULL;
  const mxArray *c13_rhs110 = NULL;
  const mxArray *c13_lhs110 = NULL;
  const mxArray *c13_rhs111 = NULL;
  const mxArray *c13_lhs111 = NULL;
  const mxArray *c13_rhs112 = NULL;
  const mxArray *c13_lhs112 = NULL;
  const mxArray *c13_rhs113 = NULL;
  const mxArray *c13_lhs113 = NULL;
  const mxArray *c13_rhs114 = NULL;
  const mxArray *c13_lhs114 = NULL;
  const mxArray *c13_rhs115 = NULL;
  const mxArray *c13_lhs115 = NULL;
  const mxArray *c13_rhs116 = NULL;
  const mxArray *c13_lhs116 = NULL;
  const mxArray *c13_rhs117 = NULL;
  const mxArray *c13_lhs117 = NULL;
  const mxArray *c13_rhs118 = NULL;
  const mxArray *c13_lhs118 = NULL;
  const mxArray *c13_rhs119 = NULL;
  const mxArray *c13_lhs119 = NULL;
  const mxArray *c13_rhs120 = NULL;
  const mxArray *c13_lhs120 = NULL;
  const mxArray *c13_rhs121 = NULL;
  const mxArray *c13_lhs121 = NULL;
  const mxArray *c13_rhs122 = NULL;
  const mxArray *c13_lhs122 = NULL;
  const mxArray *c13_rhs123 = NULL;
  const mxArray *c13_lhs123 = NULL;
  const mxArray *c13_rhs124 = NULL;
  const mxArray *c13_lhs124 = NULL;
  const mxArray *c13_rhs125 = NULL;
  const mxArray *c13_lhs125 = NULL;
  const mxArray *c13_rhs126 = NULL;
  const mxArray *c13_lhs126 = NULL;
  const mxArray *c13_rhs127 = NULL;
  const mxArray *c13_lhs127 = NULL;
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 64);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 64);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1346542740U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c13_rhs64, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs64, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 65);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 65);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 65);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c13_rhs65, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs65, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 66);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 66);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 66);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c13_rhs66, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs66, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 67);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 67);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 67);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 67);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851180U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c13_rhs67, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs67, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 68);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 68);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 68);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 68);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c13_rhs68, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs68, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 69);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 69);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 69);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 69);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c13_rhs69, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs69, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 70);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 70);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 70);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c13_rhs70, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs70, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 71);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_xdotc"), "name", "name",
                  71);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"),
                  "resolved", "resolved", 71);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742670U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c13_rhs71, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs71, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"), "context",
                  "context", 72);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 72);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 72);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 72);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1299109168U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c13_rhs72, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs72, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"), "context",
                  "context", 73);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_xdot"), "name", "name",
                  73);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 73);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m"), "resolved",
                  "resolved", 73);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742668U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c13_rhs73, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs73, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m"), "context",
                  "context", 74);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 74);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 74);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 74);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1299109168U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c13_rhs74, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs74, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs74), "lhs", "lhs",
                  74);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m!below_threshold"),
                  "context", "context", 75);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("length"), "name", "name", 75);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 75);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 75);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1303178606U), "fileTimeLo",
                  "fileTimeLo", 75);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 75);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 75);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 75);
  sf_mex_assign(&c13_rhs75, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs75, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs75), "rhs", "rhs",
                  75);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs75), "lhs", "lhs",
                  75);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m!intlength"),
                  "context", "context", 76);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 76);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 76);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 76);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 76);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 76);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 76);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 76);
  sf_mex_assign(&c13_rhs76, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs76, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs76), "rhs", "rhs",
                  76);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs76), "lhs", "lhs",
                  76);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"),
                  "context", "context", 77);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 77);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 77);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 77);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 77);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 77);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 77);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 77);
  sf_mex_assign(&c13_rhs77, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs77, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs77), "rhs", "rhs",
                  77);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs77), "lhs", "lhs",
                  77);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"),
                  "context", "context", 78);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_refblas_xdot"), "name",
                  "name", 78);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 78);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m"),
                  "resolved", "resolved", 78);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1299109172U), "fileTimeLo",
                  "fileTimeLo", 78);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 78);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 78);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 78);
  sf_mex_assign(&c13_rhs78, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs78, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs78), "rhs", "rhs",
                  78);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs78), "lhs", "lhs",
                  78);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m"),
                  "context", "context", 79);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_refblas_xdotx"), "name",
                  "name", 79);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 79);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "resolved", "resolved", 79);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1360314750U), "fileTimeLo",
                  "fileTimeLo", 79);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 79);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 79);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 79);
  sf_mex_assign(&c13_rhs79, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs79, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs79), "rhs", "rhs",
                  79);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs79), "lhs", "lhs",
                  79);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 80);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 80);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 80);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 80);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 80);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 80);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 80);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 80);
  sf_mex_assign(&c13_rhs80, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs80, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs80), "rhs", "rhs",
                  80);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs80), "lhs", "lhs",
                  80);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 81);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 81);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 81);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 81);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 81);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 81);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 81);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 81);
  sf_mex_assign(&c13_rhs81, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs81, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs81), "rhs", "rhs",
                  81);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs81), "lhs", "lhs",
                  81);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 82);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 82);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 82);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 82);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1346542740U), "fileTimeLo",
                  "fileTimeLo", 82);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 82);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 82);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 82);
  sf_mex_assign(&c13_rhs82, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs82, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs82), "rhs", "rhs",
                  82);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs82), "lhs", "lhs",
                  82);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 83);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 83);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 83);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 83);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 83);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 83);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 83);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 83);
  sf_mex_assign(&c13_rhs83, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs83, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs83), "rhs", "rhs",
                  83);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs83), "lhs", "lhs",
                  83);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 84);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 84);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 84);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 84);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 84);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 84);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 84);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 84);
  sf_mex_assign(&c13_rhs84, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs84, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs84), "rhs", "rhs",
                  84);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs84), "lhs", "lhs",
                  84);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 85);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 85);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 85);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 85);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 85);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 85);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 85);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 85);
  sf_mex_assign(&c13_rhs85, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs85, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs85), "rhs", "rhs",
                  85);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs85), "lhs", "lhs",
                  85);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 86);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_xgemv"), "name", "name",
                  86);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 86);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemv.m"),
                  "resolved", "resolved", 86);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742670U), "fileTimeLo",
                  "fileTimeLo", 86);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 86);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 86);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 86);
  sf_mex_assign(&c13_rhs86, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs86, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs86), "rhs", "rhs",
                  86);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs86), "lhs", "lhs",
                  86);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemv.m"), "context",
                  "context", 87);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 87);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 87);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 87);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1299109168U), "fileTimeLo",
                  "fileTimeLo", 87);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 87);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 87);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 87);
  sf_mex_assign(&c13_rhs87, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs87, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs87), "rhs", "rhs",
                  87);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs87), "lhs", "lhs",
                  87);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m!below_threshold"),
                  "context", "context", 88);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("length"), "name", "name", 88);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 88);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 88);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1303178606U), "fileTimeLo",
                  "fileTimeLo", 88);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 88);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 88);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 88);
  sf_mex_assign(&c13_rhs88, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs88, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs88), "rhs", "rhs",
                  88);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs88), "lhs", "lhs",
                  88);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m!below_threshold"),
                  "context", "context", 89);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("intmax"), "name", "name", 89);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 89);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 89);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1362294282U), "fileTimeLo",
                  "fileTimeLo", 89);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 89);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 89);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 89);
  sf_mex_assign(&c13_rhs89, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs89, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs89), "rhs", "rhs",
                  89);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs89), "lhs", "lhs",
                  89);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m!below_threshold"),
                  "context", "context", 90);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mtimes"), "name", "name", 90);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 90);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 90);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 90);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 90);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 90);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 90);
  sf_mex_assign(&c13_rhs90, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs90, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs90), "rhs", "rhs",
                  90);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs90), "lhs", "lhs",
                  90);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m"),
                  "context", "context", 91);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 91);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 91);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 91);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 91);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 91);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 91);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 91);
  sf_mex_assign(&c13_rhs91, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs91, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs91), "rhs", "rhs",
                  91);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs91), "lhs", "lhs",
                  91);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m"),
                  "context", "context", 92);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 92);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 92);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 92);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 92);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 92);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 92);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 92);
  sf_mex_assign(&c13_rhs92, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs92, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs92), "rhs", "rhs",
                  92);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs92), "lhs", "lhs",
                  92);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m"),
                  "context", "context", 93);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_refblas_xgemv"), "name",
                  "name", 93);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 93);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "resolved", "resolved", 93);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1360314752U), "fileTimeLo",
                  "fileTimeLo", 93);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 93);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 93);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 93);
  sf_mex_assign(&c13_rhs93, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs93, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs93), "rhs", "rhs",
                  93);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs93), "lhs", "lhs",
                  93);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 94);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 94);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 94);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 94);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 94);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 94);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 94);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 94);
  sf_mex_assign(&c13_rhs94, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs94, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs94), "rhs", "rhs",
                  94);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs94), "lhs", "lhs",
                  94);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 95);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 95);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 95);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 95);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 95);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 95);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 95);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 95);
  sf_mex_assign(&c13_rhs95, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs95, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs95), "rhs", "rhs",
                  95);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs95), "lhs", "lhs",
                  95);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 96);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 96);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 96);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 96);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851180U), "fileTimeLo",
                  "fileTimeLo", 96);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 96);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 96);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 96);
  sf_mex_assign(&c13_rhs96, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs96, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs96), "rhs", "rhs",
                  96);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs96), "lhs", "lhs",
                  96);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 97);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 97);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 97);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 97);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 97);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 97);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 97);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 97);
  sf_mex_assign(&c13_rhs97, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs97, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs97), "rhs", "rhs",
                  97);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs97), "lhs", "lhs",
                  97);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 98);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 98);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 98);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 98);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1346542740U), "fileTimeLo",
                  "fileTimeLo", 98);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 98);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 98);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 98);
  sf_mex_assign(&c13_rhs98, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs98, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs98), "rhs", "rhs",
                  98);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs98), "lhs", "lhs",
                  98);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 99);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_div"), "name", "name",
                  99);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 99);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 99);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742666U), "fileTimeLo",
                  "fileTimeLo", 99);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 99);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 99);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 99);
  sf_mex_assign(&c13_rhs99, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs99, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs99), "rhs", "rhs",
                  99);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs99), "lhs", "lhs",
                  99);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 100);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_xscal"), "name", "name",
                  100);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 100);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"),
                  "resolved", "resolved", 100);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742672U), "fileTimeLo",
                  "fileTimeLo", 100);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 100);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 100);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 100);
  sf_mex_assign(&c13_rhs100, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs100, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs100), "rhs",
                  "rhs", 100);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs100), "lhs",
                  "lhs", 100);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"), "context",
                  "context", 101);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 101);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 101);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 101);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1299109168U), "fileTimeLo",
                  "fileTimeLo", 101);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 101);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 101);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 101);
  sf_mex_assign(&c13_rhs101, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs101, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs101), "rhs",
                  "rhs", 101);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs101), "lhs",
                  "lhs", 101);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xscal.m!below_threshold"),
                  "context", "context", 102);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("length"), "name", "name",
                  102);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 102);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 102);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1303178606U), "fileTimeLo",
                  "fileTimeLo", 102);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 102);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 102);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 102);
  sf_mex_assign(&c13_rhs102, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs102, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs102), "rhs",
                  "rhs", 102);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs102), "lhs",
                  "lhs", 102);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xscal.m"),
                  "context", "context", 103);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 103);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 103);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 103);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 103);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 103);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 103);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 103);
  sf_mex_assign(&c13_rhs103, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs103, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs103), "rhs",
                  "rhs", 103);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs103), "lhs",
                  "lhs", 103);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xscal.m"),
                  "context", "context", 104);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 104);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 104);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 104);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 104);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 104);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 104);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 104);
  sf_mex_assign(&c13_rhs104, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs104, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs104), "rhs",
                  "rhs", 104);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs104), "lhs",
                  "lhs", 104);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xscal.m"),
                  "context", "context", 105);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_refblas_xscal"), "name",
                  "name", 105);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 105);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "resolved", "resolved", 105);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1299109184U), "fileTimeLo",
                  "fileTimeLo", 105);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 105);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 105);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 105);
  sf_mex_assign(&c13_rhs105, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs105, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs105), "rhs",
                  "rhs", 105);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs105), "lhs",
                  "lhs", 105);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 106);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 106);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 106);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 106);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 106);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 106);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 106);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 106);
  sf_mex_assign(&c13_rhs106, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs106, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs106), "rhs",
                  "rhs", 106);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs106), "lhs",
                  "lhs", 106);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 107);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 107);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 107);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 107);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 107);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 107);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 107);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 107);
  sf_mex_assign(&c13_rhs107, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs107, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs107), "rhs",
                  "rhs", 107);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs107), "lhs",
                  "lhs", 107);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 108);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 108);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 108);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 108);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851180U), "fileTimeLo",
                  "fileTimeLo", 108);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 108);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 108);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 108);
  sf_mex_assign(&c13_rhs108, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs108, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs108), "rhs",
                  "rhs", 108);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs108), "lhs",
                  "lhs", 108);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 109);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 109);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 109);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 109);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 109);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 109);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 109);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 109);
  sf_mex_assign(&c13_rhs109, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs109, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs109), "rhs",
                  "rhs", 109);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs109), "lhs",
                  "lhs", 109);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 110);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 110);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 110);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 110);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1346542740U), "fileTimeLo",
                  "fileTimeLo", 110);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 110);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 110);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 110);
  sf_mex_assign(&c13_rhs110, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs110, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs110), "rhs",
                  "rhs", 110);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs110), "lhs",
                  "lhs", 110);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 111);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 111);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 111);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 111);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 111);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 111);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 111);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 111);
  sf_mex_assign(&c13_rhs111, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs111, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs111), "rhs",
                  "rhs", 111);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs111), "lhs",
                  "lhs", 111);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "context", "context",
                  112);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("sin"), "name", "name", 112);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 112);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 112);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1343862786U), "fileTimeLo",
                  "fileTimeLo", 112);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 112);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 112);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 112);
  sf_mex_assign(&c13_rhs112, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs112, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs112), "rhs",
                  "rhs", 112);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs112), "lhs",
                  "lhs", 112);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 113);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 113);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 113);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 113);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851136U), "fileTimeLo",
                  "fileTimeLo", 113);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 113);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 113);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 113);
  sf_mex_assign(&c13_rhs113, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs113, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs113), "rhs",
                  "rhs", 113);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs113), "lhs",
                  "lhs", 113);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "context", "context",
                  114);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("cos"), "name", "name", 114);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 114);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 114);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1343862772U), "fileTimeLo",
                  "fileTimeLo", 114);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 114);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 114);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 114);
  sf_mex_assign(&c13_rhs114, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs114, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs114), "rhs",
                  "rhs", 114);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs114), "lhs",
                  "lhs", 114);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 115);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 115);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 115);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 115);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851122U), "fileTimeLo",
                  "fileTimeLo", 115);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 115);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 115);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 115);
  sf_mex_assign(&c13_rhs115, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs115, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs115), "rhs",
                  "rhs", 115);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs115), "lhs",
                  "lhs", 115);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "context", "context",
                  116);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("inv"), "name", "name", 116);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 116);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m"), "resolved",
                  "resolved", 116);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1305350400U), "fileTimeLo",
                  "fileTimeLo", 116);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 116);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 116);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 116);
  sf_mex_assign(&c13_rhs116, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs116, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs116), "rhs",
                  "rhs", 116);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs116), "lhs",
                  "lhs", 116);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv2x2"), "context",
                  "context", 117);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_xcabs1"), "name", "name",
                  117);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 117);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "resolved", "resolved", 117);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851106U), "fileTimeLo",
                  "fileTimeLo", 117);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 117);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 117);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 117);
  sf_mex_assign(&c13_rhs117, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs117, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs117), "rhs",
                  "rhs", 117);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs117), "lhs",
                  "lhs", 117);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "context", "context", 118);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("abs"), "name", "name", 118);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 118);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 118);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742652U), "fileTimeLo",
                  "fileTimeLo", 118);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 118);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 118);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 118);
  sf_mex_assign(&c13_rhs118, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs118, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs118), "rhs",
                  "rhs", 118);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs118), "lhs",
                  "lhs", 118);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 119);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 119);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 119);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 119);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 119);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 119);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 119);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 119);
  sf_mex_assign(&c13_rhs119, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs119, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs119), "rhs",
                  "rhs", 119);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs119), "lhs",
                  "lhs", 119);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 120);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 120);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 120);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 120);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851112U), "fileTimeLo",
                  "fileTimeLo", 120);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 120);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 120);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 120);
  sf_mex_assign(&c13_rhs120, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs120, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs120), "rhs",
                  "rhs", 120);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs120), "lhs",
                  "lhs", 120);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv2x2"), "context",
                  "context", 121);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mrdivide"), "name", "name",
                  121);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 121);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 121);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1373338908U), "fileTimeLo",
                  "fileTimeLo", 121);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 121);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1319762366U), "mFileTimeLo",
                  "mFileTimeLo", 121);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 121);
  sf_mex_assign(&c13_rhs121, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs121, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs121), "rhs",
                  "rhs", 121);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs121), "lhs",
                  "lhs", 121);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv2x2"), "context",
                  "context", 122);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mtimes"), "name", "name",
                  122);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 122);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 122);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 122);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 122);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 122);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 122);
  sf_mex_assign(&c13_rhs122, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs122, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs122), "rhs",
                  "rhs", 122);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs122), "lhs",
                  "lhs", 122);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 123);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("norm"), "name", "name", 123);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 123);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 123);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742668U), "fileTimeLo",
                  "fileTimeLo", 123);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 123);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 123);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 123);
  sf_mex_assign(&c13_rhs123, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs123, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs123), "rhs",
                  "rhs", 123);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs123), "lhs",
                  "lhs", 123);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "context",
                  "context", 124);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 124);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 124);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 124);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 124);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 124);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 124);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 124);
  sf_mex_assign(&c13_rhs124, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs124, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs124), "rhs",
                  "rhs", 124);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs124), "lhs",
                  "lhs", 124);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 125);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("abs"), "name", "name", 125);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 125);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 125);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742652U), "fileTimeLo",
                  "fileTimeLo", 125);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 125);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 125);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 125);
  sf_mex_assign(&c13_rhs125, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs125, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs125), "rhs",
                  "rhs", 125);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs125), "lhs",
                  "lhs", 125);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 126);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("isnan"), "name", "name", 126);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 126);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 126);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742658U), "fileTimeLo",
                  "fileTimeLo", 126);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 126);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 126);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 126);
  sf_mex_assign(&c13_rhs126, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs126, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs126), "rhs",
                  "rhs", 126);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs126), "lhs",
                  "lhs", 126);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 127);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 127);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 127);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 127);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 127);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 127);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 127);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 127);
  sf_mex_assign(&c13_rhs127, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs127, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs127), "rhs",
                  "rhs", 127);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs127), "lhs",
                  "lhs", 127);
  sf_mex_destroy(&c13_rhs64);
  sf_mex_destroy(&c13_lhs64);
  sf_mex_destroy(&c13_rhs65);
  sf_mex_destroy(&c13_lhs65);
  sf_mex_destroy(&c13_rhs66);
  sf_mex_destroy(&c13_lhs66);
  sf_mex_destroy(&c13_rhs67);
  sf_mex_destroy(&c13_lhs67);
  sf_mex_destroy(&c13_rhs68);
  sf_mex_destroy(&c13_lhs68);
  sf_mex_destroy(&c13_rhs69);
  sf_mex_destroy(&c13_lhs69);
  sf_mex_destroy(&c13_rhs70);
  sf_mex_destroy(&c13_lhs70);
  sf_mex_destroy(&c13_rhs71);
  sf_mex_destroy(&c13_lhs71);
  sf_mex_destroy(&c13_rhs72);
  sf_mex_destroy(&c13_lhs72);
  sf_mex_destroy(&c13_rhs73);
  sf_mex_destroy(&c13_lhs73);
  sf_mex_destroy(&c13_rhs74);
  sf_mex_destroy(&c13_lhs74);
  sf_mex_destroy(&c13_rhs75);
  sf_mex_destroy(&c13_lhs75);
  sf_mex_destroy(&c13_rhs76);
  sf_mex_destroy(&c13_lhs76);
  sf_mex_destroy(&c13_rhs77);
  sf_mex_destroy(&c13_lhs77);
  sf_mex_destroy(&c13_rhs78);
  sf_mex_destroy(&c13_lhs78);
  sf_mex_destroy(&c13_rhs79);
  sf_mex_destroy(&c13_lhs79);
  sf_mex_destroy(&c13_rhs80);
  sf_mex_destroy(&c13_lhs80);
  sf_mex_destroy(&c13_rhs81);
  sf_mex_destroy(&c13_lhs81);
  sf_mex_destroy(&c13_rhs82);
  sf_mex_destroy(&c13_lhs82);
  sf_mex_destroy(&c13_rhs83);
  sf_mex_destroy(&c13_lhs83);
  sf_mex_destroy(&c13_rhs84);
  sf_mex_destroy(&c13_lhs84);
  sf_mex_destroy(&c13_rhs85);
  sf_mex_destroy(&c13_lhs85);
  sf_mex_destroy(&c13_rhs86);
  sf_mex_destroy(&c13_lhs86);
  sf_mex_destroy(&c13_rhs87);
  sf_mex_destroy(&c13_lhs87);
  sf_mex_destroy(&c13_rhs88);
  sf_mex_destroy(&c13_lhs88);
  sf_mex_destroy(&c13_rhs89);
  sf_mex_destroy(&c13_lhs89);
  sf_mex_destroy(&c13_rhs90);
  sf_mex_destroy(&c13_lhs90);
  sf_mex_destroy(&c13_rhs91);
  sf_mex_destroy(&c13_lhs91);
  sf_mex_destroy(&c13_rhs92);
  sf_mex_destroy(&c13_lhs92);
  sf_mex_destroy(&c13_rhs93);
  sf_mex_destroy(&c13_lhs93);
  sf_mex_destroy(&c13_rhs94);
  sf_mex_destroy(&c13_lhs94);
  sf_mex_destroy(&c13_rhs95);
  sf_mex_destroy(&c13_lhs95);
  sf_mex_destroy(&c13_rhs96);
  sf_mex_destroy(&c13_lhs96);
  sf_mex_destroy(&c13_rhs97);
  sf_mex_destroy(&c13_lhs97);
  sf_mex_destroy(&c13_rhs98);
  sf_mex_destroy(&c13_lhs98);
  sf_mex_destroy(&c13_rhs99);
  sf_mex_destroy(&c13_lhs99);
  sf_mex_destroy(&c13_rhs100);
  sf_mex_destroy(&c13_lhs100);
  sf_mex_destroy(&c13_rhs101);
  sf_mex_destroy(&c13_lhs101);
  sf_mex_destroy(&c13_rhs102);
  sf_mex_destroy(&c13_lhs102);
  sf_mex_destroy(&c13_rhs103);
  sf_mex_destroy(&c13_lhs103);
  sf_mex_destroy(&c13_rhs104);
  sf_mex_destroy(&c13_lhs104);
  sf_mex_destroy(&c13_rhs105);
  sf_mex_destroy(&c13_lhs105);
  sf_mex_destroy(&c13_rhs106);
  sf_mex_destroy(&c13_lhs106);
  sf_mex_destroy(&c13_rhs107);
  sf_mex_destroy(&c13_lhs107);
  sf_mex_destroy(&c13_rhs108);
  sf_mex_destroy(&c13_lhs108);
  sf_mex_destroy(&c13_rhs109);
  sf_mex_destroy(&c13_lhs109);
  sf_mex_destroy(&c13_rhs110);
  sf_mex_destroy(&c13_lhs110);
  sf_mex_destroy(&c13_rhs111);
  sf_mex_destroy(&c13_lhs111);
  sf_mex_destroy(&c13_rhs112);
  sf_mex_destroy(&c13_lhs112);
  sf_mex_destroy(&c13_rhs113);
  sf_mex_destroy(&c13_lhs113);
  sf_mex_destroy(&c13_rhs114);
  sf_mex_destroy(&c13_lhs114);
  sf_mex_destroy(&c13_rhs115);
  sf_mex_destroy(&c13_lhs115);
  sf_mex_destroy(&c13_rhs116);
  sf_mex_destroy(&c13_lhs116);
  sf_mex_destroy(&c13_rhs117);
  sf_mex_destroy(&c13_lhs117);
  sf_mex_destroy(&c13_rhs118);
  sf_mex_destroy(&c13_lhs118);
  sf_mex_destroy(&c13_rhs119);
  sf_mex_destroy(&c13_lhs119);
  sf_mex_destroy(&c13_rhs120);
  sf_mex_destroy(&c13_lhs120);
  sf_mex_destroy(&c13_rhs121);
  sf_mex_destroy(&c13_lhs121);
  sf_mex_destroy(&c13_rhs122);
  sf_mex_destroy(&c13_lhs122);
  sf_mex_destroy(&c13_rhs123);
  sf_mex_destroy(&c13_lhs123);
  sf_mex_destroy(&c13_rhs124);
  sf_mex_destroy(&c13_lhs124);
  sf_mex_destroy(&c13_rhs125);
  sf_mex_destroy(&c13_lhs125);
  sf_mex_destroy(&c13_rhs126);
  sf_mex_destroy(&c13_lhs126);
  sf_mex_destroy(&c13_rhs127);
  sf_mex_destroy(&c13_lhs127);
}

static void c13_c_info_helper(const mxArray **c13_info)
{
  const mxArray *c13_rhs128 = NULL;
  const mxArray *c13_lhs128 = NULL;
  const mxArray *c13_rhs129 = NULL;
  const mxArray *c13_lhs129 = NULL;
  const mxArray *c13_rhs130 = NULL;
  const mxArray *c13_lhs130 = NULL;
  const mxArray *c13_rhs131 = NULL;
  const mxArray *c13_lhs131 = NULL;
  const mxArray *c13_rhs132 = NULL;
  const mxArray *c13_lhs132 = NULL;
  const mxArray *c13_rhs133 = NULL;
  const mxArray *c13_lhs133 = NULL;
  const mxArray *c13_rhs134 = NULL;
  const mxArray *c13_lhs134 = NULL;
  const mxArray *c13_rhs135 = NULL;
  const mxArray *c13_lhs135 = NULL;
  const mxArray *c13_rhs136 = NULL;
  const mxArray *c13_lhs136 = NULL;
  const mxArray *c13_rhs137 = NULL;
  const mxArray *c13_lhs137 = NULL;
  const mxArray *c13_rhs138 = NULL;
  const mxArray *c13_lhs138 = NULL;
  const mxArray *c13_rhs139 = NULL;
  const mxArray *c13_lhs139 = NULL;
  const mxArray *c13_rhs140 = NULL;
  const mxArray *c13_lhs140 = NULL;
  const mxArray *c13_rhs141 = NULL;
  const mxArray *c13_lhs141 = NULL;
  const mxArray *c13_rhs142 = NULL;
  const mxArray *c13_lhs142 = NULL;
  const mxArray *c13_rhs143 = NULL;
  const mxArray *c13_lhs143 = NULL;
  const mxArray *c13_rhs144 = NULL;
  const mxArray *c13_lhs144 = NULL;
  const mxArray *c13_rhs145 = NULL;
  const mxArray *c13_lhs145 = NULL;
  const mxArray *c13_rhs146 = NULL;
  const mxArray *c13_lhs146 = NULL;
  const mxArray *c13_rhs147 = NULL;
  const mxArray *c13_lhs147 = NULL;
  const mxArray *c13_rhs148 = NULL;
  const mxArray *c13_lhs148 = NULL;
  const mxArray *c13_rhs149 = NULL;
  const mxArray *c13_lhs149 = NULL;
  const mxArray *c13_rhs150 = NULL;
  const mxArray *c13_lhs150 = NULL;
  const mxArray *c13_rhs151 = NULL;
  const mxArray *c13_lhs151 = NULL;
  const mxArray *c13_rhs152 = NULL;
  const mxArray *c13_lhs152 = NULL;
  const mxArray *c13_rhs153 = NULL;
  const mxArray *c13_lhs153 = NULL;
  const mxArray *c13_rhs154 = NULL;
  const mxArray *c13_lhs154 = NULL;
  const mxArray *c13_rhs155 = NULL;
  const mxArray *c13_lhs155 = NULL;
  const mxArray *c13_rhs156 = NULL;
  const mxArray *c13_lhs156 = NULL;
  const mxArray *c13_rhs157 = NULL;
  const mxArray *c13_lhs157 = NULL;
  const mxArray *c13_rhs158 = NULL;
  const mxArray *c13_lhs158 = NULL;
  const mxArray *c13_rhs159 = NULL;
  const mxArray *c13_lhs159 = NULL;
  const mxArray *c13_rhs160 = NULL;
  const mxArray *c13_lhs160 = NULL;
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 128);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 128);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 128);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 128);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851176U), "fileTimeLo",
                  "fileTimeLo", 128);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 128);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 128);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 128);
  sf_mex_assign(&c13_rhs128, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs128, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs128), "rhs",
                  "rhs", 128);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs128), "lhs",
                  "lhs", 128);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "context", "context", 129);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 129);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 129);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 129);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851182U), "fileTimeLo",
                  "fileTimeLo", 129);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 129);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 129);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 129);
  sf_mex_assign(&c13_rhs129, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs129, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs129), "rhs",
                  "rhs", 129);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs129), "lhs",
                  "lhs", 129);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 130);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mtimes"), "name", "name",
                  130);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 130);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 130);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 130);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 130);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 130);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 130);
  sf_mex_assign(&c13_rhs130, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs130, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs130), "rhs",
                  "rhs", 130);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs130), "lhs",
                  "lhs", 130);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 131);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_warning"), "name",
                  "name", 131);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 131);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 131);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851202U), "fileTimeLo",
                  "fileTimeLo", 131);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 131);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 131);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 131);
  sf_mex_assign(&c13_rhs131, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs131, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs131), "rhs",
                  "rhs", 131);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs131), "lhs",
                  "lhs", 131);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 132);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("isnan"), "name", "name", 132);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 132);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 132);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742658U), "fileTimeLo",
                  "fileTimeLo", 132);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 132);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 132);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 132);
  sf_mex_assign(&c13_rhs132, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs132, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs132), "rhs",
                  "rhs", 132);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs132), "lhs",
                  "lhs", 132);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 133);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eps"), "name", "name", 133);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 133);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 133);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1326760396U), "fileTimeLo",
                  "fileTimeLo", 133);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 133);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 133);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 133);
  sf_mex_assign(&c13_rhs133, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs133, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs133), "rhs",
                  "rhs", 133);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs133), "lhs",
                  "lhs", 133);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 134);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 134);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 134);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 134);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851182U), "fileTimeLo",
                  "fileTimeLo", 134);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 134);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 134);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 134);
  sf_mex_assign(&c13_rhs134, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs134, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs134), "rhs",
                  "rhs", 134);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs134), "lhs",
                  "lhs", 134);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 135);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_eps"), "name", "name",
                  135);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 135);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 135);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1326760396U), "fileTimeLo",
                  "fileTimeLo", 135);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 135);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 135);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 135);
  sf_mex_assign(&c13_rhs135, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs135, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs135), "rhs",
                  "rhs", 135);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs135), "lhs",
                  "lhs", 135);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 136);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 136);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 136);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 136);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1326760396U), "fileTimeLo",
                  "fileTimeLo", 136);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 136);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 136);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 136);
  sf_mex_assign(&c13_rhs136, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs136, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs136), "rhs",
                  "rhs", 136);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs136), "lhs",
                  "lhs", 136);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 137);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_flt2str"), "name",
                  "name", 137);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 137);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "resolved",
                  "resolved", 137);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1360314750U), "fileTimeLo",
                  "fileTimeLo", 137);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 137);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 137);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 137);
  sf_mex_assign(&c13_rhs137, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs137, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs137), "rhs",
                  "rhs", 137);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs137), "lhs",
                  "lhs", 137);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "context",
                  "context", 138);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "name", "name", 138);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 138);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/char.m"), "resolved",
                  "resolved", 138);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1319762368U), "fileTimeLo",
                  "fileTimeLo", 138);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 138);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 138);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 138);
  sf_mex_assign(&c13_rhs138, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs138, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs138), "rhs",
                  "rhs", 138);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs138), "lhs",
                  "lhs", 138);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 139);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_xdotu"), "name", "name",
                  139);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 139);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m"),
                  "resolved", "resolved", 139);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742670U), "fileTimeLo",
                  "fileTimeLo", 139);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 139);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 139);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 139);
  sf_mex_assign(&c13_rhs139, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs139, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs139), "rhs",
                  "rhs", 139);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs139), "lhs",
                  "lhs", 139);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m"), "context",
                  "context", 140);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 140);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 140);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 140);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1299109168U), "fileTimeLo",
                  "fileTimeLo", 140);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 140);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 140);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 140);
  sf_mex_assign(&c13_rhs140, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs140, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs140), "rhs",
                  "rhs", 140);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs140), "lhs",
                  "lhs", 140);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m"), "context",
                  "context", 141);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_xdot"), "name", "name",
                  141);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 141);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m"), "resolved",
                  "resolved", 141);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742668U), "fileTimeLo",
                  "fileTimeLo", 141);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 141);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 141);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 141);
  sf_mex_assign(&c13_rhs141, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs141, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs141), "rhs",
                  "rhs", 141);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs141), "lhs",
                  "lhs", 141);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 142);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 142);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 142);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 142);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 142);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 142);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 142);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 142);
  sf_mex_assign(&c13_rhs142, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs142, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs142), "rhs",
                  "rhs", 142);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs142), "lhs",
                  "lhs", 142);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 143);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 143);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 143);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 143);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851180U), "fileTimeLo",
                  "fileTimeLo", 143);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 143);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 143);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 143);
  sf_mex_assign(&c13_rhs143, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs143, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs143), "rhs",
                  "rhs", 143);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs143), "lhs",
                  "lhs", 143);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "context", "context",
                  144);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("expm"), "name", "name", 144);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 144);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "resolved",
                  "resolved", 144);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851224U), "fileTimeLo",
                  "fileTimeLo", 144);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 144);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 144);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 144);
  sf_mex_assign(&c13_rhs144, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs144, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs144), "rhs",
                  "rhs", 144);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs144), "lhs",
                  "lhs", 144);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 145);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("norm"), "name", "name", 145);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 145);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 145);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742668U), "fileTimeLo",
                  "fileTimeLo", 145);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 145);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 145);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 145);
  sf_mex_assign(&c13_rhs145, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs145, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs145), "rhs",
                  "rhs", 145);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs145), "lhs",
                  "lhs", 145);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "context",
                  "context", 146);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("abs"), "name", "name", 146);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 146);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 146);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742652U), "fileTimeLo",
                  "fileTimeLo", 146);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 146);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 146);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 146);
  sf_mex_assign(&c13_rhs146, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs146, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs146), "rhs",
                  "rhs", 146);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs146), "lhs",
                  "lhs", 146);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m!PadeApproximantOfDegree"),
                  "context", "context", 147);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mtimes"), "name", "name",
                  147);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 147);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 147);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 147);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 147);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 147);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 147);
  sf_mex_assign(&c13_rhs147, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs147, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs147), "rhs",
                  "rhs", 147);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs147), "lhs",
                  "lhs", 147);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m!PadeApproximantOfDegree"),
                  "context", "context", 148);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mldivide"), "name", "name",
                  148);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 148);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "resolved",
                  "resolved", 148);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1373338908U), "fileTimeLo",
                  "fileTimeLo", 148);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 148);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1319762366U), "mFileTimeLo",
                  "mFileTimeLo", 148);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 148);
  sf_mex_assign(&c13_rhs148, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs148, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs148), "rhs",
                  "rhs", 148);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs148), "lhs",
                  "lhs", 148);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "context",
                  "context", 149);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("rdivide"), "name", "name",
                  149);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 149);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 149);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742680U), "fileTimeLo",
                  "fileTimeLo", 149);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 149);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 149);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 149);
  sf_mex_assign(&c13_rhs149, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs149, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs149), "rhs",
                  "rhs", 149);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs149), "lhs",
                  "lhs", 149);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 150);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("log2"), "name", "name", 150);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 150);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log2.m"), "resolved",
                  "resolved", 150);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1343862782U), "fileTimeLo",
                  "fileTimeLo", 150);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 150);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 150);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 150);
  sf_mex_assign(&c13_rhs150, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs150, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs150), "rhs",
                  "rhs", 150);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs150), "lhs",
                  "lhs", 150);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/log2.m!scalar_frexp"),
                  "context", "context", 151);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("isfinite"), "name", "name",
                  151);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 151);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "resolved",
                  "resolved", 151);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742656U), "fileTimeLo",
                  "fileTimeLo", 151);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 151);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 151);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 151);
  sf_mex_assign(&c13_rhs151, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs151, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs151), "rhs",
                  "rhs", 151);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs151), "lhs",
                  "lhs", 151);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 152);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 152);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 152);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 152);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 152);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 152);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 152);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 152);
  sf_mex_assign(&c13_rhs152, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs152, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs152), "rhs",
                  "rhs", 152);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs152), "lhs",
                  "lhs", 152);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 153);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("isinf"), "name", "name", 153);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 153);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 153);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742656U), "fileTimeLo",
                  "fileTimeLo", 153);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 153);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 153);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 153);
  sf_mex_assign(&c13_rhs153, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs153, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs153), "rhs",
                  "rhs", 153);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs153), "lhs",
                  "lhs", 153);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isfinite.m"), "context",
                  "context", 154);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("isnan"), "name", "name", 154);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 154);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 154);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742658U), "fileTimeLo",
                  "fileTimeLo", 154);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 154);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 154);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 154);
  sf_mex_assign(&c13_rhs154, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs154, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs154), "rhs",
                  "rhs", 154);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs154), "lhs",
                  "lhs", 154);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 155);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("pow2"), "name", "name", 155);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 155);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/pow2.m"), "resolved",
                  "resolved", 155);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1343862782U), "fileTimeLo",
                  "fileTimeLo", 155);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 155);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 155);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 155);
  sf_mex_assign(&c13_rhs155, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs155, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs155), "rhs",
                  "rhs", 155);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs155), "lhs",
                  "lhs", 155);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/pow2.m"), "context",
                  "context", 156);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalar_pow2"), "name",
                  "name", 156);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 156);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_pow2.m"),
                  "resolved", "resolved", 156);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286851132U), "fileTimeLo",
                  "fileTimeLo", 156);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 156);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 156);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 156);
  sf_mex_assign(&c13_rhs156, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs156, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs156), "rhs",
                  "rhs", 156);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs156), "lhs",
                  "lhs", 156);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_pow2.m"),
                  "context", "context", 157);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("power"), "name", "name", 157);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 157);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 157);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742680U), "fileTimeLo",
                  "fileTimeLo", 157);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 157);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 157);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 157);
  sf_mex_assign(&c13_rhs157, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs157, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs157), "rhs",
                  "rhs", 157);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs157), "lhs",
                  "lhs", 157);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 158);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_error"), "name", "name",
                  158);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 158);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 158);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1343862758U), "fileTimeLo",
                  "fileTimeLo", 158);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 158);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 158);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 158);
  sf_mex_assign(&c13_rhs158, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs158, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs158), "rhs",
                  "rhs", 158);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs158), "lhs",
                  "lhs", 158);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 159);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_div"), "name", "name",
                  159);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 159);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 159);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742666U), "fileTimeLo",
                  "fileTimeLo", 159);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 159);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 159);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 159);
  sf_mex_assign(&c13_rhs159, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs159, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs159), "rhs",
                  "rhs", 159);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs159), "lhs",
                  "lhs", 159);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/expm.m"), "context",
                  "context", 160);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mtimes"), "name", "name",
                  160);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 160);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 160);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 160);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 160);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 160);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 160);
  sf_mex_assign(&c13_rhs160, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c13_lhs160, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs160), "rhs",
                  "rhs", 160);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs160), "lhs",
                  "lhs", 160);
  sf_mex_destroy(&c13_rhs128);
  sf_mex_destroy(&c13_lhs128);
  sf_mex_destroy(&c13_rhs129);
  sf_mex_destroy(&c13_lhs129);
  sf_mex_destroy(&c13_rhs130);
  sf_mex_destroy(&c13_lhs130);
  sf_mex_destroy(&c13_rhs131);
  sf_mex_destroy(&c13_lhs131);
  sf_mex_destroy(&c13_rhs132);
  sf_mex_destroy(&c13_lhs132);
  sf_mex_destroy(&c13_rhs133);
  sf_mex_destroy(&c13_lhs133);
  sf_mex_destroy(&c13_rhs134);
  sf_mex_destroy(&c13_lhs134);
  sf_mex_destroy(&c13_rhs135);
  sf_mex_destroy(&c13_lhs135);
  sf_mex_destroy(&c13_rhs136);
  sf_mex_destroy(&c13_lhs136);
  sf_mex_destroy(&c13_rhs137);
  sf_mex_destroy(&c13_lhs137);
  sf_mex_destroy(&c13_rhs138);
  sf_mex_destroy(&c13_lhs138);
  sf_mex_destroy(&c13_rhs139);
  sf_mex_destroy(&c13_lhs139);
  sf_mex_destroy(&c13_rhs140);
  sf_mex_destroy(&c13_lhs140);
  sf_mex_destroy(&c13_rhs141);
  sf_mex_destroy(&c13_lhs141);
  sf_mex_destroy(&c13_rhs142);
  sf_mex_destroy(&c13_lhs142);
  sf_mex_destroy(&c13_rhs143);
  sf_mex_destroy(&c13_lhs143);
  sf_mex_destroy(&c13_rhs144);
  sf_mex_destroy(&c13_lhs144);
  sf_mex_destroy(&c13_rhs145);
  sf_mex_destroy(&c13_lhs145);
  sf_mex_destroy(&c13_rhs146);
  sf_mex_destroy(&c13_lhs146);
  sf_mex_destroy(&c13_rhs147);
  sf_mex_destroy(&c13_lhs147);
  sf_mex_destroy(&c13_rhs148);
  sf_mex_destroy(&c13_lhs148);
  sf_mex_destroy(&c13_rhs149);
  sf_mex_destroy(&c13_lhs149);
  sf_mex_destroy(&c13_rhs150);
  sf_mex_destroy(&c13_lhs150);
  sf_mex_destroy(&c13_rhs151);
  sf_mex_destroy(&c13_lhs151);
  sf_mex_destroy(&c13_rhs152);
  sf_mex_destroy(&c13_lhs152);
  sf_mex_destroy(&c13_rhs153);
  sf_mex_destroy(&c13_lhs153);
  sf_mex_destroy(&c13_rhs154);
  sf_mex_destroy(&c13_lhs154);
  sf_mex_destroy(&c13_rhs155);
  sf_mex_destroy(&c13_lhs155);
  sf_mex_destroy(&c13_rhs156);
  sf_mex_destroy(&c13_lhs156);
  sf_mex_destroy(&c13_rhs157);
  sf_mex_destroy(&c13_lhs157);
  sf_mex_destroy(&c13_rhs158);
  sf_mex_destroy(&c13_lhs158);
  sf_mex_destroy(&c13_rhs159);
  sf_mex_destroy(&c13_lhs159);
  sf_mex_destroy(&c13_rhs160);
  sf_mex_destroy(&c13_lhs160);
}

static real_T c13_mpower(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T
  c13_a)
{
  real_T c13_b_a;
  real_T c13_c_a;
  real_T c13_ak;
  real_T c13_d_a;
  real_T c13_e_a;
  real_T c13_b;
  c13_b_a = c13_a;
  c13_c_a = c13_b_a;
  c13_eml_scalar_eg(chartInstance);
  c13_ak = c13_c_a;
  c13_d_a = c13_ak;
  c13_eml_scalar_eg(chartInstance);
  c13_e_a = c13_d_a;
  c13_b = c13_d_a;
  return c13_e_a * c13_b;
}

static void c13_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static real_T c13_sqrt(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T c13_x)
{
  real_T c13_b_x;
  c13_b_x = c13_x;
  c13_b_sqrt(chartInstance, &c13_b_x);
  return c13_b_x;
}

static void c13_eml_error(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
  int32_T c13_i268;
  static char_T c13_cv0[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c13_u[30];
  const mxArray *c13_y = NULL;
  int32_T c13_i269;
  static char_T c13_cv1[4] = { 's', 'q', 'r', 't' };

  char_T c13_b_u[4];
  const mxArray *c13_b_y = NULL;
  for (c13_i268 = 0; c13_i268 < 30; c13_i268++) {
    c13_u[c13_i268] = c13_cv0[c13_i268];
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 10, 0U, 1U, 0U, 2, 1, 30),
                FALSE);
  for (c13_i269 = 0; c13_i269 < 4; c13_i269++) {
    c13_b_u[c13_i269] = c13_cv1[c13_i269];
  }

  c13_b_y = NULL;
  sf_mex_assign(&c13_b_y, sf_mex_create("y", c13_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                FALSE);
  sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U, 14,
    c13_y, 14, c13_b_y));
}

static void c13_b_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static void c13_c_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static void c13_d_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static void c13_eye(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T c13_I[25])
{
  int32_T c13_i270;
  int32_T c13_k;
  int32_T c13_b_k;
  for (c13_i270 = 0; c13_i270 < 25; c13_i270++) {
    c13_I[c13_i270] = 0.0;
  }

  for (c13_k = 1; c13_k < 6; c13_k++) {
    c13_b_k = c13_k;
    c13_I[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
             (real_T)c13_b_k), 1, 5, 1, 0) + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
             (int32_T)_SFD_INTEGER_CHECK("", (real_T)c13_b_k), 1, 5, 2, 0) - 1))
      - 1] = 1.0;
  }
}

static void c13_b_eml_error(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
  int32_T c13_i271;
  static char_T c13_cv2[48] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'c', 'h', 'o', 'l', '_', 'm', 'a', 't', 'r', 'i', 'x', 'M',
    'u', 's', 't', 'B', 'e', 'P', 'o', 's', 'D', 'e', 'f', 'W', 'i', 't', 'h',
    'R', 'e', 'a', 'l', 'D', 'i', 'a', 'g' };

  char_T c13_u[48];
  const mxArray *c13_y = NULL;
  for (c13_i271 = 0; c13_i271 < 48; c13_i271++) {
    c13_u[c13_i271] = c13_cv2[c13_i271];
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 10, 0U, 1U, 0U, 2, 1, 48),
                FALSE);
  sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U, 14,
    c13_y));
}

static void c13_eml_matlab_zpotrf(SFc13_IMM_UKFInstanceStruct *chartInstance,
  real_T c13_A[25], real_T c13_b_A[25], int32_T *c13_info)
{
  int32_T c13_i272;
  for (c13_i272 = 0; c13_i272 < 25; c13_i272++) {
    c13_b_A[c13_i272] = c13_A[c13_i272];
  }

  *c13_info = c13_b_eml_matlab_zpotrf(chartInstance, c13_b_A);
}

static void c13_check_forloop_overflow_error(SFc13_IMM_UKFInstanceStruct
  *chartInstance, boolean_T c13_overflow)
{
  int32_T c13_i273;
  static char_T c13_cv3[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c13_u[34];
  const mxArray *c13_y = NULL;
  int32_T c13_i274;
  static char_T c13_cv4[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c13_b_u[23];
  const mxArray *c13_b_y = NULL;
  if (!c13_overflow) {
  } else {
    for (c13_i273 = 0; c13_i273 < 34; c13_i273++) {
      c13_u[c13_i273] = c13_cv3[c13_i273];
    }

    c13_y = NULL;
    sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 10, 0U, 1U, 0U, 2, 1, 34),
                  FALSE);
    for (c13_i274 = 0; c13_i274 < 23; c13_i274++) {
      c13_b_u[c13_i274] = c13_cv4[c13_i274];
    }

    c13_b_y = NULL;
    sf_mex_assign(&c13_b_y, sf_mex_create("y", c13_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                  FALSE);
    sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U,
      14, c13_y, 14, c13_b_y));
  }
}

static void c13_eml_xgemv(SFc13_IMM_UKFInstanceStruct *chartInstance, int32_T
  c13_m, int32_T c13_n, int32_T c13_ia0, int32_T c13_ix0, real_T c13_y[25],
  int32_T c13_iy0, real_T c13_b_y[25])
{
  int32_T c13_i275;
  for (c13_i275 = 0; c13_i275 < 25; c13_i275++) {
    c13_b_y[c13_i275] = c13_y[c13_i275];
  }

  c13_b_eml_xgemv(chartInstance, c13_m, c13_n, c13_ia0, c13_ix0, c13_b_y,
                  c13_iy0);
}

static void c13_c_eml_error(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
  int32_T c13_i276;
  static char_T c13_cv5[19] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'p', 'o', 's', 'd', 'e', 'f' };

  char_T c13_u[19];
  const mxArray *c13_y = NULL;
  for (c13_i276 = 0; c13_i276 < 19; c13_i276++) {
    c13_u[c13_i276] = c13_cv5[c13_i276];
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 10, 0U, 1U, 0U, 2, 1, 19),
                FALSE);
  sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U, 14,
    c13_y));
}

static void c13_e_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static void c13_f_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static void c13_g_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static void c13_inv(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T c13_x[4],
                    real_T c13_y[4])
{
  int32_T c13_i277;
  real_T c13_b_x[4];
  int32_T c13_i278;
  real_T c13_c_x[4];
  real_T c13_n1x;
  int32_T c13_i279;
  real_T c13_b_y[4];
  real_T c13_n1xinv;
  real_T c13_a;
  real_T c13_b;
  real_T c13_c_y;
  real_T c13_rc;
  real_T c13_d_x;
  boolean_T c13_b_b;
  real_T c13_e_x;
  int32_T c13_i280;
  static char_T c13_cv6[8] = { '%', '%', '%', 'd', '.', '%', 'd', 'e' };

  char_T c13_u[8];
  const mxArray *c13_d_y = NULL;
  real_T c13_b_u;
  const mxArray *c13_e_y = NULL;
  real_T c13_c_u;
  const mxArray *c13_f_y = NULL;
  real_T c13_d_u;
  const mxArray *c13_g_y = NULL;
  char_T c13_str[14];
  int32_T c13_i281;
  char_T c13_b_str[14];
  boolean_T guard1 = FALSE;
  boolean_T guard2 = FALSE;
  boolean_T guard3 = FALSE;
  for (c13_i277 = 0; c13_i277 < 4; c13_i277++) {
    c13_b_x[c13_i277] = c13_x[c13_i277];
  }

  c13_inv2x2(chartInstance, c13_b_x, c13_y);
  for (c13_i278 = 0; c13_i278 < 4; c13_i278++) {
    c13_c_x[c13_i278] = c13_x[c13_i278];
  }

  c13_n1x = c13_norm(chartInstance, c13_c_x);
  for (c13_i279 = 0; c13_i279 < 4; c13_i279++) {
    c13_b_y[c13_i279] = c13_y[c13_i279];
  }

  c13_n1xinv = c13_norm(chartInstance, c13_b_y);
  c13_a = c13_n1x;
  c13_b = c13_n1xinv;
  c13_c_y = c13_a * c13_b;
  c13_rc = 1.0 / c13_c_y;
  guard1 = FALSE;
  guard2 = FALSE;
  if (c13_n1x == 0.0) {
    guard2 = TRUE;
  } else if (c13_n1xinv == 0.0) {
    guard2 = TRUE;
  } else if (c13_rc == 0.0) {
    guard1 = TRUE;
  } else {
    c13_d_x = c13_rc;
    c13_b_b = muDoubleScalarIsNaN(c13_d_x);
    guard3 = FALSE;
    if (c13_b_b) {
      guard3 = TRUE;
    } else {
      if (c13_rc < 2.2204460492503131E-16) {
        guard3 = TRUE;
      }
    }

    if (guard3 == TRUE) {
      c13_e_x = c13_rc;
      for (c13_i280 = 0; c13_i280 < 8; c13_i280++) {
        c13_u[c13_i280] = c13_cv6[c13_i280];
      }

      c13_d_y = NULL;
      sf_mex_assign(&c13_d_y, sf_mex_create("y", c13_u, 10, 0U, 1U, 0U, 2, 1, 8),
                    FALSE);
      c13_b_u = 14.0;
      c13_e_y = NULL;
      sf_mex_assign(&c13_e_y, sf_mex_create("y", &c13_b_u, 0, 0U, 0U, 0U, 0),
                    FALSE);
      c13_c_u = 6.0;
      c13_f_y = NULL;
      sf_mex_assign(&c13_f_y, sf_mex_create("y", &c13_c_u, 0, 0U, 0U, 0U, 0),
                    FALSE);
      c13_d_u = c13_e_x;
      c13_g_y = NULL;
      sf_mex_assign(&c13_g_y, sf_mex_create("y", &c13_d_u, 0, 0U, 0U, 0U, 0),
                    FALSE);
      c13_s_emlrt_marshallIn(chartInstance, sf_mex_call_debug("sprintf", 1U, 2U,
        14, sf_mex_call_debug("sprintf", 1U, 3U, 14, c13_d_y, 14, c13_e_y, 14,
        c13_f_y), 14, c13_g_y), "sprintf", c13_str);
      for (c13_i281 = 0; c13_i281 < 14; c13_i281++) {
        c13_b_str[c13_i281] = c13_str[c13_i281];
      }

      c13_b_eml_warning(chartInstance, c13_b_str);
    }
  }

  if (guard2 == TRUE) {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    c13_eml_warning(chartInstance);
  }
}

static void c13_inv2x2(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T c13_x
  [4], real_T c13_y[4])
{
  real_T c13_b_x;
  real_T c13_c_x;
  real_T c13_d_x;
  real_T c13_b_y;
  real_T c13_e_x;
  real_T c13_f_x;
  real_T c13_c_y;
  real_T c13_d;
  real_T c13_g_x;
  real_T c13_h_x;
  real_T c13_i_x;
  real_T c13_d_y;
  real_T c13_j_x;
  real_T c13_k_x;
  real_T c13_e_y;
  real_T c13_b_d;
  real_T c13_A;
  real_T c13_B;
  real_T c13_l_x;
  real_T c13_f_y;
  real_T c13_m_x;
  real_T c13_g_y;
  real_T c13_r;
  real_T c13_a;
  real_T c13_b;
  real_T c13_h_y;
  real_T c13_b_B;
  real_T c13_i_y;
  real_T c13_j_y;
  real_T c13_t;
  real_T c13_b_A;
  real_T c13_c_B;
  real_T c13_n_x;
  real_T c13_k_y;
  real_T c13_o_x;
  real_T c13_l_y;
  real_T c13_m_y;
  real_T c13_b_a;
  real_T c13_b_b;
  real_T c13_n_y;
  real_T c13_c_A;
  real_T c13_d_B;
  real_T c13_p_x;
  real_T c13_o_y;
  real_T c13_q_x;
  real_T c13_p_y;
  real_T c13_q_y;
  real_T c13_c_a;
  real_T c13_c_b;
  real_T c13_r_y;
  real_T c13_d_a;
  real_T c13_d_b;
  real_T c13_s_y;
  real_T c13_d_A;
  real_T c13_e_B;
  real_T c13_r_x;
  real_T c13_t_y;
  real_T c13_s_x;
  real_T c13_u_y;
  real_T c13_e_a;
  real_T c13_e_b;
  real_T c13_v_y;
  real_T c13_f_B;
  real_T c13_w_y;
  real_T c13_x_y;
  real_T c13_e_A;
  real_T c13_g_B;
  real_T c13_t_x;
  real_T c13_y_y;
  real_T c13_u_x;
  real_T c13_ab_y;
  real_T c13_bb_y;
  real_T c13_f_a;
  real_T c13_f_b;
  real_T c13_cb_y;
  real_T c13_g_a;
  real_T c13_g_b;
  real_T c13_db_y;
  real_T c13_f_A;
  real_T c13_h_B;
  real_T c13_v_x;
  real_T c13_eb_y;
  real_T c13_w_x;
  real_T c13_fb_y;
  real_T c13_gb_y;
  real_T c13_h_a;
  real_T c13_h_b;
  real_T c13_hb_y;
  c13_b_x = c13_x[1];
  c13_c_x = c13_b_x;
  c13_d_x = c13_c_x;
  c13_b_y = muDoubleScalarAbs(c13_d_x);
  c13_e_x = 0.0;
  c13_f_x = c13_e_x;
  c13_c_y = muDoubleScalarAbs(c13_f_x);
  c13_d = c13_b_y + c13_c_y;
  c13_g_x = c13_x[0];
  c13_h_x = c13_g_x;
  c13_i_x = c13_h_x;
  c13_d_y = muDoubleScalarAbs(c13_i_x);
  c13_j_x = 0.0;
  c13_k_x = c13_j_x;
  c13_e_y = muDoubleScalarAbs(c13_k_x);
  c13_b_d = c13_d_y + c13_e_y;
  if (c13_d > c13_b_d) {
    c13_A = c13_x[0];
    c13_B = c13_x[1];
    c13_l_x = c13_A;
    c13_f_y = c13_B;
    c13_m_x = c13_l_x;
    c13_g_y = c13_f_y;
    c13_r = c13_m_x / c13_g_y;
    c13_a = c13_r;
    c13_b = c13_x[3];
    c13_h_y = c13_a * c13_b;
    c13_b_B = c13_h_y - c13_x[2];
    c13_i_y = c13_b_B;
    c13_j_y = c13_i_y;
    c13_t = 1.0 / c13_j_y;
    c13_b_A = c13_x[3];
    c13_c_B = c13_x[1];
    c13_n_x = c13_b_A;
    c13_k_y = c13_c_B;
    c13_o_x = c13_n_x;
    c13_l_y = c13_k_y;
    c13_m_y = c13_o_x / c13_l_y;
    c13_b_a = c13_m_y;
    c13_b_b = c13_t;
    c13_n_y = c13_b_a * c13_b_b;
    c13_y[0] = c13_n_y;
    c13_y[1] = -c13_t;
    c13_c_A = -c13_x[2];
    c13_d_B = c13_x[1];
    c13_p_x = c13_c_A;
    c13_o_y = c13_d_B;
    c13_q_x = c13_p_x;
    c13_p_y = c13_o_y;
    c13_q_y = c13_q_x / c13_p_y;
    c13_c_a = c13_q_y;
    c13_c_b = c13_t;
    c13_r_y = c13_c_a * c13_c_b;
    c13_y[2] = c13_r_y;
    c13_d_a = c13_r;
    c13_d_b = c13_t;
    c13_s_y = c13_d_a * c13_d_b;
    c13_y[3] = c13_s_y;
  } else {
    c13_d_A = c13_x[1];
    c13_e_B = c13_x[0];
    c13_r_x = c13_d_A;
    c13_t_y = c13_e_B;
    c13_s_x = c13_r_x;
    c13_u_y = c13_t_y;
    c13_r = c13_s_x / c13_u_y;
    c13_e_a = c13_r;
    c13_e_b = c13_x[2];
    c13_v_y = c13_e_a * c13_e_b;
    c13_f_B = c13_x[3] - c13_v_y;
    c13_w_y = c13_f_B;
    c13_x_y = c13_w_y;
    c13_t = 1.0 / c13_x_y;
    c13_e_A = c13_x[3];
    c13_g_B = c13_x[0];
    c13_t_x = c13_e_A;
    c13_y_y = c13_g_B;
    c13_u_x = c13_t_x;
    c13_ab_y = c13_y_y;
    c13_bb_y = c13_u_x / c13_ab_y;
    c13_f_a = c13_bb_y;
    c13_f_b = c13_t;
    c13_cb_y = c13_f_a * c13_f_b;
    c13_y[0] = c13_cb_y;
    c13_g_a = -c13_r;
    c13_g_b = c13_t;
    c13_db_y = c13_g_a * c13_g_b;
    c13_y[1] = c13_db_y;
    c13_f_A = -c13_x[2];
    c13_h_B = c13_x[0];
    c13_v_x = c13_f_A;
    c13_eb_y = c13_h_B;
    c13_w_x = c13_v_x;
    c13_fb_y = c13_eb_y;
    c13_gb_y = c13_w_x / c13_fb_y;
    c13_h_a = c13_gb_y;
    c13_h_b = c13_t;
    c13_hb_y = c13_h_a * c13_h_b;
    c13_y[2] = c13_hb_y;
    c13_y[3] = c13_t;
  }
}

static real_T c13_norm(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T c13_x
  [4])
{
  real_T c13_y;
  int32_T c13_j;
  real_T c13_b_j;
  real_T c13_s;
  int32_T c13_i;
  real_T c13_b_i;
  real_T c13_b_x;
  real_T c13_c_x;
  real_T c13_b_y;
  real_T c13_d_x;
  boolean_T c13_b;
  boolean_T exitg1;
  c13_y = 0.0;
  c13_j = 0;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (c13_j < 2)) {
    c13_b_j = 1.0 + (real_T)c13_j;
    c13_s = 0.0;
    for (c13_i = 0; c13_i < 2; c13_i++) {
      c13_b_i = 1.0 + (real_T)c13_i;
      c13_b_x = c13_x[((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c13_b_i), 1, 2, 1, 0) + (((int32_T)(real_T)
        _SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", c13_b_j),
        1, 2, 2, 0) - 1) << 1)) - 1];
      c13_c_x = c13_b_x;
      c13_b_y = muDoubleScalarAbs(c13_c_x);
      c13_s += c13_b_y;
    }

    c13_d_x = c13_s;
    c13_b = muDoubleScalarIsNaN(c13_d_x);
    if (c13_b) {
      c13_y = rtNaN;
      exitg1 = TRUE;
    } else {
      if (c13_s > c13_y) {
        c13_y = c13_s;
      }

      c13_j++;
    }
  }

  return c13_y;
}

static void c13_eml_warning(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
  int32_T c13_i282;
  static char_T c13_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c13_u[27];
  const mxArray *c13_y = NULL;
  for (c13_i282 = 0; c13_i282 < 27; c13_i282++) {
    c13_u[c13_i282] = c13_varargin_1[c13_i282];
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 10, 0U, 1U, 0U, 2, 1, 27),
                FALSE);
  sf_mex_call_debug("warning", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U,
    14, c13_y));
}

static void c13_b_eml_warning(SFc13_IMM_UKFInstanceStruct *chartInstance, char_T
  c13_varargin_2[14])
{
  int32_T c13_i283;
  static char_T c13_varargin_1[33] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 'i', 'l', 'l', 'C', 'o', 'n', 'd', 'i', 't', 'i',
    'o', 'n', 'e', 'd', 'M', 'a', 't', 'r', 'i', 'x' };

  char_T c13_u[33];
  const mxArray *c13_y = NULL;
  int32_T c13_i284;
  char_T c13_b_u[14];
  const mxArray *c13_b_y = NULL;
  for (c13_i283 = 0; c13_i283 < 33; c13_i283++) {
    c13_u[c13_i283] = c13_varargin_1[c13_i283];
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 10, 0U, 1U, 0U, 2, 1, 33),
                FALSE);
  for (c13_i284 = 0; c13_i284 < 14; c13_i284++) {
    c13_b_u[c13_i284] = c13_varargin_2[c13_i284];
  }

  c13_b_y = NULL;
  sf_mex_assign(&c13_b_y, sf_mex_create("y", c13_b_u, 10, 0U, 1U, 0U, 2, 1, 14),
                FALSE);
  sf_mex_call_debug("warning", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U,
    14, c13_y, 14, c13_b_y));
}

static void c13_h_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static void c13_i_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static void c13_j_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static void c13_k_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static void c13_l_eml_scalar_eg(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

static real_T c13_expm(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T c13_A)
{
  real_T c13_F;
  real_T c13_x;
  real_T c13_b_x;
  real_T c13_c_x;
  real_T c13_normA;
  int32_T c13_i;
  real_T c13_b_i;
  static real_T c13_theta[5] = { 0.01495585217958292, 0.253939833006323,
    0.95041789961629319, 2.097847961257068, 5.3719203511481517 };

  static real_T c13_dv19[5] = { 3.0, 5.0, 7.0, 9.0, 13.0 };

  real_T c13_d_x;
  real_T c13_e_x;
  real_T c13_f_x;
  real_T c13_g_x;
  boolean_T c13_b;
  boolean_T c13_b0;
  real_T c13_h_x;
  boolean_T c13_b_b;
  boolean_T c13_b1;
  boolean_T c13_c_b;
  int32_T c13_eint;
  real_T c13_fdbl;
  int32_T c13_b_eint;
  real_T c13_b_fdbl;
  int32_T c13_c_eint;
  real_T c13_d1;
  real_T c13_d2;
  real_T c13_t;
  real_T c13_s;
  real_T c13_b_t;
  real_T c13_b_s;
  real_T c13_a;
  real_T c13_b_a;
  real_T c13_d_b;
  real_T c13_e_b;
  real_T c13_bk;
  real_T c13_f_b;
  real_T c13_br;
  real_T c13_y;
  real_T c13_i_x;
  real_T c13_b_y;
  real_T c13_c_s;
  int32_T c13_i285;
  int32_T c13_j;
  real_T c13_c_a;
  real_T c13_g_b;
  boolean_T exitg1;
  c13_x = c13_A;
  c13_b_x = c13_x;
  c13_c_x = c13_b_x;
  c13_normA = muDoubleScalarAbs(c13_c_x);
  if (c13_normA <= 5.3719203511481517) {
    c13_i = 0;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (c13_i < 5)) {
      c13_b_i = 1.0 + (real_T)c13_i;
      if (c13_normA <= c13_theta[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("",
           (int32_T)_SFD_INTEGER_CHECK("", c13_b_i), 1, 5, 1, 0) - 1]) {
        c13_F = c13_PadeApproximantOfDegree(chartInstance, c13_A, c13_dv19
          [(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", c13_b_i), 1, 5, 1, 0) - 1]);
        exitg1 = TRUE;
      } else {
        c13_i++;
      }
    }
  } else {
    c13_d_x = c13_normA / 5.3719203511481517;
    c13_e_x = c13_d_x;
    c13_f_x = c13_e_x;
    c13_g_x = c13_f_x;
    c13_b = muDoubleScalarIsInf(c13_g_x);
    c13_b0 = !c13_b;
    c13_h_x = c13_f_x;
    c13_b_b = muDoubleScalarIsNaN(c13_h_x);
    c13_b1 = !c13_b_b;
    c13_c_b = (c13_b0 && c13_b1);
    if (c13_c_b) {
      c13_fdbl = frexp(c13_e_x, &c13_eint);
      c13_b_eint = c13_eint;
      c13_b_fdbl = c13_fdbl;
      c13_c_eint = c13_b_eint;
      c13_d1 = c13_b_fdbl;
      c13_d2 = (real_T)c13_c_eint;
    } else {
      c13_d1 = c13_e_x;
      c13_d2 = 0.0;
    }

    c13_t = c13_d1;
    c13_s = c13_d2;
    c13_b_t = c13_t;
    c13_b_s = c13_s;
    if (c13_b_t == 0.5) {
      c13_b_s--;
    }

    c13_a = c13_b_s;
    c13_b_a = c13_a;
    c13_d_b = c13_b_a;
    c13_e_b = c13_d_b;
    c13_eml_scalar_eg(chartInstance);
    c13_bk = c13_e_b;
    c13_f_b = c13_bk;
    c13_eml_scalar_eg(chartInstance);
    c13_br = c13_f_b;
    c13_y = muDoubleScalarPower(2.0, c13_br);
    c13_i_x = c13_A;
    c13_b_y = c13_y;
    c13_A = c13_i_x / c13_b_y;
    c13_F = c13_PadeApproximantOfDegree(chartInstance, c13_A, 13.0);
    c13_c_s = c13_b_s;
    c13_i285 = (int32_T)c13_c_s;
    _SFD_FOR_LOOP_VECTOR_CHECK(1.0, 1.0, c13_c_s, mxDOUBLE_CLASS, c13_i285);
    for (c13_j = 0; c13_j < c13_i285; c13_j++) {
      c13_c_a = c13_F;
      c13_g_b = c13_F;
      c13_F = c13_c_a * c13_g_b;
    }
  }

  return c13_F;
}

static real_T c13_PadeApproximantOfDegree(SFc13_IMM_UKFInstanceStruct
  *chartInstance, real_T c13_A, real_T c13_m)
{
  real_T c13_a;
  real_T c13_b;
  real_T c13_A2;
  real_T c13_U;
  real_T c13_b_a;
  real_T c13_b_b;
  real_T c13_c_b;
  real_T c13_V;
  real_T c13_d;
  real_T c13_c_a;
  real_T c13_d_b;
  real_T c13_A3;
  real_T c13_e_b;
  real_T c13_y;
  real_T c13_d_a;
  real_T c13_f_b;
  real_T c13_g_b;
  real_T c13_b_y;
  real_T c13_h_b;
  real_T c13_c_y;
  real_T c13_e_a;
  real_T c13_i_b;
  real_T c13_A4;
  real_T c13_j_b;
  real_T c13_d_y;
  real_T c13_k_b;
  real_T c13_e_y;
  real_T c13_f_a;
  real_T c13_l_b;
  real_T c13_m_b;
  real_T c13_f_y;
  real_T c13_n_b;
  real_T c13_g_y;
  real_T c13_o_b;
  real_T c13_h_y;
  real_T c13_g_a;
  real_T c13_p_b;
  real_T c13_q_b;
  real_T c13_i_y;
  real_T c13_r_b;
  real_T c13_j_y;
  real_T c13_s_b;
  real_T c13_k_y;
  real_T c13_h_a;
  real_T c13_t_b;
  real_T c13_u_b;
  real_T c13_l_y;
  real_T c13_v_b;
  real_T c13_m_y;
  real_T c13_w_b;
  real_T c13_n_y;
  real_T c13_x_b;
  real_T c13_o_y;
  real_T c13_y_b;
  real_T c13_p_y;
  real_T c13_ab_b;
  real_T c13_q_y;
  real_T c13_bb_b;
  real_T c13_r_y;
  real_T c13_cb_b;
  real_T c13_s_y;
  real_T c13_db_b;
  real_T c13_t_y;
  real_T c13_i_a;
  real_T c13_eb_b;
  real_T c13_u_y;
  real_T c13_j_a;
  real_T c13_fb_b;
  real_T c13_gb_b;
  real_T c13_v_y;
  real_T c13_hb_b;
  real_T c13_w_y;
  real_T c13_ib_b;
  real_T c13_x_y;
  real_T c13_k_a;
  real_T c13_jb_b;
  real_T c13_y_y;
  real_T c13_kb_b;
  real_T c13_ab_y;
  real_T c13_lb_b;
  real_T c13_bb_y;
  real_T c13_mb_b;
  real_T c13_cb_y;
  real_T c13_uk;
  real_T c13_b_A;
  real_T c13_B;
  real_T c13_x;
  real_T c13_db_y;
  real_T c13_b_x;
  real_T c13_eb_y;
  c13_a = c13_A;
  c13_b = c13_A;
  c13_A2 = c13_a * c13_b;
  if (c13_m == 3.0) {
    c13_U = c13_A2;
    c13_U += 60.0;
    c13_b_a = c13_A;
    c13_b_b = c13_U;
    c13_U = c13_b_a * c13_b_b;
    c13_c_b = c13_A2;
    c13_V = 12.0 * c13_c_b;
    c13_d = 120.0;
  } else {
    c13_c_a = c13_A2;
    c13_d_b = c13_A2;
    c13_A3 = c13_c_a * c13_d_b;
    if (c13_m == 5.0) {
      c13_e_b = c13_A2;
      c13_y = 420.0 * c13_e_b;
      c13_U = c13_A3 + c13_y;
      c13_U += 15120.0;
      c13_d_a = c13_A;
      c13_f_b = c13_U;
      c13_U = c13_d_a * c13_f_b;
      c13_g_b = c13_A3;
      c13_b_y = 30.0 * c13_g_b;
      c13_h_b = c13_A2;
      c13_c_y = 3360.0 * c13_h_b;
      c13_V = c13_b_y + c13_c_y;
      c13_d = 30240.0;
    } else {
      c13_e_a = c13_A3;
      c13_i_b = c13_A2;
      c13_A4 = c13_e_a * c13_i_b;
      if (c13_m == 7.0) {
        c13_j_b = c13_A3;
        c13_d_y = 1512.0 * c13_j_b;
        c13_k_b = c13_A2;
        c13_e_y = 277200.0 * c13_k_b;
        c13_U = (c13_A4 + c13_d_y) + c13_e_y;
        c13_U += 8.64864E+6;
        c13_f_a = c13_A;
        c13_l_b = c13_U;
        c13_U = c13_f_a * c13_l_b;
        c13_m_b = c13_A4;
        c13_f_y = 56.0 * c13_m_b;
        c13_n_b = c13_A3;
        c13_g_y = 25200.0 * c13_n_b;
        c13_o_b = c13_A2;
        c13_h_y = 1.99584E+6 * c13_o_b;
        c13_V = (c13_f_y + c13_g_y) + c13_h_y;
        c13_d = 1.729728E+7;
      } else if (c13_m == 9.0) {
        c13_g_a = c13_A4;
        c13_p_b = c13_A2;
        c13_V = c13_g_a * c13_p_b;
        c13_q_b = c13_A4;
        c13_i_y = 3960.0 * c13_q_b;
        c13_r_b = c13_A3;
        c13_j_y = 2.16216E+6 * c13_r_b;
        c13_s_b = c13_A2;
        c13_k_y = 3.027024E+8 * c13_s_b;
        c13_U = ((c13_V + c13_i_y) + c13_j_y) + c13_k_y;
        c13_U += 8.8216128E+9;
        c13_h_a = c13_A;
        c13_t_b = c13_U;
        c13_U = c13_h_a * c13_t_b;
        c13_u_b = c13_V;
        c13_l_y = 90.0 * c13_u_b;
        c13_v_b = c13_A4;
        c13_m_y = 110880.0 * c13_v_b;
        c13_w_b = c13_A3;
        c13_n_y = 3.027024E+7 * c13_w_b;
        c13_x_b = c13_A2;
        c13_o_y = 2.0756736E+9 * c13_x_b;
        c13_V = ((c13_l_y + c13_m_y) + c13_n_y) + c13_o_y;
        c13_d = 1.76432256E+10;
      } else {
        c13_y_b = c13_A4;
        c13_p_y = 3.352212864E+10 * c13_y_b;
        c13_ab_b = c13_A3;
        c13_q_y = 1.05594705216E+13 * c13_ab_b;
        c13_bb_b = c13_A2;
        c13_r_y = 1.1873537964288E+15 * c13_bb_b;
        c13_U = (c13_p_y + c13_q_y) + c13_r_y;
        c13_U += 3.238237626624E+16;
        c13_cb_b = c13_A3;
        c13_s_y = 16380.0 * c13_cb_b;
        c13_db_b = c13_A2;
        c13_t_y = 4.08408E+7 * c13_db_b;
        c13_i_a = c13_A4;
        c13_eb_b = (c13_A4 + c13_s_y) + c13_t_y;
        c13_u_y = c13_i_a * c13_eb_b;
        c13_j_a = c13_A;
        c13_fb_b = c13_u_y + c13_U;
        c13_U = c13_j_a * c13_fb_b;
        c13_gb_b = c13_A4;
        c13_v_y = 182.0 * c13_gb_b;
        c13_hb_b = c13_A3;
        c13_w_y = 960960.0 * c13_hb_b;
        c13_ib_b = c13_A2;
        c13_x_y = 1.32324192E+9 * c13_ib_b;
        c13_k_a = c13_A4;
        c13_jb_b = (c13_v_y + c13_w_y) + c13_x_y;
        c13_y_y = c13_k_a * c13_jb_b;
        c13_kb_b = c13_A4;
        c13_ab_y = 6.704425728E+11 * c13_kb_b;
        c13_lb_b = c13_A3;
        c13_bb_y = 1.29060195264E+14 * c13_lb_b;
        c13_mb_b = c13_A2;
        c13_cb_y = 7.7717703038976E+15 * c13_mb_b;
        c13_V = ((c13_y_y + c13_ab_y) + c13_bb_y) + c13_cb_y;
        c13_d = 6.476475253248E+16;
      }
    }
  }

  c13_V += c13_d;
  c13_uk = c13_U;
  c13_U = c13_V - c13_uk;
  c13_V += c13_uk;
  c13_b_A = c13_U;
  c13_B = c13_V;
  c13_x = c13_B;
  c13_db_y = c13_b_A;
  c13_b_x = c13_x;
  c13_eb_y = c13_db_y;
  return c13_b_x / c13_eb_y;
}

static void c13_s_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_sprintf, const char_T *c13_identifier, char_T c13_y[14])
{
  emlrtMsgIdentifier c13_thisId;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_t_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_sprintf), &c13_thisId,
    c13_y);
  sf_mex_destroy(&c13_sprintf);
}

static void c13_t_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId, char_T c13_y[14])
{
  char_T c13_cv7[14];
  int32_T c13_i286;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), c13_cv7, 1, 10, 0U, 1, 0U, 2, 1,
                14);
  for (c13_i286 = 0; c13_i286 < 14; c13_i286++) {
    c13_y[c13_i286] = c13_cv7[c13_i286];
  }

  sf_mex_destroy(&c13_u);
}

static const mxArray *c13_p_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_u;
  const mxArray *c13_y = NULL;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_u = *(int32_T *)c13_inData;
  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", &c13_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, FALSE);
  return c13_mxArrayOutData;
}

static int32_T c13_u_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId)
{
  int32_T c13_y;
  int32_T c13_i287;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), &c13_i287, 1, 6, 0U, 0, 0U, 0);
  c13_y = c13_i287;
  sf_mex_destroy(&c13_u);
  return c13_y;
}

static void c13_m_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_b_sfEvent;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  int32_T c13_y;
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c13_b_sfEvent = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_y = c13_u_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_b_sfEvent),
    &c13_thisId);
  sf_mex_destroy(&c13_b_sfEvent);
  *(int32_T *)c13_outData = c13_y;
  sf_mex_destroy(&c13_mxArrayInData);
}

static uint8_T c13_v_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_b_is_active_c13_IMM_UKF, const char_T *c13_identifier)
{
  uint8_T c13_y;
  emlrtMsgIdentifier c13_thisId;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_y = c13_w_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c13_b_is_active_c13_IMM_UKF), &c13_thisId);
  sf_mex_destroy(&c13_b_is_active_c13_IMM_UKF);
  return c13_y;
}

static uint8_T c13_w_emlrt_marshallIn(SFc13_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId)
{
  uint8_T c13_y;
  uint8_T c13_u0;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), &c13_u0, 1, 3, 0U, 0, 0U, 0);
  c13_y = c13_u0;
  sf_mex_destroy(&c13_u);
  return c13_y;
}

static void c13_b_sqrt(SFc13_IMM_UKFInstanceStruct *chartInstance, real_T *c13_x)
{
  if (*c13_x < 0.0) {
    c13_eml_error(chartInstance);
  }

  *c13_x = muDoubleScalarSqrt(*c13_x);
}

static int32_T c13_b_eml_matlab_zpotrf(SFc13_IMM_UKFInstanceStruct
  *chartInstance, real_T c13_A[25])
{
  int32_T c13_info;
  int32_T c13_j;
  int32_T c13_b_j;
  int32_T c13_a;
  int32_T c13_jm1;
  int32_T c13_b_a;
  int32_T c13_c;
  int32_T c13_c_a;
  int32_T c13_b;
  int32_T c13_jj;
  int32_T c13_n;
  int32_T c13_ix0;
  int32_T c13_iy0;
  int32_T c13_b_n;
  int32_T c13_b_ix0;
  int32_T c13_b_iy0;
  int32_T c13_c_n;
  int32_T c13_c_ix0;
  int32_T c13_c_iy0;
  int32_T c13_d_n;
  int32_T c13_d_ix0;
  int32_T c13_d_iy0;
  int32_T c13_e_n;
  int32_T c13_e_ix0;
  int32_T c13_e_iy0;
  real_T c13_d;
  int32_T c13_ix;
  int32_T c13_iy;
  int32_T c13_f_n;
  int32_T c13_b_b;
  int32_T c13_c_b;
  boolean_T c13_overflow;
  int32_T c13_k;
  int32_T c13_d_a;
  int32_T c13_e_a;
  real_T c13_ajj;
  int32_T c13_d_b;
  int32_T c13_nmj;
  int32_T c13_f_a;
  int32_T c13_jp1;
  int32_T c13_g_a;
  int32_T c13_jp1j;
  int32_T c13_b_jm1;
  int32_T c13_e_b;
  int32_T c13_f_b;
  boolean_T c13_b_overflow;
  int32_T c13_b_k;
  int32_T c13_c_k;
  int32_T c13_c_jm1;
  int32_T c13_g_b;
  int32_T c13_h_b;
  boolean_T c13_c_overflow;
  int32_T c13_d_k;
  real_T c13_y;
  real_T c13_z;
  int32_T c13_g_n;
  real_T c13_h_a;
  int32_T c13_f_ix0;
  int32_T c13_h_n;
  real_T c13_i_a;
  int32_T c13_g_ix0;
  int32_T c13_i_n;
  real_T c13_j_a;
  int32_T c13_h_ix0;
  int32_T c13_i_ix0;
  int32_T c13_k_a;
  int32_T c13_b_c;
  int32_T c13_i_b;
  int32_T c13_c_c;
  int32_T c13_l_a;
  int32_T c13_j_b;
  int32_T c13_i288;
  int32_T c13_m_a;
  int32_T c13_k_b;
  int32_T c13_n_a;
  int32_T c13_l_b;
  boolean_T c13_d_overflow;
  int32_T c13_e_k;
  int32_T c13_f_k;
  boolean_T exitg1;
  c13_info = 0;
  c13_c_eml_scalar_eg(chartInstance);
  c13_j = 1;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (c13_j < 6)) {
    c13_b_j = c13_j;
    c13_a = c13_b_j - 1;
    c13_jm1 = c13_a;
    c13_b_a = c13_jm1;
    c13_c = c13_b_a * 5;
    c13_c_a = c13_b_j;
    c13_b = c13_c;
    c13_jj = c13_c_a + c13_b;
    c13_n = c13_jm1;
    c13_ix0 = c13_b_j;
    c13_iy0 = c13_b_j;
    c13_b_n = c13_n;
    c13_b_ix0 = c13_ix0;
    c13_b_iy0 = c13_iy0;
    c13_c_n = c13_b_n;
    c13_c_ix0 = c13_b_ix0;
    c13_c_iy0 = c13_b_iy0;
    c13_d_n = c13_c_n;
    c13_d_ix0 = c13_c_ix0;
    c13_d_iy0 = c13_c_iy0;
    c13_e_n = c13_d_n;
    c13_e_ix0 = c13_d_ix0;
    c13_e_iy0 = c13_d_iy0;
    c13_d = 0.0;
    if (c13_e_n < 1) {
    } else {
      c13_ix = c13_e_ix0;
      c13_iy = c13_e_iy0;
      c13_f_n = c13_e_n;
      c13_b_b = c13_f_n;
      c13_c_b = c13_b_b;
      if (1 > c13_c_b) {
        c13_overflow = FALSE;
      } else {
        c13_overflow = (c13_c_b > 2147483646);
      }

      if (c13_overflow) {
        c13_check_forloop_overflow_error(chartInstance, c13_overflow);
      }

      for (c13_k = 1; c13_k <= c13_f_n; c13_k++) {
        c13_d += c13_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c13_ix), 1, 25, 1, 0) - 1] *
          c13_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c13_iy), 1, 25, 1, 0) - 1];
        c13_d_a = c13_ix + 5;
        c13_ix = c13_d_a;
        c13_e_a = c13_iy + 5;
        c13_iy = c13_e_a;
      }
    }

    c13_ajj = c13_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c13_jj), 1, 25, 1, 0) - 1] - c13_d;
    if (c13_ajj > 0.0) {
      c13_ajj = muDoubleScalarSqrt(c13_ajj);
      c13_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c13_jj), 1, 25, 1, 0) - 1] = c13_ajj;
      if (c13_b_j < 5) {
        c13_d_b = c13_b_j;
        c13_nmj = 5 - c13_d_b;
        c13_f_a = c13_b_j;
        c13_jp1 = c13_f_a;
        c13_g_a = c13_jj + 1;
        c13_jp1j = c13_g_a;
        c13_b_jm1 = c13_jm1;
        c13_e_b = c13_b_jm1;
        c13_f_b = c13_e_b;
        if (1 > c13_f_b) {
          c13_b_overflow = FALSE;
        } else {
          c13_b_overflow = (c13_f_b > 2147483646);
        }

        if (c13_b_overflow) {
          c13_check_forloop_overflow_error(chartInstance, c13_b_overflow);
        }

        for (c13_b_k = 1; c13_b_k <= c13_b_jm1; c13_b_k++) {
          c13_c_k = c13_b_k;
          c13_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c13_b_j), 1, 5, 1, 0) + 5 *
                 (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c13_c_k), 1, 5, 2, 0) - 1)) - 1] = c13_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c13_b_j), 1, 5, 1, 0) + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c13_c_k), 1, 5, 2, 0)
               - 1)) - 1];
        }

        c13_b_eml_xgemv(chartInstance, c13_nmj, c13_jm1, c13_jp1 + 1, c13_b_j,
                        c13_A, c13_jp1j);
        c13_c_jm1 = c13_jm1;
        c13_g_b = c13_c_jm1;
        c13_h_b = c13_g_b;
        if (1 > c13_h_b) {
          c13_c_overflow = FALSE;
        } else {
          c13_c_overflow = (c13_h_b > 2147483646);
        }

        if (c13_c_overflow) {
          c13_check_forloop_overflow_error(chartInstance, c13_c_overflow);
        }

        for (c13_d_k = 1; c13_d_k <= c13_c_jm1; c13_d_k++) {
          c13_c_k = c13_d_k;
          c13_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c13_b_j), 1, 5, 1, 0) + 5 *
                 (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c13_c_k), 1, 5, 2, 0) - 1)) - 1] = c13_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c13_b_j), 1, 5, 1, 0) + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c13_c_k), 1, 5, 2, 0)
               - 1)) - 1];
        }

        c13_y = c13_ajj;
        c13_z = 1.0 / c13_y;
        c13_g_n = c13_nmj;
        c13_h_a = c13_z;
        c13_f_ix0 = c13_jp1j;
        c13_h_n = c13_g_n;
        c13_i_a = c13_h_a;
        c13_g_ix0 = c13_f_ix0;
        c13_i_n = c13_h_n;
        c13_j_a = c13_i_a;
        c13_h_ix0 = c13_g_ix0;
        c13_i_ix0 = c13_h_ix0;
        c13_k_a = c13_i_n;
        c13_b_c = c13_k_a;
        c13_i_b = c13_b_c - 1;
        c13_c_c = c13_i_b;
        c13_l_a = c13_h_ix0;
        c13_j_b = c13_c_c;
        c13_i288 = c13_l_a + c13_j_b;
        c13_m_a = c13_i_ix0;
        c13_k_b = c13_i288;
        c13_n_a = c13_m_a;
        c13_l_b = c13_k_b;
        if (c13_n_a > c13_l_b) {
          c13_d_overflow = FALSE;
        } else {
          c13_d_overflow = (c13_l_b > 2147483646);
        }

        if (c13_d_overflow) {
          c13_check_forloop_overflow_error(chartInstance, c13_d_overflow);
        }

        for (c13_e_k = c13_i_ix0; c13_e_k <= c13_i288; c13_e_k++) {
          c13_f_k = c13_e_k;
          c13_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c13_f_k), 1, 25, 1, 0) - 1] = c13_j_a *
            c13_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c13_f_k), 1, 25, 1, 0) - 1];
        }
      }

      c13_j++;
    } else {
      c13_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c13_jj), 1, 25, 1, 0) - 1] = c13_ajj;
      c13_info = c13_b_j;
      exitg1 = TRUE;
    }
  }

  return c13_info;
}

static void c13_b_eml_xgemv(SFc13_IMM_UKFInstanceStruct *chartInstance, int32_T
  c13_m, int32_T c13_n, int32_T c13_ia0, int32_T c13_ix0, real_T c13_y[25],
  int32_T c13_iy0)
{
  int32_T c13_b_m;
  int32_T c13_b_n;
  int32_T c13_b_ia0;
  int32_T c13_b_ix0;
  int32_T c13_b_iy0;
  int32_T c13_c_m;
  int32_T c13_c_n;
  real_T c13_alpha1;
  int32_T c13_c_ia0;
  int32_T c13_c_ix0;
  real_T c13_beta1;
  int32_T c13_c_iy0;
  char_T c13_TRANSA;
  int32_T c13_var;
  ptrdiff_t c13_m_t;
  int32_T c13_b_var;
  ptrdiff_t c13_n_t;
  ptrdiff_t c13_lda_t;
  ptrdiff_t c13_incx_t;
  ptrdiff_t c13_incy_t;
  double * c13_alpha1_t;
  double * c13_beta1_t;
  double * c13_yiy0_t;
  double * c13_yix0_t;
  double * c13_yia0_t;
  c13_b_m = c13_m;
  c13_b_n = c13_n;
  c13_b_ia0 = c13_ia0;
  c13_b_ix0 = c13_ix0;
  c13_b_iy0 = c13_iy0;
  if (c13_b_m < 1) {
  } else if (c13_b_n < 1) {
  } else {
    c13_c_m = c13_b_m;
    c13_c_n = c13_b_n;
    c13_alpha1 = -1.0;
    c13_c_ia0 = c13_b_ia0;
    c13_c_ix0 = c13_b_ix0;
    c13_beta1 = 1.0;
    c13_c_iy0 = c13_b_iy0;
    c13_TRANSA = 'N';
    c13_var = c13_c_m;
    c13_m_t = (ptrdiff_t)(c13_var);
    c13_b_var = c13_c_n;
    c13_n_t = (ptrdiff_t)(c13_b_var);
    c13_lda_t = (ptrdiff_t)(5);
    c13_incx_t = (ptrdiff_t)(5);
    c13_incy_t = (ptrdiff_t)(1);
    c13_alpha1_t = (double *)(&c13_alpha1);
    c13_beta1_t = (double *)(&c13_beta1);
    c13_yiy0_t = (double *)(&c13_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c13_c_iy0), 1, 25, 1, 0) - 1]);
    c13_yix0_t = (double *)(&c13_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c13_c_ix0), 1, 25, 1, 0) - 1]);
    c13_yia0_t = (double *)(&c13_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c13_c_ia0), 1, 25, 1, 0) - 1]);
    dgemv(&c13_TRANSA, &c13_m_t, &c13_n_t, c13_alpha1_t, c13_yia0_t, &c13_lda_t,
          c13_yix0_t, &c13_incx_t, c13_beta1_t, c13_yiy0_t, &c13_incy_t);
  }
}

static void init_dsm_address_info(SFc13_IMM_UKFInstanceStruct *chartInstance)
{
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c13_IMM_UKF_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2039094429U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3591841653U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3072182226U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(646454664U);
}

mxArray *sf_c13_IMM_UKF_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("yrEnNCvpLcvz83wAzYO5rF");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(5);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(5);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(5);
      pr[1] = (double)(5);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c13_IMM_UKF_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c13_IMM_UKF_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c13_IMM_UKF(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x7'type','srcId','name','auxInfo'{{M[1],M[10],T\"Lambda2\",},{M[1],M[16],T\"P_out\",},{M[1],M[9],T\"xhat2\",},{M[4],M[0],T\"P2\",S'l','i','p'{{M1x2[825 827],M[0],}}},{M[4],M[0],T\"V2\",S'l','i','p'{{M1x2[831 833],M[0],}}},{M[4],M[0],T\"W2\",S'l','i','p'{{M1x2[828 830],M[0],}}},{M[8],M[0],T\"is_active_c13_IMM_UKF\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 7, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c13_IMM_UKF_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc13_IMM_UKFInstanceStruct *chartInstance;
    chartInstance = (SFc13_IMM_UKFInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _IMM_UKFMachineNumber_,
           13,
           1,
           1,
           5,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           ssGetPath(S),
           (void *)S);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          init_script_number_translation(_IMM_UKFMachineNumber_,
            chartInstance->chartNumber);
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_IMM_UKFMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _IMM_UKFMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"xhat20");
          _SFD_SET_DATA_PROPS(1,1,1,0,"meas");
          _SFD_SET_DATA_PROPS(2,2,0,1,"xhat2");
          _SFD_SET_DATA_PROPS(3,2,0,1,"Lambda2");
          _SFD_SET_DATA_PROPS(4,2,0,1,"P_out");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,1,0,0,0,3,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,2825);
        _SFD_CV_INIT_EML_IF(0,1,0,834,848,-1,1448);
        _SFD_CV_INIT_EML_FOR(0,1,0,1615,1629,2072);
        _SFD_CV_INIT_EML_FOR(0,1,1,2123,2137,2204);
        _SFD_CV_INIT_EML_FOR(0,1,2,2300,2314,2441);

        {
          unsigned int dimVector[1];
          dimVector[0]= 5;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c13_f_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c13_g_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 5;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c13_f_sf_marshallOut,(MexInFcnForType)
            c13_f_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c13_e_sf_marshallOut,(MexInFcnForType)
          c13_e_sf_marshallIn);

        {
          unsigned int dimVector[2];
          dimVector[0]= 5;
          dimVector[1]= 5;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c13_d_sf_marshallOut,(MexInFcnForType)
            c13_d_sf_marshallIn);
        }

        {
          real_T *c13_Lambda2;
          real_T (*c13_xhat20)[5];
          real_T (*c13_meas)[2];
          real_T (*c13_xhat2)[5];
          real_T (*c13_P_out)[25];
          c13_P_out = (real_T (*)[25])ssGetOutputPortSignal(chartInstance->S, 3);
          c13_Lambda2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
          c13_xhat2 = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
          c13_meas = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 1);
          c13_xhat20 = (real_T (*)[5])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c13_xhat20);
          _SFD_SET_DATA_VALUE_PTR(1U, *c13_meas);
          _SFD_SET_DATA_VALUE_PTR(2U, *c13_xhat2);
          _SFD_SET_DATA_VALUE_PTR(3U, c13_Lambda2);
          _SFD_SET_DATA_VALUE_PTR(4U, *c13_P_out);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _IMM_UKFMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "3o002jr5JvntFmTCrpbBUH";
}

static void sf_opaque_initialize_c13_IMM_UKF(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc13_IMM_UKFInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c13_IMM_UKF((SFc13_IMM_UKFInstanceStruct*) chartInstanceVar);
  initialize_c13_IMM_UKF((SFc13_IMM_UKFInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c13_IMM_UKF(void *chartInstanceVar)
{
  enable_c13_IMM_UKF((SFc13_IMM_UKFInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c13_IMM_UKF(void *chartInstanceVar)
{
  disable_c13_IMM_UKF((SFc13_IMM_UKFInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c13_IMM_UKF(void *chartInstanceVar)
{
  sf_c13_IMM_UKF((SFc13_IMM_UKFInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c13_IMM_UKF(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c13_IMM_UKF((SFc13_IMM_UKFInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c13_IMM_UKF();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c13_IMM_UKF(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c13_IMM_UKF();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c13_IMM_UKF((SFc13_IMM_UKFInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c13_IMM_UKF(SimStruct* S)
{
  return sf_internal_get_sim_state_c13_IMM_UKF(S);
}

static void sf_opaque_set_sim_state_c13_IMM_UKF(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c13_IMM_UKF(S, st);
}

static void sf_opaque_terminate_c13_IMM_UKF(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc13_IMM_UKFInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_IMM_UKF_optimization_info();
    }

    finalize_c13_IMM_UKF((SFc13_IMM_UKFInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc13_IMM_UKF((SFc13_IMM_UKFInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c13_IMM_UKF(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c13_IMM_UKF((SFc13_IMM_UKFInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c13_IMM_UKF(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_IMM_UKF_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      13);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,13,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,13,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,13);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,13,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,13,3);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=3; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 2; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,13);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1292288340U));
  ssSetChecksum1(S,(733839598U));
  ssSetChecksum2(S,(2001602410U));
  ssSetChecksum3(S,(1566984553U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c13_IMM_UKF(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c13_IMM_UKF(SimStruct *S)
{
  SFc13_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc13_IMM_UKFInstanceStruct *)utMalloc(sizeof
    (SFc13_IMM_UKFInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc13_IMM_UKFInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c13_IMM_UKF;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c13_IMM_UKF;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c13_IMM_UKF;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c13_IMM_UKF;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c13_IMM_UKF;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c13_IMM_UKF;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c13_IMM_UKF;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c13_IMM_UKF;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c13_IMM_UKF;
  chartInstance->chartInfo.mdlStart = mdlStart_c13_IMM_UKF;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c13_IMM_UKF;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->S = S;
  ssSetUserData(S,(void *)(&(chartInstance->chartInfo)));/* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c13_IMM_UKF_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c13_IMM_UKF(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c13_IMM_UKF(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c13_IMM_UKF(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c13_IMM_UKF_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
