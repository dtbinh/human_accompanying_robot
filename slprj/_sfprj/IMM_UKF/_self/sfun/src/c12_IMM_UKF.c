/* Include files */

#include <stddef.h>
#include "blas.h"
#include "IMM_UKF_sfun.h"
#include "c12_IMM_UKF.h"
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
static const char * c12_debug_family_names[18] = { "L", "alpha", "kappa", "beta",
  "lambda", "c", "Wm", "Wc", "sP", "chi_p", "chi_m_p", "i", "nargin", "nargout",
  "x", "P", "T1", "x_m_p" };

/* Function Declarations */
static void initialize_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance);
static void initialize_params_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct
  *chartInstance);
static void enable_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance);
static void disable_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance);
static void c12_update_debugger_state_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct
  *chartInstance);
static void set_sim_state_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_st);
static void finalize_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance);
static void sf_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance);
static void c12_chartstep_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance);
static void initSimStructsc12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c12_machineNumber, uint32_T
  c12_chartNumber);
static const mxArray *c12_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData);
static void c12_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_x_m_p, const char_T *c12_identifier, real_T c12_y[5]);
static void c12_b_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId, real_T c12_y[5]);
static void c12_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData);
static const mxArray *c12_b_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData);
static const mxArray *c12_c_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData);
static real_T c12_c_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId);
static void c12_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData);
static const mxArray *c12_d_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData);
static void c12_d_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId, real_T c12_y[55]);
static void c12_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData);
static void c12_e_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId, real_T c12_y[25]);
static void c12_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData);
static const mxArray *c12_e_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData);
static void c12_f_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId, real_T c12_y[11]);
static void c12_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData);
static void c12_info_helper(const mxArray **c12_info);
static const mxArray *c12_emlrt_marshallOut(char * c12_u);
static const mxArray *c12_b_emlrt_marshallOut(uint32_T c12_u);
static void c12_b_info_helper(const mxArray **c12_info);
static real_T c12_mpower(SFc12_IMM_UKFInstanceStruct *chartInstance, real_T
  c12_a);
static void c12_eml_scalar_eg(SFc12_IMM_UKFInstanceStruct *chartInstance);
static void c12_eml_error(SFc12_IMM_UKFInstanceStruct *chartInstance);
static void c12_b_eml_scalar_eg(SFc12_IMM_UKFInstanceStruct *chartInstance);
static void c12_eml_matlab_zpotrf(SFc12_IMM_UKFInstanceStruct *chartInstance,
  real_T c12_A[25], real_T c12_b_A[25], int32_T *c12_info);
static void c12_check_forloop_overflow_error(SFc12_IMM_UKFInstanceStruct
  *chartInstance, boolean_T c12_overflow);
static void c12_eml_xgemv(SFc12_IMM_UKFInstanceStruct *chartInstance, int32_T
  c12_m, int32_T c12_n, int32_T c12_ia0, int32_T c12_ix0, real_T c12_y[25],
  int32_T c12_iy0, real_T c12_b_y[25]);
static void c12_b_eml_error(SFc12_IMM_UKFInstanceStruct *chartInstance);
static void c12_c_eml_scalar_eg(SFc12_IMM_UKFInstanceStruct *chartInstance);
static const mxArray *c12_f_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData);
static int32_T c12_g_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId);
static void c12_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData);
static uint8_T c12_h_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_b_is_active_c12_IMM_UKF, const char_T *c12_identifier);
static uint8_T c12_i_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId);
static int32_T c12_b_eml_matlab_zpotrf(SFc12_IMM_UKFInstanceStruct
  *chartInstance, real_T c12_A[25]);
static void c12_b_eml_xgemv(SFc12_IMM_UKFInstanceStruct *chartInstance, int32_T
  c12_m, int32_T c12_n, int32_T c12_ia0, int32_T c12_ix0, real_T c12_y[25],
  int32_T c12_iy0);
static void init_dsm_address_info(SFc12_IMM_UKFInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance)
{
  chartInstance->c12_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c12_is_active_c12_IMM_UKF = 0U;
}

static void initialize_params_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct
  *chartInstance)
{
}

static void enable_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c12_update_debugger_state_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct
  *chartInstance)
{
}

static const mxArray *get_sim_state_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct
  *chartInstance)
{
  const mxArray *c12_st;
  const mxArray *c12_y = NULL;
  int32_T c12_i0;
  real_T c12_u[5];
  const mxArray *c12_b_y = NULL;
  uint8_T c12_hoistedGlobal;
  uint8_T c12_b_u;
  const mxArray *c12_c_y = NULL;
  real_T (*c12_x_m_p)[5];
  c12_x_m_p = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
  c12_st = NULL;
  c12_st = NULL;
  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_createcellarray(2), FALSE);
  for (c12_i0 = 0; c12_i0 < 5; c12_i0++) {
    c12_u[c12_i0] = (*c12_x_m_p)[c12_i0];
  }

  c12_b_y = NULL;
  sf_mex_assign(&c12_b_y, sf_mex_create("y", c12_u, 0, 0U, 1U, 0U, 1, 5), FALSE);
  sf_mex_setcell(c12_y, 0, c12_b_y);
  c12_hoistedGlobal = chartInstance->c12_is_active_c12_IMM_UKF;
  c12_b_u = c12_hoistedGlobal;
  c12_c_y = NULL;
  sf_mex_assign(&c12_c_y, sf_mex_create("y", &c12_b_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c12_y, 1, c12_c_y);
  sf_mex_assign(&c12_st, c12_y, FALSE);
  return c12_st;
}

static void set_sim_state_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_st)
{
  const mxArray *c12_u;
  real_T c12_dv0[5];
  int32_T c12_i1;
  real_T (*c12_x_m_p)[5];
  c12_x_m_p = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c12_doneDoubleBufferReInit = TRUE;
  c12_u = sf_mex_dup(c12_st);
  c12_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c12_u, 0)),
                       "x_m_p", c12_dv0);
  for (c12_i1 = 0; c12_i1 < 5; c12_i1++) {
    (*c12_x_m_p)[c12_i1] = c12_dv0[c12_i1];
  }

  chartInstance->c12_is_active_c12_IMM_UKF = c12_h_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c12_u, 1)),
     "is_active_c12_IMM_UKF");
  sf_mex_destroy(&c12_u);
  c12_update_debugger_state_c12_IMM_UKF(chartInstance);
  sf_mex_destroy(&c12_st);
}

static void finalize_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance)
{
}

static void sf_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance)
{
  int32_T c12_i2;
  int32_T c12_i3;
  int32_T c12_i4;
  real_T *c12_T1;
  real_T (*c12_x_m_p)[5];
  real_T (*c12_P)[25];
  real_T (*c12_x)[5];
  c12_x_m_p = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
  c12_T1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c12_P = (real_T (*)[25])ssGetInputPortSignal(chartInstance->S, 1);
  c12_x = (real_T (*)[5])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 5U, chartInstance->c12_sfEvent);
  for (c12_i2 = 0; c12_i2 < 5; c12_i2++) {
    _SFD_DATA_RANGE_CHECK((*c12_x)[c12_i2], 0U);
  }

  for (c12_i3 = 0; c12_i3 < 25; c12_i3++) {
    _SFD_DATA_RANGE_CHECK((*c12_P)[c12_i3], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c12_T1, 2U);
  for (c12_i4 = 0; c12_i4 < 5; c12_i4++) {
    _SFD_DATA_RANGE_CHECK((*c12_x_m_p)[c12_i4], 3U);
  }

  chartInstance->c12_sfEvent = CALL_EVENT;
  c12_chartstep_c12_IMM_UKF(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_IMM_UKFMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c12_chartstep_c12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance)
{
  real_T c12_hoistedGlobal;
  int32_T c12_i5;
  real_T c12_x[5];
  int32_T c12_i6;
  real_T c12_P[25];
  real_T c12_T1;
  uint32_T c12_debug_family_var_map[18];
  real_T c12_L;
  real_T c12_alpha;
  real_T c12_kappa;
  real_T c12_beta;
  real_T c12_lambda;
  real_T c12_c;
  real_T c12_Wm[11];
  real_T c12_Wc[11];
  real_T c12_sP[25];
  real_T c12_chi_p[55];
  real_T c12_chi_m_p[55];
  real_T c12_i;
  real_T c12_nargin = 3.0;
  real_T c12_nargout = 1.0;
  real_T c12_x_m_p[5];
  int32_T c12_i7;
  int32_T c12_i8;
  int32_T c12_i9;
  real_T c12_A[25];
  int32_T c12_j;
  int32_T c12_b_j;
  int32_T c12_info;
  int32_T c12_b_info;
  int32_T c12_c_info;
  int32_T c12_d_info;
  int32_T c12_jmax;
  int32_T c12_a;
  int32_T c12_b_jmax;
  int32_T c12_b;
  int32_T c12_b_b;
  boolean_T c12_overflow;
  int32_T c12_c_j;
  int32_T c12_b_a;
  int32_T c12_i10;
  int32_T c12_c_b;
  int32_T c12_d_b;
  boolean_T c12_b_overflow;
  int32_T c12_b_i;
  int32_T c12_c_i;
  int32_T c12_i11;
  int32_T c12_i12;
  real_T c12_c_a[5];
  int32_T c12_i13;
  int32_T c12_i14;
  int32_T c12_i15;
  int32_T c12_i16;
  real_T c12_e_b[25];
  int32_T c12_i17;
  int32_T c12_i18;
  int32_T c12_i19;
  int32_T c12_i20;
  int32_T c12_i21;
  real_T c12_y[25];
  int32_T c12_i22;
  real_T c12_f_b[25];
  int32_T c12_i23;
  int32_T c12_i24;
  int32_T c12_i25;
  int32_T c12_i26;
  int32_T c12_i27;
  int32_T c12_i28;
  int32_T c12_i29;
  int32_T c12_i30;
  int32_T c12_i31;
  int32_T c12_d_i;
  real_T c12_d_a;
  real_T c12_g_b;
  real_T c12_b_y;
  real_T c12_b_x;
  real_T c12_c_x;
  real_T c12_e_a;
  real_T c12_h_b;
  real_T c12_c_y;
  real_T c12_b_A;
  real_T c12_B;
  real_T c12_d_x;
  real_T c12_d_y;
  real_T c12_e_x;
  real_T c12_e_y;
  real_T c12_f_y;
  real_T c12_f_a;
  real_T c12_i_b;
  real_T c12_g_y;
  real_T c12_f_x;
  real_T c12_g_x;
  real_T c12_g_a;
  real_T c12_j_b;
  real_T c12_h_y;
  real_T c12_c_A;
  real_T c12_b_B;
  real_T c12_h_x;
  real_T c12_i_y;
  real_T c12_i_x;
  real_T c12_j_y;
  real_T c12_k_y;
  real_T c12_h_a;
  real_T c12_k_b;
  real_T c12_l_y;
  real_T c12_j_x;
  real_T c12_k_x;
  real_T c12_i_a;
  real_T c12_l_b;
  real_T c12_m_y;
  real_T c12_j_a;
  real_T c12_m_b;
  real_T c12_n_y;
  real_T c12_l_x;
  real_T c12_m_x;
  real_T c12_k_a;
  real_T c12_n_b;
  real_T c12_o_y;
  real_T c12_l_a;
  real_T c12_o_b;
  real_T c12_p_y;
  real_T c12_n_x;
  real_T c12_o_x;
  real_T c12_m_a;
  real_T c12_p_b;
  real_T c12_q_y;
  real_T c12_d_A;
  real_T c12_c_B;
  real_T c12_p_x;
  real_T c12_r_y;
  real_T c12_q_x;
  real_T c12_s_y;
  real_T c12_t_y;
  real_T c12_n_a;
  real_T c12_q_b;
  real_T c12_u_y;
  real_T c12_r_x;
  real_T c12_s_x;
  real_T c12_o_a;
  real_T c12_r_b;
  real_T c12_v_y;
  real_T c12_e_A;
  real_T c12_d_B;
  real_T c12_t_x;
  real_T c12_w_y;
  real_T c12_u_x;
  real_T c12_x_y;
  real_T c12_y_y;
  real_T c12_p_a;
  real_T c12_s_b;
  real_T c12_ab_y;
  real_T c12_v_x;
  real_T c12_w_x;
  real_T c12_q_a;
  real_T c12_t_b;
  real_T c12_bb_y;
  real_T c12_r_a;
  real_T c12_u_b;
  real_T c12_cb_y;
  real_T c12_x_x;
  real_T c12_y_x;
  real_T c12_s_a;
  real_T c12_v_b;
  real_T c12_db_y;
  int32_T c12_i32;
  real_T c12_t_a[55];
  int32_T c12_i33;
  real_T c12_w_b[11];
  int32_T c12_i34;
  int32_T c12_i35;
  int32_T c12_i36;
  int32_T c12_i37;
  int32_T c12_i38;
  int32_T c12_i39;
  int32_T c12_i40;
  int32_T c12_i41;
  int32_T c12_i42;
  int32_T c12_i43;
  real_T (*c12_b_x_m_p)[5];
  real_T *c12_b_T1;
  real_T (*c12_b_P)[25];
  real_T (*c12_ab_x)[5];
  c12_b_x_m_p = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
  c12_b_T1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c12_b_P = (real_T (*)[25])ssGetInputPortSignal(chartInstance->S, 1);
  c12_ab_x = (real_T (*)[5])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 5U, chartInstance->c12_sfEvent);
  c12_hoistedGlobal = *c12_b_T1;
  for (c12_i5 = 0; c12_i5 < 5; c12_i5++) {
    c12_x[c12_i5] = (*c12_ab_x)[c12_i5];
  }

  for (c12_i6 = 0; c12_i6 < 25; c12_i6++) {
    c12_P[c12_i6] = (*c12_b_P)[c12_i6];
  }

  c12_T1 = c12_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 18U, 18U, c12_debug_family_names,
    c12_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c12_L, 0U, c12_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c12_alpha, 1U, c12_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c12_kappa, 2U, c12_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c12_beta, 3U, c12_b_sf_marshallOut,
    c12_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c12_lambda, 4U, c12_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c12_c, 5U, c12_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c12_Wm, 6U, c12_e_sf_marshallOut,
    c12_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c12_Wc, 7U, c12_e_sf_marshallOut,
    c12_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c12_sP, 8U, c12_c_sf_marshallOut,
    c12_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c12_chi_p, 9U, c12_d_sf_marshallOut,
    c12_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c12_chi_m_p, 10U, c12_d_sf_marshallOut,
    c12_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c12_i, 11U, c12_b_sf_marshallOut,
    c12_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c12_nargin, 12U, c12_b_sf_marshallOut,
    c12_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c12_nargout, 13U, c12_b_sf_marshallOut,
    c12_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c12_x, 14U, c12_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c12_P, 15U, c12_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c12_T1, 16U, c12_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c12_x_m_p, 17U, c12_sf_marshallOut,
    c12_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 2);
  c12_L = 5.0;
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 4);
  c12_alpha = 0.001;
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 5);
  c12_kappa = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 6);
  c12_beta = 2.0;
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 8);
  c12_lambda = -4.999995;
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 9);
  c12_c = 4.9999999998107114E-6;
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 14);
  for (c12_i7 = 0; c12_i7 < 11; c12_i7++) {
    c12_Wm[c12_i7] = 100000.00000378577;
  }

  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 15);
  for (c12_i8 = 0; c12_i8 < 11; c12_i8++) {
    c12_Wc[c12_i8] = c12_Wm[c12_i8];
  }

  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 16);
  c12_Wm[0] = -999999.00003785768;
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 17);
  c12_Wc[0] = (-999998.00003785768 - c12_mpower(chartInstance, 0.001)) +
    c12_beta;
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 18);
  c12_c = 0.0022360679774574635;
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 20);
  for (c12_i9 = 0; c12_i9 < 25; c12_i9++) {
    c12_A[c12_i9] = c12_P[c12_i9];
  }

  for (c12_j = 1; c12_j < 6; c12_j++) {
    c12_b_j = c12_j;
  }

  c12_info = c12_b_eml_matlab_zpotrf(chartInstance, c12_A);
  c12_b_info = c12_info;
  c12_c_info = c12_b_info;
  c12_d_info = c12_c_info;
  if (c12_d_info == 0) {
    c12_jmax = 6;
  } else {
    c12_b_eml_error(chartInstance);
    c12_a = c12_d_info;
    c12_jmax = c12_a;
  }

  c12_b_jmax = c12_jmax - 1;
  c12_b = c12_b_jmax;
  c12_b_b = c12_b;
  if (2 > c12_b_b) {
    c12_overflow = FALSE;
  } else {
    c12_overflow = (c12_b_b > 2147483646);
  }

  if (c12_overflow) {
    c12_check_forloop_overflow_error(chartInstance, c12_overflow);
  }

  for (c12_c_j = 2; c12_c_j <= c12_b_jmax; c12_c_j++) {
    c12_b_j = c12_c_j;
    c12_b_a = c12_b_j - 1;
    c12_i10 = c12_b_a;
    c12_c_b = c12_i10;
    c12_d_b = c12_c_b;
    if (1 > c12_d_b) {
      c12_b_overflow = FALSE;
    } else {
      c12_b_overflow = (c12_d_b > 2147483646);
    }

    if (c12_b_overflow) {
      c12_check_forloop_overflow_error(chartInstance, c12_b_overflow);
    }

    for (c12_b_i = 1; c12_b_i <= c12_i10; c12_b_i++) {
      c12_c_i = c12_b_i;
      c12_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
               (real_T)c12_c_i), 1, 5, 1, 0) + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
               "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c12_b_j), 1, 5, 2, 0)
              - 1)) - 1] = 0.0;
    }
  }

  for (c12_i11 = 0; c12_i11 < 25; c12_i11++) {
    c12_sP[c12_i11] = c12_A[c12_i11];
  }

  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 21);
  for (c12_i12 = 0; c12_i12 < 5; c12_i12++) {
    c12_c_a[c12_i12] = c12_x[c12_i12];
  }

  for (c12_i13 = 0; c12_i13 < 5; c12_i13++) {
    c12_i14 = 0;
    for (c12_i15 = 0; c12_i15 < 5; c12_i15++) {
      c12_A[c12_i14 + c12_i13] = c12_c_a[c12_i13];
      c12_i14 += 5;
    }
  }

  for (c12_i16 = 0; c12_i16 < 25; c12_i16++) {
    c12_e_b[c12_i16] = c12_sP[c12_i16];
  }

  for (c12_i17 = 0; c12_i17 < 25; c12_i17++) {
    c12_e_b[c12_i17] *= 0.0022360679774574635;
  }

  for (c12_i18 = 0; c12_i18 < 5; c12_i18++) {
    c12_c_a[c12_i18] = c12_x[c12_i18];
  }

  for (c12_i19 = 0; c12_i19 < 5; c12_i19++) {
    c12_i20 = 0;
    for (c12_i21 = 0; c12_i21 < 5; c12_i21++) {
      c12_y[c12_i20 + c12_i19] = c12_c_a[c12_i19];
      c12_i20 += 5;
    }
  }

  for (c12_i22 = 0; c12_i22 < 25; c12_i22++) {
    c12_f_b[c12_i22] = c12_sP[c12_i22];
  }

  for (c12_i23 = 0; c12_i23 < 25; c12_i23++) {
    c12_f_b[c12_i23] *= 0.0022360679774574635;
  }

  for (c12_i24 = 0; c12_i24 < 5; c12_i24++) {
    c12_chi_p[c12_i24] = c12_x[c12_i24];
  }

  c12_i25 = 0;
  for (c12_i26 = 0; c12_i26 < 5; c12_i26++) {
    for (c12_i27 = 0; c12_i27 < 5; c12_i27++) {
      c12_chi_p[(c12_i27 + c12_i25) + 5] = c12_A[c12_i27 + c12_i25] +
        c12_e_b[c12_i27 + c12_i25];
    }

    c12_i25 += 5;
  }

  c12_i28 = 0;
  for (c12_i29 = 0; c12_i29 < 5; c12_i29++) {
    for (c12_i30 = 0; c12_i30 < 5; c12_i30++) {
      c12_chi_p[(c12_i30 + c12_i28) + 30] = c12_y[c12_i30 + c12_i28] -
        c12_f_b[c12_i30 + c12_i28];
    }

    c12_i28 += 5;
  }

  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 24);
  for (c12_i31 = 0; c12_i31 < 55; c12_i31++) {
    c12_chi_m_p[c12_i31] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 25);
  c12_i = 1.0;
  c12_d_i = 0;
  while (c12_d_i < 11) {
    c12_i = 1.0 + (real_T)c12_d_i;
    CV_EML_FOR(0, 1, 0, 1);
    _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 26);
    c12_d_a = c12_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_g_b = c12_T1;
    c12_b_y = c12_d_a * c12_g_b;
    c12_b_x = c12_b_y;
    c12_c_x = c12_b_x;
    c12_c_x = muDoubleScalarSin(c12_c_x);
    c12_e_a = c12_chi_p[1 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_h_b = c12_c_x;
    c12_c_y = c12_e_a * c12_h_b;
    c12_b_A = c12_c_y;
    c12_B = c12_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_d_x = c12_b_A;
    c12_d_y = c12_B;
    c12_e_x = c12_d_x;
    c12_e_y = c12_d_y;
    c12_f_y = c12_e_x / c12_e_y;
    c12_f_a = c12_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_i_b = c12_T1;
    c12_g_y = c12_f_a * c12_i_b;
    c12_f_x = c12_g_y;
    c12_g_x = c12_f_x;
    c12_g_x = muDoubleScalarCos(c12_g_x);
    c12_g_a = c12_chi_p[3 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_j_b = 1.0 - c12_g_x;
    c12_h_y = c12_g_a * c12_j_b;
    c12_c_A = c12_h_y;
    c12_b_B = c12_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_h_x = c12_c_A;
    c12_i_y = c12_b_B;
    c12_i_x = c12_h_x;
    c12_j_y = c12_i_y;
    c12_k_y = c12_i_x / c12_j_y;
    c12_chi_m_p[5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m_p",
      (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)] = (c12_chi_p[5
      * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)] + c12_f_y) - c12_k_y;
    _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 27);
    c12_h_a = c12_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_k_b = c12_T1;
    c12_l_y = c12_h_a * c12_k_b;
    c12_j_x = c12_l_y;
    c12_k_x = c12_j_x;
    c12_k_x = muDoubleScalarCos(c12_k_x);
    c12_i_a = c12_k_x;
    c12_l_b = c12_chi_p[1 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_m_y = c12_i_a * c12_l_b;
    c12_j_a = c12_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_m_b = c12_T1;
    c12_n_y = c12_j_a * c12_m_b;
    c12_l_x = c12_n_y;
    c12_m_x = c12_l_x;
    c12_m_x = muDoubleScalarSin(c12_m_x);
    c12_k_a = c12_m_x;
    c12_n_b = c12_chi_p[3 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_o_y = c12_k_a * c12_n_b;
    c12_chi_m_p[1 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m_p",
                          (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0)
                         - 1)] = c12_m_y - c12_o_y;
    _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 28);
    c12_l_a = c12_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_o_b = c12_T1;
    c12_p_y = c12_l_a * c12_o_b;
    c12_n_x = c12_p_y;
    c12_o_x = c12_n_x;
    c12_o_x = muDoubleScalarCos(c12_o_x);
    c12_m_a = c12_chi_p[1 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_p_b = 1.0 - c12_o_x;
    c12_q_y = c12_m_a * c12_p_b;
    c12_d_A = c12_q_y;
    c12_c_B = c12_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_p_x = c12_d_A;
    c12_r_y = c12_c_B;
    c12_q_x = c12_p_x;
    c12_s_y = c12_r_y;
    c12_t_y = c12_q_x / c12_s_y;
    c12_n_a = c12_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_q_b = c12_T1;
    c12_u_y = c12_n_a * c12_q_b;
    c12_r_x = c12_u_y;
    c12_s_x = c12_r_x;
    c12_s_x = muDoubleScalarSin(c12_s_x);
    c12_o_a = c12_chi_p[3 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_r_b = c12_s_x;
    c12_v_y = c12_o_a * c12_r_b;
    c12_e_A = c12_v_y;
    c12_d_B = c12_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_t_x = c12_e_A;
    c12_w_y = c12_d_B;
    c12_u_x = c12_t_x;
    c12_x_y = c12_w_y;
    c12_y_y = c12_u_x / c12_x_y;
    c12_chi_m_p[2 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m_p",
                          (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0)
                         - 1)] = (c12_t_y + c12_chi_p[2 + 5 * ((int32_T)(real_T)
      _SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)_SFD_INTEGER_CHECK("i",
      c12_i), 1, 11, 2, 0) - 1)]) + c12_y_y;
    _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 29);
    c12_p_a = c12_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_s_b = c12_T1;
    c12_ab_y = c12_p_a * c12_s_b;
    c12_v_x = c12_ab_y;
    c12_w_x = c12_v_x;
    c12_w_x = muDoubleScalarSin(c12_w_x);
    c12_q_a = c12_chi_p[1 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_t_b = c12_w_x;
    c12_bb_y = c12_q_a * c12_t_b;
    c12_r_a = c12_chi_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_u_b = c12_T1;
    c12_cb_y = c12_r_a * c12_u_b;
    c12_x_x = c12_cb_y;
    c12_y_x = c12_x_x;
    c12_y_x = muDoubleScalarCos(c12_y_x);
    c12_s_a = c12_chi_p[3 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
      "chi_p", (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0) - 1)];
    c12_v_b = c12_y_x;
    c12_db_y = c12_s_a * c12_v_b;
    c12_chi_m_p[3 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m_p",
                          (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0)
                         - 1)] = c12_bb_y + c12_db_y;
    _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 30);
    c12_chi_m_p[4 + 5 * ((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m_p",
                          (int32_T)_SFD_INTEGER_CHECK("i", c12_i), 1, 11, 2, 0)
                         - 1)] = c12_chi_p[4 + 5 * ((int32_T)(real_T)
      _SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)_SFD_INTEGER_CHECK("i",
      c12_i), 1, 11, 2, 0) - 1)];
    c12_d_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_EML_FOR(0, 1, 0, 0);
  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, 32);
  for (c12_i32 = 0; c12_i32 < 55; c12_i32++) {
    c12_t_a[c12_i32] = c12_chi_m_p[c12_i32];
  }

  for (c12_i33 = 0; c12_i33 < 11; c12_i33++) {
    c12_w_b[c12_i33] = c12_Wm[c12_i33];
  }

  c12_c_eml_scalar_eg(chartInstance);
  c12_c_eml_scalar_eg(chartInstance);
  for (c12_i34 = 0; c12_i34 < 5; c12_i34++) {
    c12_x_m_p[c12_i34] = 0.0;
  }

  for (c12_i35 = 0; c12_i35 < 5; c12_i35++) {
    c12_x_m_p[c12_i35] = 0.0;
  }

  for (c12_i36 = 0; c12_i36 < 5; c12_i36++) {
    c12_c_a[c12_i36] = c12_x_m_p[c12_i36];
  }

  for (c12_i37 = 0; c12_i37 < 5; c12_i37++) {
    c12_x_m_p[c12_i37] = c12_c_a[c12_i37];
  }

  for (c12_i38 = 0; c12_i38 < 5; c12_i38++) {
    c12_c_a[c12_i38] = c12_x_m_p[c12_i38];
  }

  for (c12_i39 = 0; c12_i39 < 5; c12_i39++) {
    c12_x_m_p[c12_i39] = c12_c_a[c12_i39];
  }

  for (c12_i40 = 0; c12_i40 < 5; c12_i40++) {
    c12_x_m_p[c12_i40] = 0.0;
    c12_i41 = 0;
    for (c12_i42 = 0; c12_i42 < 11; c12_i42++) {
      c12_x_m_p[c12_i40] += c12_t_a[c12_i41 + c12_i40] * c12_w_b[c12_i42];
      c12_i41 += 5;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c12_sfEvent, -32);
  _SFD_SYMBOL_SCOPE_POP();
  for (c12_i43 = 0; c12_i43 < 5; c12_i43++) {
    (*c12_b_x_m_p)[c12_i43] = c12_x_m_p[c12_i43];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 5U, chartInstance->c12_sfEvent);
}

static void initSimStructsc12_IMM_UKF(SFc12_IMM_UKFInstanceStruct *chartInstance)
{
}

static void init_script_number_translation(uint32_T c12_machineNumber, uint32_T
  c12_chartNumber)
{
}

static const mxArray *c12_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData)
{
  const mxArray *c12_mxArrayOutData = NULL;
  int32_T c12_i44;
  real_T c12_b_inData[5];
  int32_T c12_i45;
  real_T c12_u[5];
  const mxArray *c12_y = NULL;
  SFc12_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc12_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c12_mxArrayOutData = NULL;
  for (c12_i44 = 0; c12_i44 < 5; c12_i44++) {
    c12_b_inData[c12_i44] = (*(real_T (*)[5])c12_inData)[c12_i44];
  }

  for (c12_i45 = 0; c12_i45 < 5; c12_i45++) {
    c12_u[c12_i45] = c12_b_inData[c12_i45];
  }

  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", c12_u, 0, 0U, 1U, 0U, 1, 5), FALSE);
  sf_mex_assign(&c12_mxArrayOutData, c12_y, FALSE);
  return c12_mxArrayOutData;
}

static void c12_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_x_m_p, const char_T *c12_identifier, real_T c12_y[5])
{
  emlrtMsgIdentifier c12_thisId;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c12_x_m_p), &c12_thisId,
    c12_y);
  sf_mex_destroy(&c12_x_m_p);
}

static void c12_b_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId, real_T c12_y[5])
{
  real_T c12_dv1[5];
  int32_T c12_i46;
  sf_mex_import(c12_parentId, sf_mex_dup(c12_u), c12_dv1, 1, 0, 0U, 1, 0U, 1, 5);
  for (c12_i46 = 0; c12_i46 < 5; c12_i46++) {
    c12_y[c12_i46] = c12_dv1[c12_i46];
  }

  sf_mex_destroy(&c12_u);
}

static void c12_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData)
{
  const mxArray *c12_x_m_p;
  const char_T *c12_identifier;
  emlrtMsgIdentifier c12_thisId;
  real_T c12_y[5];
  int32_T c12_i47;
  SFc12_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc12_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c12_x_m_p = sf_mex_dup(c12_mxArrayInData);
  c12_identifier = c12_varName;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c12_x_m_p), &c12_thisId,
    c12_y);
  sf_mex_destroy(&c12_x_m_p);
  for (c12_i47 = 0; c12_i47 < 5; c12_i47++) {
    (*(real_T (*)[5])c12_outData)[c12_i47] = c12_y[c12_i47];
  }

  sf_mex_destroy(&c12_mxArrayInData);
}

static const mxArray *c12_b_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData)
{
  const mxArray *c12_mxArrayOutData = NULL;
  real_T c12_u;
  const mxArray *c12_y = NULL;
  SFc12_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc12_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c12_mxArrayOutData = NULL;
  c12_u = *(real_T *)c12_inData;
  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", &c12_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c12_mxArrayOutData, c12_y, FALSE);
  return c12_mxArrayOutData;
}

static const mxArray *c12_c_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData)
{
  const mxArray *c12_mxArrayOutData = NULL;
  int32_T c12_i48;
  int32_T c12_i49;
  int32_T c12_i50;
  real_T c12_b_inData[25];
  int32_T c12_i51;
  int32_T c12_i52;
  int32_T c12_i53;
  real_T c12_u[25];
  const mxArray *c12_y = NULL;
  SFc12_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc12_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c12_mxArrayOutData = NULL;
  c12_i48 = 0;
  for (c12_i49 = 0; c12_i49 < 5; c12_i49++) {
    for (c12_i50 = 0; c12_i50 < 5; c12_i50++) {
      c12_b_inData[c12_i50 + c12_i48] = (*(real_T (*)[25])c12_inData)[c12_i50 +
        c12_i48];
    }

    c12_i48 += 5;
  }

  c12_i51 = 0;
  for (c12_i52 = 0; c12_i52 < 5; c12_i52++) {
    for (c12_i53 = 0; c12_i53 < 5; c12_i53++) {
      c12_u[c12_i53 + c12_i51] = c12_b_inData[c12_i53 + c12_i51];
    }

    c12_i51 += 5;
  }

  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", c12_u, 0, 0U, 1U, 0U, 2, 5, 5), FALSE);
  sf_mex_assign(&c12_mxArrayOutData, c12_y, FALSE);
  return c12_mxArrayOutData;
}

static real_T c12_c_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId)
{
  real_T c12_y;
  real_T c12_d0;
  sf_mex_import(c12_parentId, sf_mex_dup(c12_u), &c12_d0, 1, 0, 0U, 0, 0U, 0);
  c12_y = c12_d0;
  sf_mex_destroy(&c12_u);
  return c12_y;
}

static void c12_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData)
{
  const mxArray *c12_nargout;
  const char_T *c12_identifier;
  emlrtMsgIdentifier c12_thisId;
  real_T c12_y;
  SFc12_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc12_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c12_nargout = sf_mex_dup(c12_mxArrayInData);
  c12_identifier = c12_varName;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_y = c12_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c12_nargout),
    &c12_thisId);
  sf_mex_destroy(&c12_nargout);
  *(real_T *)c12_outData = c12_y;
  sf_mex_destroy(&c12_mxArrayInData);
}

static const mxArray *c12_d_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData)
{
  const mxArray *c12_mxArrayOutData = NULL;
  int32_T c12_i54;
  int32_T c12_i55;
  int32_T c12_i56;
  real_T c12_b_inData[55];
  int32_T c12_i57;
  int32_T c12_i58;
  int32_T c12_i59;
  real_T c12_u[55];
  const mxArray *c12_y = NULL;
  SFc12_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc12_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c12_mxArrayOutData = NULL;
  c12_i54 = 0;
  for (c12_i55 = 0; c12_i55 < 11; c12_i55++) {
    for (c12_i56 = 0; c12_i56 < 5; c12_i56++) {
      c12_b_inData[c12_i56 + c12_i54] = (*(real_T (*)[55])c12_inData)[c12_i56 +
        c12_i54];
    }

    c12_i54 += 5;
  }

  c12_i57 = 0;
  for (c12_i58 = 0; c12_i58 < 11; c12_i58++) {
    for (c12_i59 = 0; c12_i59 < 5; c12_i59++) {
      c12_u[c12_i59 + c12_i57] = c12_b_inData[c12_i59 + c12_i57];
    }

    c12_i57 += 5;
  }

  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", c12_u, 0, 0U, 1U, 0U, 2, 5, 11),
                FALSE);
  sf_mex_assign(&c12_mxArrayOutData, c12_y, FALSE);
  return c12_mxArrayOutData;
}

static void c12_d_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId, real_T c12_y[55])
{
  real_T c12_dv2[55];
  int32_T c12_i60;
  sf_mex_import(c12_parentId, sf_mex_dup(c12_u), c12_dv2, 1, 0, 0U, 1, 0U, 2, 5,
                11);
  for (c12_i60 = 0; c12_i60 < 55; c12_i60++) {
    c12_y[c12_i60] = c12_dv2[c12_i60];
  }

  sf_mex_destroy(&c12_u);
}

static void c12_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData)
{
  const mxArray *c12_chi_m_p;
  const char_T *c12_identifier;
  emlrtMsgIdentifier c12_thisId;
  real_T c12_y[55];
  int32_T c12_i61;
  int32_T c12_i62;
  int32_T c12_i63;
  SFc12_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc12_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c12_chi_m_p = sf_mex_dup(c12_mxArrayInData);
  c12_identifier = c12_varName;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c12_chi_m_p), &c12_thisId,
    c12_y);
  sf_mex_destroy(&c12_chi_m_p);
  c12_i61 = 0;
  for (c12_i62 = 0; c12_i62 < 11; c12_i62++) {
    for (c12_i63 = 0; c12_i63 < 5; c12_i63++) {
      (*(real_T (*)[55])c12_outData)[c12_i63 + c12_i61] = c12_y[c12_i63 +
        c12_i61];
    }

    c12_i61 += 5;
  }

  sf_mex_destroy(&c12_mxArrayInData);
}

static void c12_e_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId, real_T c12_y[25])
{
  real_T c12_dv3[25];
  int32_T c12_i64;
  sf_mex_import(c12_parentId, sf_mex_dup(c12_u), c12_dv3, 1, 0, 0U, 1, 0U, 2, 5,
                5);
  for (c12_i64 = 0; c12_i64 < 25; c12_i64++) {
    c12_y[c12_i64] = c12_dv3[c12_i64];
  }

  sf_mex_destroy(&c12_u);
}

static void c12_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData)
{
  const mxArray *c12_sP;
  const char_T *c12_identifier;
  emlrtMsgIdentifier c12_thisId;
  real_T c12_y[25];
  int32_T c12_i65;
  int32_T c12_i66;
  int32_T c12_i67;
  SFc12_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc12_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c12_sP = sf_mex_dup(c12_mxArrayInData);
  c12_identifier = c12_varName;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c12_sP), &c12_thisId, c12_y);
  sf_mex_destroy(&c12_sP);
  c12_i65 = 0;
  for (c12_i66 = 0; c12_i66 < 5; c12_i66++) {
    for (c12_i67 = 0; c12_i67 < 5; c12_i67++) {
      (*(real_T (*)[25])c12_outData)[c12_i67 + c12_i65] = c12_y[c12_i67 +
        c12_i65];
    }

    c12_i65 += 5;
  }

  sf_mex_destroy(&c12_mxArrayInData);
}

static const mxArray *c12_e_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData)
{
  const mxArray *c12_mxArrayOutData = NULL;
  int32_T c12_i68;
  real_T c12_b_inData[11];
  int32_T c12_i69;
  real_T c12_u[11];
  const mxArray *c12_y = NULL;
  SFc12_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc12_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c12_mxArrayOutData = NULL;
  for (c12_i68 = 0; c12_i68 < 11; c12_i68++) {
    c12_b_inData[c12_i68] = (*(real_T (*)[11])c12_inData)[c12_i68];
  }

  for (c12_i69 = 0; c12_i69 < 11; c12_i69++) {
    c12_u[c12_i69] = c12_b_inData[c12_i69];
  }

  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", c12_u, 0, 0U, 1U, 0U, 1, 11), FALSE);
  sf_mex_assign(&c12_mxArrayOutData, c12_y, FALSE);
  return c12_mxArrayOutData;
}

static void c12_f_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId, real_T c12_y[11])
{
  real_T c12_dv4[11];
  int32_T c12_i70;
  sf_mex_import(c12_parentId, sf_mex_dup(c12_u), c12_dv4, 1, 0, 0U, 1, 0U, 1, 11);
  for (c12_i70 = 0; c12_i70 < 11; c12_i70++) {
    c12_y[c12_i70] = c12_dv4[c12_i70];
  }

  sf_mex_destroy(&c12_u);
}

static void c12_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData)
{
  const mxArray *c12_Wc;
  const char_T *c12_identifier;
  emlrtMsgIdentifier c12_thisId;
  real_T c12_y[11];
  int32_T c12_i71;
  SFc12_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc12_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c12_Wc = sf_mex_dup(c12_mxArrayInData);
  c12_identifier = c12_varName;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c12_Wc), &c12_thisId, c12_y);
  sf_mex_destroy(&c12_Wc);
  for (c12_i71 = 0; c12_i71 < 11; c12_i71++) {
    (*(real_T (*)[11])c12_outData)[c12_i71] = c12_y[c12_i71];
  }

  sf_mex_destroy(&c12_mxArrayInData);
}

const mxArray *sf_c12_IMM_UKF_get_eml_resolved_functions_info(void)
{
  const mxArray *c12_nameCaptureInfo = NULL;
  c12_nameCaptureInfo = NULL;
  sf_mex_assign(&c12_nameCaptureInfo, sf_mex_createstruct("structure", 2, 94, 1),
                FALSE);
  c12_info_helper(&c12_nameCaptureInfo);
  c12_b_info_helper(&c12_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c12_nameCaptureInfo);
  return c12_nameCaptureInfo;
}

static void c12_info_helper(const mxArray **c12_info)
{
  const mxArray *c12_rhs0 = NULL;
  const mxArray *c12_lhs0 = NULL;
  const mxArray *c12_rhs1 = NULL;
  const mxArray *c12_lhs1 = NULL;
  const mxArray *c12_rhs2 = NULL;
  const mxArray *c12_lhs2 = NULL;
  const mxArray *c12_rhs3 = NULL;
  const mxArray *c12_lhs3 = NULL;
  const mxArray *c12_rhs4 = NULL;
  const mxArray *c12_lhs4 = NULL;
  const mxArray *c12_rhs5 = NULL;
  const mxArray *c12_lhs5 = NULL;
  const mxArray *c12_rhs6 = NULL;
  const mxArray *c12_lhs6 = NULL;
  const mxArray *c12_rhs7 = NULL;
  const mxArray *c12_lhs7 = NULL;
  const mxArray *c12_rhs8 = NULL;
  const mxArray *c12_lhs8 = NULL;
  const mxArray *c12_rhs9 = NULL;
  const mxArray *c12_lhs9 = NULL;
  const mxArray *c12_rhs10 = NULL;
  const mxArray *c12_lhs10 = NULL;
  const mxArray *c12_rhs11 = NULL;
  const mxArray *c12_lhs11 = NULL;
  const mxArray *c12_rhs12 = NULL;
  const mxArray *c12_lhs12 = NULL;
  const mxArray *c12_rhs13 = NULL;
  const mxArray *c12_lhs13 = NULL;
  const mxArray *c12_rhs14 = NULL;
  const mxArray *c12_lhs14 = NULL;
  const mxArray *c12_rhs15 = NULL;
  const mxArray *c12_lhs15 = NULL;
  const mxArray *c12_rhs16 = NULL;
  const mxArray *c12_lhs16 = NULL;
  const mxArray *c12_rhs17 = NULL;
  const mxArray *c12_lhs17 = NULL;
  const mxArray *c12_rhs18 = NULL;
  const mxArray *c12_lhs18 = NULL;
  const mxArray *c12_rhs19 = NULL;
  const mxArray *c12_lhs19 = NULL;
  const mxArray *c12_rhs20 = NULL;
  const mxArray *c12_lhs20 = NULL;
  const mxArray *c12_rhs21 = NULL;
  const mxArray *c12_lhs21 = NULL;
  const mxArray *c12_rhs22 = NULL;
  const mxArray *c12_lhs22 = NULL;
  const mxArray *c12_rhs23 = NULL;
  const mxArray *c12_lhs23 = NULL;
  const mxArray *c12_rhs24 = NULL;
  const mxArray *c12_lhs24 = NULL;
  const mxArray *c12_rhs25 = NULL;
  const mxArray *c12_lhs25 = NULL;
  const mxArray *c12_rhs26 = NULL;
  const mxArray *c12_lhs26 = NULL;
  const mxArray *c12_rhs27 = NULL;
  const mxArray *c12_lhs27 = NULL;
  const mxArray *c12_rhs28 = NULL;
  const mxArray *c12_lhs28 = NULL;
  const mxArray *c12_rhs29 = NULL;
  const mxArray *c12_lhs29 = NULL;
  const mxArray *c12_rhs30 = NULL;
  const mxArray *c12_lhs30 = NULL;
  const mxArray *c12_rhs31 = NULL;
  const mxArray *c12_lhs31 = NULL;
  const mxArray *c12_rhs32 = NULL;
  const mxArray *c12_lhs32 = NULL;
  const mxArray *c12_rhs33 = NULL;
  const mxArray *c12_lhs33 = NULL;
  const mxArray *c12_rhs34 = NULL;
  const mxArray *c12_lhs34 = NULL;
  const mxArray *c12_rhs35 = NULL;
  const mxArray *c12_lhs35 = NULL;
  const mxArray *c12_rhs36 = NULL;
  const mxArray *c12_lhs36 = NULL;
  const mxArray *c12_rhs37 = NULL;
  const mxArray *c12_lhs37 = NULL;
  const mxArray *c12_rhs38 = NULL;
  const mxArray *c12_lhs38 = NULL;
  const mxArray *c12_rhs39 = NULL;
  const mxArray *c12_lhs39 = NULL;
  const mxArray *c12_rhs40 = NULL;
  const mxArray *c12_lhs40 = NULL;
  const mxArray *c12_rhs41 = NULL;
  const mxArray *c12_lhs41 = NULL;
  const mxArray *c12_rhs42 = NULL;
  const mxArray *c12_lhs42 = NULL;
  const mxArray *c12_rhs43 = NULL;
  const mxArray *c12_lhs43 = NULL;
  const mxArray *c12_rhs44 = NULL;
  const mxArray *c12_lhs44 = NULL;
  const mxArray *c12_rhs45 = NULL;
  const mxArray *c12_lhs45 = NULL;
  const mxArray *c12_rhs46 = NULL;
  const mxArray *c12_lhs46 = NULL;
  const mxArray *c12_rhs47 = NULL;
  const mxArray *c12_lhs47 = NULL;
  const mxArray *c12_rhs48 = NULL;
  const mxArray *c12_lhs48 = NULL;
  const mxArray *c12_rhs49 = NULL;
  const mxArray *c12_lhs49 = NULL;
  const mxArray *c12_rhs50 = NULL;
  const mxArray *c12_lhs50 = NULL;
  const mxArray *c12_rhs51 = NULL;
  const mxArray *c12_lhs51 = NULL;
  const mxArray *c12_rhs52 = NULL;
  const mxArray *c12_lhs52 = NULL;
  const mxArray *c12_rhs53 = NULL;
  const mxArray *c12_lhs53 = NULL;
  const mxArray *c12_rhs54 = NULL;
  const mxArray *c12_lhs54 = NULL;
  const mxArray *c12_rhs55 = NULL;
  const mxArray *c12_lhs55 = NULL;
  const mxArray *c12_rhs56 = NULL;
  const mxArray *c12_lhs56 = NULL;
  const mxArray *c12_rhs57 = NULL;
  const mxArray *c12_lhs57 = NULL;
  const mxArray *c12_rhs58 = NULL;
  const mxArray *c12_lhs58 = NULL;
  const mxArray *c12_rhs59 = NULL;
  const mxArray *c12_lhs59 = NULL;
  const mxArray *c12_rhs60 = NULL;
  const mxArray *c12_lhs60 = NULL;
  const mxArray *c12_rhs61 = NULL;
  const mxArray *c12_lhs61 = NULL;
  const mxArray *c12_rhs62 = NULL;
  const mxArray *c12_lhs62 = NULL;
  const mxArray *c12_rhs63 = NULL;
  const mxArray *c12_lhs63 = NULL;
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("mpower"), "name", "name", 0);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c12_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c12_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("ismatrix"), "name", "name",
                  2);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1331337258U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c12_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("power"), "name", "name", 3);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742680U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c12_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c12_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 5);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 5);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c12_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 6);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 6);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1358218540U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c12_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 7);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("floor"), "name", "name", 7);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742654U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c12_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 8);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c12_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 9);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851126U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c12_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 10);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 10);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c12_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 11);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("mtimes"), "name", "name", 11);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c12_rhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 12);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 12);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c12_rhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "context", "context", 13);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("mtimes"), "name", "name", 13);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c12_rhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "context", "context", 14);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("mrdivide"), "name", "name",
                  14);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1373338908U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1319762366U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c12_rhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 15);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("rdivide"), "name", "name",
                  15);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742680U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c12_rhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 16);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c12_rhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 17);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c12_rhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_div"), "name", "name",
                  18);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742666U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c12_rhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "context", "context", 19);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("sqrt"), "name", "name", 19);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1343862786U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c12_rhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_error"), "name", "name",
                  20);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 20);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1343862758U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c12_rhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 21);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 21);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851138U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c12_rhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "context", "context", 22);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("chol"), "name", "name", 22);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m"), "resolved",
                  "resolved", 22);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1344504434U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c12_rhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m"), "context",
                  "context", 23);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_tolower"), "name",
                  "name", 23);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_tolower.m"), "resolved",
                  "resolved", 23);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742670U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c12_rhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 24);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 24);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c12_rhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 25);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("ismatrix"), "name", "name",
                  25);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1331337258U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c12_rhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 26);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 26);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1346542740U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c12_rhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 27);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("intmax"), "name", "name", 27);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1362294282U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c12_rhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 28);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_error"), "name", "name",
                  28);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1343862758U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c12_rhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 29);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_xpotrf"), "name", "name",
                  29);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xpotrf.m"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851208U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c12_rhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xpotrf.m"),
                  "context", "context", 30);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_lapack_xpotrf"), "name",
                  "name", 30);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xpotrf.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851212U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c12_rhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xpotrf.m"),
                  "context", "context", 31);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_matlab_zpotrf"), "name",
                  "name", 31);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851224U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c12_rhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 32);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 32);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c12_rhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 33);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 33);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 33);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c12_rhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 34);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 34);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1346542740U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c12_rhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 35);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 35);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c12_rhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 36);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 36);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c12_rhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 37);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 37);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 37);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851180U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c12_rhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 38);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 38);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c12_rhs38, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs38, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 39);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 39);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 39);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c12_rhs39, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs39, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 40);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 40);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c12_rhs40, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs40, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 41);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_xdotc"), "name", "name",
                  41);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742670U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c12_rhs41, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs41, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"), "context",
                  "context", 42);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 42);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1299109168U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c12_rhs42, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs42, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"), "context",
                  "context", 43);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_xdot"), "name", "name",
                  43);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m"), "resolved",
                  "resolved", 43);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742668U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c12_rhs43, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs43, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m"), "context",
                  "context", 44);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 44);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 44);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1299109168U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c12_rhs44, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs44, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m!below_threshold"),
                  "context", "context", 45);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("length"), "name", "name", 45);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 45);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1303178606U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c12_rhs45, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs45, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m!intlength"),
                  "context", "context", 46);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 46);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 46);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c12_rhs46, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs46, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"),
                  "context", "context", 47);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 47);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 47);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c12_rhs47, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs47, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"),
                  "context", "context", 48);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_refblas_xdot"), "name",
                  "name", 48);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m"),
                  "resolved", "resolved", 48);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1299109172U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c12_rhs48, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs48, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m"),
                  "context", "context", 49);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_refblas_xdotx"), "name",
                  "name", 49);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1360314750U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c12_rhs49, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs49, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 50);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 50);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 50);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c12_rhs50, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs50, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 51);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 51);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 51);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c12_rhs51, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs51, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 52);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 52);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1346542740U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c12_rhs52, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs52, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 53);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 53);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 53);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 53);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c12_rhs53, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs53, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 54);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 54);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 54);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c12_rhs54, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs54, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 55);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 55);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c12_rhs55, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs55, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 56);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_xgemv"), "name", "name",
                  56);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemv.m"),
                  "resolved", "resolved", 56);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742670U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c12_rhs56, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs56, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemv.m"), "context",
                  "context", 57);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 57);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 57);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1299109168U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c12_rhs57, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs57, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m!below_threshold"),
                  "context", "context", 58);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("length"), "name", "name", 58);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 58);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1303178606U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c12_rhs58, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs58, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m!below_threshold"),
                  "context", "context", 59);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("intmax"), "name", "name", 59);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 59);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1362294282U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c12_rhs59, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs59, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m!below_threshold"),
                  "context", "context", 60);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("mtimes"), "name", "name", 60);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 60);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 60);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c12_rhs60, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs60, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m"),
                  "context", "context", 61);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 61);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 61);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c12_rhs61, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs61, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m"),
                  "context", "context", 62);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 62);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 62);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c12_rhs62, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs62, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemv.m"),
                  "context", "context", 63);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_refblas_xgemv"), "name",
                  "name", 63);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "resolved", "resolved", 63);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1360314752U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c12_rhs63, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs63, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c12_rhs0);
  sf_mex_destroy(&c12_lhs0);
  sf_mex_destroy(&c12_rhs1);
  sf_mex_destroy(&c12_lhs1);
  sf_mex_destroy(&c12_rhs2);
  sf_mex_destroy(&c12_lhs2);
  sf_mex_destroy(&c12_rhs3);
  sf_mex_destroy(&c12_lhs3);
  sf_mex_destroy(&c12_rhs4);
  sf_mex_destroy(&c12_lhs4);
  sf_mex_destroy(&c12_rhs5);
  sf_mex_destroy(&c12_lhs5);
  sf_mex_destroy(&c12_rhs6);
  sf_mex_destroy(&c12_lhs6);
  sf_mex_destroy(&c12_rhs7);
  sf_mex_destroy(&c12_lhs7);
  sf_mex_destroy(&c12_rhs8);
  sf_mex_destroy(&c12_lhs8);
  sf_mex_destroy(&c12_rhs9);
  sf_mex_destroy(&c12_lhs9);
  sf_mex_destroy(&c12_rhs10);
  sf_mex_destroy(&c12_lhs10);
  sf_mex_destroy(&c12_rhs11);
  sf_mex_destroy(&c12_lhs11);
  sf_mex_destroy(&c12_rhs12);
  sf_mex_destroy(&c12_lhs12);
  sf_mex_destroy(&c12_rhs13);
  sf_mex_destroy(&c12_lhs13);
  sf_mex_destroy(&c12_rhs14);
  sf_mex_destroy(&c12_lhs14);
  sf_mex_destroy(&c12_rhs15);
  sf_mex_destroy(&c12_lhs15);
  sf_mex_destroy(&c12_rhs16);
  sf_mex_destroy(&c12_lhs16);
  sf_mex_destroy(&c12_rhs17);
  sf_mex_destroy(&c12_lhs17);
  sf_mex_destroy(&c12_rhs18);
  sf_mex_destroy(&c12_lhs18);
  sf_mex_destroy(&c12_rhs19);
  sf_mex_destroy(&c12_lhs19);
  sf_mex_destroy(&c12_rhs20);
  sf_mex_destroy(&c12_lhs20);
  sf_mex_destroy(&c12_rhs21);
  sf_mex_destroy(&c12_lhs21);
  sf_mex_destroy(&c12_rhs22);
  sf_mex_destroy(&c12_lhs22);
  sf_mex_destroy(&c12_rhs23);
  sf_mex_destroy(&c12_lhs23);
  sf_mex_destroy(&c12_rhs24);
  sf_mex_destroy(&c12_lhs24);
  sf_mex_destroy(&c12_rhs25);
  sf_mex_destroy(&c12_lhs25);
  sf_mex_destroy(&c12_rhs26);
  sf_mex_destroy(&c12_lhs26);
  sf_mex_destroy(&c12_rhs27);
  sf_mex_destroy(&c12_lhs27);
  sf_mex_destroy(&c12_rhs28);
  sf_mex_destroy(&c12_lhs28);
  sf_mex_destroy(&c12_rhs29);
  sf_mex_destroy(&c12_lhs29);
  sf_mex_destroy(&c12_rhs30);
  sf_mex_destroy(&c12_lhs30);
  sf_mex_destroy(&c12_rhs31);
  sf_mex_destroy(&c12_lhs31);
  sf_mex_destroy(&c12_rhs32);
  sf_mex_destroy(&c12_lhs32);
  sf_mex_destroy(&c12_rhs33);
  sf_mex_destroy(&c12_lhs33);
  sf_mex_destroy(&c12_rhs34);
  sf_mex_destroy(&c12_lhs34);
  sf_mex_destroy(&c12_rhs35);
  sf_mex_destroy(&c12_lhs35);
  sf_mex_destroy(&c12_rhs36);
  sf_mex_destroy(&c12_lhs36);
  sf_mex_destroy(&c12_rhs37);
  sf_mex_destroy(&c12_lhs37);
  sf_mex_destroy(&c12_rhs38);
  sf_mex_destroy(&c12_lhs38);
  sf_mex_destroy(&c12_rhs39);
  sf_mex_destroy(&c12_lhs39);
  sf_mex_destroy(&c12_rhs40);
  sf_mex_destroy(&c12_lhs40);
  sf_mex_destroy(&c12_rhs41);
  sf_mex_destroy(&c12_lhs41);
  sf_mex_destroy(&c12_rhs42);
  sf_mex_destroy(&c12_lhs42);
  sf_mex_destroy(&c12_rhs43);
  sf_mex_destroy(&c12_lhs43);
  sf_mex_destroy(&c12_rhs44);
  sf_mex_destroy(&c12_lhs44);
  sf_mex_destroy(&c12_rhs45);
  sf_mex_destroy(&c12_lhs45);
  sf_mex_destroy(&c12_rhs46);
  sf_mex_destroy(&c12_lhs46);
  sf_mex_destroy(&c12_rhs47);
  sf_mex_destroy(&c12_lhs47);
  sf_mex_destroy(&c12_rhs48);
  sf_mex_destroy(&c12_lhs48);
  sf_mex_destroy(&c12_rhs49);
  sf_mex_destroy(&c12_lhs49);
  sf_mex_destroy(&c12_rhs50);
  sf_mex_destroy(&c12_lhs50);
  sf_mex_destroy(&c12_rhs51);
  sf_mex_destroy(&c12_lhs51);
  sf_mex_destroy(&c12_rhs52);
  sf_mex_destroy(&c12_lhs52);
  sf_mex_destroy(&c12_rhs53);
  sf_mex_destroy(&c12_lhs53);
  sf_mex_destroy(&c12_rhs54);
  sf_mex_destroy(&c12_lhs54);
  sf_mex_destroy(&c12_rhs55);
  sf_mex_destroy(&c12_lhs55);
  sf_mex_destroy(&c12_rhs56);
  sf_mex_destroy(&c12_lhs56);
  sf_mex_destroy(&c12_rhs57);
  sf_mex_destroy(&c12_lhs57);
  sf_mex_destroy(&c12_rhs58);
  sf_mex_destroy(&c12_lhs58);
  sf_mex_destroy(&c12_rhs59);
  sf_mex_destroy(&c12_lhs59);
  sf_mex_destroy(&c12_rhs60);
  sf_mex_destroy(&c12_lhs60);
  sf_mex_destroy(&c12_rhs61);
  sf_mex_destroy(&c12_lhs61);
  sf_mex_destroy(&c12_rhs62);
  sf_mex_destroy(&c12_lhs62);
  sf_mex_destroy(&c12_rhs63);
  sf_mex_destroy(&c12_lhs63);
}

static const mxArray *c12_emlrt_marshallOut(char * c12_u)
{
  const mxArray *c12_y = NULL;
  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", c12_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c12_u)), FALSE);
  return c12_y;
}

static const mxArray *c12_b_emlrt_marshallOut(uint32_T c12_u)
{
  const mxArray *c12_y = NULL;
  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", &c12_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c12_y;
}

static void c12_b_info_helper(const mxArray **c12_info)
{
  const mxArray *c12_rhs64 = NULL;
  const mxArray *c12_lhs64 = NULL;
  const mxArray *c12_rhs65 = NULL;
  const mxArray *c12_lhs65 = NULL;
  const mxArray *c12_rhs66 = NULL;
  const mxArray *c12_lhs66 = NULL;
  const mxArray *c12_rhs67 = NULL;
  const mxArray *c12_lhs67 = NULL;
  const mxArray *c12_rhs68 = NULL;
  const mxArray *c12_lhs68 = NULL;
  const mxArray *c12_rhs69 = NULL;
  const mxArray *c12_lhs69 = NULL;
  const mxArray *c12_rhs70 = NULL;
  const mxArray *c12_lhs70 = NULL;
  const mxArray *c12_rhs71 = NULL;
  const mxArray *c12_lhs71 = NULL;
  const mxArray *c12_rhs72 = NULL;
  const mxArray *c12_lhs72 = NULL;
  const mxArray *c12_rhs73 = NULL;
  const mxArray *c12_lhs73 = NULL;
  const mxArray *c12_rhs74 = NULL;
  const mxArray *c12_lhs74 = NULL;
  const mxArray *c12_rhs75 = NULL;
  const mxArray *c12_lhs75 = NULL;
  const mxArray *c12_rhs76 = NULL;
  const mxArray *c12_lhs76 = NULL;
  const mxArray *c12_rhs77 = NULL;
  const mxArray *c12_lhs77 = NULL;
  const mxArray *c12_rhs78 = NULL;
  const mxArray *c12_lhs78 = NULL;
  const mxArray *c12_rhs79 = NULL;
  const mxArray *c12_lhs79 = NULL;
  const mxArray *c12_rhs80 = NULL;
  const mxArray *c12_lhs80 = NULL;
  const mxArray *c12_rhs81 = NULL;
  const mxArray *c12_lhs81 = NULL;
  const mxArray *c12_rhs82 = NULL;
  const mxArray *c12_lhs82 = NULL;
  const mxArray *c12_rhs83 = NULL;
  const mxArray *c12_lhs83 = NULL;
  const mxArray *c12_rhs84 = NULL;
  const mxArray *c12_lhs84 = NULL;
  const mxArray *c12_rhs85 = NULL;
  const mxArray *c12_lhs85 = NULL;
  const mxArray *c12_rhs86 = NULL;
  const mxArray *c12_lhs86 = NULL;
  const mxArray *c12_rhs87 = NULL;
  const mxArray *c12_lhs87 = NULL;
  const mxArray *c12_rhs88 = NULL;
  const mxArray *c12_lhs88 = NULL;
  const mxArray *c12_rhs89 = NULL;
  const mxArray *c12_lhs89 = NULL;
  const mxArray *c12_rhs90 = NULL;
  const mxArray *c12_lhs90 = NULL;
  const mxArray *c12_rhs91 = NULL;
  const mxArray *c12_lhs91 = NULL;
  const mxArray *c12_rhs92 = NULL;
  const mxArray *c12_lhs92 = NULL;
  const mxArray *c12_rhs93 = NULL;
  const mxArray *c12_lhs93 = NULL;
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 64);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 64);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c12_rhs64, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs64, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 65);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 65);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 65);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c12_rhs65, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs65, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 66);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 66);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 66);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 66);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851180U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c12_rhs66, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs66, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 67);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 67);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 67);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 67);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c12_rhs67, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs67, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemv.m"),
                  "context", "context", 68);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 68);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 68);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 68);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1346542740U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c12_rhs68, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs68, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 69);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_div"), "name", "name",
                  69);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 69);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 69);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742666U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c12_rhs69, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs69, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 70);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_xscal"), "name", "name",
                  70);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"),
                  "resolved", "resolved", 70);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742672U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c12_rhs70, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs70, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"), "context",
                  "context", 71);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 71);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 71);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1299109168U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c12_rhs71, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs71, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xscal.m!below_threshold"),
                  "context", "context", 72);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("length"), "name", "name", 72);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 72);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 72);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1303178606U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c12_rhs72, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs72, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xscal.m"),
                  "context", "context", 73);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 73);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 73);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 73);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c12_rhs73, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs73, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xscal.m"),
                  "context", "context", 74);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 74);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 74);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 74);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c12_rhs74, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs74, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs74), "lhs", "lhs",
                  74);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xscal.m"),
                  "context", "context", 75);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_refblas_xscal"), "name",
                  "name", 75);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 75);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "resolved", "resolved", 75);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1299109184U), "fileTimeLo",
                  "fileTimeLo", 75);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 75);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 75);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 75);
  sf_mex_assign(&c12_rhs75, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs75, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs75), "rhs", "rhs",
                  75);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs75), "lhs", "lhs",
                  75);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 76);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 76);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 76);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 76);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 76);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 76);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 76);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 76);
  sf_mex_assign(&c12_rhs76, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs76, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs76), "rhs", "rhs",
                  76);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs76), "lhs", "lhs",
                  76);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 77);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 77);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 77);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 77);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 77);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 77);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 77);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 77);
  sf_mex_assign(&c12_rhs77, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs77, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs77), "rhs", "rhs",
                  77);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs77), "lhs", "lhs",
                  77);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 78);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 78);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 78);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 78);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851180U), "fileTimeLo",
                  "fileTimeLo", 78);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 78);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 78);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 78);
  sf_mex_assign(&c12_rhs78, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs78, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs78), "rhs", "rhs",
                  78);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs78), "lhs", "lhs",
                  78);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 79);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 79);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 79);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 79);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 79);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 79);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 79);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 79);
  sf_mex_assign(&c12_rhs79, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs79, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs79), "rhs", "rhs",
                  79);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs79), "lhs", "lhs",
                  79);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xscal.m"),
                  "context", "context", 80);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 80);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 80);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 80);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1346542740U), "fileTimeLo",
                  "fileTimeLo", 80);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 80);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 80);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 80);
  sf_mex_assign(&c12_rhs80, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs80, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs80), "rhs", "rhs",
                  80);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs80), "lhs", "lhs",
                  80);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 81);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 81);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 81);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 81);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851178U), "fileTimeLo",
                  "fileTimeLo", 81);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 81);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 81);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 81);
  sf_mex_assign(&c12_rhs81, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs81, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs81), "rhs", "rhs",
                  81);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs81), "lhs", "lhs",
                  81);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "context", "context", 82);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("sin"), "name", "name", 82);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 82);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 82);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1343862786U), "fileTimeLo",
                  "fileTimeLo", 82);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 82);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 82);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 82);
  sf_mex_assign(&c12_rhs82, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs82, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs82), "rhs", "rhs",
                  82);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs82), "lhs", "lhs",
                  82);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 83);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 83);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 83);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 83);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851136U), "fileTimeLo",
                  "fileTimeLo", 83);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 83);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 83);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 83);
  sf_mex_assign(&c12_rhs83, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs83, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs83), "rhs", "rhs",
                  83);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs83), "lhs", "lhs",
                  83);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "context", "context", 84);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("cos"), "name", "name", 84);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 84);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 84);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1343862772U), "fileTimeLo",
                  "fileTimeLo", 84);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 84);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 84);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 84);
  sf_mex_assign(&c12_rhs84, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs84, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs84), "rhs", "rhs",
                  84);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs84), "lhs", "lhs",
                  84);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 85);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 85);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 85);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 85);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851122U), "fileTimeLo",
                  "fileTimeLo", 85);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 85);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 85);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 85);
  sf_mex_assign(&c12_rhs85, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs85, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs85), "rhs", "rhs",
                  85);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs85), "lhs", "lhs",
                  85);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 86);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 86);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 86);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 86);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 86);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 86);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 86);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 86);
  sf_mex_assign(&c12_rhs86, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs86, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs86), "rhs", "rhs",
                  86);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs86), "lhs", "lhs",
                  86);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 87);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 87);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 87);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 87);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 87);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 87);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 87);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 87);
  sf_mex_assign(&c12_rhs87, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs87, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs87), "rhs", "rhs",
                  87);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs87), "lhs", "lhs",
                  87);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 88);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  88);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 88);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 88);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742670U), "fileTimeLo",
                  "fileTimeLo", 88);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 88);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 88);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 88);
  sf_mex_assign(&c12_rhs88, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs88, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs88), "rhs", "rhs",
                  88);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs88), "lhs", "lhs",
                  88);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 89);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 89);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 89);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 89);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1299109168U), "fileTimeLo",
                  "fileTimeLo", 89);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 89);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 89);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 89);
  sf_mex_assign(&c12_rhs89, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs89, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs89), "rhs", "rhs",
                  89);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs89), "lhs", "lhs",
                  89);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold"),
                  "context", "context", 90);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("mtimes"), "name", "name", 90);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 90);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 90);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1363742678U), "fileTimeLo",
                  "fileTimeLo", 90);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 90);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 90);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 90);
  sf_mex_assign(&c12_rhs90, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs90, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs90), "rhs", "rhs",
                  90);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs90), "lhs", "lhs",
                  90);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 91);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 91);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 91);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 91);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 91);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 91);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 91);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 91);
  sf_mex_assign(&c12_rhs91, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs91, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs91), "rhs", "rhs",
                  91);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs91), "lhs", "lhs",
                  91);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 92);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 92);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 92);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 92);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 92);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 92);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 92);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 92);
  sf_mex_assign(&c12_rhs92, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs92, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs92), "rhs", "rhs",
                  92);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs92), "lhs", "lhs",
                  92);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 93);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("eml_refblas_xgemm"), "name",
                  "name", 93);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 93);
  sf_mex_addfield(*c12_info, c12_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "resolved", "resolved", 93);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(1360314750U), "fileTimeLo",
                  "fileTimeLo", 93);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 93);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 93);
  sf_mex_addfield(*c12_info, c12_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 93);
  sf_mex_assign(&c12_rhs93, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c12_lhs93, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_rhs93), "rhs", "rhs",
                  93);
  sf_mex_addfield(*c12_info, sf_mex_duplicatearraysafe(&c12_lhs93), "lhs", "lhs",
                  93);
  sf_mex_destroy(&c12_rhs64);
  sf_mex_destroy(&c12_lhs64);
  sf_mex_destroy(&c12_rhs65);
  sf_mex_destroy(&c12_lhs65);
  sf_mex_destroy(&c12_rhs66);
  sf_mex_destroy(&c12_lhs66);
  sf_mex_destroy(&c12_rhs67);
  sf_mex_destroy(&c12_lhs67);
  sf_mex_destroy(&c12_rhs68);
  sf_mex_destroy(&c12_lhs68);
  sf_mex_destroy(&c12_rhs69);
  sf_mex_destroy(&c12_lhs69);
  sf_mex_destroy(&c12_rhs70);
  sf_mex_destroy(&c12_lhs70);
  sf_mex_destroy(&c12_rhs71);
  sf_mex_destroy(&c12_lhs71);
  sf_mex_destroy(&c12_rhs72);
  sf_mex_destroy(&c12_lhs72);
  sf_mex_destroy(&c12_rhs73);
  sf_mex_destroy(&c12_lhs73);
  sf_mex_destroy(&c12_rhs74);
  sf_mex_destroy(&c12_lhs74);
  sf_mex_destroy(&c12_rhs75);
  sf_mex_destroy(&c12_lhs75);
  sf_mex_destroy(&c12_rhs76);
  sf_mex_destroy(&c12_lhs76);
  sf_mex_destroy(&c12_rhs77);
  sf_mex_destroy(&c12_lhs77);
  sf_mex_destroy(&c12_rhs78);
  sf_mex_destroy(&c12_lhs78);
  sf_mex_destroy(&c12_rhs79);
  sf_mex_destroy(&c12_lhs79);
  sf_mex_destroy(&c12_rhs80);
  sf_mex_destroy(&c12_lhs80);
  sf_mex_destroy(&c12_rhs81);
  sf_mex_destroy(&c12_lhs81);
  sf_mex_destroy(&c12_rhs82);
  sf_mex_destroy(&c12_lhs82);
  sf_mex_destroy(&c12_rhs83);
  sf_mex_destroy(&c12_lhs83);
  sf_mex_destroy(&c12_rhs84);
  sf_mex_destroy(&c12_lhs84);
  sf_mex_destroy(&c12_rhs85);
  sf_mex_destroy(&c12_lhs85);
  sf_mex_destroy(&c12_rhs86);
  sf_mex_destroy(&c12_lhs86);
  sf_mex_destroy(&c12_rhs87);
  sf_mex_destroy(&c12_lhs87);
  sf_mex_destroy(&c12_rhs88);
  sf_mex_destroy(&c12_lhs88);
  sf_mex_destroy(&c12_rhs89);
  sf_mex_destroy(&c12_lhs89);
  sf_mex_destroy(&c12_rhs90);
  sf_mex_destroy(&c12_lhs90);
  sf_mex_destroy(&c12_rhs91);
  sf_mex_destroy(&c12_lhs91);
  sf_mex_destroy(&c12_rhs92);
  sf_mex_destroy(&c12_lhs92);
  sf_mex_destroy(&c12_rhs93);
  sf_mex_destroy(&c12_lhs93);
}

static real_T c12_mpower(SFc12_IMM_UKFInstanceStruct *chartInstance, real_T
  c12_a)
{
  real_T c12_b_a;
  real_T c12_c_a;
  real_T c12_ak;
  real_T c12_d_a;
  real_T c12_e_a;
  real_T c12_b;
  c12_b_a = c12_a;
  c12_c_a = c12_b_a;
  c12_eml_scalar_eg(chartInstance);
  c12_ak = c12_c_a;
  c12_d_a = c12_ak;
  c12_eml_scalar_eg(chartInstance);
  c12_e_a = c12_d_a;
  c12_b = c12_d_a;
  return c12_e_a * c12_b;
}

static void c12_eml_scalar_eg(SFc12_IMM_UKFInstanceStruct *chartInstance)
{
}

static void c12_eml_error(SFc12_IMM_UKFInstanceStruct *chartInstance)
{
  int32_T c12_i72;
  static char_T c12_cv0[48] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'c', 'h', 'o', 'l', '_', 'm', 'a', 't', 'r', 'i', 'x', 'M',
    'u', 's', 't', 'B', 'e', 'P', 'o', 's', 'D', 'e', 'f', 'W', 'i', 't', 'h',
    'R', 'e', 'a', 'l', 'D', 'i', 'a', 'g' };

  char_T c12_u[48];
  const mxArray *c12_y = NULL;
  for (c12_i72 = 0; c12_i72 < 48; c12_i72++) {
    c12_u[c12_i72] = c12_cv0[c12_i72];
  }

  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", c12_u, 10, 0U, 1U, 0U, 2, 1, 48),
                FALSE);
  sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U, 14,
    c12_y));
}

static void c12_b_eml_scalar_eg(SFc12_IMM_UKFInstanceStruct *chartInstance)
{
}

static void c12_eml_matlab_zpotrf(SFc12_IMM_UKFInstanceStruct *chartInstance,
  real_T c12_A[25], real_T c12_b_A[25], int32_T *c12_info)
{
  int32_T c12_i73;
  for (c12_i73 = 0; c12_i73 < 25; c12_i73++) {
    c12_b_A[c12_i73] = c12_A[c12_i73];
  }

  *c12_info = c12_b_eml_matlab_zpotrf(chartInstance, c12_b_A);
}

static void c12_check_forloop_overflow_error(SFc12_IMM_UKFInstanceStruct
  *chartInstance, boolean_T c12_overflow)
{
  int32_T c12_i74;
  static char_T c12_cv1[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c12_u[34];
  const mxArray *c12_y = NULL;
  int32_T c12_i75;
  static char_T c12_cv2[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c12_b_u[23];
  const mxArray *c12_b_y = NULL;
  if (!c12_overflow) {
  } else {
    for (c12_i74 = 0; c12_i74 < 34; c12_i74++) {
      c12_u[c12_i74] = c12_cv1[c12_i74];
    }

    c12_y = NULL;
    sf_mex_assign(&c12_y, sf_mex_create("y", c12_u, 10, 0U, 1U, 0U, 2, 1, 34),
                  FALSE);
    for (c12_i75 = 0; c12_i75 < 23; c12_i75++) {
      c12_b_u[c12_i75] = c12_cv2[c12_i75];
    }

    c12_b_y = NULL;
    sf_mex_assign(&c12_b_y, sf_mex_create("y", c12_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                  FALSE);
    sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U,
      14, c12_y, 14, c12_b_y));
  }
}

static void c12_eml_xgemv(SFc12_IMM_UKFInstanceStruct *chartInstance, int32_T
  c12_m, int32_T c12_n, int32_T c12_ia0, int32_T c12_ix0, real_T c12_y[25],
  int32_T c12_iy0, real_T c12_b_y[25])
{
  int32_T c12_i76;
  for (c12_i76 = 0; c12_i76 < 25; c12_i76++) {
    c12_b_y[c12_i76] = c12_y[c12_i76];
  }

  c12_b_eml_xgemv(chartInstance, c12_m, c12_n, c12_ia0, c12_ix0, c12_b_y,
                  c12_iy0);
}

static void c12_b_eml_error(SFc12_IMM_UKFInstanceStruct *chartInstance)
{
  int32_T c12_i77;
  static char_T c12_cv3[19] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'p', 'o', 's', 'd', 'e', 'f' };

  char_T c12_u[19];
  const mxArray *c12_y = NULL;
  for (c12_i77 = 0; c12_i77 < 19; c12_i77++) {
    c12_u[c12_i77] = c12_cv3[c12_i77];
  }

  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", c12_u, 10, 0U, 1U, 0U, 2, 1, 19),
                FALSE);
  sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U, 14,
    c12_y));
}

static void c12_c_eml_scalar_eg(SFc12_IMM_UKFInstanceStruct *chartInstance)
{
}

static const mxArray *c12_f_sf_marshallOut(void *chartInstanceVoid, void
  *c12_inData)
{
  const mxArray *c12_mxArrayOutData = NULL;
  int32_T c12_u;
  const mxArray *c12_y = NULL;
  SFc12_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc12_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c12_mxArrayOutData = NULL;
  c12_u = *(int32_T *)c12_inData;
  c12_y = NULL;
  sf_mex_assign(&c12_y, sf_mex_create("y", &c12_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c12_mxArrayOutData, c12_y, FALSE);
  return c12_mxArrayOutData;
}

static int32_T c12_g_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId)
{
  int32_T c12_y;
  int32_T c12_i78;
  sf_mex_import(c12_parentId, sf_mex_dup(c12_u), &c12_i78, 1, 6, 0U, 0, 0U, 0);
  c12_y = c12_i78;
  sf_mex_destroy(&c12_u);
  return c12_y;
}

static void c12_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c12_mxArrayInData, const char_T *c12_varName, void *c12_outData)
{
  const mxArray *c12_b_sfEvent;
  const char_T *c12_identifier;
  emlrtMsgIdentifier c12_thisId;
  int32_T c12_y;
  SFc12_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc12_IMM_UKFInstanceStruct *)chartInstanceVoid;
  c12_b_sfEvent = sf_mex_dup(c12_mxArrayInData);
  c12_identifier = c12_varName;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_y = c12_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c12_b_sfEvent),
    &c12_thisId);
  sf_mex_destroy(&c12_b_sfEvent);
  *(int32_T *)c12_outData = c12_y;
  sf_mex_destroy(&c12_mxArrayInData);
}

static uint8_T c12_h_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_b_is_active_c12_IMM_UKF, const char_T *c12_identifier)
{
  uint8_T c12_y;
  emlrtMsgIdentifier c12_thisId;
  c12_thisId.fIdentifier = c12_identifier;
  c12_thisId.fParent = NULL;
  c12_y = c12_i_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c12_b_is_active_c12_IMM_UKF), &c12_thisId);
  sf_mex_destroy(&c12_b_is_active_c12_IMM_UKF);
  return c12_y;
}

static uint8_T c12_i_emlrt_marshallIn(SFc12_IMM_UKFInstanceStruct *chartInstance,
  const mxArray *c12_u, const emlrtMsgIdentifier *c12_parentId)
{
  uint8_T c12_y;
  uint8_T c12_u0;
  sf_mex_import(c12_parentId, sf_mex_dup(c12_u), &c12_u0, 1, 3, 0U, 0, 0U, 0);
  c12_y = c12_u0;
  sf_mex_destroy(&c12_u);
  return c12_y;
}

static int32_T c12_b_eml_matlab_zpotrf(SFc12_IMM_UKFInstanceStruct
  *chartInstance, real_T c12_A[25])
{
  int32_T c12_info;
  int32_T c12_j;
  int32_T c12_b_j;
  int32_T c12_a;
  int32_T c12_jm1;
  int32_T c12_b_a;
  int32_T c12_c;
  int32_T c12_c_a;
  int32_T c12_b;
  int32_T c12_jj;
  int32_T c12_n;
  int32_T c12_ix0;
  int32_T c12_iy0;
  int32_T c12_b_n;
  int32_T c12_b_ix0;
  int32_T c12_b_iy0;
  int32_T c12_c_n;
  int32_T c12_c_ix0;
  int32_T c12_c_iy0;
  int32_T c12_d_n;
  int32_T c12_d_ix0;
  int32_T c12_d_iy0;
  int32_T c12_e_n;
  int32_T c12_e_ix0;
  int32_T c12_e_iy0;
  real_T c12_d;
  int32_T c12_ix;
  int32_T c12_iy;
  int32_T c12_f_n;
  int32_T c12_b_b;
  int32_T c12_c_b;
  boolean_T c12_overflow;
  int32_T c12_k;
  int32_T c12_d_a;
  int32_T c12_e_a;
  real_T c12_ajj;
  int32_T c12_d_b;
  int32_T c12_nmj;
  int32_T c12_f_a;
  int32_T c12_jp1;
  int32_T c12_g_a;
  int32_T c12_jp1j;
  int32_T c12_b_jm1;
  int32_T c12_e_b;
  int32_T c12_f_b;
  boolean_T c12_b_overflow;
  int32_T c12_b_k;
  int32_T c12_c_k;
  int32_T c12_c_jm1;
  int32_T c12_g_b;
  int32_T c12_h_b;
  boolean_T c12_c_overflow;
  int32_T c12_d_k;
  real_T c12_y;
  real_T c12_z;
  int32_T c12_g_n;
  real_T c12_h_a;
  int32_T c12_f_ix0;
  int32_T c12_h_n;
  real_T c12_i_a;
  int32_T c12_g_ix0;
  int32_T c12_i_n;
  real_T c12_j_a;
  int32_T c12_h_ix0;
  int32_T c12_i_ix0;
  int32_T c12_k_a;
  int32_T c12_b_c;
  int32_T c12_i_b;
  int32_T c12_c_c;
  int32_T c12_l_a;
  int32_T c12_j_b;
  int32_T c12_i79;
  int32_T c12_m_a;
  int32_T c12_k_b;
  int32_T c12_n_a;
  int32_T c12_l_b;
  boolean_T c12_d_overflow;
  int32_T c12_e_k;
  int32_T c12_f_k;
  boolean_T exitg1;
  c12_info = 0;
  c12_b_eml_scalar_eg(chartInstance);
  c12_j = 1;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (c12_j < 6)) {
    c12_b_j = c12_j;
    c12_a = c12_b_j - 1;
    c12_jm1 = c12_a;
    c12_b_a = c12_jm1;
    c12_c = c12_b_a * 5;
    c12_c_a = c12_b_j;
    c12_b = c12_c;
    c12_jj = c12_c_a + c12_b;
    c12_n = c12_jm1;
    c12_ix0 = c12_b_j;
    c12_iy0 = c12_b_j;
    c12_b_n = c12_n;
    c12_b_ix0 = c12_ix0;
    c12_b_iy0 = c12_iy0;
    c12_c_n = c12_b_n;
    c12_c_ix0 = c12_b_ix0;
    c12_c_iy0 = c12_b_iy0;
    c12_d_n = c12_c_n;
    c12_d_ix0 = c12_c_ix0;
    c12_d_iy0 = c12_c_iy0;
    c12_e_n = c12_d_n;
    c12_e_ix0 = c12_d_ix0;
    c12_e_iy0 = c12_d_iy0;
    c12_d = 0.0;
    if (c12_e_n < 1) {
    } else {
      c12_ix = c12_e_ix0;
      c12_iy = c12_e_iy0;
      c12_f_n = c12_e_n;
      c12_b_b = c12_f_n;
      c12_c_b = c12_b_b;
      if (1 > c12_c_b) {
        c12_overflow = FALSE;
      } else {
        c12_overflow = (c12_c_b > 2147483646);
      }

      if (c12_overflow) {
        c12_check_forloop_overflow_error(chartInstance, c12_overflow);
      }

      for (c12_k = 1; c12_k <= c12_f_n; c12_k++) {
        c12_d += c12_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c12_ix), 1, 25, 1, 0) - 1] *
          c12_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c12_iy), 1, 25, 1, 0) - 1];
        c12_d_a = c12_ix + 5;
        c12_ix = c12_d_a;
        c12_e_a = c12_iy + 5;
        c12_iy = c12_e_a;
      }
    }

    c12_ajj = c12_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c12_jj), 1, 25, 1, 0) - 1] - c12_d;
    if (c12_ajj > 0.0) {
      c12_ajj = muDoubleScalarSqrt(c12_ajj);
      c12_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c12_jj), 1, 25, 1, 0) - 1] = c12_ajj;
      if (c12_b_j < 5) {
        c12_d_b = c12_b_j;
        c12_nmj = 5 - c12_d_b;
        c12_f_a = c12_b_j;
        c12_jp1 = c12_f_a;
        c12_g_a = c12_jj + 1;
        c12_jp1j = c12_g_a;
        c12_b_jm1 = c12_jm1;
        c12_e_b = c12_b_jm1;
        c12_f_b = c12_e_b;
        if (1 > c12_f_b) {
          c12_b_overflow = FALSE;
        } else {
          c12_b_overflow = (c12_f_b > 2147483646);
        }

        if (c12_b_overflow) {
          c12_check_forloop_overflow_error(chartInstance, c12_b_overflow);
        }

        for (c12_b_k = 1; c12_b_k <= c12_b_jm1; c12_b_k++) {
          c12_c_k = c12_b_k;
          c12_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c12_b_j), 1, 5, 1, 0) + 5 *
                 (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c12_c_k), 1, 5, 2, 0) - 1)) - 1] = c12_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c12_b_j), 1, 5, 1, 0) + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c12_c_k), 1, 5, 2, 0)
               - 1)) - 1];
        }

        c12_b_eml_xgemv(chartInstance, c12_nmj, c12_jm1, c12_jp1 + 1, c12_b_j,
                        c12_A, c12_jp1j);
        c12_c_jm1 = c12_jm1;
        c12_g_b = c12_c_jm1;
        c12_h_b = c12_g_b;
        if (1 > c12_h_b) {
          c12_c_overflow = FALSE;
        } else {
          c12_c_overflow = (c12_h_b > 2147483646);
        }

        if (c12_c_overflow) {
          c12_check_forloop_overflow_error(chartInstance, c12_c_overflow);
        }

        for (c12_d_k = 1; c12_d_k <= c12_c_jm1; c12_d_k++) {
          c12_c_k = c12_d_k;
          c12_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c12_b_j), 1, 5, 1, 0) + 5 *
                 (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                    (real_T)c12_c_k), 1, 5, 2, 0) - 1)) - 1] = c12_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c12_b_j), 1, 5, 1, 0) + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK
               ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c12_c_k), 1, 5, 2, 0)
               - 1)) - 1];
        }

        c12_y = c12_ajj;
        c12_z = 1.0 / c12_y;
        c12_g_n = c12_nmj;
        c12_h_a = c12_z;
        c12_f_ix0 = c12_jp1j;
        c12_h_n = c12_g_n;
        c12_i_a = c12_h_a;
        c12_g_ix0 = c12_f_ix0;
        c12_i_n = c12_h_n;
        c12_j_a = c12_i_a;
        c12_h_ix0 = c12_g_ix0;
        c12_i_ix0 = c12_h_ix0;
        c12_k_a = c12_i_n;
        c12_b_c = c12_k_a;
        c12_i_b = c12_b_c - 1;
        c12_c_c = c12_i_b;
        c12_l_a = c12_h_ix0;
        c12_j_b = c12_c_c;
        c12_i79 = c12_l_a + c12_j_b;
        c12_m_a = c12_i_ix0;
        c12_k_b = c12_i79;
        c12_n_a = c12_m_a;
        c12_l_b = c12_k_b;
        if (c12_n_a > c12_l_b) {
          c12_d_overflow = FALSE;
        } else {
          c12_d_overflow = (c12_l_b > 2147483646);
        }

        if (c12_d_overflow) {
          c12_check_forloop_overflow_error(chartInstance, c12_d_overflow);
        }

        for (c12_e_k = c12_i_ix0; c12_e_k <= c12_i79; c12_e_k++) {
          c12_f_k = c12_e_k;
          c12_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c12_f_k), 1, 25, 1, 0) - 1] = c12_j_a *
            c12_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c12_f_k), 1, 25, 1, 0) - 1];
        }
      }

      c12_j++;
    } else {
      c12_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c12_jj), 1, 25, 1, 0) - 1] = c12_ajj;
      c12_info = c12_b_j;
      exitg1 = TRUE;
    }
  }

  return c12_info;
}

static void c12_b_eml_xgemv(SFc12_IMM_UKFInstanceStruct *chartInstance, int32_T
  c12_m, int32_T c12_n, int32_T c12_ia0, int32_T c12_ix0, real_T c12_y[25],
  int32_T c12_iy0)
{
  int32_T c12_b_m;
  int32_T c12_b_n;
  int32_T c12_b_ia0;
  int32_T c12_b_ix0;
  int32_T c12_b_iy0;
  int32_T c12_c_m;
  int32_T c12_c_n;
  real_T c12_alpha1;
  int32_T c12_c_ia0;
  int32_T c12_c_ix0;
  real_T c12_beta1;
  int32_T c12_c_iy0;
  char_T c12_TRANSA;
  int32_T c12_var;
  ptrdiff_t c12_m_t;
  int32_T c12_b_var;
  ptrdiff_t c12_n_t;
  ptrdiff_t c12_lda_t;
  ptrdiff_t c12_incx_t;
  ptrdiff_t c12_incy_t;
  double * c12_alpha1_t;
  double * c12_beta1_t;
  double * c12_yiy0_t;
  double * c12_yix0_t;
  double * c12_yia0_t;
  c12_b_m = c12_m;
  c12_b_n = c12_n;
  c12_b_ia0 = c12_ia0;
  c12_b_ix0 = c12_ix0;
  c12_b_iy0 = c12_iy0;
  if (c12_b_m < 1) {
  } else if (c12_b_n < 1) {
  } else {
    c12_c_m = c12_b_m;
    c12_c_n = c12_b_n;
    c12_alpha1 = -1.0;
    c12_c_ia0 = c12_b_ia0;
    c12_c_ix0 = c12_b_ix0;
    c12_beta1 = 1.0;
    c12_c_iy0 = c12_b_iy0;
    c12_TRANSA = 'N';
    c12_var = c12_c_m;
    c12_m_t = (ptrdiff_t)(c12_var);
    c12_b_var = c12_c_n;
    c12_n_t = (ptrdiff_t)(c12_b_var);
    c12_lda_t = (ptrdiff_t)(5);
    c12_incx_t = (ptrdiff_t)(5);
    c12_incy_t = (ptrdiff_t)(1);
    c12_alpha1_t = (double *)(&c12_alpha1);
    c12_beta1_t = (double *)(&c12_beta1);
    c12_yiy0_t = (double *)(&c12_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c12_c_iy0), 1, 25, 1, 0) - 1]);
    c12_yix0_t = (double *)(&c12_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c12_c_ix0), 1, 25, 1, 0) - 1]);
    c12_yia0_t = (double *)(&c12_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c12_c_ia0), 1, 25, 1, 0) - 1]);
    dgemv(&c12_TRANSA, &c12_m_t, &c12_n_t, c12_alpha1_t, c12_yia0_t, &c12_lda_t,
          c12_yix0_t, &c12_incx_t, c12_beta1_t, c12_yiy0_t, &c12_incy_t);
  }
}

static void init_dsm_address_info(SFc12_IMM_UKFInstanceStruct *chartInstance)
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

void sf_c12_IMM_UKF_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2458350672U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1339977128U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1228857506U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1732967944U);
}

mxArray *sf_c12_IMM_UKF_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("KW1AE7NhFgv790h5fIyQK");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
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
      pr[0] = (double)(5);
      pr[1] = (double)(5);
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
      pr[0] = (double)(1);
      pr[1] = (double)(1);
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
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c12_IMM_UKF_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c12_IMM_UKF_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c12_IMM_UKF(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[9],T\"x_m_p\",},{M[8],M[0],T\"is_active_c12_IMM_UKF\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c12_IMM_UKF_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc12_IMM_UKFInstanceStruct *chartInstance;
    chartInstance = (SFc12_IMM_UKFInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _IMM_UKFMachineNumber_,
           12,
           1,
           1,
           4,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"x");
          _SFD_SET_DATA_PROPS(1,1,1,0,"P");
          _SFD_SET_DATA_PROPS(2,1,1,0,"T1");
          _SFD_SET_DATA_PROPS(3,2,0,1,"x_m_p");
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
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,1,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,1335);
        _SFD_CV_INIT_EML_FOR(0,1,0,823,837,1298);

        {
          unsigned int dimVector[1];
          dimVector[0]= 5;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c12_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 5;
          dimVector[1]= 5;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c12_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c12_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 5;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c12_sf_marshallOut,(MexInFcnForType)
            c12_sf_marshallIn);
        }

        {
          real_T *c12_T1;
          real_T (*c12_x)[5];
          real_T (*c12_P)[25];
          real_T (*c12_x_m_p)[5];
          c12_x_m_p = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
          c12_T1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c12_P = (real_T (*)[25])ssGetInputPortSignal(chartInstance->S, 1);
          c12_x = (real_T (*)[5])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c12_x);
          _SFD_SET_DATA_VALUE_PTR(1U, *c12_P);
          _SFD_SET_DATA_VALUE_PTR(2U, c12_T1);
          _SFD_SET_DATA_VALUE_PTR(3U, *c12_x_m_p);
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
  return "P1t2z1hzLvOYcOWXlUGyb";
}

static void sf_opaque_initialize_c12_IMM_UKF(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc12_IMM_UKFInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c12_IMM_UKF((SFc12_IMM_UKFInstanceStruct*) chartInstanceVar);
  initialize_c12_IMM_UKF((SFc12_IMM_UKFInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c12_IMM_UKF(void *chartInstanceVar)
{
  enable_c12_IMM_UKF((SFc12_IMM_UKFInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c12_IMM_UKF(void *chartInstanceVar)
{
  disable_c12_IMM_UKF((SFc12_IMM_UKFInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c12_IMM_UKF(void *chartInstanceVar)
{
  sf_c12_IMM_UKF((SFc12_IMM_UKFInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c12_IMM_UKF(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c12_IMM_UKF((SFc12_IMM_UKFInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c12_IMM_UKF();/* state var info */
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

extern void sf_internal_set_sim_state_c12_IMM_UKF(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c12_IMM_UKF();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c12_IMM_UKF((SFc12_IMM_UKFInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c12_IMM_UKF(SimStruct* S)
{
  return sf_internal_get_sim_state_c12_IMM_UKF(S);
}

static void sf_opaque_set_sim_state_c12_IMM_UKF(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c12_IMM_UKF(S, st);
}

static void sf_opaque_terminate_c12_IMM_UKF(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc12_IMM_UKFInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_IMM_UKF_optimization_info();
    }

    finalize_c12_IMM_UKF((SFc12_IMM_UKFInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc12_IMM_UKF((SFc12_IMM_UKFInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c12_IMM_UKF(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c12_IMM_UKF((SFc12_IMM_UKFInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c12_IMM_UKF(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_IMM_UKF_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      12);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,12,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,12,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,12);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,12,3);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,12,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 3; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,12);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1791623954U));
  ssSetChecksum1(S,(2849634082U));
  ssSetChecksum2(S,(2026043968U));
  ssSetChecksum3(S,(3813902845U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c12_IMM_UKF(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c12_IMM_UKF(SimStruct *S)
{
  SFc12_IMM_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc12_IMM_UKFInstanceStruct *)utMalloc(sizeof
    (SFc12_IMM_UKFInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc12_IMM_UKFInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c12_IMM_UKF;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c12_IMM_UKF;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c12_IMM_UKF;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c12_IMM_UKF;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c12_IMM_UKF;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c12_IMM_UKF;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c12_IMM_UKF;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c12_IMM_UKF;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c12_IMM_UKF;
  chartInstance->chartInfo.mdlStart = mdlStart_c12_IMM_UKF;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c12_IMM_UKF;
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

void c12_IMM_UKF_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c12_IMM_UKF(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c12_IMM_UKF(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c12_IMM_UKF(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c12_IMM_UKF_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
