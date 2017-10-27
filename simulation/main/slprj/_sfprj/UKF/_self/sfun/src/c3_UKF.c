/* Include files */

#include <stddef.h>
#include "blas.h"
#include "UKF_sfun.h"
#include "c3_UKF.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "UKF_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c3_debug_family_names[13] = { "L", "sP", "chi_p", "chi_m_p",
  "i", "nargin", "nargout", "T1", "x", "P", "c", "Wm", "x_m_p" };

/* Function Declarations */
static void initialize_c3_UKF(SFc3_UKFInstanceStruct *chartInstance);
static void initialize_params_c3_UKF(SFc3_UKFInstanceStruct *chartInstance);
static void enable_c3_UKF(SFc3_UKFInstanceStruct *chartInstance);
static void disable_c3_UKF(SFc3_UKFInstanceStruct *chartInstance);
static void c3_update_debugger_state_c3_UKF(SFc3_UKFInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c3_UKF(SFc3_UKFInstanceStruct *chartInstance);
static void set_sim_state_c3_UKF(SFc3_UKFInstanceStruct *chartInstance, const
  mxArray *c3_st);
static void finalize_c3_UKF(SFc3_UKFInstanceStruct *chartInstance);
static void sf_gateway_c3_UKF(SFc3_UKFInstanceStruct *chartInstance);
static void c3_chartstep_c3_UKF(SFc3_UKFInstanceStruct *chartInstance);
static void initSimStructsc3_UKF(SFc3_UKFInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber, uint32_T c3_instanceNumber);
static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData);
static void c3_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance, const
  mxArray *c3_x_m_p, const char_T *c3_identifier, real_T c3_y[5]);
static void c3_b_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[5]);
static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static real_T c3_c_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_e_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_d_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[55]);
static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static void c3_e_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[25]);
static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static void c3_info_helper(const mxArray **c3_info);
static const mxArray *c3_emlrt_marshallOut(const char * c3_u);
static const mxArray *c3_b_emlrt_marshallOut(const uint32_T c3_u);
static void c3_b_info_helper(const mxArray **c3_info);
static int32_T c3_intmax(SFc3_UKFInstanceStruct *chartInstance);
static void c3_eml_error(SFc3_UKFInstanceStruct *chartInstance);
static void c3_eml_matlab_zpotrf(SFc3_UKFInstanceStruct *chartInstance, real_T
  c3_A[25], real_T c3_b_A[25], int32_T *c3_info);
static real_T c3_eml_xdotc(SFc3_UKFInstanceStruct *chartInstance, int32_T c3_n,
  real_T c3_x[25], int32_T c3_ix0, real_T c3_y[25], int32_T c3_iy0);
static void c3_check_forloop_overflow_error(SFc3_UKFInstanceStruct
  *chartInstance, boolean_T c3_overflow);
static void c3_eml_xgemv(SFc3_UKFInstanceStruct *chartInstance, int32_T c3_m,
  int32_T c3_n, int32_T c3_ia0, int32_T c3_ix0, real_T c3_y[25], int32_T c3_iy0,
  real_T c3_b_y[25]);
static void c3_below_threshold(SFc3_UKFInstanceStruct *chartInstance);
static void c3_b_eml_error(SFc3_UKFInstanceStruct *chartInstance);
static void c3_eml_scalar_eg(SFc3_UKFInstanceStruct *chartInstance);
static void c3_eml_xgemm(SFc3_UKFInstanceStruct *chartInstance, real_T c3_A[55],
  real_T c3_B[11], real_T c3_C[5], real_T c3_b_C[5]);
static const mxArray *c3_f_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static int32_T c3_f_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static uint8_T c3_g_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance,
  const mxArray *c3_b_is_active_c3_UKF, const char_T *c3_identifier);
static uint8_T c3_h_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static int32_T c3_b_eml_matlab_zpotrf(SFc3_UKFInstanceStruct *chartInstance,
  real_T c3_A[25]);
static void c3_b_eml_xgemv(SFc3_UKFInstanceStruct *chartInstance, int32_T c3_m,
  int32_T c3_n, int32_T c3_ia0, int32_T c3_ix0, real_T c3_y[25], int32_T c3_iy0);
static void c3_b_eml_xgemm(SFc3_UKFInstanceStruct *chartInstance, real_T c3_A[55],
  real_T c3_B[11], real_T c3_C[5]);
static void init_dsm_address_info(SFc3_UKFInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c3_UKF(SFc3_UKFInstanceStruct *chartInstance)
{
  chartInstance->c3_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c3_is_active_c3_UKF = 0U;
}

static void initialize_params_c3_UKF(SFc3_UKFInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void enable_c3_UKF(SFc3_UKFInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c3_UKF(SFc3_UKFInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c3_update_debugger_state_c3_UKF(SFc3_UKFInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c3_UKF(SFc3_UKFInstanceStruct *chartInstance)
{
  const mxArray *c3_st;
  const mxArray *c3_y = NULL;
  int32_T c3_i0;
  real_T c3_u[5];
  const mxArray *c3_b_y = NULL;
  uint8_T c3_hoistedGlobal;
  uint8_T c3_b_u;
  const mxArray *c3_c_y = NULL;
  real_T (*c3_x_m_p)[5];
  c3_x_m_p = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
  c3_st = NULL;
  c3_st = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_createcellmatrix(2, 1), false);
  for (c3_i0 = 0; c3_i0 < 5; c3_i0++) {
    c3_u[c3_i0] = (*c3_x_m_p)[c3_i0];
  }

  c3_b_y = NULL;
  sf_mex_assign(&c3_b_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 5), false);
  sf_mex_setcell(c3_y, 0, c3_b_y);
  c3_hoistedGlobal = chartInstance->c3_is_active_c3_UKF;
  c3_b_u = c3_hoistedGlobal;
  c3_c_y = NULL;
  sf_mex_assign(&c3_c_y, sf_mex_create("y", &c3_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c3_y, 1, c3_c_y);
  sf_mex_assign(&c3_st, c3_y, false);
  return c3_st;
}

static void set_sim_state_c3_UKF(SFc3_UKFInstanceStruct *chartInstance, const
  mxArray *c3_st)
{
  const mxArray *c3_u;
  real_T c3_dv0[5];
  int32_T c3_i1;
  real_T (*c3_x_m_p)[5];
  c3_x_m_p = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c3_doneDoubleBufferReInit = true;
  c3_u = sf_mex_dup(c3_st);
  c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 0)),
                      "x_m_p", c3_dv0);
  for (c3_i1 = 0; c3_i1 < 5; c3_i1++) {
    (*c3_x_m_p)[c3_i1] = c3_dv0[c3_i1];
  }

  chartInstance->c3_is_active_c3_UKF = c3_g_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c3_u, 1)), "is_active_c3_UKF");
  sf_mex_destroy(&c3_u);
  c3_update_debugger_state_c3_UKF(chartInstance);
  sf_mex_destroy(&c3_st);
}

static void finalize_c3_UKF(SFc3_UKFInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c3_UKF(SFc3_UKFInstanceStruct *chartInstance)
{
  int32_T c3_i2;
  int32_T c3_i3;
  int32_T c3_i4;
  int32_T c3_i5;
  real_T *c3_T1;
  real_T *c3_c;
  real_T (*c3_Wm)[11];
  real_T (*c3_P)[25];
  real_T (*c3_x_m_p)[5];
  real_T (*c3_x)[5];
  c3_Wm = (real_T (*)[11])ssGetInputPortSignal(chartInstance->S, 4);
  c3_c = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c3_P = (real_T (*)[25])ssGetInputPortSignal(chartInstance->S, 2);
  c3_x_m_p = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
  c3_x = (real_T (*)[5])ssGetInputPortSignal(chartInstance->S, 1);
  c3_T1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 1U, chartInstance->c3_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c3_T1, 0U);
  for (c3_i2 = 0; c3_i2 < 5; c3_i2++) {
    _SFD_DATA_RANGE_CHECK((*c3_x)[c3_i2], 1U);
  }

  chartInstance->c3_sfEvent = CALL_EVENT;
  c3_chartstep_c3_UKF(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_UKFMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c3_i3 = 0; c3_i3 < 5; c3_i3++) {
    _SFD_DATA_RANGE_CHECK((*c3_x_m_p)[c3_i3], 2U);
  }

  for (c3_i4 = 0; c3_i4 < 25; c3_i4++) {
    _SFD_DATA_RANGE_CHECK((*c3_P)[c3_i4], 3U);
  }

  _SFD_DATA_RANGE_CHECK(*c3_c, 4U);
  for (c3_i5 = 0; c3_i5 < 11; c3_i5++) {
    _SFD_DATA_RANGE_CHECK((*c3_Wm)[c3_i5], 5U);
  }
}

static void c3_chartstep_c3_UKF(SFc3_UKFInstanceStruct *chartInstance)
{
  real_T c3_hoistedGlobal;
  real_T c3_b_hoistedGlobal;
  real_T c3_T1;
  int32_T c3_i6;
  real_T c3_x[5];
  int32_T c3_i7;
  real_T c3_P[25];
  real_T c3_c;
  int32_T c3_i8;
  real_T c3_Wm[11];
  uint32_T c3_debug_family_var_map[13];
  real_T c3_L;
  real_T c3_sP[25];
  real_T c3_chi_p[55];
  real_T c3_chi_m_p[55];
  real_T c3_i;
  real_T c3_nargin = 5.0;
  real_T c3_nargout = 1.0;
  real_T c3_x_m_p[5];
  int32_T c3_i9;
  real_T c3_A[25];
  int32_T c3_j;
  int32_T c3_b_j;
  int32_T c3_info;
  int32_T c3_b_info;
  int32_T c3_c_info;
  int32_T c3_d_info;
  int32_T c3_jmax;
  int32_T c3_a;
  int32_T c3_b_a;
  int32_T c3_b_jmax;
  int32_T c3_b;
  int32_T c3_b_b;
  boolean_T c3_overflow;
  int32_T c3_c_j;
  int32_T c3_c_a;
  int32_T c3_d_a;
  int32_T c3_i10;
  int32_T c3_c_b;
  int32_T c3_d_b;
  boolean_T c3_b_overflow;
  int32_T c3_b_i;
  int32_T c3_c_i;
  int32_T c3_i11;
  int32_T c3_i12;
  real_T c3_e_a[5];
  int32_T c3_i13;
  int32_T c3_i14;
  int32_T c3_i15;
  real_T c3_f_a;
  int32_T c3_i16;
  real_T c3_e_b[25];
  int32_T c3_i17;
  int32_T c3_i18;
  int32_T c3_i19;
  int32_T c3_i20;
  int32_T c3_i21;
  real_T c3_y[25];
  real_T c3_g_a;
  int32_T c3_i22;
  real_T c3_f_b[25];
  int32_T c3_i23;
  int32_T c3_i24;
  int32_T c3_i25;
  int32_T c3_i26;
  int32_T c3_i27;
  int32_T c3_i28;
  int32_T c3_i29;
  int32_T c3_i30;
  int32_T c3_i31;
  int32_T c3_d_i;
  real_T c3_b_x;
  real_T c3_c_x;
  real_T c3_b_A;
  real_T c3_B;
  real_T c3_d_x;
  real_T c3_b_y;
  real_T c3_e_x;
  real_T c3_c_y;
  real_T c3_f_x;
  real_T c3_d_y;
  real_T c3_e_y;
  real_T c3_g_x;
  real_T c3_h_x;
  real_T c3_c_A;
  real_T c3_b_B;
  real_T c3_i_x;
  real_T c3_f_y;
  real_T c3_j_x;
  real_T c3_g_y;
  real_T c3_k_x;
  real_T c3_h_y;
  real_T c3_i_y;
  real_T c3_l_x;
  real_T c3_m_x;
  real_T c3_n_x;
  real_T c3_o_x;
  real_T c3_p_x;
  real_T c3_q_x;
  real_T c3_d_A;
  real_T c3_c_B;
  real_T c3_r_x;
  real_T c3_j_y;
  real_T c3_s_x;
  real_T c3_k_y;
  real_T c3_t_x;
  real_T c3_l_y;
  real_T c3_m_y;
  real_T c3_u_x;
  real_T c3_v_x;
  real_T c3_e_A;
  real_T c3_d_B;
  real_T c3_w_x;
  real_T c3_n_y;
  real_T c3_x_x;
  real_T c3_o_y;
  real_T c3_y_x;
  real_T c3_p_y;
  real_T c3_q_y;
  real_T c3_ab_x;
  real_T c3_bb_x;
  real_T c3_cb_x;
  real_T c3_db_x;
  int32_T c3_i32;
  real_T c3_h_a[55];
  int32_T c3_i33;
  real_T c3_g_b[11];
  int32_T c3_i34;
  int32_T c3_i35;
  int32_T c3_i36;
  real_T c3_dv1[55];
  int32_T c3_i37;
  real_T c3_dv2[11];
  int32_T c3_i38;
  real_T c3_dv3[55];
  int32_T c3_i39;
  real_T c3_dv4[11];
  int32_T c3_i40;
  real_T (*c3_b_x_m_p)[5];
  real_T *c3_b_T1;
  real_T *c3_b_c;
  real_T (*c3_b_Wm)[11];
  real_T (*c3_b_P)[25];
  real_T (*c3_eb_x)[5];
  c3_b_Wm = (real_T (*)[11])ssGetInputPortSignal(chartInstance->S, 4);
  c3_b_c = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c3_b_P = (real_T (*)[25])ssGetInputPortSignal(chartInstance->S, 2);
  c3_b_x_m_p = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
  c3_eb_x = (real_T (*)[5])ssGetInputPortSignal(chartInstance->S, 1);
  c3_b_T1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c3_sfEvent);
  c3_hoistedGlobal = *c3_b_T1;
  c3_b_hoistedGlobal = *c3_b_c;
  c3_T1 = c3_hoistedGlobal;
  for (c3_i6 = 0; c3_i6 < 5; c3_i6++) {
    c3_x[c3_i6] = (*c3_eb_x)[c3_i6];
  }

  for (c3_i7 = 0; c3_i7 < 25; c3_i7++) {
    c3_P[c3_i7] = (*c3_b_P)[c3_i7];
  }

  c3_c = c3_b_hoistedGlobal;
  for (c3_i8 = 0; c3_i8 < 11; c3_i8++) {
    c3_Wm[c3_i8] = (*c3_b_Wm)[c3_i8];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 13U, 13U, c3_debug_family_names,
    c3_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_L, 0U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_sP, 1U, c3_d_sf_marshallOut,
    c3_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_chi_p, 2U, c3_e_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_chi_m_p, 3U, c3_e_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i, 4U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 5U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 6U, c3_c_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_T1, 7U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_x, 8U, c3_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_P, 9U, c3_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_c, 10U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_Wm, 11U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_x_m_p, 12U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 2);
  c3_L = 5.0;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 3);
  for (c3_i9 = 0; c3_i9 < 25; c3_i9++) {
    c3_A[c3_i9] = c3_P[c3_i9];
  }

  for (c3_j = 1; c3_j < 6; c3_j++) {
    c3_b_j = c3_j;
  }

  c3_info = c3_b_eml_matlab_zpotrf(chartInstance, c3_A);
  c3_b_info = c3_info;
  c3_c_info = c3_b_info;
  c3_d_info = c3_c_info;
  if (c3_d_info == 0) {
    c3_jmax = 5;
  } else {
    c3_b_eml_error(chartInstance);
    c3_a = c3_d_info;
    c3_b_a = c3_a - 1;
    c3_jmax = c3_b_a;
  }

  c3_b_jmax = c3_jmax;
  c3_b = c3_b_jmax;
  c3_b_b = c3_b;
  if (2 > c3_b_b) {
    c3_overflow = false;
  } else {
    c3_overflow = (c3_b_b > c3_intmax(chartInstance) - 1);
  }

  if (c3_overflow) {
    c3_check_forloop_overflow_error(chartInstance, c3_overflow);
  }

  for (c3_c_j = 2; c3_c_j <= c3_b_jmax; c3_c_j++) {
    c3_b_j = c3_c_j;
    c3_c_a = c3_b_j;
    c3_d_a = c3_c_a - 1;
    c3_i10 = c3_d_a;
    c3_c_b = c3_i10;
    c3_d_b = c3_c_b;
    if (1 > c3_d_b) {
      c3_b_overflow = false;
    } else {
      c3_b_overflow = (c3_d_b > c3_intmax(chartInstance) - 1);
    }

    if (c3_b_overflow) {
      c3_check_forloop_overflow_error(chartInstance, c3_b_overflow);
    }

    for (c3_b_i = 1; c3_b_i <= c3_i10; c3_b_i++) {
      c3_c_i = c3_b_i;
      c3_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
              (real_T)c3_c_i), 1, 5, 1, 0) + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
              (int32_T)_SFD_INTEGER_CHECK("", (real_T)c3_b_j), 1, 5, 2, 0) - 1))
        - 1] = 0.0;
    }
  }

  for (c3_i11 = 0; c3_i11 < 25; c3_i11++) {
    c3_sP[c3_i11] = c3_A[c3_i11];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 4);
  for (c3_i12 = 0; c3_i12 < 5; c3_i12++) {
    c3_e_a[c3_i12] = c3_x[c3_i12];
  }

  for (c3_i13 = 0; c3_i13 < 5; c3_i13++) {
    c3_i14 = 0;
    for (c3_i15 = 0; c3_i15 < 5; c3_i15++) {
      c3_A[c3_i14 + c3_i13] = c3_e_a[c3_i13];
      c3_i14 += 5;
    }
  }

  c3_f_a = c3_c;
  for (c3_i16 = 0; c3_i16 < 25; c3_i16++) {
    c3_e_b[c3_i16] = c3_sP[c3_i16];
  }

  for (c3_i17 = 0; c3_i17 < 25; c3_i17++) {
    c3_e_b[c3_i17] *= c3_f_a;
  }

  for (c3_i18 = 0; c3_i18 < 5; c3_i18++) {
    c3_e_a[c3_i18] = c3_x[c3_i18];
  }

  for (c3_i19 = 0; c3_i19 < 5; c3_i19++) {
    c3_i20 = 0;
    for (c3_i21 = 0; c3_i21 < 5; c3_i21++) {
      c3_y[c3_i20 + c3_i19] = c3_e_a[c3_i19];
      c3_i20 += 5;
    }
  }

  c3_g_a = c3_c;
  for (c3_i22 = 0; c3_i22 < 25; c3_i22++) {
    c3_f_b[c3_i22] = c3_sP[c3_i22];
  }

  for (c3_i23 = 0; c3_i23 < 25; c3_i23++) {
    c3_f_b[c3_i23] *= c3_g_a;
  }

  for (c3_i24 = 0; c3_i24 < 5; c3_i24++) {
    c3_chi_p[c3_i24] = c3_x[c3_i24];
  }

  c3_i25 = 0;
  for (c3_i26 = 0; c3_i26 < 5; c3_i26++) {
    for (c3_i27 = 0; c3_i27 < 5; c3_i27++) {
      c3_chi_p[(c3_i27 + c3_i25) + 5] = c3_A[c3_i27 + c3_i25] + c3_e_b[c3_i27 +
        c3_i25];
    }

    c3_i25 += 5;
  }

  c3_i28 = 0;
  for (c3_i29 = 0; c3_i29 < 5; c3_i29++) {
    for (c3_i30 = 0; c3_i30 < 5; c3_i30++) {
      c3_chi_p[(c3_i30 + c3_i28) + 30] = c3_y[c3_i30 + c3_i28] - c3_f_b[c3_i30 +
        c3_i28];
    }

    c3_i28 += 5;
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 7);
  for (c3_i31 = 0; c3_i31 < 55; c3_i31++) {
    c3_chi_m_p[c3_i31] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 8);
  c3_i = 1.0;
  c3_d_i = 0;
  while (c3_d_i < 11) {
    c3_i = 1.0 + (real_T)c3_d_i;
    CV_EML_FOR(0, 1, 0, 1);
    _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 9);
    c3_b_x = c3_chi_p[4 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] * c3_T1;
    c3_c_x = c3_b_x;
    c3_c_x = muDoubleScalarSin(c3_c_x);
    c3_b_A = c3_chi_p[1 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] * c3_c_x;
    c3_B = c3_chi_p[4 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)];
    c3_d_x = c3_b_A;
    c3_b_y = c3_B;
    c3_e_x = c3_d_x;
    c3_c_y = c3_b_y;
    c3_f_x = c3_e_x;
    c3_d_y = c3_c_y;
    c3_e_y = c3_f_x / c3_d_y;
    c3_g_x = c3_chi_p[4 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] * c3_T1;
    c3_h_x = c3_g_x;
    c3_h_x = muDoubleScalarCos(c3_h_x);
    c3_c_A = c3_chi_p[3 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] * (1.0 - c3_h_x);
    c3_b_B = c3_chi_p[4 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)];
    c3_i_x = c3_c_A;
    c3_f_y = c3_b_B;
    c3_j_x = c3_i_x;
    c3_g_y = c3_f_y;
    c3_k_x = c3_j_x;
    c3_h_y = c3_g_y;
    c3_i_y = c3_k_x / c3_h_y;
    c3_chi_m_p[5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] = (c3_chi_p[5 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)_SFD_INTEGER_CHECK("i",
      c3_i), 1, 11, 2, 0) - 1)] + c3_e_y) - c3_i_y;
    _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 10);
    c3_l_x = c3_chi_p[4 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] * c3_T1;
    c3_m_x = c3_l_x;
    c3_m_x = muDoubleScalarCos(c3_m_x);
    c3_n_x = c3_chi_p[4 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] * c3_T1;
    c3_o_x = c3_n_x;
    c3_o_x = muDoubleScalarSin(c3_o_x);
    c3_chi_m_p[1 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] = c3_m_x * c3_chi_p[1 +
      5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)_SFD_INTEGER_CHECK("i",
             c3_i), 1, 11, 2, 0) - 1)] - c3_o_x * c3_chi_p[3 + 5 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)_SFD_INTEGER_CHECK("i",
         c3_i), 1, 11, 2, 0) - 1)];
    _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 11);
    c3_p_x = c3_chi_p[4 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] * c3_T1;
    c3_q_x = c3_p_x;
    c3_q_x = muDoubleScalarCos(c3_q_x);
    c3_d_A = c3_chi_p[1 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] * (1.0 - c3_q_x);
    c3_c_B = c3_chi_p[4 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)];
    c3_r_x = c3_d_A;
    c3_j_y = c3_c_B;
    c3_s_x = c3_r_x;
    c3_k_y = c3_j_y;
    c3_t_x = c3_s_x;
    c3_l_y = c3_k_y;
    c3_m_y = c3_t_x / c3_l_y;
    c3_u_x = c3_chi_p[4 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] * c3_T1;
    c3_v_x = c3_u_x;
    c3_v_x = muDoubleScalarSin(c3_v_x);
    c3_e_A = c3_chi_p[3 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] * c3_v_x;
    c3_d_B = c3_chi_p[4 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)];
    c3_w_x = c3_e_A;
    c3_n_y = c3_d_B;
    c3_x_x = c3_w_x;
    c3_o_y = c3_n_y;
    c3_y_x = c3_x_x;
    c3_p_y = c3_o_y;
    c3_q_y = c3_y_x / c3_p_y;
    c3_chi_m_p[2 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] = (c3_m_y + c3_chi_p[2 +
      5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)_SFD_INTEGER_CHECK("i",
      c3_i), 1, 11, 2, 0) - 1)]) + c3_q_y;
    _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 12);
    c3_ab_x = c3_chi_p[4 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] * c3_T1;
    c3_bb_x = c3_ab_x;
    c3_bb_x = muDoubleScalarSin(c3_bb_x);
    c3_cb_x = c3_chi_p[4 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] * c3_T1;
    c3_db_x = c3_cb_x;
    c3_db_x = muDoubleScalarCos(c3_db_x);
    c3_chi_m_p[3 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] = c3_chi_p[1 + 5 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)_SFD_INTEGER_CHECK("i",
         c3_i), 1, 11, 2, 0) - 1)] * c3_bb_x + c3_chi_p[3 + 5 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)_SFD_INTEGER_CHECK("i",
         c3_i), 1, 11, 2, 0) - 1)] * c3_db_x;
    _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 13);
    c3_chi_m_p[4 + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_m_p", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 11, 2, 0) - 1)] = c3_chi_p[4 + 5 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("chi_p", (int32_T)_SFD_INTEGER_CHECK("i",
         c3_i), 1, 11, 2, 0) - 1)];
    c3_d_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_EML_FOR(0, 1, 0, 0);
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 15);
  for (c3_i32 = 0; c3_i32 < 55; c3_i32++) {
    c3_h_a[c3_i32] = c3_chi_m_p[c3_i32];
  }

  for (c3_i33 = 0; c3_i33 < 11; c3_i33++) {
    c3_g_b[c3_i33] = c3_Wm[c3_i33];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i34 = 0; c3_i34 < 5; c3_i34++) {
    c3_x_m_p[c3_i34] = 0.0;
  }

  for (c3_i35 = 0; c3_i35 < 5; c3_i35++) {
    c3_x_m_p[c3_i35] = 0.0;
  }

  for (c3_i36 = 0; c3_i36 < 55; c3_i36++) {
    c3_dv1[c3_i36] = c3_h_a[c3_i36];
  }

  for (c3_i37 = 0; c3_i37 < 11; c3_i37++) {
    c3_dv2[c3_i37] = c3_g_b[c3_i37];
  }

  for (c3_i38 = 0; c3_i38 < 55; c3_i38++) {
    c3_dv3[c3_i38] = c3_dv1[c3_i38];
  }

  for (c3_i39 = 0; c3_i39 < 11; c3_i39++) {
    c3_dv4[c3_i39] = c3_dv2[c3_i39];
  }

  c3_b_eml_xgemm(chartInstance, c3_dv3, c3_dv4, c3_x_m_p);
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, -15);
  _SFD_SYMBOL_SCOPE_POP();
  for (c3_i40 = 0; c3_i40 < 5; c3_i40++) {
    (*c3_b_x_m_p)[c3_i40] = c3_x_m_p[c3_i40];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c3_sfEvent);
}

static void initSimStructsc3_UKF(SFc3_UKFInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber, uint32_T c3_instanceNumber)
{
  (void)c3_machineNumber;
  (void)c3_chartNumber;
  (void)c3_instanceNumber;
}

static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i41;
  real_T c3_b_inData[5];
  int32_T c3_i42;
  real_T c3_u[5];
  const mxArray *c3_y = NULL;
  SFc3_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc3_UKFInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i41 = 0; c3_i41 < 5; c3_i41++) {
    c3_b_inData[c3_i41] = (*(real_T (*)[5])c3_inData)[c3_i41];
  }

  for (c3_i42 = 0; c3_i42 < 5; c3_i42++) {
    c3_u[c3_i42] = c3_b_inData[c3_i42];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 5), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance, const
  mxArray *c3_x_m_p, const char_T *c3_identifier, real_T c3_y[5])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_x_m_p), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_x_m_p);
}

static void c3_b_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[5])
{
  real_T c3_dv5[5];
  int32_T c3_i43;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv5, 1, 0, 0U, 1, 0U, 1, 5);
  for (c3_i43 = 0; c3_i43 < 5; c3_i43++) {
    c3_y[c3_i43] = c3_dv5[c3_i43];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_x_m_p;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[5];
  int32_T c3_i44;
  SFc3_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc3_UKFInstanceStruct *)chartInstanceVoid;
  c3_x_m_p = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_x_m_p), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_x_m_p);
  for (c3_i44 = 0; c3_i44 < 5; c3_i44++) {
    (*(real_T (*)[5])c3_outData)[c3_i44] = c3_y[c3_i44];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i45;
  real_T c3_b_inData[11];
  int32_T c3_i46;
  real_T c3_u[11];
  const mxArray *c3_y = NULL;
  SFc3_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc3_UKFInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i45 = 0; c3_i45 < 11; c3_i45++) {
    c3_b_inData[c3_i45] = (*(real_T (*)[11])c3_inData)[c3_i45];
  }

  for (c3_i46 = 0; c3_i46 < 11; c3_i46++) {
    c3_u[c3_i46] = c3_b_inData[c3_i46];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 11), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  real_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc3_UKFInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(real_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i47;
  int32_T c3_i48;
  int32_T c3_i49;
  real_T c3_b_inData[25];
  int32_T c3_i50;
  int32_T c3_i51;
  int32_T c3_i52;
  real_T c3_u[25];
  const mxArray *c3_y = NULL;
  SFc3_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc3_UKFInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i47 = 0;
  for (c3_i48 = 0; c3_i48 < 5; c3_i48++) {
    for (c3_i49 = 0; c3_i49 < 5; c3_i49++) {
      c3_b_inData[c3_i49 + c3_i47] = (*(real_T (*)[25])c3_inData)[c3_i49 +
        c3_i47];
    }

    c3_i47 += 5;
  }

  c3_i50 = 0;
  for (c3_i51 = 0; c3_i51 < 5; c3_i51++) {
    for (c3_i52 = 0; c3_i52 < 5; c3_i52++) {
      c3_u[c3_i52 + c3_i50] = c3_b_inData[c3_i52 + c3_i50];
    }

    c3_i50 += 5;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 5, 5), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static real_T c3_c_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  real_T c3_y;
  real_T c3_d0;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_d0, 1, 0, 0U, 0, 0U, 0);
  c3_y = c3_d0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_nargout;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y;
  SFc3_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc3_UKFInstanceStruct *)chartInstanceVoid;
  c3_nargout = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_nargout), &c3_thisId);
  sf_mex_destroy(&c3_nargout);
  *(real_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_e_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i53;
  int32_T c3_i54;
  int32_T c3_i55;
  real_T c3_b_inData[55];
  int32_T c3_i56;
  int32_T c3_i57;
  int32_T c3_i58;
  real_T c3_u[55];
  const mxArray *c3_y = NULL;
  SFc3_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc3_UKFInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i53 = 0;
  for (c3_i54 = 0; c3_i54 < 11; c3_i54++) {
    for (c3_i55 = 0; c3_i55 < 5; c3_i55++) {
      c3_b_inData[c3_i55 + c3_i53] = (*(real_T (*)[55])c3_inData)[c3_i55 +
        c3_i53];
    }

    c3_i53 += 5;
  }

  c3_i56 = 0;
  for (c3_i57 = 0; c3_i57 < 11; c3_i57++) {
    for (c3_i58 = 0; c3_i58 < 5; c3_i58++) {
      c3_u[c3_i58 + c3_i56] = c3_b_inData[c3_i58 + c3_i56];
    }

    c3_i56 += 5;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 5, 11), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_d_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[55])
{
  real_T c3_dv6[55];
  int32_T c3_i59;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv6, 1, 0, 0U, 1, 0U, 2, 5, 11);
  for (c3_i59 = 0; c3_i59 < 55; c3_i59++) {
    c3_y[c3_i59] = c3_dv6[c3_i59];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_chi_m_p;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[55];
  int32_T c3_i60;
  int32_T c3_i61;
  int32_T c3_i62;
  SFc3_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc3_UKFInstanceStruct *)chartInstanceVoid;
  c3_chi_m_p = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_chi_m_p), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_chi_m_p);
  c3_i60 = 0;
  for (c3_i61 = 0; c3_i61 < 11; c3_i61++) {
    for (c3_i62 = 0; c3_i62 < 5; c3_i62++) {
      (*(real_T (*)[55])c3_outData)[c3_i62 + c3_i60] = c3_y[c3_i62 + c3_i60];
    }

    c3_i60 += 5;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static void c3_e_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[25])
{
  real_T c3_dv7[25];
  int32_T c3_i63;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv7, 1, 0, 0U, 1, 0U, 2, 5, 5);
  for (c3_i63 = 0; c3_i63 < 25; c3_i63++) {
    c3_y[c3_i63] = c3_dv7[c3_i63];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_sP;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[25];
  int32_T c3_i64;
  int32_T c3_i65;
  int32_T c3_i66;
  SFc3_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc3_UKFInstanceStruct *)chartInstanceVoid;
  c3_sP = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_sP), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_sP);
  c3_i64 = 0;
  for (c3_i65 = 0; c3_i65 < 5; c3_i65++) {
    for (c3_i66 = 0; c3_i66 < 5; c3_i66++) {
      (*(real_T (*)[25])c3_outData)[c3_i66 + c3_i64] = c3_y[c3_i66 + c3_i64];
    }

    c3_i64 += 5;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

const mxArray *sf_c3_UKF_get_eml_resolved_functions_info(void)
{
  const mxArray *c3_nameCaptureInfo = NULL;
  c3_nameCaptureInfo = NULL;
  sf_mex_assign(&c3_nameCaptureInfo, sf_mex_createstruct("structure", 2, 92, 1),
                false);
  c3_info_helper(&c3_nameCaptureInfo);
  c3_b_info_helper(&c3_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c3_nameCaptureInfo);
  return c3_nameCaptureInfo;
}

static void c3_info_helper(const mxArray **c3_info)
{
  const mxArray *c3_rhs0 = NULL;
  const mxArray *c3_lhs0 = NULL;
  const mxArray *c3_rhs1 = NULL;
  const mxArray *c3_lhs1 = NULL;
  const mxArray *c3_rhs2 = NULL;
  const mxArray *c3_lhs2 = NULL;
  const mxArray *c3_rhs3 = NULL;
  const mxArray *c3_lhs3 = NULL;
  const mxArray *c3_rhs4 = NULL;
  const mxArray *c3_lhs4 = NULL;
  const mxArray *c3_rhs5 = NULL;
  const mxArray *c3_lhs5 = NULL;
  const mxArray *c3_rhs6 = NULL;
  const mxArray *c3_lhs6 = NULL;
  const mxArray *c3_rhs7 = NULL;
  const mxArray *c3_lhs7 = NULL;
  const mxArray *c3_rhs8 = NULL;
  const mxArray *c3_lhs8 = NULL;
  const mxArray *c3_rhs9 = NULL;
  const mxArray *c3_lhs9 = NULL;
  const mxArray *c3_rhs10 = NULL;
  const mxArray *c3_lhs10 = NULL;
  const mxArray *c3_rhs11 = NULL;
  const mxArray *c3_lhs11 = NULL;
  const mxArray *c3_rhs12 = NULL;
  const mxArray *c3_lhs12 = NULL;
  const mxArray *c3_rhs13 = NULL;
  const mxArray *c3_lhs13 = NULL;
  const mxArray *c3_rhs14 = NULL;
  const mxArray *c3_lhs14 = NULL;
  const mxArray *c3_rhs15 = NULL;
  const mxArray *c3_lhs15 = NULL;
  const mxArray *c3_rhs16 = NULL;
  const mxArray *c3_lhs16 = NULL;
  const mxArray *c3_rhs17 = NULL;
  const mxArray *c3_lhs17 = NULL;
  const mxArray *c3_rhs18 = NULL;
  const mxArray *c3_lhs18 = NULL;
  const mxArray *c3_rhs19 = NULL;
  const mxArray *c3_lhs19 = NULL;
  const mxArray *c3_rhs20 = NULL;
  const mxArray *c3_lhs20 = NULL;
  const mxArray *c3_rhs21 = NULL;
  const mxArray *c3_lhs21 = NULL;
  const mxArray *c3_rhs22 = NULL;
  const mxArray *c3_lhs22 = NULL;
  const mxArray *c3_rhs23 = NULL;
  const mxArray *c3_lhs23 = NULL;
  const mxArray *c3_rhs24 = NULL;
  const mxArray *c3_lhs24 = NULL;
  const mxArray *c3_rhs25 = NULL;
  const mxArray *c3_lhs25 = NULL;
  const mxArray *c3_rhs26 = NULL;
  const mxArray *c3_lhs26 = NULL;
  const mxArray *c3_rhs27 = NULL;
  const mxArray *c3_lhs27 = NULL;
  const mxArray *c3_rhs28 = NULL;
  const mxArray *c3_lhs28 = NULL;
  const mxArray *c3_rhs29 = NULL;
  const mxArray *c3_lhs29 = NULL;
  const mxArray *c3_rhs30 = NULL;
  const mxArray *c3_lhs30 = NULL;
  const mxArray *c3_rhs31 = NULL;
  const mxArray *c3_lhs31 = NULL;
  const mxArray *c3_rhs32 = NULL;
  const mxArray *c3_lhs32 = NULL;
  const mxArray *c3_rhs33 = NULL;
  const mxArray *c3_lhs33 = NULL;
  const mxArray *c3_rhs34 = NULL;
  const mxArray *c3_lhs34 = NULL;
  const mxArray *c3_rhs35 = NULL;
  const mxArray *c3_lhs35 = NULL;
  const mxArray *c3_rhs36 = NULL;
  const mxArray *c3_lhs36 = NULL;
  const mxArray *c3_rhs37 = NULL;
  const mxArray *c3_lhs37 = NULL;
  const mxArray *c3_rhs38 = NULL;
  const mxArray *c3_lhs38 = NULL;
  const mxArray *c3_rhs39 = NULL;
  const mxArray *c3_lhs39 = NULL;
  const mxArray *c3_rhs40 = NULL;
  const mxArray *c3_lhs40 = NULL;
  const mxArray *c3_rhs41 = NULL;
  const mxArray *c3_lhs41 = NULL;
  const mxArray *c3_rhs42 = NULL;
  const mxArray *c3_lhs42 = NULL;
  const mxArray *c3_rhs43 = NULL;
  const mxArray *c3_lhs43 = NULL;
  const mxArray *c3_rhs44 = NULL;
  const mxArray *c3_lhs44 = NULL;
  const mxArray *c3_rhs45 = NULL;
  const mxArray *c3_lhs45 = NULL;
  const mxArray *c3_rhs46 = NULL;
  const mxArray *c3_lhs46 = NULL;
  const mxArray *c3_rhs47 = NULL;
  const mxArray *c3_lhs47 = NULL;
  const mxArray *c3_rhs48 = NULL;
  const mxArray *c3_lhs48 = NULL;
  const mxArray *c3_rhs49 = NULL;
  const mxArray *c3_lhs49 = NULL;
  const mxArray *c3_rhs50 = NULL;
  const mxArray *c3_lhs50 = NULL;
  const mxArray *c3_rhs51 = NULL;
  const mxArray *c3_lhs51 = NULL;
  const mxArray *c3_rhs52 = NULL;
  const mxArray *c3_lhs52 = NULL;
  const mxArray *c3_rhs53 = NULL;
  const mxArray *c3_lhs53 = NULL;
  const mxArray *c3_rhs54 = NULL;
  const mxArray *c3_lhs54 = NULL;
  const mxArray *c3_rhs55 = NULL;
  const mxArray *c3_lhs55 = NULL;
  const mxArray *c3_rhs56 = NULL;
  const mxArray *c3_lhs56 = NULL;
  const mxArray *c3_rhs57 = NULL;
  const mxArray *c3_lhs57 = NULL;
  const mxArray *c3_rhs58 = NULL;
  const mxArray *c3_lhs58 = NULL;
  const mxArray *c3_rhs59 = NULL;
  const mxArray *c3_lhs59 = NULL;
  const mxArray *c3_rhs60 = NULL;
  const mxArray *c3_lhs60 = NULL;
  const mxArray *c3_rhs61 = NULL;
  const mxArray *c3_lhs61 = NULL;
  const mxArray *c3_rhs62 = NULL;
  const mxArray *c3_lhs62 = NULL;
  const mxArray *c3_rhs63 = NULL;
  const mxArray *c3_lhs63 = NULL;
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("chol"), "name", "name", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1344504434U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c3_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_tolower"), "name", "name",
                  1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_tolower.m"), "resolved",
                  "resolved", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363742670U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c3_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c3_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("ismatrix"), "name", "name", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1331337258U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c3_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376013088U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c3_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmax"), "name", "name", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362294282U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c3_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1381882700U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c3_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_error"), "name", "name", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1343862758U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c3_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xpotrf"), "name", "name",
                  8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xpotrf.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286851208U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c3_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xpotrf.m"),
                  "context", "context", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_lapack_xpotrf"), "name",
                  "name", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xpotrf.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286851212U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c3_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xpotrf.m"),
                  "context", "context", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_matlab_zpotrf"), "name",
                  "name", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286851224U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c3_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c3_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376013088U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c3_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340320U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c3_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376013088U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c3_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372614816U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c3_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372615560U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c3_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372614816U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c3_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372615560U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c3_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372614816U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c3_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372615560U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c3_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xdotc"), "name", "name",
                  21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376013090U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c3_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c3_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"), "context",
                  "context", 23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xdotc"),
                  "name", "name", 23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdotc.p"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c3_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdotc.p"),
                  "context", "context", 24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xdot"),
                  "name", "name", 24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c3_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p"),
                  "context", "context", 25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c3_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p!below_threshold"),
                  "context", "context", 26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c3_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1381882700U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c3_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p!below_threshold"),
                  "context", "context", 28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("length"), "name", "name", 28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1303178606U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c3_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m!intlength"),
                  "context", "context", 29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c3_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p"),
                  "context", "context", 30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xdot"),
                  "name", "name", 30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdot.p"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c3_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdot.p"),
                  "context", "context", 31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xdotx"),
                  "name", "name", 31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdotx.p"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c3_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdotx.p"),
                  "context", "context", 32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340320U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c3_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdotx.p"),
                  "context", "context", 33);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 33);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376013088U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c3_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdotx.p"),
                  "context", "context", 34);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 34);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 34);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372615560U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c3_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 35);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 35);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 35);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372614816U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c3_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 36);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 36);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 36);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372615560U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c3_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 37);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 37);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372614816U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c3_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 38);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 38);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372615560U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c3_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 39);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xgemv"), "name", "name",
                  39);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemv.m"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376013090U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c3_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemv.m"), "context",
                  "context", 40);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 40);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c3_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemv.m"), "context",
                  "context", 41);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xgemv"),
                  "name", "name", 41);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c3_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p"),
                  "context", "context", 42);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 42);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c3_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p!below_threshold"),
                  "context", "context", 43);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 43);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c3_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p!below_threshold"),
                  "context", "context", 44);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("length"), "name", "name", 44);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 44);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1303178606U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c3_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p!below_threshold"),
                  "context", "context", 45);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.int"),
                  "name", "name", 45);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/int.p"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c3_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p!below_threshold"),
                  "context", "context", 46);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmax"), "name", "name", 46);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 46);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362294282U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c3_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p"),
                  "context", "context", 47);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 47);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 47);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340320U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c3_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p"),
                  "context", "context", 48);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xgemv"),
                  "name", "name", 48);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemv.p"),
                  "resolved", "resolved", 48);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c3_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemv.p"),
                  "context", "context", 49);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 49);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372615560U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c3_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemv.p"),
                  "context", "context", 50);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 50);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 50);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372615560U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c3_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemv.p"),
                  "context", "context", 51);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 51);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 51);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 51);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372615560U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c3_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemv.p"),
                  "context", "context", 52);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 52);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376013088U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c3_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p"),
                  "context", "context", 53);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.int"),
                  "name", "name", 53);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/int.p"),
                  "resolved", "resolved", 53);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c3_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p!ceval_xgemv"),
                  "context", "context", 54);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.size_ptr"),
                  "name", "name", 54);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/size_ptr.p"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c3_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p!ceval_xgemv"),
                  "context", "context", 55);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.c_cast"),
                  "name", "name", 55);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("int32"), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/c_cast.p"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340320U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c3_rhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 56);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_div"), "name", "name", 56);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 56);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376013088U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c3_rhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 57);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 57);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 57);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340320U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c3_rhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 58);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xscal"), "name", "name",
                  58);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"),
                  "resolved", "resolved", 58);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376013092U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c3_rhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"), "context",
                  "context", 59);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 59);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 59);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c3_rhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"), "context",
                  "context", 60);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xscal"),
                  "name", "name", 60);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 60);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p"),
                  "resolved", "resolved", 60);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c3_rhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p"),
                  "context", "context", 61);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 61);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 61);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c3_rhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p!below_threshold"),
                  "context", "context", 62);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 62);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 62);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c3_rhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p!below_threshold"),
                  "context", "context", 63);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("length"), "name", "name", 63);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 63);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1303178606U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c3_rhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c3_rhs0);
  sf_mex_destroy(&c3_lhs0);
  sf_mex_destroy(&c3_rhs1);
  sf_mex_destroy(&c3_lhs1);
  sf_mex_destroy(&c3_rhs2);
  sf_mex_destroy(&c3_lhs2);
  sf_mex_destroy(&c3_rhs3);
  sf_mex_destroy(&c3_lhs3);
  sf_mex_destroy(&c3_rhs4);
  sf_mex_destroy(&c3_lhs4);
  sf_mex_destroy(&c3_rhs5);
  sf_mex_destroy(&c3_lhs5);
  sf_mex_destroy(&c3_rhs6);
  sf_mex_destroy(&c3_lhs6);
  sf_mex_destroy(&c3_rhs7);
  sf_mex_destroy(&c3_lhs7);
  sf_mex_destroy(&c3_rhs8);
  sf_mex_destroy(&c3_lhs8);
  sf_mex_destroy(&c3_rhs9);
  sf_mex_destroy(&c3_lhs9);
  sf_mex_destroy(&c3_rhs10);
  sf_mex_destroy(&c3_lhs10);
  sf_mex_destroy(&c3_rhs11);
  sf_mex_destroy(&c3_lhs11);
  sf_mex_destroy(&c3_rhs12);
  sf_mex_destroy(&c3_lhs12);
  sf_mex_destroy(&c3_rhs13);
  sf_mex_destroy(&c3_lhs13);
  sf_mex_destroy(&c3_rhs14);
  sf_mex_destroy(&c3_lhs14);
  sf_mex_destroy(&c3_rhs15);
  sf_mex_destroy(&c3_lhs15);
  sf_mex_destroy(&c3_rhs16);
  sf_mex_destroy(&c3_lhs16);
  sf_mex_destroy(&c3_rhs17);
  sf_mex_destroy(&c3_lhs17);
  sf_mex_destroy(&c3_rhs18);
  sf_mex_destroy(&c3_lhs18);
  sf_mex_destroy(&c3_rhs19);
  sf_mex_destroy(&c3_lhs19);
  sf_mex_destroy(&c3_rhs20);
  sf_mex_destroy(&c3_lhs20);
  sf_mex_destroy(&c3_rhs21);
  sf_mex_destroy(&c3_lhs21);
  sf_mex_destroy(&c3_rhs22);
  sf_mex_destroy(&c3_lhs22);
  sf_mex_destroy(&c3_rhs23);
  sf_mex_destroy(&c3_lhs23);
  sf_mex_destroy(&c3_rhs24);
  sf_mex_destroy(&c3_lhs24);
  sf_mex_destroy(&c3_rhs25);
  sf_mex_destroy(&c3_lhs25);
  sf_mex_destroy(&c3_rhs26);
  sf_mex_destroy(&c3_lhs26);
  sf_mex_destroy(&c3_rhs27);
  sf_mex_destroy(&c3_lhs27);
  sf_mex_destroy(&c3_rhs28);
  sf_mex_destroy(&c3_lhs28);
  sf_mex_destroy(&c3_rhs29);
  sf_mex_destroy(&c3_lhs29);
  sf_mex_destroy(&c3_rhs30);
  sf_mex_destroy(&c3_lhs30);
  sf_mex_destroy(&c3_rhs31);
  sf_mex_destroy(&c3_lhs31);
  sf_mex_destroy(&c3_rhs32);
  sf_mex_destroy(&c3_lhs32);
  sf_mex_destroy(&c3_rhs33);
  sf_mex_destroy(&c3_lhs33);
  sf_mex_destroy(&c3_rhs34);
  sf_mex_destroy(&c3_lhs34);
  sf_mex_destroy(&c3_rhs35);
  sf_mex_destroy(&c3_lhs35);
  sf_mex_destroy(&c3_rhs36);
  sf_mex_destroy(&c3_lhs36);
  sf_mex_destroy(&c3_rhs37);
  sf_mex_destroy(&c3_lhs37);
  sf_mex_destroy(&c3_rhs38);
  sf_mex_destroy(&c3_lhs38);
  sf_mex_destroy(&c3_rhs39);
  sf_mex_destroy(&c3_lhs39);
  sf_mex_destroy(&c3_rhs40);
  sf_mex_destroy(&c3_lhs40);
  sf_mex_destroy(&c3_rhs41);
  sf_mex_destroy(&c3_lhs41);
  sf_mex_destroy(&c3_rhs42);
  sf_mex_destroy(&c3_lhs42);
  sf_mex_destroy(&c3_rhs43);
  sf_mex_destroy(&c3_lhs43);
  sf_mex_destroy(&c3_rhs44);
  sf_mex_destroy(&c3_lhs44);
  sf_mex_destroy(&c3_rhs45);
  sf_mex_destroy(&c3_lhs45);
  sf_mex_destroy(&c3_rhs46);
  sf_mex_destroy(&c3_lhs46);
  sf_mex_destroy(&c3_rhs47);
  sf_mex_destroy(&c3_lhs47);
  sf_mex_destroy(&c3_rhs48);
  sf_mex_destroy(&c3_lhs48);
  sf_mex_destroy(&c3_rhs49);
  sf_mex_destroy(&c3_lhs49);
  sf_mex_destroy(&c3_rhs50);
  sf_mex_destroy(&c3_lhs50);
  sf_mex_destroy(&c3_rhs51);
  sf_mex_destroy(&c3_lhs51);
  sf_mex_destroy(&c3_rhs52);
  sf_mex_destroy(&c3_lhs52);
  sf_mex_destroy(&c3_rhs53);
  sf_mex_destroy(&c3_lhs53);
  sf_mex_destroy(&c3_rhs54);
  sf_mex_destroy(&c3_lhs54);
  sf_mex_destroy(&c3_rhs55);
  sf_mex_destroy(&c3_lhs55);
  sf_mex_destroy(&c3_rhs56);
  sf_mex_destroy(&c3_lhs56);
  sf_mex_destroy(&c3_rhs57);
  sf_mex_destroy(&c3_lhs57);
  sf_mex_destroy(&c3_rhs58);
  sf_mex_destroy(&c3_lhs58);
  sf_mex_destroy(&c3_rhs59);
  sf_mex_destroy(&c3_lhs59);
  sf_mex_destroy(&c3_rhs60);
  sf_mex_destroy(&c3_lhs60);
  sf_mex_destroy(&c3_rhs61);
  sf_mex_destroy(&c3_lhs61);
  sf_mex_destroy(&c3_rhs62);
  sf_mex_destroy(&c3_lhs62);
  sf_mex_destroy(&c3_rhs63);
  sf_mex_destroy(&c3_lhs63);
}

static const mxArray *c3_emlrt_marshallOut(const char * c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c3_u)), false);
  return c3_y;
}

static const mxArray *c3_b_emlrt_marshallOut(const uint32_T c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 7, 0U, 0U, 0U, 0), false);
  return c3_y;
}

static void c3_b_info_helper(const mxArray **c3_info)
{
  const mxArray *c3_rhs64 = NULL;
  const mxArray *c3_lhs64 = NULL;
  const mxArray *c3_rhs65 = NULL;
  const mxArray *c3_lhs65 = NULL;
  const mxArray *c3_rhs66 = NULL;
  const mxArray *c3_lhs66 = NULL;
  const mxArray *c3_rhs67 = NULL;
  const mxArray *c3_lhs67 = NULL;
  const mxArray *c3_rhs68 = NULL;
  const mxArray *c3_lhs68 = NULL;
  const mxArray *c3_rhs69 = NULL;
  const mxArray *c3_lhs69 = NULL;
  const mxArray *c3_rhs70 = NULL;
  const mxArray *c3_lhs70 = NULL;
  const mxArray *c3_rhs71 = NULL;
  const mxArray *c3_lhs71 = NULL;
  const mxArray *c3_rhs72 = NULL;
  const mxArray *c3_lhs72 = NULL;
  const mxArray *c3_rhs73 = NULL;
  const mxArray *c3_lhs73 = NULL;
  const mxArray *c3_rhs74 = NULL;
  const mxArray *c3_lhs74 = NULL;
  const mxArray *c3_rhs75 = NULL;
  const mxArray *c3_lhs75 = NULL;
  const mxArray *c3_rhs76 = NULL;
  const mxArray *c3_lhs76 = NULL;
  const mxArray *c3_rhs77 = NULL;
  const mxArray *c3_lhs77 = NULL;
  const mxArray *c3_rhs78 = NULL;
  const mxArray *c3_lhs78 = NULL;
  const mxArray *c3_rhs79 = NULL;
  const mxArray *c3_lhs79 = NULL;
  const mxArray *c3_rhs80 = NULL;
  const mxArray *c3_lhs80 = NULL;
  const mxArray *c3_rhs81 = NULL;
  const mxArray *c3_lhs81 = NULL;
  const mxArray *c3_rhs82 = NULL;
  const mxArray *c3_lhs82 = NULL;
  const mxArray *c3_rhs83 = NULL;
  const mxArray *c3_lhs83 = NULL;
  const mxArray *c3_rhs84 = NULL;
  const mxArray *c3_lhs84 = NULL;
  const mxArray *c3_rhs85 = NULL;
  const mxArray *c3_lhs85 = NULL;
  const mxArray *c3_rhs86 = NULL;
  const mxArray *c3_lhs86 = NULL;
  const mxArray *c3_rhs87 = NULL;
  const mxArray *c3_lhs87 = NULL;
  const mxArray *c3_rhs88 = NULL;
  const mxArray *c3_lhs88 = NULL;
  const mxArray *c3_rhs89 = NULL;
  const mxArray *c3_lhs89 = NULL;
  const mxArray *c3_rhs90 = NULL;
  const mxArray *c3_lhs90 = NULL;
  const mxArray *c3_rhs91 = NULL;
  const mxArray *c3_lhs91 = NULL;
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p"),
                  "context", "context", 64);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 64);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340320U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c3_rhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p"),
                  "context", "context", 65);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xscal"),
                  "name", "name", 65);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "resolved", "resolved", 65);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c3_rhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "context", "context", 66);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 66);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 66);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372615560U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c3_rhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "context", "context", 67);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 67);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 67);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 67);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372615560U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c3_rhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "context", "context", 68);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 68);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 68);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 68);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372615560U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c3_rhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "context", "context", 69);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 69);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 69);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 69);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376013088U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c3_rhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 70);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 70);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 70);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372614816U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c3_rhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 71);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 71);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 71);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1383909694U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c3_rhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 72);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 72);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 72);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 72);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c3_rhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 73);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("sin"), "name", "name", 73);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 73);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 73);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1343862786U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c3_rhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 74);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 74);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 74);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 74);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286851136U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c3_rhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs74), "lhs", "lhs",
                  74);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 75);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("mrdivide"), "name", "name", 75);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 75);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 75);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1388492496U), "fileTimeLo",
                  "fileTimeLo", 75);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 75);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1370042286U), "mFileTimeLo",
                  "mFileTimeLo", 75);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 75);
  sf_mex_assign(&c3_rhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs75), "rhs", "rhs",
                  75);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs75), "lhs", "lhs",
                  75);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 76);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 76);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 76);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 76);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 76);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 76);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 76);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 76);
  sf_mex_assign(&c3_rhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs76), "rhs", "rhs",
                  76);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs76), "lhs", "lhs",
                  76);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 77);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("rdivide"), "name", "name", 77);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 77);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 77);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363742680U), "fileTimeLo",
                  "fileTimeLo", 77);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 77);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 77);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 77);
  sf_mex_assign(&c3_rhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs77), "rhs", "rhs",
                  77);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs77), "lhs", "lhs",
                  77);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 78);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 78);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 78);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 78);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363743356U), "fileTimeLo",
                  "fileTimeLo", 78);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 78);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 78);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 78);
  sf_mex_assign(&c3_rhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs78), "rhs", "rhs",
                  78);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs78), "lhs", "lhs",
                  78);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 79);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 79);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 79);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 79);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286851196U), "fileTimeLo",
                  "fileTimeLo", 79);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 79);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 79);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 79);
  sf_mex_assign(&c3_rhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs79), "rhs", "rhs",
                  79);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs79), "lhs", "lhs",
                  79);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 80);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_div"), "name", "name", 80);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 80);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 80);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376013088U), "fileTimeLo",
                  "fileTimeLo", 80);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 80);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 80);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 80);
  sf_mex_assign(&c3_rhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs80), "rhs", "rhs",
                  80);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs80), "lhs", "lhs",
                  80);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 81);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("cos"), "name", "name", 81);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 81);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 81);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1343862772U), "fileTimeLo",
                  "fileTimeLo", 81);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 81);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 81);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 81);
  sf_mex_assign(&c3_rhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs81), "rhs", "rhs",
                  81);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs81), "lhs", "lhs",
                  81);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 82);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 82);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 82);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 82);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286851122U), "fileTimeLo",
                  "fileTimeLo", 82);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 82);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 82);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 82);
  sf_mex_assign(&c3_rhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs82), "rhs", "rhs",
                  82);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs82), "lhs", "lhs",
                  82);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 83);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 83);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 83);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 83);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323202978U), "fileTimeLo",
                  "fileTimeLo", 83);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 83);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 83);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 83);
  sf_mex_assign(&c3_rhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs83), "rhs", "rhs",
                  83);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs83), "lhs", "lhs",
                  83);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 84);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 84);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 84);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 84);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376013088U), "fileTimeLo",
                  "fileTimeLo", 84);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 84);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 84);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 84);
  sf_mex_assign(&c3_rhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs84), "rhs", "rhs",
                  84);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs84), "lhs", "lhs",
                  84);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 85);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  85);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 85);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 85);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1376013090U), "fileTimeLo",
                  "fileTimeLo", 85);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 85);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 85);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 85);
  sf_mex_assign(&c3_rhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs85), "rhs", "rhs",
                  85);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs85), "lhs", "lhs",
                  85);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 86);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 86);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 86);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 86);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 86);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 86);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 86);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 86);
  sf_mex_assign(&c3_rhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs86), "rhs", "rhs",
                  86);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs86), "lhs", "lhs",
                  86);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 87);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 87);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 87);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 87);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 87);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 87);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 87);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 87);
  sf_mex_assign(&c3_rhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs87), "rhs", "rhs",
                  87);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs87), "lhs", "lhs",
                  87);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 88);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 88);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 88);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 88);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 88);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 88);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 88);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 88);
  sf_mex_assign(&c3_rhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs88), "rhs", "rhs",
                  88);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs88), "lhs", "lhs",
                  88);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 89);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 89);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 89);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 89);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 89);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 89);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 89);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 89);
  sf_mex_assign(&c3_rhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs89), "rhs", "rhs",
                  89);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs89), "lhs", "lhs",
                  89);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 90);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 90);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 90);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 90);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340320U), "fileTimeLo",
                  "fileTimeLo", 90);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 90);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 90);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 90);
  sf_mex_assign(&c3_rhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs90), "rhs", "rhs",
                  90);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs90), "lhs", "lhs",
                  90);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 91);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 91);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 91);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 91);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389340322U), "fileTimeLo",
                  "fileTimeLo", 91);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 91);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 91);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 91);
  sf_mex_assign(&c3_rhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs91), "rhs", "rhs",
                  91);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs91), "lhs", "lhs",
                  91);
  sf_mex_destroy(&c3_rhs64);
  sf_mex_destroy(&c3_lhs64);
  sf_mex_destroy(&c3_rhs65);
  sf_mex_destroy(&c3_lhs65);
  sf_mex_destroy(&c3_rhs66);
  sf_mex_destroy(&c3_lhs66);
  sf_mex_destroy(&c3_rhs67);
  sf_mex_destroy(&c3_lhs67);
  sf_mex_destroy(&c3_rhs68);
  sf_mex_destroy(&c3_lhs68);
  sf_mex_destroy(&c3_rhs69);
  sf_mex_destroy(&c3_lhs69);
  sf_mex_destroy(&c3_rhs70);
  sf_mex_destroy(&c3_lhs70);
  sf_mex_destroy(&c3_rhs71);
  sf_mex_destroy(&c3_lhs71);
  sf_mex_destroy(&c3_rhs72);
  sf_mex_destroy(&c3_lhs72);
  sf_mex_destroy(&c3_rhs73);
  sf_mex_destroy(&c3_lhs73);
  sf_mex_destroy(&c3_rhs74);
  sf_mex_destroy(&c3_lhs74);
  sf_mex_destroy(&c3_rhs75);
  sf_mex_destroy(&c3_lhs75);
  sf_mex_destroy(&c3_rhs76);
  sf_mex_destroy(&c3_lhs76);
  sf_mex_destroy(&c3_rhs77);
  sf_mex_destroy(&c3_lhs77);
  sf_mex_destroy(&c3_rhs78);
  sf_mex_destroy(&c3_lhs78);
  sf_mex_destroy(&c3_rhs79);
  sf_mex_destroy(&c3_lhs79);
  sf_mex_destroy(&c3_rhs80);
  sf_mex_destroy(&c3_lhs80);
  sf_mex_destroy(&c3_rhs81);
  sf_mex_destroy(&c3_lhs81);
  sf_mex_destroy(&c3_rhs82);
  sf_mex_destroy(&c3_lhs82);
  sf_mex_destroy(&c3_rhs83);
  sf_mex_destroy(&c3_lhs83);
  sf_mex_destroy(&c3_rhs84);
  sf_mex_destroy(&c3_lhs84);
  sf_mex_destroy(&c3_rhs85);
  sf_mex_destroy(&c3_lhs85);
  sf_mex_destroy(&c3_rhs86);
  sf_mex_destroy(&c3_lhs86);
  sf_mex_destroy(&c3_rhs87);
  sf_mex_destroy(&c3_lhs87);
  sf_mex_destroy(&c3_rhs88);
  sf_mex_destroy(&c3_lhs88);
  sf_mex_destroy(&c3_rhs89);
  sf_mex_destroy(&c3_lhs89);
  sf_mex_destroy(&c3_rhs90);
  sf_mex_destroy(&c3_lhs90);
  sf_mex_destroy(&c3_rhs91);
  sf_mex_destroy(&c3_lhs91);
}

static int32_T c3_intmax(SFc3_UKFInstanceStruct *chartInstance)
{
  (void)chartInstance;
  return MAX_int32_T;
}

static void c3_eml_error(SFc3_UKFInstanceStruct *chartInstance)
{
  int32_T c3_i67;
  static char_T c3_cv0[48] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'c', 'h', 'o', 'l', '_', 'm', 'a', 't', 'r', 'i', 'x', 'M',
    'u', 's', 't', 'B', 'e', 'P', 'o', 's', 'D', 'e', 'f', 'W', 'i', 't', 'h',
    'R', 'e', 'a', 'l', 'D', 'i', 'a', 'g' };

  char_T c3_u[48];
  const mxArray *c3_y = NULL;
  (void)chartInstance;
  for (c3_i67 = 0; c3_i67 < 48; c3_i67++) {
    c3_u[c3_i67] = c3_cv0[c3_i67];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 48), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c3_y));
}

static void c3_eml_matlab_zpotrf(SFc3_UKFInstanceStruct *chartInstance, real_T
  c3_A[25], real_T c3_b_A[25], int32_T *c3_info)
{
  int32_T c3_i68;
  for (c3_i68 = 0; c3_i68 < 25; c3_i68++) {
    c3_b_A[c3_i68] = c3_A[c3_i68];
  }

  *c3_info = c3_b_eml_matlab_zpotrf(chartInstance, c3_b_A);
}

static real_T c3_eml_xdotc(SFc3_UKFInstanceStruct *chartInstance, int32_T c3_n,
  real_T c3_x[25], int32_T c3_ix0, real_T c3_y[25], int32_T c3_iy0)
{
  real_T c3_d;
  int32_T c3_b_n;
  int32_T c3_b_ix0;
  int32_T c3_b_iy0;
  int32_T c3_c_n;
  int32_T c3_c_ix0;
  int32_T c3_c_iy0;
  int32_T c3_d_n;
  int32_T c3_d_ix0;
  int32_T c3_d_iy0;
  int32_T c3_e_n;
  int32_T c3_e_ix0;
  int32_T c3_e_iy0;
  int32_T c3_ix;
  int32_T c3_iy;
  int32_T c3_f_n;
  int32_T c3_b;
  int32_T c3_b_b;
  boolean_T c3_overflow;
  int32_T c3_k;
  int32_T c3_a;
  int32_T c3_b_a;
  c3_b_n = c3_n;
  c3_b_ix0 = c3_ix0;
  c3_b_iy0 = c3_iy0;
  c3_c_n = c3_b_n;
  c3_c_ix0 = c3_b_ix0;
  c3_c_iy0 = c3_b_iy0;
  c3_d_n = c3_c_n;
  c3_d_ix0 = c3_c_ix0;
  c3_d_iy0 = c3_c_iy0;
  c3_e_n = c3_d_n;
  c3_e_ix0 = c3_d_ix0;
  c3_e_iy0 = c3_d_iy0;
  c3_d = 0.0;
  if (c3_e_n < 1) {
  } else {
    c3_ix = c3_e_ix0;
    c3_iy = c3_e_iy0;
    c3_f_n = c3_e_n;
    c3_b = c3_f_n;
    c3_b_b = c3_b;
    if (1 > c3_b_b) {
      c3_overflow = false;
    } else {
      c3_overflow = (c3_b_b > c3_intmax(chartInstance) - 1);
    }

    if (c3_overflow) {
      c3_check_forloop_overflow_error(chartInstance, c3_overflow);
    }

    for (c3_k = 1; c3_k <= c3_f_n; c3_k++) {
      c3_d += c3_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c3_ix), 1, 25, 1, 0) - 1] * c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK
        ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c3_iy), 1, 25, 1, 0) - 1];
      c3_a = c3_ix + 5;
      c3_ix = c3_a;
      c3_b_a = c3_iy + 5;
      c3_iy = c3_b_a;
    }
  }

  return c3_d;
}

static void c3_check_forloop_overflow_error(SFc3_UKFInstanceStruct
  *chartInstance, boolean_T c3_overflow)
{
  int32_T c3_i69;
  static char_T c3_cv1[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c3_u[34];
  const mxArray *c3_y = NULL;
  int32_T c3_i70;
  static char_T c3_cv2[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c3_b_u[23];
  const mxArray *c3_b_y = NULL;
  (void)chartInstance;
  if (!c3_overflow) {
  } else {
    for (c3_i69 = 0; c3_i69 < 34; c3_i69++) {
      c3_u[c3_i69] = c3_cv1[c3_i69];
    }

    c3_y = NULL;
    sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 34),
                  false);
    for (c3_i70 = 0; c3_i70 < 23; c3_i70++) {
      c3_b_u[c3_i70] = c3_cv2[c3_i70];
    }

    c3_b_y = NULL;
    sf_mex_assign(&c3_b_y, sf_mex_create("y", c3_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                  false);
    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
      1U, 2U, 14, c3_y, 14, c3_b_y));
  }
}

static void c3_eml_xgemv(SFc3_UKFInstanceStruct *chartInstance, int32_T c3_m,
  int32_T c3_n, int32_T c3_ia0, int32_T c3_ix0, real_T c3_y[25], int32_T c3_iy0,
  real_T c3_b_y[25])
{
  int32_T c3_i71;
  for (c3_i71 = 0; c3_i71 < 25; c3_i71++) {
    c3_b_y[c3_i71] = c3_y[c3_i71];
  }

  c3_b_eml_xgemv(chartInstance, c3_m, c3_n, c3_ia0, c3_ix0, c3_b_y, c3_iy0);
}

static void c3_below_threshold(SFc3_UKFInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c3_b_eml_error(SFc3_UKFInstanceStruct *chartInstance)
{
  int32_T c3_i72;
  static char_T c3_cv3[19] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'p', 'o', 's', 'd', 'e', 'f' };

  char_T c3_u[19];
  const mxArray *c3_y = NULL;
  (void)chartInstance;
  for (c3_i72 = 0; c3_i72 < 19; c3_i72++) {
    c3_u[c3_i72] = c3_cv3[c3_i72];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 19), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c3_y));
}

static void c3_eml_scalar_eg(SFc3_UKFInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c3_eml_xgemm(SFc3_UKFInstanceStruct *chartInstance, real_T c3_A[55],
  real_T c3_B[11], real_T c3_C[5], real_T c3_b_C[5])
{
  int32_T c3_i73;
  int32_T c3_i74;
  real_T c3_b_A[55];
  int32_T c3_i75;
  real_T c3_b_B[11];
  for (c3_i73 = 0; c3_i73 < 5; c3_i73++) {
    c3_b_C[c3_i73] = c3_C[c3_i73];
  }

  for (c3_i74 = 0; c3_i74 < 55; c3_i74++) {
    c3_b_A[c3_i74] = c3_A[c3_i74];
  }

  for (c3_i75 = 0; c3_i75 < 11; c3_i75++) {
    c3_b_B[c3_i75] = c3_B[c3_i75];
  }

  c3_b_eml_xgemm(chartInstance, c3_b_A, c3_b_B, c3_b_C);
}

static const mxArray *c3_f_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc3_UKFInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(int32_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static int32_T c3_f_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  int32_T c3_y;
  int32_T c3_i76;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_i76, 1, 6, 0U, 0, 0U, 0);
  c3_y = c3_i76;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_sfEvent;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  int32_T c3_y;
  SFc3_UKFInstanceStruct *chartInstance;
  chartInstance = (SFc3_UKFInstanceStruct *)chartInstanceVoid;
  c3_b_sfEvent = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_sfEvent),
    &c3_thisId);
  sf_mex_destroy(&c3_b_sfEvent);
  *(int32_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static uint8_T c3_g_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance,
  const mxArray *c3_b_is_active_c3_UKF, const char_T *c3_identifier)
{
  uint8_T c3_y;
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_is_active_c3_UKF),
    &c3_thisId);
  sf_mex_destroy(&c3_b_is_active_c3_UKF);
  return c3_y;
}

static uint8_T c3_h_emlrt_marshallIn(SFc3_UKFInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  uint8_T c3_y;
  uint8_T c3_u0;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_u0, 1, 3, 0U, 0, 0U, 0);
  c3_y = c3_u0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static int32_T c3_b_eml_matlab_zpotrf(SFc3_UKFInstanceStruct *chartInstance,
  real_T c3_A[25])
{
  int32_T c3_info;
  int32_T c3_j;
  int32_T c3_b_j;
  int32_T c3_a;
  int32_T c3_b_a;
  int32_T c3_jm1;
  int32_T c3_c_a;
  int32_T c3_d_a;
  int32_T c3_c;
  int32_T c3_e_a;
  int32_T c3_b;
  int32_T c3_f_a;
  int32_T c3_b_b;
  int32_T c3_jj;
  int32_T c3_i77;
  int32_T c3_i78;
  int32_T c3_i79;
  real_T c3_b_A[25];
  int32_T c3_i80;
  int32_T c3_i81;
  int32_T c3_i82;
  real_T c3_c_A[25];
  real_T c3_ajj;
  int32_T c3_c_b;
  int32_T c3_d_b;
  int32_T c3_nmj;
  int32_T c3_g_a;
  int32_T c3_h_a;
  int32_T c3_jp1;
  int32_T c3_i_a;
  int32_T c3_j_a;
  int32_T c3_jp1j;
  int32_T c3_b_jm1;
  int32_T c3_e_b;
  int32_T c3_f_b;
  boolean_T c3_overflow;
  int32_T c3_k;
  int32_T c3_b_k;
  int32_T c3_c_jm1;
  int32_T c3_g_b;
  int32_T c3_h_b;
  boolean_T c3_b_overflow;
  int32_T c3_c_k;
  real_T c3_y;
  real_T c3_b_y;
  real_T c3_z;
  int32_T c3_n;
  real_T c3_k_a;
  int32_T c3_ix0;
  int32_T c3_b_n;
  real_T c3_l_a;
  int32_T c3_b_ix0;
  int32_T c3_c_n;
  real_T c3_m_a;
  int32_T c3_c_ix0;
  int32_T c3_d_ix0;
  int32_T c3_n_a;
  int32_T c3_b_c;
  int32_T c3_i_b;
  int32_T c3_c_c;
  int32_T c3_o_a;
  int32_T c3_j_b;
  int32_T c3_i83;
  int32_T c3_p_a;
  int32_T c3_k_b;
  int32_T c3_q_a;
  int32_T c3_l_b;
  boolean_T c3_c_overflow;
  int32_T c3_d_k;
  int32_T c3_e_k;
  boolean_T exitg1;
  c3_info = 0;
  c3_j = 1;
  exitg1 = false;
  while ((exitg1 == false) && (c3_j < 6)) {
    c3_b_j = c3_j;
    c3_a = c3_b_j;
    c3_b_a = c3_a - 1;
    c3_jm1 = c3_b_a;
    c3_c_a = c3_jm1;
    c3_d_a = c3_c_a;
    c3_c = c3_d_a * 5;
    c3_e_a = c3_b_j;
    c3_b = c3_c;
    c3_f_a = c3_e_a;
    c3_b_b = c3_b;
    c3_jj = c3_f_a + c3_b_b;
    c3_i77 = 0;
    for (c3_i78 = 0; c3_i78 < 5; c3_i78++) {
      for (c3_i79 = 0; c3_i79 < 5; c3_i79++) {
        c3_b_A[c3_i79 + c3_i77] = c3_A[c3_i79 + c3_i77];
      }

      c3_i77 += 5;
    }

    c3_i80 = 0;
    for (c3_i81 = 0; c3_i81 < 5; c3_i81++) {
      for (c3_i82 = 0; c3_i82 < 5; c3_i82++) {
        c3_c_A[c3_i82 + c3_i80] = c3_A[c3_i82 + c3_i80];
      }

      c3_i80 += 5;
    }

    c3_ajj = c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c3_jj), 1, 25, 1, 0) - 1] - c3_eml_xdotc(chartInstance, c3_jm1,
      c3_b_A, c3_b_j, c3_c_A, c3_b_j);
    if (c3_ajj > 0.0) {
      c3_ajj = muDoubleScalarSqrt(c3_ajj);
      c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c3_jj), 1, 25, 1, 0) - 1] = c3_ajj;
      if (c3_b_j < 5) {
        c3_c_b = c3_b_j;
        c3_d_b = c3_c_b;
        c3_nmj = 5 - c3_d_b;
        c3_g_a = c3_b_j;
        c3_h_a = c3_g_a;
        c3_jp1 = c3_h_a;
        c3_i_a = c3_jj;
        c3_j_a = c3_i_a + 1;
        c3_jp1j = c3_j_a;
        c3_b_jm1 = c3_jm1;
        c3_e_b = c3_b_jm1;
        c3_f_b = c3_e_b;
        if (1 > c3_f_b) {
          c3_overflow = false;
        } else {
          c3_overflow = (c3_f_b > c3_intmax(chartInstance) - 1);
        }

        if (c3_overflow) {
          c3_check_forloop_overflow_error(chartInstance, c3_overflow);
        }

        for (c3_k = 1; c3_k <= c3_b_jm1; c3_k++) {
          c3_b_k = c3_k;
          c3_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c3_b_j), 1, 5, 1, 0) + 5 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c3_b_k), 1, 5, 2, 0) - 1)) - 1] = c3_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c3_b_j), 1, 5, 1, 0) + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c3_b_k), 1, 5, 2, 0)
               - 1)) - 1];
        }

        c3_b_eml_xgemv(chartInstance, c3_nmj, c3_jm1, c3_jp1 + 1, c3_b_j, c3_A,
                       c3_jp1j);
        c3_c_jm1 = c3_jm1;
        c3_g_b = c3_c_jm1;
        c3_h_b = c3_g_b;
        if (1 > c3_h_b) {
          c3_b_overflow = false;
        } else {
          c3_b_overflow = (c3_h_b > c3_intmax(chartInstance) - 1);
        }

        if (c3_b_overflow) {
          c3_check_forloop_overflow_error(chartInstance, c3_b_overflow);
        }

        for (c3_c_k = 1; c3_c_k <= c3_c_jm1; c3_c_k++) {
          c3_b_k = c3_c_k;
          c3_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c3_b_j), 1, 5, 1, 0) + 5 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c3_b_k), 1, 5, 2, 0) - 1)) - 1] = c3_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c3_b_j), 1, 5, 1, 0) + 5 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c3_b_k), 1, 5, 2, 0)
               - 1)) - 1];
        }

        c3_y = c3_ajj;
        c3_b_y = c3_y;
        c3_z = 1.0 / c3_b_y;
        c3_n = c3_nmj;
        c3_k_a = c3_z;
        c3_ix0 = c3_jp1j;
        c3_b_n = c3_n;
        c3_l_a = c3_k_a;
        c3_b_ix0 = c3_ix0;
        c3_below_threshold(chartInstance);
        c3_c_n = c3_b_n;
        c3_m_a = c3_l_a;
        c3_c_ix0 = c3_b_ix0;
        c3_d_ix0 = c3_c_ix0;
        c3_n_a = c3_c_n;
        c3_b_c = c3_n_a;
        c3_i_b = c3_b_c - 1;
        c3_c_c = c3_i_b;
        c3_o_a = c3_c_ix0;
        c3_j_b = c3_c_c;
        c3_i83 = c3_o_a + c3_j_b;
        c3_p_a = c3_d_ix0;
        c3_k_b = c3_i83;
        c3_q_a = c3_p_a;
        c3_l_b = c3_k_b;
        if (c3_q_a > c3_l_b) {
          c3_c_overflow = false;
        } else {
          c3_c_overflow = (c3_l_b > c3_intmax(chartInstance) - 1);
        }

        if (c3_c_overflow) {
          c3_check_forloop_overflow_error(chartInstance, c3_c_overflow);
        }

        for (c3_d_k = c3_d_ix0; c3_d_k <= c3_i83; c3_d_k++) {
          c3_e_k = c3_d_k;
          c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c3_e_k), 1, 25, 1, 0) - 1] = c3_m_a *
            c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c3_e_k), 1, 25, 1, 0) - 1];
        }
      }

      c3_j++;
    } else {
      c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c3_jj), 1, 25, 1, 0) - 1] = c3_ajj;
      c3_info = c3_b_j;
      exitg1 = true;
    }
  }

  return c3_info;
}

static void c3_b_eml_xgemv(SFc3_UKFInstanceStruct *chartInstance, int32_T c3_m,
  int32_T c3_n, int32_T c3_ia0, int32_T c3_ix0, real_T c3_y[25], int32_T c3_iy0)
{
  int32_T c3_b_m;
  int32_T c3_b_n;
  int32_T c3_b_ia0;
  int32_T c3_b_ix0;
  int32_T c3_b_iy0;
  int32_T c3_c_m;
  int32_T c3_c_n;
  real_T c3_alpha1;
  int32_T c3_c_ia0;
  int32_T c3_c_ix0;
  real_T c3_beta1;
  int32_T c3_c_iy0;
  char_T c3_TRANSA;
  int32_T c3_var;
  ptrdiff_t c3_m_t;
  int32_T c3_b_var;
  ptrdiff_t c3_n_t;
  ptrdiff_t c3_lda_t;
  ptrdiff_t c3_incx_t;
  ptrdiff_t c3_incy_t;
  double * c3_alpha1_t;
  double * c3_beta1_t;
  double * c3_yiy0_t;
  double * c3_yix0_t;
  double * c3_yia0_t;
  (void)chartInstance;
  c3_b_m = c3_m;
  c3_b_n = c3_n;
  c3_b_ia0 = c3_ia0;
  c3_b_ix0 = c3_ix0;
  c3_b_iy0 = c3_iy0;
  if (c3_b_m < 1) {
  } else if (c3_b_n < 1) {
  } else {
    c3_c_m = c3_b_m;
    c3_c_n = c3_b_n;
    c3_alpha1 = -1.0;
    c3_c_ia0 = c3_b_ia0;
    c3_c_ix0 = c3_b_ix0;
    c3_beta1 = 1.0;
    c3_c_iy0 = c3_b_iy0;
    c3_TRANSA = 'N';
    c3_var = c3_c_m;
    c3_m_t = (ptrdiff_t)(c3_var);
    c3_b_var = c3_c_n;
    c3_n_t = (ptrdiff_t)(c3_b_var);
    c3_lda_t = (ptrdiff_t)(5);
    c3_incx_t = (ptrdiff_t)(5);
    c3_incy_t = (ptrdiff_t)(1);
    c3_alpha1_t = (double *)(&c3_alpha1);
    c3_beta1_t = (double *)(&c3_beta1);
    c3_yiy0_t = (double *)(&c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c3_c_iy0), 1, 25, 1, 0) - 1]);
    c3_yix0_t = (double *)(&c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c3_c_ix0), 1, 25, 1, 0) - 1]);
    c3_yia0_t = (double *)(&c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c3_c_ia0), 1, 25, 1, 0) - 1]);
    dgemv(&c3_TRANSA, &c3_m_t, &c3_n_t, c3_alpha1_t, c3_yia0_t, &c3_lda_t,
          c3_yix0_t, &c3_incx_t, c3_beta1_t, c3_yiy0_t, &c3_incy_t);
  }
}

static void c3_b_eml_xgemm(SFc3_UKFInstanceStruct *chartInstance, real_T c3_A[55],
  real_T c3_B[11], real_T c3_C[5])
{
  int32_T c3_i84;
  int32_T c3_i85;
  int32_T c3_i86;
  (void)chartInstance;
  for (c3_i84 = 0; c3_i84 < 5; c3_i84++) {
    c3_C[c3_i84] = 0.0;
    c3_i85 = 0;
    for (c3_i86 = 0; c3_i86 < 11; c3_i86++) {
      c3_C[c3_i84] += c3_A[c3_i85 + c3_i84] * c3_B[c3_i86];
      c3_i85 += 5;
    }
  }
}

static void init_dsm_address_info(SFc3_UKFInstanceStruct *chartInstance)
{
  (void)chartInstance;
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

void sf_c3_UKF_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2769288693U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(884653730U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(530377896U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1763521019U);
}

mxArray *sf_c3_UKF_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("NnWxH2s1a8ym0gqUZCFC6C");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,5,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(11);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));
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

mxArray *sf_c3_UKF_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c3_UKF_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c3_UKF(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[9],T\"x_m_p\",},{M[8],M[0],T\"is_active_c3_UKF\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c3_UKF_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc3_UKFInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc3_UKFInstanceStruct *) chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _UKFMachineNumber_,
           3,
           1,
           1,
           0,
           6,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize ist own list of scripts */
        init_script_number_translation(_UKFMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_UKFMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _UKFMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"T1");
          _SFD_SET_DATA_PROPS(1,1,1,0,"x");
          _SFD_SET_DATA_PROPS(2,2,0,1,"x_m_p");
          _SFD_SET_DATA_PROPS(3,1,1,0,"P");
          _SFD_SET_DATA_PROPS(4,1,1,0,"c");
          _SFD_SET_DATA_PROPS(5,1,1,0,"Wm");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,690);
        _SFD_CV_INIT_EML_FOR(0,1,0,178,192,653);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_c_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 5;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 5;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)
            c3_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 5;
          dimVector[1]= 5;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_c_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 11;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T *c3_T1;
          real_T *c3_c;
          real_T (*c3_x)[5];
          real_T (*c3_x_m_p)[5];
          real_T (*c3_P)[25];
          real_T (*c3_Wm)[11];
          c3_Wm = (real_T (*)[11])ssGetInputPortSignal(chartInstance->S, 4);
          c3_c = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
          c3_P = (real_T (*)[25])ssGetInputPortSignal(chartInstance->S, 2);
          c3_x_m_p = (real_T (*)[5])ssGetOutputPortSignal(chartInstance->S, 1);
          c3_x = (real_T (*)[5])ssGetInputPortSignal(chartInstance->S, 1);
          c3_T1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, c3_T1);
          _SFD_SET_DATA_VALUE_PTR(1U, *c3_x);
          _SFD_SET_DATA_VALUE_PTR(2U, *c3_x_m_p);
          _SFD_SET_DATA_VALUE_PTR(3U, *c3_P);
          _SFD_SET_DATA_VALUE_PTR(4U, c3_c);
          _SFD_SET_DATA_VALUE_PTR(5U, *c3_Wm);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _UKFMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "sDKT22WCGOSegTvLM86zGC";
}

static void sf_opaque_initialize_c3_UKF(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc3_UKFInstanceStruct*) chartInstanceVar)->S,0);
  initialize_params_c3_UKF((SFc3_UKFInstanceStruct*) chartInstanceVar);
  initialize_c3_UKF((SFc3_UKFInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c3_UKF(void *chartInstanceVar)
{
  enable_c3_UKF((SFc3_UKFInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c3_UKF(void *chartInstanceVar)
{
  disable_c3_UKF((SFc3_UKFInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c3_UKF(void *chartInstanceVar)
{
  sf_gateway_c3_UKF((SFc3_UKFInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c3_UKF(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c3_UKF((SFc3_UKFInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c3_UKF();/* state var info */
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

extern void sf_internal_set_sim_state_c3_UKF(SimStruct* S, const mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c3_UKF();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c3_UKF((SFc3_UKFInstanceStruct*)chartInfo->chartInstance,
                       mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c3_UKF(SimStruct* S)
{
  return sf_internal_get_sim_state_c3_UKF(S);
}

static void sf_opaque_set_sim_state_c3_UKF(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c3_UKF(S, st);
}

static void sf_opaque_terminate_c3_UKF(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc3_UKFInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_UKF_optimization_info();
    }

    finalize_c3_UKF((SFc3_UKFInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc3_UKF((SFc3_UKFInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c3_UKF(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    initialize_params_c3_UKF((SFc3_UKFInstanceStruct*)(chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c3_UKF(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_UKF_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,3);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,3,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,3,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,3);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,3,5);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,3,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 5; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,3);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1028761753U));
  ssSetChecksum1(S,(3169988463U));
  ssSetChecksum2(S,(4011782949U));
  ssSetChecksum3(S,(1103032927U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c3_UKF(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c3_UKF(SimStruct *S)
{
  SFc3_UKFInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc3_UKFInstanceStruct *)utMalloc(sizeof
    (SFc3_UKFInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc3_UKFInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c3_UKF;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c3_UKF;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c3_UKF;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c3_UKF;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c3_UKF;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c3_UKF;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c3_UKF;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c3_UKF;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c3_UKF;
  chartInstance->chartInfo.mdlStart = mdlStart_c3_UKF;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c3_UKF;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  crtInfo->instanceInfo = (&(chartInstance->chartInfo));
  crtInfo->isJITEnabled = false;
  ssSetUserData(S,(void *)(crtInfo));  /* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c3_UKF_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c3_UKF(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c3_UKF(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c3_UKF(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c3_UKF_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}