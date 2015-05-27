/* Include files */

#include "IMM_UKF_sfun.h"
#include "IMM_UKF_sfun_debug_macros.h"
#include "c2_IMM_UKF.h"
#include "c8_IMM_UKF.h"
#include "c9_IMM_UKF.h"
#include "c10_IMM_UKF.h"
#include "c11_IMM_UKF.h"
#include "c12_IMM_UKF.h"
#include "c13_IMM_UKF.h"
#include "c14_IMM_UKF.h"
#include "c15_IMM_UKF.h"
#include "c16_IMM_UKF.h"
#include "c17_IMM_UKF.h"
#include "c18_IMM_UKF.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
uint32_T _IMM_UKFMachineNumber_;
real_T _sfTime_;

/* Function Declarations */

/* Function Definitions */
void IMM_UKF_initializer(void)
{
}

void IMM_UKF_terminator(void)
{
}

/* SFunction Glue Code */
unsigned int sf_IMM_UKF_method_dispatcher(SimStruct *simstructPtr, unsigned int
  chartFileNumber, const char* specsCksum, int_T method, void *data)
{
  if (chartFileNumber==2) {
    c2_IMM_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==8) {
    c8_IMM_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==9) {
    c9_IMM_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==10) {
    c10_IMM_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==11) {
    c11_IMM_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==12) {
    c12_IMM_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==13) {
    c13_IMM_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==14) {
    c14_IMM_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==15) {
    c15_IMM_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==16) {
    c16_IMM_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==17) {
    c17_IMM_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==18) {
    c18_IMM_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  return 0;
}

unsigned int sf_IMM_UKF_process_check_sum_call( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[20];
  if (nrhs<1 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the checksum */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"sf_get_check_sum"))
    return 0;
  plhs[0] = mxCreateDoubleMatrix( 1,4,mxREAL);
  if (nrhs>1 && mxIsChar(prhs[1])) {
    mxGetString(prhs[1], commandName,sizeof(commandName)/sizeof(char));
    commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
    if (!strcmp(commandName,"machine")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(951495287U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3848696068U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(73294267U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1321277321U);
    } else if (!strcmp(commandName,"exportedFcn")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0U);
    } else if (!strcmp(commandName,"makefile")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1102473787U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3909254941U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(236694278U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1501946943U);
    } else if (nrhs==3 && !strcmp(commandName,"chart")) {
      unsigned int chartFileNumber;
      chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
      switch (chartFileNumber) {
       case 2:
        {
          extern void sf_c2_IMM_UKF_get_check_sum(mxArray *plhs[]);
          sf_c2_IMM_UKF_get_check_sum(plhs);
          break;
        }

       case 8:
        {
          extern void sf_c8_IMM_UKF_get_check_sum(mxArray *plhs[]);
          sf_c8_IMM_UKF_get_check_sum(plhs);
          break;
        }

       case 9:
        {
          extern void sf_c9_IMM_UKF_get_check_sum(mxArray *plhs[]);
          sf_c9_IMM_UKF_get_check_sum(plhs);
          break;
        }

       case 10:
        {
          extern void sf_c10_IMM_UKF_get_check_sum(mxArray *plhs[]);
          sf_c10_IMM_UKF_get_check_sum(plhs);
          break;
        }

       case 11:
        {
          extern void sf_c11_IMM_UKF_get_check_sum(mxArray *plhs[]);
          sf_c11_IMM_UKF_get_check_sum(plhs);
          break;
        }

       case 12:
        {
          extern void sf_c12_IMM_UKF_get_check_sum(mxArray *plhs[]);
          sf_c12_IMM_UKF_get_check_sum(plhs);
          break;
        }

       case 13:
        {
          extern void sf_c13_IMM_UKF_get_check_sum(mxArray *plhs[]);
          sf_c13_IMM_UKF_get_check_sum(plhs);
          break;
        }

       case 14:
        {
          extern void sf_c14_IMM_UKF_get_check_sum(mxArray *plhs[]);
          sf_c14_IMM_UKF_get_check_sum(plhs);
          break;
        }

       case 15:
        {
          extern void sf_c15_IMM_UKF_get_check_sum(mxArray *plhs[]);
          sf_c15_IMM_UKF_get_check_sum(plhs);
          break;
        }

       case 16:
        {
          extern void sf_c16_IMM_UKF_get_check_sum(mxArray *plhs[]);
          sf_c16_IMM_UKF_get_check_sum(plhs);
          break;
        }

       case 17:
        {
          extern void sf_c17_IMM_UKF_get_check_sum(mxArray *plhs[]);
          sf_c17_IMM_UKF_get_check_sum(plhs);
          break;
        }

       case 18:
        {
          extern void sf_c18_IMM_UKF_get_check_sum(mxArray *plhs[]);
          sf_c18_IMM_UKF_get_check_sum(plhs);
          break;
        }

       default:
        ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0.0);
      }
    } else if (!strcmp(commandName,"target")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2083502392U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1110276785U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3258378658U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3926592909U);
    } else {
      return 0;
    }
  } else {
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2827851971U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3215987912U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(36673005U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(887203722U);
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_IMM_UKF_autoinheritance_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[32];
  char aiChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the autoinheritance_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_autoinheritance_info"))
    return 0;
  mxGetString(prhs[2], aiChksum,sizeof(aiChksum)/sizeof(char));
  aiChksum[(sizeof(aiChksum)/sizeof(char)-1)] = '\0';

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 2:
      {
        if (strcmp(aiChksum, "Q3QtH4hS6Xkbjurw6cr5pB") == 0) {
          extern mxArray *sf_c2_IMM_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c2_IMM_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 8:
      {
        if (strcmp(aiChksum, "KW1AE7NhFgv790h5fIyQK") == 0) {
          extern mxArray *sf_c8_IMM_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c8_IMM_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 9:
      {
        if (strcmp(aiChksum, "KW1AE7NhFgv790h5fIyQK") == 0) {
          extern mxArray *sf_c9_IMM_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c9_IMM_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 10:
      {
        if (strcmp(aiChksum, "KW1AE7NhFgv790h5fIyQK") == 0) {
          extern mxArray *sf_c10_IMM_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c10_IMM_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 11:
      {
        if (strcmp(aiChksum, "KW1AE7NhFgv790h5fIyQK") == 0) {
          extern mxArray *sf_c11_IMM_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c11_IMM_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 12:
      {
        if (strcmp(aiChksum, "KW1AE7NhFgv790h5fIyQK") == 0) {
          extern mxArray *sf_c12_IMM_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c12_IMM_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 13:
      {
        if (strcmp(aiChksum, "yrEnNCvpLcvz83wAzYO5rF") == 0) {
          extern mxArray *sf_c13_IMM_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c13_IMM_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 14:
      {
        if (strcmp(aiChksum, "wKnG7el3bKfdknIUb6ofXG") == 0) {
          extern mxArray *sf_c14_IMM_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c14_IMM_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 15:
      {
        if (strcmp(aiChksum, "wKnG7el3bKfdknIUb6ofXG") == 0) {
          extern mxArray *sf_c15_IMM_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c15_IMM_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 16:
      {
        if (strcmp(aiChksum, "wKnG7el3bKfdknIUb6ofXG") == 0) {
          extern mxArray *sf_c16_IMM_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c16_IMM_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 17:
      {
        if (strcmp(aiChksum, "wKnG7el3bKfdknIUb6ofXG") == 0) {
          extern mxArray *sf_c17_IMM_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c17_IMM_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 18:
      {
        if (strcmp(aiChksum, "wKnG7el3bKfdknIUb6ofXG") == 0) {
          extern mxArray *sf_c18_IMM_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c18_IMM_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_IMM_UKF_get_eml_resolved_functions_info( int nlhs, mxArray *
  plhs[], int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[64];
  if (nrhs<2 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the get_eml_resolved_functions_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_eml_resolved_functions_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 2:
      {
        extern const mxArray *sf_c2_IMM_UKF_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c2_IMM_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 8:
      {
        extern const mxArray *sf_c8_IMM_UKF_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c8_IMM_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 9:
      {
        extern const mxArray *sf_c9_IMM_UKF_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c9_IMM_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 10:
      {
        extern const mxArray *sf_c10_IMM_UKF_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c10_IMM_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 11:
      {
        extern const mxArray *sf_c11_IMM_UKF_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c11_IMM_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 12:
      {
        extern const mxArray *sf_c12_IMM_UKF_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c12_IMM_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 13:
      {
        extern const mxArray *sf_c13_IMM_UKF_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c13_IMM_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 14:
      {
        extern const mxArray *sf_c14_IMM_UKF_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c14_IMM_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 15:
      {
        extern const mxArray *sf_c15_IMM_UKF_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c15_IMM_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 16:
      {
        extern const mxArray *sf_c16_IMM_UKF_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c16_IMM_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 17:
      {
        extern const mxArray *sf_c17_IMM_UKF_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c17_IMM_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 18:
      {
        extern const mxArray *sf_c18_IMM_UKF_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c18_IMM_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_IMM_UKF_third_party_uses_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the third_party_uses_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_third_party_uses_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 2:
      {
        if (strcmp(tpChksum, "t1N8OdSFit5mbzHJx6JtUF") == 0) {
          extern mxArray *sf_c2_IMM_UKF_third_party_uses_info(void);
          plhs[0] = sf_c2_IMM_UKF_third_party_uses_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "P1t2z1hzLvOYcOWXlUGyb") == 0) {
          extern mxArray *sf_c8_IMM_UKF_third_party_uses_info(void);
          plhs[0] = sf_c8_IMM_UKF_third_party_uses_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "P1t2z1hzLvOYcOWXlUGyb") == 0) {
          extern mxArray *sf_c9_IMM_UKF_third_party_uses_info(void);
          plhs[0] = sf_c9_IMM_UKF_third_party_uses_info();
          break;
        }
      }

     case 10:
      {
        if (strcmp(tpChksum, "P1t2z1hzLvOYcOWXlUGyb") == 0) {
          extern mxArray *sf_c10_IMM_UKF_third_party_uses_info(void);
          plhs[0] = sf_c10_IMM_UKF_third_party_uses_info();
          break;
        }
      }

     case 11:
      {
        if (strcmp(tpChksum, "P1t2z1hzLvOYcOWXlUGyb") == 0) {
          extern mxArray *sf_c11_IMM_UKF_third_party_uses_info(void);
          plhs[0] = sf_c11_IMM_UKF_third_party_uses_info();
          break;
        }
      }

     case 12:
      {
        if (strcmp(tpChksum, "P1t2z1hzLvOYcOWXlUGyb") == 0) {
          extern mxArray *sf_c12_IMM_UKF_third_party_uses_info(void);
          plhs[0] = sf_c12_IMM_UKF_third_party_uses_info();
          break;
        }
      }

     case 13:
      {
        if (strcmp(tpChksum, "3o002jr5JvntFmTCrpbBUH") == 0) {
          extern mxArray *sf_c13_IMM_UKF_third_party_uses_info(void);
          plhs[0] = sf_c13_IMM_UKF_third_party_uses_info();
          break;
        }
      }

     case 14:
      {
        if (strcmp(tpChksum, "M76TmyorBx1UKnd1ZOlkbF") == 0) {
          extern mxArray *sf_c14_IMM_UKF_third_party_uses_info(void);
          plhs[0] = sf_c14_IMM_UKF_third_party_uses_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "M76TmyorBx1UKnd1ZOlkbF") == 0) {
          extern mxArray *sf_c15_IMM_UKF_third_party_uses_info(void);
          plhs[0] = sf_c15_IMM_UKF_third_party_uses_info();
          break;
        }
      }

     case 16:
      {
        if (strcmp(tpChksum, "M76TmyorBx1UKnd1ZOlkbF") == 0) {
          extern mxArray *sf_c16_IMM_UKF_third_party_uses_info(void);
          plhs[0] = sf_c16_IMM_UKF_third_party_uses_info();
          break;
        }
      }

     case 17:
      {
        if (strcmp(tpChksum, "M76TmyorBx1UKnd1ZOlkbF") == 0) {
          extern mxArray *sf_c17_IMM_UKF_third_party_uses_info(void);
          plhs[0] = sf_c17_IMM_UKF_third_party_uses_info();
          break;
        }
      }

     case 18:
      {
        if (strcmp(tpChksum, "M76TmyorBx1UKnd1ZOlkbF") == 0) {
          extern mxArray *sf_c18_IMM_UKF_third_party_uses_info(void);
          plhs[0] = sf_c18_IMM_UKF_third_party_uses_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

unsigned int sf_IMM_UKF_updateBuildInfo_args_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the updateBuildInfo_args_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_updateBuildInfo_args_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 2:
      {
        if (strcmp(tpChksum, "t1N8OdSFit5mbzHJx6JtUF") == 0) {
          extern mxArray *sf_c2_IMM_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c2_IMM_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "P1t2z1hzLvOYcOWXlUGyb") == 0) {
          extern mxArray *sf_c8_IMM_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c8_IMM_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "P1t2z1hzLvOYcOWXlUGyb") == 0) {
          extern mxArray *sf_c9_IMM_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c9_IMM_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 10:
      {
        if (strcmp(tpChksum, "P1t2z1hzLvOYcOWXlUGyb") == 0) {
          extern mxArray *sf_c10_IMM_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c10_IMM_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 11:
      {
        if (strcmp(tpChksum, "P1t2z1hzLvOYcOWXlUGyb") == 0) {
          extern mxArray *sf_c11_IMM_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c11_IMM_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 12:
      {
        if (strcmp(tpChksum, "P1t2z1hzLvOYcOWXlUGyb") == 0) {
          extern mxArray *sf_c12_IMM_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c12_IMM_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 13:
      {
        if (strcmp(tpChksum, "3o002jr5JvntFmTCrpbBUH") == 0) {
          extern mxArray *sf_c13_IMM_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c13_IMM_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 14:
      {
        if (strcmp(tpChksum, "M76TmyorBx1UKnd1ZOlkbF") == 0) {
          extern mxArray *sf_c14_IMM_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c14_IMM_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "M76TmyorBx1UKnd1ZOlkbF") == 0) {
          extern mxArray *sf_c15_IMM_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c15_IMM_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 16:
      {
        if (strcmp(tpChksum, "M76TmyorBx1UKnd1ZOlkbF") == 0) {
          extern mxArray *sf_c16_IMM_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c16_IMM_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 17:
      {
        if (strcmp(tpChksum, "M76TmyorBx1UKnd1ZOlkbF") == 0) {
          extern mxArray *sf_c17_IMM_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c17_IMM_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 18:
      {
        if (strcmp(tpChksum, "M76TmyorBx1UKnd1ZOlkbF") == 0) {
          extern mxArray *sf_c18_IMM_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c18_IMM_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

void IMM_UKF_debug_initialize(struct SfDebugInstanceStruct* debugInstance)
{
  _IMM_UKFMachineNumber_ = sf_debug_initialize_machine(debugInstance,"IMM_UKF",
    "sfun",0,12,0,0,0);
  sf_debug_set_machine_event_thresholds(debugInstance,_IMM_UKFMachineNumber_,0,0);
  sf_debug_set_machine_data_thresholds(debugInstance,_IMM_UKFMachineNumber_,0);
}

void IMM_UKF_register_exported_symbols(SimStruct* S)
{
}

static mxArray* sRtwOptimizationInfoStruct= NULL;
mxArray* load_IMM_UKF_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct==NULL) {
    sRtwOptimizationInfoStruct = sf_load_rtw_optimization_info("IMM_UKF",
      "IMM_UKF");
    mexMakeArrayPersistent(sRtwOptimizationInfoStruct);
  }

  return(sRtwOptimizationInfoStruct);
}

void unload_IMM_UKF_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct!=NULL) {
    mxDestroyArray(sRtwOptimizationInfoStruct);
    sRtwOptimizationInfoStruct = NULL;
  }
}
