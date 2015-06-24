/* Include files */

#include "UKF_sfun.h"
#include "UKF_sfun_debug_macros.h"
#include "c2_UKF.h"
#include "c3_UKF.h"
#include "c4_UKF.h"
#include "c5_UKF.h"
#include "c6_UKF.h"
#include "c7_UKF.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
uint32_T _UKFMachineNumber_;

/* Function Declarations */

/* Function Definitions */
void UKF_initializer(void)
{
}

void UKF_terminator(void)
{
}

/* SFunction Glue Code */
unsigned int sf_UKF_method_dispatcher(SimStruct *simstructPtr, unsigned int
  chartFileNumber, const char* specsCksum, int_T method, void *data)
{
  if (chartFileNumber==2) {
    c2_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==3) {
    c3_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==4) {
    c4_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==5) {
    c5_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==6) {
    c6_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==7) {
    c7_UKF_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  return 0;
}

unsigned int sf_UKF_process_check_sum_call( int nlhs, mxArray * plhs[], int nrhs,
  const mxArray * prhs[] )
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
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(4021368973U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1756586097U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3193271797U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(998575220U);
    } else if (!strcmp(commandName,"exportedFcn")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0U);
    } else if (!strcmp(commandName,"makefile")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3494672649U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(4068180583U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(517701289U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(4145179517U);
    } else if (nrhs==3 && !strcmp(commandName,"chart")) {
      unsigned int chartFileNumber;
      chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
      switch (chartFileNumber) {
       case 2:
        {
          extern void sf_c2_UKF_get_check_sum(mxArray *plhs[]);
          sf_c2_UKF_get_check_sum(plhs);
          break;
        }

       case 3:
        {
          extern void sf_c3_UKF_get_check_sum(mxArray *plhs[]);
          sf_c3_UKF_get_check_sum(plhs);
          break;
        }

       case 4:
        {
          extern void sf_c4_UKF_get_check_sum(mxArray *plhs[]);
          sf_c4_UKF_get_check_sum(plhs);
          break;
        }

       case 5:
        {
          extern void sf_c5_UKF_get_check_sum(mxArray *plhs[]);
          sf_c5_UKF_get_check_sum(plhs);
          break;
        }

       case 6:
        {
          extern void sf_c6_UKF_get_check_sum(mxArray *plhs[]);
          sf_c6_UKF_get_check_sum(plhs);
          break;
        }

       case 7:
        {
          extern void sf_c7_UKF_get_check_sum(mxArray *plhs[]);
          sf_c7_UKF_get_check_sum(plhs);
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
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(406623480U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2102367403U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3126301061U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3209822257U);
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_UKF_autoinheritance_info( int nlhs, mxArray * plhs[], int nrhs,
  const mxArray * prhs[] )
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
        if (strcmp(aiChksum, "4cgweEmHYemb3D1CEpNovG") == 0) {
          extern mxArray *sf_c2_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c2_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 3:
      {
        if (strcmp(aiChksum, "NnWxH2s1a8ym0gqUZCFC6C") == 0) {
          extern mxArray *sf_c3_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c3_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 4:
      {
        if (strcmp(aiChksum, "NnWxH2s1a8ym0gqUZCFC6C") == 0) {
          extern mxArray *sf_c4_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c4_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 5:
      {
        if (strcmp(aiChksum, "NnWxH2s1a8ym0gqUZCFC6C") == 0) {
          extern mxArray *sf_c5_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c5_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 6:
      {
        if (strcmp(aiChksum, "NnWxH2s1a8ym0gqUZCFC6C") == 0) {
          extern mxArray *sf_c6_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c6_UKF_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 7:
      {
        if (strcmp(aiChksum, "NnWxH2s1a8ym0gqUZCFC6C") == 0) {
          extern mxArray *sf_c7_UKF_get_autoinheritance_info(void);
          plhs[0] = sf_c7_UKF_get_autoinheritance_info();
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

unsigned int sf_UKF_get_eml_resolved_functions_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
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
        extern const mxArray *sf_c2_UKF_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c2_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 3:
      {
        extern const mxArray *sf_c3_UKF_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c3_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 4:
      {
        extern const mxArray *sf_c4_UKF_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c4_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 5:
      {
        extern const mxArray *sf_c5_UKF_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c5_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 6:
      {
        extern const mxArray *sf_c6_UKF_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c6_UKF_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 7:
      {
        extern const mxArray *sf_c7_UKF_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c7_UKF_get_eml_resolved_functions_info();
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

unsigned int sf_UKF_third_party_uses_info( int nlhs, mxArray * plhs[], int nrhs,
  const mxArray * prhs[] )
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
        if (strcmp(tpChksum, "9L0J2D7PrXRe5rd100cCcF") == 0) {
          extern mxArray *sf_c2_UKF_third_party_uses_info(void);
          plhs[0] = sf_c2_UKF_third_party_uses_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "sDKT22WCGOSegTvLM86zGC") == 0) {
          extern mxArray *sf_c3_UKF_third_party_uses_info(void);
          plhs[0] = sf_c3_UKF_third_party_uses_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "sDKT22WCGOSegTvLM86zGC") == 0) {
          extern mxArray *sf_c4_UKF_third_party_uses_info(void);
          plhs[0] = sf_c4_UKF_third_party_uses_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "sDKT22WCGOSegTvLM86zGC") == 0) {
          extern mxArray *sf_c5_UKF_third_party_uses_info(void);
          plhs[0] = sf_c5_UKF_third_party_uses_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "sDKT22WCGOSegTvLM86zGC") == 0) {
          extern mxArray *sf_c6_UKF_third_party_uses_info(void);
          plhs[0] = sf_c6_UKF_third_party_uses_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "sDKT22WCGOSegTvLM86zGC") == 0) {
          extern mxArray *sf_c7_UKF_third_party_uses_info(void);
          plhs[0] = sf_c7_UKF_third_party_uses_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

unsigned int sf_UKF_updateBuildInfo_args_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
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
        if (strcmp(tpChksum, "9L0J2D7PrXRe5rd100cCcF") == 0) {
          extern mxArray *sf_c2_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c2_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "sDKT22WCGOSegTvLM86zGC") == 0) {
          extern mxArray *sf_c3_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c3_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "sDKT22WCGOSegTvLM86zGC") == 0) {
          extern mxArray *sf_c4_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c4_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "sDKT22WCGOSegTvLM86zGC") == 0) {
          extern mxArray *sf_c5_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c5_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "sDKT22WCGOSegTvLM86zGC") == 0) {
          extern mxArray *sf_c6_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c6_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "sDKT22WCGOSegTvLM86zGC") == 0) {
          extern mxArray *sf_c7_UKF_updateBuildInfo_args_info(void);
          plhs[0] = sf_c7_UKF_updateBuildInfo_args_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

void UKF_debug_initialize(struct SfDebugInstanceStruct* debugInstance)
{
  _UKFMachineNumber_ = sf_debug_initialize_machine(debugInstance,"UKF","sfun",0,
    6,0,0,0);
  sf_debug_set_machine_event_thresholds(debugInstance,_UKFMachineNumber_,0,0);
  sf_debug_set_machine_data_thresholds(debugInstance,_UKFMachineNumber_,0);
}

void UKF_register_exported_symbols(SimStruct* S)
{
}

static mxArray* sRtwOptimizationInfoStruct= NULL;
mxArray* load_UKF_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct==NULL) {
    sRtwOptimizationInfoStruct = sf_load_rtw_optimization_info("UKF", "UKF");
    mexMakeArrayPersistent(sRtwOptimizationInfoStruct);
  }

  return(sRtwOptimizationInfoStruct);
}

void unload_UKF_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct!=NULL) {
    mxDestroyArray(sRtwOptimizationInfoStruct);
    sRtwOptimizationInfoStruct = NULL;
  }
}
