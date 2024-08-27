/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: VandD.h
 *
 * Code generated for Simulink model 'VandD'.
 *
 * Model version                  : 1.68
 * Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
 * C/C++ source code generated on : Sun Jan 28 01:50:34 2024
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: ARM Compatible->ARM Cortex
 * Emulation hardware selection:
 *    Differs from embedded hardware (MATLAB Host)
 * Code generation objectives:
 *    1. Execution efficiency
 *    2. RAM efficiency
 * Validation result: Not run
 */

#ifndef RTW_HEADER_VandD_h_
#define RTW_HEADER_VandD_h_
#ifndef VandD_COMMON_INCLUDES_
#define VandD_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* VandD_COMMON_INCLUDES_ */

#include <string.h>
#include <stddef.h>

/* Model Code Variants */

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

#define VandD_M                        (rtM)

/* Forward declaration for rtModel */
typedef struct tag_RTM RT_MODEL;

/* Custom Type definition for MATLAB Function: '<S5>/FWS Controller' */
#ifndef struct_tag_RLzPOlY5WWc0uSZcYYgODD
#define struct_tag_RLzPOlY5WWc0uSZcYYgODD

struct tag_RLzPOlY5WWc0uSZcYYgODD
{
  real_T f1[2];
};

#endif                                 /* struct_tag_RLzPOlY5WWc0uSZcYYgODD */

#ifndef typedef_cell_wrap_0
#define typedef_cell_wrap_0

typedef struct tag_RLzPOlY5WWc0uSZcYYgODD cell_wrap_0;

#endif                                 /* typedef_cell_wrap_0 */

/* Block signals and states (default storage) for system '<Root>/VandD' */
typedef struct {
  cell_wrap_0 Traj_loc[801];
  cell_wrap_0 PP_loc[5];
  cell_wrap_0 PP[5];
  cell_wrap_0 b;
  cell_wrap_0 c;
  cell_wrap_0 d;
  cell_wrap_0 e;
  cell_wrap_0 f;
  real_T ManualSwitch[4];              /* '<S6>/Manual Switch' */
  real_T Product9[3];                  /* '<S4>/Product9' */
  real_T Phi[4];                       /* '<S113>/MATLAB Function' */
  real_T B_c[16];
  real_T B_c_m[16];
  real_T A2[16];
  real_T A4[16];
  real_T A6[16];
  real_T Q[16];
  real_T T[16];
  real_T dv[16];
  real_T b_a[16];
  real_T cBuffer[16];
  real_T aBuffer[16];
  real_T cBuffer_c[16];
  real_T b_x[16];
  real_T x[9];
  real_T b_x_k[9];
  real_T y[9];                         /* '<S4>/MATLAB Function3' */
  real_T x_ld[5];
  real_T y_ld[5];
  real_T Ft[4];                        /* '<S4>/MATLAB Function4' */
  real_T rtb_Ft_c[4];
  real_T rtb_y_b[3];
  real_T rtb_y_p[3];
  real_T rtb_y_c[3];
  real_T rtb_y_f[3];
  real_T V[16];
  real_T A6_g[16];
  real_T work[4];
  real_T tau[3];
  real_T v[3];
  real_T work_g[4];
  real_T vx;                           /* '<S4>/Integrator' */
  real_T psi_dot;                      /* '<S4>/Integrator2' */
  real_T vy;                           /* '<S4>/Integrator1' */
  real_T zs_dot;                       /* '<S4>/Integrator14' */
  real_T Integrator10;                 /* '<S4>/Integrator10' */
  real_T Integrator12;                 /* '<S4>/Integrator12' */
  real_T Integrator3;                  /* '<S4>/Integrator3' */
  real_T Integrator8;                  /* '<S4>/Integrator8' */
  real_T NProdOut;                     /* '<S51>/NProd Out' */
  real_T FilterCoefficient;            /* '<S99>/Filter Coefficient' */
  real_T Gain6;                        /* '<S4>/Gain6' */
  real_T Integrator19;                 /* '<S4>/Integrator19' */
  real_T Gain11;                       /* '<S4>/Gain11' */
  real_T Gain3;                        /* '<S4>/Gain3' */
  real_T Integrator16;                 /* '<S4>/Integrator16' */
  real_T Gain12;                       /* '<S4>/Gain12' */
  real_T Gain5;                        /* '<S4>/Gain5' */
  real_T IProdOut;                     /* '<S45>/IProd Out' */
  real_T IntegralGain;                 /* '<S93>/Integral Gain' */
  real_T p;                            /* '<S4>/Integrator5' */
  real_T p_o;                          /* '<S4>/Integrator6' */
  real_T TransferFcn;                  /* '<S3>/Transfer Fcn' */
  real_T zu2_ddot;                     /* '<S4>/1//mu2' */
  real_T zu3_ddot;                     /* '<S4>/1//mu3' */
  real_T zu4_ddot;                     /* '<S4>/1//mu4' */
  real_T Gain10;                       /* '<S4>/Gain10' */
  real_T zu1_ddot;                     /* '<S4>/`1//mu1' */
  real_T radnoise;                     /* '<S3>/Gain' */
  real_T Sum;                          /* '<S3>/Sum' */
  real_T euler_dot1;                   /* '<S4>/MATLAB Function5' */
  real_T euler_dot2;                   /* '<S4>/MATLAB Function5' */
  real_T euler_dot3;                   /* '<S4>/MATLAB Function5' */
  real_T Dz1;                          /* '<S4>/MATLAB Function' */
  real_T Dz2;                          /* '<S4>/MATLAB Function' */
  real_T Dz3;                          /* '<S4>/MATLAB Function' */
  real_T Dz4;                          /* '<S4>/MATLAB Function' */
  real_T Z_b1;                         /* '<S113>/MATLAB Function1' */
  real_T Z_b2;                         /* '<S113>/MATLAB Function1' */
  real_T Qbk;                          /* '<S113>/MATLAB Function' */
  real_T Qkk;                          /* '<S113>/MATLAB Function' */
  real_T Qnk;                          /* '<S113>/MATLAB Function' */
  real_T TimeStampA;                   /* '<S4>/Derivative' */
  real_T LastUAtTimeA;                 /* '<S4>/Derivative' */
  real_T TimeStampB;                   /* '<S4>/Derivative' */
  real_T LastUAtTimeB;                 /* '<S4>/Derivative' */
  real_T TimeStampA_l;                 /* '<S4>/Derivative1' */
  real_T LastUAtTimeA_k;               /* '<S4>/Derivative1' */
  real_T TimeStampB_o;                 /* '<S4>/Derivative1' */
  real_T LastUAtTimeB_a;               /* '<S4>/Derivative1' */
  real_T TimeStampA_p;                 /* '<S4>/Derivative2' */
  real_T LastUAtTimeA_n;               /* '<S4>/Derivative2' */
  real_T TimeStampB_n;                 /* '<S4>/Derivative2' */
  real_T LastUAtTimeB_e;               /* '<S4>/Derivative2' */
  real_T TimeStampA_o;                 /* '<S4>/Derivative3' */
  real_T LastUAtTimeA_e;               /* '<S4>/Derivative3' */
  real_T TimeStampB_l;                 /* '<S4>/Derivative3' */
  real_T LastUAtTimeB_h;               /* '<S4>/Derivative3' */
  real_T Memory_PreviousInput;         /* '<S113>/Memory' */
  real_T Memory1_PreviousInput;        /* '<S113>/Memory1' */
  real_T absx11;
  real_T absx21;
  real_T absx31;
  real_T Derivative2;                  /* '<S4>/Derivative2' */
  real_T Derivative3;                  /* '<S4>/Derivative3' */
  real_T Fz;
  real_T rtb_Fs_idx_0;
  real_T rtb_Fs_idx_1;
  real_T rtb_Fs_idx_2;
  real_T rtb_alpha_idx_0;
  real_T rtb_alpha_idx_1;
  real_T rtb_alpha_idx_2;
  real_T rtb_alpha_idx_3;
  real_T BCD_idx_0;
  real_T BCD_idx_1;
  real_T BCD_idx_2;
  real_T Fz_idx_0;
  real_T D_idx_0;
  real_T B_idx_1;
  real_T Fz_idx_1;
  real_T D_idx_1;
  real_T B_idx_2;
  real_T Fz_idx_2;
  real_T D_idx_2;
  real_T B_idx_3;
  real_T D_tmp;
  real_T exptj;
  real_T d6;
  real_T eta1;
  real_T b_varargin_1;
  real_T e_m;
  real_T ed2;
  real_T smax;
  real_T s;
  real_T d_n;
  real_T A6_p;
  real_T A4_l;
  real_T tst;
  real_T htmp1;
  real_T ab;
  real_T ba;
  real_T aa;
  real_T h12;
  real_T sn;
  real_T tst_tmp;
  real_T tst_tmp_tmp;
  real_T temp;
  real_T p_j;
  real_T bcmax;
  real_T bcmis;
  real_T scale;
  real_T z;
  real_T a;
  real_T alpha1;
  real_T xnorm;
  real_T c_d;
  real_T xnorm_g;
  real_T a_l;
  real_T scale_d;
  real_T absxk;
  real_T t;
  real_T x_d;
  real_T c_u;
  real_T temp_l;
  real_T temp_o;
  real_T temp_b;
  int32_T ipiv[4];
  int32_T blockFormat[3];
  int32_T ipiv_n[4];
  int32_T p1;
  int32_T p2;
  int32_T p3;
  int32_T itmp;
  int32_T b_s;
  int32_T e_j;
  int32_T e_i;
  int32_T A2_tmp;
  int32_T i;
  int32_T n;
  int32_T nb;
  int32_T nbitson;
  int32_T b_n;
  int32_T i1;
  int32_T c_tmp;
  int32_T c_tmp_tmp;
  int32_T j;
  int32_T i_b;
  int32_T iaii;
  int32_T knt;
  int32_T lastc;
  int32_T iy;
  int32_T f_l;
  int32_T i_h;
  uint32_T state[2];                   /* '<S113>/MATLAB Function1' */
  uint32_T state_o[2];                 /* '<S113>/MATLAB Function1' */
  uint32_T state_d[625];               /* '<S113>/MATLAB Function1' */
  uint32_T u32[2];
  uint32_T method;                     /* '<S113>/MATLAB Function1' */
  uint32_T method_b;                   /* '<S113>/MATLAB Function1' */
  uint32_T state_g;                    /* '<S113>/MATLAB Function1' */
  int8_T p_b[4];
} DW_VandD;

/* Continuous states for system '<Root>/VandD' */
typedef struct {
  real_T Integrator20_CSTATE;          /* '<S4>/Integrator20' */
  real_T Integrator21_CSTATE;          /* '<S4>/Integrator21' */
  real_T Integrator22_CSTATE;          /* '<S4>/Integrator22' */
  real_T Integrator_CSTATE;            /* '<S4>/Integrator' */
  real_T Integrator2_CSTATE;           /* '<S4>/Integrator2' */
  real_T Integrator1_CSTATE;           /* '<S4>/Integrator1' */
  real_T Integrator15_CSTATE;          /* '<S4>/Integrator15' */
  real_T Integrator14_CSTATE;          /* '<S4>/Integrator14' */
  real_T Integrator11_CSTATE;          /* '<S4>/Integrator11' */
  real_T Integrator10_CSTATE;          /* '<S4>/Integrator10' */
  real_T Integrator13_CSTATE;          /* '<S4>/Integrator13' */
  real_T Integrator12_CSTATE;          /* '<S4>/Integrator12' */
  real_T Integrator7_CSTATE;           /* '<S4>/Integrator7' */
  real_T Integrator3_CSTATE;           /* '<S4>/Integrator3' */
  real_T Integrator9_CSTATE;           /* '<S4>/Integrator9' */
  real_T Integrator8_CSTATE;           /* '<S4>/Integrator8' */
  real_T Integrator_CSTATE_p;          /* '<S48>/Integrator' */
  real_T Filter_CSTATE;                /* '<S43>/Filter' */
  real_T Integrator_CSTATE_b;          /* '<S96>/Integrator' */
  real_T Filter_CSTATE_n;              /* '<S91>/Filter' */
  real_T Integrator19_CSTATE;          /* '<S4>/Integrator19' */
  real_T Integrator17_CSTATE;          /* '<S4>/Integrator17' */
  real_T Integrator16_CSTATE;          /* '<S4>/Integrator16' */
  real_T Integrator18_CSTATE;          /* '<S4>/Integrator18' */
  real_T Integrator5_CSTATE;           /* '<S4>/Integrator5' */
  real_T Integrator6_CSTATE;           /* '<S4>/Integrator6' */
  real_T TransferFcn_CSTATE;           /* '<S3>/Transfer Fcn' */
  real_T Integrator4_CSTATE;           /* '<S4>/Integrator4' */
  real_T Integrator23_CSTATE;          /* '<S4>/Integrator23' */
} X_VandD;

/* State derivatives for system '<Root>/VandD' */
typedef struct {
  real_T Integrator20_CSTATE;          /* '<S4>/Integrator20' */
  real_T Integrator21_CSTATE;          /* '<S4>/Integrator21' */
  real_T Integrator22_CSTATE;          /* '<S4>/Integrator22' */
  real_T Integrator_CSTATE;            /* '<S4>/Integrator' */
  real_T Integrator2_CSTATE;           /* '<S4>/Integrator2' */
  real_T Integrator1_CSTATE;           /* '<S4>/Integrator1' */
  real_T Integrator15_CSTATE;          /* '<S4>/Integrator15' */
  real_T Integrator14_CSTATE;          /* '<S4>/Integrator14' */
  real_T Integrator11_CSTATE;          /* '<S4>/Integrator11' */
  real_T Integrator10_CSTATE;          /* '<S4>/Integrator10' */
  real_T Integrator13_CSTATE;          /* '<S4>/Integrator13' */
  real_T Integrator12_CSTATE;          /* '<S4>/Integrator12' */
  real_T Integrator7_CSTATE;           /* '<S4>/Integrator7' */
  real_T Integrator3_CSTATE;           /* '<S4>/Integrator3' */
  real_T Integrator9_CSTATE;           /* '<S4>/Integrator9' */
  real_T Integrator8_CSTATE;           /* '<S4>/Integrator8' */
  real_T Integrator_CSTATE_p;          /* '<S48>/Integrator' */
  real_T Filter_CSTATE;                /* '<S43>/Filter' */
  real_T Integrator_CSTATE_b;          /* '<S96>/Integrator' */
  real_T Filter_CSTATE_n;              /* '<S91>/Filter' */
  real_T Integrator19_CSTATE;          /* '<S4>/Integrator19' */
  real_T Integrator17_CSTATE;          /* '<S4>/Integrator17' */
  real_T Integrator16_CSTATE;          /* '<S4>/Integrator16' */
  real_T Integrator18_CSTATE;          /* '<S4>/Integrator18' */
  real_T Integrator5_CSTATE;           /* '<S4>/Integrator5' */
  real_T Integrator6_CSTATE;           /* '<S4>/Integrator6' */
  real_T TransferFcn_CSTATE;           /* '<S3>/Transfer Fcn' */
  real_T Integrator4_CSTATE;           /* '<S4>/Integrator4' */
  real_T Integrator23_CSTATE;          /* '<S4>/Integrator23' */
} XDot_VandD;

/* State Disabled for system '<Root>/VandD' */
typedef struct {
  boolean_T Integrator20_CSTATE;       /* '<S4>/Integrator20' */
  boolean_T Integrator21_CSTATE;       /* '<S4>/Integrator21' */
  boolean_T Integrator22_CSTATE;       /* '<S4>/Integrator22' */
  boolean_T Integrator_CSTATE;         /* '<S4>/Integrator' */
  boolean_T Integrator2_CSTATE;        /* '<S4>/Integrator2' */
  boolean_T Integrator1_CSTATE;        /* '<S4>/Integrator1' */
  boolean_T Integrator15_CSTATE;       /* '<S4>/Integrator15' */
  boolean_T Integrator14_CSTATE;       /* '<S4>/Integrator14' */
  boolean_T Integrator11_CSTATE;       /* '<S4>/Integrator11' */
  boolean_T Integrator10_CSTATE;       /* '<S4>/Integrator10' */
  boolean_T Integrator13_CSTATE;       /* '<S4>/Integrator13' */
  boolean_T Integrator12_CSTATE;       /* '<S4>/Integrator12' */
  boolean_T Integrator7_CSTATE;        /* '<S4>/Integrator7' */
  boolean_T Integrator3_CSTATE;        /* '<S4>/Integrator3' */
  boolean_T Integrator9_CSTATE;        /* '<S4>/Integrator9' */
  boolean_T Integrator8_CSTATE;        /* '<S4>/Integrator8' */
  boolean_T Integrator_CSTATE_p;       /* '<S48>/Integrator' */
  boolean_T Filter_CSTATE;             /* '<S43>/Filter' */
  boolean_T Integrator_CSTATE_b;       /* '<S96>/Integrator' */
  boolean_T Filter_CSTATE_n;           /* '<S91>/Filter' */
  boolean_T Integrator19_CSTATE;       /* '<S4>/Integrator19' */
  boolean_T Integrator17_CSTATE;       /* '<S4>/Integrator17' */
  boolean_T Integrator16_CSTATE;       /* '<S4>/Integrator16' */
  boolean_T Integrator18_CSTATE;       /* '<S4>/Integrator18' */
  boolean_T Integrator5_CSTATE;        /* '<S4>/Integrator5' */
  boolean_T Integrator6_CSTATE;        /* '<S4>/Integrator6' */
  boolean_T TransferFcn_CSTATE;        /* '<S3>/Transfer Fcn' */
  boolean_T Integrator4_CSTATE;        /* '<S4>/Integrator4' */
  boolean_T Integrator23_CSTATE;       /* '<S4>/Integrator23' */
} XDis_VandD;

/* Block signals and states (default storage) for system '<Root>' */
typedef struct {
  DW_VandD VandD_gs;                   /* '<Root>/VandD' */
  real_T Fy[4];                        /* '<S4>/MATLAB Function2' */
} DW;

/* Continuous states (default storage) */
typedef struct {
  X_VandD VandD_gs;                    /* '<Root>/VandD' */
} X;

/* State derivatives (default storage) */
typedef struct {
  XDot_VandD VandD_gs;                 /* '<Root>/VandD' */
} XDot;

/* State disabled  */
typedef struct {
  XDis_VandD VandD_gs;                 /* '<Root>/VandD' */
} XDis;

#ifndef ODE4_INTG
#define ODE4_INTG

/* ODE4 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[4];                        /* derivatives */
} ODE4_IntgData;

#endif

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T lamda[4];                     /* '<Root>/lamda' */
} ExtU;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T Fy1;                          /* '<Root>/Fy1' */
  real_T Fy2;                          /* '<Root>/Fy2' */
  real_T Fy3;                          /* '<Root>/Fy3' */
  real_T Fy4;                          /* '<Root>/Fy4' */
  real_T throttleforce[4];             /* '<Root>/throttle force' */
  real_T Steeringwheelangle;           /* '<Root>/Steering wheel angle' */
  real_T ActualVelocity;               /* '<Root>/Actual Velocity' */
  real_T LateralAcceleration;          /* '<Root>/Lateral Acceleration' */
  real_T AngularVelocity;              /* '<Root>/Angular Velocity' */
  real_T LateralVelocity;              /* '<Root>/Lateral Velocity' */
  real_T x_position;                   /* '<Root>/x_position' */
  real_T y_position;                   /* '<Root>/y_position' */
  real_T orientation;                  /* '<Root>/orientation' */
} ExtY;

/* Parameters for system: '<Root>/VandD' */
struct P_VandD_ {
  real_T PIDController1_D;             /* Mask Parameter: PIDController1_D
                                        * Referenced by: '<S90>/Derivative Gain'
                                        */
  real_T PIDController1_I;             /* Mask Parameter: PIDController1_I
                                        * Referenced by: '<S93>/Integral Gain'
                                        */
  real_T PIDController_InitialConditionForFilter;
                      /* Mask Parameter: PIDController_InitialConditionForFilter
                       * Referenced by: '<S43>/Filter'
                       */
  real_T PIDController1_InitialConditionForFilter;
                     /* Mask Parameter: PIDController1_InitialConditionForFilter
                      * Referenced by: '<S91>/Filter'
                      */
  real_T PIDController_InitialConditionForIntegrator;
                  /* Mask Parameter: PIDController_InitialConditionForIntegrator
                   * Referenced by: '<S48>/Integrator'
                   */
  real_T PIDController1_InitialConditionForIntegrator;
                 /* Mask Parameter: PIDController1_InitialConditionForIntegrator
                  * Referenced by: '<S96>/Integrator'
                  */
  real_T PIDController1_N;             /* Mask Parameter: PIDController1_N
                                        * Referenced by: '<S99>/Filter Coefficient'
                                        */
  real_T PIDController1_P;             /* Mask Parameter: PIDController1_P
                                        * Referenced by: '<S101>/Proportional Gain'
                                        */
  real_T Out1_Y0;                      /* Computed Parameter: Out1_Y0
                                        * Referenced by: '<S10>/Out1'
                                        */
  real_T Out1_Y0_j;                    /* Computed Parameter: Out1_Y0_j
                                        * Referenced by: '<S11>/Out1'
                                        */
  real_T Out1_Y0_c;                    /* Computed Parameter: Out1_Y0_c
                                        * Referenced by: '<S12>/Out1'
                                        */
  real_T Out1_Y0_b;                    /* Computed Parameter: Out1_Y0_b
                                        * Referenced by: '<S13>/Out1'
                                        */
  real_T Integrator20_IC;              /* Expression: 0
                                        * Referenced by: '<S4>/Integrator20'
                                        */
  real_T Integrator21_IC;              /* Expression: 0
                                        * Referenced by: '<S4>/Integrator21'
                                        */
  real_T Integrator22_IC;              /* Expression: 0
                                        * Referenced by: '<S4>/Integrator22'
                                        */
  real_T hzpchzpchzpchzpc_Value[4];    /* Expression: [hzpc,hzpc,hzpc,hzpc]
                                        * Referenced by: '<S4>/[hzpc,hzpc,hzpc,hzpc]'
                                        */
  real_T l1l1l2l2_Value[4];            /* Expression: [+l1, +l1, -l2, -l2]
                                        * Referenced by: '<S4>/[l1,l1,-l2,-l2]'
                                        */
  real_T hzrc_Value[4];                /* Expression: [hzrc,hzrc,hzrc,hzrc]
                                        * Referenced by: '<S4>/hzrc'
                                        */
  real_T b1b2b3b4_Value[4];            /* Expression: [b1, -b2, b3, -b4]
                                        * Referenced by: '<S4>/b1, -b2, b3, -b4'
                                        */
  real_T b1b2b3b4_Value_f[4];          /* Expression: [+b1, -b2, +b3, -b4]
                                        * Referenced by: '<S4>/[b1, -b2, b3, -b4]'
                                        */
  real_T l1l1l2l2_Value_m[4];          /* Expression: [l1, l1, -l2, -l2]
                                        * Referenced by: '<S4>/[l1, l1, -l2, -l2]'
                                        */
  real_T b1b2b3b4_Value_j[4];          /* Expression: [-b1, +b2, -b3, +b4]
                                        * Referenced by: '<S4>/[-b1, b2, -b3, b4]'
                                        */
  real_T Integrator15_IC;              /* Expression: -0.058
                                        * Referenced by: '<S4>/Integrator15'
                                        */
  real_T Integrator14_IC;              /* Expression: 0
                                        * Referenced by: '<S4>/Integrator14'
                                        */
  real_T Integrator10_IC;              /* Expression: 0
                                        * Referenced by: '<S4>/Integrator10'
                                        */
  real_T Integrator12_IC;              /* Expression: 0
                                        * Referenced by: '<S4>/Integrator12'
                                        */
  real_T Integrator3_IC;               /* Expression: 0
                                        * Referenced by: '<S4>/Integrator3'
                                        */
  real_T Integrator8_IC;               /* Expression: 0
                                        * Referenced by: '<S4>/Integrator8'
                                        */
  real_T Constant_Value;               /* Expression: 10
                                        * Referenced by: '<S6>/Constant'
                                        */
  real_T Constant1_Value;              /* Expression: 10
                                        * Referenced by: '<S6>/Constant1'
                                        */
  real_T Integrator19_IC;              /* Expression: 0
                                        * Referenced by: '<S4>/Integrator19'
                                        */
  real_T Integrator17_IC;              /* Expression: 0
                                        * Referenced by: '<S4>/Integrator17'
                                        */
  real_T Integrator16_IC;              /* Expression: 0
                                        * Referenced by: '<S4>/Integrator16'
                                        */
  real_T Integrator18_IC;              /* Expression: 0
                                        * Referenced by: '<S4>/Integrator18'
                                        */
  real_T TransferFcn_A;                /* Computed Parameter: TransferFcn_A
                                        * Referenced by: '<S3>/Transfer Fcn'
                                        */
  real_T TransferFcn_C;                /* Computed Parameter: TransferFcn_C
                                        * Referenced by: '<S3>/Transfer Fcn'
                                        */
  real_T Memory_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S113>/Memory'
                                        */
  real_T Memory1_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S113>/Memory1'
                                        */
  real_T Gain_Gain;                    /* Expression: pi/180
                                        * Referenced by: '<S113>/Gain'
                                        */
  real_T Gain_Gain_h;                  /* Expression: 1
                                        * Referenced by: '<S3>/Gain'
                                        */
  real_T Constant_Value_l;             /* Expression: 1
                                        * Referenced by: '<S8>/Constant'
                                        */
  real_T Constant1_Value_m;            /* Expression: 1
                                        * Referenced by: '<S8>/Constant1'
                                        */
  real_T Constant_Value_c;             /* Expression: 1
                                        * Referenced by: '<S9>/Constant'
                                        */
  real_T Constant1_Value_g;            /* Expression: 1
                                        * Referenced by: '<S9>/Constant1'
                                        */
  uint8_T ManualSwitch_CurrentSetting;
                              /* Computed Parameter: ManualSwitch_CurrentSetting
                               * Referenced by: '<S6>/Manual Switch'
                               */
};

/* Parameters for system: '<Root>/VandD' */
typedef struct P_VandD_ P_VandD;

/* Parameters (default storage) */
struct P_ {
  real_T B;                            /* Variable: B
                                        * Referenced by: '<S113>/Constant3'
                                        */
  real_T Fs;                           /* Variable: Fs
                                        * Referenced by: '<S113>/Constant'
                                        */
  real_T Ix;                           /* Variable: Ix
                                        * Referenced by: '<S4>/Gain12'
                                        */
  real_T Ix_Iy;                        /* Variable: Ix_Iy
                                        * Referenced by: '<S4>/Ix-Iy'
                                        */
  real_T Iy;                           /* Variable: Iy
                                        * Referenced by: '<S4>/Gain11'
                                        */
  real_T Iy_Iz;                        /* Variable: Iy_Iz
                                        * Referenced by: '<S4>/Iy-Iz'
                                        */
  real_T Iz;                           /* Variable: Iz
                                        * Referenced by: '<S4>/Gain3'
                                        */
  real_T Iz_Ix;                        /* Variable: Iz_Ix
                                        * Referenced by: '<S4>/Iz-Ix'
                                        */
  real_T K;                            /* Variable: K
                                        * Referenced by: '<S113>/Constant2'
                                        */
  real_T Model;                        /* Variable: Model
                                        * Referenced by: '<S4>/MATLAB Function2'
                                        */
  real_T N;                            /* Variable: N
                                        * Referenced by: '<S113>/Constant1'
                                        */
  real_T THETAo;                       /* Variable: THETAo
                                        * Referenced by: '<S4>/Integrator4'
                                        */
  real_T Tc;                           /* Variable: Tc
                                        * Referenced by: '<S113>/Constant4'
                                        */
  real_T Xo;                           /* Variable: Xo
                                        * Referenced by: '<S4>/Integrator5'
                                        */
  real_T Yo;                           /* Variable: Yo
                                        * Referenced by:
                                        *   '<S4>/Integrator23'
                                        *   '<S4>/Integrator6'
                                        */
  real_T b1;                           /* Variable: b1
                                        * Referenced by: '<S4>/MATLAB Function'
                                        */
  real_T b2;                           /* Variable: b2
                                        * Referenced by: '<S4>/MATLAB Function'
                                        */
  real_T b3;                           /* Variable: b3
                                        * Referenced by: '<S4>/MATLAB Function'
                                        */
  real_T b4;                           /* Variable: b4
                                        * Referenced by: '<S4>/MATLAB Function'
                                        */
  real_T dpc;                          /* Variable: dpc
                                        * Referenced by:
                                        *   '<S4>/drc'
                                        *   '<S4>/m*g*dpc'
                                        */
  real_T drc;                          /* Variable: drc
                                        * Referenced by:
                                        *   '<S4>/dpc'
                                        *   '<S4>/m*g*drc'
                                        */
  real_T g;                            /* Variable: g
                                        * Referenced by:
                                        *   '<S4>/hzp4'
                                        *   '<S4>/hzpc1'
                                        *   '<S4>/hzpc2'
                                        *   '<S4>/hzpc3'
                                        *   '<S4>/m*g'
                                        *   '<S4>/m*g*dpc'
                                        *   '<S4>/m*g*drc'
                                        */
  real_T l1;                           /* Variable: l1
                                        * Referenced by:
                                        *   '<S4>/MATLAB Function'
                                        *   '<S4>/MATLAB Function1'
                                        *   '<S5>/FWS Controller'
                                        */
  real_T l2;                           /* Variable: l2
                                        * Referenced by:
                                        *   '<S4>/MATLAB Function'
                                        *   '<S4>/MATLAB Function1'
                                        *   '<S5>/FWS Controller'
                                        */
  real_T m;                            /* Variable: m
                                        * Referenced by:
                                        *   '<S4>/m*g'
                                        *   '<S4>/m*g*dpc'
                                        *   '<S4>/m*g*drc'
                                        *   '<S4>/Gain1'
                                        *   '<S4>/Gain10'
                                        *   '<S4>/Gain2'
                                        *   '<S4>/Gain5'
                                        *   '<S4>/Gain6'
                                        *   '<S4>/Gain8'
                                        *   '<S5>/FWS Controller'
                                        */
  real_T mu1;                          /* Variable: mu1
                                        * Referenced by:
                                        *   '<S4>/hzpc1'
                                        *   '<S4>/`1//mu1'
                                        */
  real_T mu2;                          /* Variable: mu2
                                        * Referenced by:
                                        *   '<S4>/hzpc2'
                                        *   '<S4>/1//mu2'
                                        */
  real_T mu3;                          /* Variable: mu3
                                        * Referenced by:
                                        *   '<S4>/hzpc3'
                                        *   '<S4>/1//mu3'
                                        */
  real_T mu4;                          /* Variable: mu4
                                        * Referenced by:
                                        *   '<S4>/hzp4'
                                        *   '<S4>/1//mu4'
                                        */
  real_T radius;                       /* Variable: radius
                                        * Referenced by: '<S5>/FWS Controller'
                                        */
  real_T ro;                           /* Variable: ro
                                        * Referenced by: '<S4>/Integrator2'
                                        */
  real_T sw_delta_f;                   /* Variable: sw_delta_f
                                        * Referenced by:
                                        *   '<S8>/Constant2'
                                        *   '<S9>/Constant2'
                                        */
  real_T sw_phi;                       /* Variable: sw_phi
                                        * Referenced by: '<S4>/Gain12'
                                        */
  real_T sw_psi;                       /* Variable: sw_psi
                                        * Referenced by: '<S4>/Gain3'
                                        */
  real_T sw_theta;                     /* Variable: sw_theta
                                        * Referenced by: '<S4>/Gain11'
                                        */
  real_T sw_x;                         /* Variable: sw_x
                                        * Referenced by: '<S4>/Gain6'
                                        */
  real_T sw_y;                         /* Variable: sw_y
                                        * Referenced by: '<S4>/Gain5'
                                        */
  real_T sw_zs;                        /* Variable: sw_zs
                                        * Referenced by: '<S4>/Gain10'
                                        */
  real_T sw_zu1;                       /* Variable: sw_zu1
                                        * Referenced by: '<S4>/`1//mu1'
                                        */
  real_T sw_zu2;                       /* Variable: sw_zu2
                                        * Referenced by: '<S4>/1//mu2'
                                        */
  real_T sw_zu3;                       /* Variable: sw_zu3
                                        * Referenced by: '<S4>/1//mu3'
                                        */
  real_T sw_zu4;                       /* Variable: sw_zu4
                                        * Referenced by: '<S4>/1//mu4'
                                        */
  real_T vxo;                          /* Variable: vxo
                                        * Referenced by: '<S4>/Integrator'
                                        */
  real_T vyo;                          /* Variable: vyo
                                        * Referenced by: '<S4>/Integrator1'
                                        */
  real_T x[801];                       /* Variable: x
                                        * Referenced by: '<S5>/FWS Controller'
                                        */
  real_T y[801];                       /* Variable: y
                                        * Referenced by: '<S5>/FWS Controller'
                                        */
  real_T zu1_0;                        /* Variable: zu1_0
                                        * Referenced by: '<S4>/Integrator7'
                                        */
  real_T zu2_0;                        /* Variable: zu2_0
                                        * Referenced by: '<S4>/Integrator9'
                                        */
  real_T zu3_0;                        /* Variable: zu3_0
                                        * Referenced by: '<S4>/Integrator11'
                                        */
  real_T zu4_0;                        /* Variable: zu4_0
                                        * Referenced by: '<S4>/Integrator13'
                                        */
  P_VandD VandD_gs;                    /* '<Root>/VandD' */
};

/* Parameters (default storage) */
typedef struct P_ P;

/* Real-time Model Data Structure */
struct tag_RTM {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[29];
  real_T odeF[4][29];
  ODE4_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    time_T stepSize0;
    uint32_T clockTick1;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block parameters (default storage) */
extern P rtP;

/* Continuous states (default storage) */
extern X rtX;

/* Block signals and states (default storage) */
extern DW rtDW;

/* External inputs (root inport signals with default storage) */
extern ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY rtY;

/* Model entry point functions */
extern void VandD_initialize(void);
extern void VandD_step(void);

/* Real-time Model object */
extern RT_MODEL *const rtM;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S5>/From Workspace' : Unused code path elimination
 * Block '<S8>/Product' : Unused code path elimination
 * Block '<S8>/Product1' : Unused code path elimination
 * Block '<S8>/Sum' : Unused code path elimination
 * Block '<S9>/Product' : Unused code path elimination
 * Block '<S9>/Product1' : Unused code path elimination
 * Block '<S9>/Sum' : Unused code path elimination
 * Block '<S6>/Constant2' : Unused code path elimination
 * Block '<S3>/Scope' : Unused code path elimination
 * Block '<S4>/Fs1' : Unused code path elimination
 * Block '<S4>/Fs10' : Unused code path elimination
 * Block '<S4>/Fs11' : Unused code path elimination
 * Block '<S4>/Fs12' : Unused code path elimination
 * Block '<S4>/Fs17' : Unused code path elimination
 * Block '<S4>/Fs2' : Unused code path elimination
 * Block '<S4>/Fs3' : Unused code path elimination
 * Block '<S4>/Fs4' : Unused code path elimination
 * Block '<S4>/Fs5' : Unused code path elimination
 * Block '<S4>/Fs6' : Unused code path elimination
 * Block '<S4>/Fs7' : Unused code path elimination
 * Block '<S4>/Scope' : Unused code path elimination
 * Block '<S4>/Scope1' : Unused code path elimination
 * Block '<S4>/Scope2' : Unused code path elimination
 * Block '<S4>/Scope3' : Unused code path elimination
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Note that this particular code originates from a subsystem build,
 * and has its own system numbers different from the parent model.
 * Refer to the system hierarchy for this subsystem below, and use the
 * MATLAB hilite_system command to trace the generated code back
 * to the parent model.  For example,
 *
 * hilite_system('MalaksModelSub/VandD')    - opens subsystem MalaksModelSub/VandD
 * hilite_system('MalaksModelSub/VandD/Kp') - opens and selects block Kp
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'MalaksModelSub'
 * '<S1>'   : 'MalaksModelSub/VandD'
 * '<S2>'   : 'MalaksModelSub/VandD/Driver Subsystem'
 * '<S3>'   : 'MalaksModelSub/VandD/Orientation Gyro'
 * '<S4>'   : 'MalaksModelSub/VandD/Vehicle model'
 * '<S5>'   : 'MalaksModelSub/VandD/Driver Subsystem/Driver Steering Control'
 * '<S6>'   : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control'
 * '<S7>'   : 'MalaksModelSub/VandD/Driver Subsystem/Driver Steering Control/FWS Controller'
 * '<S8>'   : 'MalaksModelSub/VandD/Driver Subsystem/Driver Steering Control/Subsystem'
 * '<S9>'   : 'MalaksModelSub/VandD/Driver Subsystem/Driver Steering Control/Subsystem1'
 * '<S10>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Steering Control/Subsystem/If Action Subsystem'
 * '<S11>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Steering Control/Subsystem/If Action Subsystem1'
 * '<S12>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Steering Control/Subsystem1/If Action Subsystem'
 * '<S13>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Steering Control/Subsystem1/If Action Subsystem1'
 * '<S14>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/MATLAB Function'
 * '<S15>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller'
 * '<S16>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1'
 * '<S17>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Anti-windup'
 * '<S18>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/D Gain'
 * '<S19>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Filter'
 * '<S20>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Filter ICs'
 * '<S21>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/I Gain'
 * '<S22>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Ideal P Gain'
 * '<S23>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Ideal P Gain Fdbk'
 * '<S24>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Integrator'
 * '<S25>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Integrator ICs'
 * '<S26>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/N Copy'
 * '<S27>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/N Gain'
 * '<S28>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/P Copy'
 * '<S29>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Parallel P Gain'
 * '<S30>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Reset Signal'
 * '<S31>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Saturation'
 * '<S32>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Saturation Fdbk'
 * '<S33>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Sum'
 * '<S34>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Sum Fdbk'
 * '<S35>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Tracking Mode'
 * '<S36>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Tracking Mode Sum'
 * '<S37>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Tsamp - Integral'
 * '<S38>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Tsamp - Ngain'
 * '<S39>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/postSat Signal'
 * '<S40>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/preSat Signal'
 * '<S41>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Anti-windup/Passthrough'
 * '<S42>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/D Gain/External Parameters'
 * '<S43>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Filter/Cont. Filter'
 * '<S44>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Filter ICs/Internal IC - Filter'
 * '<S45>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/I Gain/External Parameters'
 * '<S46>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Ideal P Gain/Passthrough'
 * '<S47>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S48>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Integrator/Continuous'
 * '<S49>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Integrator ICs/Internal IC'
 * '<S50>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/N Copy/Disabled'
 * '<S51>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/N Gain/External Parameters'
 * '<S52>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/P Copy/Disabled'
 * '<S53>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Parallel P Gain/External Parameters'
 * '<S54>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Reset Signal/Disabled'
 * '<S55>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Saturation/Passthrough'
 * '<S56>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Saturation Fdbk/Disabled'
 * '<S57>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Sum/Sum_PID'
 * '<S58>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Sum Fdbk/Disabled'
 * '<S59>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Tracking Mode/Disabled'
 * '<S60>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Tracking Mode Sum/Passthrough'
 * '<S61>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Tsamp - Integral/Passthrough'
 * '<S62>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/Tsamp - Ngain/Passthrough'
 * '<S63>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/postSat Signal/Forward_Path'
 * '<S64>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller/preSat Signal/Forward_Path'
 * '<S65>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Anti-windup'
 * '<S66>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/D Gain'
 * '<S67>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Filter'
 * '<S68>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Filter ICs'
 * '<S69>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/I Gain'
 * '<S70>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Ideal P Gain'
 * '<S71>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Ideal P Gain Fdbk'
 * '<S72>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Integrator'
 * '<S73>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Integrator ICs'
 * '<S74>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/N Copy'
 * '<S75>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/N Gain'
 * '<S76>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/P Copy'
 * '<S77>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Parallel P Gain'
 * '<S78>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Reset Signal'
 * '<S79>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Saturation'
 * '<S80>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Saturation Fdbk'
 * '<S81>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Sum'
 * '<S82>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Sum Fdbk'
 * '<S83>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Tracking Mode'
 * '<S84>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Tracking Mode Sum'
 * '<S85>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Tsamp - Integral'
 * '<S86>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Tsamp - Ngain'
 * '<S87>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/postSat Signal'
 * '<S88>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/preSat Signal'
 * '<S89>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Anti-windup/Passthrough'
 * '<S90>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/D Gain/Internal Parameters'
 * '<S91>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Filter/Cont. Filter'
 * '<S92>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Filter ICs/Internal IC - Filter'
 * '<S93>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/I Gain/Internal Parameters'
 * '<S94>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Ideal P Gain/Passthrough'
 * '<S95>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Ideal P Gain Fdbk/Disabled'
 * '<S96>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Integrator/Continuous'
 * '<S97>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Integrator ICs/Internal IC'
 * '<S98>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/N Copy/Disabled'
 * '<S99>'  : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/N Gain/Internal Parameters'
 * '<S100>' : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/P Copy/Disabled'
 * '<S101>' : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Parallel P Gain/Internal Parameters'
 * '<S102>' : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Reset Signal/Disabled'
 * '<S103>' : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Saturation/Passthrough'
 * '<S104>' : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Saturation Fdbk/Disabled'
 * '<S105>' : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Sum/Sum_PID'
 * '<S106>' : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Sum Fdbk/Disabled'
 * '<S107>' : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Tracking Mode/Disabled'
 * '<S108>' : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Tracking Mode Sum/Passthrough'
 * '<S109>' : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Tsamp - Integral/Passthrough'
 * '<S110>' : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/Tsamp - Ngain/Passthrough'
 * '<S111>' : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/postSat Signal/Forward_Path'
 * '<S112>' : 'MalaksModelSub/VandD/Driver Subsystem/Driver Throttle Control/PID Controller1/preSat Signal/Forward_Path'
 * '<S113>' : 'MalaksModelSub/VandD/Orientation Gyro/Allan Noise'
 * '<S114>' : 'MalaksModelSub/VandD/Orientation Gyro/Allan Noise/MATLAB Function'
 * '<S115>' : 'MalaksModelSub/VandD/Orientation Gyro/Allan Noise/MATLAB Function1'
 * '<S116>' : 'MalaksModelSub/VandD/Vehicle model/MATLAB Function'
 * '<S117>' : 'MalaksModelSub/VandD/Vehicle model/MATLAB Function1'
 * '<S118>' : 'MalaksModelSub/VandD/Vehicle model/MATLAB Function2'
 * '<S119>' : 'MalaksModelSub/VandD/Vehicle model/MATLAB Function3'
 * '<S120>' : 'MalaksModelSub/VandD/Vehicle model/MATLAB Function4'
 * '<S121>' : 'MalaksModelSub/VandD/Vehicle model/MATLAB Function5'
 */
#endif                                 /* RTW_HEADER_VandD_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
