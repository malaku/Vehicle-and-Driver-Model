/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: VandD.c
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

#include "VandD.h"
#include "rtwtypes.h"
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <net/if.h>
#include <unistd.h> // Add this line for write() and close()
#include <linux/can.h>
#include <linux/can/raw.h>
#include <arpa/inet.h>

#define NumBitsPerChar                 8U
#define CAN_INTERFACE "can1"

int openCANSocket() {
    int p;
    struct sockaddr_can addr;
    struct ifreq ifr;

    p = socket(PF_CAN, SOCK_RAW, CAN_RAW);
    if (p == -1) {
        perror("socket");
        exit(EXIT_FAILURE);
    }

    strcpy(ifr.ifr_name, CAN_INTERFACE);
    ioctl(p, SIOCGIFINDEX, &ifr);

    addr.can_family = AF_CAN;
    addr.can_ifindex = ifr.ifr_ifindex;

    if (bind(p, (struct sockaddr *)&addr, sizeof(addr)) < 0) {
        perror("bind");
        exit(EXIT_FAILURE);
    }

    return p;
}

void closeCANSocket(int socketDescriptor) {
    close(socketDescriptor);
}

void sendCANMessage(int socketDescriptor, uint32_t messageID, float value) {
    struct can_frame frame;

    frame.can_id = messageID;
    frame.can_id &= CAN_SFF_MASK; // Ensure it's a standard identifier
    frame.can_dlc = sizeof(float);

    // Serialize the float value into network byte order (big-endian)
    memcpy(frame.data, &value, sizeof(float));
    //printf("Value sent: %f\n", value);

    if (write(socketDescriptor, &frame, sizeof(struct can_frame)) != sizeof(struct can_frame)) {
        perror("Error sending CAN message");
    }
}

void receiveCANMessage(int socketDescriptor, int targetID, float *tempvar) {
    struct can_frame frame;

    //printf("Listening for CAN messages...\n");

    while (1) {
        ssize_t nbytes = read(socketDescriptor, &frame, sizeof(struct can_frame));
        if (nbytes < 0) {
            perror("read");
            break;
        }

        if ((frame.can_id & CAN_ERR_FLAG) == 0 && (frame.can_id & CAN_RTR_FLAG) == 0 &&
            (frame.can_id & CAN_SFF_MASK) == targetID && frame.can_dlc == sizeof(float)) {
            // Deserialize the float value from network byte order and copy it directly into *tempvar
            memcpy(tempvar, frame.data, sizeof(float));
            //printf("ID: %x value: %f\n", frame.can_id, *tempvar);
            break;
        }
    }
}
void receiveCAN() {
    int socketDescriptor = openCANSocket();

    int targetMessageID1 = 0x101; // Desired FW1
    int targetMessageID2 = 0x201; // Desired FW1
    int targetMessageID3 = 0x301; // Desired FW1
    int targetMessageID4 = 0x401; // Desired FW1

    uint8_t receivedData[4]; // Buffer to store received data
    float lamda[4];

    receiveCANMessage(socketDescriptor, targetMessageID1, &lamda[0]);
    receiveCANMessage(socketDescriptor, targetMessageID2, &lamda[1]);
    receiveCANMessage(socketDescriptor, targetMessageID3, &lamda[2]);
    receiveCANMessage(socketDescriptor, targetMessageID4, &lamda[3]);

    rtU.lamda[0] = lamda[0];
    rtU.lamda[1] = lamda[1];
    rtU.lamda[2] = lamda[2];
    rtU.lamda[3] = lamda[3];

    closeCANSocket(socketDescriptor);
}

/* Private macros used by the generated code to access rtModel */
#ifndef rtmIsMajorTimeStep
#define rtmIsMajorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
#define rtmIsMinorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
#define rtmSetTPtr(rtm, val)           ((rtm)->Timing.t = (val))
#endif

/* Continuous states */
X rtX;

/* Block signals and states (default storage) */
DW rtDW;

/* External inputs (root inport signals with default storage) */
ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
ExtY rtY;

/* Real-time model */
static RT_MODEL rtM_;
RT_MODEL *const rtM = &rtM_;
extern real_T rt_powd_snf(real_T u0, real_T u1);
static void VandD_Init(DW_VandD *localDW, P_VandD *localP, X_VandD *localX);
static void VandD_Deriv(DW_VandD *localDW, P_VandD *localP, X_VandD *localX,
  XDot_VandD *localXdot);
static void VandD_Update(RT_MODEL * const rtM, DW_VandD *localDW);
static void VandD_g(RT_MODEL * const rtM, const real_T rtu_lamda[4], real_T
                    rtyyyy_Fy1Fy2Fy3Fy4[4], real_T *rty_Steeringwheelangle,
                    DW_VandD *localDW, P_VandD *localP, X_VandD *localX);

/* private model entry point functions */
extern void VandD_derivatives(void);

/* Forward declaration for local functions */
static boolean_T anyNonFinite(const real_T x[16]);
static real_T xnrm2(int32_T n, const real_T x[16], int32_T ix0);
static real_T rt_hypotd_snf_m(real_T u0, real_T u1, DW_VandD *localDW);
static void xzlarf(int32_T m, int32_T n, int32_T iv0, real_T tau, real_T C[16],
                   int32_T ic0, real_T work[4], DW_VandD *localDW);
static void xgehrd(real_T a[16], real_T tau[3], DW_VandD *localDW);
static void xdlanv2(real_T *a, real_T *b, real_T *c, real_T *d, real_T *rt1r,
                    real_T *rt1i, real_T *rt2r, real_T *rt2i, real_T *cs, real_T
                    *sn, DW_VandD *localDW);
static void xrot(int32_T n, real_T x[16], int32_T ix0, int32_T iy0, real_T c,
                 real_T s, DW_VandD *localDW);
static void xrot_g(int32_T n, real_T x[16], int32_T ix0, int32_T iy0, real_T c,
                   real_T s, DW_VandD *localDW);
static void xrot_gn(real_T x[16], int32_T ix0, int32_T iy0, real_T c, real_T s,
                    DW_VandD *localDW);
static real_T xnrm2_n(int32_T n, const real_T x[3], DW_VandD *localDW);
static real_T xzlarfg(int32_T n, real_T *alpha1, real_T x[3], DW_VandD *localDW);
static int32_T eml_dlahqr(real_T h[16], real_T z[16], DW_VandD *localDW);
static void schur(const real_T A[16], real_T V[16], real_T T[16], DW_VandD
                  *localDW);
static real_T norm(const real_T x[16]);
static void xzgetrf(real_T A[16], int32_T ipiv[4], int32_T *info, DW_VandD
                    *localDW);
static void xtrsm(const real_T A[16], real_T B_0[16]);
static void inv(const real_T x[16], real_T y[16], DW_VandD *localDW);
static void mpower(const real_T a[16], real_T b, real_T c[16], DW_VandD *localDW);
static real_T log2_p(real_T x);
static void padeApproximation(const real_T A[16], const real_T A2[16], const
  real_T A4[16], const real_T A6[16], int32_T m, real_T F[16], DW_VandD *localDW);
static void recomputeBlockDiag(const real_T A[16], real_T F[16], const int32_T
  blockFormat[3]);
static void expm(real_T A[16], real_T F[16], DW_VandD *localDW);
static real_T eml_rand_shr3cong(uint32_T state[2]);
static void genrandu(uint32_T s, uint32_T *state, real_T *r);
static void genrand_uint32_vector(uint32_T mt[625], uint32_T u[2]);
static real_T genrandu_d(uint32_T mt[625]);
static real_T randn(DW_VandD *localDW);
static real_T rtGetInf(void);
static real32_T rtGetInfF(void);
static real_T rtGetMinusInf(void);
static real32_T rtGetMinusInfF(void);
static real_T rtGetNaN(void);
static real32_T rtGetNaNF(void);
extern real_T rtInf;
extern real_T rtMinusInf;
extern real_T rtNaN;
extern real32_T rtInfF;
extern real32_T rtMinusInfF;
extern real32_T rtNaNF;
static void rt_InitInfAndNaN(size_t realSize);
static boolean_T rtIsInf(real_T value);
static boolean_T rtIsInfF(real32_T value);
static boolean_T rtIsNaN(real_T value);
static boolean_T rtIsNaNF(real32_T value);
typedef struct {
  struct {
    uint32_T wordH;
    uint32_T wordL;
  } words;
} BigEndianIEEEDouble;

typedef struct {
  struct {
    uint32_T wordL;
    uint32_T wordH;
  } words;
} LittleEndianIEEEDouble;

typedef struct {
  union {
    real32_T wordLreal;
    uint32_T wordLuint;
  } wordL;
} IEEESingle;

real_T rtInf;
real_T rtMinusInf;
real_T rtNaN;
real32_T rtInfF;
real32_T rtMinusInfF;
real32_T rtNaNF;

/*
 * Initialize rtInf needed by the generated code.
 * Inf is initialized as non-signaling. Assumes IEEE.
 */
static real_T rtGetInf(void)
{
  size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
  real_T inf = 0.0;
  if (bitsPerReal == 32U) {
    inf = rtGetInfF();
  } else {
    union {
      LittleEndianIEEEDouble bitVal;
      real_T fltVal;
    } tmpVal;

    tmpVal.bitVal.words.wordH = 0x7FF00000U;
    tmpVal.bitVal.words.wordL = 0x00000000U;
    inf = tmpVal.fltVal;
  }

  return inf;
}

/*
 * Initialize rtInfF needed by the generated code.
 * Inf is initialized as non-signaling. Assumes IEEE.
 */
static real32_T rtGetInfF(void)
{
  IEEESingle infF;
  infF.wordL.wordLuint = 0x7F800000U;
  return infF.wordL.wordLreal;
}

/*
 * Initialize rtMinusInf needed by the generated code.
 * Inf is initialized as non-signaling. Assumes IEEE.
 */
static real_T rtGetMinusInf(void)
{
  size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
  real_T minf = 0.0;
  if (bitsPerReal == 32U) {
    minf = rtGetMinusInfF();
  } else {
    union {
      LittleEndianIEEEDouble bitVal;
      real_T fltVal;
    } tmpVal;

    tmpVal.bitVal.words.wordH = 0xFFF00000U;
    tmpVal.bitVal.words.wordL = 0x00000000U;
    minf = tmpVal.fltVal;
  }

  return minf;
}

/*
 * Initialize rtMinusInfF needed by the generated code.
 * Inf is initialized as non-signaling. Assumes IEEE.
 */
static real32_T rtGetMinusInfF(void)
{
  IEEESingle minfF;
  minfF.wordL.wordLuint = 0xFF800000U;
  return minfF.wordL.wordLreal;
}

/*
 * Initialize rtNaN needed by the generated code.
 * NaN is initialized as non-signaling. Assumes IEEE.
 */
static real_T rtGetNaN(void)
{
  size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
  real_T nan = 0.0;
  if (bitsPerReal == 32U) {
    nan = rtGetNaNF();
  } else {
    union {
      LittleEndianIEEEDouble bitVal;
      real_T fltVal;
    } tmpVal;

    tmpVal.bitVal.words.wordH = 0xFFF80000U;
    tmpVal.bitVal.words.wordL = 0x00000000U;
    nan = tmpVal.fltVal;
  }

  return nan;
}

/*
 * Initialize rtNaNF needed by the generated code.
 * NaN is initialized as non-signaling. Assumes IEEE.
 */
static real32_T rtGetNaNF(void)
{
  IEEESingle nanF = { { 0.0F } };

  nanF.wordL.wordLuint = 0xFFC00000U;
  return nanF.wordL.wordLreal;
}

/*
 * Initialize the rtInf, rtMinusInf, and rtNaN needed by the
 * generated code. NaN is initialized as non-signaling. Assumes IEEE.
 */
static void rt_InitInfAndNaN(size_t realSize)
{
  (void) (realSize);
  rtNaN = rtGetNaN();
  rtNaNF = rtGetNaNF();
  rtInf = rtGetInf();
  rtInfF = rtGetInfF();
  rtMinusInf = rtGetMinusInf();
  rtMinusInfF = rtGetMinusInfF();
}

/* Test if value is infinite */
static boolean_T rtIsInf(real_T value)
{
  return (boolean_T)((value==rtInf || value==rtMinusInf) ? 1U : 0U);
}

/* Test if single-precision value is infinite */
static boolean_T rtIsInfF(real32_T value)
{
  return (boolean_T)(((value)==rtInfF || (value)==rtMinusInfF) ? 1U : 0U);
}

/* Test if value is not a number */
static boolean_T rtIsNaN(real_T value)
{
  boolean_T result = (boolean_T) 0;
  size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
  if (bitsPerReal == 32U) {
    result = rtIsNaNF((real32_T)value);
  } else {
    union {
      LittleEndianIEEEDouble bitVal;
      real_T fltVal;
    } tmpVal;

    tmpVal.fltVal = value;
    result = (boolean_T)((tmpVal.bitVal.words.wordH & 0x7FF00000) == 0x7FF00000 &&
                         ( (tmpVal.bitVal.words.wordH & 0x000FFFFF) != 0 ||
                          (tmpVal.bitVal.words.wordL != 0) ));
  }

  return result;
}

/* Test if single-precision value is not a number */
static boolean_T rtIsNaNF(real32_T value)
{
  IEEESingle tmp;
  tmp.wordL.wordLreal = value;
  return (boolean_T)( (tmp.wordL.wordLuint & 0x7F800000) == 0x7F800000 &&
                     (tmp.wordL.wordLuint & 0x007FFFFF) != 0 );
}

/*
 * This function updates continuous states using the ODE4 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE4_IntgData *id = (ODE4_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T *f3 = id->f[3];
  real_T temp;
  int_T i;
  int_T nXc = 29;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  VandD_derivatives();

  /* f1 = f(t + (h/2), y + (h/2)*f0) */
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  VandD_step();
  VandD_derivatives();

  /* f2 = f(t + (h/2), y + (h/2)*f1) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  VandD_step();
  VandD_derivatives();

  /* f3 = f(t + h, y + h*f2) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  VandD_step();
  VandD_derivatives();

  /* tnew = t + h
     ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3) */
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i]);
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static boolean_T anyNonFinite(const real_T x[16])
{
  real_T x_0;
  int32_T k;
  boolean_T b_p;
  b_p = true;
  for (k = 0; k < 16; k++) {
    x_0 = x[k];
    if (b_p && (rtIsInf(x_0) || rtIsNaN(x_0))) {
      b_p = false;
    }
  }

  return !b_p;
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static real_T xnrm2(int32_T n, const real_T x[16], int32_T ix0)
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  int32_T k;
  int32_T kend;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

static real_T rt_hypotd_snf_m(real_T u0, real_T u1, DW_VandD *localDW)
{
  real_T y;
  localDW->a = fabs(u0);
  y = fabs(u1);
  if (localDW->a < y) {
    localDW->a /= y;
    y *= sqrt(localDW->a * localDW->a + 1.0);
  } else if (localDW->a > y) {
    y /= localDW->a;
    y = sqrt(y * y + 1.0) * localDW->a;
  } else if (!rtIsNaN(y)) {
    y = localDW->a * 1.4142135623730951;
  }

  return y;
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void xzlarf(int32_T m, int32_T n, int32_T iv0, real_T tau, real_T C[16],
                   int32_T ic0, real_T work[4], DW_VandD *localDW)
{
  int32_T coltop;
  int32_T d;
  int32_T exitg1;
  int32_T ia;
  int32_T iac;
  int32_T ix;
  int32_T jy;
  int32_T lastc;
  int32_T lastv;
  boolean_T exitg2;
  if (tau != 0.0) {
    lastv = m;
    lastc = iv0 + m;
    while ((lastv > 0) && (C[lastc - 2] == 0.0)) {
      lastv--;
      lastc--;
    }

    lastc = n - 1;
    exitg2 = false;
    while ((!exitg2) && (lastc + 1 > 0)) {
      coltop = (lastc << 2) + ic0;
      jy = coltop;
      do {
        exitg1 = 0;
        if (jy <= (coltop + lastv) - 1) {
          if (C[jy - 1] != 0.0) {
            exitg1 = 1;
          } else {
            jy++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = -1;
  }

  if (lastv > 0) {
    if (lastc + 1 != 0) {
      for (coltop = 0; coltop <= lastc; coltop++) {
        work[coltop] = 0.0;
      }

      coltop = 0;
      jy = (lastc << 2) + ic0;
      for (iac = ic0; iac <= jy; iac += 4) {
        ix = iv0;
        localDW->c_d = 0.0;
        d = (iac + lastv) - 1;
        for (ia = iac; ia <= d; ia++) {
          localDW->c_d += C[ia - 1] * C[ix - 1];
          ix++;
        }

        work[coltop] += localDW->c_d;
        coltop++;
      }
    }

    if (!(-tau == 0.0)) {
      coltop = ic0 - 1;
      jy = 0;
      for (iac = 0; iac <= lastc; iac++) {
        if (work[jy] != 0.0) {
          localDW->c_d = work[jy] * -tau;
          ix = iv0;
          d = coltop;
          ia = lastv + coltop;
          while (d + 1 <= ia) {
            C[d] += C[ix - 1] * localDW->c_d;
            ix++;
            d++;
          }
        }

        jy++;
        coltop += 4;
      }
    }
  }
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void xgehrd(real_T a[16], real_T tau[3], DW_VandD *localDW)
{
  int32_T b_ia;
  int32_T b_ix;
  int32_T exitg1;
  int32_T ia;
  int32_T jy;
  int32_T rowleft;
  boolean_T exitg2;
  localDW->work_g[0] = 0.0;
  localDW->work_g[1] = 0.0;
  localDW->work_g[2] = 0.0;
  localDW->work_g[3] = 0.0;
  localDW->alpha1 = a[1];
  localDW->lastc = 3;
  tau[0] = 0.0;
  localDW->xnorm = xnrm2(2, a, 3);
  if (localDW->xnorm != 0.0) {
    localDW->xnorm = rt_hypotd_snf_m(a[1], localDW->xnorm, localDW);
    if (a[1] >= 0.0) {
      localDW->xnorm = -localDW->xnorm;
    }

    if (fabs(localDW->xnorm) < 1.0020841800044864E-292) {
      localDW->knt = 0;
      do {
        localDW->knt++;
        for (rowleft = 3; rowleft < 5; rowleft++) {
          a[rowleft - 1] *= 9.9792015476736E+291;
        }

        localDW->xnorm *= 9.9792015476736E+291;
        localDW->alpha1 *= 9.9792015476736E+291;
      } while ((fabs(localDW->xnorm) < 1.0020841800044864E-292) && (localDW->knt
                < 20));

      localDW->xnorm = rt_hypotd_snf_m(localDW->alpha1, xnrm2(2, a, 3), localDW);
      if (localDW->alpha1 >= 0.0) {
        localDW->xnorm = -localDW->xnorm;
      }

      tau[0] = (localDW->xnorm - localDW->alpha1) / localDW->xnorm;
      localDW->alpha1 = 1.0 / (localDW->alpha1 - localDW->xnorm);
      while (localDW->lastc <= 4) {
        a[localDW->lastc - 1] *= localDW->alpha1;
        localDW->lastc++;
      }

      localDW->lastc = 0;
      while (localDW->lastc <= localDW->knt - 1) {
        localDW->xnorm *= 1.0020841800044864E-292;
        localDW->lastc++;
      }

      localDW->alpha1 = localDW->xnorm;
    } else {
      tau[0] = (localDW->xnorm - a[1]) / localDW->xnorm;
      localDW->alpha1 = 1.0 / (a[1] - localDW->xnorm);
      while (localDW->lastc <= 4) {
        a[localDW->lastc - 1] *= localDW->alpha1;
        localDW->lastc++;
      }

      localDW->alpha1 = localDW->xnorm;
    }
  }

  a[1] = 1.0;
  jy = 1;
  if (tau[0] != 0.0) {
    localDW->knt = 2;
    localDW->lastc = 3;
    while ((localDW->knt + 1 > 0) && (a[localDW->lastc] == 0.0)) {
      localDW->knt--;
      localDW->lastc--;
    }

    localDW->lastc = 4;
    exitg2 = false;
    while ((!exitg2) && (localDW->lastc > 0)) {
      b_ia = localDW->lastc + 4;
      do {
        exitg1 = 0;
        if (b_ia <= ((localDW->knt << 2) + localDW->lastc) + 4) {
          if (a[b_ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            b_ia += 4;
          }
        } else {
          localDW->lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    localDW->knt = -1;
    localDW->lastc = 0;
  }

  if (localDW->knt + 1 > 0) {
    if (localDW->lastc != 0) {
      for (rowleft = 0; rowleft < localDW->lastc; rowleft++) {
        localDW->work_g[rowleft] = 0.0;
      }

      rowleft = 1;
      b_ia = (localDW->knt << 2) + 5;
      for (b_ix = 5; b_ix <= b_ia; b_ix += 4) {
        localDW->iy = 0;
        localDW->f_l = (b_ix + localDW->lastc) - 1;
        for (ia = b_ix; ia <= localDW->f_l; ia++) {
          localDW->work_g[localDW->iy] += a[ia - 1] * a[rowleft];
          localDW->iy++;
        }

        rowleft++;
      }
    }

    if (!(-tau[0] == 0.0)) {
      rowleft = 4;
      for (b_ia = 0; b_ia <= localDW->knt; b_ia++) {
        if (a[jy] != 0.0) {
          localDW->xnorm = a[jy] * -tau[0];
          b_ix = 0;
          localDW->iy = rowleft;
          localDW->f_l = localDW->lastc + rowleft;
          while (localDW->iy + 1 <= localDW->f_l) {
            a[localDW->iy] += localDW->work_g[b_ix] * localDW->xnorm;
            b_ix++;
            localDW->iy++;
          }
        }

        jy++;
        rowleft += 4;
      }
    }
  }

  xzlarf(3, 3, 2, tau[0], a, 6, localDW->work_g, localDW);
  a[1] = localDW->alpha1;
  localDW->alpha1 = a[6];
  tau[1] = 0.0;
  localDW->xnorm = xnrm2(1, a, 8);
  if (localDW->xnorm != 0.0) {
    localDW->xnorm = rt_hypotd_snf_m(a[6], localDW->xnorm, localDW);
    if (a[6] >= 0.0) {
      localDW->xnorm = -localDW->xnorm;
    }

    if (fabs(localDW->xnorm) < 1.0020841800044864E-292) {
      localDW->knt = 0;
      do {
        localDW->knt++;
        a[7] *= 9.9792015476736E+291;
        localDW->xnorm *= 9.9792015476736E+291;
        localDW->alpha1 *= 9.9792015476736E+291;
      } while ((fabs(localDW->xnorm) < 1.0020841800044864E-292) && (localDW->knt
                < 20));

      localDW->xnorm = rt_hypotd_snf_m(localDW->alpha1, xnrm2(1, a, 8), localDW);
      if (localDW->alpha1 >= 0.0) {
        localDW->xnorm = -localDW->xnorm;
      }

      tau[1] = (localDW->xnorm - localDW->alpha1) / localDW->xnorm;
      a[7] *= 1.0 / (localDW->alpha1 - localDW->xnorm);
      localDW->lastc = 0;
      while (localDW->lastc <= localDW->knt - 1) {
        localDW->xnorm *= 1.0020841800044864E-292;
        localDW->lastc++;
      }

      localDW->alpha1 = localDW->xnorm;
    } else {
      tau[1] = (localDW->xnorm - a[6]) / localDW->xnorm;
      a[7] *= 1.0 / (a[6] - localDW->xnorm);
      localDW->alpha1 = localDW->xnorm;
    }
  }

  a[6] = 1.0;
  jy = 6;
  if (tau[1] != 0.0) {
    localDW->knt = 1;
    localDW->lastc = 7;
    while ((localDW->knt + 1 > 0) && (a[localDW->lastc] == 0.0)) {
      localDW->knt--;
      localDW->lastc--;
    }

    localDW->lastc = 4;
    exitg2 = false;
    while ((!exitg2) && (localDW->lastc > 0)) {
      b_ia = localDW->lastc + 8;
      do {
        exitg1 = 0;
        if (b_ia <= ((localDW->knt << 2) + localDW->lastc) + 8) {
          if (a[b_ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            b_ia += 4;
          }
        } else {
          localDW->lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    localDW->knt = -1;
    localDW->lastc = 0;
  }

  if (localDW->knt + 1 > 0) {
    if (localDW->lastc != 0) {
      for (rowleft = 0; rowleft < localDW->lastc; rowleft++) {
        localDW->work_g[rowleft] = 0.0;
      }

      rowleft = 6;
      b_ia = (localDW->knt << 2) + 9;
      for (b_ix = 9; b_ix <= b_ia; b_ix += 4) {
        localDW->iy = 0;
        localDW->f_l = (b_ix + localDW->lastc) - 1;
        for (ia = b_ix; ia <= localDW->f_l; ia++) {
          localDW->work_g[localDW->iy] += a[ia - 1] * a[rowleft];
          localDW->iy++;
        }

        rowleft++;
      }
    }

    if (!(-tau[1] == 0.0)) {
      rowleft = 8;
      for (b_ia = 0; b_ia <= localDW->knt; b_ia++) {
        if (a[jy] != 0.0) {
          localDW->xnorm = a[jy] * -tau[1];
          b_ix = 0;
          localDW->iy = rowleft;
          localDW->f_l = localDW->lastc + rowleft;
          while (localDW->iy + 1 <= localDW->f_l) {
            a[localDW->iy] += localDW->work_g[b_ix] * localDW->xnorm;
            b_ix++;
            localDW->iy++;
          }
        }

        jy++;
        rowleft += 4;
      }
    }
  }

  xzlarf(2, 2, 7, tau[1], a, 11, localDW->work_g, localDW);
  a[6] = localDW->alpha1;
  localDW->alpha1 = a[11];
  tau[2] = 0.0;
  localDW->xnorm = xnrm2(0, a, 12);
  if (localDW->xnorm != 0.0) {
    localDW->xnorm = rt_hypotd_snf_m(a[11], localDW->xnorm, localDW);
    if (a[11] >= 0.0) {
      localDW->xnorm = -localDW->xnorm;
    }

    if (fabs(localDW->xnorm) < 1.0020841800044864E-292) {
      localDW->knt = 0;
      do {
        localDW->knt++;
        localDW->xnorm *= 9.9792015476736E+291;
        localDW->alpha1 *= 9.9792015476736E+291;
      } while ((fabs(localDW->xnorm) < 1.0020841800044864E-292) && (localDW->knt
                < 20));

      localDW->xnorm = rt_hypotd_snf_m(localDW->alpha1, xnrm2(0, a, 12), localDW);
      if (localDW->alpha1 >= 0.0) {
        localDW->xnorm = -localDW->xnorm;
      }

      tau[2] = (localDW->xnorm - localDW->alpha1) / localDW->xnorm;
      localDW->lastc = 0;
      while (localDW->lastc <= localDW->knt - 1) {
        localDW->xnorm *= 1.0020841800044864E-292;
        localDW->lastc++;
      }

      localDW->alpha1 = localDW->xnorm;
    } else {
      tau[2] = (localDW->xnorm - a[11]) / localDW->xnorm;
      localDW->alpha1 = localDW->xnorm;
    }
  }

  a[11] = 1.0;
  jy = 11;
  if (tau[2] != 0.0) {
    localDW->knt = 0;
    localDW->lastc = 11;
    while ((localDW->knt + 1 > 0) && (a[localDW->lastc] == 0.0)) {
      localDW->knt--;
      localDW->lastc--;
    }

    localDW->lastc = 4;
    exitg2 = false;
    while ((!exitg2) && (localDW->lastc > 0)) {
      b_ia = localDW->lastc + 12;
      do {
        exitg1 = 0;
        if (b_ia <= ((localDW->knt << 2) + localDW->lastc) + 12) {
          if (a[b_ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            b_ia += 4;
          }
        } else {
          localDW->lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    localDW->knt = -1;
    localDW->lastc = 0;
  }

  if (localDW->knt + 1 > 0) {
    if (localDW->lastc != 0) {
      for (rowleft = 0; rowleft < localDW->lastc; rowleft++) {
        localDW->work_g[rowleft] = 0.0;
      }

      rowleft = 11;
      b_ia = (localDW->knt << 2) + 13;
      for (b_ix = 13; b_ix <= b_ia; b_ix += 4) {
        localDW->iy = 0;
        localDW->f_l = (b_ix + localDW->lastc) - 1;
        for (ia = b_ix; ia <= localDW->f_l; ia++) {
          localDW->work_g[localDW->iy] += a[ia - 1] * a[rowleft];
          localDW->iy++;
        }

        rowleft++;
      }
    }

    if (!(-tau[2] == 0.0)) {
      rowleft = 12;
      for (b_ia = 0; b_ia <= localDW->knt; b_ia++) {
        if (a[jy] != 0.0) {
          localDW->xnorm = a[jy] * -tau[2];
          b_ix = 0;
          localDW->iy = rowleft;
          localDW->f_l = localDW->lastc + rowleft;
          while (localDW->iy + 1 <= localDW->f_l) {
            a[localDW->iy] += localDW->work_g[b_ix] * localDW->xnorm;
            b_ix++;
            localDW->iy++;
          }
        }

        jy++;
        rowleft += 4;
      }
    }
  }

  xzlarf(1, 1, 12, tau[2], a, 16, localDW->work_g, localDW);
  a[11] = localDW->alpha1;
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void xdlanv2(real_T *a, real_T *b, real_T *c, real_T *d, real_T *rt1r,
                    real_T *rt1i, real_T *rt2r, real_T *rt2i, real_T *cs, real_T
                    *sn, DW_VandD *localDW)
{
  int32_T b_0;
  int32_T c_0;
  boolean_T tmp;
  if (*c == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (*b == 0.0) {
    *cs = 0.0;
    *sn = 1.0;
    localDW->temp = *d;
    *d = *a;
    *a = localDW->temp;
    *b = -*c;
    *c = 0.0;
  } else {
    localDW->temp = *a - *d;
    if ((localDW->temp == 0.0) && ((*b < 0.0) != (*c < 0.0))) {
      *cs = 1.0;
      *sn = 0.0;
    } else {
      localDW->p_j = 0.5 * localDW->temp;
      localDW->bcmis = fabs(*b);
      localDW->z = fabs(*c);
      tmp = rtIsNaN(localDW->z);
      if ((localDW->bcmis >= localDW->z) || tmp) {
        localDW->bcmax = localDW->bcmis;
      } else {
        localDW->bcmax = localDW->z;
      }

      if ((localDW->bcmis <= localDW->z) || tmp) {
        localDW->z = localDW->bcmis;
      }

      if (!(*b < 0.0)) {
        b_0 = 1;
      } else {
        b_0 = -1;
      }

      if (!(*c < 0.0)) {
        c_0 = 1;
      } else {
        c_0 = -1;
      }

      localDW->bcmis = localDW->z * (real_T)b_0 * (real_T)c_0;
      localDW->scale = fabs(localDW->p_j);
      if ((!(localDW->scale >= localDW->bcmax)) && (!rtIsNaN(localDW->bcmax))) {
        localDW->scale = localDW->bcmax;
      }

      localDW->z = localDW->p_j / localDW->scale * localDW->p_j + localDW->bcmax
        / localDW->scale * localDW->bcmis;
      if (localDW->z >= 8.8817841970012523E-16) {
        if (!(localDW->p_j < 0.0)) {
          localDW->z = sqrt(localDW->scale) * sqrt(localDW->z);
        } else {
          localDW->z = -(sqrt(localDW->scale) * sqrt(localDW->z));
        }

        localDW->z += localDW->p_j;
        *a = *d + localDW->z;
        *d -= localDW->bcmax / localDW->z * localDW->bcmis;
        localDW->bcmax = rt_hypotd_snf_m(*c, localDW->z, localDW);
        *cs = localDW->z / localDW->bcmax;
        *sn = *c / localDW->bcmax;
        *b -= *c;
        *c = 0.0;
      } else {
        localDW->bcmis = *b + *c;
        localDW->bcmax = rt_hypotd_snf_m(localDW->bcmis, localDW->temp, localDW);
        *cs = sqrt((fabs(localDW->bcmis) / localDW->bcmax + 1.0) * 0.5);
        if (!(localDW->bcmis < 0.0)) {
          b_0 = 1;
        } else {
          b_0 = -1;
        }

        *sn = -(localDW->p_j / (localDW->bcmax * *cs)) * (real_T)b_0;
        localDW->temp = *a * *cs + *b * *sn;
        localDW->p_j = -*a * *sn + *b * *cs;
        localDW->bcmax = *c * *cs + *d * *sn;
        localDW->bcmis = -*c * *sn + *d * *cs;
        *b = localDW->p_j * *cs + localDW->bcmis * *sn;
        *c = -localDW->temp * *sn + localDW->bcmax * *cs;
        localDW->temp = ((localDW->temp * *cs + localDW->bcmax * *sn) +
                         (-localDW->p_j * *sn + localDW->bcmis * *cs)) * 0.5;
        *a = localDW->temp;
        *d = localDW->temp;
        if (*c != 0.0) {
          if (*b != 0.0) {
            if ((*b < 0.0) == (*c < 0.0)) {
              localDW->z = sqrt(fabs(*b));
              localDW->bcmis = sqrt(fabs(*c));
              localDW->p_j = localDW->z * localDW->bcmis;
              if (*c < 0.0) {
                localDW->p_j = -localDW->p_j;
              }

              localDW->bcmax = 1.0 / sqrt(fabs(*b + *c));
              *a = localDW->temp + localDW->p_j;
              *d = localDW->temp - localDW->p_j;
              *b -= *c;
              *c = 0.0;
              localDW->p_j = localDW->z * localDW->bcmax;
              localDW->bcmax *= localDW->bcmis;
              localDW->temp = *cs * localDW->p_j - *sn * localDW->bcmax;
              *sn = *cs * localDW->bcmax + *sn * localDW->p_j;
              *cs = localDW->temp;
            }
          } else {
            *b = -*c;
            *c = 0.0;
            localDW->temp = *cs;
            *cs = -*sn;
            *sn = localDW->temp;
          }
        }
      }
    }
  }

  *rt1r = *a;
  *rt2r = *d;
  if (*c == 0.0) {
    *rt1i = 0.0;
    *rt2i = 0.0;
  } else {
    *rt1i = sqrt(fabs(*b)) * sqrt(fabs(*c));
    *rt2i = -*rt1i;
  }
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void xrot(int32_T n, real_T x[16], int32_T ix0, int32_T iy0, real_T c,
                 real_T s, DW_VandD *localDW)
{
  int32_T ix;
  int32_T iy;
  int32_T k;
  if (n >= 1) {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      localDW->temp_o = c * x[ix] + s * x[iy];
      x[iy] = c * x[iy] - s * x[ix];
      x[ix] = localDW->temp_o;
      iy += 4;
      ix += 4;
    }
  }
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void xrot_g(int32_T n, real_T x[16], int32_T ix0, int32_T iy0, real_T c,
                   real_T s, DW_VandD *localDW)
{
  int32_T ix;
  int32_T iy;
  int32_T k;
  if (n >= 1) {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      localDW->temp_l = c * x[ix] + s * x[iy];
      x[iy] = c * x[iy] - s * x[ix];
      x[ix] = localDW->temp_l;
      iy++;
      ix++;
    }
  }
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void xrot_gn(real_T x[16], int32_T ix0, int32_T iy0, real_T c, real_T s,
                    DW_VandD *localDW)
{
  real_T temp_tmp;
  localDW->temp_b = x[iy0 - 1];
  temp_tmp = x[ix0 - 1];
  x[iy0 - 1] = localDW->temp_b * c - temp_tmp * s;
  x[ix0 - 1] = temp_tmp * c + localDW->temp_b * s;
  localDW->temp_b = x[ix0] * c + x[iy0] * s;
  x[iy0] = x[iy0] * c - x[ix0] * s;
  x[ix0] = localDW->temp_b;
  localDW->temp_b = x[iy0 + 1];
  temp_tmp = x[ix0 + 1];
  x[iy0 + 1] = localDW->temp_b * c - temp_tmp * s;
  x[ix0 + 1] = temp_tmp * c + localDW->temp_b * s;
  localDW->temp_b = x[iy0 + 2];
  temp_tmp = x[ix0 + 2];
  x[iy0 + 2] = localDW->temp_b * c - temp_tmp * s;
  x[ix0 + 2] = temp_tmp * c + localDW->temp_b * s;
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static real_T xnrm2_n(int32_T n, const real_T x[3], DW_VandD *localDW)
{
  real_T y;
  int32_T k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[1]);
    } else {
      localDW->scale_d = 3.3121686421112381E-170;
      for (k = 2; k <= n + 1; k++) {
        localDW->absxk = fabs(x[k - 1]);
        if (localDW->absxk > localDW->scale_d) {
          localDW->t = localDW->scale_d / localDW->absxk;
          y = y * localDW->t * localDW->t + 1.0;
          localDW->scale_d = localDW->absxk;
        } else {
          localDW->t = localDW->absxk / localDW->scale_d;
          y += localDW->t * localDW->t;
        }
      }

      y = localDW->scale_d * sqrt(y);
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static real_T xzlarfg(int32_T n, real_T *alpha1, real_T x[3], DW_VandD *localDW)
{
  real_T tau;
  int32_T b_k;
  int32_T knt;
  tau = 0.0;
  if (n > 0) {
    localDW->xnorm_g = xnrm2_n(n - 1, x, localDW);
    if (localDW->xnorm_g != 0.0) {
      localDW->xnorm_g = rt_hypotd_snf_m(*alpha1, localDW->xnorm_g, localDW);
      if (*alpha1 >= 0.0) {
        localDW->xnorm_g = -localDW->xnorm_g;
      }

      if (fabs(localDW->xnorm_g) < 1.0020841800044864E-292) {
        knt = 0;
        do {
          knt++;
          for (b_k = 1; b_k < n; b_k++) {
            x[b_k] *= 9.9792015476736E+291;
          }

          localDW->xnorm_g *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while ((fabs(localDW->xnorm_g) < 1.0020841800044864E-292) && (knt < 20));

        localDW->xnorm_g = rt_hypotd_snf_m(*alpha1, xnrm2_n(n - 1, x, localDW),
          localDW);
        if (*alpha1 >= 0.0) {
          localDW->xnorm_g = -localDW->xnorm_g;
        }

        tau = (localDW->xnorm_g - *alpha1) / localDW->xnorm_g;
        localDW->a_l = 1.0 / (*alpha1 - localDW->xnorm_g);
        for (b_k = 1; b_k < n; b_k++) {
          x[b_k] *= localDW->a_l;
        }

        for (b_k = 0; b_k < knt; b_k++) {
          localDW->xnorm_g *= 1.0020841800044864E-292;
        }

        *alpha1 = localDW->xnorm_g;
      } else {
        tau = (localDW->xnorm_g - *alpha1) / localDW->xnorm_g;
        localDW->a_l = 1.0 / (*alpha1 - localDW->xnorm_g);
        for (knt = 1; knt < n; knt++) {
          x[knt] *= localDW->a_l;
        }

        *alpha1 = localDW->xnorm_g;
      }
    }
  }

  return tau;
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static int32_T eml_dlahqr(real_T h[16], real_T z[16], DW_VandD *localDW)
{
  int32_T L;
  int32_T b_k;
  int32_T h12_tmp;
  int32_T h12_tmp_0;
  int32_T hoffset;
  int32_T info;
  int32_T its;
  int32_T j;
  int32_T k;
  int32_T m;
  int32_T nr;
  boolean_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  boolean_T goto150;
  info = 0;
  localDW->v[0] = 0.0;
  localDW->v[1] = 0.0;
  localDW->v[2] = 0.0;
  h[2] = 0.0;
  h[3] = 0.0;
  h[7] = 0.0;
  localDW->i_h = 3;
  exitg1 = false;
  while ((!exitg1) && (localDW->i_h + 1 >= 1)) {
    L = 1;
    goto150 = false;
    its = 0;
    exitg2 = false;
    while ((!exitg2) && (its < 301)) {
      k = localDW->i_h;
      exitg3 = false;
      while ((!exitg3) && (k + 1 > L)) {
        hoffset = ((k - 1) << 2) + k;
        if (fabs(h[hoffset]) <= 4.0083367200179456E-292) {
          exitg3 = true;
        } else {
          nr = (k << 2) + k;
          localDW->tst = fabs(h[hoffset - 1]) + fabs(h[nr]);
          if (localDW->tst == 0.0) {
            if (k - 1 >= 1) {
              localDW->tst = fabs(h[(((k - 2) << 2) + k) - 1]);
            }

            if (k + 2 <= 4) {
              localDW->tst += fabs(h[nr + 1]);
            }
          }

          if (fabs(h[hoffset]) <= 2.2204460492503131E-16 * localDW->tst) {
            localDW->htmp1 = fabs(h[hoffset]);
            localDW->tst = fabs(h[nr - 1]);
            if (localDW->htmp1 > localDW->tst) {
              localDW->ab = localDW->htmp1;
              localDW->ba = localDW->tst;
            } else {
              localDW->ab = localDW->tst;
              localDW->ba = localDW->htmp1;
            }

            localDW->tst = h[nr];
            localDW->htmp1 = fabs(localDW->tst);
            localDW->tst = fabs(h[hoffset - 1] - localDW->tst);
            if (localDW->htmp1 > localDW->tst) {
              localDW->aa = localDW->htmp1;
              localDW->htmp1 = localDW->tst;
            } else {
              localDW->aa = localDW->tst;
            }

            localDW->tst = localDW->aa + localDW->ab;
            localDW->htmp1 = localDW->aa / localDW->tst * localDW->htmp1 *
              2.2204460492503131E-16;
            if ((localDW->htmp1 <= 4.0083367200179456E-292) || rtIsNaN
                (localDW->htmp1)) {
              localDW->htmp1 = 4.0083367200179456E-292;
            }

            if (localDW->ab / localDW->tst * localDW->ba <= localDW->htmp1) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }
      }

      L = k + 1;
      if (k + 1 > 1) {
        h[k + ((k - 1) << 2)] = 0.0;
      }

      if (k + 1 >= localDW->i_h) {
        goto150 = true;
        exitg2 = true;
      } else {
        switch (its) {
         case 10:
          nr = (k << 2) + k;
          localDW->tst = fabs(h[(((k + 1) << 2) + k) + 2]) + fabs(h[nr + 1]);
          localDW->ba = 0.75 * localDW->tst + h[nr];
          localDW->h12 = -0.4375 * localDW->tst;
          localDW->aa = localDW->tst;
          localDW->htmp1 = localDW->ba;
          break;

         case 20:
          localDW->tst = fabs(h[(((localDW->i_h - 2) << 2) + localDW->i_h) - 1])
            + fabs(h[((localDW->i_h - 1) << 2) + localDW->i_h]);
          localDW->ba = h[(localDW->i_h << 2) + localDW->i_h] + 0.75 *
            localDW->tst;
          localDW->h12 = -0.4375 * localDW->tst;
          localDW->aa = localDW->tst;
          localDW->htmp1 = localDW->ba;
          break;

         default:
          m = ((localDW->i_h - 1) << 2) + localDW->i_h;
          localDW->ba = h[m - 1];
          localDW->aa = h[m];
          localDW->h12 = h[((localDW->i_h << 2) + localDW->i_h) - 1];
          localDW->htmp1 = h[(localDW->i_h << 2) + localDW->i_h];
          break;
        }

        localDW->tst = ((fabs(localDW->ba) + fabs(localDW->h12)) + fabs
                        (localDW->aa)) + fabs(localDW->htmp1);
        if (localDW->tst == 0.0) {
          localDW->ba = 0.0;
          localDW->htmp1 = 0.0;
          localDW->ab = 0.0;
          localDW->aa = 0.0;
        } else {
          localDW->ba /= localDW->tst;
          localDW->htmp1 /= localDW->tst;
          localDW->ab = (localDW->ba + localDW->htmp1) / 2.0;
          localDW->ba = (localDW->ba - localDW->ab) * (localDW->htmp1 -
            localDW->ab) - localDW->h12 / localDW->tst * (localDW->aa /
            localDW->tst);
          localDW->aa = sqrt(fabs(localDW->ba));
          if (localDW->ba >= 0.0) {
            localDW->ba = localDW->ab * localDW->tst;
            localDW->ab = localDW->ba;
            localDW->htmp1 = localDW->aa * localDW->tst;
            localDW->aa = -localDW->htmp1;
          } else {
            localDW->ba = localDW->ab + localDW->aa;
            localDW->ab -= localDW->aa;
            if (fabs(localDW->ba - localDW->htmp1) <= fabs(localDW->ab -
                 localDW->htmp1)) {
              localDW->ba *= localDW->tst;
              localDW->ab = localDW->ba;
            } else {
              localDW->ab *= localDW->tst;
              localDW->ba = localDW->ab;
            }

            localDW->htmp1 = 0.0;
            localDW->aa = 0.0;
          }
        }

        m = localDW->i_h - 1;
        exitg3 = false;
        while ((!exitg3) && (m >= k + 1)) {
          nr = ((m - 1) << 2) + m;
          localDW->h12 = h[nr];
          localDW->tst_tmp_tmp = h[nr - 1];
          localDW->tst_tmp = localDW->tst_tmp_tmp - localDW->ab;
          localDW->tst = (fabs(localDW->tst_tmp) + fabs(localDW->aa)) + fabs
            (localDW->h12);
          localDW->h12 /= localDW->tst;
          nr = (m << 2) + m;
          localDW->v[0] = (localDW->tst_tmp / localDW->tst *
                           (localDW->tst_tmp_tmp - localDW->ba) + h[nr - 1] *
                           localDW->h12) - localDW->aa / localDW->tst *
            localDW->htmp1;
          localDW->tst_tmp = h[nr];
          localDW->v[1] = (((localDW->tst_tmp_tmp + localDW->tst_tmp) -
                            localDW->ba) - localDW->ab) * localDW->h12;
          localDW->v[2] = h[nr + 1] * localDW->h12;
          localDW->tst = (fabs(localDW->v[0]) + fabs(localDW->v[1])) + fabs
            (localDW->v[2]);
          localDW->v[0] /= localDW->tst;
          localDW->v[1] /= localDW->tst;
          localDW->v[2] /= localDW->tst;
          if (k + 1 == m) {
            exitg3 = true;
          } else {
            hoffset = ((m - 2) << 2) + m;
            if (fabs(h[hoffset - 1]) * (fabs(localDW->v[1]) + fabs(localDW->v[2]))
                <= ((fabs(h[hoffset - 2]) + fabs(localDW->tst_tmp_tmp)) + fabs
                    (localDW->tst_tmp)) * (2.2204460492503131E-16 * fabs
                 (localDW->v[0]))) {
              exitg3 = true;
            } else {
              m--;
            }
          }
        }

        for (b_k = m; b_k <= localDW->i_h; b_k++) {
          nr = (localDW->i_h - b_k) + 2;
          if (nr >= 3) {
            nr = 3;
          }

          if (b_k > m) {
            hoffset = ((b_k - 2) << 2) + b_k;
            for (j = 0; j < nr; j++) {
              localDW->v[j] = h[(j + hoffset) - 1];
            }
          }

          localDW->htmp1 = localDW->v[0];
          localDW->tst = xzlarfg(nr, &localDW->htmp1, localDW->v, localDW);
          localDW->v[0] = localDW->htmp1;
          if (b_k > m) {
            h[(b_k + ((b_k - 2) << 2)) - 1] = localDW->htmp1;
            h[b_k + ((b_k - 2) << 2)] = 0.0;
            if (b_k < localDW->i_h) {
              h[(b_k + ((b_k - 2) << 2)) + 1] = 0.0;
            }
          } else if (m > k + 1) {
            hoffset = (((b_k - 2) << 2) + b_k) - 1;
            h[hoffset] *= 1.0 - localDW->tst;
          }

          localDW->htmp1 = localDW->v[1];
          localDW->ab = localDW->tst * localDW->v[1];
          switch (nr) {
           case 3:
            localDW->ba = localDW->v[2];
            localDW->aa = localDW->tst * localDW->v[2];
            for (nr = b_k - 1; nr + 1 < 5; nr++) {
              h12_tmp = (nr << 2) + b_k;
              localDW->h12 = (h[h12_tmp - 1] + h[h12_tmp] * localDW->htmp1) +
                h[h12_tmp + 1] * localDW->ba;
              h[h12_tmp - 1] -= localDW->h12 * localDW->tst;
              h[h12_tmp] -= localDW->h12 * localDW->ab;
              h[h12_tmp + 1] -= localDW->h12 * localDW->aa;
            }

            if (b_k + 3 <= localDW->i_h + 1) {
              nr = b_k + 3;
            } else {
              nr = localDW->i_h + 1;
            }

            for (hoffset = 0; hoffset < nr; hoffset++) {
              h12_tmp = ((b_k - 1) << 2) + hoffset;
              j = (b_k << 2) + hoffset;
              h12_tmp_0 = ((b_k + 1) << 2) + hoffset;
              localDW->h12 = (h[j] * localDW->htmp1 + h[h12_tmp]) + h[h12_tmp_0]
                * localDW->ba;
              h[h12_tmp] -= localDW->h12 * localDW->tst;
              h[j] -= localDW->h12 * localDW->ab;
              h[h12_tmp_0] -= localDW->h12 * localDW->aa;
            }

            for (nr = 0; nr < 4; nr++) {
              h12_tmp = ((b_k - 1) << 2) + nr;
              j = (b_k << 2) + nr;
              h12_tmp_0 = ((b_k + 1) << 2) + nr;
              localDW->h12 = (z[j] * localDW->htmp1 + z[h12_tmp]) + z[h12_tmp_0]
                * localDW->ba;
              z[h12_tmp] -= localDW->h12 * localDW->tst;
              z[j] -= localDW->h12 * localDW->ab;
              z[h12_tmp_0] -= localDW->h12 * localDW->aa;
            }
            break;

           case 2:
            for (nr = b_k - 1; nr + 1 < 5; nr++) {
              h12_tmp = (nr << 2) + b_k;
              localDW->ba = h[h12_tmp - 1];
              localDW->h12 = h[h12_tmp] * localDW->htmp1 + localDW->ba;
              h[h12_tmp - 1] = localDW->ba - localDW->h12 * localDW->tst;
              h[h12_tmp] -= localDW->h12 * localDW->ab;
            }

            for (nr = 0; nr <= localDW->i_h; nr++) {
              h12_tmp = ((b_k - 1) << 2) + nr;
              j = (b_k << 2) + nr;
              localDW->h12 = h[j] * localDW->htmp1 + h[h12_tmp];
              h[h12_tmp] -= localDW->h12 * localDW->tst;
              h[j] -= localDW->h12 * localDW->ab;
            }

            for (nr = 0; nr < 4; nr++) {
              h12_tmp = (b_k - 1) << 2;
              j = b_k << 2;
              localDW->h12 = z[j + nr] * localDW->htmp1 + z[h12_tmp + nr];
              hoffset = h12_tmp + nr;
              z[hoffset] -= localDW->h12 * localDW->tst;
              hoffset = j + nr;
              z[hoffset] -= localDW->h12 * localDW->ab;
            }
            break;
          }
        }

        its++;
      }
    }

    if (!goto150) {
      info = localDW->i_h + 1;
      exitg1 = true;
    } else {
      if ((localDW->i_h + 1 != L) && (L == localDW->i_h)) {
        nr = (localDW->i_h << 2) + localDW->i_h;
        localDW->tst = h[nr - 1];
        its = ((localDW->i_h - 1) << 2) + localDW->i_h;
        localDW->htmp1 = h[its];
        localDW->ab = h[nr];
        xdlanv2(&h[(localDW->i_h + ((localDW->i_h - 1) << 2)) - 1],
                &localDW->tst, &localDW->htmp1, &localDW->ab, &localDW->ba,
                &localDW->aa, &localDW->h12, &localDW->tst_tmp_tmp,
                &localDW->tst_tmp, &localDW->sn, localDW);
        h[nr - 1] = localDW->tst;
        h[its] = localDW->htmp1;
        h[nr] = localDW->ab;
        if (localDW->i_h + 1 < 4) {
          xrot(3 - localDW->i_h, h, localDW->i_h + ((localDW->i_h + 1) << 2),
               (localDW->i_h + ((localDW->i_h + 1) << 2)) + 1, localDW->tst_tmp,
               localDW->sn, localDW);
        }

        xrot_g(localDW->i_h - 1, h, ((localDW->i_h - 1) << 2) + 1, (localDW->i_h
                << 2) + 1, localDW->tst_tmp, localDW->sn, localDW);
        xrot_gn(z, ((localDW->i_h - 1) << 2) + 1, (localDW->i_h << 2) + 1,
                localDW->tst_tmp, localDW->sn, localDW);
      }

      localDW->i_h = L - 2;
    }
  }

  return info;
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void schur(const real_T A[16], real_T V[16], real_T T[16], DW_VandD
                  *localDW)
{
  if (anyNonFinite(A)) {
    for (localDW->iaii = 0; localDW->iaii < 16; localDW->iaii++) {
      V[localDW->iaii] = (rtNaN);
    }

    localDW->iaii = 2;
    while (localDW->iaii < 5) {
      V[localDW->iaii - 1] = 0.0;
      localDW->iaii++;
    }

    localDW->iaii = 3;
    while (localDW->iaii < 5) {
      V[localDW->iaii + 3] = 0.0;
      localDW->iaii++;
    }

    V[11] = 0.0;
    for (localDW->iaii = 0; localDW->iaii < 16; localDW->iaii++) {
      T[localDW->iaii] = (rtNaN);
    }
  } else {
    memcpy(&T[0], &A[0], sizeof(real_T) << 4U);
    xgehrd(T, localDW->tau, localDW);
    memcpy(&V[0], &T[0], sizeof(real_T) << 4U);
    localDW->iaii = 0;
    while (localDW->iaii <= 2) {
      V[localDW->iaii + 12] = 0.0;
      localDW->iaii++;
    }

    localDW->iaii = 0;
    while (localDW->iaii <= 1) {
      V[localDW->iaii + 8] = 0.0;
      localDW->iaii++;
    }

    localDW->iaii = 1;
    while (localDW->iaii + 3 <= 4) {
      V[localDW->iaii + 10] = V[localDW->iaii + 6];
      localDW->iaii++;
    }

    V[4] = 0.0;
    localDW->iaii = 0;
    while (localDW->iaii + 3 <= 4) {
      V[localDW->iaii + 6] = V[localDW->iaii + 2];
      localDW->iaii++;
    }

    localDW->work[0] = 0.0;
    V[1] = 0.0;
    localDW->work[1] = 0.0;
    V[2] = 0.0;
    localDW->work[2] = 0.0;
    V[3] = 0.0;
    localDW->work[3] = 0.0;
    V[0] = 1.0;
    V[15] = 1.0 - localDW->tau[2];
    localDW->iaii = 0;
    while (localDW->iaii <= 1) {
      V[14 - localDW->iaii] = 0.0;
      localDW->iaii++;
    }

    V[10] = 1.0;
    xzlarf(2, 1, 11, localDW->tau[1], V, 15, localDW->work, localDW);
    localDW->iaii = 11;
    while (localDW->iaii + 1 <= 12) {
      V[localDW->iaii] *= -localDW->tau[1];
      localDW->iaii++;
    }

    V[10] = 1.0 - localDW->tau[1];
    V[9] = 0.0;
    V[5] = 1.0;
    xzlarf(3, 2, 6, localDW->tau[0], V, 10, localDW->work, localDW);
    localDW->iaii = 6;
    while (localDW->iaii + 1 <= 8) {
      V[localDW->iaii] *= -localDW->tau[0];
      localDW->iaii++;
    }

    V[5] = 1.0 - localDW->tau[0];
    eml_dlahqr(T, V, localDW);
    T[3] = 0.0;
  }
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static real_T norm(const real_T x[16])
{
  real_T s;
  real_T y;
  int32_T j;
  int32_T s_tmp;
  boolean_T exitg1;
  y = 0.0;
  j = 0;
  exitg1 = false;
  while ((!exitg1) && (j < 4)) {
    s_tmp = j << 2;
    s = ((fabs(x[s_tmp + 1]) + fabs(x[s_tmp])) + fabs(x[s_tmp + 2])) + fabs
      (x[s_tmp + 3]);
    if (rtIsNaN(s)) {
      y = (rtNaN);
      exitg1 = true;
    } else {
      if (s > y) {
        y = s;
      }

      j++;
    }
  }

  return y;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T tmp;
  real_T tmp_0;
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void xzgetrf(real_T A[16], int32_T ipiv[4], int32_T *info, DW_VandD
                    *localDW)
{
  int32_T b_ix;
  int32_T c;
  int32_T c_ix;
  int32_T ijA;
  int32_T ix;
  int32_T j;
  int32_T jA;
  int32_T jj;
  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  ipiv[3] = 4;
  *info = 0;
  for (j = 0; j < 3; j++) {
    jj = j * 5;
    jA = 0;
    ix = jj;
    localDW->smax = fabs(A[jj]);
    for (b_ix = 2; b_ix <= 4 - j; b_ix++) {
      ix++;
      localDW->s = fabs(A[ix]);
      if (localDW->s > localDW->smax) {
        jA = b_ix - 1;
        localDW->smax = localDW->s;
      }
    }

    if (A[jj + jA] != 0.0) {
      if (jA != 0) {
        jA += j;
        ipiv[j] = jA + 1;
        localDW->smax = A[j];
        A[j] = A[jA];
        A[jA] = localDW->smax;
        localDW->smax = A[j + 4];
        A[j + 4] = A[jA + 4];
        A[jA + 4] = localDW->smax;
        localDW->smax = A[j + 8];
        A[j + 8] = A[jA + 8];
        A[jA + 8] = localDW->smax;
        localDW->smax = A[j + 12];
        A[j + 12] = A[jA + 12];
        A[jA + 12] = localDW->smax;
      }

      jA = (jj - j) + 4;
      for (ix = jj + 1; ix < jA; ix++) {
        A[ix] /= A[jj];
      }
    } else {
      *info = j + 1;
    }

    jA = jj;
    ix = jj + 4;
    for (b_ix = 0; b_ix <= 2 - j; b_ix++) {
      if (A[ix] != 0.0) {
        localDW->smax = -A[ix];
        c_ix = jj + 1;
        ijA = jA + 5;
        c = (jA - j) + 8;
        while (ijA + 1 <= c) {
          A[ijA] += A[c_ix] * localDW->smax;
          c_ix++;
          ijA++;
        }
      }

      ix += 4;
      jA += 4;
    }
  }

  if ((*info == 0) && (!(A[15] != 0.0))) {
    *info = 4;
  }
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void xtrsm(const real_T A[16], real_T B_0[16])
{
  real_T tmp;
  int32_T i;
  int32_T j;
  int32_T jBcol;
  int32_T tmp_0;
  for (j = 0; j < 4; j++) {
    jBcol = j << 2;
    tmp = B_0[jBcol + 3];
    if (tmp != 0.0) {
      B_0[jBcol + 3] = tmp / A[15];
      for (i = 0; i < 3; i++) {
        tmp_0 = i + jBcol;
        B_0[tmp_0] -= B_0[jBcol + 3] * A[i + 12];
      }
    }

    tmp = B_0[jBcol + 2];
    if (tmp != 0.0) {
      B_0[jBcol + 2] = tmp / A[10];
      for (i = 0; i < 2; i++) {
        tmp_0 = i + jBcol;
        B_0[tmp_0] -= B_0[jBcol + 2] * A[i + 8];
      }
    }

    tmp = B_0[jBcol + 1];
    if (tmp != 0.0) {
      B_0[jBcol + 1] = tmp / A[5];
      B_0[jBcol] -= B_0[jBcol + 1] * A[4];
    }

    if (B_0[jBcol] != 0.0) {
      B_0[jBcol] /= A[0];
    }
  }
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void inv(const real_T x[16], real_T y[16], DW_VandD *localDW)
{
  int32_T info;
  int32_T pipk;
  int32_T y_tmp;
  int32_T y_tmp_0;
  for (info = 0; info < 16; info++) {
    y[info] = 0.0;
    localDW->b_x[info] = x[info];
  }

  xzgetrf(localDW->b_x, localDW->ipiv, &info, localDW);
  localDW->p_b[0] = 1;
  localDW->p_b[1] = 2;
  localDW->p_b[2] = 3;
  localDW->p_b[3] = 4;
  if (localDW->ipiv[0] > 1) {
    pipk = localDW->p_b[localDW->ipiv[0] - 1];
    localDW->p_b[localDW->ipiv[0] - 1] = 1;
    localDW->p_b[0] = (int8_T)pipk;
  }

  if (localDW->ipiv[1] > 2) {
    pipk = localDW->p_b[localDW->ipiv[1] - 1];
    localDW->p_b[localDW->ipiv[1] - 1] = localDW->p_b[1];
    localDW->p_b[1] = (int8_T)pipk;
  }

  if (localDW->ipiv[2] > 3) {
    pipk = localDW->p_b[localDW->ipiv[2] - 1];
    localDW->p_b[localDW->ipiv[2] - 1] = localDW->p_b[2];
    localDW->p_b[2] = (int8_T)pipk;
  }

  pipk = localDW->p_b[0] - 1;
  y[(localDW->p_b[0] - 1) << 2] = 1.0;
  localDW->j = 0;
  while (localDW->j + 1 < 5) {
    info = pipk << 2;
    y_tmp_0 = info + localDW->j;
    if (y[y_tmp_0] != 0.0) {
      localDW->i_b = localDW->j + 1;
      while (localDW->i_b + 1 < 5) {
        y_tmp = info + localDW->i_b;
        y[y_tmp] -= localDW->b_x[(localDW->j << 2) + localDW->i_b] * y[y_tmp_0];
        localDW->i_b++;
      }
    }

    localDW->j++;
  }

  pipk = localDW->p_b[1] - 1;
  y[((localDW->p_b[1] - 1) << 2) + 1] = 1.0;
  localDW->j = 1;
  while (localDW->j + 1 < 5) {
    info = pipk << 2;
    y_tmp_0 = info + localDW->j;
    if (y[y_tmp_0] != 0.0) {
      localDW->i_b = localDW->j + 1;
      while (localDW->i_b + 1 < 5) {
        y_tmp = info + localDW->i_b;
        y[y_tmp] -= localDW->b_x[(localDW->j << 2) + localDW->i_b] * y[y_tmp_0];
        localDW->i_b++;
      }
    }

    localDW->j++;
  }

  info = (localDW->p_b[2] - 1) << 2;
  y[info + 2] = 1.0;
  localDW->j = 2;
  while (localDW->j + 1 < 5) {
    y_tmp_0 = info + localDW->j;
    if (y[y_tmp_0] != 0.0) {
      localDW->i_b = localDW->j + 1;
      while (localDW->i_b + 1 < 5) {
        y_tmp = info + localDW->i_b;
        y[y_tmp] -= localDW->b_x[(localDW->j << 2) + localDW->i_b] * y[y_tmp_0];
        localDW->i_b++;
      }
    }

    localDW->j++;
  }

  y_tmp = (localDW->p_b[3] - 1) << 2;
  y[y_tmp + 3] = 1.0;
  localDW->j = 3;
  while (localDW->j + 1 < 5) {
    info = y_tmp + localDW->j;
    if (y[info] != 0.0) {
      localDW->i_b = localDW->j + 1;
      while (localDW->i_b + 1 < 5) {
        y_tmp_0 = y_tmp + localDW->i_b;
        y[y_tmp_0] -= localDW->b_x[(localDW->j << 2) + localDW->i_b] * y[info];
        localDW->i_b++;
      }
    }

    localDW->j++;
  }

  xtrsm(localDW->b_x, y);
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void mpower(const real_T a[16], real_T b, real_T c[16], DW_VandD *localDW)
{
  int32_T exitg1;
  boolean_T aBufferInUse;
  boolean_T firstmult;
  boolean_T lsb;
  if (floor(b) == b) {
    localDW->e_m = fabs(b);
    if (localDW->e_m <= 2.147483647E+9) {
      memcpy(&localDW->b_a[0], &a[0], sizeof(real_T) << 4U);
      localDW->n = (int32_T)localDW->e_m;
      localDW->b_n = localDW->n;
      localDW->nbitson = 0;
      localDW->nb = -2;
      while (localDW->b_n > 0) {
        localDW->nb++;
        if ((localDW->b_n & 1U) != 0U) {
          localDW->nbitson++;
        }

        localDW->b_n >>= 1;
      }

      if (localDW->e_m <= 2.0) {
        if (b == 2.0) {
          for (localDW->nbitson = 0; localDW->nbitson < 4; localDW->nbitson++) {
            for (localDW->i1 = 0; localDW->i1 < 4; localDW->i1++) {
              localDW->c_tmp_tmp = localDW->nbitson << 2;
              localDW->c_tmp = localDW->i1 + localDW->c_tmp_tmp;
              c[localDW->c_tmp] = 0.0;
              c[localDW->c_tmp] += a[localDW->c_tmp_tmp] * a[localDW->i1];
              c[localDW->c_tmp] += a[localDW->c_tmp_tmp + 1] * a[localDW->i1 + 4];
              c[localDW->c_tmp] += a[localDW->c_tmp_tmp + 2] * a[localDW->i1 + 8];
              c[localDW->c_tmp] += a[localDW->c_tmp_tmp + 3] * a[localDW->i1 +
                12];
            }
          }
        } else if (b == 1.0) {
          memcpy(&c[0], &a[0], sizeof(real_T) << 4U);
        } else if (b == -1.0) {
          inv(a, c, localDW);
        } else if (b == -2.0) {
          for (localDW->nbitson = 0; localDW->nbitson < 4; localDW->nbitson++) {
            for (localDW->i1 = 0; localDW->i1 < 4; localDW->i1++) {
              localDW->n = localDW->nbitson << 2;
              localDW->b_n = localDW->i1 + localDW->n;
              localDW->b_a[localDW->b_n] = 0.0;
              localDW->b_a[localDW->b_n] += a[localDW->n] * a[localDW->i1];
              localDW->b_a[localDW->b_n] += a[localDW->n + 1] * a[localDW->i1 +
                4];
              localDW->b_a[localDW->b_n] += a[localDW->n + 2] * a[localDW->i1 +
                8];
              localDW->b_a[localDW->b_n] += a[localDW->n + 3] * a[localDW->i1 +
                12];
            }
          }

          inv(localDW->b_a, c, localDW);
        } else {
          firstmult = false;
          for (localDW->n = 0; localDW->n < 16; localDW->n++) {
            if (!firstmult) {
              firstmult = rtIsNaN(a[localDW->n]);
            }
          }

          if (firstmult) {
            for (localDW->nbitson = 0; localDW->nbitson < 16; localDW->nbitson++)
            {
              c[localDW->nbitson] = (rtNaN);
            }
          } else {
            memset(&c[0], 0, sizeof(real_T) << 4U);
            c[0] = 1.0;
            c[5] = 1.0;
            c[10] = 1.0;
            c[15] = 1.0;
          }
        }
      } else {
        firstmult = true;
        aBufferInUse = false;
        lsb = ((localDW->nbitson & 1U) != 0U);
        lsb = ((lsb && (b < 0.0)) || ((!lsb) && (b >= 0.0)));
        while (localDW->b_n <= localDW->nb) {
          if ((localDW->n & 1U) != 0U) {
            if (firstmult) {
              firstmult = false;
              if (lsb) {
                if (aBufferInUse) {
                  memcpy(&localDW->cBuffer[0], &localDW->aBuffer[0], sizeof
                         (real_T) << 4U);
                } else {
                  memcpy(&localDW->cBuffer[0], &localDW->b_a[0], sizeof(real_T) <<
                         4U);
                }
              } else if (aBufferInUse) {
                memcpy(&c[0], &localDW->aBuffer[0], sizeof(real_T) << 4U);
              } else {
                memcpy(&c[0], &localDW->b_a[0], sizeof(real_T) << 4U);
              }
            } else {
              if (aBufferInUse) {
                if (lsb) {
                  for (localDW->nbitson = 0; localDW->nbitson < 4;
                       localDW->nbitson++) {
                    for (localDW->i1 = 0; localDW->i1 < 4; localDW->i1++) {
                      localDW->c_tmp_tmp = localDW->i1 << 2;
                      localDW->c_tmp = localDW->nbitson + localDW->c_tmp_tmp;
                      c[localDW->c_tmp] = 0.0;
                      c[localDW->c_tmp] += localDW->aBuffer[localDW->c_tmp_tmp] *
                        localDW->cBuffer[localDW->nbitson];
                      c[localDW->c_tmp] += localDW->aBuffer[localDW->c_tmp_tmp +
                        1] * localDW->cBuffer[localDW->nbitson + 4];
                      c[localDW->c_tmp] += localDW->aBuffer[localDW->c_tmp_tmp +
                        2] * localDW->cBuffer[localDW->nbitson + 8];
                      c[localDW->c_tmp] += localDW->aBuffer[localDW->c_tmp_tmp +
                        3] * localDW->cBuffer[localDW->nbitson + 12];
                    }
                  }
                } else {
                  for (localDW->nbitson = 0; localDW->nbitson < 4;
                       localDW->nbitson++) {
                    for (localDW->i1 = 0; localDW->i1 < 4; localDW->i1++) {
                      localDW->c_tmp_tmp = localDW->i1 << 2;
                      localDW->c_tmp = localDW->nbitson + localDW->c_tmp_tmp;
                      localDW->cBuffer[localDW->c_tmp] = 0.0;
                      localDW->cBuffer[localDW->c_tmp] += localDW->
                        aBuffer[localDW->c_tmp_tmp] * c[localDW->nbitson];
                      localDW->cBuffer[localDW->c_tmp] += localDW->
                        aBuffer[localDW->c_tmp_tmp + 1] * c[localDW->nbitson + 4];
                      localDW->cBuffer[localDW->c_tmp] += localDW->
                        aBuffer[localDW->c_tmp_tmp + 2] * c[localDW->nbitson + 8];
                      localDW->cBuffer[localDW->c_tmp] += localDW->
                        aBuffer[localDW->c_tmp_tmp + 3] * c[localDW->nbitson +
                        12];
                    }
                  }
                }
              } else if (lsb) {
                for (localDW->nbitson = 0; localDW->nbitson < 4;
                     localDW->nbitson++) {
                  for (localDW->i1 = 0; localDW->i1 < 4; localDW->i1++) {
                    localDW->c_tmp_tmp = localDW->i1 << 2;
                    localDW->c_tmp = localDW->nbitson + localDW->c_tmp_tmp;
                    c[localDW->c_tmp] = 0.0;
                    c[localDW->c_tmp] += localDW->b_a[localDW->c_tmp_tmp] *
                      localDW->cBuffer[localDW->nbitson];
                    c[localDW->c_tmp] += localDW->b_a[localDW->c_tmp_tmp + 1] *
                      localDW->cBuffer[localDW->nbitson + 4];
                    c[localDW->c_tmp] += localDW->b_a[localDW->c_tmp_tmp + 2] *
                      localDW->cBuffer[localDW->nbitson + 8];
                    c[localDW->c_tmp] += localDW->b_a[localDW->c_tmp_tmp + 3] *
                      localDW->cBuffer[localDW->nbitson + 12];
                  }
                }
              } else {
                for (localDW->nbitson = 0; localDW->nbitson < 4;
                     localDW->nbitson++) {
                  for (localDW->i1 = 0; localDW->i1 < 4; localDW->i1++) {
                    localDW->c_tmp_tmp = localDW->i1 << 2;
                    localDW->c_tmp = localDW->nbitson + localDW->c_tmp_tmp;
                    localDW->cBuffer[localDW->c_tmp] = 0.0;
                    localDW->cBuffer[localDW->c_tmp] += localDW->b_a
                      [localDW->c_tmp_tmp] * c[localDW->nbitson];
                    localDW->cBuffer[localDW->c_tmp] += localDW->b_a
                      [localDW->c_tmp_tmp + 1] * c[localDW->nbitson + 4];
                    localDW->cBuffer[localDW->c_tmp] += localDW->b_a
                      [localDW->c_tmp_tmp + 2] * c[localDW->nbitson + 8];
                    localDW->cBuffer[localDW->c_tmp] += localDW->b_a
                      [localDW->c_tmp_tmp + 3] * c[localDW->nbitson + 12];
                  }
                }
              }

              lsb = !lsb;
            }
          }

          localDW->n >>= 1;
          if (aBufferInUse) {
            for (localDW->nbitson = 0; localDW->nbitson < 4; localDW->nbitson++)
            {
              for (localDW->i1 = 0; localDW->i1 < 4; localDW->i1++) {
                localDW->c_tmp_tmp = localDW->i1 << 2;
                localDW->c_tmp = localDW->nbitson + localDW->c_tmp_tmp;
                localDW->b_a[localDW->c_tmp] = 0.0;
                localDW->b_a[localDW->c_tmp] += localDW->aBuffer
                  [localDW->c_tmp_tmp] * localDW->aBuffer[localDW->nbitson];
                localDW->b_a[localDW->c_tmp] += localDW->aBuffer
                  [localDW->c_tmp_tmp + 1] * localDW->aBuffer[localDW->nbitson +
                  4];
                localDW->b_a[localDW->c_tmp] += localDW->aBuffer
                  [localDW->c_tmp_tmp + 2] * localDW->aBuffer[localDW->nbitson +
                  8];
                localDW->b_a[localDW->c_tmp] += localDW->aBuffer
                  [localDW->c_tmp_tmp + 3] * localDW->aBuffer[localDW->nbitson +
                  12];
              }
            }
          } else {
            for (localDW->nbitson = 0; localDW->nbitson < 4; localDW->nbitson++)
            {
              for (localDW->i1 = 0; localDW->i1 < 4; localDW->i1++) {
                localDW->c_tmp_tmp = localDW->i1 << 2;
                localDW->c_tmp = localDW->nbitson + localDW->c_tmp_tmp;
                localDW->aBuffer[localDW->c_tmp] = 0.0;
                localDW->aBuffer[localDW->c_tmp] += localDW->b_a
                  [localDW->c_tmp_tmp] * localDW->b_a[localDW->nbitson];
                localDW->aBuffer[localDW->c_tmp] += localDW->b_a
                  [localDW->c_tmp_tmp + 1] * localDW->b_a[localDW->nbitson + 4];
                localDW->aBuffer[localDW->c_tmp] += localDW->b_a
                  [localDW->c_tmp_tmp + 2] * localDW->b_a[localDW->nbitson + 8];
                localDW->aBuffer[localDW->c_tmp] += localDW->b_a
                  [localDW->c_tmp_tmp + 3] * localDW->b_a[localDW->nbitson + 12];
              }
            }
          }

          aBufferInUse = !aBufferInUse;
          localDW->b_n++;
        }

        if (firstmult) {
          if (b < 0.0) {
            if (aBufferInUse) {
              inv(localDW->aBuffer, c, localDW);
            } else {
              inv(localDW->b_a, c, localDW);
            }
          } else if (aBufferInUse) {
            memcpy(&c[0], &localDW->aBuffer[0], sizeof(real_T) << 4U);
          } else {
            memcpy(&c[0], &localDW->b_a[0], sizeof(real_T) << 4U);
          }
        } else if (b < 0.0) {
          for (localDW->nbitson = 0; localDW->nbitson < 4; localDW->nbitson++) {
            for (localDW->i1 = 0; localDW->i1 < 4; localDW->i1++) {
              localDW->c_tmp_tmp = localDW->nbitson << 2;
              localDW->c_tmp = localDW->i1 + localDW->c_tmp_tmp;
              localDW->cBuffer[localDW->c_tmp] = 0.0;
              localDW->cBuffer_c[localDW->c_tmp] = 0.0;
              localDW->e_m = c[localDW->i1];
              localDW->cBuffer[localDW->c_tmp] += localDW->aBuffer
                [localDW->c_tmp_tmp] * localDW->e_m;
              localDW->cBuffer_c[localDW->c_tmp] += localDW->b_a
                [localDW->c_tmp_tmp] * localDW->e_m;
              localDW->e_m = c[localDW->i1 + 4];
              localDW->cBuffer[localDW->c_tmp] += localDW->aBuffer
                [localDW->c_tmp_tmp + 1] * localDW->e_m;
              localDW->cBuffer_c[localDW->c_tmp] += localDW->b_a
                [localDW->c_tmp_tmp + 1] * localDW->e_m;
              localDW->e_m = c[localDW->i1 + 8];
              localDW->cBuffer[localDW->c_tmp] += localDW->aBuffer
                [localDW->c_tmp_tmp + 2] * localDW->e_m;
              localDW->cBuffer_c[localDW->c_tmp] += localDW->b_a
                [localDW->c_tmp_tmp + 2] * localDW->e_m;
              localDW->e_m = c[localDW->i1 + 12];
              localDW->cBuffer[localDW->c_tmp] += localDW->aBuffer
                [localDW->c_tmp_tmp + 3] * localDW->e_m;
              localDW->cBuffer_c[localDW->c_tmp] += localDW->b_a
                [localDW->c_tmp_tmp + 3] * localDW->e_m;
            }
          }

          if (aBufferInUse) {
            memcpy(&localDW->b_a[0], &localDW->cBuffer[0], sizeof(real_T) << 4U);
          } else {
            memcpy(&localDW->b_a[0], &localDW->cBuffer_c[0], sizeof(real_T) <<
                   4U);
          }

          inv(localDW->b_a, c, localDW);
        } else {
          for (localDW->nbitson = 0; localDW->nbitson < 4; localDW->nbitson++) {
            for (localDW->i1 = 0; localDW->i1 < 4; localDW->i1++) {
              localDW->c_tmp_tmp = localDW->nbitson << 2;
              localDW->c_tmp = localDW->i1 + localDW->c_tmp_tmp;
              localDW->cBuffer_c[localDW->c_tmp] = 0.0;
              c[localDW->c_tmp] = 0.0;
              localDW->e_m = localDW->cBuffer[localDW->i1];
              localDW->cBuffer_c[localDW->c_tmp] += localDW->aBuffer
                [localDW->c_tmp_tmp] * localDW->e_m;
              c[localDW->c_tmp] += localDW->b_a[localDW->c_tmp_tmp] *
                localDW->e_m;
              localDW->e_m = localDW->cBuffer[localDW->i1 + 4];
              localDW->cBuffer_c[localDW->c_tmp] += localDW->aBuffer
                [localDW->c_tmp_tmp + 1] * localDW->e_m;
              c[localDW->c_tmp] += localDW->b_a[localDW->c_tmp_tmp + 1] *
                localDW->e_m;
              localDW->e_m = localDW->cBuffer[localDW->i1 + 8];
              localDW->cBuffer_c[localDW->c_tmp] += localDW->aBuffer
                [localDW->c_tmp_tmp + 2] * localDW->e_m;
              c[localDW->c_tmp] += localDW->b_a[localDW->c_tmp_tmp + 2] *
                localDW->e_m;
              localDW->e_m = localDW->cBuffer[localDW->i1 + 12];
              localDW->cBuffer_c[localDW->c_tmp] += localDW->aBuffer
                [localDW->c_tmp_tmp + 3] * localDW->e_m;
              c[localDW->c_tmp] += localDW->b_a[localDW->c_tmp_tmp + 3] *
                localDW->e_m;
            }
          }

          if (aBufferInUse) {
            memcpy(&c[0], &localDW->cBuffer_c[0], sizeof(real_T) << 4U);
          }
        }
      }
    } else {
      memcpy(&localDW->b_a[0], &a[0], sizeof(real_T) << 4U);
      if (!rtIsInf(b)) {
        firstmult = true;
        do {
          exitg1 = 0;
          localDW->ed2 = floor(localDW->e_m / 2.0);
          if (2.0 * localDW->ed2 != localDW->e_m) {
            if (firstmult) {
              memcpy(&c[0], &localDW->b_a[0], sizeof(real_T) << 4U);
              firstmult = false;
            } else {
              for (localDW->nbitson = 0; localDW->nbitson < 4; localDW->nbitson
                   ++) {
                for (localDW->i1 = 0; localDW->i1 < 4; localDW->i1++) {
                  localDW->c_tmp_tmp = localDW->i1 << 2;
                  localDW->c_tmp = localDW->nbitson + localDW->c_tmp_tmp;
                  localDW->cBuffer[localDW->c_tmp] = 0.0;
                  localDW->cBuffer[localDW->c_tmp] += localDW->b_a
                    [localDW->c_tmp_tmp] * c[localDW->nbitson];
                  localDW->cBuffer[localDW->c_tmp] += localDW->b_a
                    [localDW->c_tmp_tmp + 1] * c[localDW->nbitson + 4];
                  localDW->cBuffer[localDW->c_tmp] += localDW->b_a
                    [localDW->c_tmp_tmp + 2] * c[localDW->nbitson + 8];
                  localDW->cBuffer[localDW->c_tmp] += localDW->b_a
                    [localDW->c_tmp_tmp + 3] * c[localDW->nbitson + 12];
                }
              }

              memcpy(&c[0], &localDW->cBuffer[0], sizeof(real_T) << 4U);
            }
          }

          if (localDW->ed2 == 0.0) {
            exitg1 = 1;
          } else {
            localDW->e_m = localDW->ed2;
            for (localDW->nbitson = 0; localDW->nbitson < 4; localDW->nbitson++)
            {
              for (localDW->i1 = 0; localDW->i1 < 4; localDW->i1++) {
                localDW->c_tmp_tmp = localDW->i1 << 2;
                localDW->c_tmp = localDW->nbitson + localDW->c_tmp_tmp;
                localDW->aBuffer[localDW->c_tmp] = 0.0;
                localDW->aBuffer[localDW->c_tmp] += localDW->b_a
                  [localDW->c_tmp_tmp] * localDW->b_a[localDW->nbitson];
                localDW->aBuffer[localDW->c_tmp] += localDW->b_a
                  [localDW->c_tmp_tmp + 1] * localDW->b_a[localDW->nbitson + 4];
                localDW->aBuffer[localDW->c_tmp] += localDW->b_a
                  [localDW->c_tmp_tmp + 2] * localDW->b_a[localDW->nbitson + 8];
                localDW->aBuffer[localDW->c_tmp] += localDW->b_a
                  [localDW->c_tmp_tmp + 3] * localDW->b_a[localDW->nbitson + 12];
              }
            }

            memcpy(&localDW->b_a[0], &localDW->aBuffer[0], sizeof(real_T) << 4U);
          }
        } while (exitg1 == 0);

        if (b < 0.0) {
          memcpy(&localDW->b_a[0], &c[0], sizeof(real_T) << 4);
          inv(localDW->b_a, c, localDW);
        }
      } else {
        for (localDW->nbitson = 0; localDW->nbitson < 16; localDW->nbitson++) {
          c[localDW->nbitson] = (rtNaN);
        }
      }
    }
  } else {
    for (localDW->nbitson = 0; localDW->nbitson < 16; localDW->nbitson++) {
      c[localDW->nbitson] = (rtNaN);
    }
  }
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static real_T log2_p(real_T x)
{
  real_T f;
  real_T t;
  int32_T inte;
  if (x == 0.0) {
    f = (rtMinusInf);
  } else if (x < 0.0) {
    f = (rtNaN);
  } else if ((!rtIsInf(x)) && (!rtIsNaN(x))) {
    t = frexp(x, &inte);
    if (t == 0.5) {
      f = (real_T)inte - 1.0;
    } else if ((inte == 1) && (t < 0.75)) {
      f = log(2.0 * t) / 0.69314718055994529;
    } else {
      f = log(t) / 0.69314718055994529 + (real_T)inte;
    }
  } else {
    f = x;
  }

  return f;
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void padeApproximation(const real_T A[16], const real_T A2[16], const
  real_T A4[16], const real_T A6[16], int32_T m, real_T F[16], DW_VandD *localDW)
{
  int32_T A6_tmp;
  int32_T A6_tmp_0;
  int32_T e_k;
  int32_T ip;
  switch (m) {
   case 3:
    memcpy(&F[0], &A2[0], sizeof(real_T) << 4U);
    F[0] += 60.0;
    F[5] += 60.0;
    F[10] += 60.0;
    F[15] += 60.0;
    for (e_k = 0; e_k < 4; e_k++) {
      for (ip = 0; ip < 4; ip++) {
        A6_tmp = e_k << 2;
        A6_tmp_0 = ip + A6_tmp;
        localDW->A6_g[A6_tmp_0] = 0.0;
        localDW->A6_g[A6_tmp_0] += F[A6_tmp] * A[ip];
        localDW->A6_g[A6_tmp_0] += F[A6_tmp + 1] * A[ip + 4];
        localDW->A6_g[A6_tmp_0] += F[A6_tmp + 2] * A[ip + 8];
        localDW->A6_g[A6_tmp_0] += F[A6_tmp + 3] * A[ip + 12];
      }
    }

    for (e_k = 0; e_k < 16; e_k++) {
      F[e_k] = localDW->A6_g[e_k];
      localDW->V[e_k] = 12.0 * A2[e_k];
    }

    localDW->d_n = 120.0;
    break;

   case 5:
    for (e_k = 0; e_k < 16; e_k++) {
      F[e_k] = 420.0 * A2[e_k] + A4[e_k];
    }

    F[0] += 15120.0;
    F[5] += 15120.0;
    F[10] += 15120.0;
    F[15] += 15120.0;
    for (e_k = 0; e_k < 4; e_k++) {
      for (ip = 0; ip < 4; ip++) {
        A6_tmp = e_k << 2;
        A6_tmp_0 = ip + A6_tmp;
        localDW->A6_g[A6_tmp_0] = 0.0;
        localDW->A6_g[A6_tmp_0] += F[A6_tmp] * A[ip];
        localDW->A6_g[A6_tmp_0] += F[A6_tmp + 1] * A[ip + 4];
        localDW->A6_g[A6_tmp_0] += F[A6_tmp + 2] * A[ip + 8];
        localDW->A6_g[A6_tmp_0] += F[A6_tmp + 3] * A[ip + 12];
      }
    }

    for (e_k = 0; e_k < 16; e_k++) {
      F[e_k] = localDW->A6_g[e_k];
      localDW->V[e_k] = 30.0 * A4[e_k] + 3360.0 * A2[e_k];
    }

    localDW->d_n = 30240.0;
    break;

   case 7:
    for (e_k = 0; e_k < 16; e_k++) {
      F[e_k] = (1512.0 * A4[e_k] + A6[e_k]) + 277200.0 * A2[e_k];
    }

    F[0] += 8.64864E+6;
    F[5] += 8.64864E+6;
    F[10] += 8.64864E+6;
    F[15] += 8.64864E+6;
    for (e_k = 0; e_k < 4; e_k++) {
      for (ip = 0; ip < 4; ip++) {
        A6_tmp = e_k << 2;
        A6_tmp_0 = ip + A6_tmp;
        localDW->A6_g[A6_tmp_0] = 0.0;
        localDW->A6_g[A6_tmp_0] += F[A6_tmp] * A[ip];
        localDW->A6_g[A6_tmp_0] += F[A6_tmp + 1] * A[ip + 4];
        localDW->A6_g[A6_tmp_0] += F[A6_tmp + 2] * A[ip + 8];
        localDW->A6_g[A6_tmp_0] += F[A6_tmp + 3] * A[ip + 12];
      }
    }

    for (e_k = 0; e_k < 16; e_k++) {
      F[e_k] = localDW->A6_g[e_k];
      localDW->V[e_k] = (56.0 * A6[e_k] + 25200.0 * A4[e_k]) + 1.99584E+6 *
        A2[e_k];
    }

    localDW->d_n = 1.729728E+7;
    break;

   case 9:
    for (e_k = 0; e_k < 4; e_k++) {
      for (ip = 0; ip < 4; ip++) {
        A6_tmp_0 = e_k << 2;
        A6_tmp = ip + A6_tmp_0;
        localDW->V[A6_tmp] = 0.0;
        localDW->V[A6_tmp] += A2[A6_tmp_0] * A6[ip];
        localDW->V[A6_tmp] += A2[A6_tmp_0 + 1] * A6[ip + 4];
        localDW->V[A6_tmp] += A2[A6_tmp_0 + 2] * A6[ip + 8];
        localDW->V[A6_tmp] += A2[A6_tmp_0 + 3] * A6[ip + 12];
      }
    }

    for (e_k = 0; e_k < 16; e_k++) {
      F[e_k] = ((3960.0 * A6[e_k] + localDW->V[e_k]) + 2.16216E+6 * A4[e_k]) +
        3.027024E+8 * A2[e_k];
    }

    F[0] += 8.8216128E+9;
    F[5] += 8.8216128E+9;
    F[10] += 8.8216128E+9;
    F[15] += 8.8216128E+9;
    for (e_k = 0; e_k < 4; e_k++) {
      for (ip = 0; ip < 4; ip++) {
        A6_tmp = e_k << 2;
        A6_tmp_0 = ip + A6_tmp;
        localDW->A6_g[A6_tmp_0] = 0.0;
        localDW->A6_g[A6_tmp_0] += F[A6_tmp] * A[ip];
        localDW->A6_g[A6_tmp_0] += F[A6_tmp + 1] * A[ip + 4];
        localDW->A6_g[A6_tmp_0] += F[A6_tmp + 2] * A[ip + 8];
        localDW->A6_g[A6_tmp_0] += F[A6_tmp + 3] * A[ip + 12];
      }
    }

    for (e_k = 0; e_k < 16; e_k++) {
      F[e_k] = localDW->A6_g[e_k];
      localDW->V[e_k] = ((90.0 * localDW->V[e_k] + 110880.0 * A6[e_k]) +
                         3.027024E+7 * A4[e_k]) + 2.0756736E+9 * A2[e_k];
    }

    localDW->d_n = 1.76432256E+10;
    break;

   default:
    for (e_k = 0; e_k < 16; e_k++) {
      localDW->d_n = A2[e_k];
      localDW->A4_l = A4[e_k];
      localDW->A6_p = A6[e_k];
      F[e_k] = (3.352212864E+10 * localDW->A6_p + 1.05594705216E+13 *
                localDW->A4_l) + 1.1873537964288E+15 * localDW->d_n;
      localDW->V[e_k] = (16380.0 * localDW->A4_l + localDW->A6_p) + 4.08408E+7 *
        localDW->d_n;
    }

    F[0] += 3.238237626624E+16;
    F[5] += 3.238237626624E+16;
    F[10] += 3.238237626624E+16;
    F[15] += 3.238237626624E+16;
    for (e_k = 0; e_k < 4; e_k++) {
      for (ip = 0; ip < 4; ip++) {
        A6_tmp = ip << 2;
        A6_tmp_0 = A6_tmp + e_k;
        localDW->A6_g[A6_tmp_0] = (((localDW->V[A6_tmp + 1] * A6[e_k + 4] +
          localDW->V[A6_tmp] * A6[e_k]) + localDW->V[A6_tmp + 2] * A6[e_k + 8])
          + localDW->V[A6_tmp + 3] * A6[e_k + 12]) + F[A6_tmp_0];
      }
    }

    for (e_k = 0; e_k < 4; e_k++) {
      for (ip = 0; ip < 4; ip++) {
        A6_tmp_0 = e_k << 2;
        A6_tmp = ip + A6_tmp_0;
        F[A6_tmp] = 0.0;
        F[A6_tmp] += localDW->A6_g[A6_tmp_0] * A[ip];
        F[A6_tmp] += localDW->A6_g[A6_tmp_0 + 1] * A[ip + 4];
        F[A6_tmp] += localDW->A6_g[A6_tmp_0 + 2] * A[ip + 8];
        F[A6_tmp] += localDW->A6_g[A6_tmp_0 + 3] * A[ip + 12];
      }
    }

    for (e_k = 0; e_k < 16; e_k++) {
      localDW->A6_g[e_k] = (182.0 * A6[e_k] + 960960.0 * A4[e_k]) +
        1.32324192E+9 * A2[e_k];
    }

    for (e_k = 0; e_k < 4; e_k++) {
      for (ip = 0; ip < 4; ip++) {
        A6_tmp = ip << 2;
        A6_tmp_0 = A6_tmp + e_k;
        localDW->V[A6_tmp_0] = (((((localDW->A6_g[A6_tmp + 1] * A6[e_k + 4] +
          localDW->A6_g[A6_tmp] * A6[e_k]) + localDW->A6_g[A6_tmp + 2] * A6[e_k
          + 8]) + localDW->A6_g[A6_tmp + 3] * A6[e_k + 12]) + A6[A6_tmp_0] *
          6.704425728E+11) + A4[A6_tmp_0] * 1.29060195264E+14) + A2[A6_tmp_0] *
          7.7717703038976E+15;
      }
    }

    localDW->d_n = 6.476475253248E+16;
    break;
  }

  localDW->V[0] += localDW->d_n;
  localDW->V[5] += localDW->d_n;
  localDW->V[10] += localDW->d_n;
  localDW->V[15] += localDW->d_n;
  for (e_k = 0; e_k < 16; e_k++) {
    localDW->d_n = F[e_k];
    localDW->V[e_k] -= localDW->d_n;
    F[e_k] = 2.0 * localDW->d_n;
  }

  xzgetrf(localDW->V, localDW->ipiv_n, &e_k, localDW);
  for (e_k = 0; e_k < 3; e_k++) {
    ip = localDW->ipiv_n[e_k];
    if (e_k + 1 != ip) {
      localDW->d_n = F[e_k];
      F[e_k] = F[ip - 1];
      F[ip - 1] = localDW->d_n;
      localDW->d_n = F[e_k + 4];
      F[e_k + 4] = F[ip + 3];
      F[ip + 3] = localDW->d_n;
      localDW->d_n = F[e_k + 8];
      F[e_k + 8] = F[ip + 7];
      F[ip + 7] = localDW->d_n;
      localDW->d_n = F[e_k + 12];
      F[e_k + 12] = F[ip + 11];
      F[ip + 11] = localDW->d_n;
    }
  }

  for (e_k = 0; e_k < 4; e_k++) {
    ip = e_k << 2;
    if (F[ip] != 0.0) {
      for (A6_tmp_0 = 2; A6_tmp_0 < 5; A6_tmp_0++) {
        A6_tmp = (A6_tmp_0 + ip) - 1;
        F[A6_tmp] -= localDW->V[A6_tmp_0 - 1] * F[ip];
      }
    }

    if (F[ip + 1] != 0.0) {
      for (A6_tmp_0 = 3; A6_tmp_0 < 5; A6_tmp_0++) {
        A6_tmp = (A6_tmp_0 + ip) - 1;
        F[A6_tmp] -= F[ip + 1] * localDW->V[A6_tmp_0 + 3];
      }
    }

    localDW->d_n = F[ip + 2];
    if (localDW->d_n != 0.0) {
      F[ip + 3] -= localDW->d_n * localDW->V[11];
    }
  }

  xtrsm(localDW->V, F);
  F[0]++;
  F[5]++;
  F[10]++;
  F[15]++;
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void recomputeBlockDiag(const real_T A[16], real_T F[16], const int32_T
  blockFormat[3])
{
  real_T avg;
  real_T expa11;
  real_T expa22;
  real_T u1;
  switch (blockFormat[0]) {
   case 0:
    break;

   case 1:
    expa11 = exp(A[0]);
    expa22 = exp(A[5]);
    avg = (A[0] + A[5]) / 2.0;
    u1 = fabs(A[0] - A[5]) / 2.0;
    if ((avg >= u1) || rtIsNaN(u1)) {
      u1 = avg;
    }

    if (u1 < 709.782712893384) {
      u1 = (A[5] - A[0]) / 2.0;
      if (u1 == 0.0) {
        u1 = 1.0;
      } else {
        u1 = sinh(u1) / u1;
      }

      avg = A[4] * exp(avg) * u1;
    } else {
      avg = (expa22 - expa11) * A[4] / (A[5] - A[0]);
    }

    F[0] = expa11;
    F[4] = avg;
    F[5] = expa22;
    break;

   case 2:
    expa11 = sqrt(fabs(A[1] * A[4]));
    expa22 = exp(A[0]);
    if (expa11 == 0.0) {
      avg = 1.0;
    } else {
      avg = sin(expa11) / expa11;
    }

    F[0] = expa22 * cos(expa11);
    F[1] = expa22 * A[1] * avg;
    F[4] = expa22 * A[4] * avg;
    F[5] = F[0];
    break;
  }

  switch (blockFormat[1]) {
   case 0:
    break;

   case 1:
    expa11 = exp(A[5]);
    expa22 = exp(A[10]);
    avg = (A[5] + A[10]) / 2.0;
    u1 = fabs(A[5] - A[10]) / 2.0;
    if ((avg >= u1) || rtIsNaN(u1)) {
      u1 = avg;
    }

    if (u1 < 709.782712893384) {
      u1 = (A[10] - A[5]) / 2.0;
      if (u1 == 0.0) {
        u1 = 1.0;
      } else {
        u1 = sinh(u1) / u1;
      }

      avg = A[9] * exp(avg) * u1;
    } else {
      avg = (expa22 - expa11) * A[9] / (A[10] - A[5]);
    }

    F[5] = expa11;
    F[9] = avg;
    F[10] = expa22;
    break;

   case 2:
    expa11 = sqrt(fabs(A[6] * A[9]));
    expa22 = exp(A[5]);
    if (expa11 == 0.0) {
      avg = 1.0;
    } else {
      avg = sin(expa11) / expa11;
    }

    F[5] = expa22 * cos(expa11);
    F[6] = expa22 * A[6] * avg;
    F[9] = expa22 * A[9] * avg;
    F[10] = F[5];
    break;
  }

  switch (blockFormat[2]) {
   case 0:
    break;

   case 1:
    expa11 = exp(A[10]);
    expa22 = exp(A[15]);
    avg = (A[10] + A[15]) / 2.0;
    u1 = fabs(A[10] - A[15]) / 2.0;
    if ((avg >= u1) || rtIsNaN(u1)) {
      u1 = avg;
    }

    if (u1 < 709.782712893384) {
      u1 = (A[15] - A[10]) / 2.0;
      if (u1 == 0.0) {
        u1 = 1.0;
      } else {
        u1 = sinh(u1) / u1;
      }

      avg = A[14] * exp(avg) * u1;
    } else {
      avg = (expa22 - expa11) * A[14] / (A[15] - A[10]);
    }

    F[10] = expa11;
    F[14] = avg;
    F[15] = expa22;
    break;

   case 2:
    expa11 = sqrt(fabs(A[11] * A[14]));
    expa22 = exp(A[10]);
    if (expa11 == 0.0) {
      avg = 1.0;
    } else {
      avg = sin(expa11) / expa11;
    }

    F[10] = expa22 * cos(expa11);
    F[11] = expa22 * A[11] * avg;
    F[14] = expa22 * A[14] * avg;
    F[15] = F[10];
    break;
  }

  if (blockFormat[2] == 0) {
    F[15] = exp(A[15]);
  }
}

/* Function for MATLAB Function: '<S113>/MATLAB Function' */
static void expm(real_T A[16], real_T F[16], DW_VandD *localDW)
{
  int32_T exitg1;
  boolean_T exitg2;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  boolean_T recomputeDiags;
  if (anyNonFinite(A)) {
    for (localDW->e_i = 0; localDW->e_i < 16; localDW->e_i++) {
      F[localDW->e_i] = (rtNaN);
    }
  } else {
    recomputeDiags = true;
    localDW->e_j = 0;
    exitg2 = false;
    while ((!exitg2) && (localDW->e_j < 4)) {
      localDW->e_i = 0;
      do {
        exitg1 = 0;
        if (localDW->e_i < 4) {
          if ((localDW->e_i != localDW->e_j) && (!(A[(localDW->e_j << 2) +
                localDW->e_i] == 0.0))) {
            recomputeDiags = false;
            exitg1 = 1;
          } else {
            localDW->e_i++;
          }
        } else {
          localDW->e_j++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }

    if (recomputeDiags) {
      memset(&F[0], 0, sizeof(real_T) << 4U);
      F[0] = exp(A[0]);
      F[5] = exp(A[5]);
      F[10] = exp(A[10]);
      F[15] = exp(A[15]);
    } else {
      recomputeDiags = true;
      localDW->e_j = 0;
      exitg2 = false;
      while ((!exitg2) && (localDW->e_j < 4)) {
        localDW->e_i = 0;
        do {
          exitg1 = 0;
          if (localDW->e_i <= localDW->e_j) {
            if (!(A[(localDW->e_j << 2) + localDW->e_i] == A[(localDW->e_i << 2)
                  + localDW->e_j])) {
              recomputeDiags = false;
              exitg1 = 1;
            } else {
              localDW->e_i++;
            }
          } else {
            localDW->e_j++;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }

      if (recomputeDiags) {
        schur(A, localDW->Q, localDW->T, localDW);
        for (localDW->b_s = 0; localDW->b_s < 4; localDW->b_s++) {
          localDW->e_i = localDW->b_s << 2;
          localDW->exptj = exp(localDW->T[localDW->e_i + localDW->b_s]);
          localDW->A2[localDW->e_i] = localDW->Q[localDW->e_i] * localDW->exptj;
          localDW->A2[localDW->e_i + 1] = localDW->Q[localDW->e_i + 1] *
            localDW->exptj;
          localDW->A2[localDW->e_i + 2] = localDW->Q[localDW->e_i + 2] *
            localDW->exptj;
          localDW->A2[localDW->e_i + 3] = localDW->Q[localDW->e_i + 3] *
            localDW->exptj;
        }

        for (localDW->e_i = 0; localDW->e_i < 4; localDW->e_i++) {
          for (localDW->i = 0; localDW->i < 4; localDW->i++) {
            localDW->A2_tmp = (localDW->i << 2) + localDW->e_i;
            F[localDW->A2_tmp] = 0.0;
            F[localDW->A2_tmp] += localDW->A2[localDW->e_i] * localDW->Q
              [localDW->i];
            F[localDW->A2_tmp] += localDW->A2[localDW->e_i + 4] * localDW->
              Q[localDW->i + 4];
            F[localDW->A2_tmp] += localDW->A2[localDW->e_i + 8] * localDW->
              Q[localDW->i + 8];
            F[localDW->A2_tmp] += localDW->A2[localDW->e_i + 12] * localDW->
              Q[localDW->i + 12];
          }
        }

        for (localDW->e_i = 0; localDW->e_i < 4; localDW->e_i++) {
          localDW->A2_tmp = localDW->e_i << 2;
          localDW->A2[localDW->A2_tmp] = (F[localDW->A2_tmp] + F[localDW->e_i]) /
            2.0;
          localDW->A2[localDW->A2_tmp + 1] = (F[localDW->A2_tmp + 1] + F
            [localDW->e_i + 4]) / 2.0;
          localDW->A2[localDW->A2_tmp + 2] = (F[localDW->A2_tmp + 2] + F
            [localDW->e_i + 8]) / 2.0;
          localDW->A2[localDW->A2_tmp + 3] = (F[localDW->A2_tmp + 3] + F
            [localDW->e_i + 12]) / 2.0;
        }

        memcpy(&F[0], &localDW->A2[0], sizeof(real_T) << 4U);
      } else {
        recomputeDiags = true;
        localDW->e_j = 3;
        while (recomputeDiags && (localDW->e_j <= 4)) {
          localDW->e_i = localDW->e_j;
          while (recomputeDiags && (localDW->e_i <= 4)) {
            recomputeDiags = (A[(((localDW->e_j - 3) << 2) + localDW->e_i) - 1] ==
                              0.0);
            localDW->e_i++;
          }

          localDW->e_j++;
        }

        if (recomputeDiags) {
          localDW->e_j = 1;
          exitg2 = false;
          while ((!exitg2) && (localDW->e_j - 1 < 3)) {
            localDW->e_i = ((localDW->e_j - 1) << 2) + localDW->e_j;
            localDW->exptj = A[localDW->e_i];
            if (localDW->exptj != 0.0) {
              if ((localDW->e_j != 3) && (A[((localDW->e_j << 2) + localDW->e_j)
                   + 1] != 0.0)) {
                recomputeDiags = false;
                exitg2 = true;
              } else {
                localDW->i = (localDW->e_j << 2) + localDW->e_j;
                if (A[localDW->e_i - 1] != A[localDW->i]) {
                  recomputeDiags = false;
                  exitg2 = true;
                } else {
                  localDW->d6 = A[localDW->i - 1];
                  if (!rtIsNaN(localDW->exptj)) {
                    if (localDW->exptj < 0.0) {
                      localDW->exptj = -1.0;
                    } else {
                      localDW->exptj = (localDW->exptj > 0.0);
                    }
                  }

                  if (!rtIsNaN(localDW->d6)) {
                    if (localDW->d6 < 0.0) {
                      localDW->d6 = -1.0;
                    } else {
                      localDW->d6 = (localDW->d6 > 0.0);
                    }
                  }

                  if (localDW->exptj * localDW->d6 != -1.0) {
                    recomputeDiags = false;
                    exitg2 = true;
                  } else {
                    localDW->e_j++;
                  }
                }
              }
            } else {
              localDW->e_j++;
            }
          }
        }

        localDW->exptj = 0.0;
        for (localDW->e_i = 0; localDW->e_i < 4; localDW->e_i++) {
          for (localDW->i = 0; localDW->i < 4; localDW->i++) {
            localDW->e_j = localDW->e_i << 2;
            localDW->A2_tmp = localDW->i + localDW->e_j;
            localDW->A2[localDW->A2_tmp] = 0.0;
            localDW->A2[localDW->A2_tmp] += A[localDW->e_j] * A[localDW->i];
            localDW->A2[localDW->A2_tmp] += A[localDW->e_j + 1] * A[localDW->i +
              4];
            localDW->A2[localDW->A2_tmp] += A[localDW->e_j + 2] * A[localDW->i +
              8];
            localDW->A2[localDW->A2_tmp] += A[localDW->e_j + 3] * A[localDW->i +
              12];
          }
        }

        for (localDW->e_i = 0; localDW->e_i < 4; localDW->e_i++) {
          for (localDW->i = 0; localDW->i < 4; localDW->i++) {
            localDW->e_j = localDW->i << 2;
            localDW->A2_tmp = localDW->e_i + localDW->e_j;
            localDW->A4[localDW->A2_tmp] = 0.0;
            localDW->A4[localDW->A2_tmp] += localDW->A2[localDW->e_j] *
              localDW->A2[localDW->e_i];
            localDW->A4[localDW->A2_tmp] += localDW->A2[localDW->e_j + 1] *
              localDW->A2[localDW->e_i + 4];
            localDW->A4[localDW->A2_tmp] += localDW->A2[localDW->e_j + 2] *
              localDW->A2[localDW->e_i + 8];
            localDW->A4[localDW->A2_tmp] += localDW->A2[localDW->e_j + 3] *
              localDW->A2[localDW->e_i + 12];
          }

          for (localDW->i = 0; localDW->i < 4; localDW->i++) {
            localDW->e_j = localDW->i << 2;
            localDW->A2_tmp = localDW->e_i + localDW->e_j;
            localDW->A6[localDW->A2_tmp] = 0.0;
            localDW->A6[localDW->A2_tmp] += localDW->A2[localDW->e_j] *
              localDW->A4[localDW->e_i];
            localDW->A6[localDW->A2_tmp] += localDW->A2[localDW->e_j + 1] *
              localDW->A4[localDW->e_i + 4];
            localDW->A6[localDW->A2_tmp] += localDW->A2[localDW->e_j + 2] *
              localDW->A4[localDW->e_i + 8];
            localDW->A6[localDW->A2_tmp] += localDW->A2[localDW->e_j + 3] *
              localDW->A4[localDW->e_i + 12];
          }
        }

        localDW->d6 = rt_powd_snf(norm(localDW->A6), 0.16666666666666666);
        localDW->eta1 = rt_powd_snf(norm(localDW->A4), 0.25);
        if ((!(localDW->eta1 >= localDW->d6)) && (!rtIsNaN(localDW->d6))) {
          localDW->eta1 = localDW->d6;
        }

        guard1 = false;
        guard2 = false;
        guard3 = false;
        guard4 = false;
        if (localDW->eta1 <= 0.01495585217958292) {
          for (localDW->e_j = 0; localDW->e_j < 16; localDW->e_j++) {
            localDW->dv[localDW->e_j] = 0.19285012468241128 * fabs(A
              [localDW->e_j]);
          }

          mpower(localDW->dv, 7.0, localDW->Q, localDW);
          localDW->b_varargin_1 = log2_p(norm(localDW->Q) / norm(A) * 2.0 /
            2.2204460492503131E-16) / 6.0;
          localDW->b_varargin_1 = ceil(localDW->b_varargin_1);
          if (!(localDW->b_varargin_1 >= 0.0)) {
            localDW->b_varargin_1 = 0.0;
          }

          if (localDW->b_varargin_1 == 0.0) {
            localDW->b_s = 3;
          } else {
            guard4 = true;
          }
        } else {
          guard4 = true;
        }

        if (guard4) {
          if (localDW->eta1 <= 0.253939833006323) {
            for (localDW->e_j = 0; localDW->e_j < 16; localDW->e_j++) {
              localDW->dv[localDW->e_j] = 0.12321872304378752 * fabs(A
                [localDW->e_j]);
            }

            mpower(localDW->dv, 11.0, localDW->Q, localDW);
            localDW->eta1 = log2_p(norm(localDW->Q) / norm(A) * 2.0 /
              2.2204460492503131E-16) / 10.0;
            localDW->b_varargin_1 = ceil(localDW->eta1);
            if (!(localDW->b_varargin_1 >= 0.0)) {
              localDW->b_varargin_1 = 0.0;
            }

            if (localDW->b_varargin_1 == 0.0) {
              localDW->b_s = 5;
            } else {
              guard3 = true;
            }
          } else {
            guard3 = true;
          }
        }

        if (guard3) {
          mpower(localDW->A4, 2.0, localDW->Q, localDW);
          localDW->eta1 = rt_powd_snf(norm(localDW->Q), 0.125);
          if ((!(localDW->d6 >= localDW->eta1)) && (!rtIsNaN(localDW->eta1))) {
            localDW->d6 = localDW->eta1;
          }

          if (localDW->d6 <= 0.95041789961629319) {
            for (localDW->e_j = 0; localDW->e_j < 16; localDW->e_j++) {
              localDW->dv[localDW->e_j] = 0.090475336558796943 * fabs(A
                [localDW->e_j]);
            }

            mpower(localDW->dv, 15.0, localDW->Q, localDW);
            localDW->b_varargin_1 = log2_p(norm(localDW->Q) / norm(A) * 2.0 /
              2.2204460492503131E-16) / 14.0;
            localDW->b_varargin_1 = ceil(localDW->b_varargin_1);
            if (!(localDW->b_varargin_1 >= 0.0)) {
              localDW->b_varargin_1 = 0.0;
            }

            if (localDW->b_varargin_1 == 0.0) {
              localDW->b_s = 7;
            } else {
              guard2 = true;
            }
          } else {
            guard2 = true;
          }
        }

        if (guard2) {
          if (localDW->d6 <= 2.097847961257068) {
            for (localDW->e_j = 0; localDW->e_j < 16; localDW->e_j++) {
              localDW->dv[localDW->e_j] = 0.071467735648795785 * fabs(A
                [localDW->e_j]);
            }

            mpower(localDW->dv, 19.0, localDW->Q, localDW);
            localDW->b_varargin_1 = log2_p(norm(localDW->Q) / norm(A) * 2.0 /
              2.2204460492503131E-16) / 18.0;
            localDW->b_varargin_1 = ceil(localDW->b_varargin_1);
            if (!(localDW->b_varargin_1 >= 0.0)) {
              localDW->b_varargin_1 = 0.0;
            }

            if (localDW->b_varargin_1 == 0.0) {
              localDW->b_s = 9;
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }
        }

        if (guard1) {
          for (localDW->e_i = 0; localDW->e_i < 4; localDW->e_i++) {
            for (localDW->i = 0; localDW->i < 4; localDW->i++) {
              localDW->e_j = localDW->e_i << 2;
              localDW->A2_tmp = localDW->i + localDW->e_j;
              localDW->Q[localDW->A2_tmp] = 0.0;
              localDW->Q[localDW->A2_tmp] += localDW->A6[localDW->e_j] *
                localDW->A4[localDW->i];
              localDW->Q[localDW->A2_tmp] += localDW->A6[localDW->e_j + 1] *
                localDW->A4[localDW->i + 4];
              localDW->Q[localDW->A2_tmp] += localDW->A6[localDW->e_j + 2] *
                localDW->A4[localDW->i + 8];
              localDW->Q[localDW->A2_tmp] += localDW->A6[localDW->e_j + 3] *
                localDW->A4[localDW->i + 12];
            }
          }

          localDW->exptj = rt_powd_snf(norm(localDW->Q), 0.1);
          if ((localDW->eta1 >= localDW->exptj) || rtIsNaN(localDW->exptj)) {
            localDW->exptj = localDW->eta1;
          }

          if ((localDW->d6 <= localDW->exptj) || rtIsNaN(localDW->exptj)) {
            localDW->exptj = localDW->d6;
          }

          localDW->exptj = log2_p(localDW->exptj / 5.3719203511481517);
          localDW->exptj = ceil(localDW->exptj);
          if (!(localDW->exptj >= 0.0)) {
            localDW->exptj = 0.0;
          }

          localDW->d6 = rt_powd_snf(2.0, localDW->exptj);
          for (localDW->e_j = 0; localDW->e_j < 16; localDW->e_j++) {
            localDW->eta1 = A[localDW->e_j] / localDW->d6;
            localDW->dv[localDW->e_j] = 0.05031554467093536 * fabs(localDW->eta1);
            localDW->T[localDW->e_j] = localDW->eta1;
          }

          mpower(localDW->dv, 27.0, localDW->Q, localDW);
          localDW->d6 = log2_p(norm(localDW->Q) / norm(localDW->T) * 2.0 /
                               2.2204460492503131E-16) / 26.0;
          localDW->b_varargin_1 = ceil(localDW->d6);
          if (!(localDW->b_varargin_1 >= 0.0)) {
            localDW->b_varargin_1 = 0.0;
          }

          localDW->exptj += localDW->b_varargin_1;
          if (rtIsInf(localDW->exptj)) {
            localDW->d6 = norm(A) / 5.3719203511481517;
            if ((!rtIsInf(localDW->d6)) && (!rtIsNaN(localDW->d6))) {
              localDW->d6 = frexp(localDW->d6, &localDW->b_s);
            } else {
              localDW->b_s = 0;
            }

            localDW->exptj = localDW->b_s;
            if (localDW->d6 == 0.5) {
              localDW->exptj = (real_T)localDW->b_s - 1.0;
            }
          }

          localDW->b_s = 13;
        }

        if (localDW->exptj != 0.0) {
          localDW->d6 = rt_powd_snf(2.0, localDW->exptj);
          for (localDW->e_i = 0; localDW->e_i < 16; localDW->e_i++) {
            A[localDW->e_i] /= localDW->d6;
          }

          localDW->d6 = rt_powd_snf(2.0, 2.0 * localDW->exptj);
          for (localDW->e_i = 0; localDW->e_i < 16; localDW->e_i++) {
            localDW->A2[localDW->e_i] /= localDW->d6;
          }

          localDW->d6 = rt_powd_snf(2.0, 4.0 * localDW->exptj);
          for (localDW->e_i = 0; localDW->e_i < 16; localDW->e_i++) {
            localDW->A4[localDW->e_i] /= localDW->d6;
          }

          localDW->d6 = rt_powd_snf(2.0, 6.0 * localDW->exptj);
          for (localDW->e_i = 0; localDW->e_i < 16; localDW->e_i++) {
            localDW->A6[localDW->e_i] /= localDW->d6;
          }
        }

        if (recomputeDiags) {
          localDW->blockFormat[0] = 0;
          localDW->blockFormat[1] = 0;
          localDW->blockFormat[2] = 0;
          localDW->e_j = 1;
          while (localDW->e_j < 3) {
            localDW->d6 = A[((localDW->e_j - 1) << 2) + localDW->e_j];
            if (localDW->d6 != 0.0) {
              localDW->blockFormat[localDW->e_j - 1] = 2;
              localDW->blockFormat[localDW->e_j] = 0;
              localDW->e_j += 2;
            } else if ((localDW->d6 == 0.0) && (A[((localDW->e_j << 2) +
                         localDW->e_j) + 1] == 0.0)) {
              localDW->blockFormat[localDW->e_j - 1] = 1;
              localDW->e_j++;
            } else {
              localDW->blockFormat[localDW->e_j - 1] = 0;
              localDW->e_j++;
            }
          }

          if (A[11] != 0.0) {
            localDW->blockFormat[2] = 2;
          } else {
            switch (localDW->blockFormat[1]) {
             case 0:
              localDW->blockFormat[2] = 1;
              break;

             case 1:
              localDW->blockFormat[2] = 1;
              break;
            }
          }
        }

        padeApproximation(A, localDW->A2, localDW->A4, localDW->A6, localDW->b_s,
                          F, localDW);
        if (recomputeDiags) {
          recomputeBlockDiag(A, F, localDW->blockFormat);
        }

        localDW->b_s = 0;
        while (localDW->b_s <= (int32_T)localDW->exptj - 1) {
          for (localDW->e_i = 0; localDW->e_i < 4; localDW->e_i++) {
            for (localDW->i = 0; localDW->i < 4; localDW->i++) {
              localDW->e_j = localDW->i << 2;
              localDW->A2_tmp = localDW->e_i + localDW->e_j;
              localDW->A2[localDW->A2_tmp] = 0.0;
              localDW->A2[localDW->A2_tmp] += F[localDW->e_j] * F[localDW->e_i];
              localDW->A2[localDW->A2_tmp] += F[localDW->e_j + 1] * F
                [localDW->e_i + 4];
              localDW->A2[localDW->A2_tmp] += F[localDW->e_j + 2] * F
                [localDW->e_i + 8];
              localDW->A2[localDW->A2_tmp] += F[localDW->e_j + 3] * F
                [localDW->e_i + 12];
            }
          }

          memcpy(&F[0], &localDW->A2[0], sizeof(real_T) << 4U);
          if (recomputeDiags) {
            for (localDW->e_i = 0; localDW->e_i < 16; localDW->e_i++) {
              A[localDW->e_i] *= 2.0;
            }

            recomputeBlockDiag(A, F, localDW->blockFormat);
          }

          localDW->b_s++;
        }
      }
    }
  }
}

/* Function for MATLAB Function: '<S113>/MATLAB Function1' */
static real_T eml_rand_shr3cong(uint32_T state[2])
{
  real_T r;
  real_T s;
  real_T x;
  real_T y;
  int32_T j;
  uint32_T icng;
  uint32_T jsr;
  uint32_T ui;
  static const real_T b[65] = { 0.340945, 0.4573146, 0.5397793, 0.6062427,
    0.6631691, 0.7136975, 0.7596125, 0.8020356, 0.8417227, 0.8792102, 0.9148948,
    0.9490791, 0.9820005, 1.0138492, 1.044781, 1.0749254, 1.1043917, 1.1332738,
    1.161653, 1.189601, 1.2171815, 1.2444516, 1.2714635, 1.298265, 1.3249008,
    1.3514125, 1.3778399, 1.4042211, 1.4305929, 1.4569915, 1.4834527, 1.5100122,
    1.5367061, 1.5635712, 1.5906454, 1.617968, 1.6455802, 1.6735255, 1.7018503,
    1.7306045, 1.7598422, 1.7896223, 1.8200099, 1.851077, 1.8829044, 1.9155831,
    1.9492166, 1.9839239, 2.0198431, 2.0571356, 2.095993, 2.136645, 2.1793713,
    2.2245175, 2.2725186, 2.3239338, 2.3795008, 2.4402218, 2.5075117, 2.5834658,
    2.6713916, 2.7769942, 2.7769942, 2.7769942, 2.7769942 };

  icng = 69069U * state[0] + 1234567U;
  jsr = state[1] << 13 ^ state[1];
  jsr ^= jsr >> 17;
  jsr ^= jsr << 5;
  ui = icng + jsr;
  j = (int32_T)((ui & 63U) + 1U);
  r = (real_T)(int32_T)ui * 4.6566128730773926E-10 * b[j];
  x = fabs(r);
  y = b[j - 1];
  if (!(x <= y)) {
    x = (x - y) / (b[j] - y);
    icng = 69069U * icng + 1234567U;
    jsr ^= jsr << 13;
    jsr ^= jsr >> 17;
    jsr ^= jsr << 5;
    y = (real_T)(int32_T)(icng + jsr) * 2.328306436538696E-10 + 0.5;
    s = x + y;
    if (s > 1.301198) {
      if (r < 0.0) {
        r = 0.4878992 * x - 0.4878992;
      } else {
        r = 0.4878992 - 0.4878992 * x;
      }
    } else if (!(s <= 0.9689279)) {
      x = 0.4878992 - 0.4878992 * x;
      if (y > 12.67706 - exp(-0.5 * x * x) * 12.37586) {
        if (r < 0.0) {
          r = -x;
        } else {
          r = x;
        }
      } else if (!(exp(-0.5 * b[j] * b[j]) + y * 0.01958303 / b[j] <= exp(-0.5 *
        r * r))) {
        do {
          icng = 69069U * icng + 1234567U;
          jsr ^= jsr << 13;
          jsr ^= jsr >> 17;
          jsr ^= jsr << 5;
          x = log((real_T)(int32_T)(icng + jsr) * 2.328306436538696E-10 + 0.5) /
            2.776994;
          icng = 69069U * icng + 1234567U;
          jsr ^= jsr << 13;
          jsr ^= jsr >> 17;
          jsr ^= jsr << 5;
        } while (!(log((real_T)(int32_T)(icng + jsr) * 2.328306436538696E-10 +
                       0.5) * -2.0 > x * x));

        if (r < 0.0) {
          r = x - 2.776994;
        } else {
          r = 2.776994 - x;
        }
      }
    }
  }

  state[0] = icng;
  state[1] = jsr;
  return r;
}

/* Function for MATLAB Function: '<S113>/MATLAB Function1' */
static void genrandu(uint32_T s, uint32_T *state, real_T *r)
{
  int32_T hi;
  uint32_T a;
  uint32_T b;
  hi = (int32_T)(s / 127773U);
  a = (s - hi * 127773U) * 16807U;
  b = 2836U * hi;
  if (a < b) {
    *state = ~(b - a) & 2147483647U;
  } else {
    *state = a - b;
  }

  *r = (real_T)*state * 4.6566128752457969E-10;
}

/* Function for MATLAB Function: '<S113>/MATLAB Function1' */
static void genrand_uint32_vector(uint32_T mt[625], uint32_T u[2])
{
  int32_T j;
  int32_T kk;
  uint32_T mti;
  uint32_T y;
  for (j = 0; j < 2; j++) {
    mti = mt[624] + 1U;
    if (mt[624] + 1U >= 625U) {
      for (kk = 0; kk < 227; kk++) {
        mti = (mt[kk + 1] & 2147483647U) | (mt[kk] & 2147483648U);
        if ((mti & 1U) == 0U) {
          mti >>= 1U;
        } else {
          mti = mti >> 1U ^ 2567483615U;
        }

        mt[kk] = mt[kk + 397] ^ mti;
      }

      for (kk = 0; kk < 396; kk++) {
        mti = (mt[kk + 227] & 2147483648U) | (mt[kk + 228] & 2147483647U);
        if ((mti & 1U) == 0U) {
          mti >>= 1U;
        } else {
          mti = mti >> 1U ^ 2567483615U;
        }

        mt[kk + 227] = mt[kk] ^ mti;
      }

      mti = (mt[623] & 2147483648U) | (mt[0] & 2147483647U);
      if ((mti & 1U) == 0U) {
        mti >>= 1U;
      } else {
        mti = mti >> 1U ^ 2567483615U;
      }

      mt[623] = mt[396] ^ mti;
      mti = 1U;
    }

    y = mt[(int32_T)mti - 1];
    mt[624] = mti;
    y ^= y >> 11U;
    y ^= y << 7U & 2636928640U;
    y ^= y << 15U & 4022730752U;
    u[j] = y >> 18U ^ y;
  }
}

/* Function for MATLAB Function: '<S113>/MATLAB Function1' */
static real_T genrandu_d(uint32_T mt[625])
{
  real_T r;
  int32_T exitg1;
  int32_T k;
  uint32_T u[2];
  uint32_T b_r;
  boolean_T b_isvalid;
  boolean_T exitg2;

  /* ========================= COPYRIGHT NOTICE ============================ */
  /*  This is a uniform (0,1) pseudorandom number generator based on:        */
  /*                                                                         */
  /*  A C-program for MT19937, with initialization improved 2002/1/26.       */
  /*  Coded by Takuji Nishimura and Makoto Matsumoto.                        */
  /*                                                                         */
  /*  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,      */
  /*  All rights reserved.                                                   */
  /*                                                                         */
  /*  Redistribution and use in source and binary forms, with or without     */
  /*  modification, are permitted provided that the following conditions     */
  /*  are met:                                                               */
  /*                                                                         */
  /*    1. Redistributions of source code must retain the above copyright    */
  /*       notice, this list of conditions and the following disclaimer.     */
  /*                                                                         */
  /*    2. Redistributions in binary form must reproduce the above copyright */
  /*       notice, this list of conditions and the following disclaimer      */
  /*       in the documentation and/or other materials provided with the     */
  /*       distribution.                                                     */
  /*                                                                         */
  /*    3. The names of its contributors may not be used to endorse or       */
  /*       promote products derived from this software without specific      */
  /*       prior written permission.                                         */
  /*                                                                         */
  /*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS    */
  /*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT      */
  /*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR  */
  /*  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT  */
  /*  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,  */
  /*  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT       */
  /*  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,  */
  /*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY  */
  /*  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT    */
  /*  (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE */
  /*  OF THIS  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  */
  /*                                                                         */
  /* =============================   END   ================================= */
  do {
    exitg1 = 0;
    genrand_uint32_vector(mt, u);
    r = ((real_T)(u[0] >> 5U) * 6.7108864E+7 + (real_T)(u[1] >> 6U)) *
      1.1102230246251565E-16;
    if (r == 0.0) {
      b_isvalid = ((mt[624] >= 1U) && (mt[624] < 625U));
      if (b_isvalid) {
        b_isvalid = false;
        k = 1;
        exitg2 = false;
        while ((!exitg2) && (k < 625)) {
          if (mt[k - 1] == 0U) {
            k++;
          } else {
            b_isvalid = true;
            exitg2 = true;
          }
        }
      }

      if (!b_isvalid) {
        b_r = 5489U;
        mt[0] = 5489U;
        for (k = 0; k < 623; k++) {
          b_r = ((b_r >> 30U ^ b_r) * 1812433253U + k) + 1U;
          mt[k + 1] = b_r;
        }

        mt[624] = 624U;
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return r;
}

/* Function for MATLAB Function: '<S113>/MATLAB Function1' */
static real_T randn(DW_VandD *localDW)
{
  real_T r;
  int32_T i;
  uint32_T c_state;
  uint32_T d_state;
  static const real_T b[257] = { 0.0, 0.215241895984875, 0.286174591792068,
    0.335737519214422, 0.375121332878378, 0.408389134611989, 0.43751840220787,
    0.46363433679088, 0.487443966139235, 0.50942332960209, 0.529909720661557,
    0.549151702327164, 0.567338257053817, 0.584616766106378, 0.601104617755991,
    0.61689699000775, 0.63207223638606, 0.646695714894993, 0.660822574244419,
    0.674499822837293, 0.687767892795788, 0.700661841106814, 0.713212285190975,
    0.725446140909999, 0.737387211434295, 0.749056662017815, 0.760473406430107,
    0.771654424224568, 0.782615023307232, 0.793369058840623, 0.80392911698997,
    0.814306670135215, 0.824512208752291, 0.834555354086381, 0.844444954909153,
    0.854189171008163, 0.863795545553308, 0.87327106808886, 0.882622229585165,
    0.891855070732941, 0.900975224461221, 0.909987953496718, 0.91889818364959,
    0.927710533401999, 0.936429340286575, 0.945058684468165, 0.953602409881086,
    0.96206414322304, 0.970447311064224, 0.978755155294224, 0.986990747099062,
    0.99515699963509, 1.00325667954467, 1.01129241744, 1.01926671746548,
    1.02718196603564, 1.03504043983344, 1.04284431314415, 1.05059566459093,
    1.05829648333067, 1.06594867476212, 1.07355406579244, 1.0811144097034,
    1.08863139065398, 1.09610662785202, 1.10354167942464, 1.11093804601357,
    1.11829717411934, 1.12562045921553, 1.13290924865253, 1.14016484436815,
    1.14738850542085, 1.15458145035993, 1.16174485944561, 1.16887987673083,
    1.17598761201545, 1.18306914268269, 1.19012551542669, 1.19715774787944,
    1.20416683014438, 1.2111537262437, 1.21811937548548, 1.22506469375653,
    1.23199057474614, 1.23889789110569, 1.24578749554863, 1.2526602218949,
    1.25951688606371, 1.26635828701823, 1.27318520766536, 1.27999841571382,
    1.28679866449324, 1.29358669373695, 1.30036323033084, 1.30712898903073,
    1.31388467315022, 1.32063097522106, 1.32736857762793, 1.33409815321936,
    1.3408203658964, 1.34753587118059, 1.35424531676263, 1.36094934303328,
    1.36764858359748, 1.37434366577317, 1.38103521107586, 1.38772383568998,
    1.39441015092814, 1.40109476367925, 1.4077782768464, 1.41446128977547,
    1.42114439867531, 1.42782819703026, 1.43451327600589, 1.44120022484872,
    1.44788963128058, 1.45458208188841, 1.46127816251028, 1.46797845861808,
    1.47468355569786, 1.48139403962819, 1.48811049705745, 1.49483351578049,
    1.50156368511546, 1.50830159628131, 1.51504784277671, 1.521803020761,
    1.52856772943771, 1.53534257144151, 1.542128153229, 1.54892508547417,
    1.55573398346918, 1.56255546753104, 1.56939016341512, 1.57623870273591,
    1.58310172339603, 1.58997987002419, 1.59687379442279, 1.60378415602609,
    1.61071162236983, 1.61765686957301, 1.62462058283303, 1.63160345693487,
    1.63860619677555, 1.64562951790478, 1.65267414708306, 1.65974082285818,
    1.66683029616166, 1.67394333092612, 1.68108070472517, 1.68824320943719,
    1.69543165193456, 1.70264685479992, 1.7098896570713, 1.71716091501782,
    1.72446150294804, 1.73179231405296, 1.73915426128591, 1.74654827828172,
    1.75397532031767, 1.76143636531891, 1.76893241491127, 1.77646449552452,
    1.78403365954944, 1.79164098655216, 1.79928758454972, 1.80697459135082,
    1.81470317596628, 1.82247454009388, 1.83028991968276, 1.83815058658281,
    1.84605785028518, 1.8540130597602, 1.86201760539967, 1.87007292107127,
    1.878180486293, 1.88634182853678, 1.8945585256707, 1.90283220855043,
    1.91116456377125, 1.91955733659319, 1.92801233405266, 1.93653142827569,
    1.94511656000868, 1.95376974238465, 1.96249306494436, 1.97128869793366,
    1.98015889690048, 1.98910600761744, 1.99813247135842, 2.00724083056053,
    2.0164337349062, 2.02571394786385, 2.03508435372962, 2.04454796521753,
    2.05410793165065, 2.06376754781173, 2.07353026351874, 2.0833996939983,
    2.09337963113879, 2.10347405571488, 2.11368715068665, 2.12402331568952,
    2.13448718284602, 2.14508363404789, 2.15581781987674, 2.16669518035431,
    2.17772146774029, 2.18890277162636, 2.20024554661128, 2.21175664288416,
    2.22344334009251, 2.23531338492992, 2.24737503294739, 2.25963709517379,
    2.27210899022838, 2.28480080272449, 2.29772334890286, 2.31088825060137,
    2.32430801887113, 2.33799614879653, 2.35196722737914, 2.36623705671729,
    2.38082279517208, 2.39574311978193, 2.41101841390112, 2.42667098493715,
    2.44272531820036, 2.4592083743347, 2.47614993967052, 2.49358304127105,
    2.51154444162669, 2.53007523215985, 2.54922155032478, 2.56903545268184,
    2.58957598670829, 2.61091051848882, 2.63311639363158, 2.65628303757674,
    2.68051464328574, 2.70593365612306, 2.73268535904401, 2.76094400527999,
    2.79092117400193, 2.82287739682644, 2.85713873087322, 2.89412105361341,
    2.93436686720889, 2.97860327988184, 3.02783779176959, 3.08352613200214,
    3.147889289518, 3.2245750520478, 3.32024473383983, 3.44927829856143,
    3.65415288536101, 3.91075795952492 };

  static const real_T c[257] = { 1.0, 0.977101701267673, 0.959879091800108,
    0.9451989534423, 0.932060075959231, 0.919991505039348, 0.908726440052131,
    0.898095921898344, 0.887984660755834, 0.878309655808918, 0.869008688036857,
    0.860033621196332, 0.851346258458678, 0.842915653112205, 0.834716292986884,
    0.826726833946222, 0.818929191603703, 0.811307874312656, 0.803849483170964,
    0.796542330422959, 0.789376143566025, 0.782341832654803, 0.775431304981187,
    0.768637315798486, 0.761953346836795, 0.755373506507096, 0.748892447219157,
    0.742505296340151, 0.736207598126863, 0.729995264561476, 0.72386453346863,
    0.717811932630722, 0.711834248878248, 0.705928501332754, 0.700091918136512,
    0.694321916126117, 0.688616083004672, 0.682972161644995, 0.677388036218774,
    0.671861719897082, 0.66639134390875, 0.660975147776663, 0.655611470579697,
    0.650298743110817, 0.645035480820822, 0.639820277453057, 0.634651799287624,
    0.629528779924837, 0.624450015547027, 0.619414360605834, 0.614420723888914,
    0.609468064925773, 0.604555390697468, 0.599681752619125, 0.594846243767987,
    0.590047996332826, 0.585286179263371, 0.580559996100791, 0.575868682972354,
    0.571211506735253, 0.566587763256165, 0.561996775814525, 0.557437893618766,
    0.552910490425833, 0.548413963255266, 0.543947731190026, 0.539511234256952,
    0.535103932380458, 0.530725304403662, 0.526374847171684, 0.522052074672322,
    0.517756517229756, 0.513487720747327, 0.509245245995748, 0.505028667943468,
    0.500837575126149, 0.49667156905249, 0.492530263643869, 0.488413284705458,
    0.484320269426683, 0.480250865909047, 0.476204732719506, 0.47218153846773,
    0.468180961405694, 0.464202689048174, 0.460246417812843, 0.456311852678716,
    0.452398706861849, 0.448506701507203, 0.444635565395739, 0.440785034665804,
    0.436954852547985, 0.433144769112652, 0.429354541029442, 0.425583931338022,
    0.421832709229496, 0.418100649837848, 0.414387534040891, 0.410693148270188,
    0.407017284329473, 0.403359739221114, 0.399720314980197, 0.396098818515832,
    0.392495061459315, 0.388908860018789, 0.385340034840077, 0.381788410873393,
    0.378253817245619, 0.374736087137891, 0.371235057668239, 0.367750569779032,
    0.364282468129004, 0.360830600989648, 0.357394820145781, 0.353974980800077,
    0.350570941481406, 0.347182563956794, 0.343809713146851, 0.340452257044522,
    0.337110066637006, 0.333783015830718, 0.330470981379163, 0.327173842813601,
    0.323891482376391, 0.320623784956905, 0.317370638029914, 0.314131931596337,
    0.310907558126286, 0.307697412504292, 0.30450139197665, 0.301319396100803,
    0.298151326696685, 0.294997087799962, 0.291856585617095, 0.288729728482183,
    0.285616426815502, 0.282516593083708, 0.279430141761638, 0.276356989295668,
    0.273297054068577, 0.270250256365875, 0.267216518343561, 0.264195763997261,
    0.261187919132721, 0.258192911337619, 0.255210669954662, 0.252241126055942,
    0.249284212418529, 0.246339863501264, 0.24340801542275, 0.240488605940501,
    0.237581574431238, 0.23468686187233, 0.231804410824339, 0.228934165414681,
    0.226076071322381, 0.223230075763918, 0.220396127480152, 0.217574176724331,
    0.214764175251174, 0.211966076307031, 0.209179834621125, 0.206405406397881,
    0.203642749310335, 0.200891822494657, 0.198152586545776, 0.195425003514135,
    0.192709036903589, 0.190004651670465, 0.187311814223801, 0.1846304924268,
    0.181960655599523, 0.179302274522848, 0.176655321443735, 0.174019770081839,
    0.171395595637506, 0.168782774801212, 0.166181285764482, 0.163591108232366,
    0.161012223437511, 0.158444614155925, 0.15588826472448, 0.153343161060263,
    0.150809290681846, 0.148286642732575, 0.145775208005994, 0.143274978973514,
    0.140785949814445, 0.138308116448551, 0.135841476571254, 0.133386029691669,
    0.130941777173644, 0.12850872228, 0.126086870220186, 0.123676228201597,
    0.12127680548479, 0.11888861344291, 0.116511665625611, 0.114145977827839,
    0.111791568163838, 0.109448457146812, 0.107116667774684, 0.104796225622487,
    0.102487158941935, 0.10018949876881, 0.0979032790388625, 0.095628536713009,
    0.093365311912691, 0.0911136480663738, 0.0888735920682759,
    0.0866451944505581, 0.0844285095703535, 0.082223595813203,
    0.0800305158146631, 0.0778493367020961, 0.0756801303589272,
    0.0735229737139814, 0.0713779490588905, 0.0692451443970068,
    0.0671246538277886, 0.065016577971243, 0.0629210244377582, 0.06083810834954,
    0.0587679529209339, 0.0567106901062031, 0.0546664613248891,
    0.0526354182767924, 0.0506177238609479, 0.0486135532158687,
    0.0466230949019305, 0.0446465522512946, 0.0426841449164746,
    0.0407361106559411, 0.0388027074045262, 0.0368842156885674,
    0.0349809414617162, 0.0330932194585786, 0.0312214171919203,
    0.0293659397581334, 0.0275272356696031, 0.0257058040085489,
    0.0239022033057959, 0.0221170627073089, 0.0203510962300445,
    0.0186051212757247, 0.0168800831525432, 0.0151770883079353,
    0.0134974506017399, 0.0118427578579079, 0.0102149714397015,
    0.00861658276939875, 0.00705087547137324, 0.00552240329925101,
    0.00403797259336304, 0.00260907274610216, 0.0012602859304986,
    0.000477467764609386 };

  int32_T exitg1;
  switch (localDW->method) {
   case 0U:
    switch (localDW->method_b) {
     case 4U:
      c_state = localDW->state_g;
      do {
        genrandu(c_state, &d_state, &localDW->c_u);
        genrandu(d_state, &c_state, &localDW->x_d);
        localDW->c_u = 2.0 * localDW->c_u - 1.0;
        localDW->x_d = 2.0 * localDW->x_d - 1.0;
        localDW->x_d = localDW->x_d * localDW->x_d + localDW->c_u * localDW->c_u;
      } while (!(localDW->x_d <= 1.0));

      r = sqrt(-2.0 * log(localDW->x_d) / localDW->x_d) * localDW->c_u;
      localDW->state_g = c_state;
      break;

     case 5U:
      r = eml_rand_shr3cong(localDW->state_o);
      break;

     default:
      do {
        exitg1 = 0;
        genrand_uint32_vector(localDW->state_d, localDW->u32);
        i = (int32_T)((localDW->u32[1] >> 24U) + 1U);
        r = (((real_T)(localDW->u32[0] >> 3U) * 1.6777216E+7 + (real_T)((int32_T)
               localDW->u32[1] & 16777215)) * 2.2204460492503131E-16 - 1.0) *
          b[i];
        if (fabs(r) <= b[i - 1]) {
          exitg1 = 1;
        } else if (i < 256) {
          localDW->x_d = genrandu_d(localDW->state_d);
          if ((c[i - 1] - c[i]) * localDW->x_d + c[i] < exp(-0.5 * r * r)) {
            exitg1 = 1;
          }
        } else {
          do {
            localDW->x_d = genrandu_d(localDW->state_d);
            localDW->x_d = log(localDW->x_d) * 0.273661237329758;
            localDW->c_u = genrandu_d(localDW->state_d);
          } while (!(-2.0 * log(localDW->c_u) > localDW->x_d * localDW->x_d));

          if (r < 0.0) {
            r = localDW->x_d - 3.65415288536101;
          } else {
            r = 3.65415288536101 - localDW->x_d;
          }

          exitg1 = 1;
        }
      } while (exitg1 == 0);
      break;
    }
    break;

   case 4U:
    c_state = localDW->state[0];
    do {
      genrandu(c_state, &d_state, &localDW->c_u);
      genrandu(d_state, &c_state, &localDW->x_d);
      localDW->c_u = 2.0 * localDW->c_u - 1.0;
      localDW->x_d = 2.0 * localDW->x_d - 1.0;
      localDW->x_d = localDW->x_d * localDW->x_d + localDW->c_u * localDW->c_u;
    } while (!(localDW->x_d <= 1.0));

    r = sqrt(-2.0 * log(localDW->x_d) / localDW->x_d) * localDW->c_u;
    localDW->state[0] = c_state;
    break;

   default:
    r = eml_rand_shr3cong(localDW->state);
    break;
  }

  return r;
}

/* System initialize for atomic system: '<Root>/VandD' */
static void VandD_Init(DW_VandD *localDW, P_VandD *localP, X_VandD *localX)
{
  static const uint32_T tmp[625] = { 5489U, 1301868182U, 2938499221U,
    2950281878U, 1875628136U, 751856242U, 944701696U, 2243192071U, 694061057U,
    219885934U, 2066767472U, 3182869408U, 485472502U, 2336857883U, 1071588843U,
    3418470598U, 951210697U, 3693558366U, 2923482051U, 1793174584U, 2982310801U,
    1586906132U, 1951078751U, 1808158765U, 1733897588U, 431328322U, 4202539044U,
    530658942U, 1714810322U, 3025256284U, 3342585396U, 1937033938U, 2640572511U,
    1654299090U, 3692403553U, 4233871309U, 3497650794U, 862629010U, 2943236032U,
    2426458545U, 1603307207U, 1133453895U, 3099196360U, 2208657629U, 2747653927U,
    931059398U, 761573964U, 3157853227U, 785880413U, 730313442U, 124945756U,
    2937117055U, 3295982469U, 1724353043U, 3021675344U, 3884886417U, 4010150098U,
    4056961966U, 699635835U, 2681338818U, 1339167484U, 720757518U, 2800161476U,
    2376097373U, 1532957371U, 3902664099U, 1238982754U, 3725394514U, 3449176889U,
    3570962471U, 4287636090U, 4087307012U, 3603343627U, 202242161U, 2995682783U,
    1620962684U, 3704723357U, 371613603U, 2814834333U, 2111005706U, 624778151U,
    2094172212U, 4284947003U, 1211977835U, 991917094U, 1570449747U, 2962370480U,
    1259410321U, 170182696U, 146300961U, 2836829791U, 619452428U, 2723670296U,
    1881399711U, 1161269684U, 1675188680U, 4132175277U, 780088327U, 3409462821U,
    1036518241U, 1834958505U, 3048448173U, 161811569U, 618488316U, 44795092U,
    3918322701U, 1924681712U, 3239478144U, 383254043U, 4042306580U, 2146983041U,
    3992780527U, 3518029708U, 3545545436U, 3901231469U, 1896136409U, 2028528556U,
    2339662006U, 501326714U, 2060962201U, 2502746480U, 561575027U, 581893337U,
    3393774360U, 1778912547U, 3626131687U, 2175155826U, 319853231U, 986875531U,
    819755096U, 2915734330U, 2688355739U, 3482074849U, 2736559U, 2296975761U,
    1029741190U, 2876812646U, 690154749U, 579200347U, 4027461746U, 1285330465U,
    2701024045U, 4117700889U, 759495121U, 3332270341U, 2313004527U, 2277067795U,
    4131855432U, 2722057515U, 1264804546U, 3848622725U, 2211267957U, 4100593547U,
    959123777U, 2130745407U, 3194437393U, 486673947U, 1377371204U, 17472727U,
    352317554U, 3955548058U, 159652094U, 1232063192U, 3835177280U, 49423123U,
    3083993636U, 733092U, 2120519771U, 2573409834U, 1112952433U, 3239502554U,
    761045320U, 1087580692U, 2540165110U, 641058802U, 1792435497U, 2261799288U,
    1579184083U, 627146892U, 2165744623U, 2200142389U, 2167590760U, 2381418376U,
    1793358889U, 3081659520U, 1663384067U, 2009658756U, 2689600308U, 739136266U,
    2304581039U, 3529067263U, 591360555U, 525209271U, 3131882996U, 294230224U,
    2076220115U, 3113580446U, 1245621585U, 1386885462U, 3203270426U, 123512128U,
    12350217U, 354956375U, 4282398238U, 3356876605U, 3888857667U, 157639694U,
    2616064085U, 1563068963U, 2762125883U, 4045394511U, 4180452559U, 3294769488U,
    1684529556U, 1002945951U, 3181438866U, 22506664U, 691783457U, 2685221343U,
    171579916U, 3878728600U, 2475806724U, 2030324028U, 3331164912U, 1708711359U,
    1970023127U, 2859691344U, 2588476477U, 2748146879U, 136111222U, 2967685492U,
    909517429U, 2835297809U, 3206906216U, 3186870716U, 341264097U, 2542035121U,
    3353277068U, 548223577U, 3170936588U, 1678403446U, 297435620U, 2337555430U,
    466603495U, 1132321815U, 1208589219U, 696392160U, 894244439U, 2562678859U,
    470224582U, 3306867480U, 201364898U, 2075966438U, 1767227936U, 2929737987U,
    3674877796U, 2654196643U, 3692734598U, 3528895099U, 2796780123U, 3048728353U,
    842329300U, 191554730U, 2922459673U, 3489020079U, 3979110629U, 1022523848U,
    2202932467U, 3583655201U, 3565113719U, 587085778U, 4176046313U, 3013713762U,
    950944241U, 396426791U, 3784844662U, 3477431613U, 3594592395U, 2782043838U,
    3392093507U, 3106564952U, 2829419931U, 1358665591U, 2206918825U, 3170783123U,
    31522386U, 2988194168U, 1782249537U, 1105080928U, 843500134U, 1225290080U,
    1521001832U, 3605886097U, 2802786495U, 2728923319U, 3996284304U, 903417639U,
    1171249804U, 1020374987U, 2824535874U, 423621996U, 1988534473U, 2493544470U,
    1008604435U, 1756003503U, 1488867287U, 1386808992U, 732088248U, 1780630732U,
    2482101014U, 976561178U, 1543448953U, 2602866064U, 2021139923U, 1952599828U,
    2360242564U, 2117959962U, 2753061860U, 2388623612U, 4138193781U, 2962920654U,
    2284970429U, 766920861U, 3457264692U, 2879611383U, 815055854U, 2332929068U,
    1254853997U, 3740375268U, 3799380844U, 4091048725U, 2006331129U, 1982546212U,
    686850534U, 1907447564U, 2682801776U, 2780821066U, 998290361U, 1342433871U,
    4195430425U, 607905174U, 3902331779U, 2454067926U, 1708133115U, 1170874362U,
    2008609376U, 3260320415U, 2211196135U, 433538229U, 2728786374U, 2189520818U,
    262554063U, 1182318347U, 3710237267U, 1221022450U, 715966018U, 2417068910U,
    2591870721U, 2870691989U, 3418190842U, 4238214053U, 1540704231U, 1575580968U,
    2095917976U, 4078310857U, 2313532447U, 2110690783U, 4056346629U, 4061784526U,
    1123218514U, 551538993U, 597148360U, 4120175196U, 3581618160U, 3181170517U,
    422862282U, 3227524138U, 1713114790U, 662317149U, 1230418732U, 928171837U,
    1324564878U, 1928816105U, 1786535431U, 2878099422U, 3290185549U, 539474248U,
    1657512683U, 552370646U, 1671741683U, 3655312128U, 1552739510U, 2605208763U,
    1441755014U, 181878989U, 3124053868U, 1447103986U, 3183906156U, 1728556020U,
    3502241336U, 3055466967U, 1013272474U, 818402132U, 1715099063U, 2900113506U,
    397254517U, 4194863039U, 1009068739U, 232864647U, 2540223708U, 2608288560U,
    2415367765U, 478404847U, 3455100648U, 3182600021U, 2115988978U, 434269567U,
    4117179324U, 3461774077U, 887256537U, 3545801025U, 286388911U, 3451742129U,
    1981164769U, 786667016U, 3310123729U, 3097811076U, 2224235657U, 2959658883U,
    3370969234U, 2514770915U, 3345656436U, 2677010851U, 2206236470U, 271648054U,
    2342188545U, 4292848611U, 3646533909U, 3754009956U, 3803931226U, 4160647125U,
    1477814055U, 4043852216U, 1876372354U, 3133294443U, 3871104810U, 3177020907U,
    2074304428U, 3479393793U, 759562891U, 164128153U, 1839069216U, 2114162633U,
    3989947309U, 3611054956U, 1333547922U, 835429831U, 494987340U, 171987910U,
    1252001001U, 370809172U, 3508925425U, 2535703112U, 1276855041U, 1922855120U,
    835673414U, 3030664304U, 613287117U, 171219893U, 3423096126U, 3376881639U,
    2287770315U, 1658692645U, 1262815245U, 3957234326U, 1168096164U, 2968737525U,
    2655813712U, 2132313144U, 3976047964U, 326516571U, 353088456U, 3679188938U,
    3205649712U, 2654036126U, 1249024881U, 880166166U, 691800469U, 2229503665U,
    1673458056U, 4032208375U, 1851778863U, 2563757330U, 376742205U, 1794655231U,
    340247333U, 1505873033U, 396524441U, 879666767U, 3335579166U, 3260764261U,
    3335999539U, 506221798U, 4214658741U, 975887814U, 2080536343U, 3360539560U,
    571586418U, 138896374U, 4234352651U, 2737620262U, 3928362291U, 1516365296U,
    38056726U, 3599462320U, 3585007266U, 3850961033U, 471667319U, 1536883193U,
    2310166751U, 1861637689U, 2530999841U, 4139843801U, 2710569485U, 827578615U,
    2012334720U, 2907369459U, 3029312804U, 2820112398U, 1965028045U, 35518606U,
    2478379033U, 643747771U, 1924139484U, 4123405127U, 3811735531U, 3429660832U,
    3285177704U, 1948416081U, 1311525291U, 1183517742U, 1739192232U, 3979815115U,
    2567840007U, 4116821529U, 213304419U, 4125718577U, 1473064925U, 2442436592U,
    1893310111U, 4195361916U, 3747569474U, 828465101U, 2991227658U, 750582866U,
    1205170309U, 1409813056U, 678418130U, 1171531016U, 3821236156U, 354504587U,
    4202874632U, 3882511497U, 1893248677U, 1903078632U, 26340130U, 2069166240U,
    3657122492U, 3725758099U, 831344905U, 811453383U, 3447711422U, 2434543565U,
    4166886888U, 3358210805U, 4142984013U, 2988152326U, 3527824853U, 982082992U,
    2809155763U, 190157081U, 3340214818U, 2365432395U, 2548636180U, 2894533366U,
    3474657421U, 2372634704U, 2845748389U, 43024175U, 2774226648U, 1987702864U,
    3186502468U, 453610222U, 4204736567U, 1392892630U, 2471323686U, 2470534280U,
    3541393095U, 4269885866U, 3909911300U, 759132955U, 1482612480U, 667715263U,
    1795580598U, 2337923983U, 3390586366U, 581426223U, 1515718634U, 476374295U,
    705213300U, 363062054U, 2084697697U, 2407503428U, 2292957699U, 2426213835U,
    2199989172U, 1987356470U, 4026755612U, 2147252133U, 270400031U, 1367820199U,
    2369854699U, 2844269403U, 79981964U, 624U };

  /* InitializeConditions for Integrator: '<S4>/Integrator20' */
  localX->Integrator20_CSTATE = localP->Integrator20_IC;

  /* InitializeConditions for Integrator: '<S4>/Integrator21' */
  localX->Integrator21_CSTATE = localP->Integrator21_IC;

  /* InitializeConditions for Integrator: '<S4>/Integrator22' */
  localX->Integrator22_CSTATE = localP->Integrator22_IC;

  /* InitializeConditions for Integrator: '<S4>/Integrator' */
  localX->Integrator_CSTATE = rtP.vxo;

  /* InitializeConditions for Integrator: '<S4>/Integrator2' */
  localX->Integrator2_CSTATE = rtP.ro;

  /* InitializeConditions for Integrator: '<S4>/Integrator1' */
  localX->Integrator1_CSTATE = rtP.vyo;

  /* InitializeConditions for Integrator: '<S4>/Integrator15' */
  localX->Integrator15_CSTATE = localP->Integrator15_IC;

  /* InitializeConditions for Integrator: '<S4>/Integrator14' */
  localX->Integrator14_CSTATE = localP->Integrator14_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  localDW->TimeStampA = (rtInf);
  localDW->TimeStampB = (rtInf);

  /* InitializeConditions for Derivative: '<S4>/Derivative1' */
  localDW->TimeStampA_l = (rtInf);
  localDW->TimeStampB_o = (rtInf);

  /* InitializeConditions for Derivative: '<S4>/Derivative2' */
  localDW->TimeStampA_p = (rtInf);
  localDW->TimeStampB_n = (rtInf);

  /* InitializeConditions for Derivative: '<S4>/Derivative3' */
  localDW->TimeStampA_o = (rtInf);
  localDW->TimeStampB_l = (rtInf);

  /* InitializeConditions for Integrator: '<S4>/Integrator11' */
  localX->Integrator11_CSTATE = rtP.zu3_0;

  /* InitializeConditions for Integrator: '<S4>/Integrator10' */
  localX->Integrator10_CSTATE = localP->Integrator10_IC;

  /* InitializeConditions for Integrator: '<S4>/Integrator13' */
  localX->Integrator13_CSTATE = rtP.zu4_0;

  /* InitializeConditions for Integrator: '<S4>/Integrator12' */
  localX->Integrator12_CSTATE = localP->Integrator12_IC;

  /* InitializeConditions for Integrator: '<S4>/Integrator7' */
  localX->Integrator7_CSTATE = rtP.zu1_0;

  /* InitializeConditions for Integrator: '<S4>/Integrator3' */
  localX->Integrator3_CSTATE = localP->Integrator3_IC;

  /* InitializeConditions for Integrator: '<S4>/Integrator9' */
  localX->Integrator9_CSTATE = rtP.zu2_0;

  /* InitializeConditions for Integrator: '<S4>/Integrator8' */
  localX->Integrator8_CSTATE = localP->Integrator8_IC;

  /* InitializeConditions for Integrator: '<S48>/Integrator' */
  localX->Integrator_CSTATE_p =
    localP->PIDController_InitialConditionForIntegrator;

  /* InitializeConditions for Integrator: '<S43>/Filter' */
  localX->Filter_CSTATE = localP->PIDController_InitialConditionForFilter;

  /* InitializeConditions for Integrator: '<S96>/Integrator' */
  localX->Integrator_CSTATE_b =
    localP->PIDController1_InitialConditionForIntegrator;

  /* InitializeConditions for Integrator: '<S91>/Filter' */
  localX->Filter_CSTATE_n = localP->PIDController1_InitialConditionForFilter;

  /* InitializeConditions for Integrator: '<S4>/Integrator19' */
  localX->Integrator19_CSTATE = localP->Integrator19_IC;

  /* InitializeConditions for Integrator: '<S4>/Integrator17' */
  localX->Integrator17_CSTATE = localP->Integrator17_IC;

  /* InitializeConditions for Integrator: '<S4>/Integrator16' */
  localX->Integrator16_CSTATE = localP->Integrator16_IC;

  /* InitializeConditions for Integrator: '<S4>/Integrator18' */
  localX->Integrator18_CSTATE = localP->Integrator18_IC;

  /* InitializeConditions for Integrator: '<S4>/Integrator5' */
  localX->Integrator5_CSTATE = rtP.Xo;

  /* InitializeConditions for Integrator: '<S4>/Integrator6' */
  localX->Integrator6_CSTATE = rtP.Yo;

  /* InitializeConditions for TransferFcn: '<S3>/Transfer Fcn' */
  localX->TransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Memory: '<S113>/Memory' */
  localDW->Memory_PreviousInput = localP->Memory_InitialCondition;

  /* InitializeConditions for Memory: '<S113>/Memory1' */
  localDW->Memory1_PreviousInput = localP->Memory1_InitialCondition;

  /* InitializeConditions for Integrator: '<S4>/Integrator4' */
  localX->Integrator4_CSTATE = rtP.THETAo;

  /* InitializeConditions for Integrator: '<S4>/Integrator23' */
  localX->Integrator23_CSTATE = rtP.Yo;

  /* SystemInitialize for MATLAB Function: '<S113>/MATLAB Function1' */
  memcpy(&localDW->state_d[0], &tmp[0], 625U * sizeof(uint32_T));
  localDW->method_b = 7U;
  localDW->state_g = 1144108930U;
  localDW->state[0] = 362436069U;
  localDW->state_o[0] = 362436069U;
  localDW->state[1] = 521288629U;
  localDW->state_o[1] = 521288629U;
}

/* Outputs for atomic system: '<Root>/VandD' */
static void VandD_g(RT_MODEL * const rtM, const real_T rtu_lamda[4], real_T
                    rtyyyy_Fy1Fy2Fy3Fy4[4], real_T *rty_Steeringwheelangle,
                    DW_VandD *localDW, P_VandD *localP, X_VandD *localX)
{
  real_T *lastU;
  boolean_T rEQ0;
  static const real_T a[5] = { 1.8, 3.6, 5.3999999999999995, 7.2, 9.0 };

  static const real_T g[5] = { 0.2, 0.6, 1.0, 0.8, 0.4 };

  /* MATLAB Function: '<S4>/MATLAB Function3' incorporates:
   *  Integrator: '<S4>/Integrator20'
   *  Integrator: '<S4>/Integrator21'
   *  Integrator: '<S4>/Integrator22'
   */
  localDW->absx11 = cos(localX->Integrator22_CSTATE);
  localDW->absx21 = sin(localX->Integrator20_CSTATE);
  localDW->absx31 = sin(localX->Integrator22_CSTATE);
  localDW->Derivative2 = cos(localX->Integrator20_CSTATE);
  localDW->Derivative3 = sin(localX->Integrator21_CSTATE);
  localDW->rtb_alpha_idx_0 = cos(localX->Integrator21_CSTATE);
  localDW->x[0] = localDW->rtb_alpha_idx_0 * localDW->absx11;
  localDW->x[3] = localDW->absx21 * localDW->Derivative3 * localDW->absx11 +
    localDW->Derivative2 * localDW->absx31;
  localDW->x[6] = localDW->absx21 * localDW->absx31 - localDW->Derivative2 *
    localDW->Derivative3 * localDW->absx11;
  localDW->x[1] = -localDW->rtb_alpha_idx_0 * localDW->absx31;
  localDW->x[4] = localDW->Derivative2 * localDW->absx11 - sin
    (localX->Integrator20_CSTATE) * sin(localX->Integrator21_CSTATE) *
    localDW->absx31;
  localDW->x[7] = cos(localX->Integrator20_CSTATE) * sin
    (localX->Integrator21_CSTATE) * localDW->absx31 + localDW->absx21 *
    localDW->absx11;
  localDW->x[2] = localDW->Derivative3;
  localDW->x[5] = -localDW->absx21 * localDW->rtb_alpha_idx_0;
  localDW->x[8] = localDW->Derivative2 * localDW->rtb_alpha_idx_0;
  memcpy(&localDW->b_x_k[0], &localDW->x[0], 9U * sizeof(real_T));
  localDW->p1 = 0;
  localDW->p2 = 3;
  localDW->p3 = 6;
  localDW->absx11 = fabs(localDW->x[0]);
  localDW->absx21 = fabs(localDW->x[1]);
  localDW->absx31 = fabs(localDW->Derivative3);
  if ((localDW->absx21 > localDW->absx11) && (localDW->absx21 > localDW->absx31))
  {
    localDW->p1 = 3;
    localDW->p2 = 0;
    localDW->b_x_k[0] = localDW->x[1];
    localDW->b_x_k[1] = localDW->x[0];
    localDW->b_x_k[3] = localDW->x[4];
    localDW->b_x_k[4] = localDW->x[3];
    localDW->b_x_k[6] = localDW->x[7];
    localDW->b_x_k[7] = localDW->x[6];
  } else if (localDW->absx31 > localDW->absx11) {
    localDW->p1 = 6;
    localDW->p3 = 0;
    localDW->b_x_k[0] = localDW->Derivative3;
    localDW->b_x_k[2] = localDW->x[0];
    localDW->b_x_k[3] = localDW->x[5];
    localDW->b_x_k[5] = localDW->x[3];
    localDW->b_x_k[6] = localDW->x[8];
    localDW->b_x_k[8] = localDW->x[6];
  }

  localDW->b_x_k[1] /= localDW->b_x_k[0];
  localDW->b_x_k[2] /= localDW->b_x_k[0];
  localDW->b_x_k[4] -= localDW->b_x_k[1] * localDW->b_x_k[3];
  localDW->b_x_k[5] -= localDW->b_x_k[2] * localDW->b_x_k[3];
  localDW->b_x_k[7] -= localDW->b_x_k[1] * localDW->b_x_k[6];
  localDW->b_x_k[8] -= localDW->b_x_k[2] * localDW->b_x_k[6];
  if (fabs(localDW->b_x_k[5]) > fabs(localDW->b_x_k[4])) {
    localDW->itmp = localDW->p2;
    localDW->p2 = localDW->p3;
    localDW->p3 = localDW->itmp;
    localDW->absx11 = localDW->b_x_k[1];
    localDW->b_x_k[1] = localDW->b_x_k[2];
    localDW->b_x_k[2] = localDW->absx11;
    localDW->absx11 = localDW->b_x_k[4];
    localDW->b_x_k[4] = localDW->b_x_k[5];
    localDW->b_x_k[5] = localDW->absx11;
    localDW->absx11 = localDW->b_x_k[7];
    localDW->b_x_k[7] = localDW->b_x_k[8];
    localDW->b_x_k[8] = localDW->absx11;
  }

  localDW->b_x_k[5] /= localDW->b_x_k[4];
  localDW->b_x_k[8] -= localDW->b_x_k[5] * localDW->b_x_k[7];
  localDW->absx11 = (localDW->b_x_k[1] * localDW->b_x_k[5] - localDW->b_x_k[2]) /
    localDW->b_x_k[8];
  localDW->absx21 = -(localDW->b_x_k[7] * localDW->absx11 + localDW->b_x_k[1]) /
    localDW->b_x_k[4];
  localDW->y[localDW->p1] = ((1.0 - localDW->b_x_k[3] * localDW->absx21) -
    localDW->b_x_k[6] * localDW->absx11) / localDW->b_x_k[0];
  localDW->y[localDW->p1 + 1] = localDW->absx21;
  localDW->y[localDW->p1 + 2] = localDW->absx11;
  localDW->absx11 = -localDW->b_x_k[5] / localDW->b_x_k[8];
  localDW->absx21 = (1.0 - localDW->b_x_k[7] * localDW->absx11) / localDW->
    b_x_k[4];
  localDW->y[localDW->p2] = -(localDW->b_x_k[3] * localDW->absx21 +
    localDW->b_x_k[6] * localDW->absx11) / localDW->b_x_k[0];
  localDW->y[localDW->p2 + 1] = localDW->absx21;
  localDW->y[localDW->p2 + 2] = localDW->absx11;
  localDW->absx11 = 1.0 / localDW->b_x_k[8];
  localDW->absx21 = -localDW->b_x_k[7] * localDW->absx11 / localDW->b_x_k[4];
  localDW->y[localDW->p3] = -(localDW->b_x_k[3] * localDW->absx21 +
    localDW->b_x_k[6] * localDW->absx11) / localDW->b_x_k[0];
  localDW->y[localDW->p3 + 1] = localDW->absx21;
  localDW->y[localDW->p3 + 2] = localDW->absx11;

  /* End of MATLAB Function: '<S4>/MATLAB Function3' */

  /* Integrator: '<S4>/Integrator' */
  localDW->vx = localX->Integrator_CSTATE;

  /* Integrator: '<S4>/Integrator2' */
  localDW->psi_dot = localX->Integrator2_CSTATE;

  /* Integrator: '<S4>/Integrator1' */
  localDW->vy = localX->Integrator1_CSTATE;

  /* MATLAB Function: '<S4>/MATLAB Function1' */
  localDW->absx11 = (rtP.l1 * localDW->psi_dot + localDW->vy) / localDW->vx;
  localDW->rtb_alpha_idx_0 = rtu_lamda[0] - localDW->absx11;
  localDW->rtb_alpha_idx_1 = rtu_lamda[1] - localDW->absx11;
  localDW->absx11 = (localDW->vy - rtP.l2 * localDW->psi_dot) / localDW->vx;
  localDW->rtb_alpha_idx_2 = rtu_lamda[2] - localDW->absx11;
  localDW->rtb_alpha_idx_3 = rtu_lamda[3] - localDW->absx11;

  /* Integrator: '<S4>/Integrator14' */
  localDW->zs_dot = localX->Integrator14_CSTATE;

  /* MATLAB Function: '<S4>/MATLAB Function' */
  for (localDW->p1 = 0; localDW->p1 < 3; localDW->p1++) {
    localDW->absx11 = localDW->y[localDW->p1] * rtP.l1;
    localDW->rtb_Fs_idx_0 = localDW->y[localDW->p1 + 3];
    localDW->rtb_Fs_idx_1 = localDW->y[localDW->p1 + 6] * 0.0;
    localDW->rtb_Fs_idx_2 = localDW->y[localDW->p1] * -rtP.l2;
    localDW->rtb_y_b[localDW->p1] = (localDW->rtb_Fs_idx_0 * -rtP.b1 +
      localDW->absx11) + localDW->rtb_Fs_idx_1;
    localDW->rtb_y_p[localDW->p1] = (localDW->rtb_Fs_idx_0 * rtP.b2 +
      localDW->absx11) + localDW->rtb_Fs_idx_1;
    localDW->rtb_y_c[localDW->p1] = (localDW->rtb_Fs_idx_0 * -rtP.b3 +
      localDW->rtb_Fs_idx_2) + localDW->rtb_Fs_idx_1;
    localDW->rtb_y_f[localDW->p1] = (localDW->rtb_Fs_idx_0 * rtP.b4 +
      localDW->rtb_Fs_idx_2) + localDW->rtb_Fs_idx_1;
  }

  localDW->Dz1 = localDW->rtb_y_b[2];
  localDW->Dz2 = localDW->rtb_y_p[2];
  localDW->Dz3 = localDW->rtb_y_c[2];
  localDW->Dz4 = localDW->rtb_y_f[2];

  /* End of MATLAB Function: '<S4>/MATLAB Function' */

  /* Derivative: '<S4>/Derivative' incorporates:
   *  Clock: '<S113>/Clock'
   *  Derivative: '<S4>/Derivative1'
   *  Derivative: '<S4>/Derivative2'
   *  Derivative: '<S4>/Derivative3'
   */
  localDW->absx11 = rtM->Timing.t[0];
  if ((localDW->TimeStampA >= localDW->absx11) && (localDW->TimeStampB >=
       localDW->absx11)) {
    localDW->absx21 = 0.0;
  } else {
    localDW->Derivative3 = localDW->TimeStampA;
    lastU = &localDW->LastUAtTimeA;
    if (localDW->TimeStampA < localDW->TimeStampB) {
      if (localDW->TimeStampB < localDW->absx11) {
        localDW->Derivative3 = localDW->TimeStampB;
        lastU = &localDW->LastUAtTimeB;
      }
    } else if (localDW->TimeStampA >= localDW->absx11) {
      localDW->Derivative3 = localDW->TimeStampB;
      lastU = &localDW->LastUAtTimeB;
    }

    localDW->absx21 = (localDW->Dz1 - *lastU) / (localDW->absx11 -
      localDW->Derivative3);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Derivative: '<S4>/Derivative1' */
  if ((localDW->TimeStampA_l >= localDW->absx11) && (localDW->TimeStampB_o >=
       localDW->absx11)) {
    localDW->absx31 = 0.0;
  } else {
    localDW->Derivative3 = localDW->TimeStampA_l;
    lastU = &localDW->LastUAtTimeA_k;
    if (localDW->TimeStampA_l < localDW->TimeStampB_o) {
      if (localDW->TimeStampB_o < localDW->absx11) {
        localDW->Derivative3 = localDW->TimeStampB_o;
        lastU = &localDW->LastUAtTimeB_a;
      }
    } else if (localDW->TimeStampA_l >= localDW->absx11) {
      localDW->Derivative3 = localDW->TimeStampB_o;
      lastU = &localDW->LastUAtTimeB_a;
    }

    localDW->absx31 = (localDW->Dz2 - *lastU) / (localDW->absx11 -
      localDW->Derivative3);
  }

  /* Derivative: '<S4>/Derivative2' */
  if ((localDW->TimeStampA_p >= localDW->absx11) && (localDW->TimeStampB_n >=
       localDW->absx11)) {
    localDW->Derivative2 = 0.0;
  } else {
    localDW->Derivative3 = localDW->TimeStampA_p;
    lastU = &localDW->LastUAtTimeA_n;
    if (localDW->TimeStampA_p < localDW->TimeStampB_n) {
      if (localDW->TimeStampB_n < localDW->absx11) {
        localDW->Derivative3 = localDW->TimeStampB_n;
        lastU = &localDW->LastUAtTimeB_e;
      }
    } else if (localDW->TimeStampA_p >= localDW->absx11) {
      localDW->Derivative3 = localDW->TimeStampB_n;
      lastU = &localDW->LastUAtTimeB_e;
    }

    localDW->Derivative2 = (localDW->Dz3 - *lastU) / (localDW->absx11 -
      localDW->Derivative3);
  }

  /* Derivative: '<S4>/Derivative3' */
  if ((localDW->TimeStampA_o >= localDW->absx11) && (localDW->TimeStampB_l >=
       localDW->absx11)) {
    localDW->Derivative3 = 0.0;
  } else {
    localDW->Derivative3 = localDW->TimeStampA_o;
    lastU = &localDW->LastUAtTimeA_e;
    if (localDW->TimeStampA_o < localDW->TimeStampB_l) {
      if (localDW->TimeStampB_l < localDW->absx11) {
        localDW->Derivative3 = localDW->TimeStampB_l;
        lastU = &localDW->LastUAtTimeB_h;
      }
    } else if (localDW->TimeStampA_o >= localDW->absx11) {
      localDW->Derivative3 = localDW->TimeStampB_l;
      lastU = &localDW->LastUAtTimeB_h;
    }

    localDW->Derivative3 = (localDW->Dz4 - *lastU) / (localDW->absx11 -
      localDW->Derivative3);
  }

  /* Integrator: '<S4>/Integrator10' */
  localDW->Integrator10 = localX->Integrator10_CSTATE;

  /* Integrator: '<S4>/Integrator12' */
  localDW->Integrator12 = localX->Integrator12_CSTATE;

  /* Integrator: '<S4>/Integrator3' */
  localDW->Integrator3 = localX->Integrator3_CSTATE;

  /* Integrator: '<S4>/Integrator8' */
  localDW->Integrator8 = localX->Integrator8_CSTATE;

  /* MATLAB Function: '<S4>/MATLAB Function4' incorporates:
   *  Integrator: '<S4>/Integrator11'
   *  Integrator: '<S4>/Integrator13'
   *  Integrator: '<S4>/Integrator15'
   *  Integrator: '<S4>/Integrator7'
   *  Integrator: '<S4>/Integrator9'
   */
  localDW->rtb_Fs_idx_0 = -(((localX->Integrator15_CSTATE - localDW->Dz1) -
    localX->Integrator7_CSTATE) * 265000.0 + ((localDW->zs_dot - localDW->absx21)
    - localDW->Integrator3) * 4500.0);
  localDW->rtb_Fs_idx_1 = -(((localX->Integrator15_CSTATE - localDW->Dz2) -
    localX->Integrator9_CSTATE) * 265000.0 + ((localDW->zs_dot - localDW->absx31)
    - localDW->Integrator8) * 4500.0);
  localDW->rtb_Fs_idx_2 = -(((localX->Integrator15_CSTATE - localDW->Dz3) -
    localX->Integrator11_CSTATE) * 52500.0 + ((localDW->zs_dot -
    localDW->Derivative2) - localDW->Integrator10) * 11207.0);
  localDW->Derivative3 = -(((localX->Integrator15_CSTATE - localDW->Dz4) -
    localX->Integrator13_CSTATE) * 52500.0 + ((localDW->zs_dot -
    localDW->Derivative3) - localDW->Integrator12) * 11207.0);
  localDW->Ft[0] = -(200000.0 * localX->Integrator7_CSTATE);
  localDW->Ft[1] = -(200000.0 * localX->Integrator9_CSTATE);
  localDW->Ft[2] = -(200000.0 * localX->Integrator11_CSTATE);
  localDW->Ft[3] = -(200000.0 * localX->Integrator13_CSTATE);

  /* MATLAB Function: '<S4>/MATLAB Function2' */
  localDW->Fz = localDW->Ft[0] / 1000.0;
  localDW->absx21 = localDW->Fz * localDW->Fz * -22.1 + 1011.0 * localDW->Fz;
  localDW->Derivative2 = sin(atan(0.208 * localDW->Fz) * 1.82) * 1078.0;
  localDW->absx31 = localDW->Derivative2 / (1.3 * localDW->absx21);
  localDW->Fz_idx_0 = localDW->Fz;
  localDW->D_idx_0 = localDW->absx21;
  localDW->BCD_idx_0 = localDW->Derivative2;
  localDW->Fz = localDW->Ft[1] / 1000.0;
  localDW->absx21 = localDW->Fz * localDW->Fz * -22.1 + 1011.0 * localDW->Fz;
  localDW->Derivative2 = sin(atan(0.208 * localDW->Fz) * 1.82) * 1078.0;
  localDW->B_idx_1 = localDW->Derivative2 / (1.3 * localDW->absx21);
  localDW->Fz_idx_1 = localDW->Fz;
  localDW->D_idx_1 = localDW->absx21;
  localDW->BCD_idx_1 = localDW->Derivative2;
  localDW->Fz = localDW->Ft[2] / 1000.0;
  localDW->absx21 = localDW->Fz * localDW->Fz * -22.1 + 1011.0 * localDW->Fz;
  localDW->Derivative2 = sin(atan(0.208 * localDW->Fz) * 1.82) * 1078.0;
  localDW->B_idx_2 = localDW->Derivative2 / (1.3 * localDW->absx21);
  localDW->Fz_idx_2 = localDW->Fz;
  localDW->D_idx_2 = localDW->absx21;
  localDW->BCD_idx_2 = localDW->Derivative2;
  localDW->Fz = localDW->Ft[3] / 1000.0;
  localDW->D_tmp = localDW->Fz * localDW->Fz;
  localDW->absx21 = localDW->D_tmp * -22.1 + 1011.0 * localDW->Fz;
  localDW->Derivative2 = sin(atan(0.208 * localDW->Fz) * 1.82) * 1078.0;
  localDW->B_idx_3 = localDW->Derivative2 / (1.3 * localDW->absx21);
  localDW->BCD_idx_0 = (localDW->BCD_idx_0 + localDW->BCD_idx_1) * 57.3;
  localDW->Fz_idx_0 = (localDW->Fz_idx_0 * localDW->Fz_idx_0 * 0.0 + -0.354 *
                       localDW->Fz_idx_0) + 0.707;
  localDW->Fz_idx_1 = (localDW->Fz_idx_1 * localDW->Fz_idx_1 * 0.0 + -0.354 *
                       localDW->Fz_idx_1) + 0.707;
  localDW->Fz_idx_2 = (localDW->Fz_idx_2 * localDW->Fz_idx_2 * 0.0 + -0.354 *
                       localDW->Fz_idx_2) + 0.707;
  localDW->Fz = (localDW->D_tmp * 0.0 + -0.354 * localDW->Fz) + 0.707;
  if (rtP.Model == 1.0) {
    rtyyyy_Fy1Fy2Fy3Fy4[0] = localDW->BCD_idx_0 * localDW->rtb_alpha_idx_0;
    rtyyyy_Fy1Fy2Fy3Fy4[1] = localDW->BCD_idx_0 * localDW->rtb_alpha_idx_1;
    localDW->rtb_alpha_idx_0 = (localDW->BCD_idx_2 + localDW->Derivative2) *
      57.3;
    rtyyyy_Fy1Fy2Fy3Fy4[2] = localDW->rtb_alpha_idx_0 * localDW->rtb_alpha_idx_2;
    rtyyyy_Fy1Fy2Fy3Fy4[3] = localDW->rtb_alpha_idx_0 * localDW->rtb_alpha_idx_3;
  } else {
    rtyyyy_Fy1Fy2Fy3Fy4[0] = sin(atan(((1.0 - localDW->Fz_idx_0) * 57.3 *
      localDW->rtb_alpha_idx_0 + atan(localDW->absx31 * 57.3 *
      localDW->rtb_alpha_idx_0) * (localDW->Fz_idx_0 / localDW->absx31)) *
      localDW->absx31) * 1.3) * localDW->D_idx_0;
    rtyyyy_Fy1Fy2Fy3Fy4[1] = sin(atan(((1.0 - localDW->Fz_idx_1) * 57.3 *
      localDW->rtb_alpha_idx_1 + atan(localDW->B_idx_1 * 57.3 *
      localDW->rtb_alpha_idx_1) * (localDW->Fz_idx_1 / localDW->B_idx_1)) *
      localDW->B_idx_1) * 1.3) * localDW->D_idx_1;
    rtyyyy_Fy1Fy2Fy3Fy4[2] = sin(atan(((1.0 - localDW->Fz_idx_2) * 57.3 *
      localDW->rtb_alpha_idx_2 + atan(localDW->B_idx_2 * 57.3 *
      localDW->rtb_alpha_idx_2) * (localDW->Fz_idx_2 / localDW->B_idx_2)) *
      localDW->B_idx_2) * 1.3) * localDW->D_idx_2;
    rtyyyy_Fy1Fy2Fy3Fy4[3] = sin(atan(((1.0 - localDW->Fz) * 57.3 *
      localDW->rtb_alpha_idx_3 + atan(localDW->B_idx_3 * 57.3 *
      localDW->rtb_alpha_idx_3) * (localDW->Fz / localDW->B_idx_3)) *
      localDW->B_idx_3) * 1.3) * localDW->absx21;
  }

  /* Sum: '<S6>/Sum' incorporates:
   *  Constant: '<S6>/Constant'
   */
  localDW->absx21 = localP->Constant_Value - localDW->vx;

  /* Product: '<S51>/NProd Out' incorporates:
   *  Constant: '<S6>/Constant1'
   *  Integrator: '<S43>/Filter'
   *  MATLAB Function: '<S6>/MATLAB Function'
   *  Product: '<S42>/DProd Out'
   *  Sum: '<S43>/SumD'
   */
  localDW->NProdOut = (0.05 * localDW->vx * localDW->absx21 -
                       localX->Filter_CSTATE) * localP->Constant1_Value;

  /* Sum: '<S57>/Sum' incorporates:
   *  Integrator: '<S48>/Integrator'
   *  MATLAB Function: '<S6>/MATLAB Function'
   *  Product: '<S53>/PProd Out'
   */
  localDW->absx31 = ((70.0 * localDW->vx + 5.0) * localDW->absx21 +
                     localX->Integrator_CSTATE_p) + localDW->NProdOut;

  /* Integrator: '<S91>/Filter' */
  *rty_Steeringwheelangle = localX->Filter_CSTATE_n;

  /* Gain: '<S99>/Filter Coefficient' incorporates:
   *  Gain: '<S90>/Derivative Gain'
   *  Sum: '<S91>/SumD'
   */
  localDW->FilterCoefficient = (localP->PIDController1_D * localDW->absx21 -
    *rty_Steeringwheelangle) * localP->PIDController1_N;

  /* Sum: '<S105>/Sum' incorporates:
   *  Gain: '<S101>/Proportional Gain'
   *  Integrator: '<S96>/Integrator'
   */
  localDW->Derivative2 = (localP->PIDController1_P * localDW->absx21 +
    localX->Integrator_CSTATE_b) + localDW->FilterCoefficient;

  /* ManualSwitch: '<S6>/Manual Switch' */
  if (localP->ManualSwitch_CurrentSetting == 1) {
    /* ManualSwitch: '<S6>/Manual Switch' */
    localDW->ManualSwitch[0] = localDW->absx31;
    localDW->ManualSwitch[1] = localDW->absx31;
    localDW->ManualSwitch[2] = localDW->absx31;
    localDW->ManualSwitch[3] = localDW->absx31;
  } else {
    /* ManualSwitch: '<S6>/Manual Switch' */
    localDW->ManualSwitch[0] = localDW->Derivative2;
    localDW->ManualSwitch[1] = localDW->Derivative2;
    localDW->ManualSwitch[2] = localDW->Derivative2;
    localDW->ManualSwitch[3] = localDW->Derivative2;
  }

  /* End of ManualSwitch: '<S6>/Manual Switch' */
  for (localDW->p1 = 0; localDW->p1 < 3; localDW->p1++) {
    /* Product: '<S4>/Product9' incorporates:
     *  SignalConversion generated from: '<S4>/Vector Concatenate'
     * */
    localDW->Product9[localDW->p1] = 0.0;
    localDW->Product9[localDW->p1] += localDW->y[localDW->p1] * localDW->vx;
    localDW->Product9[localDW->p1] += localDW->y[localDW->p1 + 3] * localDW->vy;
    localDW->Product9[localDW->p1] += localDW->y[localDW->p1 + 6] *
      localDW->zs_dot;
  }

  /* Integrator: '<S4>/Integrator19' */
  localDW->Integrator19 = localX->Integrator19_CSTATE;

  /* Integrator: '<S4>/Integrator16' */
  localDW->Integrator16 = localX->Integrator16_CSTATE;

  /* MATLAB Function: '<S4>/MATLAB Function5' incorporates:
   *  Integrator: '<S4>/Integrator21'
   *  Integrator: '<S4>/Integrator22'
   */
  localDW->rtb_alpha_idx_0 = sin(localX->Integrator22_CSTATE) *
    localDW->Integrator16;
  localDW->rtb_alpha_idx_2 = cos(localX->Integrator22_CSTATE) *
    localDW->Integrator19;
  localDW->euler_dot1 = (localDW->rtb_alpha_idx_2 - localDW->rtb_alpha_idx_0) /
    cos(localX->Integrator21_CSTATE);
  localDW->euler_dot2 = sin(localX->Integrator22_CSTATE) * localDW->Integrator19
    + cos(localX->Integrator22_CSTATE) * localDW->Integrator16;
  localDW->euler_dot3 = (localDW->rtb_alpha_idx_0 - localDW->rtb_alpha_idx_2) *
    tan(localDW->Integrator16) + localDW->psi_dot;

  /* Product: '<S4>/Product14' incorporates:
   *  Constant: '<S4>/[hzpc,hzpc,hzpc,hzpc]'
   */
  *rty_Steeringwheelangle = ((localP->hzpchzpchzpchzpc_Value[0] *
    localDW->ManualSwitch[0] + localP->hzpchzpchzpchzpc_Value[1] *
    localDW->ManualSwitch[1]) + localP->hzpchzpchzpchzpc_Value[2] *
    localDW->ManualSwitch[2]) + localP->hzpchzpchzpchzpc_Value[3] *
    localDW->ManualSwitch[3];

  /* Gain: '<S4>/Gain6' incorporates:
   *  Gain: '<S4>/Gain2'
   *  Product: '<S4>/Product'
   *  Product: '<S4>/Product17'
   *  Sum: '<S4>/Sum'
   *  Sum: '<S4>/Sum of Elements1'
   */
  localDW->Gain6 = (((((localDW->ManualSwitch[0] + localDW->ManualSwitch[1]) +
                       localDW->ManualSwitch[2]) + localDW->ManualSwitch[3]) -
                     (((rtu_lamda[0] * rtyyyy_Fy1Fy2Fy3Fy4[0] + rtu_lamda[1] *
                        rtyyyy_Fy1Fy2Fy3Fy4[1]) + rtu_lamda[2] *
                       rtyyyy_Fy1Fy2Fy3Fy4[2]) + rtu_lamda[3] *
                      rtyyyy_Fy1Fy2Fy3Fy4[3])) + localDW->vy * localDW->psi_dot *
                    rtP.m) * (rtP.sw_x / rtP.m);

  /* Gain: '<S4>/Gain11' incorporates:
   *  Constant: '<S4>/Iz-Ix'
   *  Constant: '<S4>/[l1,l1,-l2,-l2]'
   *  Constant: '<S4>/dpc'
   *  Constant: '<S4>/m*g*dpc'
   *  Gain: '<S4>/Gain8'
   *  Integrator: '<S4>/Integrator17'
   *  Product: '<S4>/Product13'
   *  Product: '<S4>/Product15'
   *  Product: '<S4>/Product18'
   *  Product: '<S4>/wywz2'
   *  Product: '<S4>/wzwx'
   *  Sum: '<S4>/Sum8'
   */
  localDW->Gain11 = (((((((localP->l1l1l2l2_Value[0] * localDW->rtb_Fs_idx_0 +
    localP->l1l1l2l2_Value[1] * localDW->rtb_Fs_idx_1) + localP->l1l1l2l2_Value
    [2] * localDW->rtb_Fs_idx_2) + localP->l1l1l2l2_Value[3] *
    localDW->Derivative3) + *rty_Steeringwheelangle) + rtP.m * localDW->Gain6 *
                       rtP.drc) + localDW->psi_dot * localDW->Integrator19 *
                      rtP.Iz_Ix) + rtP.m * rtP.g * rtP.dpc *
                     localX->Integrator17_CSTATE) * (rtP.sw_theta / rtP.Iy);

  /* Product: '<S4>/Product2' incorporates:
   *  Constant: '<S4>/[b1, -b2, b3, -b4]'
   */
  *rty_Steeringwheelangle = ((localP->b1b2b3b4_Value_f[0] *
    localDW->ManualSwitch[0] + localP->b1b2b3b4_Value_f[1] *
    localDW->ManualSwitch[1]) + localP->b1b2b3b4_Value_f[2] *
    localDW->ManualSwitch[2]) + localP->b1b2b3b4_Value_f[3] *
    localDW->ManualSwitch[3];

  /* Gain: '<S4>/Gain3' incorporates:
   *  Constant: '<S4>/Ix-Iy'
   *  Constant: '<S4>/[-b1, b2, -b3, b4]'
   *  Constant: '<S4>/[l1, l1, -l2, -l2]'
   *  Product: '<S4>/Product3'
   *  Product: '<S4>/Product4'
   *  Product: '<S4>/Product5'
   *  Product: '<S4>/Product6'
   *  Product: '<S4>/Product7'
   *  Product: '<S4>/wxwy'
   *  Product: '<S4>/wywz3'
   *  Sum: '<S4>/Sum2'
   */
  localDW->Gain3 = (((((((localP->l1l1l2l2_Value_m[0] * rtyyyy_Fy1Fy2Fy3Fy4[0] +
    localP->l1l1l2l2_Value_m[1] * rtyyyy_Fy1Fy2Fy3Fy4[1]) +
    localP->l1l1l2l2_Value_m[2] * rtyyyy_Fy1Fy2Fy3Fy4[2]) +
                        localP->l1l1l2l2_Value_m[3] * rtyyyy_Fy1Fy2Fy3Fy4[3]) + *
                       rty_Steeringwheelangle) + (((localP->l1l1l2l2_Value_m[0] *
    rtu_lamda[0] * localDW->ManualSwitch[0] + localP->l1l1l2l2_Value_m[1] *
    rtu_lamda[1] * localDW->ManualSwitch[1]) + localP->l1l1l2l2_Value_m[2] *
    rtu_lamda[2] * localDW->ManualSwitch[2]) + localP->l1l1l2l2_Value_m[3] *
    rtu_lamda[3] * localDW->ManualSwitch[3])) + (((localP->b1b2b3b4_Value_j[0] *
    rtu_lamda[0] * rtyyyy_Fy1Fy2Fy3Fy4[0] + localP->b1b2b3b4_Value_j[1] *
    rtu_lamda[1] * rtyyyy_Fy1Fy2Fy3Fy4[1]) + localP->b1b2b3b4_Value_j[2] *
    rtu_lamda[2] * rtyyyy_Fy1Fy2Fy3Fy4[2]) + localP->b1b2b3b4_Value_j[3] *
    rtu_lamda[3] * rtyyyy_Fy1Fy2Fy3Fy4[3])) + localDW->Integrator19 *
                    localDW->psi_dot * rtP.Ix_Iy) * (rtP.sw_psi / rtP.Iz);

  /* Integrator: '<S4>/Integrator18' */
  *rty_Steeringwheelangle = localX->Integrator18_CSTATE;

  /* Product: '<S4>/Product16' incorporates:
   *  Constant: '<S4>/m*g*drc'
   */
  localDW->rtb_alpha_idx_1 = rtP.m * rtP.g * rtP.drc * *rty_Steeringwheelangle;

  /* Gain: '<S4>/Gain1' incorporates:
   *  Product: '<S4>/Product8'
   */
  *rty_Steeringwheelangle = localDW->vx * localDW->psi_dot * rtP.m;

  /* Gain: '<S4>/Gain12' incorporates:
   *  Constant: '<S4>/Iy-Iz'
   *  Constant: '<S4>/b1, -b2, b3, -b4'
   *  Constant: '<S4>/drc'
   *  Constant: '<S4>/hzrc'
   *  Product: '<S4>/Product10'
   *  Product: '<S4>/Product11'
   *  Product: '<S4>/Product12'
   *  Product: '<S4>/Product13'
   *  Product: '<S4>/wywz'
   *  Product: '<S4>/wywz1'
   *  Sum: '<S4>/Sum9'
   */
  localDW->Gain12 = (((((((localP->hzrc_Value[0] * rtyyyy_Fy1Fy2Fy3Fy4[0] +
    localP->hzrc_Value[1] * rtyyyy_Fy1Fy2Fy3Fy4[1]) + localP->hzrc_Value[2] *
    rtyyyy_Fy1Fy2Fy3Fy4[2]) + localP->hzrc_Value[3] * rtyyyy_Fy1Fy2Fy3Fy4[3]) +
                        localDW->rtb_alpha_idx_1) + *rty_Steeringwheelangle *
                       rtP.dpc) + (((localP->b1b2b3b4_Value[0] *
    localDW->rtb_Fs_idx_0 + localP->b1b2b3b4_Value[1] * localDW->rtb_Fs_idx_1) +
    localP->b1b2b3b4_Value[2] * localDW->rtb_Fs_idx_2) + localP->b1b2b3b4_Value
    [3] * localDW->Derivative3)) + localDW->psi_dot * localDW->Integrator16 *
                     rtP.Iy_Iz) * (rtP.sw_phi / rtP.Ix);

  /* Gain: '<S4>/Gain5' incorporates:
   *  Product: '<S4>/Product1'
   *  Sum: '<S4>/Sum of Elements'
   *  Sum: '<S4>/Sum1'
   */
  localDW->Gain5 = (((((rtu_lamda[0] * localDW->ManualSwitch[0] + rtu_lamda[1] *
                        localDW->ManualSwitch[1]) + rtu_lamda[2] *
                       localDW->ManualSwitch[2]) + rtu_lamda[3] *
                      localDW->ManualSwitch[3]) + (((rtyyyy_Fy1Fy2Fy3Fy4[0] +
    rtyyyy_Fy1Fy2Fy3Fy4[1]) + rtyyyy_Fy1Fy2Fy3Fy4[2]) + rtyyyy_Fy1Fy2Fy3Fy4[3]))
                    - *rty_Steeringwheelangle) * (rtP.sw_y / rtP.m);

  /* Product: '<S45>/IProd Out' incorporates:
   *  MATLAB Function: '<S6>/MATLAB Function'
   */
  localDW->IProdOut = 0.05 * localDW->vx * localDW->absx21;

  /* Gain: '<S93>/Integral Gain' */
  localDW->IntegralGain = localP->PIDController1_I * localDW->absx21;

  /* Integrator: '<S4>/Integrator5' */
  localDW->p = localX->Integrator5_CSTATE;

  /* Integrator: '<S4>/Integrator6' */
  localDW->p_o = localX->Integrator6_CSTATE;

  /* TransferFcn: '<S3>/Transfer Fcn' */
  localDW->TransferFcn = localP->TransferFcn_C * localX->TransferFcn_CSTATE;

  /* MATLAB Function: '<S5>/FWS Controller' incorporates:
   *  MATLAB Function: '<S4>/MATLAB Function2'
   *  SignalConversion generated from: '<S1>/Vector Concatenate2'
   * */
  localDW->rtb_alpha_idx_2 = rtP.radius;
  localDW->rtb_alpha_idx_3 = cos(localDW->TransferFcn);
  localDW->Derivative2 = sin(localDW->TransferFcn);
  if (rtP.radius == 0.0) {
    localDW->rtb_alpha_idx_2 = 1000.0;
  }

  localDW->Traj_loc[0].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[0];
  localDW->Traj_loc[1].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[1];
  localDW->Traj_loc[2].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[2];
  localDW->Traj_loc[3].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[3];
  localDW->Traj_loc[4].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[4];
  localDW->Traj_loc[5].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[5];
  localDW->Traj_loc[6].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[6];
  localDW->Traj_loc[7].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[7];
  localDW->Traj_loc[8].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[8];
  localDW->Traj_loc[9].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[9];
  localDW->Traj_loc[10].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[10];
  localDW->Traj_loc[11].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[11];
  localDW->Traj_loc[12].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[12];
  localDW->Traj_loc[13].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[13];
  localDW->Traj_loc[14].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[14];
  localDW->Traj_loc[15].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[15];
  localDW->Traj_loc[16].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[16];
  localDW->Traj_loc[17].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[17];
  localDW->Traj_loc[18].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[18];
  localDW->Traj_loc[19].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[19];
  localDW->Traj_loc[20].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[20];
  localDW->Traj_loc[21].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[21];
  localDW->Traj_loc[22].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[22];
  localDW->Traj_loc[23].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[23];
  localDW->Traj_loc[24].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[24];
  localDW->Traj_loc[25].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[25];
  localDW->Traj_loc[26].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[26];
  localDW->Traj_loc[27].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[27];
  localDW->Traj_loc[28].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[28];
  localDW->Traj_loc[29].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[29];
  localDW->Traj_loc[30].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[30];
  localDW->Traj_loc[31].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[31];
  localDW->Traj_loc[32].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[32];
  localDW->Traj_loc[33].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[33];
  localDW->Traj_loc[34].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[34];
  localDW->Traj_loc[35].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[35];
  localDW->Traj_loc[36].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[36];
  localDW->Traj_loc[37].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[37];
  localDW->Traj_loc[38].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[38];
  localDW->Traj_loc[39].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[39];
  localDW->Traj_loc[40].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[40];
  localDW->Traj_loc[41].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[41];
  localDW->Traj_loc[42].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[42];
  localDW->Traj_loc[43].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[43];
  localDW->Traj_loc[44].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[44];
  localDW->Traj_loc[45].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[45];
  localDW->Traj_loc[46].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[46];
  localDW->Traj_loc[47].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[47];
  localDW->Traj_loc[48].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[48];
  localDW->Traj_loc[49].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[49];
  localDW->Traj_loc[50].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[50];
  localDW->Traj_loc[51].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[51];
  localDW->Traj_loc[52].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[52];
  localDW->Traj_loc[53].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[53];
  localDW->Traj_loc[54].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[54];
  localDW->Traj_loc[55].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[55];
  localDW->Traj_loc[56].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[56];
  localDW->Traj_loc[57].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[57];
  localDW->Traj_loc[58].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[58];
  localDW->Traj_loc[59].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[59];
  localDW->Traj_loc[60].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[60];
  localDW->Traj_loc[61].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[61];
  localDW->Traj_loc[62].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[62];
  localDW->Traj_loc[63].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[63];
  localDW->Traj_loc[64].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[64];
  localDW->Traj_loc[65].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[65];
  localDW->Traj_loc[66].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[66];
  localDW->Traj_loc[67].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[67];
  localDW->Traj_loc[68].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[68];
  localDW->Traj_loc[69].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[69];
  localDW->Traj_loc[70].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[70];
  localDW->Traj_loc[71].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[71];
  localDW->Traj_loc[72].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[72];
  localDW->Traj_loc[73].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[73];
  localDW->Traj_loc[74].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[74];
  localDW->Traj_loc[75].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[75];
  localDW->Traj_loc[76].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[76];
  localDW->Traj_loc[77].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[77];
  localDW->Traj_loc[78].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[78];
  localDW->Traj_loc[79].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[79];
  localDW->Traj_loc[80].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[80];
  localDW->Traj_loc[81].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[81];
  localDW->Traj_loc[82].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[82];
  localDW->Traj_loc[83].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[83];
  localDW->Traj_loc[84].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[84];
  localDW->Traj_loc[85].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[85];
  localDW->Traj_loc[86].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[86];
  localDW->Traj_loc[87].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[87];
  localDW->Traj_loc[88].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[88];
  localDW->Traj_loc[89].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[89];
  localDW->Traj_loc[90].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[90];
  localDW->Traj_loc[91].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[91];
  localDW->Traj_loc[92].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[92];
  localDW->Traj_loc[93].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[93];
  localDW->Traj_loc[94].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[94];
  localDW->Traj_loc[95].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[95];
  localDW->Traj_loc[96].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[96];
  localDW->Traj_loc[97].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[97];
  localDW->Traj_loc[98].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[98];
  localDW->Traj_loc[99].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[99];
  localDW->Traj_loc[100].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[100];
  localDW->Traj_loc[101].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[101];
  localDW->Traj_loc[102].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[102];
  localDW->Traj_loc[103].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[103];
  localDW->Traj_loc[104].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[104];
  localDW->Traj_loc[105].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[105];
  localDW->Traj_loc[106].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[106];
  localDW->Traj_loc[107].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[107];
  localDW->Traj_loc[108].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[108];
  localDW->Traj_loc[109].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[109];
  localDW->Traj_loc[110].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[110];
  localDW->Traj_loc[111].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[111];
  localDW->Traj_loc[112].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[112];
  localDW->Traj_loc[113].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[113];
  localDW->Traj_loc[114].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[114];
  localDW->Traj_loc[115].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[115];
  localDW->Traj_loc[116].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[116];
  localDW->Traj_loc[117].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[117];
  localDW->Traj_loc[118].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[118];
  localDW->Traj_loc[119].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[119];
  localDW->Traj_loc[120].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[120];
  localDW->Traj_loc[121].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[121];
  localDW->Traj_loc[122].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[122];
  localDW->Traj_loc[123].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[123];
  localDW->Traj_loc[124].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[124];
  localDW->Traj_loc[125].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[125];
  localDW->Traj_loc[126].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[126];
  localDW->Traj_loc[127].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[127];
  localDW->Traj_loc[128].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[128];
  localDW->Traj_loc[129].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[129];
  localDW->Traj_loc[130].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[130];
  localDW->Traj_loc[131].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[131];
  localDW->Traj_loc[132].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[132];
  localDW->Traj_loc[133].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[133];
  localDW->Traj_loc[134].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[134];
  localDW->Traj_loc[135].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[135];
  localDW->Traj_loc[136].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[136];
  localDW->Traj_loc[137].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[137];
  localDW->Traj_loc[138].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[138];
  localDW->Traj_loc[139].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[139];
  localDW->Traj_loc[140].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[140];
  localDW->Traj_loc[141].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[141];
  localDW->Traj_loc[142].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[142];
  localDW->Traj_loc[143].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[143];
  localDW->Traj_loc[144].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[144];
  localDW->Traj_loc[145].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[145];
  localDW->Traj_loc[146].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[146];
  localDW->Traj_loc[147].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[147];
  localDW->Traj_loc[148].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[148];
  localDW->Traj_loc[149].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[149];
  localDW->Traj_loc[150].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[150];
  localDW->Traj_loc[151].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[151];
  localDW->Traj_loc[152].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[152];
  localDW->Traj_loc[153].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[153];
  localDW->Traj_loc[154].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[154];
  localDW->Traj_loc[155].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[155];
  localDW->Traj_loc[156].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[156];
  localDW->Traj_loc[157].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[157];
  localDW->Traj_loc[158].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[158];
  localDW->Traj_loc[159].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[159];
  localDW->Traj_loc[160].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[160];
  localDW->Traj_loc[161].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[161];
  localDW->Traj_loc[162].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[162];
  localDW->Traj_loc[163].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[163];
  localDW->Traj_loc[164].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[164];
  localDW->Traj_loc[165].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[165];
  localDW->Traj_loc[166].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[166];
  localDW->Traj_loc[167].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[167];
  localDW->Traj_loc[168].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[168];
  localDW->Traj_loc[169].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[169];
  localDW->Traj_loc[170].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[170];
  localDW->Traj_loc[171].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[171];
  localDW->Traj_loc[172].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[172];
  localDW->Traj_loc[173].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[173];
  localDW->Traj_loc[174].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[174];
  localDW->Traj_loc[175].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[175];
  localDW->Traj_loc[176].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[176];
  localDW->Traj_loc[177].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[177];
  localDW->Traj_loc[178].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[178];
  localDW->Traj_loc[179].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[179];
  localDW->Traj_loc[180].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[180];
  localDW->Traj_loc[181].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[181];
  localDW->Traj_loc[182].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[182];
  localDW->Traj_loc[183].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[183];
  localDW->Traj_loc[184].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[184];
  localDW->Traj_loc[185].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[185];
  localDW->Traj_loc[186].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[186];
  localDW->Traj_loc[187].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[187];
  localDW->Traj_loc[188].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[188];
  localDW->Traj_loc[189].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[189];
  localDW->Traj_loc[190].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[190];
  localDW->Traj_loc[191].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[191];
  localDW->Traj_loc[192].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[192];
  localDW->Traj_loc[193].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[193];
  localDW->Traj_loc[194].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[194];
  localDW->Traj_loc[195].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[195];
  localDW->Traj_loc[196].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[196];
  localDW->Traj_loc[197].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[197];
  localDW->Traj_loc[198].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[198];
  localDW->Traj_loc[199].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[199];
  localDW->Traj_loc[200].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[200];
  localDW->Traj_loc[201].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[201];
  localDW->Traj_loc[202].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[202];
  localDW->Traj_loc[203].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[203];
  localDW->Traj_loc[204].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[204];
  localDW->Traj_loc[205].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[205];
  localDW->Traj_loc[206].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[206];
  localDW->Traj_loc[207].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[207];
  localDW->Traj_loc[208].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[208];
  localDW->Traj_loc[209].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[209];
  localDW->Traj_loc[210].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[210];
  localDW->Traj_loc[211].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[211];
  localDW->Traj_loc[212].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[212];
  localDW->Traj_loc[213].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[213];
  localDW->Traj_loc[214].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[214];
  localDW->Traj_loc[215].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[215];
  localDW->Traj_loc[216].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[216];
  localDW->Traj_loc[217].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[217];
  localDW->Traj_loc[218].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[218];
  localDW->Traj_loc[219].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[219];
  localDW->Traj_loc[220].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[220];
  localDW->Traj_loc[221].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[221];
  localDW->Traj_loc[222].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[222];
  localDW->Traj_loc[223].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[223];
  localDW->Traj_loc[224].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[224];
  localDW->Traj_loc[225].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[225];
  localDW->Traj_loc[226].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[226];
  localDW->Traj_loc[227].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[227];
  localDW->Traj_loc[228].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[228];
  localDW->Traj_loc[229].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[229];
  localDW->Traj_loc[230].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[230];
  localDW->Traj_loc[231].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[231];
  localDW->Traj_loc[232].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[232];
  localDW->Traj_loc[233].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[233];
  localDW->Traj_loc[234].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[234];
  localDW->Traj_loc[235].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[235];
  localDW->Traj_loc[236].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[236];
  localDW->Traj_loc[237].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[237];
  localDW->Traj_loc[238].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[238];
  localDW->Traj_loc[239].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[239];
  localDW->Traj_loc[240].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[240];
  localDW->Traj_loc[241].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[241];
  localDW->Traj_loc[242].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[242];
  localDW->Traj_loc[243].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[243];
  localDW->Traj_loc[244].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[244];
  localDW->Traj_loc[245].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[245];
  localDW->Traj_loc[246].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[246];
  localDW->Traj_loc[247].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[247];
  localDW->Traj_loc[248].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[248];
  localDW->Traj_loc[249].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[249];
  localDW->Traj_loc[250].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[250];
  localDW->Traj_loc[251].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[251];
  localDW->Traj_loc[252].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[252];
  localDW->Traj_loc[253].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[253];
  localDW->Traj_loc[254].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[254];
  localDW->Traj_loc[255].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[255];
  localDW->Traj_loc[256].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[256];
  localDW->Traj_loc[257].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[257];
  localDW->Traj_loc[258].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[258];
  localDW->Traj_loc[259].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[259];
  localDW->Traj_loc[260].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[260];
  localDW->Traj_loc[261].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[261];
  localDW->Traj_loc[262].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[262];
  localDW->Traj_loc[263].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[263];
  localDW->Traj_loc[264].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[264];
  localDW->Traj_loc[265].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[265];
  localDW->Traj_loc[266].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[266];
  localDW->Traj_loc[267].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[267];
  localDW->Traj_loc[268].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[268];
  localDW->Traj_loc[269].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[269];
  localDW->Traj_loc[270].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[270];
  localDW->Traj_loc[271].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[271];
  localDW->Traj_loc[272].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[272];
  localDW->Traj_loc[273].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[273];
  localDW->Traj_loc[274].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[274];
  localDW->Traj_loc[275].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[275];
  localDW->Traj_loc[276].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[276];
  localDW->Traj_loc[277].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[277];
  localDW->Traj_loc[278].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[278];
  localDW->Traj_loc[279].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[279];
  localDW->Traj_loc[280].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[280];
  localDW->Traj_loc[281].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[281];
  localDW->Traj_loc[282].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[282];
  localDW->Traj_loc[283].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[283];
  localDW->Traj_loc[284].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[284];
  localDW->Traj_loc[285].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[285];
  localDW->Traj_loc[286].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[286];
  localDW->Traj_loc[287].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[287];
  localDW->Traj_loc[288].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[288];
  localDW->Traj_loc[289].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[289];
  localDW->Traj_loc[290].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[290];
  localDW->Traj_loc[291].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[291];
  localDW->Traj_loc[292].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[292];
  localDW->Traj_loc[293].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[293];
  localDW->Traj_loc[294].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[294];
  localDW->Traj_loc[295].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[295];
  localDW->Traj_loc[296].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[296];
  localDW->Traj_loc[297].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[297];
  localDW->Traj_loc[298].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[298];
  localDW->Traj_loc[299].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[299];
  localDW->Traj_loc[300].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[300];
  localDW->Traj_loc[301].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[301];
  localDW->Traj_loc[302].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[302];
  localDW->Traj_loc[303].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[303];
  localDW->Traj_loc[304].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[304];
  localDW->Traj_loc[305].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[305];
  localDW->Traj_loc[306].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[306];
  localDW->Traj_loc[307].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[307];
  localDW->Traj_loc[308].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[308];
  localDW->Traj_loc[309].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[309];
  localDW->Traj_loc[310].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[310];
  localDW->Traj_loc[311].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[311];
  localDW->Traj_loc[312].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[312];
  localDW->Traj_loc[313].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[313];
  localDW->Traj_loc[314].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[314];
  localDW->Traj_loc[315].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[315];
  localDW->Traj_loc[316].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[316];
  localDW->Traj_loc[317].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[317];
  localDW->Traj_loc[318].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[318];
  localDW->Traj_loc[319].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[319];
  localDW->Traj_loc[320].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[320];
  localDW->Traj_loc[321].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[321];
  localDW->Traj_loc[322].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[322];
  localDW->Traj_loc[323].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[323];
  localDW->Traj_loc[324].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[324];
  localDW->Traj_loc[325].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[325];
  localDW->Traj_loc[326].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[326];
  localDW->Traj_loc[327].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[327];
  localDW->Traj_loc[328].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[328];
  localDW->Traj_loc[329].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[329];
  localDW->Traj_loc[330].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[330];
  localDW->Traj_loc[331].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[331];
  localDW->Traj_loc[332].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[332];
  localDW->Traj_loc[333].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[333];
  localDW->Traj_loc[334].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[334];
  localDW->Traj_loc[335].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[335];
  localDW->Traj_loc[336].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[336];
  localDW->Traj_loc[337].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[337];
  localDW->Traj_loc[338].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[338];
  localDW->Traj_loc[339].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[339];
  localDW->Traj_loc[340].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[340];
  localDW->Traj_loc[341].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[341];
  localDW->Traj_loc[342].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[342];
  localDW->Traj_loc[343].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[343];
  localDW->Traj_loc[344].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[344];
  localDW->Traj_loc[345].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[345];
  localDW->Traj_loc[346].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[346];
  localDW->Traj_loc[347].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[347];
  localDW->Traj_loc[348].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[348];
  localDW->Traj_loc[349].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[349];
  localDW->Traj_loc[350].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[350];
  localDW->Traj_loc[351].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[351];
  localDW->Traj_loc[352].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[352];
  localDW->Traj_loc[353].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[353];
  localDW->Traj_loc[354].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[354];
  localDW->Traj_loc[355].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[355];
  localDW->Traj_loc[356].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[356];
  localDW->Traj_loc[357].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[357];
  localDW->Traj_loc[358].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[358];
  localDW->Traj_loc[359].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[359];
  localDW->Traj_loc[360].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[360];
  localDW->Traj_loc[361].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[361];
  localDW->Traj_loc[362].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[362];
  localDW->Traj_loc[363].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[363];
  localDW->Traj_loc[364].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[364];
  localDW->Traj_loc[365].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[365];
  localDW->Traj_loc[366].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[366];
  localDW->Traj_loc[367].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[367];
  localDW->Traj_loc[368].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[368];
  localDW->Traj_loc[369].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[369];
  localDW->Traj_loc[370].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[370];
  localDW->Traj_loc[371].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[371];
  localDW->Traj_loc[372].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[372];
  localDW->Traj_loc[373].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[373];
  localDW->Traj_loc[374].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[374];
  localDW->Traj_loc[375].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[375];
  localDW->Traj_loc[376].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[376];
  localDW->Traj_loc[377].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[377];
  localDW->Traj_loc[378].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[378];
  localDW->Traj_loc[379].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[379];
  localDW->Traj_loc[380].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[380];
  localDW->Traj_loc[381].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[381];
  localDW->Traj_loc[382].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[382];
  localDW->Traj_loc[383].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[383];
  localDW->Traj_loc[384].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[384];
  localDW->Traj_loc[385].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[385];
  localDW->Traj_loc[386].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[386];
  localDW->Traj_loc[387].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[387];
  localDW->Traj_loc[388].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[388];
  localDW->Traj_loc[389].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[389];
  localDW->Traj_loc[390].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[390];
  localDW->Traj_loc[391].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[391];
  localDW->Traj_loc[392].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[392];
  localDW->Traj_loc[393].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[393];
  localDW->Traj_loc[394].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[394];
  localDW->Traj_loc[395].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[395];
  localDW->Traj_loc[396].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[396];
  localDW->Traj_loc[397].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[397];
  localDW->Traj_loc[398].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[398];
  localDW->Traj_loc[399].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[399];
  localDW->Traj_loc[400].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[400];
  localDW->Traj_loc[401].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[401];
  localDW->Traj_loc[402].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[402];
  localDW->Traj_loc[403].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[403];
  localDW->Traj_loc[404].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[404];
  localDW->Traj_loc[405].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[405];
  localDW->Traj_loc[406].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[406];
  localDW->Traj_loc[407].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[407];
  localDW->Traj_loc[408].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[408];
  localDW->Traj_loc[409].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[409];
  localDW->Traj_loc[410].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[410];
  localDW->Traj_loc[411].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[411];
  localDW->Traj_loc[412].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[412];
  localDW->Traj_loc[413].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[413];
  localDW->Traj_loc[414].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[414];
  localDW->Traj_loc[415].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[415];
  localDW->Traj_loc[416].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[416];
  localDW->Traj_loc[417].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[417];
  localDW->Traj_loc[418].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[418];
  localDW->Traj_loc[419].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[419];
  localDW->Traj_loc[420].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[420];
  localDW->Traj_loc[421].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[421];
  localDW->Traj_loc[422].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[422];
  localDW->Traj_loc[423].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[423];
  localDW->Traj_loc[424].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[424];
  localDW->Traj_loc[425].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[425];
  localDW->Traj_loc[426].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[426];
  localDW->Traj_loc[427].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[427];
  localDW->Traj_loc[428].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[428];
  localDW->Traj_loc[429].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[429];
  localDW->Traj_loc[430].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[430];
  localDW->Traj_loc[431].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[431];
  localDW->Traj_loc[432].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[432];
  localDW->Traj_loc[433].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[433];
  localDW->Traj_loc[434].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[434];
  localDW->Traj_loc[435].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[435];
  localDW->Traj_loc[436].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[436];
  localDW->Traj_loc[437].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[437];
  localDW->Traj_loc[438].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[438];
  localDW->Traj_loc[439].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[439];
  localDW->Traj_loc[440].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[440];
  localDW->Traj_loc[441].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[441];
  localDW->Traj_loc[442].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[442];
  localDW->Traj_loc[443].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[443];
  localDW->Traj_loc[444].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[444];
  localDW->Traj_loc[445].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[445];
  localDW->Traj_loc[446].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[446];
  localDW->Traj_loc[447].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[447];
  localDW->Traj_loc[448].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[448];
  localDW->Traj_loc[449].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[449];
  localDW->Traj_loc[450].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[450];
  localDW->Traj_loc[451].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[451];
  localDW->Traj_loc[452].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[452];
  localDW->Traj_loc[453].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[453];
  localDW->Traj_loc[454].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[454];
  localDW->Traj_loc[455].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[455];
  localDW->Traj_loc[456].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[456];
  localDW->Traj_loc[457].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[457];
  localDW->Traj_loc[458].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[458];
  localDW->Traj_loc[459].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[459];
  localDW->Traj_loc[460].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[460];
  localDW->Traj_loc[461].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[461];
  localDW->Traj_loc[462].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[462];
  localDW->Traj_loc[463].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[463];
  localDW->Traj_loc[464].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[464];
  localDW->Traj_loc[465].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[465];
  localDW->Traj_loc[466].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[466];
  localDW->Traj_loc[467].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[467];
  localDW->Traj_loc[468].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[468];
  localDW->Traj_loc[469].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[469];
  localDW->Traj_loc[470].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[470];
  localDW->Traj_loc[471].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[471];
  localDW->Traj_loc[472].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[472];
  localDW->Traj_loc[473].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[473];
  localDW->Traj_loc[474].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[474];
  localDW->Traj_loc[475].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[475];
  localDW->Traj_loc[476].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[476];
  localDW->Traj_loc[477].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[477];
  localDW->Traj_loc[478].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[478];
  localDW->Traj_loc[479].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[479];
  localDW->Traj_loc[480].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[480];
  localDW->Traj_loc[481].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[481];
  localDW->Traj_loc[482].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[482];
  localDW->Traj_loc[483].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[483];
  localDW->Traj_loc[484].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[484];
  localDW->Traj_loc[485].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[485];
  localDW->Traj_loc[486].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[486];
  localDW->Traj_loc[487].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[487];
  localDW->Traj_loc[488].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[488];
  localDW->Traj_loc[489].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[489];
  localDW->Traj_loc[490].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[490];
  localDW->Traj_loc[491].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[491];
  localDW->Traj_loc[492].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[492];
  localDW->Traj_loc[493].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[493];
  localDW->Traj_loc[494].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[494];
  localDW->Traj_loc[495].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[495];
  localDW->Traj_loc[496].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[496];
  localDW->Traj_loc[497].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[497];
  localDW->Traj_loc[498].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[498];
  localDW->Traj_loc[499].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[499];
  localDW->Traj_loc[500].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[500];
  localDW->Traj_loc[501].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[501];
  localDW->Traj_loc[502].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[502];
  localDW->Traj_loc[503].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[503];
  localDW->Traj_loc[504].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[504];
  localDW->Traj_loc[505].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[505];
  localDW->Traj_loc[506].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[506];
  localDW->Traj_loc[507].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[507];
  localDW->Traj_loc[508].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[508];
  localDW->Traj_loc[509].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[509];
  localDW->Traj_loc[510].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[510];
  localDW->Traj_loc[511].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[511];
  localDW->Traj_loc[512].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[512];
  localDW->Traj_loc[513].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[513];
  localDW->Traj_loc[514].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[514];
  localDW->Traj_loc[515].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[515];
  localDW->Traj_loc[516].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[516];
  localDW->Traj_loc[517].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[517];
  localDW->Traj_loc[518].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[518];
  localDW->Traj_loc[519].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[519];
  localDW->Traj_loc[520].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[520];
  localDW->Traj_loc[521].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[521];
  localDW->Traj_loc[522].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[522];
  localDW->Traj_loc[523].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[523];
  localDW->Traj_loc[524].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[524];
  localDW->Traj_loc[525].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[525];
  localDW->Traj_loc[526].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[526];
  localDW->Traj_loc[527].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[527];
  localDW->Traj_loc[528].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[528];
  localDW->Traj_loc[529].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[529];
  localDW->Traj_loc[530].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[530];
  localDW->Traj_loc[531].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[531];
  localDW->Traj_loc[532].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[532];
  localDW->Traj_loc[533].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[533];
  localDW->Traj_loc[534].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[534];
  localDW->Traj_loc[535].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[535];
  localDW->Traj_loc[536].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[536];
  localDW->Traj_loc[537].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[537];
  localDW->Traj_loc[538].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[538];
  localDW->Traj_loc[539].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[539];
  localDW->Traj_loc[540].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[540];
  localDW->Traj_loc[541].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[541];
  localDW->Traj_loc[542].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[542];
  localDW->Traj_loc[543].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[543];
  localDW->Traj_loc[544].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[544];
  localDW->Traj_loc[545].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[545];
  localDW->Traj_loc[546].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[546];
  localDW->Traj_loc[547].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[547];
  localDW->Traj_loc[548].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[548];
  localDW->Traj_loc[549].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[549];
  localDW->Traj_loc[550].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[550];
  localDW->Traj_loc[551].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[551];
  localDW->Traj_loc[552].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[552];
  localDW->Traj_loc[553].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[553];
  localDW->Traj_loc[554].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[554];
  localDW->Traj_loc[555].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[555];
  localDW->Traj_loc[556].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[556];
  localDW->Traj_loc[557].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[557];
  localDW->Traj_loc[558].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[558];
  localDW->Traj_loc[559].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[559];
  localDW->Traj_loc[560].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[560];
  localDW->Traj_loc[561].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[561];
  localDW->Traj_loc[562].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[562];
  localDW->Traj_loc[563].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[563];
  localDW->Traj_loc[564].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[564];
  localDW->Traj_loc[565].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[565];
  localDW->Traj_loc[566].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[566];
  localDW->Traj_loc[567].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[567];
  localDW->Traj_loc[568].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[568];
  localDW->Traj_loc[569].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[569];
  localDW->Traj_loc[570].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[570];
  localDW->Traj_loc[571].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[571];
  localDW->Traj_loc[572].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[572];
  localDW->Traj_loc[573].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[573];
  localDW->Traj_loc[574].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[574];
  localDW->Traj_loc[575].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[575];
  localDW->Traj_loc[576].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[576];
  localDW->Traj_loc[577].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[577];
  localDW->Traj_loc[578].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[578];
  localDW->Traj_loc[579].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[579];
  localDW->Traj_loc[580].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[580];
  localDW->Traj_loc[581].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[581];
  localDW->Traj_loc[582].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[582];
  localDW->Traj_loc[583].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[583];
  localDW->Traj_loc[584].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[584];
  localDW->Traj_loc[585].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[585];
  localDW->Traj_loc[586].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[586];
  localDW->Traj_loc[587].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[587];
  localDW->Traj_loc[588].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[588];
  localDW->Traj_loc[589].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[589];
  localDW->Traj_loc[590].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[590];
  localDW->Traj_loc[591].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[591];
  localDW->Traj_loc[592].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[592];
  localDW->Traj_loc[593].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[593];
  localDW->Traj_loc[594].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[594];
  localDW->Traj_loc[595].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[595];
  localDW->Traj_loc[596].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[596];
  localDW->Traj_loc[597].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[597];
  localDW->Traj_loc[598].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[598];
  localDW->Traj_loc[599].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[599];
  localDW->Traj_loc[600].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[600];
  localDW->Traj_loc[601].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[601];
  localDW->Traj_loc[602].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[602];
  localDW->Traj_loc[603].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[603];
  localDW->Traj_loc[604].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[604];
  localDW->Traj_loc[605].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[605];
  localDW->Traj_loc[606].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[606];
  localDW->Traj_loc[607].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[607];
  localDW->Traj_loc[608].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[608];
  localDW->Traj_loc[609].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[609];
  localDW->Traj_loc[610].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[610];
  localDW->Traj_loc[611].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[611];
  localDW->Traj_loc[612].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[612];
  localDW->Traj_loc[613].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[613];
  localDW->Traj_loc[614].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[614];
  localDW->Traj_loc[615].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[615];
  localDW->Traj_loc[616].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[616];
  localDW->Traj_loc[617].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[617];
  localDW->Traj_loc[618].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[618];
  localDW->Traj_loc[619].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[619];
  localDW->Traj_loc[620].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[620];
  localDW->Traj_loc[621].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[621];
  localDW->Traj_loc[622].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[622];
  localDW->Traj_loc[623].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[623];
  localDW->Traj_loc[624].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[624];
  localDW->Traj_loc[625].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[625];
  localDW->Traj_loc[626].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[626];
  localDW->Traj_loc[627].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[627];
  localDW->Traj_loc[628].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[628];
  localDW->Traj_loc[629].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[629];
  localDW->Traj_loc[630].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[630];
  localDW->Traj_loc[631].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[631];
  localDW->Traj_loc[632].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[632];
  localDW->Traj_loc[633].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[633];
  localDW->Traj_loc[634].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[634];
  localDW->Traj_loc[635].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[635];
  localDW->Traj_loc[636].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[636];
  localDW->Traj_loc[637].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[637];
  localDW->Traj_loc[638].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[638];
  localDW->Traj_loc[639].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[639];
  localDW->Traj_loc[640].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[640];
  localDW->Traj_loc[641].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[641];
  localDW->Traj_loc[642].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[642];
  localDW->Traj_loc[643].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[643];
  localDW->Traj_loc[644].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[644];
  localDW->Traj_loc[645].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[645];
  localDW->Traj_loc[646].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[646];
  localDW->Traj_loc[647].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[647];
  localDW->Traj_loc[648].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[648];
  localDW->Traj_loc[649].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[649];
  localDW->Traj_loc[650].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[650];
  localDW->Traj_loc[651].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[651];
  localDW->Traj_loc[652].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[652];
  localDW->Traj_loc[653].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[653];
  localDW->Traj_loc[654].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[654];
  localDW->Traj_loc[655].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[655];
  localDW->Traj_loc[656].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[656];
  localDW->Traj_loc[657].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[657];
  localDW->Traj_loc[658].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[658];
  localDW->Traj_loc[659].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[659];
  localDW->Traj_loc[660].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[660];
  localDW->Traj_loc[661].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[661];
  localDW->Traj_loc[662].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[662];
  localDW->Traj_loc[663].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[663];
  localDW->Traj_loc[664].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[664];
  localDW->Traj_loc[665].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[665];
  localDW->Traj_loc[666].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[666];
  localDW->Traj_loc[667].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[667];
  localDW->Traj_loc[668].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[668];
  localDW->Traj_loc[669].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[669];
  localDW->Traj_loc[670].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[670];
  localDW->Traj_loc[671].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[671];
  localDW->Traj_loc[672].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[672];
  localDW->Traj_loc[673].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[673];
  localDW->Traj_loc[674].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[674];
  localDW->Traj_loc[675].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[675];
  localDW->Traj_loc[676].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[676];
  localDW->Traj_loc[677].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[677];
  localDW->Traj_loc[678].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[678];
  localDW->Traj_loc[679].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[679];
  localDW->Traj_loc[680].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[680];
  localDW->Traj_loc[681].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[681];
  localDW->Traj_loc[682].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[682];
  localDW->Traj_loc[683].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[683];
  localDW->Traj_loc[684].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[684];
  localDW->Traj_loc[685].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[685];
  localDW->Traj_loc[686].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[686];
  localDW->Traj_loc[687].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[687];
  localDW->Traj_loc[688].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[688];
  localDW->Traj_loc[689].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[689];
  localDW->Traj_loc[690].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[690];
  localDW->Traj_loc[691].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[691];
  localDW->Traj_loc[692].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[692];
  localDW->Traj_loc[693].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[693];
  localDW->Traj_loc[694].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[694];
  localDW->Traj_loc[695].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[695];
  localDW->Traj_loc[696].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[696];
  localDW->Traj_loc[697].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[697];
  localDW->Traj_loc[698].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[698];
  localDW->Traj_loc[699].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[699];
  localDW->Traj_loc[700].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[700];
  localDW->Traj_loc[701].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[701];
  localDW->Traj_loc[702].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[702];
  localDW->Traj_loc[703].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[703];
  localDW->Traj_loc[704].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[704];
  localDW->Traj_loc[705].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[705];
  localDW->Traj_loc[706].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[706];
  localDW->Traj_loc[707].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[707];
  localDW->Traj_loc[708].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[708];
  localDW->Traj_loc[709].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[709];
  localDW->Traj_loc[710].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[710];
  localDW->Traj_loc[711].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[711];
  localDW->Traj_loc[712].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[712];
  localDW->Traj_loc[713].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[713];
  localDW->Traj_loc[714].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[714];
  localDW->Traj_loc[715].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[715];
  localDW->Traj_loc[716].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[716];
  localDW->Traj_loc[717].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[717];
  localDW->Traj_loc[718].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[718];
  localDW->Traj_loc[719].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[719];
  localDW->Traj_loc[720].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[720];
  localDW->Traj_loc[721].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[721];
  localDW->Traj_loc[722].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[722];
  localDW->Traj_loc[723].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[723];
  localDW->Traj_loc[724].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[724];
  localDW->Traj_loc[725].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[725];
  localDW->Traj_loc[726].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[726];
  localDW->Traj_loc[727].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[727];
  localDW->Traj_loc[728].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[728];
  localDW->Traj_loc[729].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[729];
  localDW->Traj_loc[730].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[730];
  localDW->Traj_loc[731].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[731];
  localDW->Traj_loc[732].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[732];
  localDW->Traj_loc[733].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[733];
  localDW->Traj_loc[734].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[734];
  localDW->Traj_loc[735].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[735];
  localDW->Traj_loc[736].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[736];
  localDW->Traj_loc[737].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[737];
  localDW->Traj_loc[738].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[738];
  localDW->Traj_loc[739].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[739];
  localDW->Traj_loc[740].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[740];
  localDW->Traj_loc[741].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[741];
  localDW->Traj_loc[742].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[742];
  localDW->Traj_loc[743].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[743];
  localDW->Traj_loc[744].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[744];
  localDW->Traj_loc[745].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[745];
  localDW->Traj_loc[746].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[746];
  localDW->Traj_loc[747].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[747];
  localDW->Traj_loc[748].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[748];
  localDW->Traj_loc[749].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[749];
  localDW->Traj_loc[750].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[750];
  localDW->Traj_loc[751].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[751];
  localDW->Traj_loc[752].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[752];
  localDW->Traj_loc[753].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[753];
  localDW->Traj_loc[754].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[754];
  localDW->Traj_loc[755].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[755];
  localDW->Traj_loc[756].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[756];
  localDW->Traj_loc[757].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[757];
  localDW->Traj_loc[758].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[758];
  localDW->Traj_loc[759].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[759];
  localDW->Traj_loc[760].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[760];
  localDW->Traj_loc[761].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[761];
  localDW->Traj_loc[762].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[762];
  localDW->Traj_loc[763].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[763];
  localDW->Traj_loc[764].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[764];
  localDW->Traj_loc[765].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[765];
  localDW->Traj_loc[766].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[766];
  localDW->Traj_loc[767].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[767];
  localDW->Traj_loc[768].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[768];
  localDW->Traj_loc[769].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[769];
  localDW->Traj_loc[770].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[770];
  localDW->Traj_loc[771].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[771];
  localDW->Traj_loc[772].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[772];
  localDW->Traj_loc[773].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[773];
  localDW->Traj_loc[774].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[774];
  localDW->Traj_loc[775].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[775];
  localDW->Traj_loc[776].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[776];
  localDW->Traj_loc[777].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[777];
  localDW->Traj_loc[778].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[778];
  localDW->Traj_loc[779].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[779];
  localDW->Traj_loc[780].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[780];
  localDW->Traj_loc[781].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[781];
  localDW->Traj_loc[782].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[782];
  localDW->Traj_loc[783].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[783];
  localDW->Traj_loc[784].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[784];
  localDW->Traj_loc[785].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[785];
  localDW->Traj_loc[786].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[786];
  localDW->Traj_loc[787].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[787];
  localDW->Traj_loc[788].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[788];
  localDW->Traj_loc[789].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[789];
  localDW->Traj_loc[790].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[790];
  localDW->Traj_loc[791].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[791];
  localDW->Traj_loc[792].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[792];
  localDW->Traj_loc[793].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[793];
  localDW->Traj_loc[794].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[794];
  localDW->Traj_loc[795].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[795];
  localDW->Traj_loc[796].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[796];
  localDW->Traj_loc[797].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[797];
  localDW->Traj_loc[798].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[798];
  localDW->Traj_loc[799].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[799];
  localDW->Traj_loc[800].f1[0] = localDW->rtb_alpha_idx_3 * rtP.x[800];
  localDW->Traj_loc[0].f1[0] += localDW->Derivative2 * rtP.y[0];
  localDW->Traj_loc[1].f1[0] += localDW->Derivative2 * rtP.y[1];
  localDW->Traj_loc[2].f1[0] += localDW->Derivative2 * rtP.y[2];
  localDW->Traj_loc[3].f1[0] += localDW->Derivative2 * rtP.y[3];
  localDW->Traj_loc[4].f1[0] += localDW->Derivative2 * rtP.y[4];
  localDW->Traj_loc[5].f1[0] += localDW->Derivative2 * rtP.y[5];
  localDW->Traj_loc[6].f1[0] += localDW->Derivative2 * rtP.y[6];
  localDW->Traj_loc[7].f1[0] += localDW->Derivative2 * rtP.y[7];
  localDW->Traj_loc[8].f1[0] += localDW->Derivative2 * rtP.y[8];
  localDW->Traj_loc[9].f1[0] += localDW->Derivative2 * rtP.y[9];
  localDW->Traj_loc[10].f1[0] += localDW->Derivative2 * rtP.y[10];
  localDW->Traj_loc[11].f1[0] += localDW->Derivative2 * rtP.y[11];
  localDW->Traj_loc[12].f1[0] += localDW->Derivative2 * rtP.y[12];
  localDW->Traj_loc[13].f1[0] += localDW->Derivative2 * rtP.y[13];
  localDW->Traj_loc[14].f1[0] += localDW->Derivative2 * rtP.y[14];
  localDW->Traj_loc[15].f1[0] += localDW->Derivative2 * rtP.y[15];
  localDW->Traj_loc[16].f1[0] += localDW->Derivative2 * rtP.y[16];
  localDW->Traj_loc[17].f1[0] += localDW->Derivative2 * rtP.y[17];
  localDW->Traj_loc[18].f1[0] += localDW->Derivative2 * rtP.y[18];
  localDW->Traj_loc[19].f1[0] += localDW->Derivative2 * rtP.y[19];
  localDW->Traj_loc[20].f1[0] += localDW->Derivative2 * rtP.y[20];
  localDW->Traj_loc[21].f1[0] += localDW->Derivative2 * rtP.y[21];
  localDW->Traj_loc[22].f1[0] += localDW->Derivative2 * rtP.y[22];
  localDW->Traj_loc[23].f1[0] += localDW->Derivative2 * rtP.y[23];
  localDW->Traj_loc[24].f1[0] += localDW->Derivative2 * rtP.y[24];
  localDW->Traj_loc[25].f1[0] += localDW->Derivative2 * rtP.y[25];
  localDW->Traj_loc[26].f1[0] += localDW->Derivative2 * rtP.y[26];
  localDW->Traj_loc[27].f1[0] += localDW->Derivative2 * rtP.y[27];
  localDW->Traj_loc[28].f1[0] += localDW->Derivative2 * rtP.y[28];
  localDW->Traj_loc[29].f1[0] += localDW->Derivative2 * rtP.y[29];
  localDW->Traj_loc[30].f1[0] += localDW->Derivative2 * rtP.y[30];
  localDW->Traj_loc[31].f1[0] += localDW->Derivative2 * rtP.y[31];
  localDW->Traj_loc[32].f1[0] += localDW->Derivative2 * rtP.y[32];
  localDW->Traj_loc[33].f1[0] += localDW->Derivative2 * rtP.y[33];
  localDW->Traj_loc[34].f1[0] += localDW->Derivative2 * rtP.y[34];
  localDW->Traj_loc[35].f1[0] += localDW->Derivative2 * rtP.y[35];
  localDW->Traj_loc[36].f1[0] += localDW->Derivative2 * rtP.y[36];
  localDW->Traj_loc[37].f1[0] += localDW->Derivative2 * rtP.y[37];
  localDW->Traj_loc[38].f1[0] += localDW->Derivative2 * rtP.y[38];
  localDW->Traj_loc[39].f1[0] += localDW->Derivative2 * rtP.y[39];
  localDW->Traj_loc[40].f1[0] += localDW->Derivative2 * rtP.y[40];
  localDW->Traj_loc[41].f1[0] += localDW->Derivative2 * rtP.y[41];
  localDW->Traj_loc[42].f1[0] += localDW->Derivative2 * rtP.y[42];
  localDW->Traj_loc[43].f1[0] += localDW->Derivative2 * rtP.y[43];
  localDW->Traj_loc[44].f1[0] += localDW->Derivative2 * rtP.y[44];
  localDW->Traj_loc[45].f1[0] += localDW->Derivative2 * rtP.y[45];
  localDW->Traj_loc[46].f1[0] += localDW->Derivative2 * rtP.y[46];
  localDW->Traj_loc[47].f1[0] += localDW->Derivative2 * rtP.y[47];
  localDW->Traj_loc[48].f1[0] += localDW->Derivative2 * rtP.y[48];
  localDW->Traj_loc[49].f1[0] += localDW->Derivative2 * rtP.y[49];
  localDW->Traj_loc[50].f1[0] += localDW->Derivative2 * rtP.y[50];
  localDW->Traj_loc[51].f1[0] += localDW->Derivative2 * rtP.y[51];
  localDW->Traj_loc[52].f1[0] += localDW->Derivative2 * rtP.y[52];
  localDW->Traj_loc[53].f1[0] += localDW->Derivative2 * rtP.y[53];
  localDW->Traj_loc[54].f1[0] += localDW->Derivative2 * rtP.y[54];
  localDW->Traj_loc[55].f1[0] += localDW->Derivative2 * rtP.y[55];
  localDW->Traj_loc[56].f1[0] += localDW->Derivative2 * rtP.y[56];
  localDW->Traj_loc[57].f1[0] += localDW->Derivative2 * rtP.y[57];
  localDW->Traj_loc[58].f1[0] += localDW->Derivative2 * rtP.y[58];
  localDW->Traj_loc[59].f1[0] += localDW->Derivative2 * rtP.y[59];
  localDW->Traj_loc[60].f1[0] += localDW->Derivative2 * rtP.y[60];
  localDW->Traj_loc[61].f1[0] += localDW->Derivative2 * rtP.y[61];
  localDW->Traj_loc[62].f1[0] += localDW->Derivative2 * rtP.y[62];
  localDW->Traj_loc[63].f1[0] += localDW->Derivative2 * rtP.y[63];
  localDW->Traj_loc[64].f1[0] += localDW->Derivative2 * rtP.y[64];
  localDW->Traj_loc[65].f1[0] += localDW->Derivative2 * rtP.y[65];
  localDW->Traj_loc[66].f1[0] += localDW->Derivative2 * rtP.y[66];
  localDW->Traj_loc[67].f1[0] += localDW->Derivative2 * rtP.y[67];
  localDW->Traj_loc[68].f1[0] += localDW->Derivative2 * rtP.y[68];
  localDW->Traj_loc[69].f1[0] += localDW->Derivative2 * rtP.y[69];
  localDW->Traj_loc[70].f1[0] += localDW->Derivative2 * rtP.y[70];
  localDW->Traj_loc[71].f1[0] += localDW->Derivative2 * rtP.y[71];
  localDW->Traj_loc[72].f1[0] += localDW->Derivative2 * rtP.y[72];
  localDW->Traj_loc[73].f1[0] += localDW->Derivative2 * rtP.y[73];
  localDW->Traj_loc[74].f1[0] += localDW->Derivative2 * rtP.y[74];
  localDW->Traj_loc[75].f1[0] += localDW->Derivative2 * rtP.y[75];
  localDW->Traj_loc[76].f1[0] += localDW->Derivative2 * rtP.y[76];
  localDW->Traj_loc[77].f1[0] += localDW->Derivative2 * rtP.y[77];
  localDW->Traj_loc[78].f1[0] += localDW->Derivative2 * rtP.y[78];
  localDW->Traj_loc[79].f1[0] += localDW->Derivative2 * rtP.y[79];
  localDW->Traj_loc[80].f1[0] += localDW->Derivative2 * rtP.y[80];
  localDW->Traj_loc[81].f1[0] += localDW->Derivative2 * rtP.y[81];
  localDW->Traj_loc[82].f1[0] += localDW->Derivative2 * rtP.y[82];
  localDW->Traj_loc[83].f1[0] += localDW->Derivative2 * rtP.y[83];
  localDW->Traj_loc[84].f1[0] += localDW->Derivative2 * rtP.y[84];
  localDW->Traj_loc[85].f1[0] += localDW->Derivative2 * rtP.y[85];
  localDW->Traj_loc[86].f1[0] += localDW->Derivative2 * rtP.y[86];
  localDW->Traj_loc[87].f1[0] += localDW->Derivative2 * rtP.y[87];
  localDW->Traj_loc[88].f1[0] += localDW->Derivative2 * rtP.y[88];
  localDW->Traj_loc[89].f1[0] += localDW->Derivative2 * rtP.y[89];
  localDW->Traj_loc[90].f1[0] += localDW->Derivative2 * rtP.y[90];
  localDW->Traj_loc[91].f1[0] += localDW->Derivative2 * rtP.y[91];
  localDW->Traj_loc[92].f1[0] += localDW->Derivative2 * rtP.y[92];
  localDW->Traj_loc[93].f1[0] += localDW->Derivative2 * rtP.y[93];
  localDW->Traj_loc[94].f1[0] += localDW->Derivative2 * rtP.y[94];
  localDW->Traj_loc[95].f1[0] += localDW->Derivative2 * rtP.y[95];
  localDW->Traj_loc[96].f1[0] += localDW->Derivative2 * rtP.y[96];
  localDW->Traj_loc[97].f1[0] += localDW->Derivative2 * rtP.y[97];
  localDW->Traj_loc[98].f1[0] += localDW->Derivative2 * rtP.y[98];
  localDW->Traj_loc[99].f1[0] += localDW->Derivative2 * rtP.y[99];
  localDW->Traj_loc[100].f1[0] += localDW->Derivative2 * rtP.y[100];
  localDW->Traj_loc[101].f1[0] += localDW->Derivative2 * rtP.y[101];
  localDW->Traj_loc[102].f1[0] += localDW->Derivative2 * rtP.y[102];
  localDW->Traj_loc[103].f1[0] += localDW->Derivative2 * rtP.y[103];
  localDW->Traj_loc[104].f1[0] += localDW->Derivative2 * rtP.y[104];
  localDW->Traj_loc[105].f1[0] += localDW->Derivative2 * rtP.y[105];
  localDW->Traj_loc[106].f1[0] += localDW->Derivative2 * rtP.y[106];
  localDW->Traj_loc[107].f1[0] += localDW->Derivative2 * rtP.y[107];
  localDW->Traj_loc[108].f1[0] += localDW->Derivative2 * rtP.y[108];
  localDW->Traj_loc[109].f1[0] += localDW->Derivative2 * rtP.y[109];
  localDW->Traj_loc[110].f1[0] += localDW->Derivative2 * rtP.y[110];
  localDW->Traj_loc[111].f1[0] += localDW->Derivative2 * rtP.y[111];
  localDW->Traj_loc[112].f1[0] += localDW->Derivative2 * rtP.y[112];
  localDW->Traj_loc[113].f1[0] += localDW->Derivative2 * rtP.y[113];
  localDW->Traj_loc[114].f1[0] += localDW->Derivative2 * rtP.y[114];
  localDW->Traj_loc[115].f1[0] += localDW->Derivative2 * rtP.y[115];
  localDW->Traj_loc[116].f1[0] += localDW->Derivative2 * rtP.y[116];
  localDW->Traj_loc[117].f1[0] += localDW->Derivative2 * rtP.y[117];
  localDW->Traj_loc[118].f1[0] += localDW->Derivative2 * rtP.y[118];
  localDW->Traj_loc[119].f1[0] += localDW->Derivative2 * rtP.y[119];
  localDW->Traj_loc[120].f1[0] += localDW->Derivative2 * rtP.y[120];
  localDW->Traj_loc[121].f1[0] += localDW->Derivative2 * rtP.y[121];
  localDW->Traj_loc[122].f1[0] += localDW->Derivative2 * rtP.y[122];
  localDW->Traj_loc[123].f1[0] += localDW->Derivative2 * rtP.y[123];
  localDW->Traj_loc[124].f1[0] += localDW->Derivative2 * rtP.y[124];
  localDW->Traj_loc[125].f1[0] += localDW->Derivative2 * rtP.y[125];
  localDW->Traj_loc[126].f1[0] += localDW->Derivative2 * rtP.y[126];
  localDW->Traj_loc[127].f1[0] += localDW->Derivative2 * rtP.y[127];
  localDW->Traj_loc[128].f1[0] += localDW->Derivative2 * rtP.y[128];
  localDW->Traj_loc[129].f1[0] += localDW->Derivative2 * rtP.y[129];
  localDW->Traj_loc[130].f1[0] += localDW->Derivative2 * rtP.y[130];
  localDW->Traj_loc[131].f1[0] += localDW->Derivative2 * rtP.y[131];
  localDW->Traj_loc[132].f1[0] += localDW->Derivative2 * rtP.y[132];
  localDW->Traj_loc[133].f1[0] += localDW->Derivative2 * rtP.y[133];
  localDW->Traj_loc[134].f1[0] += localDW->Derivative2 * rtP.y[134];
  localDW->Traj_loc[135].f1[0] += localDW->Derivative2 * rtP.y[135];
  localDW->Traj_loc[136].f1[0] += localDW->Derivative2 * rtP.y[136];
  localDW->Traj_loc[137].f1[0] += localDW->Derivative2 * rtP.y[137];
  localDW->Traj_loc[138].f1[0] += localDW->Derivative2 * rtP.y[138];
  localDW->Traj_loc[139].f1[0] += localDW->Derivative2 * rtP.y[139];
  localDW->Traj_loc[140].f1[0] += localDW->Derivative2 * rtP.y[140];
  localDW->Traj_loc[141].f1[0] += localDW->Derivative2 * rtP.y[141];
  localDW->Traj_loc[142].f1[0] += localDW->Derivative2 * rtP.y[142];
  localDW->Traj_loc[143].f1[0] += localDW->Derivative2 * rtP.y[143];
  localDW->Traj_loc[144].f1[0] += localDW->Derivative2 * rtP.y[144];
  localDW->Traj_loc[145].f1[0] += localDW->Derivative2 * rtP.y[145];
  localDW->Traj_loc[146].f1[0] += localDW->Derivative2 * rtP.y[146];
  localDW->Traj_loc[147].f1[0] += localDW->Derivative2 * rtP.y[147];
  localDW->Traj_loc[148].f1[0] += localDW->Derivative2 * rtP.y[148];
  localDW->Traj_loc[149].f1[0] += localDW->Derivative2 * rtP.y[149];
  localDW->Traj_loc[150].f1[0] += localDW->Derivative2 * rtP.y[150];
  localDW->Traj_loc[151].f1[0] += localDW->Derivative2 * rtP.y[151];
  localDW->Traj_loc[152].f1[0] += localDW->Derivative2 * rtP.y[152];
  localDW->Traj_loc[153].f1[0] += localDW->Derivative2 * rtP.y[153];
  localDW->Traj_loc[154].f1[0] += localDW->Derivative2 * rtP.y[154];
  localDW->Traj_loc[155].f1[0] += localDW->Derivative2 * rtP.y[155];
  localDW->Traj_loc[156].f1[0] += localDW->Derivative2 * rtP.y[156];
  localDW->Traj_loc[157].f1[0] += localDW->Derivative2 * rtP.y[157];
  localDW->Traj_loc[158].f1[0] += localDW->Derivative2 * rtP.y[158];
  localDW->Traj_loc[159].f1[0] += localDW->Derivative2 * rtP.y[159];
  localDW->Traj_loc[160].f1[0] += localDW->Derivative2 * rtP.y[160];
  localDW->Traj_loc[161].f1[0] += localDW->Derivative2 * rtP.y[161];
  localDW->Traj_loc[162].f1[0] += localDW->Derivative2 * rtP.y[162];
  localDW->Traj_loc[163].f1[0] += localDW->Derivative2 * rtP.y[163];
  localDW->Traj_loc[164].f1[0] += localDW->Derivative2 * rtP.y[164];
  localDW->Traj_loc[165].f1[0] += localDW->Derivative2 * rtP.y[165];
  localDW->Traj_loc[166].f1[0] += localDW->Derivative2 * rtP.y[166];
  localDW->Traj_loc[167].f1[0] += localDW->Derivative2 * rtP.y[167];
  localDW->Traj_loc[168].f1[0] += localDW->Derivative2 * rtP.y[168];
  localDW->Traj_loc[169].f1[0] += localDW->Derivative2 * rtP.y[169];
  localDW->Traj_loc[170].f1[0] += localDW->Derivative2 * rtP.y[170];
  localDW->Traj_loc[171].f1[0] += localDW->Derivative2 * rtP.y[171];
  localDW->Traj_loc[172].f1[0] += localDW->Derivative2 * rtP.y[172];
  localDW->Traj_loc[173].f1[0] += localDW->Derivative2 * rtP.y[173];
  localDW->Traj_loc[174].f1[0] += localDW->Derivative2 * rtP.y[174];
  localDW->Traj_loc[175].f1[0] += localDW->Derivative2 * rtP.y[175];
  localDW->Traj_loc[176].f1[0] += localDW->Derivative2 * rtP.y[176];
  localDW->Traj_loc[177].f1[0] += localDW->Derivative2 * rtP.y[177];
  localDW->Traj_loc[178].f1[0] += localDW->Derivative2 * rtP.y[178];
  localDW->Traj_loc[179].f1[0] += localDW->Derivative2 * rtP.y[179];
  localDW->Traj_loc[180].f1[0] += localDW->Derivative2 * rtP.y[180];
  localDW->Traj_loc[181].f1[0] += localDW->Derivative2 * rtP.y[181];
  localDW->Traj_loc[182].f1[0] += localDW->Derivative2 * rtP.y[182];
  localDW->Traj_loc[183].f1[0] += localDW->Derivative2 * rtP.y[183];
  localDW->Traj_loc[184].f1[0] += localDW->Derivative2 * rtP.y[184];
  localDW->Traj_loc[185].f1[0] += localDW->Derivative2 * rtP.y[185];
  localDW->Traj_loc[186].f1[0] += localDW->Derivative2 * rtP.y[186];
  localDW->Traj_loc[187].f1[0] += localDW->Derivative2 * rtP.y[187];
  localDW->Traj_loc[188].f1[0] += localDW->Derivative2 * rtP.y[188];
  localDW->Traj_loc[189].f1[0] += localDW->Derivative2 * rtP.y[189];
  localDW->Traj_loc[190].f1[0] += localDW->Derivative2 * rtP.y[190];
  localDW->Traj_loc[191].f1[0] += localDW->Derivative2 * rtP.y[191];
  localDW->Traj_loc[192].f1[0] += localDW->Derivative2 * rtP.y[192];
  localDW->Traj_loc[193].f1[0] += localDW->Derivative2 * rtP.y[193];
  localDW->Traj_loc[194].f1[0] += localDW->Derivative2 * rtP.y[194];
  localDW->Traj_loc[195].f1[0] += localDW->Derivative2 * rtP.y[195];
  localDW->Traj_loc[196].f1[0] += localDW->Derivative2 * rtP.y[196];
  localDW->Traj_loc[197].f1[0] += localDW->Derivative2 * rtP.y[197];
  localDW->Traj_loc[198].f1[0] += localDW->Derivative2 * rtP.y[198];
  localDW->Traj_loc[199].f1[0] += localDW->Derivative2 * rtP.y[199];
  localDW->Traj_loc[200].f1[0] += localDW->Derivative2 * rtP.y[200];
  localDW->Traj_loc[201].f1[0] += localDW->Derivative2 * rtP.y[201];
  localDW->Traj_loc[202].f1[0] += localDW->Derivative2 * rtP.y[202];
  localDW->Traj_loc[203].f1[0] += localDW->Derivative2 * rtP.y[203];
  localDW->Traj_loc[204].f1[0] += localDW->Derivative2 * rtP.y[204];
  localDW->Traj_loc[205].f1[0] += localDW->Derivative2 * rtP.y[205];
  localDW->Traj_loc[206].f1[0] += localDW->Derivative2 * rtP.y[206];
  localDW->Traj_loc[207].f1[0] += localDW->Derivative2 * rtP.y[207];
  localDW->Traj_loc[208].f1[0] += localDW->Derivative2 * rtP.y[208];
  localDW->Traj_loc[209].f1[0] += localDW->Derivative2 * rtP.y[209];
  localDW->Traj_loc[210].f1[0] += localDW->Derivative2 * rtP.y[210];
  localDW->Traj_loc[211].f1[0] += localDW->Derivative2 * rtP.y[211];
  localDW->Traj_loc[212].f1[0] += localDW->Derivative2 * rtP.y[212];
  localDW->Traj_loc[213].f1[0] += localDW->Derivative2 * rtP.y[213];
  localDW->Traj_loc[214].f1[0] += localDW->Derivative2 * rtP.y[214];
  localDW->Traj_loc[215].f1[0] += localDW->Derivative2 * rtP.y[215];
  localDW->Traj_loc[216].f1[0] += localDW->Derivative2 * rtP.y[216];
  localDW->Traj_loc[217].f1[0] += localDW->Derivative2 * rtP.y[217];
  localDW->Traj_loc[218].f1[0] += localDW->Derivative2 * rtP.y[218];
  localDW->Traj_loc[219].f1[0] += localDW->Derivative2 * rtP.y[219];
  localDW->Traj_loc[220].f1[0] += localDW->Derivative2 * rtP.y[220];
  localDW->Traj_loc[221].f1[0] += localDW->Derivative2 * rtP.y[221];
  localDW->Traj_loc[222].f1[0] += localDW->Derivative2 * rtP.y[222];
  localDW->Traj_loc[223].f1[0] += localDW->Derivative2 * rtP.y[223];
  localDW->Traj_loc[224].f1[0] += localDW->Derivative2 * rtP.y[224];
  localDW->Traj_loc[225].f1[0] += localDW->Derivative2 * rtP.y[225];
  localDW->Traj_loc[226].f1[0] += localDW->Derivative2 * rtP.y[226];
  localDW->Traj_loc[227].f1[0] += localDW->Derivative2 * rtP.y[227];
  localDW->Traj_loc[228].f1[0] += localDW->Derivative2 * rtP.y[228];
  localDW->Traj_loc[229].f1[0] += localDW->Derivative2 * rtP.y[229];
  localDW->Traj_loc[230].f1[0] += localDW->Derivative2 * rtP.y[230];
  localDW->Traj_loc[231].f1[0] += localDW->Derivative2 * rtP.y[231];
  localDW->Traj_loc[232].f1[0] += localDW->Derivative2 * rtP.y[232];
  localDW->Traj_loc[233].f1[0] += localDW->Derivative2 * rtP.y[233];
  localDW->Traj_loc[234].f1[0] += localDW->Derivative2 * rtP.y[234];
  localDW->Traj_loc[235].f1[0] += localDW->Derivative2 * rtP.y[235];
  localDW->Traj_loc[236].f1[0] += localDW->Derivative2 * rtP.y[236];
  localDW->Traj_loc[237].f1[0] += localDW->Derivative2 * rtP.y[237];
  localDW->Traj_loc[238].f1[0] += localDW->Derivative2 * rtP.y[238];
  localDW->Traj_loc[239].f1[0] += localDW->Derivative2 * rtP.y[239];
  localDW->Traj_loc[240].f1[0] += localDW->Derivative2 * rtP.y[240];
  localDW->Traj_loc[241].f1[0] += localDW->Derivative2 * rtP.y[241];
  localDW->Traj_loc[242].f1[0] += localDW->Derivative2 * rtP.y[242];
  localDW->Traj_loc[243].f1[0] += localDW->Derivative2 * rtP.y[243];
  localDW->Traj_loc[244].f1[0] += localDW->Derivative2 * rtP.y[244];
  localDW->Traj_loc[245].f1[0] += localDW->Derivative2 * rtP.y[245];
  localDW->Traj_loc[246].f1[0] += localDW->Derivative2 * rtP.y[246];
  localDW->Traj_loc[247].f1[0] += localDW->Derivative2 * rtP.y[247];
  localDW->Traj_loc[248].f1[0] += localDW->Derivative2 * rtP.y[248];
  localDW->Traj_loc[249].f1[0] += localDW->Derivative2 * rtP.y[249];
  localDW->Traj_loc[250].f1[0] += localDW->Derivative2 * rtP.y[250];
  localDW->Traj_loc[251].f1[0] += localDW->Derivative2 * rtP.y[251];
  localDW->Traj_loc[252].f1[0] += localDW->Derivative2 * rtP.y[252];
  localDW->Traj_loc[253].f1[0] += localDW->Derivative2 * rtP.y[253];
  localDW->Traj_loc[254].f1[0] += localDW->Derivative2 * rtP.y[254];
  localDW->Traj_loc[255].f1[0] += localDW->Derivative2 * rtP.y[255];
  localDW->Traj_loc[256].f1[0] += localDW->Derivative2 * rtP.y[256];
  localDW->Traj_loc[257].f1[0] += localDW->Derivative2 * rtP.y[257];
  localDW->Traj_loc[258].f1[0] += localDW->Derivative2 * rtP.y[258];
  localDW->Traj_loc[259].f1[0] += localDW->Derivative2 * rtP.y[259];
  localDW->Traj_loc[260].f1[0] += localDW->Derivative2 * rtP.y[260];
  localDW->Traj_loc[261].f1[0] += localDW->Derivative2 * rtP.y[261];
  localDW->Traj_loc[262].f1[0] += localDW->Derivative2 * rtP.y[262];
  localDW->Traj_loc[263].f1[0] += localDW->Derivative2 * rtP.y[263];
  localDW->Traj_loc[264].f1[0] += localDW->Derivative2 * rtP.y[264];
  localDW->Traj_loc[265].f1[0] += localDW->Derivative2 * rtP.y[265];
  localDW->Traj_loc[266].f1[0] += localDW->Derivative2 * rtP.y[266];
  localDW->Traj_loc[267].f1[0] += localDW->Derivative2 * rtP.y[267];
  localDW->Traj_loc[268].f1[0] += localDW->Derivative2 * rtP.y[268];
  localDW->Traj_loc[269].f1[0] += localDW->Derivative2 * rtP.y[269];
  localDW->Traj_loc[270].f1[0] += localDW->Derivative2 * rtP.y[270];
  localDW->Traj_loc[271].f1[0] += localDW->Derivative2 * rtP.y[271];
  localDW->Traj_loc[272].f1[0] += localDW->Derivative2 * rtP.y[272];
  localDW->Traj_loc[273].f1[0] += localDW->Derivative2 * rtP.y[273];
  localDW->Traj_loc[274].f1[0] += localDW->Derivative2 * rtP.y[274];
  localDW->Traj_loc[275].f1[0] += localDW->Derivative2 * rtP.y[275];
  localDW->Traj_loc[276].f1[0] += localDW->Derivative2 * rtP.y[276];
  localDW->Traj_loc[277].f1[0] += localDW->Derivative2 * rtP.y[277];
  localDW->Traj_loc[278].f1[0] += localDW->Derivative2 * rtP.y[278];
  localDW->Traj_loc[279].f1[0] += localDW->Derivative2 * rtP.y[279];
  localDW->Traj_loc[280].f1[0] += localDW->Derivative2 * rtP.y[280];
  localDW->Traj_loc[281].f1[0] += localDW->Derivative2 * rtP.y[281];
  localDW->Traj_loc[282].f1[0] += localDW->Derivative2 * rtP.y[282];
  localDW->Traj_loc[283].f1[0] += localDW->Derivative2 * rtP.y[283];
  localDW->Traj_loc[284].f1[0] += localDW->Derivative2 * rtP.y[284];
  localDW->Traj_loc[285].f1[0] += localDW->Derivative2 * rtP.y[285];
  localDW->Traj_loc[286].f1[0] += localDW->Derivative2 * rtP.y[286];
  localDW->Traj_loc[287].f1[0] += localDW->Derivative2 * rtP.y[287];
  localDW->Traj_loc[288].f1[0] += localDW->Derivative2 * rtP.y[288];
  localDW->Traj_loc[289].f1[0] += localDW->Derivative2 * rtP.y[289];
  localDW->Traj_loc[290].f1[0] += localDW->Derivative2 * rtP.y[290];
  localDW->Traj_loc[291].f1[0] += localDW->Derivative2 * rtP.y[291];
  localDW->Traj_loc[292].f1[0] += localDW->Derivative2 * rtP.y[292];
  localDW->Traj_loc[293].f1[0] += localDW->Derivative2 * rtP.y[293];
  localDW->Traj_loc[294].f1[0] += localDW->Derivative2 * rtP.y[294];
  localDW->Traj_loc[295].f1[0] += localDW->Derivative2 * rtP.y[295];
  localDW->Traj_loc[296].f1[0] += localDW->Derivative2 * rtP.y[296];
  localDW->Traj_loc[297].f1[0] += localDW->Derivative2 * rtP.y[297];
  localDW->Traj_loc[298].f1[0] += localDW->Derivative2 * rtP.y[298];
  localDW->Traj_loc[299].f1[0] += localDW->Derivative2 * rtP.y[299];
  localDW->Traj_loc[300].f1[0] += localDW->Derivative2 * rtP.y[300];
  localDW->Traj_loc[301].f1[0] += localDW->Derivative2 * rtP.y[301];
  localDW->Traj_loc[302].f1[0] += localDW->Derivative2 * rtP.y[302];
  localDW->Traj_loc[303].f1[0] += localDW->Derivative2 * rtP.y[303];
  localDW->Traj_loc[304].f1[0] += localDW->Derivative2 * rtP.y[304];
  localDW->Traj_loc[305].f1[0] += localDW->Derivative2 * rtP.y[305];
  localDW->Traj_loc[306].f1[0] += localDW->Derivative2 * rtP.y[306];
  localDW->Traj_loc[307].f1[0] += localDW->Derivative2 * rtP.y[307];
  localDW->Traj_loc[308].f1[0] += localDW->Derivative2 * rtP.y[308];
  localDW->Traj_loc[309].f1[0] += localDW->Derivative2 * rtP.y[309];
  localDW->Traj_loc[310].f1[0] += localDW->Derivative2 * rtP.y[310];
  localDW->Traj_loc[311].f1[0] += localDW->Derivative2 * rtP.y[311];
  localDW->Traj_loc[312].f1[0] += localDW->Derivative2 * rtP.y[312];
  localDW->Traj_loc[313].f1[0] += localDW->Derivative2 * rtP.y[313];
  localDW->Traj_loc[314].f1[0] += localDW->Derivative2 * rtP.y[314];
  localDW->Traj_loc[315].f1[0] += localDW->Derivative2 * rtP.y[315];
  localDW->Traj_loc[316].f1[0] += localDW->Derivative2 * rtP.y[316];
  localDW->Traj_loc[317].f1[0] += localDW->Derivative2 * rtP.y[317];
  localDW->Traj_loc[318].f1[0] += localDW->Derivative2 * rtP.y[318];
  localDW->Traj_loc[319].f1[0] += localDW->Derivative2 * rtP.y[319];
  localDW->Traj_loc[320].f1[0] += localDW->Derivative2 * rtP.y[320];
  localDW->Traj_loc[321].f1[0] += localDW->Derivative2 * rtP.y[321];
  localDW->Traj_loc[322].f1[0] += localDW->Derivative2 * rtP.y[322];
  localDW->Traj_loc[323].f1[0] += localDW->Derivative2 * rtP.y[323];
  localDW->Traj_loc[324].f1[0] += localDW->Derivative2 * rtP.y[324];
  localDW->Traj_loc[325].f1[0] += localDW->Derivative2 * rtP.y[325];
  localDW->Traj_loc[326].f1[0] += localDW->Derivative2 * rtP.y[326];
  localDW->Traj_loc[327].f1[0] += localDW->Derivative2 * rtP.y[327];
  localDW->Traj_loc[328].f1[0] += localDW->Derivative2 * rtP.y[328];
  localDW->Traj_loc[329].f1[0] += localDW->Derivative2 * rtP.y[329];
  localDW->Traj_loc[330].f1[0] += localDW->Derivative2 * rtP.y[330];
  localDW->Traj_loc[331].f1[0] += localDW->Derivative2 * rtP.y[331];
  localDW->Traj_loc[332].f1[0] += localDW->Derivative2 * rtP.y[332];
  localDW->Traj_loc[333].f1[0] += localDW->Derivative2 * rtP.y[333];
  localDW->Traj_loc[334].f1[0] += localDW->Derivative2 * rtP.y[334];
  localDW->Traj_loc[335].f1[0] += localDW->Derivative2 * rtP.y[335];
  localDW->Traj_loc[336].f1[0] += localDW->Derivative2 * rtP.y[336];
  localDW->Traj_loc[337].f1[0] += localDW->Derivative2 * rtP.y[337];
  localDW->Traj_loc[338].f1[0] += localDW->Derivative2 * rtP.y[338];
  localDW->Traj_loc[339].f1[0] += localDW->Derivative2 * rtP.y[339];
  localDW->Traj_loc[340].f1[0] += localDW->Derivative2 * rtP.y[340];
  localDW->Traj_loc[341].f1[0] += localDW->Derivative2 * rtP.y[341];
  localDW->Traj_loc[342].f1[0] += localDW->Derivative2 * rtP.y[342];
  localDW->Traj_loc[343].f1[0] += localDW->Derivative2 * rtP.y[343];
  localDW->Traj_loc[344].f1[0] += localDW->Derivative2 * rtP.y[344];
  localDW->Traj_loc[345].f1[0] += localDW->Derivative2 * rtP.y[345];
  localDW->Traj_loc[346].f1[0] += localDW->Derivative2 * rtP.y[346];
  localDW->Traj_loc[347].f1[0] += localDW->Derivative2 * rtP.y[347];
  localDW->Traj_loc[348].f1[0] += localDW->Derivative2 * rtP.y[348];
  localDW->Traj_loc[349].f1[0] += localDW->Derivative2 * rtP.y[349];
  localDW->Traj_loc[350].f1[0] += localDW->Derivative2 * rtP.y[350];
  localDW->Traj_loc[351].f1[0] += localDW->Derivative2 * rtP.y[351];
  localDW->Traj_loc[352].f1[0] += localDW->Derivative2 * rtP.y[352];
  localDW->Traj_loc[353].f1[0] += localDW->Derivative2 * rtP.y[353];
  localDW->Traj_loc[354].f1[0] += localDW->Derivative2 * rtP.y[354];
  localDW->Traj_loc[355].f1[0] += localDW->Derivative2 * rtP.y[355];
  localDW->Traj_loc[356].f1[0] += localDW->Derivative2 * rtP.y[356];
  localDW->Traj_loc[357].f1[0] += localDW->Derivative2 * rtP.y[357];
  localDW->Traj_loc[358].f1[0] += localDW->Derivative2 * rtP.y[358];
  localDW->Traj_loc[359].f1[0] += localDW->Derivative2 * rtP.y[359];
  localDW->Traj_loc[360].f1[0] += localDW->Derivative2 * rtP.y[360];
  localDW->Traj_loc[361].f1[0] += localDW->Derivative2 * rtP.y[361];
  localDW->Traj_loc[362].f1[0] += localDW->Derivative2 * rtP.y[362];
  localDW->Traj_loc[363].f1[0] += localDW->Derivative2 * rtP.y[363];
  localDW->Traj_loc[364].f1[0] += localDW->Derivative2 * rtP.y[364];
  localDW->Traj_loc[365].f1[0] += localDW->Derivative2 * rtP.y[365];
  localDW->Traj_loc[366].f1[0] += localDW->Derivative2 * rtP.y[366];
  localDW->Traj_loc[367].f1[0] += localDW->Derivative2 * rtP.y[367];
  localDW->Traj_loc[368].f1[0] += localDW->Derivative2 * rtP.y[368];
  localDW->Traj_loc[369].f1[0] += localDW->Derivative2 * rtP.y[369];
  localDW->Traj_loc[370].f1[0] += localDW->Derivative2 * rtP.y[370];
  localDW->Traj_loc[371].f1[0] += localDW->Derivative2 * rtP.y[371];
  localDW->Traj_loc[372].f1[0] += localDW->Derivative2 * rtP.y[372];
  localDW->Traj_loc[373].f1[0] += localDW->Derivative2 * rtP.y[373];
  localDW->Traj_loc[374].f1[0] += localDW->Derivative2 * rtP.y[374];
  localDW->Traj_loc[375].f1[0] += localDW->Derivative2 * rtP.y[375];
  localDW->Traj_loc[376].f1[0] += localDW->Derivative2 * rtP.y[376];
  localDW->Traj_loc[377].f1[0] += localDW->Derivative2 * rtP.y[377];
  localDW->Traj_loc[378].f1[0] += localDW->Derivative2 * rtP.y[378];
  localDW->Traj_loc[379].f1[0] += localDW->Derivative2 * rtP.y[379];
  localDW->Traj_loc[380].f1[0] += localDW->Derivative2 * rtP.y[380];
  localDW->Traj_loc[381].f1[0] += localDW->Derivative2 * rtP.y[381];
  localDW->Traj_loc[382].f1[0] += localDW->Derivative2 * rtP.y[382];
  localDW->Traj_loc[383].f1[0] += localDW->Derivative2 * rtP.y[383];
  localDW->Traj_loc[384].f1[0] += localDW->Derivative2 * rtP.y[384];
  localDW->Traj_loc[385].f1[0] += localDW->Derivative2 * rtP.y[385];
  localDW->Traj_loc[386].f1[0] += localDW->Derivative2 * rtP.y[386];
  localDW->Traj_loc[387].f1[0] += localDW->Derivative2 * rtP.y[387];
  localDW->Traj_loc[388].f1[0] += localDW->Derivative2 * rtP.y[388];
  localDW->Traj_loc[389].f1[0] += localDW->Derivative2 * rtP.y[389];
  localDW->Traj_loc[390].f1[0] += localDW->Derivative2 * rtP.y[390];
  localDW->Traj_loc[391].f1[0] += localDW->Derivative2 * rtP.y[391];
  localDW->Traj_loc[392].f1[0] += localDW->Derivative2 * rtP.y[392];
  localDW->Traj_loc[393].f1[0] += localDW->Derivative2 * rtP.y[393];
  localDW->Traj_loc[394].f1[0] += localDW->Derivative2 * rtP.y[394];
  localDW->Traj_loc[395].f1[0] += localDW->Derivative2 * rtP.y[395];
  localDW->Traj_loc[396].f1[0] += localDW->Derivative2 * rtP.y[396];
  localDW->Traj_loc[397].f1[0] += localDW->Derivative2 * rtP.y[397];
  localDW->Traj_loc[398].f1[0] += localDW->Derivative2 * rtP.y[398];
  localDW->Traj_loc[399].f1[0] += localDW->Derivative2 * rtP.y[399];
  localDW->Traj_loc[400].f1[0] += localDW->Derivative2 * rtP.y[400];
  localDW->Traj_loc[401].f1[0] += localDW->Derivative2 * rtP.y[401];
  localDW->Traj_loc[402].f1[0] += localDW->Derivative2 * rtP.y[402];
  localDW->Traj_loc[403].f1[0] += localDW->Derivative2 * rtP.y[403];
  localDW->Traj_loc[404].f1[0] += localDW->Derivative2 * rtP.y[404];
  localDW->Traj_loc[405].f1[0] += localDW->Derivative2 * rtP.y[405];
  localDW->Traj_loc[406].f1[0] += localDW->Derivative2 * rtP.y[406];
  localDW->Traj_loc[407].f1[0] += localDW->Derivative2 * rtP.y[407];
  localDW->Traj_loc[408].f1[0] += localDW->Derivative2 * rtP.y[408];
  localDW->Traj_loc[409].f1[0] += localDW->Derivative2 * rtP.y[409];
  localDW->Traj_loc[410].f1[0] += localDW->Derivative2 * rtP.y[410];
  localDW->Traj_loc[411].f1[0] += localDW->Derivative2 * rtP.y[411];
  localDW->Traj_loc[412].f1[0] += localDW->Derivative2 * rtP.y[412];
  localDW->Traj_loc[413].f1[0] += localDW->Derivative2 * rtP.y[413];
  localDW->Traj_loc[414].f1[0] += localDW->Derivative2 * rtP.y[414];
  localDW->Traj_loc[415].f1[0] += localDW->Derivative2 * rtP.y[415];
  localDW->Traj_loc[416].f1[0] += localDW->Derivative2 * rtP.y[416];
  localDW->Traj_loc[417].f1[0] += localDW->Derivative2 * rtP.y[417];
  localDW->Traj_loc[418].f1[0] += localDW->Derivative2 * rtP.y[418];
  localDW->Traj_loc[419].f1[0] += localDW->Derivative2 * rtP.y[419];
  localDW->Traj_loc[420].f1[0] += localDW->Derivative2 * rtP.y[420];
  localDW->Traj_loc[421].f1[0] += localDW->Derivative2 * rtP.y[421];
  localDW->Traj_loc[422].f1[0] += localDW->Derivative2 * rtP.y[422];
  localDW->Traj_loc[423].f1[0] += localDW->Derivative2 * rtP.y[423];
  localDW->Traj_loc[424].f1[0] += localDW->Derivative2 * rtP.y[424];
  localDW->Traj_loc[425].f1[0] += localDW->Derivative2 * rtP.y[425];
  localDW->Traj_loc[426].f1[0] += localDW->Derivative2 * rtP.y[426];
  localDW->Traj_loc[427].f1[0] += localDW->Derivative2 * rtP.y[427];
  localDW->Traj_loc[428].f1[0] += localDW->Derivative2 * rtP.y[428];
  localDW->Traj_loc[429].f1[0] += localDW->Derivative2 * rtP.y[429];
  localDW->Traj_loc[430].f1[0] += localDW->Derivative2 * rtP.y[430];
  localDW->Traj_loc[431].f1[0] += localDW->Derivative2 * rtP.y[431];
  localDW->Traj_loc[432].f1[0] += localDW->Derivative2 * rtP.y[432];
  localDW->Traj_loc[433].f1[0] += localDW->Derivative2 * rtP.y[433];
  localDW->Traj_loc[434].f1[0] += localDW->Derivative2 * rtP.y[434];
  localDW->Traj_loc[435].f1[0] += localDW->Derivative2 * rtP.y[435];
  localDW->Traj_loc[436].f1[0] += localDW->Derivative2 * rtP.y[436];
  localDW->Traj_loc[437].f1[0] += localDW->Derivative2 * rtP.y[437];
  localDW->Traj_loc[438].f1[0] += localDW->Derivative2 * rtP.y[438];
  localDW->Traj_loc[439].f1[0] += localDW->Derivative2 * rtP.y[439];
  localDW->Traj_loc[440].f1[0] += localDW->Derivative2 * rtP.y[440];
  localDW->Traj_loc[441].f1[0] += localDW->Derivative2 * rtP.y[441];
  localDW->Traj_loc[442].f1[0] += localDW->Derivative2 * rtP.y[442];
  localDW->Traj_loc[443].f1[0] += localDW->Derivative2 * rtP.y[443];
  localDW->Traj_loc[444].f1[0] += localDW->Derivative2 * rtP.y[444];
  localDW->Traj_loc[445].f1[0] += localDW->Derivative2 * rtP.y[445];
  localDW->Traj_loc[446].f1[0] += localDW->Derivative2 * rtP.y[446];
  localDW->Traj_loc[447].f1[0] += localDW->Derivative2 * rtP.y[447];
  localDW->Traj_loc[448].f1[0] += localDW->Derivative2 * rtP.y[448];
  localDW->Traj_loc[449].f1[0] += localDW->Derivative2 * rtP.y[449];
  localDW->Traj_loc[450].f1[0] += localDW->Derivative2 * rtP.y[450];
  localDW->Traj_loc[451].f1[0] += localDW->Derivative2 * rtP.y[451];
  localDW->Traj_loc[452].f1[0] += localDW->Derivative2 * rtP.y[452];
  localDW->Traj_loc[453].f1[0] += localDW->Derivative2 * rtP.y[453];
  localDW->Traj_loc[454].f1[0] += localDW->Derivative2 * rtP.y[454];
  localDW->Traj_loc[455].f1[0] += localDW->Derivative2 * rtP.y[455];
  localDW->Traj_loc[456].f1[0] += localDW->Derivative2 * rtP.y[456];
  localDW->Traj_loc[457].f1[0] += localDW->Derivative2 * rtP.y[457];
  localDW->Traj_loc[458].f1[0] += localDW->Derivative2 * rtP.y[458];
  localDW->Traj_loc[459].f1[0] += localDW->Derivative2 * rtP.y[459];
  localDW->Traj_loc[460].f1[0] += localDW->Derivative2 * rtP.y[460];
  localDW->Traj_loc[461].f1[0] += localDW->Derivative2 * rtP.y[461];
  localDW->Traj_loc[462].f1[0] += localDW->Derivative2 * rtP.y[462];
  localDW->Traj_loc[463].f1[0] += localDW->Derivative2 * rtP.y[463];
  localDW->Traj_loc[464].f1[0] += localDW->Derivative2 * rtP.y[464];
  localDW->Traj_loc[465].f1[0] += localDW->Derivative2 * rtP.y[465];
  localDW->Traj_loc[466].f1[0] += localDW->Derivative2 * rtP.y[466];
  localDW->Traj_loc[467].f1[0] += localDW->Derivative2 * rtP.y[467];
  localDW->Traj_loc[468].f1[0] += localDW->Derivative2 * rtP.y[468];
  localDW->Traj_loc[469].f1[0] += localDW->Derivative2 * rtP.y[469];
  localDW->Traj_loc[470].f1[0] += localDW->Derivative2 * rtP.y[470];
  localDW->Traj_loc[471].f1[0] += localDW->Derivative2 * rtP.y[471];
  localDW->Traj_loc[472].f1[0] += localDW->Derivative2 * rtP.y[472];
  localDW->Traj_loc[473].f1[0] += localDW->Derivative2 * rtP.y[473];
  localDW->Traj_loc[474].f1[0] += localDW->Derivative2 * rtP.y[474];
  localDW->Traj_loc[475].f1[0] += localDW->Derivative2 * rtP.y[475];
  localDW->Traj_loc[476].f1[0] += localDW->Derivative2 * rtP.y[476];
  localDW->Traj_loc[477].f1[0] += localDW->Derivative2 * rtP.y[477];
  localDW->Traj_loc[478].f1[0] += localDW->Derivative2 * rtP.y[478];
  localDW->Traj_loc[479].f1[0] += localDW->Derivative2 * rtP.y[479];
  localDW->Traj_loc[480].f1[0] += localDW->Derivative2 * rtP.y[480];
  localDW->Traj_loc[481].f1[0] += localDW->Derivative2 * rtP.y[481];
  localDW->Traj_loc[482].f1[0] += localDW->Derivative2 * rtP.y[482];
  localDW->Traj_loc[483].f1[0] += localDW->Derivative2 * rtP.y[483];
  localDW->Traj_loc[484].f1[0] += localDW->Derivative2 * rtP.y[484];
  localDW->Traj_loc[485].f1[0] += localDW->Derivative2 * rtP.y[485];
  localDW->Traj_loc[486].f1[0] += localDW->Derivative2 * rtP.y[486];
  localDW->Traj_loc[487].f1[0] += localDW->Derivative2 * rtP.y[487];
  localDW->Traj_loc[488].f1[0] += localDW->Derivative2 * rtP.y[488];
  localDW->Traj_loc[489].f1[0] += localDW->Derivative2 * rtP.y[489];
  localDW->Traj_loc[490].f1[0] += localDW->Derivative2 * rtP.y[490];
  localDW->Traj_loc[491].f1[0] += localDW->Derivative2 * rtP.y[491];
  localDW->Traj_loc[492].f1[0] += localDW->Derivative2 * rtP.y[492];
  localDW->Traj_loc[493].f1[0] += localDW->Derivative2 * rtP.y[493];
  localDW->Traj_loc[494].f1[0] += localDW->Derivative2 * rtP.y[494];
  localDW->Traj_loc[495].f1[0] += localDW->Derivative2 * rtP.y[495];
  localDW->Traj_loc[496].f1[0] += localDW->Derivative2 * rtP.y[496];
  localDW->Traj_loc[497].f1[0] += localDW->Derivative2 * rtP.y[497];
  localDW->Traj_loc[498].f1[0] += localDW->Derivative2 * rtP.y[498];
  localDW->Traj_loc[499].f1[0] += localDW->Derivative2 * rtP.y[499];
  localDW->Traj_loc[500].f1[0] += localDW->Derivative2 * rtP.y[500];
  localDW->Traj_loc[501].f1[0] += localDW->Derivative2 * rtP.y[501];
  localDW->Traj_loc[502].f1[0] += localDW->Derivative2 * rtP.y[502];
  localDW->Traj_loc[503].f1[0] += localDW->Derivative2 * rtP.y[503];
  localDW->Traj_loc[504].f1[0] += localDW->Derivative2 * rtP.y[504];
  localDW->Traj_loc[505].f1[0] += localDW->Derivative2 * rtP.y[505];
  localDW->Traj_loc[506].f1[0] += localDW->Derivative2 * rtP.y[506];
  localDW->Traj_loc[507].f1[0] += localDW->Derivative2 * rtP.y[507];
  localDW->Traj_loc[508].f1[0] += localDW->Derivative2 * rtP.y[508];
  localDW->Traj_loc[509].f1[0] += localDW->Derivative2 * rtP.y[509];
  localDW->Traj_loc[510].f1[0] += localDW->Derivative2 * rtP.y[510];
  localDW->Traj_loc[511].f1[0] += localDW->Derivative2 * rtP.y[511];
  localDW->Traj_loc[512].f1[0] += localDW->Derivative2 * rtP.y[512];
  localDW->Traj_loc[513].f1[0] += localDW->Derivative2 * rtP.y[513];
  localDW->Traj_loc[514].f1[0] += localDW->Derivative2 * rtP.y[514];
  localDW->Traj_loc[515].f1[0] += localDW->Derivative2 * rtP.y[515];
  localDW->Traj_loc[516].f1[0] += localDW->Derivative2 * rtP.y[516];
  localDW->Traj_loc[517].f1[0] += localDW->Derivative2 * rtP.y[517];
  localDW->Traj_loc[518].f1[0] += localDW->Derivative2 * rtP.y[518];
  localDW->Traj_loc[519].f1[0] += localDW->Derivative2 * rtP.y[519];
  localDW->Traj_loc[520].f1[0] += localDW->Derivative2 * rtP.y[520];
  localDW->Traj_loc[521].f1[0] += localDW->Derivative2 * rtP.y[521];
  localDW->Traj_loc[522].f1[0] += localDW->Derivative2 * rtP.y[522];
  localDW->Traj_loc[523].f1[0] += localDW->Derivative2 * rtP.y[523];
  localDW->Traj_loc[524].f1[0] += localDW->Derivative2 * rtP.y[524];
  localDW->Traj_loc[525].f1[0] += localDW->Derivative2 * rtP.y[525];
  localDW->Traj_loc[526].f1[0] += localDW->Derivative2 * rtP.y[526];
  localDW->Traj_loc[527].f1[0] += localDW->Derivative2 * rtP.y[527];
  localDW->Traj_loc[528].f1[0] += localDW->Derivative2 * rtP.y[528];
  localDW->Traj_loc[529].f1[0] += localDW->Derivative2 * rtP.y[529];
  localDW->Traj_loc[530].f1[0] += localDW->Derivative2 * rtP.y[530];
  localDW->Traj_loc[531].f1[0] += localDW->Derivative2 * rtP.y[531];
  localDW->Traj_loc[532].f1[0] += localDW->Derivative2 * rtP.y[532];
  localDW->Traj_loc[533].f1[0] += localDW->Derivative2 * rtP.y[533];
  localDW->Traj_loc[534].f1[0] += localDW->Derivative2 * rtP.y[534];
  localDW->Traj_loc[535].f1[0] += localDW->Derivative2 * rtP.y[535];
  localDW->Traj_loc[536].f1[0] += localDW->Derivative2 * rtP.y[536];
  localDW->Traj_loc[537].f1[0] += localDW->Derivative2 * rtP.y[537];
  localDW->Traj_loc[538].f1[0] += localDW->Derivative2 * rtP.y[538];
  localDW->Traj_loc[539].f1[0] += localDW->Derivative2 * rtP.y[539];
  localDW->Traj_loc[540].f1[0] += localDW->Derivative2 * rtP.y[540];
  localDW->Traj_loc[541].f1[0] += localDW->Derivative2 * rtP.y[541];
  localDW->Traj_loc[542].f1[0] += localDW->Derivative2 * rtP.y[542];
  localDW->Traj_loc[543].f1[0] += localDW->Derivative2 * rtP.y[543];
  localDW->Traj_loc[544].f1[0] += localDW->Derivative2 * rtP.y[544];
  localDW->Traj_loc[545].f1[0] += localDW->Derivative2 * rtP.y[545];
  localDW->Traj_loc[546].f1[0] += localDW->Derivative2 * rtP.y[546];
  localDW->Traj_loc[547].f1[0] += localDW->Derivative2 * rtP.y[547];
  localDW->Traj_loc[548].f1[0] += localDW->Derivative2 * rtP.y[548];
  localDW->Traj_loc[549].f1[0] += localDW->Derivative2 * rtP.y[549];
  localDW->Traj_loc[550].f1[0] += localDW->Derivative2 * rtP.y[550];
  localDW->Traj_loc[551].f1[0] += localDW->Derivative2 * rtP.y[551];
  localDW->Traj_loc[552].f1[0] += localDW->Derivative2 * rtP.y[552];
  localDW->Traj_loc[553].f1[0] += localDW->Derivative2 * rtP.y[553];
  localDW->Traj_loc[554].f1[0] += localDW->Derivative2 * rtP.y[554];
  localDW->Traj_loc[555].f1[0] += localDW->Derivative2 * rtP.y[555];
  localDW->Traj_loc[556].f1[0] += localDW->Derivative2 * rtP.y[556];
  localDW->Traj_loc[557].f1[0] += localDW->Derivative2 * rtP.y[557];
  localDW->Traj_loc[558].f1[0] += localDW->Derivative2 * rtP.y[558];
  localDW->Traj_loc[559].f1[0] += localDW->Derivative2 * rtP.y[559];
  localDW->Traj_loc[560].f1[0] += localDW->Derivative2 * rtP.y[560];
  localDW->Traj_loc[561].f1[0] += localDW->Derivative2 * rtP.y[561];
  localDW->Traj_loc[562].f1[0] += localDW->Derivative2 * rtP.y[562];
  localDW->Traj_loc[563].f1[0] += localDW->Derivative2 * rtP.y[563];
  localDW->Traj_loc[564].f1[0] += localDW->Derivative2 * rtP.y[564];
  localDW->Traj_loc[565].f1[0] += localDW->Derivative2 * rtP.y[565];
  localDW->Traj_loc[566].f1[0] += localDW->Derivative2 * rtP.y[566];
  localDW->Traj_loc[567].f1[0] += localDW->Derivative2 * rtP.y[567];
  localDW->Traj_loc[568].f1[0] += localDW->Derivative2 * rtP.y[568];
  localDW->Traj_loc[569].f1[0] += localDW->Derivative2 * rtP.y[569];
  localDW->Traj_loc[570].f1[0] += localDW->Derivative2 * rtP.y[570];
  localDW->Traj_loc[571].f1[0] += localDW->Derivative2 * rtP.y[571];
  localDW->Traj_loc[572].f1[0] += localDW->Derivative2 * rtP.y[572];
  localDW->Traj_loc[573].f1[0] += localDW->Derivative2 * rtP.y[573];
  localDW->Traj_loc[574].f1[0] += localDW->Derivative2 * rtP.y[574];
  localDW->Traj_loc[575].f1[0] += localDW->Derivative2 * rtP.y[575];
  localDW->Traj_loc[576].f1[0] += localDW->Derivative2 * rtP.y[576];
  localDW->Traj_loc[577].f1[0] += localDW->Derivative2 * rtP.y[577];
  localDW->Traj_loc[578].f1[0] += localDW->Derivative2 * rtP.y[578];
  localDW->Traj_loc[579].f1[0] += localDW->Derivative2 * rtP.y[579];
  localDW->Traj_loc[580].f1[0] += localDW->Derivative2 * rtP.y[580];
  localDW->Traj_loc[581].f1[0] += localDW->Derivative2 * rtP.y[581];
  localDW->Traj_loc[582].f1[0] += localDW->Derivative2 * rtP.y[582];
  localDW->Traj_loc[583].f1[0] += localDW->Derivative2 * rtP.y[583];
  localDW->Traj_loc[584].f1[0] += localDW->Derivative2 * rtP.y[584];
  localDW->Traj_loc[585].f1[0] += localDW->Derivative2 * rtP.y[585];
  localDW->Traj_loc[586].f1[0] += localDW->Derivative2 * rtP.y[586];
  localDW->Traj_loc[587].f1[0] += localDW->Derivative2 * rtP.y[587];
  localDW->Traj_loc[588].f1[0] += localDW->Derivative2 * rtP.y[588];
  localDW->Traj_loc[589].f1[0] += localDW->Derivative2 * rtP.y[589];
  localDW->Traj_loc[590].f1[0] += localDW->Derivative2 * rtP.y[590];
  localDW->Traj_loc[591].f1[0] += localDW->Derivative2 * rtP.y[591];
  localDW->Traj_loc[592].f1[0] += localDW->Derivative2 * rtP.y[592];
  localDW->Traj_loc[593].f1[0] += localDW->Derivative2 * rtP.y[593];
  localDW->Traj_loc[594].f1[0] += localDW->Derivative2 * rtP.y[594];
  localDW->Traj_loc[595].f1[0] += localDW->Derivative2 * rtP.y[595];
  localDW->Traj_loc[596].f1[0] += localDW->Derivative2 * rtP.y[596];
  localDW->Traj_loc[597].f1[0] += localDW->Derivative2 * rtP.y[597];
  localDW->Traj_loc[598].f1[0] += localDW->Derivative2 * rtP.y[598];
  localDW->Traj_loc[599].f1[0] += localDW->Derivative2 * rtP.y[599];
  localDW->Traj_loc[600].f1[0] += localDW->Derivative2 * rtP.y[600];
  localDW->Traj_loc[601].f1[0] += localDW->Derivative2 * rtP.y[601];
  localDW->Traj_loc[602].f1[0] += localDW->Derivative2 * rtP.y[602];
  localDW->Traj_loc[603].f1[0] += localDW->Derivative2 * rtP.y[603];
  localDW->Traj_loc[604].f1[0] += localDW->Derivative2 * rtP.y[604];
  localDW->Traj_loc[605].f1[0] += localDW->Derivative2 * rtP.y[605];
  localDW->Traj_loc[606].f1[0] += localDW->Derivative2 * rtP.y[606];
  localDW->Traj_loc[607].f1[0] += localDW->Derivative2 * rtP.y[607];
  localDW->Traj_loc[608].f1[0] += localDW->Derivative2 * rtP.y[608];
  localDW->Traj_loc[609].f1[0] += localDW->Derivative2 * rtP.y[609];
  localDW->Traj_loc[610].f1[0] += localDW->Derivative2 * rtP.y[610];
  localDW->Traj_loc[611].f1[0] += localDW->Derivative2 * rtP.y[611];
  localDW->Traj_loc[612].f1[0] += localDW->Derivative2 * rtP.y[612];
  localDW->Traj_loc[613].f1[0] += localDW->Derivative2 * rtP.y[613];
  localDW->Traj_loc[614].f1[0] += localDW->Derivative2 * rtP.y[614];
  localDW->Traj_loc[615].f1[0] += localDW->Derivative2 * rtP.y[615];
  localDW->Traj_loc[616].f1[0] += localDW->Derivative2 * rtP.y[616];
  localDW->Traj_loc[617].f1[0] += localDW->Derivative2 * rtP.y[617];
  localDW->Traj_loc[618].f1[0] += localDW->Derivative2 * rtP.y[618];
  localDW->Traj_loc[619].f1[0] += localDW->Derivative2 * rtP.y[619];
  localDW->Traj_loc[620].f1[0] += localDW->Derivative2 * rtP.y[620];
  localDW->Traj_loc[621].f1[0] += localDW->Derivative2 * rtP.y[621];
  localDW->Traj_loc[622].f1[0] += localDW->Derivative2 * rtP.y[622];
  localDW->Traj_loc[623].f1[0] += localDW->Derivative2 * rtP.y[623];
  localDW->Traj_loc[624].f1[0] += localDW->Derivative2 * rtP.y[624];
  localDW->Traj_loc[625].f1[0] += localDW->Derivative2 * rtP.y[625];
  localDW->Traj_loc[626].f1[0] += localDW->Derivative2 * rtP.y[626];
  localDW->Traj_loc[627].f1[0] += localDW->Derivative2 * rtP.y[627];
  localDW->Traj_loc[628].f1[0] += localDW->Derivative2 * rtP.y[628];
  localDW->Traj_loc[629].f1[0] += localDW->Derivative2 * rtP.y[629];
  localDW->Traj_loc[630].f1[0] += localDW->Derivative2 * rtP.y[630];
  localDW->Traj_loc[631].f1[0] += localDW->Derivative2 * rtP.y[631];
  localDW->Traj_loc[632].f1[0] += localDW->Derivative2 * rtP.y[632];
  localDW->Traj_loc[633].f1[0] += localDW->Derivative2 * rtP.y[633];
  localDW->Traj_loc[634].f1[0] += localDW->Derivative2 * rtP.y[634];
  localDW->Traj_loc[635].f1[0] += localDW->Derivative2 * rtP.y[635];
  localDW->Traj_loc[636].f1[0] += localDW->Derivative2 * rtP.y[636];
  localDW->Traj_loc[637].f1[0] += localDW->Derivative2 * rtP.y[637];
  localDW->Traj_loc[638].f1[0] += localDW->Derivative2 * rtP.y[638];
  localDW->Traj_loc[639].f1[0] += localDW->Derivative2 * rtP.y[639];
  localDW->Traj_loc[640].f1[0] += localDW->Derivative2 * rtP.y[640];
  localDW->Traj_loc[641].f1[0] += localDW->Derivative2 * rtP.y[641];
  localDW->Traj_loc[642].f1[0] += localDW->Derivative2 * rtP.y[642];
  localDW->Traj_loc[643].f1[0] += localDW->Derivative2 * rtP.y[643];
  localDW->Traj_loc[644].f1[0] += localDW->Derivative2 * rtP.y[644];
  localDW->Traj_loc[645].f1[0] += localDW->Derivative2 * rtP.y[645];
  localDW->Traj_loc[646].f1[0] += localDW->Derivative2 * rtP.y[646];
  localDW->Traj_loc[647].f1[0] += localDW->Derivative2 * rtP.y[647];
  localDW->Traj_loc[648].f1[0] += localDW->Derivative2 * rtP.y[648];
  localDW->Traj_loc[649].f1[0] += localDW->Derivative2 * rtP.y[649];
  localDW->Traj_loc[650].f1[0] += localDW->Derivative2 * rtP.y[650];
  localDW->Traj_loc[651].f1[0] += localDW->Derivative2 * rtP.y[651];
  localDW->Traj_loc[652].f1[0] += localDW->Derivative2 * rtP.y[652];
  localDW->Traj_loc[653].f1[0] += localDW->Derivative2 * rtP.y[653];
  localDW->Traj_loc[654].f1[0] += localDW->Derivative2 * rtP.y[654];
  localDW->Traj_loc[655].f1[0] += localDW->Derivative2 * rtP.y[655];
  localDW->Traj_loc[656].f1[0] += localDW->Derivative2 * rtP.y[656];
  localDW->Traj_loc[657].f1[0] += localDW->Derivative2 * rtP.y[657];
  localDW->Traj_loc[658].f1[0] += localDW->Derivative2 * rtP.y[658];
  localDW->Traj_loc[659].f1[0] += localDW->Derivative2 * rtP.y[659];
  localDW->Traj_loc[660].f1[0] += localDW->Derivative2 * rtP.y[660];
  localDW->Traj_loc[661].f1[0] += localDW->Derivative2 * rtP.y[661];
  localDW->Traj_loc[662].f1[0] += localDW->Derivative2 * rtP.y[662];
  localDW->Traj_loc[663].f1[0] += localDW->Derivative2 * rtP.y[663];
  localDW->Traj_loc[664].f1[0] += localDW->Derivative2 * rtP.y[664];
  localDW->Traj_loc[665].f1[0] += localDW->Derivative2 * rtP.y[665];
  localDW->Traj_loc[666].f1[0] += localDW->Derivative2 * rtP.y[666];
  localDW->Traj_loc[667].f1[0] += localDW->Derivative2 * rtP.y[667];
  localDW->Traj_loc[668].f1[0] += localDW->Derivative2 * rtP.y[668];
  localDW->Traj_loc[669].f1[0] += localDW->Derivative2 * rtP.y[669];
  localDW->Traj_loc[670].f1[0] += localDW->Derivative2 * rtP.y[670];
  localDW->Traj_loc[671].f1[0] += localDW->Derivative2 * rtP.y[671];
  localDW->Traj_loc[672].f1[0] += localDW->Derivative2 * rtP.y[672];
  localDW->Traj_loc[673].f1[0] += localDW->Derivative2 * rtP.y[673];
  localDW->Traj_loc[674].f1[0] += localDW->Derivative2 * rtP.y[674];
  localDW->Traj_loc[675].f1[0] += localDW->Derivative2 * rtP.y[675];
  localDW->Traj_loc[676].f1[0] += localDW->Derivative2 * rtP.y[676];
  localDW->Traj_loc[677].f1[0] += localDW->Derivative2 * rtP.y[677];
  localDW->Traj_loc[678].f1[0] += localDW->Derivative2 * rtP.y[678];
  localDW->Traj_loc[679].f1[0] += localDW->Derivative2 * rtP.y[679];
  localDW->Traj_loc[680].f1[0] += localDW->Derivative2 * rtP.y[680];
  localDW->Traj_loc[681].f1[0] += localDW->Derivative2 * rtP.y[681];
  localDW->Traj_loc[682].f1[0] += localDW->Derivative2 * rtP.y[682];
  localDW->Traj_loc[683].f1[0] += localDW->Derivative2 * rtP.y[683];
  localDW->Traj_loc[684].f1[0] += localDW->Derivative2 * rtP.y[684];
  localDW->Traj_loc[685].f1[0] += localDW->Derivative2 * rtP.y[685];
  localDW->Traj_loc[686].f1[0] += localDW->Derivative2 * rtP.y[686];
  localDW->Traj_loc[687].f1[0] += localDW->Derivative2 * rtP.y[687];
  localDW->Traj_loc[688].f1[0] += localDW->Derivative2 * rtP.y[688];
  localDW->Traj_loc[689].f1[0] += localDW->Derivative2 * rtP.y[689];
  localDW->Traj_loc[690].f1[0] += localDW->Derivative2 * rtP.y[690];
  localDW->Traj_loc[691].f1[0] += localDW->Derivative2 * rtP.y[691];
  localDW->Traj_loc[692].f1[0] += localDW->Derivative2 * rtP.y[692];
  localDW->Traj_loc[693].f1[0] += localDW->Derivative2 * rtP.y[693];
  localDW->Traj_loc[694].f1[0] += localDW->Derivative2 * rtP.y[694];
  localDW->Traj_loc[695].f1[0] += localDW->Derivative2 * rtP.y[695];
  localDW->Traj_loc[696].f1[0] += localDW->Derivative2 * rtP.y[696];
  localDW->Traj_loc[697].f1[0] += localDW->Derivative2 * rtP.y[697];
  localDW->Traj_loc[698].f1[0] += localDW->Derivative2 * rtP.y[698];
  localDW->Traj_loc[699].f1[0] += localDW->Derivative2 * rtP.y[699];
  localDW->Traj_loc[700].f1[0] += localDW->Derivative2 * rtP.y[700];
  localDW->Traj_loc[701].f1[0] += localDW->Derivative2 * rtP.y[701];
  localDW->Traj_loc[702].f1[0] += localDW->Derivative2 * rtP.y[702];
  localDW->Traj_loc[703].f1[0] += localDW->Derivative2 * rtP.y[703];
  localDW->Traj_loc[704].f1[0] += localDW->Derivative2 * rtP.y[704];
  localDW->Traj_loc[705].f1[0] += localDW->Derivative2 * rtP.y[705];
  localDW->Traj_loc[706].f1[0] += localDW->Derivative2 * rtP.y[706];
  localDW->Traj_loc[707].f1[0] += localDW->Derivative2 * rtP.y[707];
  localDW->Traj_loc[708].f1[0] += localDW->Derivative2 * rtP.y[708];
  localDW->Traj_loc[709].f1[0] += localDW->Derivative2 * rtP.y[709];
  localDW->Traj_loc[710].f1[0] += localDW->Derivative2 * rtP.y[710];
  localDW->Traj_loc[711].f1[0] += localDW->Derivative2 * rtP.y[711];
  localDW->Traj_loc[712].f1[0] += localDW->Derivative2 * rtP.y[712];
  localDW->Traj_loc[713].f1[0] += localDW->Derivative2 * rtP.y[713];
  localDW->Traj_loc[714].f1[0] += localDW->Derivative2 * rtP.y[714];
  localDW->Traj_loc[715].f1[0] += localDW->Derivative2 * rtP.y[715];
  localDW->Traj_loc[716].f1[0] += localDW->Derivative2 * rtP.y[716];
  localDW->Traj_loc[717].f1[0] += localDW->Derivative2 * rtP.y[717];
  localDW->Traj_loc[718].f1[0] += localDW->Derivative2 * rtP.y[718];
  localDW->Traj_loc[719].f1[0] += localDW->Derivative2 * rtP.y[719];
  localDW->Traj_loc[720].f1[0] += localDW->Derivative2 * rtP.y[720];
  localDW->Traj_loc[721].f1[0] += localDW->Derivative2 * rtP.y[721];
  localDW->Traj_loc[722].f1[0] += localDW->Derivative2 * rtP.y[722];
  localDW->Traj_loc[723].f1[0] += localDW->Derivative2 * rtP.y[723];
  localDW->Traj_loc[724].f1[0] += localDW->Derivative2 * rtP.y[724];
  localDW->Traj_loc[725].f1[0] += localDW->Derivative2 * rtP.y[725];
  localDW->Traj_loc[726].f1[0] += localDW->Derivative2 * rtP.y[726];
  localDW->Traj_loc[727].f1[0] += localDW->Derivative2 * rtP.y[727];
  localDW->Traj_loc[728].f1[0] += localDW->Derivative2 * rtP.y[728];
  localDW->Traj_loc[729].f1[0] += localDW->Derivative2 * rtP.y[729];
  localDW->Traj_loc[730].f1[0] += localDW->Derivative2 * rtP.y[730];
  localDW->Traj_loc[731].f1[0] += localDW->Derivative2 * rtP.y[731];
  localDW->Traj_loc[732].f1[0] += localDW->Derivative2 * rtP.y[732];
  localDW->Traj_loc[733].f1[0] += localDW->Derivative2 * rtP.y[733];
  localDW->Traj_loc[734].f1[0] += localDW->Derivative2 * rtP.y[734];
  localDW->Traj_loc[735].f1[0] += localDW->Derivative2 * rtP.y[735];
  localDW->Traj_loc[736].f1[0] += localDW->Derivative2 * rtP.y[736];
  localDW->Traj_loc[737].f1[0] += localDW->Derivative2 * rtP.y[737];
  localDW->Traj_loc[738].f1[0] += localDW->Derivative2 * rtP.y[738];
  localDW->Traj_loc[739].f1[0] += localDW->Derivative2 * rtP.y[739];
  localDW->Traj_loc[740].f1[0] += localDW->Derivative2 * rtP.y[740];
  localDW->Traj_loc[741].f1[0] += localDW->Derivative2 * rtP.y[741];
  localDW->Traj_loc[742].f1[0] += localDW->Derivative2 * rtP.y[742];
  localDW->Traj_loc[743].f1[0] += localDW->Derivative2 * rtP.y[743];
  localDW->Traj_loc[744].f1[0] += localDW->Derivative2 * rtP.y[744];
  localDW->Traj_loc[745].f1[0] += localDW->Derivative2 * rtP.y[745];
  localDW->Traj_loc[746].f1[0] += localDW->Derivative2 * rtP.y[746];
  localDW->Traj_loc[747].f1[0] += localDW->Derivative2 * rtP.y[747];
  localDW->Traj_loc[748].f1[0] += localDW->Derivative2 * rtP.y[748];
  localDW->Traj_loc[749].f1[0] += localDW->Derivative2 * rtP.y[749];
  localDW->Traj_loc[750].f1[0] += localDW->Derivative2 * rtP.y[750];
  localDW->Traj_loc[751].f1[0] += localDW->Derivative2 * rtP.y[751];
  localDW->Traj_loc[752].f1[0] += localDW->Derivative2 * rtP.y[752];
  localDW->Traj_loc[753].f1[0] += localDW->Derivative2 * rtP.y[753];
  localDW->Traj_loc[754].f1[0] += localDW->Derivative2 * rtP.y[754];
  localDW->Traj_loc[755].f1[0] += localDW->Derivative2 * rtP.y[755];
  localDW->Traj_loc[756].f1[0] += localDW->Derivative2 * rtP.y[756];
  localDW->Traj_loc[757].f1[0] += localDW->Derivative2 * rtP.y[757];
  localDW->Traj_loc[758].f1[0] += localDW->Derivative2 * rtP.y[758];
  localDW->Traj_loc[759].f1[0] += localDW->Derivative2 * rtP.y[759];
  localDW->Traj_loc[760].f1[0] += localDW->Derivative2 * rtP.y[760];
  localDW->Traj_loc[761].f1[0] += localDW->Derivative2 * rtP.y[761];
  localDW->Traj_loc[762].f1[0] += localDW->Derivative2 * rtP.y[762];
  localDW->Traj_loc[763].f1[0] += localDW->Derivative2 * rtP.y[763];
  localDW->Traj_loc[764].f1[0] += localDW->Derivative2 * rtP.y[764];
  localDW->Traj_loc[765].f1[0] += localDW->Derivative2 * rtP.y[765];
  localDW->Traj_loc[766].f1[0] += localDW->Derivative2 * rtP.y[766];
  localDW->Traj_loc[767].f1[0] += localDW->Derivative2 * rtP.y[767];
  localDW->Traj_loc[768].f1[0] += localDW->Derivative2 * rtP.y[768];
  localDW->Traj_loc[769].f1[0] += localDW->Derivative2 * rtP.y[769];
  localDW->Traj_loc[770].f1[0] += localDW->Derivative2 * rtP.y[770];
  localDW->Traj_loc[771].f1[0] += localDW->Derivative2 * rtP.y[771];
  localDW->Traj_loc[772].f1[0] += localDW->Derivative2 * rtP.y[772];
  localDW->Traj_loc[773].f1[0] += localDW->Derivative2 * rtP.y[773];
  localDW->Traj_loc[774].f1[0] += localDW->Derivative2 * rtP.y[774];
  localDW->Traj_loc[775].f1[0] += localDW->Derivative2 * rtP.y[775];
  localDW->Traj_loc[776].f1[0] += localDW->Derivative2 * rtP.y[776];
  localDW->Traj_loc[777].f1[0] += localDW->Derivative2 * rtP.y[777];
  localDW->Traj_loc[778].f1[0] += localDW->Derivative2 * rtP.y[778];
  localDW->Traj_loc[779].f1[0] += localDW->Derivative2 * rtP.y[779];
  localDW->Traj_loc[780].f1[0] += localDW->Derivative2 * rtP.y[780];
  localDW->Traj_loc[781].f1[0] += localDW->Derivative2 * rtP.y[781];
  localDW->Traj_loc[782].f1[0] += localDW->Derivative2 * rtP.y[782];
  localDW->Traj_loc[783].f1[0] += localDW->Derivative2 * rtP.y[783];
  localDW->Traj_loc[784].f1[0] += localDW->Derivative2 * rtP.y[784];
  localDW->Traj_loc[785].f1[0] += localDW->Derivative2 * rtP.y[785];
  localDW->Traj_loc[786].f1[0] += localDW->Derivative2 * rtP.y[786];
  localDW->Traj_loc[787].f1[0] += localDW->Derivative2 * rtP.y[787];
  localDW->Traj_loc[788].f1[0] += localDW->Derivative2 * rtP.y[788];
  localDW->Traj_loc[789].f1[0] += localDW->Derivative2 * rtP.y[789];
  localDW->Traj_loc[790].f1[0] += localDW->Derivative2 * rtP.y[790];
  localDW->Traj_loc[791].f1[0] += localDW->Derivative2 * rtP.y[791];
  localDW->Traj_loc[792].f1[0] += localDW->Derivative2 * rtP.y[792];
  localDW->Traj_loc[793].f1[0] += localDW->Derivative2 * rtP.y[793];
  localDW->Traj_loc[794].f1[0] += localDW->Derivative2 * rtP.y[794];
  localDW->Traj_loc[795].f1[0] += localDW->Derivative2 * rtP.y[795];
  localDW->Traj_loc[796].f1[0] += localDW->Derivative2 * rtP.y[796];
  localDW->Traj_loc[797].f1[0] += localDW->Derivative2 * rtP.y[797];
  localDW->Traj_loc[798].f1[0] += localDW->Derivative2 * rtP.y[798];
  localDW->Traj_loc[799].f1[0] += localDW->Derivative2 * rtP.y[799];
  localDW->Traj_loc[800].f1[0] += localDW->Derivative2 * rtP.y[800];
  localDW->Traj_loc[0].f1[1] = localDW->Derivative2 * rtP.x[0];
  localDW->Traj_loc[1].f1[1] = localDW->Derivative2 * rtP.x[1];
  localDW->Traj_loc[2].f1[1] = localDW->Derivative2 * rtP.x[2];
  localDW->Traj_loc[3].f1[1] = localDW->Derivative2 * rtP.x[3];
  localDW->Traj_loc[4].f1[1] = localDW->Derivative2 * rtP.x[4];
  localDW->Traj_loc[5].f1[1] = localDW->Derivative2 * rtP.x[5];
  localDW->Traj_loc[6].f1[1] = localDW->Derivative2 * rtP.x[6];
  localDW->Traj_loc[7].f1[1] = localDW->Derivative2 * rtP.x[7];
  localDW->Traj_loc[8].f1[1] = localDW->Derivative2 * rtP.x[8];
  localDW->Traj_loc[9].f1[1] = localDW->Derivative2 * rtP.x[9];
  localDW->Traj_loc[10].f1[1] = localDW->Derivative2 * rtP.x[10];
  localDW->Traj_loc[11].f1[1] = localDW->Derivative2 * rtP.x[11];
  localDW->Traj_loc[12].f1[1] = localDW->Derivative2 * rtP.x[12];
  localDW->Traj_loc[13].f1[1] = localDW->Derivative2 * rtP.x[13];
  localDW->Traj_loc[14].f1[1] = localDW->Derivative2 * rtP.x[14];
  localDW->Traj_loc[15].f1[1] = localDW->Derivative2 * rtP.x[15];
  localDW->Traj_loc[16].f1[1] = localDW->Derivative2 * rtP.x[16];
  localDW->Traj_loc[17].f1[1] = localDW->Derivative2 * rtP.x[17];
  localDW->Traj_loc[18].f1[1] = localDW->Derivative2 * rtP.x[18];
  localDW->Traj_loc[19].f1[1] = localDW->Derivative2 * rtP.x[19];
  localDW->Traj_loc[20].f1[1] = localDW->Derivative2 * rtP.x[20];
  localDW->Traj_loc[21].f1[1] = localDW->Derivative2 * rtP.x[21];
  localDW->Traj_loc[22].f1[1] = localDW->Derivative2 * rtP.x[22];
  localDW->Traj_loc[23].f1[1] = localDW->Derivative2 * rtP.x[23];
  localDW->Traj_loc[24].f1[1] = localDW->Derivative2 * rtP.x[24];
  localDW->Traj_loc[25].f1[1] = localDW->Derivative2 * rtP.x[25];
  localDW->Traj_loc[26].f1[1] = localDW->Derivative2 * rtP.x[26];
  localDW->Traj_loc[27].f1[1] = localDW->Derivative2 * rtP.x[27];
  localDW->Traj_loc[28].f1[1] = localDW->Derivative2 * rtP.x[28];
  localDW->Traj_loc[29].f1[1] = localDW->Derivative2 * rtP.x[29];
  localDW->Traj_loc[30].f1[1] = localDW->Derivative2 * rtP.x[30];
  localDW->Traj_loc[31].f1[1] = localDW->Derivative2 * rtP.x[31];
  localDW->Traj_loc[32].f1[1] = localDW->Derivative2 * rtP.x[32];
  localDW->Traj_loc[33].f1[1] = localDW->Derivative2 * rtP.x[33];
  localDW->Traj_loc[34].f1[1] = localDW->Derivative2 * rtP.x[34];
  localDW->Traj_loc[35].f1[1] = localDW->Derivative2 * rtP.x[35];
  localDW->Traj_loc[36].f1[1] = localDW->Derivative2 * rtP.x[36];
  localDW->Traj_loc[37].f1[1] = localDW->Derivative2 * rtP.x[37];
  localDW->Traj_loc[38].f1[1] = localDW->Derivative2 * rtP.x[38];
  localDW->Traj_loc[39].f1[1] = localDW->Derivative2 * rtP.x[39];
  localDW->Traj_loc[40].f1[1] = localDW->Derivative2 * rtP.x[40];
  localDW->Traj_loc[41].f1[1] = localDW->Derivative2 * rtP.x[41];
  localDW->Traj_loc[42].f1[1] = localDW->Derivative2 * rtP.x[42];
  localDW->Traj_loc[43].f1[1] = localDW->Derivative2 * rtP.x[43];
  localDW->Traj_loc[44].f1[1] = localDW->Derivative2 * rtP.x[44];
  localDW->Traj_loc[45].f1[1] = localDW->Derivative2 * rtP.x[45];
  localDW->Traj_loc[46].f1[1] = localDW->Derivative2 * rtP.x[46];
  localDW->Traj_loc[47].f1[1] = localDW->Derivative2 * rtP.x[47];
  localDW->Traj_loc[48].f1[1] = localDW->Derivative2 * rtP.x[48];
  localDW->Traj_loc[49].f1[1] = localDW->Derivative2 * rtP.x[49];
  localDW->Traj_loc[50].f1[1] = localDW->Derivative2 * rtP.x[50];
  localDW->Traj_loc[51].f1[1] = localDW->Derivative2 * rtP.x[51];
  localDW->Traj_loc[52].f1[1] = localDW->Derivative2 * rtP.x[52];
  localDW->Traj_loc[53].f1[1] = localDW->Derivative2 * rtP.x[53];
  localDW->Traj_loc[54].f1[1] = localDW->Derivative2 * rtP.x[54];
  localDW->Traj_loc[55].f1[1] = localDW->Derivative2 * rtP.x[55];
  localDW->Traj_loc[56].f1[1] = localDW->Derivative2 * rtP.x[56];
  localDW->Traj_loc[57].f1[1] = localDW->Derivative2 * rtP.x[57];
  localDW->Traj_loc[58].f1[1] = localDW->Derivative2 * rtP.x[58];
  localDW->Traj_loc[59].f1[1] = localDW->Derivative2 * rtP.x[59];
  localDW->Traj_loc[60].f1[1] = localDW->Derivative2 * rtP.x[60];
  localDW->Traj_loc[61].f1[1] = localDW->Derivative2 * rtP.x[61];
  localDW->Traj_loc[62].f1[1] = localDW->Derivative2 * rtP.x[62];
  localDW->Traj_loc[63].f1[1] = localDW->Derivative2 * rtP.x[63];
  localDW->Traj_loc[64].f1[1] = localDW->Derivative2 * rtP.x[64];
  localDW->Traj_loc[65].f1[1] = localDW->Derivative2 * rtP.x[65];
  localDW->Traj_loc[66].f1[1] = localDW->Derivative2 * rtP.x[66];
  localDW->Traj_loc[67].f1[1] = localDW->Derivative2 * rtP.x[67];
  localDW->Traj_loc[68].f1[1] = localDW->Derivative2 * rtP.x[68];
  localDW->Traj_loc[69].f1[1] = localDW->Derivative2 * rtP.x[69];
  localDW->Traj_loc[70].f1[1] = localDW->Derivative2 * rtP.x[70];
  localDW->Traj_loc[71].f1[1] = localDW->Derivative2 * rtP.x[71];
  localDW->Traj_loc[72].f1[1] = localDW->Derivative2 * rtP.x[72];
  localDW->Traj_loc[73].f1[1] = localDW->Derivative2 * rtP.x[73];
  localDW->Traj_loc[74].f1[1] = localDW->Derivative2 * rtP.x[74];
  localDW->Traj_loc[75].f1[1] = localDW->Derivative2 * rtP.x[75];
  localDW->Traj_loc[76].f1[1] = localDW->Derivative2 * rtP.x[76];
  localDW->Traj_loc[77].f1[1] = localDW->Derivative2 * rtP.x[77];
  localDW->Traj_loc[78].f1[1] = localDW->Derivative2 * rtP.x[78];
  localDW->Traj_loc[79].f1[1] = localDW->Derivative2 * rtP.x[79];
  localDW->Traj_loc[80].f1[1] = localDW->Derivative2 * rtP.x[80];
  localDW->Traj_loc[81].f1[1] = localDW->Derivative2 * rtP.x[81];
  localDW->Traj_loc[82].f1[1] = localDW->Derivative2 * rtP.x[82];
  localDW->Traj_loc[83].f1[1] = localDW->Derivative2 * rtP.x[83];
  localDW->Traj_loc[84].f1[1] = localDW->Derivative2 * rtP.x[84];
  localDW->Traj_loc[85].f1[1] = localDW->Derivative2 * rtP.x[85];
  localDW->Traj_loc[86].f1[1] = localDW->Derivative2 * rtP.x[86];
  localDW->Traj_loc[87].f1[1] = localDW->Derivative2 * rtP.x[87];
  localDW->Traj_loc[88].f1[1] = localDW->Derivative2 * rtP.x[88];
  localDW->Traj_loc[89].f1[1] = localDW->Derivative2 * rtP.x[89];
  localDW->Traj_loc[90].f1[1] = localDW->Derivative2 * rtP.x[90];
  localDW->Traj_loc[91].f1[1] = localDW->Derivative2 * rtP.x[91];
  localDW->Traj_loc[92].f1[1] = localDW->Derivative2 * rtP.x[92];
  localDW->Traj_loc[93].f1[1] = localDW->Derivative2 * rtP.x[93];
  localDW->Traj_loc[94].f1[1] = localDW->Derivative2 * rtP.x[94];
  localDW->Traj_loc[95].f1[1] = localDW->Derivative2 * rtP.x[95];
  localDW->Traj_loc[96].f1[1] = localDW->Derivative2 * rtP.x[96];
  localDW->Traj_loc[97].f1[1] = localDW->Derivative2 * rtP.x[97];
  localDW->Traj_loc[98].f1[1] = localDW->Derivative2 * rtP.x[98];
  localDW->Traj_loc[99].f1[1] = localDW->Derivative2 * rtP.x[99];
  localDW->Traj_loc[100].f1[1] = localDW->Derivative2 * rtP.x[100];
  localDW->Traj_loc[101].f1[1] = localDW->Derivative2 * rtP.x[101];
  localDW->Traj_loc[102].f1[1] = localDW->Derivative2 * rtP.x[102];
  localDW->Traj_loc[103].f1[1] = localDW->Derivative2 * rtP.x[103];
  localDW->Traj_loc[104].f1[1] = localDW->Derivative2 * rtP.x[104];
  localDW->Traj_loc[105].f1[1] = localDW->Derivative2 * rtP.x[105];
  localDW->Traj_loc[106].f1[1] = localDW->Derivative2 * rtP.x[106];
  localDW->Traj_loc[107].f1[1] = localDW->Derivative2 * rtP.x[107];
  localDW->Traj_loc[108].f1[1] = localDW->Derivative2 * rtP.x[108];
  localDW->Traj_loc[109].f1[1] = localDW->Derivative2 * rtP.x[109];
  localDW->Traj_loc[110].f1[1] = localDW->Derivative2 * rtP.x[110];
  localDW->Traj_loc[111].f1[1] = localDW->Derivative2 * rtP.x[111];
  localDW->Traj_loc[112].f1[1] = localDW->Derivative2 * rtP.x[112];
  localDW->Traj_loc[113].f1[1] = localDW->Derivative2 * rtP.x[113];
  localDW->Traj_loc[114].f1[1] = localDW->Derivative2 * rtP.x[114];
  localDW->Traj_loc[115].f1[1] = localDW->Derivative2 * rtP.x[115];
  localDW->Traj_loc[116].f1[1] = localDW->Derivative2 * rtP.x[116];
  localDW->Traj_loc[117].f1[1] = localDW->Derivative2 * rtP.x[117];
  localDW->Traj_loc[118].f1[1] = localDW->Derivative2 * rtP.x[118];
  localDW->Traj_loc[119].f1[1] = localDW->Derivative2 * rtP.x[119];
  localDW->Traj_loc[120].f1[1] = localDW->Derivative2 * rtP.x[120];
  localDW->Traj_loc[121].f1[1] = localDW->Derivative2 * rtP.x[121];
  localDW->Traj_loc[122].f1[1] = localDW->Derivative2 * rtP.x[122];
  localDW->Traj_loc[123].f1[1] = localDW->Derivative2 * rtP.x[123];
  localDW->Traj_loc[124].f1[1] = localDW->Derivative2 * rtP.x[124];
  localDW->Traj_loc[125].f1[1] = localDW->Derivative2 * rtP.x[125];
  localDW->Traj_loc[126].f1[1] = localDW->Derivative2 * rtP.x[126];
  localDW->Traj_loc[127].f1[1] = localDW->Derivative2 * rtP.x[127];
  localDW->Traj_loc[128].f1[1] = localDW->Derivative2 * rtP.x[128];
  localDW->Traj_loc[129].f1[1] = localDW->Derivative2 * rtP.x[129];
  localDW->Traj_loc[130].f1[1] = localDW->Derivative2 * rtP.x[130];
  localDW->Traj_loc[131].f1[1] = localDW->Derivative2 * rtP.x[131];
  localDW->Traj_loc[132].f1[1] = localDW->Derivative2 * rtP.x[132];
  localDW->Traj_loc[133].f1[1] = localDW->Derivative2 * rtP.x[133];
  localDW->Traj_loc[134].f1[1] = localDW->Derivative2 * rtP.x[134];
  localDW->Traj_loc[135].f1[1] = localDW->Derivative2 * rtP.x[135];
  localDW->Traj_loc[136].f1[1] = localDW->Derivative2 * rtP.x[136];
  localDW->Traj_loc[137].f1[1] = localDW->Derivative2 * rtP.x[137];
  localDW->Traj_loc[138].f1[1] = localDW->Derivative2 * rtP.x[138];
  localDW->Traj_loc[139].f1[1] = localDW->Derivative2 * rtP.x[139];
  localDW->Traj_loc[140].f1[1] = localDW->Derivative2 * rtP.x[140];
  localDW->Traj_loc[141].f1[1] = localDW->Derivative2 * rtP.x[141];
  localDW->Traj_loc[142].f1[1] = localDW->Derivative2 * rtP.x[142];
  localDW->Traj_loc[143].f1[1] = localDW->Derivative2 * rtP.x[143];
  localDW->Traj_loc[144].f1[1] = localDW->Derivative2 * rtP.x[144];
  localDW->Traj_loc[145].f1[1] = localDW->Derivative2 * rtP.x[145];
  localDW->Traj_loc[146].f1[1] = localDW->Derivative2 * rtP.x[146];
  localDW->Traj_loc[147].f1[1] = localDW->Derivative2 * rtP.x[147];
  localDW->Traj_loc[148].f1[1] = localDW->Derivative2 * rtP.x[148];
  localDW->Traj_loc[149].f1[1] = localDW->Derivative2 * rtP.x[149];
  localDW->Traj_loc[150].f1[1] = localDW->Derivative2 * rtP.x[150];
  localDW->Traj_loc[151].f1[1] = localDW->Derivative2 * rtP.x[151];
  localDW->Traj_loc[152].f1[1] = localDW->Derivative2 * rtP.x[152];
  localDW->Traj_loc[153].f1[1] = localDW->Derivative2 * rtP.x[153];
  localDW->Traj_loc[154].f1[1] = localDW->Derivative2 * rtP.x[154];
  localDW->Traj_loc[155].f1[1] = localDW->Derivative2 * rtP.x[155];
  localDW->Traj_loc[156].f1[1] = localDW->Derivative2 * rtP.x[156];
  localDW->Traj_loc[157].f1[1] = localDW->Derivative2 * rtP.x[157];
  localDW->Traj_loc[158].f1[1] = localDW->Derivative2 * rtP.x[158];
  localDW->Traj_loc[159].f1[1] = localDW->Derivative2 * rtP.x[159];
  localDW->Traj_loc[160].f1[1] = localDW->Derivative2 * rtP.x[160];
  localDW->Traj_loc[161].f1[1] = localDW->Derivative2 * rtP.x[161];
  localDW->Traj_loc[162].f1[1] = localDW->Derivative2 * rtP.x[162];
  localDW->Traj_loc[163].f1[1] = localDW->Derivative2 * rtP.x[163];
  localDW->Traj_loc[164].f1[1] = localDW->Derivative2 * rtP.x[164];
  localDW->Traj_loc[165].f1[1] = localDW->Derivative2 * rtP.x[165];
  localDW->Traj_loc[166].f1[1] = localDW->Derivative2 * rtP.x[166];
  localDW->Traj_loc[167].f1[1] = localDW->Derivative2 * rtP.x[167];
  localDW->Traj_loc[168].f1[1] = localDW->Derivative2 * rtP.x[168];
  localDW->Traj_loc[169].f1[1] = localDW->Derivative2 * rtP.x[169];
  localDW->Traj_loc[170].f1[1] = localDW->Derivative2 * rtP.x[170];
  localDW->Traj_loc[171].f1[1] = localDW->Derivative2 * rtP.x[171];
  localDW->Traj_loc[172].f1[1] = localDW->Derivative2 * rtP.x[172];
  localDW->Traj_loc[173].f1[1] = localDW->Derivative2 * rtP.x[173];
  localDW->Traj_loc[174].f1[1] = localDW->Derivative2 * rtP.x[174];
  localDW->Traj_loc[175].f1[1] = localDW->Derivative2 * rtP.x[175];
  localDW->Traj_loc[176].f1[1] = localDW->Derivative2 * rtP.x[176];
  localDW->Traj_loc[177].f1[1] = localDW->Derivative2 * rtP.x[177];
  localDW->Traj_loc[178].f1[1] = localDW->Derivative2 * rtP.x[178];
  localDW->Traj_loc[179].f1[1] = localDW->Derivative2 * rtP.x[179];
  localDW->Traj_loc[180].f1[1] = localDW->Derivative2 * rtP.x[180];
  localDW->Traj_loc[181].f1[1] = localDW->Derivative2 * rtP.x[181];
  localDW->Traj_loc[182].f1[1] = localDW->Derivative2 * rtP.x[182];
  localDW->Traj_loc[183].f1[1] = localDW->Derivative2 * rtP.x[183];
  localDW->Traj_loc[184].f1[1] = localDW->Derivative2 * rtP.x[184];
  localDW->Traj_loc[185].f1[1] = localDW->Derivative2 * rtP.x[185];
  localDW->Traj_loc[186].f1[1] = localDW->Derivative2 * rtP.x[186];
  localDW->Traj_loc[187].f1[1] = localDW->Derivative2 * rtP.x[187];
  localDW->Traj_loc[188].f1[1] = localDW->Derivative2 * rtP.x[188];
  localDW->Traj_loc[189].f1[1] = localDW->Derivative2 * rtP.x[189];
  localDW->Traj_loc[190].f1[1] = localDW->Derivative2 * rtP.x[190];
  localDW->Traj_loc[191].f1[1] = localDW->Derivative2 * rtP.x[191];
  localDW->Traj_loc[192].f1[1] = localDW->Derivative2 * rtP.x[192];
  localDW->Traj_loc[193].f1[1] = localDW->Derivative2 * rtP.x[193];
  localDW->Traj_loc[194].f1[1] = localDW->Derivative2 * rtP.x[194];
  localDW->Traj_loc[195].f1[1] = localDW->Derivative2 * rtP.x[195];
  localDW->Traj_loc[196].f1[1] = localDW->Derivative2 * rtP.x[196];
  localDW->Traj_loc[197].f1[1] = localDW->Derivative2 * rtP.x[197];
  localDW->Traj_loc[198].f1[1] = localDW->Derivative2 * rtP.x[198];
  localDW->Traj_loc[199].f1[1] = localDW->Derivative2 * rtP.x[199];
  localDW->Traj_loc[200].f1[1] = localDW->Derivative2 * rtP.x[200];
  localDW->Traj_loc[201].f1[1] = localDW->Derivative2 * rtP.x[201];
  localDW->Traj_loc[202].f1[1] = localDW->Derivative2 * rtP.x[202];
  localDW->Traj_loc[203].f1[1] = localDW->Derivative2 * rtP.x[203];
  localDW->Traj_loc[204].f1[1] = localDW->Derivative2 * rtP.x[204];
  localDW->Traj_loc[205].f1[1] = localDW->Derivative2 * rtP.x[205];
  localDW->Traj_loc[206].f1[1] = localDW->Derivative2 * rtP.x[206];
  localDW->Traj_loc[207].f1[1] = localDW->Derivative2 * rtP.x[207];
  localDW->Traj_loc[208].f1[1] = localDW->Derivative2 * rtP.x[208];
  localDW->Traj_loc[209].f1[1] = localDW->Derivative2 * rtP.x[209];
  localDW->Traj_loc[210].f1[1] = localDW->Derivative2 * rtP.x[210];
  localDW->Traj_loc[211].f1[1] = localDW->Derivative2 * rtP.x[211];
  localDW->Traj_loc[212].f1[1] = localDW->Derivative2 * rtP.x[212];
  localDW->Traj_loc[213].f1[1] = localDW->Derivative2 * rtP.x[213];
  localDW->Traj_loc[214].f1[1] = localDW->Derivative2 * rtP.x[214];
  localDW->Traj_loc[215].f1[1] = localDW->Derivative2 * rtP.x[215];
  localDW->Traj_loc[216].f1[1] = localDW->Derivative2 * rtP.x[216];
  localDW->Traj_loc[217].f1[1] = localDW->Derivative2 * rtP.x[217];
  localDW->Traj_loc[218].f1[1] = localDW->Derivative2 * rtP.x[218];
  localDW->Traj_loc[219].f1[1] = localDW->Derivative2 * rtP.x[219];
  localDW->Traj_loc[220].f1[1] = localDW->Derivative2 * rtP.x[220];
  localDW->Traj_loc[221].f1[1] = localDW->Derivative2 * rtP.x[221];
  localDW->Traj_loc[222].f1[1] = localDW->Derivative2 * rtP.x[222];
  localDW->Traj_loc[223].f1[1] = localDW->Derivative2 * rtP.x[223];
  localDW->Traj_loc[224].f1[1] = localDW->Derivative2 * rtP.x[224];
  localDW->Traj_loc[225].f1[1] = localDW->Derivative2 * rtP.x[225];
  localDW->Traj_loc[226].f1[1] = localDW->Derivative2 * rtP.x[226];
  localDW->Traj_loc[227].f1[1] = localDW->Derivative2 * rtP.x[227];
  localDW->Traj_loc[228].f1[1] = localDW->Derivative2 * rtP.x[228];
  localDW->Traj_loc[229].f1[1] = localDW->Derivative2 * rtP.x[229];
  localDW->Traj_loc[230].f1[1] = localDW->Derivative2 * rtP.x[230];
  localDW->Traj_loc[231].f1[1] = localDW->Derivative2 * rtP.x[231];
  localDW->Traj_loc[232].f1[1] = localDW->Derivative2 * rtP.x[232];
  localDW->Traj_loc[233].f1[1] = localDW->Derivative2 * rtP.x[233];
  localDW->Traj_loc[234].f1[1] = localDW->Derivative2 * rtP.x[234];
  localDW->Traj_loc[235].f1[1] = localDW->Derivative2 * rtP.x[235];
  localDW->Traj_loc[236].f1[1] = localDW->Derivative2 * rtP.x[236];
  localDW->Traj_loc[237].f1[1] = localDW->Derivative2 * rtP.x[237];
  localDW->Traj_loc[238].f1[1] = localDW->Derivative2 * rtP.x[238];
  localDW->Traj_loc[239].f1[1] = localDW->Derivative2 * rtP.x[239];
  localDW->Traj_loc[240].f1[1] = localDW->Derivative2 * rtP.x[240];
  localDW->Traj_loc[241].f1[1] = localDW->Derivative2 * rtP.x[241];
  localDW->Traj_loc[242].f1[1] = localDW->Derivative2 * rtP.x[242];
  localDW->Traj_loc[243].f1[1] = localDW->Derivative2 * rtP.x[243];
  localDW->Traj_loc[244].f1[1] = localDW->Derivative2 * rtP.x[244];
  localDW->Traj_loc[245].f1[1] = localDW->Derivative2 * rtP.x[245];
  localDW->Traj_loc[246].f1[1] = localDW->Derivative2 * rtP.x[246];
  localDW->Traj_loc[247].f1[1] = localDW->Derivative2 * rtP.x[247];
  localDW->Traj_loc[248].f1[1] = localDW->Derivative2 * rtP.x[248];
  localDW->Traj_loc[249].f1[1] = localDW->Derivative2 * rtP.x[249];
  localDW->Traj_loc[250].f1[1] = localDW->Derivative2 * rtP.x[250];
  localDW->Traj_loc[251].f1[1] = localDW->Derivative2 * rtP.x[251];
  localDW->Traj_loc[252].f1[1] = localDW->Derivative2 * rtP.x[252];
  localDW->Traj_loc[253].f1[1] = localDW->Derivative2 * rtP.x[253];
  localDW->Traj_loc[254].f1[1] = localDW->Derivative2 * rtP.x[254];
  localDW->Traj_loc[255].f1[1] = localDW->Derivative2 * rtP.x[255];
  localDW->Traj_loc[256].f1[1] = localDW->Derivative2 * rtP.x[256];
  localDW->Traj_loc[257].f1[1] = localDW->Derivative2 * rtP.x[257];
  localDW->Traj_loc[258].f1[1] = localDW->Derivative2 * rtP.x[258];
  localDW->Traj_loc[259].f1[1] = localDW->Derivative2 * rtP.x[259];
  localDW->Traj_loc[260].f1[1] = localDW->Derivative2 * rtP.x[260];
  localDW->Traj_loc[261].f1[1] = localDW->Derivative2 * rtP.x[261];
  localDW->Traj_loc[262].f1[1] = localDW->Derivative2 * rtP.x[262];
  localDW->Traj_loc[263].f1[1] = localDW->Derivative2 * rtP.x[263];
  localDW->Traj_loc[264].f1[1] = localDW->Derivative2 * rtP.x[264];
  localDW->Traj_loc[265].f1[1] = localDW->Derivative2 * rtP.x[265];
  localDW->Traj_loc[266].f1[1] = localDW->Derivative2 * rtP.x[266];
  localDW->Traj_loc[267].f1[1] = localDW->Derivative2 * rtP.x[267];
  localDW->Traj_loc[268].f1[1] = localDW->Derivative2 * rtP.x[268];
  localDW->Traj_loc[269].f1[1] = localDW->Derivative2 * rtP.x[269];
  localDW->Traj_loc[270].f1[1] = localDW->Derivative2 * rtP.x[270];
  localDW->Traj_loc[271].f1[1] = localDW->Derivative2 * rtP.x[271];
  localDW->Traj_loc[272].f1[1] = localDW->Derivative2 * rtP.x[272];
  localDW->Traj_loc[273].f1[1] = localDW->Derivative2 * rtP.x[273];
  localDW->Traj_loc[274].f1[1] = localDW->Derivative2 * rtP.x[274];
  localDW->Traj_loc[275].f1[1] = localDW->Derivative2 * rtP.x[275];
  localDW->Traj_loc[276].f1[1] = localDW->Derivative2 * rtP.x[276];
  localDW->Traj_loc[277].f1[1] = localDW->Derivative2 * rtP.x[277];
  localDW->Traj_loc[278].f1[1] = localDW->Derivative2 * rtP.x[278];
  localDW->Traj_loc[279].f1[1] = localDW->Derivative2 * rtP.x[279];
  localDW->Traj_loc[280].f1[1] = localDW->Derivative2 * rtP.x[280];
  localDW->Traj_loc[281].f1[1] = localDW->Derivative2 * rtP.x[281];
  localDW->Traj_loc[282].f1[1] = localDW->Derivative2 * rtP.x[282];
  localDW->Traj_loc[283].f1[1] = localDW->Derivative2 * rtP.x[283];
  localDW->Traj_loc[284].f1[1] = localDW->Derivative2 * rtP.x[284];
  localDW->Traj_loc[285].f1[1] = localDW->Derivative2 * rtP.x[285];
  localDW->Traj_loc[286].f1[1] = localDW->Derivative2 * rtP.x[286];
  localDW->Traj_loc[287].f1[1] = localDW->Derivative2 * rtP.x[287];
  localDW->Traj_loc[288].f1[1] = localDW->Derivative2 * rtP.x[288];
  localDW->Traj_loc[289].f1[1] = localDW->Derivative2 * rtP.x[289];
  localDW->Traj_loc[290].f1[1] = localDW->Derivative2 * rtP.x[290];
  localDW->Traj_loc[291].f1[1] = localDW->Derivative2 * rtP.x[291];
  localDW->Traj_loc[292].f1[1] = localDW->Derivative2 * rtP.x[292];
  localDW->Traj_loc[293].f1[1] = localDW->Derivative2 * rtP.x[293];
  localDW->Traj_loc[294].f1[1] = localDW->Derivative2 * rtP.x[294];
  localDW->Traj_loc[295].f1[1] = localDW->Derivative2 * rtP.x[295];
  localDW->Traj_loc[296].f1[1] = localDW->Derivative2 * rtP.x[296];
  localDW->Traj_loc[297].f1[1] = localDW->Derivative2 * rtP.x[297];
  localDW->Traj_loc[298].f1[1] = localDW->Derivative2 * rtP.x[298];
  localDW->Traj_loc[299].f1[1] = localDW->Derivative2 * rtP.x[299];
  localDW->Traj_loc[300].f1[1] = localDW->Derivative2 * rtP.x[300];
  localDW->Traj_loc[301].f1[1] = localDW->Derivative2 * rtP.x[301];
  localDW->Traj_loc[302].f1[1] = localDW->Derivative2 * rtP.x[302];
  localDW->Traj_loc[303].f1[1] = localDW->Derivative2 * rtP.x[303];
  localDW->Traj_loc[304].f1[1] = localDW->Derivative2 * rtP.x[304];
  localDW->Traj_loc[305].f1[1] = localDW->Derivative2 * rtP.x[305];
  localDW->Traj_loc[306].f1[1] = localDW->Derivative2 * rtP.x[306];
  localDW->Traj_loc[307].f1[1] = localDW->Derivative2 * rtP.x[307];
  localDW->Traj_loc[308].f1[1] = localDW->Derivative2 * rtP.x[308];
  localDW->Traj_loc[309].f1[1] = localDW->Derivative2 * rtP.x[309];
  localDW->Traj_loc[310].f1[1] = localDW->Derivative2 * rtP.x[310];
  localDW->Traj_loc[311].f1[1] = localDW->Derivative2 * rtP.x[311];
  localDW->Traj_loc[312].f1[1] = localDW->Derivative2 * rtP.x[312];
  localDW->Traj_loc[313].f1[1] = localDW->Derivative2 * rtP.x[313];
  localDW->Traj_loc[314].f1[1] = localDW->Derivative2 * rtP.x[314];
  localDW->Traj_loc[315].f1[1] = localDW->Derivative2 * rtP.x[315];
  localDW->Traj_loc[316].f1[1] = localDW->Derivative2 * rtP.x[316];
  localDW->Traj_loc[317].f1[1] = localDW->Derivative2 * rtP.x[317];
  localDW->Traj_loc[318].f1[1] = localDW->Derivative2 * rtP.x[318];
  localDW->Traj_loc[319].f1[1] = localDW->Derivative2 * rtP.x[319];
  localDW->Traj_loc[320].f1[1] = localDW->Derivative2 * rtP.x[320];
  localDW->Traj_loc[321].f1[1] = localDW->Derivative2 * rtP.x[321];
  localDW->Traj_loc[322].f1[1] = localDW->Derivative2 * rtP.x[322];
  localDW->Traj_loc[323].f1[1] = localDW->Derivative2 * rtP.x[323];
  localDW->Traj_loc[324].f1[1] = localDW->Derivative2 * rtP.x[324];
  localDW->Traj_loc[325].f1[1] = localDW->Derivative2 * rtP.x[325];
  localDW->Traj_loc[326].f1[1] = localDW->Derivative2 * rtP.x[326];
  localDW->Traj_loc[327].f1[1] = localDW->Derivative2 * rtP.x[327];
  localDW->Traj_loc[328].f1[1] = localDW->Derivative2 * rtP.x[328];
  localDW->Traj_loc[329].f1[1] = localDW->Derivative2 * rtP.x[329];
  localDW->Traj_loc[330].f1[1] = localDW->Derivative2 * rtP.x[330];
  localDW->Traj_loc[331].f1[1] = localDW->Derivative2 * rtP.x[331];
  localDW->Traj_loc[332].f1[1] = localDW->Derivative2 * rtP.x[332];
  localDW->Traj_loc[333].f1[1] = localDW->Derivative2 * rtP.x[333];
  localDW->Traj_loc[334].f1[1] = localDW->Derivative2 * rtP.x[334];
  localDW->Traj_loc[335].f1[1] = localDW->Derivative2 * rtP.x[335];
  localDW->Traj_loc[336].f1[1] = localDW->Derivative2 * rtP.x[336];
  localDW->Traj_loc[337].f1[1] = localDW->Derivative2 * rtP.x[337];
  localDW->Traj_loc[338].f1[1] = localDW->Derivative2 * rtP.x[338];
  localDW->Traj_loc[339].f1[1] = localDW->Derivative2 * rtP.x[339];
  localDW->Traj_loc[340].f1[1] = localDW->Derivative2 * rtP.x[340];
  localDW->Traj_loc[341].f1[1] = localDW->Derivative2 * rtP.x[341];
  localDW->Traj_loc[342].f1[1] = localDW->Derivative2 * rtP.x[342];
  localDW->Traj_loc[343].f1[1] = localDW->Derivative2 * rtP.x[343];
  localDW->Traj_loc[344].f1[1] = localDW->Derivative2 * rtP.x[344];
  localDW->Traj_loc[345].f1[1] = localDW->Derivative2 * rtP.x[345];
  localDW->Traj_loc[346].f1[1] = localDW->Derivative2 * rtP.x[346];
  localDW->Traj_loc[347].f1[1] = localDW->Derivative2 * rtP.x[347];
  localDW->Traj_loc[348].f1[1] = localDW->Derivative2 * rtP.x[348];
  localDW->Traj_loc[349].f1[1] = localDW->Derivative2 * rtP.x[349];
  localDW->Traj_loc[350].f1[1] = localDW->Derivative2 * rtP.x[350];
  localDW->Traj_loc[351].f1[1] = localDW->Derivative2 * rtP.x[351];
  localDW->Traj_loc[352].f1[1] = localDW->Derivative2 * rtP.x[352];
  localDW->Traj_loc[353].f1[1] = localDW->Derivative2 * rtP.x[353];
  localDW->Traj_loc[354].f1[1] = localDW->Derivative2 * rtP.x[354];
  localDW->Traj_loc[355].f1[1] = localDW->Derivative2 * rtP.x[355];
  localDW->Traj_loc[356].f1[1] = localDW->Derivative2 * rtP.x[356];
  localDW->Traj_loc[357].f1[1] = localDW->Derivative2 * rtP.x[357];
  localDW->Traj_loc[358].f1[1] = localDW->Derivative2 * rtP.x[358];
  localDW->Traj_loc[359].f1[1] = localDW->Derivative2 * rtP.x[359];
  localDW->Traj_loc[360].f1[1] = localDW->Derivative2 * rtP.x[360];
  localDW->Traj_loc[361].f1[1] = localDW->Derivative2 * rtP.x[361];
  localDW->Traj_loc[362].f1[1] = localDW->Derivative2 * rtP.x[362];
  localDW->Traj_loc[363].f1[1] = localDW->Derivative2 * rtP.x[363];
  localDW->Traj_loc[364].f1[1] = localDW->Derivative2 * rtP.x[364];
  localDW->Traj_loc[365].f1[1] = localDW->Derivative2 * rtP.x[365];
  localDW->Traj_loc[366].f1[1] = localDW->Derivative2 * rtP.x[366];
  localDW->Traj_loc[367].f1[1] = localDW->Derivative2 * rtP.x[367];
  localDW->Traj_loc[368].f1[1] = localDW->Derivative2 * rtP.x[368];
  localDW->Traj_loc[369].f1[1] = localDW->Derivative2 * rtP.x[369];
  localDW->Traj_loc[370].f1[1] = localDW->Derivative2 * rtP.x[370];
  localDW->Traj_loc[371].f1[1] = localDW->Derivative2 * rtP.x[371];
  localDW->Traj_loc[372].f1[1] = localDW->Derivative2 * rtP.x[372];
  localDW->Traj_loc[373].f1[1] = localDW->Derivative2 * rtP.x[373];
  localDW->Traj_loc[374].f1[1] = localDW->Derivative2 * rtP.x[374];
  localDW->Traj_loc[375].f1[1] = localDW->Derivative2 * rtP.x[375];
  localDW->Traj_loc[376].f1[1] = localDW->Derivative2 * rtP.x[376];
  localDW->Traj_loc[377].f1[1] = localDW->Derivative2 * rtP.x[377];
  localDW->Traj_loc[378].f1[1] = localDW->Derivative2 * rtP.x[378];
  localDW->Traj_loc[379].f1[1] = localDW->Derivative2 * rtP.x[379];
  localDW->Traj_loc[380].f1[1] = localDW->Derivative2 * rtP.x[380];
  localDW->Traj_loc[381].f1[1] = localDW->Derivative2 * rtP.x[381];
  localDW->Traj_loc[382].f1[1] = localDW->Derivative2 * rtP.x[382];
  localDW->Traj_loc[383].f1[1] = localDW->Derivative2 * rtP.x[383];
  localDW->Traj_loc[384].f1[1] = localDW->Derivative2 * rtP.x[384];
  localDW->Traj_loc[385].f1[1] = localDW->Derivative2 * rtP.x[385];
  localDW->Traj_loc[386].f1[1] = localDW->Derivative2 * rtP.x[386];
  localDW->Traj_loc[387].f1[1] = localDW->Derivative2 * rtP.x[387];
  localDW->Traj_loc[388].f1[1] = localDW->Derivative2 * rtP.x[388];
  localDW->Traj_loc[389].f1[1] = localDW->Derivative2 * rtP.x[389];
  localDW->Traj_loc[390].f1[1] = localDW->Derivative2 * rtP.x[390];
  localDW->Traj_loc[391].f1[1] = localDW->Derivative2 * rtP.x[391];
  localDW->Traj_loc[392].f1[1] = localDW->Derivative2 * rtP.x[392];
  localDW->Traj_loc[393].f1[1] = localDW->Derivative2 * rtP.x[393];
  localDW->Traj_loc[394].f1[1] = localDW->Derivative2 * rtP.x[394];
  localDW->Traj_loc[395].f1[1] = localDW->Derivative2 * rtP.x[395];
  localDW->Traj_loc[396].f1[1] = localDW->Derivative2 * rtP.x[396];
  localDW->Traj_loc[397].f1[1] = localDW->Derivative2 * rtP.x[397];
  localDW->Traj_loc[398].f1[1] = localDW->Derivative2 * rtP.x[398];
  localDW->Traj_loc[399].f1[1] = localDW->Derivative2 * rtP.x[399];
  localDW->Traj_loc[400].f1[1] = localDW->Derivative2 * rtP.x[400];
  localDW->Traj_loc[401].f1[1] = localDW->Derivative2 * rtP.x[401];
  localDW->Traj_loc[402].f1[1] = localDW->Derivative2 * rtP.x[402];
  localDW->Traj_loc[403].f1[1] = localDW->Derivative2 * rtP.x[403];
  localDW->Traj_loc[404].f1[1] = localDW->Derivative2 * rtP.x[404];
  localDW->Traj_loc[405].f1[1] = localDW->Derivative2 * rtP.x[405];
  localDW->Traj_loc[406].f1[1] = localDW->Derivative2 * rtP.x[406];
  localDW->Traj_loc[407].f1[1] = localDW->Derivative2 * rtP.x[407];
  localDW->Traj_loc[408].f1[1] = localDW->Derivative2 * rtP.x[408];
  localDW->Traj_loc[409].f1[1] = localDW->Derivative2 * rtP.x[409];
  localDW->Traj_loc[410].f1[1] = localDW->Derivative2 * rtP.x[410];
  localDW->Traj_loc[411].f1[1] = localDW->Derivative2 * rtP.x[411];
  localDW->Traj_loc[412].f1[1] = localDW->Derivative2 * rtP.x[412];
  localDW->Traj_loc[413].f1[1] = localDW->Derivative2 * rtP.x[413];
  localDW->Traj_loc[414].f1[1] = localDW->Derivative2 * rtP.x[414];
  localDW->Traj_loc[415].f1[1] = localDW->Derivative2 * rtP.x[415];
  localDW->Traj_loc[416].f1[1] = localDW->Derivative2 * rtP.x[416];
  localDW->Traj_loc[417].f1[1] = localDW->Derivative2 * rtP.x[417];
  localDW->Traj_loc[418].f1[1] = localDW->Derivative2 * rtP.x[418];
  localDW->Traj_loc[419].f1[1] = localDW->Derivative2 * rtP.x[419];
  localDW->Traj_loc[420].f1[1] = localDW->Derivative2 * rtP.x[420];
  localDW->Traj_loc[421].f1[1] = localDW->Derivative2 * rtP.x[421];
  localDW->Traj_loc[422].f1[1] = localDW->Derivative2 * rtP.x[422];
  localDW->Traj_loc[423].f1[1] = localDW->Derivative2 * rtP.x[423];
  localDW->Traj_loc[424].f1[1] = localDW->Derivative2 * rtP.x[424];
  localDW->Traj_loc[425].f1[1] = localDW->Derivative2 * rtP.x[425];
  localDW->Traj_loc[426].f1[1] = localDW->Derivative2 * rtP.x[426];
  localDW->Traj_loc[427].f1[1] = localDW->Derivative2 * rtP.x[427];
  localDW->Traj_loc[428].f1[1] = localDW->Derivative2 * rtP.x[428];
  localDW->Traj_loc[429].f1[1] = localDW->Derivative2 * rtP.x[429];
  localDW->Traj_loc[430].f1[1] = localDW->Derivative2 * rtP.x[430];
  localDW->Traj_loc[431].f1[1] = localDW->Derivative2 * rtP.x[431];
  localDW->Traj_loc[432].f1[1] = localDW->Derivative2 * rtP.x[432];
  localDW->Traj_loc[433].f1[1] = localDW->Derivative2 * rtP.x[433];
  localDW->Traj_loc[434].f1[1] = localDW->Derivative2 * rtP.x[434];
  localDW->Traj_loc[435].f1[1] = localDW->Derivative2 * rtP.x[435];
  localDW->Traj_loc[436].f1[1] = localDW->Derivative2 * rtP.x[436];
  localDW->Traj_loc[437].f1[1] = localDW->Derivative2 * rtP.x[437];
  localDW->Traj_loc[438].f1[1] = localDW->Derivative2 * rtP.x[438];
  localDW->Traj_loc[439].f1[1] = localDW->Derivative2 * rtP.x[439];
  localDW->Traj_loc[440].f1[1] = localDW->Derivative2 * rtP.x[440];
  localDW->Traj_loc[441].f1[1] = localDW->Derivative2 * rtP.x[441];
  localDW->Traj_loc[442].f1[1] = localDW->Derivative2 * rtP.x[442];
  localDW->Traj_loc[443].f1[1] = localDW->Derivative2 * rtP.x[443];
  localDW->Traj_loc[444].f1[1] = localDW->Derivative2 * rtP.x[444];
  localDW->Traj_loc[445].f1[1] = localDW->Derivative2 * rtP.x[445];
  localDW->Traj_loc[446].f1[1] = localDW->Derivative2 * rtP.x[446];
  localDW->Traj_loc[447].f1[1] = localDW->Derivative2 * rtP.x[447];
  localDW->Traj_loc[448].f1[1] = localDW->Derivative2 * rtP.x[448];
  localDW->Traj_loc[449].f1[1] = localDW->Derivative2 * rtP.x[449];
  localDW->Traj_loc[450].f1[1] = localDW->Derivative2 * rtP.x[450];
  localDW->Traj_loc[451].f1[1] = localDW->Derivative2 * rtP.x[451];
  localDW->Traj_loc[452].f1[1] = localDW->Derivative2 * rtP.x[452];
  localDW->Traj_loc[453].f1[1] = localDW->Derivative2 * rtP.x[453];
  localDW->Traj_loc[454].f1[1] = localDW->Derivative2 * rtP.x[454];
  localDW->Traj_loc[455].f1[1] = localDW->Derivative2 * rtP.x[455];
  localDW->Traj_loc[456].f1[1] = localDW->Derivative2 * rtP.x[456];
  localDW->Traj_loc[457].f1[1] = localDW->Derivative2 * rtP.x[457];
  localDW->Traj_loc[458].f1[1] = localDW->Derivative2 * rtP.x[458];
  localDW->Traj_loc[459].f1[1] = localDW->Derivative2 * rtP.x[459];
  localDW->Traj_loc[460].f1[1] = localDW->Derivative2 * rtP.x[460];
  localDW->Traj_loc[461].f1[1] = localDW->Derivative2 * rtP.x[461];
  localDW->Traj_loc[462].f1[1] = localDW->Derivative2 * rtP.x[462];
  localDW->Traj_loc[463].f1[1] = localDW->Derivative2 * rtP.x[463];
  localDW->Traj_loc[464].f1[1] = localDW->Derivative2 * rtP.x[464];
  localDW->Traj_loc[465].f1[1] = localDW->Derivative2 * rtP.x[465];
  localDW->Traj_loc[466].f1[1] = localDW->Derivative2 * rtP.x[466];
  localDW->Traj_loc[467].f1[1] = localDW->Derivative2 * rtP.x[467];
  localDW->Traj_loc[468].f1[1] = localDW->Derivative2 * rtP.x[468];
  localDW->Traj_loc[469].f1[1] = localDW->Derivative2 * rtP.x[469];
  localDW->Traj_loc[470].f1[1] = localDW->Derivative2 * rtP.x[470];
  localDW->Traj_loc[471].f1[1] = localDW->Derivative2 * rtP.x[471];
  localDW->Traj_loc[472].f1[1] = localDW->Derivative2 * rtP.x[472];
  localDW->Traj_loc[473].f1[1] = localDW->Derivative2 * rtP.x[473];
  localDW->Traj_loc[474].f1[1] = localDW->Derivative2 * rtP.x[474];
  localDW->Traj_loc[475].f1[1] = localDW->Derivative2 * rtP.x[475];
  localDW->Traj_loc[476].f1[1] = localDW->Derivative2 * rtP.x[476];
  localDW->Traj_loc[477].f1[1] = localDW->Derivative2 * rtP.x[477];
  localDW->Traj_loc[478].f1[1] = localDW->Derivative2 * rtP.x[478];
  localDW->Traj_loc[479].f1[1] = localDW->Derivative2 * rtP.x[479];
  localDW->Traj_loc[480].f1[1] = localDW->Derivative2 * rtP.x[480];
  localDW->Traj_loc[481].f1[1] = localDW->Derivative2 * rtP.x[481];
  localDW->Traj_loc[482].f1[1] = localDW->Derivative2 * rtP.x[482];
  localDW->Traj_loc[483].f1[1] = localDW->Derivative2 * rtP.x[483];
  localDW->Traj_loc[484].f1[1] = localDW->Derivative2 * rtP.x[484];
  localDW->Traj_loc[485].f1[1] = localDW->Derivative2 * rtP.x[485];
  localDW->Traj_loc[486].f1[1] = localDW->Derivative2 * rtP.x[486];
  localDW->Traj_loc[487].f1[1] = localDW->Derivative2 * rtP.x[487];
  localDW->Traj_loc[488].f1[1] = localDW->Derivative2 * rtP.x[488];
  localDW->Traj_loc[489].f1[1] = localDW->Derivative2 * rtP.x[489];
  localDW->Traj_loc[490].f1[1] = localDW->Derivative2 * rtP.x[490];
  localDW->Traj_loc[491].f1[1] = localDW->Derivative2 * rtP.x[491];
  localDW->Traj_loc[492].f1[1] = localDW->Derivative2 * rtP.x[492];
  localDW->Traj_loc[493].f1[1] = localDW->Derivative2 * rtP.x[493];
  localDW->Traj_loc[494].f1[1] = localDW->Derivative2 * rtP.x[494];
  localDW->Traj_loc[495].f1[1] = localDW->Derivative2 * rtP.x[495];
  localDW->Traj_loc[496].f1[1] = localDW->Derivative2 * rtP.x[496];
  localDW->Traj_loc[497].f1[1] = localDW->Derivative2 * rtP.x[497];
  localDW->Traj_loc[498].f1[1] = localDW->Derivative2 * rtP.x[498];
  localDW->Traj_loc[499].f1[1] = localDW->Derivative2 * rtP.x[499];
  localDW->Traj_loc[500].f1[1] = localDW->Derivative2 * rtP.x[500];
  localDW->Traj_loc[501].f1[1] = localDW->Derivative2 * rtP.x[501];
  localDW->Traj_loc[502].f1[1] = localDW->Derivative2 * rtP.x[502];
  localDW->Traj_loc[503].f1[1] = localDW->Derivative2 * rtP.x[503];
  localDW->Traj_loc[504].f1[1] = localDW->Derivative2 * rtP.x[504];
  localDW->Traj_loc[505].f1[1] = localDW->Derivative2 * rtP.x[505];
  localDW->Traj_loc[506].f1[1] = localDW->Derivative2 * rtP.x[506];
  localDW->Traj_loc[507].f1[1] = localDW->Derivative2 * rtP.x[507];
  localDW->Traj_loc[508].f1[1] = localDW->Derivative2 * rtP.x[508];
  localDW->Traj_loc[509].f1[1] = localDW->Derivative2 * rtP.x[509];
  localDW->Traj_loc[510].f1[1] = localDW->Derivative2 * rtP.x[510];
  localDW->Traj_loc[511].f1[1] = localDW->Derivative2 * rtP.x[511];
  localDW->Traj_loc[512].f1[1] = localDW->Derivative2 * rtP.x[512];
  localDW->Traj_loc[513].f1[1] = localDW->Derivative2 * rtP.x[513];
  localDW->Traj_loc[514].f1[1] = localDW->Derivative2 * rtP.x[514];
  localDW->Traj_loc[515].f1[1] = localDW->Derivative2 * rtP.x[515];
  localDW->Traj_loc[516].f1[1] = localDW->Derivative2 * rtP.x[516];
  localDW->Traj_loc[517].f1[1] = localDW->Derivative2 * rtP.x[517];
  localDW->Traj_loc[518].f1[1] = localDW->Derivative2 * rtP.x[518];
  localDW->Traj_loc[519].f1[1] = localDW->Derivative2 * rtP.x[519];
  localDW->Traj_loc[520].f1[1] = localDW->Derivative2 * rtP.x[520];
  localDW->Traj_loc[521].f1[1] = localDW->Derivative2 * rtP.x[521];
  localDW->Traj_loc[522].f1[1] = localDW->Derivative2 * rtP.x[522];
  localDW->Traj_loc[523].f1[1] = localDW->Derivative2 * rtP.x[523];
  localDW->Traj_loc[524].f1[1] = localDW->Derivative2 * rtP.x[524];
  localDW->Traj_loc[525].f1[1] = localDW->Derivative2 * rtP.x[525];
  localDW->Traj_loc[526].f1[1] = localDW->Derivative2 * rtP.x[526];
  localDW->Traj_loc[527].f1[1] = localDW->Derivative2 * rtP.x[527];
  localDW->Traj_loc[528].f1[1] = localDW->Derivative2 * rtP.x[528];
  localDW->Traj_loc[529].f1[1] = localDW->Derivative2 * rtP.x[529];
  localDW->Traj_loc[530].f1[1] = localDW->Derivative2 * rtP.x[530];
  localDW->Traj_loc[531].f1[1] = localDW->Derivative2 * rtP.x[531];
  localDW->Traj_loc[532].f1[1] = localDW->Derivative2 * rtP.x[532];
  localDW->Traj_loc[533].f1[1] = localDW->Derivative2 * rtP.x[533];
  localDW->Traj_loc[534].f1[1] = localDW->Derivative2 * rtP.x[534];
  localDW->Traj_loc[535].f1[1] = localDW->Derivative2 * rtP.x[535];
  localDW->Traj_loc[536].f1[1] = localDW->Derivative2 * rtP.x[536];
  localDW->Traj_loc[537].f1[1] = localDW->Derivative2 * rtP.x[537];
  localDW->Traj_loc[538].f1[1] = localDW->Derivative2 * rtP.x[538];
  localDW->Traj_loc[539].f1[1] = localDW->Derivative2 * rtP.x[539];
  localDW->Traj_loc[540].f1[1] = localDW->Derivative2 * rtP.x[540];
  localDW->Traj_loc[541].f1[1] = localDW->Derivative2 * rtP.x[541];
  localDW->Traj_loc[542].f1[1] = localDW->Derivative2 * rtP.x[542];
  localDW->Traj_loc[543].f1[1] = localDW->Derivative2 * rtP.x[543];
  localDW->Traj_loc[544].f1[1] = localDW->Derivative2 * rtP.x[544];
  localDW->Traj_loc[545].f1[1] = localDW->Derivative2 * rtP.x[545];
  localDW->Traj_loc[546].f1[1] = localDW->Derivative2 * rtP.x[546];
  localDW->Traj_loc[547].f1[1] = localDW->Derivative2 * rtP.x[547];
  localDW->Traj_loc[548].f1[1] = localDW->Derivative2 * rtP.x[548];
  localDW->Traj_loc[549].f1[1] = localDW->Derivative2 * rtP.x[549];
  localDW->Traj_loc[550].f1[1] = localDW->Derivative2 * rtP.x[550];
  localDW->Traj_loc[551].f1[1] = localDW->Derivative2 * rtP.x[551];
  localDW->Traj_loc[552].f1[1] = localDW->Derivative2 * rtP.x[552];
  localDW->Traj_loc[553].f1[1] = localDW->Derivative2 * rtP.x[553];
  localDW->Traj_loc[554].f1[1] = localDW->Derivative2 * rtP.x[554];
  localDW->Traj_loc[555].f1[1] = localDW->Derivative2 * rtP.x[555];
  localDW->Traj_loc[556].f1[1] = localDW->Derivative2 * rtP.x[556];
  localDW->Traj_loc[557].f1[1] = localDW->Derivative2 * rtP.x[557];
  localDW->Traj_loc[558].f1[1] = localDW->Derivative2 * rtP.x[558];
  localDW->Traj_loc[559].f1[1] = localDW->Derivative2 * rtP.x[559];
  localDW->Traj_loc[560].f1[1] = localDW->Derivative2 * rtP.x[560];
  localDW->Traj_loc[561].f1[1] = localDW->Derivative2 * rtP.x[561];
  localDW->Traj_loc[562].f1[1] = localDW->Derivative2 * rtP.x[562];
  localDW->Traj_loc[563].f1[1] = localDW->Derivative2 * rtP.x[563];
  localDW->Traj_loc[564].f1[1] = localDW->Derivative2 * rtP.x[564];
  localDW->Traj_loc[565].f1[1] = localDW->Derivative2 * rtP.x[565];
  localDW->Traj_loc[566].f1[1] = localDW->Derivative2 * rtP.x[566];
  localDW->Traj_loc[567].f1[1] = localDW->Derivative2 * rtP.x[567];
  localDW->Traj_loc[568].f1[1] = localDW->Derivative2 * rtP.x[568];
  localDW->Traj_loc[569].f1[1] = localDW->Derivative2 * rtP.x[569];
  localDW->Traj_loc[570].f1[1] = localDW->Derivative2 * rtP.x[570];
  localDW->Traj_loc[571].f1[1] = localDW->Derivative2 * rtP.x[571];
  localDW->Traj_loc[572].f1[1] = localDW->Derivative2 * rtP.x[572];
  localDW->Traj_loc[573].f1[1] = localDW->Derivative2 * rtP.x[573];
  localDW->Traj_loc[574].f1[1] = localDW->Derivative2 * rtP.x[574];
  localDW->Traj_loc[575].f1[1] = localDW->Derivative2 * rtP.x[575];
  localDW->Traj_loc[576].f1[1] = localDW->Derivative2 * rtP.x[576];
  localDW->Traj_loc[577].f1[1] = localDW->Derivative2 * rtP.x[577];
  localDW->Traj_loc[578].f1[1] = localDW->Derivative2 * rtP.x[578];
  localDW->Traj_loc[579].f1[1] = localDW->Derivative2 * rtP.x[579];
  localDW->Traj_loc[580].f1[1] = localDW->Derivative2 * rtP.x[580];
  localDW->Traj_loc[581].f1[1] = localDW->Derivative2 * rtP.x[581];
  localDW->Traj_loc[582].f1[1] = localDW->Derivative2 * rtP.x[582];
  localDW->Traj_loc[583].f1[1] = localDW->Derivative2 * rtP.x[583];
  localDW->Traj_loc[584].f1[1] = localDW->Derivative2 * rtP.x[584];
  localDW->Traj_loc[585].f1[1] = localDW->Derivative2 * rtP.x[585];
  localDW->Traj_loc[586].f1[1] = localDW->Derivative2 * rtP.x[586];
  localDW->Traj_loc[587].f1[1] = localDW->Derivative2 * rtP.x[587];
  localDW->Traj_loc[588].f1[1] = localDW->Derivative2 * rtP.x[588];
  localDW->Traj_loc[589].f1[1] = localDW->Derivative2 * rtP.x[589];
  localDW->Traj_loc[590].f1[1] = localDW->Derivative2 * rtP.x[590];
  localDW->Traj_loc[591].f1[1] = localDW->Derivative2 * rtP.x[591];
  localDW->Traj_loc[592].f1[1] = localDW->Derivative2 * rtP.x[592];
  localDW->Traj_loc[593].f1[1] = localDW->Derivative2 * rtP.x[593];
  localDW->Traj_loc[594].f1[1] = localDW->Derivative2 * rtP.x[594];
  localDW->Traj_loc[595].f1[1] = localDW->Derivative2 * rtP.x[595];
  localDW->Traj_loc[596].f1[1] = localDW->Derivative2 * rtP.x[596];
  localDW->Traj_loc[597].f1[1] = localDW->Derivative2 * rtP.x[597];
  localDW->Traj_loc[598].f1[1] = localDW->Derivative2 * rtP.x[598];
  localDW->Traj_loc[599].f1[1] = localDW->Derivative2 * rtP.x[599];
  localDW->Traj_loc[600].f1[1] = localDW->Derivative2 * rtP.x[600];
  localDW->Traj_loc[601].f1[1] = localDW->Derivative2 * rtP.x[601];
  localDW->Traj_loc[602].f1[1] = localDW->Derivative2 * rtP.x[602];
  localDW->Traj_loc[603].f1[1] = localDW->Derivative2 * rtP.x[603];
  localDW->Traj_loc[604].f1[1] = localDW->Derivative2 * rtP.x[604];
  localDW->Traj_loc[605].f1[1] = localDW->Derivative2 * rtP.x[605];
  localDW->Traj_loc[606].f1[1] = localDW->Derivative2 * rtP.x[606];
  localDW->Traj_loc[607].f1[1] = localDW->Derivative2 * rtP.x[607];
  localDW->Traj_loc[608].f1[1] = localDW->Derivative2 * rtP.x[608];
  localDW->Traj_loc[609].f1[1] = localDW->Derivative2 * rtP.x[609];
  localDW->Traj_loc[610].f1[1] = localDW->Derivative2 * rtP.x[610];
  localDW->Traj_loc[611].f1[1] = localDW->Derivative2 * rtP.x[611];
  localDW->Traj_loc[612].f1[1] = localDW->Derivative2 * rtP.x[612];
  localDW->Traj_loc[613].f1[1] = localDW->Derivative2 * rtP.x[613];
  localDW->Traj_loc[614].f1[1] = localDW->Derivative2 * rtP.x[614];
  localDW->Traj_loc[615].f1[1] = localDW->Derivative2 * rtP.x[615];
  localDW->Traj_loc[616].f1[1] = localDW->Derivative2 * rtP.x[616];
  localDW->Traj_loc[617].f1[1] = localDW->Derivative2 * rtP.x[617];
  localDW->Traj_loc[618].f1[1] = localDW->Derivative2 * rtP.x[618];
  localDW->Traj_loc[619].f1[1] = localDW->Derivative2 * rtP.x[619];
  localDW->Traj_loc[620].f1[1] = localDW->Derivative2 * rtP.x[620];
  localDW->Traj_loc[621].f1[1] = localDW->Derivative2 * rtP.x[621];
  localDW->Traj_loc[622].f1[1] = localDW->Derivative2 * rtP.x[622];
  localDW->Traj_loc[623].f1[1] = localDW->Derivative2 * rtP.x[623];
  localDW->Traj_loc[624].f1[1] = localDW->Derivative2 * rtP.x[624];
  localDW->Traj_loc[625].f1[1] = localDW->Derivative2 * rtP.x[625];
  localDW->Traj_loc[626].f1[1] = localDW->Derivative2 * rtP.x[626];
  localDW->Traj_loc[627].f1[1] = localDW->Derivative2 * rtP.x[627];
  localDW->Traj_loc[628].f1[1] = localDW->Derivative2 * rtP.x[628];
  localDW->Traj_loc[629].f1[1] = localDW->Derivative2 * rtP.x[629];
  localDW->Traj_loc[630].f1[1] = localDW->Derivative2 * rtP.x[630];
  localDW->Traj_loc[631].f1[1] = localDW->Derivative2 * rtP.x[631];
  localDW->Traj_loc[632].f1[1] = localDW->Derivative2 * rtP.x[632];
  localDW->Traj_loc[633].f1[1] = localDW->Derivative2 * rtP.x[633];
  localDW->Traj_loc[634].f1[1] = localDW->Derivative2 * rtP.x[634];
  localDW->Traj_loc[635].f1[1] = localDW->Derivative2 * rtP.x[635];
  localDW->Traj_loc[636].f1[1] = localDW->Derivative2 * rtP.x[636];
  localDW->Traj_loc[637].f1[1] = localDW->Derivative2 * rtP.x[637];
  localDW->Traj_loc[638].f1[1] = localDW->Derivative2 * rtP.x[638];
  localDW->Traj_loc[639].f1[1] = localDW->Derivative2 * rtP.x[639];
  localDW->Traj_loc[640].f1[1] = localDW->Derivative2 * rtP.x[640];
  localDW->Traj_loc[641].f1[1] = localDW->Derivative2 * rtP.x[641];
  localDW->Traj_loc[642].f1[1] = localDW->Derivative2 * rtP.x[642];
  localDW->Traj_loc[643].f1[1] = localDW->Derivative2 * rtP.x[643];
  localDW->Traj_loc[644].f1[1] = localDW->Derivative2 * rtP.x[644];
  localDW->Traj_loc[645].f1[1] = localDW->Derivative2 * rtP.x[645];
  localDW->Traj_loc[646].f1[1] = localDW->Derivative2 * rtP.x[646];
  localDW->Traj_loc[647].f1[1] = localDW->Derivative2 * rtP.x[647];
  localDW->Traj_loc[648].f1[1] = localDW->Derivative2 * rtP.x[648];
  localDW->Traj_loc[649].f1[1] = localDW->Derivative2 * rtP.x[649];
  localDW->Traj_loc[650].f1[1] = localDW->Derivative2 * rtP.x[650];
  localDW->Traj_loc[651].f1[1] = localDW->Derivative2 * rtP.x[651];
  localDW->Traj_loc[652].f1[1] = localDW->Derivative2 * rtP.x[652];
  localDW->Traj_loc[653].f1[1] = localDW->Derivative2 * rtP.x[653];
  localDW->Traj_loc[654].f1[1] = localDW->Derivative2 * rtP.x[654];
  localDW->Traj_loc[655].f1[1] = localDW->Derivative2 * rtP.x[655];
  localDW->Traj_loc[656].f1[1] = localDW->Derivative2 * rtP.x[656];
  localDW->Traj_loc[657].f1[1] = localDW->Derivative2 * rtP.x[657];
  localDW->Traj_loc[658].f1[1] = localDW->Derivative2 * rtP.x[658];
  localDW->Traj_loc[659].f1[1] = localDW->Derivative2 * rtP.x[659];
  localDW->Traj_loc[660].f1[1] = localDW->Derivative2 * rtP.x[660];
  localDW->Traj_loc[661].f1[1] = localDW->Derivative2 * rtP.x[661];
  localDW->Traj_loc[662].f1[1] = localDW->Derivative2 * rtP.x[662];
  localDW->Traj_loc[663].f1[1] = localDW->Derivative2 * rtP.x[663];
  localDW->Traj_loc[664].f1[1] = localDW->Derivative2 * rtP.x[664];
  localDW->Traj_loc[665].f1[1] = localDW->Derivative2 * rtP.x[665];
  localDW->Traj_loc[666].f1[1] = localDW->Derivative2 * rtP.x[666];
  localDW->Traj_loc[667].f1[1] = localDW->Derivative2 * rtP.x[667];
  localDW->Traj_loc[668].f1[1] = localDW->Derivative2 * rtP.x[668];
  localDW->Traj_loc[669].f1[1] = localDW->Derivative2 * rtP.x[669];
  localDW->Traj_loc[670].f1[1] = localDW->Derivative2 * rtP.x[670];
  localDW->Traj_loc[671].f1[1] = localDW->Derivative2 * rtP.x[671];
  localDW->Traj_loc[672].f1[1] = localDW->Derivative2 * rtP.x[672];
  localDW->Traj_loc[673].f1[1] = localDW->Derivative2 * rtP.x[673];
  localDW->Traj_loc[674].f1[1] = localDW->Derivative2 * rtP.x[674];
  localDW->Traj_loc[675].f1[1] = localDW->Derivative2 * rtP.x[675];
  localDW->Traj_loc[676].f1[1] = localDW->Derivative2 * rtP.x[676];
  localDW->Traj_loc[677].f1[1] = localDW->Derivative2 * rtP.x[677];
  localDW->Traj_loc[678].f1[1] = localDW->Derivative2 * rtP.x[678];
  localDW->Traj_loc[679].f1[1] = localDW->Derivative2 * rtP.x[679];
  localDW->Traj_loc[680].f1[1] = localDW->Derivative2 * rtP.x[680];
  localDW->Traj_loc[681].f1[1] = localDW->Derivative2 * rtP.x[681];
  localDW->Traj_loc[682].f1[1] = localDW->Derivative2 * rtP.x[682];
  localDW->Traj_loc[683].f1[1] = localDW->Derivative2 * rtP.x[683];
  localDW->Traj_loc[684].f1[1] = localDW->Derivative2 * rtP.x[684];
  localDW->Traj_loc[685].f1[1] = localDW->Derivative2 * rtP.x[685];
  localDW->Traj_loc[686].f1[1] = localDW->Derivative2 * rtP.x[686];
  localDW->Traj_loc[687].f1[1] = localDW->Derivative2 * rtP.x[687];
  localDW->Traj_loc[688].f1[1] = localDW->Derivative2 * rtP.x[688];
  localDW->Traj_loc[689].f1[1] = localDW->Derivative2 * rtP.x[689];
  localDW->Traj_loc[690].f1[1] = localDW->Derivative2 * rtP.x[690];
  localDW->Traj_loc[691].f1[1] = localDW->Derivative2 * rtP.x[691];
  localDW->Traj_loc[692].f1[1] = localDW->Derivative2 * rtP.x[692];
  localDW->Traj_loc[693].f1[1] = localDW->Derivative2 * rtP.x[693];
  localDW->Traj_loc[694].f1[1] = localDW->Derivative2 * rtP.x[694];
  localDW->Traj_loc[695].f1[1] = localDW->Derivative2 * rtP.x[695];
  localDW->Traj_loc[696].f1[1] = localDW->Derivative2 * rtP.x[696];
  localDW->Traj_loc[697].f1[1] = localDW->Derivative2 * rtP.x[697];
  localDW->Traj_loc[698].f1[1] = localDW->Derivative2 * rtP.x[698];
  localDW->Traj_loc[699].f1[1] = localDW->Derivative2 * rtP.x[699];
  localDW->Traj_loc[700].f1[1] = localDW->Derivative2 * rtP.x[700];
  localDW->Traj_loc[701].f1[1] = localDW->Derivative2 * rtP.x[701];
  localDW->Traj_loc[702].f1[1] = localDW->Derivative2 * rtP.x[702];
  localDW->Traj_loc[703].f1[1] = localDW->Derivative2 * rtP.x[703];
  localDW->Traj_loc[704].f1[1] = localDW->Derivative2 * rtP.x[704];
  localDW->Traj_loc[705].f1[1] = localDW->Derivative2 * rtP.x[705];
  localDW->Traj_loc[706].f1[1] = localDW->Derivative2 * rtP.x[706];
  localDW->Traj_loc[707].f1[1] = localDW->Derivative2 * rtP.x[707];
  localDW->Traj_loc[708].f1[1] = localDW->Derivative2 * rtP.x[708];
  localDW->Traj_loc[709].f1[1] = localDW->Derivative2 * rtP.x[709];
  localDW->Traj_loc[710].f1[1] = localDW->Derivative2 * rtP.x[710];
  localDW->Traj_loc[711].f1[1] = localDW->Derivative2 * rtP.x[711];
  localDW->Traj_loc[712].f1[1] = localDW->Derivative2 * rtP.x[712];
  localDW->Traj_loc[713].f1[1] = localDW->Derivative2 * rtP.x[713];
  localDW->Traj_loc[714].f1[1] = localDW->Derivative2 * rtP.x[714];
  localDW->Traj_loc[715].f1[1] = localDW->Derivative2 * rtP.x[715];
  localDW->Traj_loc[716].f1[1] = localDW->Derivative2 * rtP.x[716];
  localDW->Traj_loc[717].f1[1] = localDW->Derivative2 * rtP.x[717];
  localDW->Traj_loc[718].f1[1] = localDW->Derivative2 * rtP.x[718];
  localDW->Traj_loc[719].f1[1] = localDW->Derivative2 * rtP.x[719];
  localDW->Traj_loc[720].f1[1] = localDW->Derivative2 * rtP.x[720];
  localDW->Traj_loc[721].f1[1] = localDW->Derivative2 * rtP.x[721];
  localDW->Traj_loc[722].f1[1] = localDW->Derivative2 * rtP.x[722];
  localDW->Traj_loc[723].f1[1] = localDW->Derivative2 * rtP.x[723];
  localDW->Traj_loc[724].f1[1] = localDW->Derivative2 * rtP.x[724];
  localDW->Traj_loc[725].f1[1] = localDW->Derivative2 * rtP.x[725];
  localDW->Traj_loc[726].f1[1] = localDW->Derivative2 * rtP.x[726];
  localDW->Traj_loc[727].f1[1] = localDW->Derivative2 * rtP.x[727];
  localDW->Traj_loc[728].f1[1] = localDW->Derivative2 * rtP.x[728];
  localDW->Traj_loc[729].f1[1] = localDW->Derivative2 * rtP.x[729];
  localDW->Traj_loc[730].f1[1] = localDW->Derivative2 * rtP.x[730];
  localDW->Traj_loc[731].f1[1] = localDW->Derivative2 * rtP.x[731];
  localDW->Traj_loc[732].f1[1] = localDW->Derivative2 * rtP.x[732];
  localDW->Traj_loc[733].f1[1] = localDW->Derivative2 * rtP.x[733];
  localDW->Traj_loc[734].f1[1] = localDW->Derivative2 * rtP.x[734];
  localDW->Traj_loc[735].f1[1] = localDW->Derivative2 * rtP.x[735];
  localDW->Traj_loc[736].f1[1] = localDW->Derivative2 * rtP.x[736];
  localDW->Traj_loc[737].f1[1] = localDW->Derivative2 * rtP.x[737];
  localDW->Traj_loc[738].f1[1] = localDW->Derivative2 * rtP.x[738];
  localDW->Traj_loc[739].f1[1] = localDW->Derivative2 * rtP.x[739];
  localDW->Traj_loc[740].f1[1] = localDW->Derivative2 * rtP.x[740];
  localDW->Traj_loc[741].f1[1] = localDW->Derivative2 * rtP.x[741];
  localDW->Traj_loc[742].f1[1] = localDW->Derivative2 * rtP.x[742];
  localDW->Traj_loc[743].f1[1] = localDW->Derivative2 * rtP.x[743];
  localDW->Traj_loc[744].f1[1] = localDW->Derivative2 * rtP.x[744];
  localDW->Traj_loc[745].f1[1] = localDW->Derivative2 * rtP.x[745];
  localDW->Traj_loc[746].f1[1] = localDW->Derivative2 * rtP.x[746];
  localDW->Traj_loc[747].f1[1] = localDW->Derivative2 * rtP.x[747];
  localDW->Traj_loc[748].f1[1] = localDW->Derivative2 * rtP.x[748];
  localDW->Traj_loc[749].f1[1] = localDW->Derivative2 * rtP.x[749];
  localDW->Traj_loc[750].f1[1] = localDW->Derivative2 * rtP.x[750];
  localDW->Traj_loc[751].f1[1] = localDW->Derivative2 * rtP.x[751];
  localDW->Traj_loc[752].f1[1] = localDW->Derivative2 * rtP.x[752];
  localDW->Traj_loc[753].f1[1] = localDW->Derivative2 * rtP.x[753];
  localDW->Traj_loc[754].f1[1] = localDW->Derivative2 * rtP.x[754];
  localDW->Traj_loc[755].f1[1] = localDW->Derivative2 * rtP.x[755];
  localDW->Traj_loc[756].f1[1] = localDW->Derivative2 * rtP.x[756];
  localDW->Traj_loc[757].f1[1] = localDW->Derivative2 * rtP.x[757];
  localDW->Traj_loc[758].f1[1] = localDW->Derivative2 * rtP.x[758];
  localDW->Traj_loc[759].f1[1] = localDW->Derivative2 * rtP.x[759];
  localDW->Traj_loc[760].f1[1] = localDW->Derivative2 * rtP.x[760];
  localDW->Traj_loc[761].f1[1] = localDW->Derivative2 * rtP.x[761];
  localDW->Traj_loc[762].f1[1] = localDW->Derivative2 * rtP.x[762];
  localDW->Traj_loc[763].f1[1] = localDW->Derivative2 * rtP.x[763];
  localDW->Traj_loc[764].f1[1] = localDW->Derivative2 * rtP.x[764];
  localDW->Traj_loc[765].f1[1] = localDW->Derivative2 * rtP.x[765];
  localDW->Traj_loc[766].f1[1] = localDW->Derivative2 * rtP.x[766];
  localDW->Traj_loc[767].f1[1] = localDW->Derivative2 * rtP.x[767];
  localDW->Traj_loc[768].f1[1] = localDW->Derivative2 * rtP.x[768];
  localDW->Traj_loc[769].f1[1] = localDW->Derivative2 * rtP.x[769];
  localDW->Traj_loc[770].f1[1] = localDW->Derivative2 * rtP.x[770];
  localDW->Traj_loc[771].f1[1] = localDW->Derivative2 * rtP.x[771];
  localDW->Traj_loc[772].f1[1] = localDW->Derivative2 * rtP.x[772];
  localDW->Traj_loc[773].f1[1] = localDW->Derivative2 * rtP.x[773];
  localDW->Traj_loc[774].f1[1] = localDW->Derivative2 * rtP.x[774];
  localDW->Traj_loc[775].f1[1] = localDW->Derivative2 * rtP.x[775];
  localDW->Traj_loc[776].f1[1] = localDW->Derivative2 * rtP.x[776];
  localDW->Traj_loc[777].f1[1] = localDW->Derivative2 * rtP.x[777];
  localDW->Traj_loc[778].f1[1] = localDW->Derivative2 * rtP.x[778];
  localDW->Traj_loc[779].f1[1] = localDW->Derivative2 * rtP.x[779];
  localDW->Traj_loc[780].f1[1] = localDW->Derivative2 * rtP.x[780];
  localDW->Traj_loc[781].f1[1] = localDW->Derivative2 * rtP.x[781];
  localDW->Traj_loc[782].f1[1] = localDW->Derivative2 * rtP.x[782];
  localDW->Traj_loc[783].f1[1] = localDW->Derivative2 * rtP.x[783];
  localDW->Traj_loc[784].f1[1] = localDW->Derivative2 * rtP.x[784];
  localDW->Traj_loc[785].f1[1] = localDW->Derivative2 * rtP.x[785];
  localDW->Traj_loc[786].f1[1] = localDW->Derivative2 * rtP.x[786];
  localDW->Traj_loc[787].f1[1] = localDW->Derivative2 * rtP.x[787];
  localDW->Traj_loc[788].f1[1] = localDW->Derivative2 * rtP.x[788];
  localDW->Traj_loc[789].f1[1] = localDW->Derivative2 * rtP.x[789];
  localDW->Traj_loc[790].f1[1] = localDW->Derivative2 * rtP.x[790];
  localDW->Traj_loc[791].f1[1] = localDW->Derivative2 * rtP.x[791];
  localDW->Traj_loc[792].f1[1] = localDW->Derivative2 * rtP.x[792];
  localDW->Traj_loc[793].f1[1] = localDW->Derivative2 * rtP.x[793];
  localDW->Traj_loc[794].f1[1] = localDW->Derivative2 * rtP.x[794];
  localDW->Traj_loc[795].f1[1] = localDW->Derivative2 * rtP.x[795];
  localDW->Traj_loc[796].f1[1] = localDW->Derivative2 * rtP.x[796];
  localDW->Traj_loc[797].f1[1] = localDW->Derivative2 * rtP.x[797];
  localDW->Traj_loc[798].f1[1] = localDW->Derivative2 * rtP.x[798];
  localDW->Traj_loc[799].f1[1] = localDW->Derivative2 * rtP.x[799];
  localDW->Traj_loc[800].f1[1] = localDW->Derivative2 * rtP.x[800];
  localDW->Traj_loc[0].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[0];
  localDW->Traj_loc[1].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[1];
  localDW->Traj_loc[2].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[2];
  localDW->Traj_loc[3].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[3];
  localDW->Traj_loc[4].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[4];
  localDW->Traj_loc[5].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[5];
  localDW->Traj_loc[6].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[6];
  localDW->Traj_loc[7].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[7];
  localDW->Traj_loc[8].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[8];
  localDW->Traj_loc[9].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[9];
  localDW->Traj_loc[10].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[10];
  localDW->Traj_loc[11].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[11];
  localDW->Traj_loc[12].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[12];
  localDW->Traj_loc[13].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[13];
  localDW->Traj_loc[14].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[14];
  localDW->Traj_loc[15].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[15];
  localDW->Traj_loc[16].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[16];
  localDW->Traj_loc[17].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[17];
  localDW->Traj_loc[18].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[18];
  localDW->Traj_loc[19].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[19];
  localDW->Traj_loc[20].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[20];
  localDW->Traj_loc[21].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[21];
  localDW->Traj_loc[22].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[22];
  localDW->Traj_loc[23].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[23];
  localDW->Traj_loc[24].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[24];
  localDW->Traj_loc[25].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[25];
  localDW->Traj_loc[26].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[26];
  localDW->Traj_loc[27].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[27];
  localDW->Traj_loc[28].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[28];
  localDW->Traj_loc[29].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[29];
  localDW->Traj_loc[30].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[30];
  localDW->Traj_loc[31].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[31];
  localDW->Traj_loc[32].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[32];
  localDW->Traj_loc[33].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[33];
  localDW->Traj_loc[34].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[34];
  localDW->Traj_loc[35].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[35];
  localDW->Traj_loc[36].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[36];
  localDW->Traj_loc[37].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[37];
  localDW->Traj_loc[38].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[38];
  localDW->Traj_loc[39].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[39];
  localDW->Traj_loc[40].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[40];
  localDW->Traj_loc[41].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[41];
  localDW->Traj_loc[42].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[42];
  localDW->Traj_loc[43].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[43];
  localDW->Traj_loc[44].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[44];
  localDW->Traj_loc[45].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[45];
  localDW->Traj_loc[46].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[46];
  localDW->Traj_loc[47].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[47];
  localDW->Traj_loc[48].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[48];
  localDW->Traj_loc[49].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[49];
  localDW->Traj_loc[50].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[50];
  localDW->Traj_loc[51].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[51];
  localDW->Traj_loc[52].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[52];
  localDW->Traj_loc[53].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[53];
  localDW->Traj_loc[54].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[54];
  localDW->Traj_loc[55].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[55];
  localDW->Traj_loc[56].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[56];
  localDW->Traj_loc[57].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[57];
  localDW->Traj_loc[58].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[58];
  localDW->Traj_loc[59].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[59];
  localDW->Traj_loc[60].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[60];
  localDW->Traj_loc[61].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[61];
  localDW->Traj_loc[62].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[62];
  localDW->Traj_loc[63].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[63];
  localDW->Traj_loc[64].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[64];
  localDW->Traj_loc[65].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[65];
  localDW->Traj_loc[66].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[66];
  localDW->Traj_loc[67].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[67];
  localDW->Traj_loc[68].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[68];
  localDW->Traj_loc[69].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[69];
  localDW->Traj_loc[70].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[70];
  localDW->Traj_loc[71].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[71];
  localDW->Traj_loc[72].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[72];
  localDW->Traj_loc[73].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[73];
  localDW->Traj_loc[74].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[74];
  localDW->Traj_loc[75].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[75];
  localDW->Traj_loc[76].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[76];
  localDW->Traj_loc[77].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[77];
  localDW->Traj_loc[78].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[78];
  localDW->Traj_loc[79].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[79];
  localDW->Traj_loc[80].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[80];
  localDW->Traj_loc[81].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[81];
  localDW->Traj_loc[82].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[82];
  localDW->Traj_loc[83].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[83];
  localDW->Traj_loc[84].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[84];
  localDW->Traj_loc[85].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[85];
  localDW->Traj_loc[86].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[86];
  localDW->Traj_loc[87].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[87];
  localDW->Traj_loc[88].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[88];
  localDW->Traj_loc[89].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[89];
  localDW->Traj_loc[90].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[90];
  localDW->Traj_loc[91].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[91];
  localDW->Traj_loc[92].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[92];
  localDW->Traj_loc[93].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[93];
  localDW->Traj_loc[94].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[94];
  localDW->Traj_loc[95].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[95];
  localDW->Traj_loc[96].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[96];
  localDW->Traj_loc[97].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[97];
  localDW->Traj_loc[98].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[98];
  localDW->Traj_loc[99].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[99];
  localDW->Traj_loc[100].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[100];
  localDW->Traj_loc[101].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[101];
  localDW->Traj_loc[102].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[102];
  localDW->Traj_loc[103].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[103];
  localDW->Traj_loc[104].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[104];
  localDW->Traj_loc[105].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[105];
  localDW->Traj_loc[106].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[106];
  localDW->Traj_loc[107].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[107];
  localDW->Traj_loc[108].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[108];
  localDW->Traj_loc[109].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[109];
  localDW->Traj_loc[110].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[110];
  localDW->Traj_loc[111].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[111];
  localDW->Traj_loc[112].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[112];
  localDW->Traj_loc[113].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[113];
  localDW->Traj_loc[114].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[114];
  localDW->Traj_loc[115].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[115];
  localDW->Traj_loc[116].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[116];
  localDW->Traj_loc[117].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[117];
  localDW->Traj_loc[118].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[118];
  localDW->Traj_loc[119].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[119];
  localDW->Traj_loc[120].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[120];
  localDW->Traj_loc[121].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[121];
  localDW->Traj_loc[122].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[122];
  localDW->Traj_loc[123].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[123];
  localDW->Traj_loc[124].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[124];
  localDW->Traj_loc[125].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[125];
  localDW->Traj_loc[126].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[126];
  localDW->Traj_loc[127].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[127];
  localDW->Traj_loc[128].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[128];
  localDW->Traj_loc[129].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[129];
  localDW->Traj_loc[130].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[130];
  localDW->Traj_loc[131].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[131];
  localDW->Traj_loc[132].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[132];
  localDW->Traj_loc[133].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[133];
  localDW->Traj_loc[134].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[134];
  localDW->Traj_loc[135].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[135];
  localDW->Traj_loc[136].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[136];
  localDW->Traj_loc[137].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[137];
  localDW->Traj_loc[138].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[138];
  localDW->Traj_loc[139].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[139];
  localDW->Traj_loc[140].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[140];
  localDW->Traj_loc[141].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[141];
  localDW->Traj_loc[142].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[142];
  localDW->Traj_loc[143].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[143];
  localDW->Traj_loc[144].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[144];
  localDW->Traj_loc[145].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[145];
  localDW->Traj_loc[146].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[146];
  localDW->Traj_loc[147].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[147];
  localDW->Traj_loc[148].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[148];
  localDW->Traj_loc[149].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[149];
  localDW->Traj_loc[150].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[150];
  localDW->Traj_loc[151].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[151];
  localDW->Traj_loc[152].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[152];
  localDW->Traj_loc[153].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[153];
  localDW->Traj_loc[154].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[154];
  localDW->Traj_loc[155].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[155];
  localDW->Traj_loc[156].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[156];
  localDW->Traj_loc[157].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[157];
  localDW->Traj_loc[158].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[158];
  localDW->Traj_loc[159].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[159];
  localDW->Traj_loc[160].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[160];
  localDW->Traj_loc[161].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[161];
  localDW->Traj_loc[162].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[162];
  localDW->Traj_loc[163].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[163];
  localDW->Traj_loc[164].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[164];
  localDW->Traj_loc[165].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[165];
  localDW->Traj_loc[166].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[166];
  localDW->Traj_loc[167].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[167];
  localDW->Traj_loc[168].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[168];
  localDW->Traj_loc[169].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[169];
  localDW->Traj_loc[170].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[170];
  localDW->Traj_loc[171].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[171];
  localDW->Traj_loc[172].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[172];
  localDW->Traj_loc[173].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[173];
  localDW->Traj_loc[174].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[174];
  localDW->Traj_loc[175].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[175];
  localDW->Traj_loc[176].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[176];
  localDW->Traj_loc[177].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[177];
  localDW->Traj_loc[178].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[178];
  localDW->Traj_loc[179].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[179];
  localDW->Traj_loc[180].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[180];
  localDW->Traj_loc[181].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[181];
  localDW->Traj_loc[182].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[182];
  localDW->Traj_loc[183].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[183];
  localDW->Traj_loc[184].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[184];
  localDW->Traj_loc[185].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[185];
  localDW->Traj_loc[186].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[186];
  localDW->Traj_loc[187].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[187];
  localDW->Traj_loc[188].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[188];
  localDW->Traj_loc[189].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[189];
  localDW->Traj_loc[190].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[190];
  localDW->Traj_loc[191].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[191];
  localDW->Traj_loc[192].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[192];
  localDW->Traj_loc[193].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[193];
  localDW->Traj_loc[194].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[194];
  localDW->Traj_loc[195].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[195];
  localDW->Traj_loc[196].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[196];
  localDW->Traj_loc[197].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[197];
  localDW->Traj_loc[198].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[198];
  localDW->Traj_loc[199].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[199];
  localDW->Traj_loc[200].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[200];
  localDW->Traj_loc[201].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[201];
  localDW->Traj_loc[202].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[202];
  localDW->Traj_loc[203].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[203];
  localDW->Traj_loc[204].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[204];
  localDW->Traj_loc[205].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[205];
  localDW->Traj_loc[206].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[206];
  localDW->Traj_loc[207].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[207];
  localDW->Traj_loc[208].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[208];
  localDW->Traj_loc[209].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[209];
  localDW->Traj_loc[210].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[210];
  localDW->Traj_loc[211].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[211];
  localDW->Traj_loc[212].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[212];
  localDW->Traj_loc[213].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[213];
  localDW->Traj_loc[214].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[214];
  localDW->Traj_loc[215].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[215];
  localDW->Traj_loc[216].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[216];
  localDW->Traj_loc[217].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[217];
  localDW->Traj_loc[218].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[218];
  localDW->Traj_loc[219].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[219];
  localDW->Traj_loc[220].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[220];
  localDW->Traj_loc[221].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[221];
  localDW->Traj_loc[222].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[222];
  localDW->Traj_loc[223].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[223];
  localDW->Traj_loc[224].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[224];
  localDW->Traj_loc[225].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[225];
  localDW->Traj_loc[226].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[226];
  localDW->Traj_loc[227].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[227];
  localDW->Traj_loc[228].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[228];
  localDW->Traj_loc[229].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[229];
  localDW->Traj_loc[230].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[230];
  localDW->Traj_loc[231].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[231];
  localDW->Traj_loc[232].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[232];
  localDW->Traj_loc[233].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[233];
  localDW->Traj_loc[234].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[234];
  localDW->Traj_loc[235].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[235];
  localDW->Traj_loc[236].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[236];
  localDW->Traj_loc[237].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[237];
  localDW->Traj_loc[238].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[238];
  localDW->Traj_loc[239].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[239];
  localDW->Traj_loc[240].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[240];
  localDW->Traj_loc[241].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[241];
  localDW->Traj_loc[242].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[242];
  localDW->Traj_loc[243].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[243];
  localDW->Traj_loc[244].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[244];
  localDW->Traj_loc[245].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[245];
  localDW->Traj_loc[246].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[246];
  localDW->Traj_loc[247].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[247];
  localDW->Traj_loc[248].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[248];
  localDW->Traj_loc[249].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[249];
  localDW->Traj_loc[250].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[250];
  localDW->Traj_loc[251].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[251];
  localDW->Traj_loc[252].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[252];
  localDW->Traj_loc[253].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[253];
  localDW->Traj_loc[254].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[254];
  localDW->Traj_loc[255].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[255];
  localDW->Traj_loc[256].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[256];
  localDW->Traj_loc[257].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[257];
  localDW->Traj_loc[258].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[258];
  localDW->Traj_loc[259].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[259];
  localDW->Traj_loc[260].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[260];
  localDW->Traj_loc[261].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[261];
  localDW->Traj_loc[262].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[262];
  localDW->Traj_loc[263].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[263];
  localDW->Traj_loc[264].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[264];
  localDW->Traj_loc[265].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[265];
  localDW->Traj_loc[266].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[266];
  localDW->Traj_loc[267].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[267];
  localDW->Traj_loc[268].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[268];
  localDW->Traj_loc[269].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[269];
  localDW->Traj_loc[270].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[270];
  localDW->Traj_loc[271].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[271];
  localDW->Traj_loc[272].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[272];
  localDW->Traj_loc[273].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[273];
  localDW->Traj_loc[274].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[274];
  localDW->Traj_loc[275].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[275];
  localDW->Traj_loc[276].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[276];
  localDW->Traj_loc[277].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[277];
  localDW->Traj_loc[278].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[278];
  localDW->Traj_loc[279].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[279];
  localDW->Traj_loc[280].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[280];
  localDW->Traj_loc[281].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[281];
  localDW->Traj_loc[282].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[282];
  localDW->Traj_loc[283].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[283];
  localDW->Traj_loc[284].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[284];
  localDW->Traj_loc[285].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[285];
  localDW->Traj_loc[286].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[286];
  localDW->Traj_loc[287].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[287];
  localDW->Traj_loc[288].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[288];
  localDW->Traj_loc[289].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[289];
  localDW->Traj_loc[290].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[290];
  localDW->Traj_loc[291].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[291];
  localDW->Traj_loc[292].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[292];
  localDW->Traj_loc[293].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[293];
  localDW->Traj_loc[294].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[294];
  localDW->Traj_loc[295].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[295];
  localDW->Traj_loc[296].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[296];
  localDW->Traj_loc[297].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[297];
  localDW->Traj_loc[298].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[298];
  localDW->Traj_loc[299].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[299];
  localDW->Traj_loc[300].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[300];
  localDW->Traj_loc[301].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[301];
  localDW->Traj_loc[302].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[302];
  localDW->Traj_loc[303].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[303];
  localDW->Traj_loc[304].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[304];
  localDW->Traj_loc[305].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[305];
  localDW->Traj_loc[306].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[306];
  localDW->Traj_loc[307].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[307];
  localDW->Traj_loc[308].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[308];
  localDW->Traj_loc[309].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[309];
  localDW->Traj_loc[310].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[310];
  localDW->Traj_loc[311].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[311];
  localDW->Traj_loc[312].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[312];
  localDW->Traj_loc[313].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[313];
  localDW->Traj_loc[314].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[314];
  localDW->Traj_loc[315].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[315];
  localDW->Traj_loc[316].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[316];
  localDW->Traj_loc[317].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[317];
  localDW->Traj_loc[318].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[318];
  localDW->Traj_loc[319].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[319];
  localDW->Traj_loc[320].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[320];
  localDW->Traj_loc[321].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[321];
  localDW->Traj_loc[322].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[322];
  localDW->Traj_loc[323].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[323];
  localDW->Traj_loc[324].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[324];
  localDW->Traj_loc[325].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[325];
  localDW->Traj_loc[326].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[326];
  localDW->Traj_loc[327].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[327];
  localDW->Traj_loc[328].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[328];
  localDW->Traj_loc[329].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[329];
  localDW->Traj_loc[330].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[330];
  localDW->Traj_loc[331].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[331];
  localDW->Traj_loc[332].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[332];
  localDW->Traj_loc[333].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[333];
  localDW->Traj_loc[334].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[334];
  localDW->Traj_loc[335].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[335];
  localDW->Traj_loc[336].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[336];
  localDW->Traj_loc[337].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[337];
  localDW->Traj_loc[338].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[338];
  localDW->Traj_loc[339].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[339];
  localDW->Traj_loc[340].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[340];
  localDW->Traj_loc[341].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[341];
  localDW->Traj_loc[342].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[342];
  localDW->Traj_loc[343].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[343];
  localDW->Traj_loc[344].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[344];
  localDW->Traj_loc[345].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[345];
  localDW->Traj_loc[346].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[346];
  localDW->Traj_loc[347].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[347];
  localDW->Traj_loc[348].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[348];
  localDW->Traj_loc[349].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[349];
  localDW->Traj_loc[350].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[350];
  localDW->Traj_loc[351].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[351];
  localDW->Traj_loc[352].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[352];
  localDW->Traj_loc[353].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[353];
  localDW->Traj_loc[354].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[354];
  localDW->Traj_loc[355].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[355];
  localDW->Traj_loc[356].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[356];
  localDW->Traj_loc[357].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[357];
  localDW->Traj_loc[358].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[358];
  localDW->Traj_loc[359].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[359];
  localDW->Traj_loc[360].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[360];
  localDW->Traj_loc[361].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[361];
  localDW->Traj_loc[362].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[362];
  localDW->Traj_loc[363].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[363];
  localDW->Traj_loc[364].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[364];
  localDW->Traj_loc[365].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[365];
  localDW->Traj_loc[366].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[366];
  localDW->Traj_loc[367].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[367];
  localDW->Traj_loc[368].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[368];
  localDW->Traj_loc[369].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[369];
  localDW->Traj_loc[370].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[370];
  localDW->Traj_loc[371].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[371];
  localDW->Traj_loc[372].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[372];
  localDW->Traj_loc[373].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[373];
  localDW->Traj_loc[374].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[374];
  localDW->Traj_loc[375].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[375];
  localDW->Traj_loc[376].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[376];
  localDW->Traj_loc[377].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[377];
  localDW->Traj_loc[378].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[378];
  localDW->Traj_loc[379].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[379];
  localDW->Traj_loc[380].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[380];
  localDW->Traj_loc[381].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[381];
  localDW->Traj_loc[382].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[382];
  localDW->Traj_loc[383].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[383];
  localDW->Traj_loc[384].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[384];
  localDW->Traj_loc[385].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[385];
  localDW->Traj_loc[386].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[386];
  localDW->Traj_loc[387].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[387];
  localDW->Traj_loc[388].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[388];
  localDW->Traj_loc[389].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[389];
  localDW->Traj_loc[390].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[390];
  localDW->Traj_loc[391].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[391];
  localDW->Traj_loc[392].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[392];
  localDW->Traj_loc[393].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[393];
  localDW->Traj_loc[394].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[394];
  localDW->Traj_loc[395].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[395];
  localDW->Traj_loc[396].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[396];
  localDW->Traj_loc[397].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[397];
  localDW->Traj_loc[398].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[398];
  localDW->Traj_loc[399].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[399];
  localDW->Traj_loc[400].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[400];
  localDW->Traj_loc[401].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[401];
  localDW->Traj_loc[402].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[402];
  localDW->Traj_loc[403].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[403];
  localDW->Traj_loc[404].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[404];
  localDW->Traj_loc[405].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[405];
  localDW->Traj_loc[406].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[406];
  localDW->Traj_loc[407].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[407];
  localDW->Traj_loc[408].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[408];
  localDW->Traj_loc[409].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[409];
  localDW->Traj_loc[410].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[410];
  localDW->Traj_loc[411].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[411];
  localDW->Traj_loc[412].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[412];
  localDW->Traj_loc[413].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[413];
  localDW->Traj_loc[414].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[414];
  localDW->Traj_loc[415].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[415];
  localDW->Traj_loc[416].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[416];
  localDW->Traj_loc[417].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[417];
  localDW->Traj_loc[418].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[418];
  localDW->Traj_loc[419].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[419];
  localDW->Traj_loc[420].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[420];
  localDW->Traj_loc[421].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[421];
  localDW->Traj_loc[422].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[422];
  localDW->Traj_loc[423].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[423];
  localDW->Traj_loc[424].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[424];
  localDW->Traj_loc[425].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[425];
  localDW->Traj_loc[426].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[426];
  localDW->Traj_loc[427].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[427];
  localDW->Traj_loc[428].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[428];
  localDW->Traj_loc[429].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[429];
  localDW->Traj_loc[430].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[430];
  localDW->Traj_loc[431].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[431];
  localDW->Traj_loc[432].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[432];
  localDW->Traj_loc[433].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[433];
  localDW->Traj_loc[434].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[434];
  localDW->Traj_loc[435].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[435];
  localDW->Traj_loc[436].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[436];
  localDW->Traj_loc[437].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[437];
  localDW->Traj_loc[438].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[438];
  localDW->Traj_loc[439].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[439];
  localDW->Traj_loc[440].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[440];
  localDW->Traj_loc[441].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[441];
  localDW->Traj_loc[442].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[442];
  localDW->Traj_loc[443].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[443];
  localDW->Traj_loc[444].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[444];
  localDW->Traj_loc[445].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[445];
  localDW->Traj_loc[446].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[446];
  localDW->Traj_loc[447].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[447];
  localDW->Traj_loc[448].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[448];
  localDW->Traj_loc[449].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[449];
  localDW->Traj_loc[450].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[450];
  localDW->Traj_loc[451].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[451];
  localDW->Traj_loc[452].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[452];
  localDW->Traj_loc[453].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[453];
  localDW->Traj_loc[454].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[454];
  localDW->Traj_loc[455].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[455];
  localDW->Traj_loc[456].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[456];
  localDW->Traj_loc[457].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[457];
  localDW->Traj_loc[458].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[458];
  localDW->Traj_loc[459].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[459];
  localDW->Traj_loc[460].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[460];
  localDW->Traj_loc[461].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[461];
  localDW->Traj_loc[462].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[462];
  localDW->Traj_loc[463].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[463];
  localDW->Traj_loc[464].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[464];
  localDW->Traj_loc[465].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[465];
  localDW->Traj_loc[466].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[466];
  localDW->Traj_loc[467].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[467];
  localDW->Traj_loc[468].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[468];
  localDW->Traj_loc[469].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[469];
  localDW->Traj_loc[470].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[470];
  localDW->Traj_loc[471].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[471];
  localDW->Traj_loc[472].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[472];
  localDW->Traj_loc[473].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[473];
  localDW->Traj_loc[474].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[474];
  localDW->Traj_loc[475].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[475];
  localDW->Traj_loc[476].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[476];
  localDW->Traj_loc[477].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[477];
  localDW->Traj_loc[478].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[478];
  localDW->Traj_loc[479].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[479];
  localDW->Traj_loc[480].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[480];
  localDW->Traj_loc[481].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[481];
  localDW->Traj_loc[482].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[482];
  localDW->Traj_loc[483].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[483];
  localDW->Traj_loc[484].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[484];
  localDW->Traj_loc[485].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[485];
  localDW->Traj_loc[486].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[486];
  localDW->Traj_loc[487].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[487];
  localDW->Traj_loc[488].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[488];
  localDW->Traj_loc[489].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[489];
  localDW->Traj_loc[490].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[490];
  localDW->Traj_loc[491].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[491];
  localDW->Traj_loc[492].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[492];
  localDW->Traj_loc[493].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[493];
  localDW->Traj_loc[494].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[494];
  localDW->Traj_loc[495].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[495];
  localDW->Traj_loc[496].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[496];
  localDW->Traj_loc[497].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[497];
  localDW->Traj_loc[498].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[498];
  localDW->Traj_loc[499].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[499];
  localDW->Traj_loc[500].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[500];
  localDW->Traj_loc[501].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[501];
  localDW->Traj_loc[502].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[502];
  localDW->Traj_loc[503].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[503];
  localDW->Traj_loc[504].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[504];
  localDW->Traj_loc[505].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[505];
  localDW->Traj_loc[506].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[506];
  localDW->Traj_loc[507].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[507];
  localDW->Traj_loc[508].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[508];
  localDW->Traj_loc[509].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[509];
  localDW->Traj_loc[510].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[510];
  localDW->Traj_loc[511].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[511];
  localDW->Traj_loc[512].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[512];
  localDW->Traj_loc[513].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[513];
  localDW->Traj_loc[514].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[514];
  localDW->Traj_loc[515].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[515];
  localDW->Traj_loc[516].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[516];
  localDW->Traj_loc[517].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[517];
  localDW->Traj_loc[518].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[518];
  localDW->Traj_loc[519].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[519];
  localDW->Traj_loc[520].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[520];
  localDW->Traj_loc[521].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[521];
  localDW->Traj_loc[522].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[522];
  localDW->Traj_loc[523].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[523];
  localDW->Traj_loc[524].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[524];
  localDW->Traj_loc[525].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[525];
  localDW->Traj_loc[526].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[526];
  localDW->Traj_loc[527].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[527];
  localDW->Traj_loc[528].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[528];
  localDW->Traj_loc[529].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[529];
  localDW->Traj_loc[530].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[530];
  localDW->Traj_loc[531].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[531];
  localDW->Traj_loc[532].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[532];
  localDW->Traj_loc[533].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[533];
  localDW->Traj_loc[534].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[534];
  localDW->Traj_loc[535].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[535];
  localDW->Traj_loc[536].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[536];
  localDW->Traj_loc[537].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[537];
  localDW->Traj_loc[538].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[538];
  localDW->Traj_loc[539].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[539];
  localDW->Traj_loc[540].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[540];
  localDW->Traj_loc[541].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[541];
  localDW->Traj_loc[542].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[542];
  localDW->Traj_loc[543].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[543];
  localDW->Traj_loc[544].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[544];
  localDW->Traj_loc[545].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[545];
  localDW->Traj_loc[546].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[546];
  localDW->Traj_loc[547].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[547];
  localDW->Traj_loc[548].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[548];
  localDW->Traj_loc[549].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[549];
  localDW->Traj_loc[550].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[550];
  localDW->Traj_loc[551].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[551];
  localDW->Traj_loc[552].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[552];
  localDW->Traj_loc[553].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[553];
  localDW->Traj_loc[554].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[554];
  localDW->Traj_loc[555].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[555];
  localDW->Traj_loc[556].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[556];
  localDW->Traj_loc[557].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[557];
  localDW->Traj_loc[558].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[558];
  localDW->Traj_loc[559].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[559];
  localDW->Traj_loc[560].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[560];
  localDW->Traj_loc[561].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[561];
  localDW->Traj_loc[562].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[562];
  localDW->Traj_loc[563].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[563];
  localDW->Traj_loc[564].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[564];
  localDW->Traj_loc[565].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[565];
  localDW->Traj_loc[566].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[566];
  localDW->Traj_loc[567].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[567];
  localDW->Traj_loc[568].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[568];
  localDW->Traj_loc[569].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[569];
  localDW->Traj_loc[570].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[570];
  localDW->Traj_loc[571].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[571];
  localDW->Traj_loc[572].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[572];
  localDW->Traj_loc[573].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[573];
  localDW->Traj_loc[574].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[574];
  localDW->Traj_loc[575].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[575];
  localDW->Traj_loc[576].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[576];
  localDW->Traj_loc[577].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[577];
  localDW->Traj_loc[578].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[578];
  localDW->Traj_loc[579].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[579];
  localDW->Traj_loc[580].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[580];
  localDW->Traj_loc[581].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[581];
  localDW->Traj_loc[582].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[582];
  localDW->Traj_loc[583].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[583];
  localDW->Traj_loc[584].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[584];
  localDW->Traj_loc[585].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[585];
  localDW->Traj_loc[586].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[586];
  localDW->Traj_loc[587].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[587];
  localDW->Traj_loc[588].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[588];
  localDW->Traj_loc[589].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[589];
  localDW->Traj_loc[590].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[590];
  localDW->Traj_loc[591].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[591];
  localDW->Traj_loc[592].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[592];
  localDW->Traj_loc[593].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[593];
  localDW->Traj_loc[594].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[594];
  localDW->Traj_loc[595].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[595];
  localDW->Traj_loc[596].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[596];
  localDW->Traj_loc[597].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[597];
  localDW->Traj_loc[598].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[598];
  localDW->Traj_loc[599].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[599];
  localDW->Traj_loc[600].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[600];
  localDW->Traj_loc[601].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[601];
  localDW->Traj_loc[602].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[602];
  localDW->Traj_loc[603].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[603];
  localDW->Traj_loc[604].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[604];
  localDW->Traj_loc[605].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[605];
  localDW->Traj_loc[606].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[606];
  localDW->Traj_loc[607].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[607];
  localDW->Traj_loc[608].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[608];
  localDW->Traj_loc[609].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[609];
  localDW->Traj_loc[610].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[610];
  localDW->Traj_loc[611].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[611];
  localDW->Traj_loc[612].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[612];
  localDW->Traj_loc[613].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[613];
  localDW->Traj_loc[614].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[614];
  localDW->Traj_loc[615].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[615];
  localDW->Traj_loc[616].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[616];
  localDW->Traj_loc[617].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[617];
  localDW->Traj_loc[618].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[618];
  localDW->Traj_loc[619].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[619];
  localDW->Traj_loc[620].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[620];
  localDW->Traj_loc[621].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[621];
  localDW->Traj_loc[622].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[622];
  localDW->Traj_loc[623].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[623];
  localDW->Traj_loc[624].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[624];
  localDW->Traj_loc[625].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[625];
  localDW->Traj_loc[626].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[626];
  localDW->Traj_loc[627].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[627];
  localDW->Traj_loc[628].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[628];
  localDW->Traj_loc[629].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[629];
  localDW->Traj_loc[630].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[630];
  localDW->Traj_loc[631].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[631];
  localDW->Traj_loc[632].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[632];
  localDW->Traj_loc[633].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[633];
  localDW->Traj_loc[634].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[634];
  localDW->Traj_loc[635].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[635];
  localDW->Traj_loc[636].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[636];
  localDW->Traj_loc[637].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[637];
  localDW->Traj_loc[638].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[638];
  localDW->Traj_loc[639].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[639];
  localDW->Traj_loc[640].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[640];
  localDW->Traj_loc[641].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[641];
  localDW->Traj_loc[642].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[642];
  localDW->Traj_loc[643].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[643];
  localDW->Traj_loc[644].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[644];
  localDW->Traj_loc[645].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[645];
  localDW->Traj_loc[646].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[646];
  localDW->Traj_loc[647].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[647];
  localDW->Traj_loc[648].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[648];
  localDW->Traj_loc[649].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[649];
  localDW->Traj_loc[650].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[650];
  localDW->Traj_loc[651].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[651];
  localDW->Traj_loc[652].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[652];
  localDW->Traj_loc[653].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[653];
  localDW->Traj_loc[654].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[654];
  localDW->Traj_loc[655].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[655];
  localDW->Traj_loc[656].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[656];
  localDW->Traj_loc[657].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[657];
  localDW->Traj_loc[658].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[658];
  localDW->Traj_loc[659].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[659];
  localDW->Traj_loc[660].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[660];
  localDW->Traj_loc[661].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[661];
  localDW->Traj_loc[662].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[662];
  localDW->Traj_loc[663].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[663];
  localDW->Traj_loc[664].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[664];
  localDW->Traj_loc[665].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[665];
  localDW->Traj_loc[666].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[666];
  localDW->Traj_loc[667].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[667];
  localDW->Traj_loc[668].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[668];
  localDW->Traj_loc[669].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[669];
  localDW->Traj_loc[670].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[670];
  localDW->Traj_loc[671].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[671];
  localDW->Traj_loc[672].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[672];
  localDW->Traj_loc[673].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[673];
  localDW->Traj_loc[674].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[674];
  localDW->Traj_loc[675].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[675];
  localDW->Traj_loc[676].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[676];
  localDW->Traj_loc[677].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[677];
  localDW->Traj_loc[678].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[678];
  localDW->Traj_loc[679].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[679];
  localDW->Traj_loc[680].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[680];
  localDW->Traj_loc[681].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[681];
  localDW->Traj_loc[682].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[682];
  localDW->Traj_loc[683].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[683];
  localDW->Traj_loc[684].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[684];
  localDW->Traj_loc[685].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[685];
  localDW->Traj_loc[686].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[686];
  localDW->Traj_loc[687].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[687];
  localDW->Traj_loc[688].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[688];
  localDW->Traj_loc[689].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[689];
  localDW->Traj_loc[690].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[690];
  localDW->Traj_loc[691].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[691];
  localDW->Traj_loc[692].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[692];
  localDW->Traj_loc[693].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[693];
  localDW->Traj_loc[694].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[694];
  localDW->Traj_loc[695].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[695];
  localDW->Traj_loc[696].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[696];
  localDW->Traj_loc[697].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[697];
  localDW->Traj_loc[698].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[698];
  localDW->Traj_loc[699].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[699];
  localDW->Traj_loc[700].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[700];
  localDW->Traj_loc[701].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[701];
  localDW->Traj_loc[702].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[702];
  localDW->Traj_loc[703].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[703];
  localDW->Traj_loc[704].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[704];
  localDW->Traj_loc[705].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[705];
  localDW->Traj_loc[706].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[706];
  localDW->Traj_loc[707].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[707];
  localDW->Traj_loc[708].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[708];
  localDW->Traj_loc[709].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[709];
  localDW->Traj_loc[710].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[710];
  localDW->Traj_loc[711].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[711];
  localDW->Traj_loc[712].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[712];
  localDW->Traj_loc[713].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[713];
  localDW->Traj_loc[714].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[714];
  localDW->Traj_loc[715].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[715];
  localDW->Traj_loc[716].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[716];
  localDW->Traj_loc[717].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[717];
  localDW->Traj_loc[718].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[718];
  localDW->Traj_loc[719].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[719];
  localDW->Traj_loc[720].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[720];
  localDW->Traj_loc[721].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[721];
  localDW->Traj_loc[722].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[722];
  localDW->Traj_loc[723].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[723];
  localDW->Traj_loc[724].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[724];
  localDW->Traj_loc[725].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[725];
  localDW->Traj_loc[726].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[726];
  localDW->Traj_loc[727].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[727];
  localDW->Traj_loc[728].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[728];
  localDW->Traj_loc[729].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[729];
  localDW->Traj_loc[730].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[730];
  localDW->Traj_loc[731].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[731];
  localDW->Traj_loc[732].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[732];
  localDW->Traj_loc[733].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[733];
  localDW->Traj_loc[734].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[734];
  localDW->Traj_loc[735].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[735];
  localDW->Traj_loc[736].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[736];
  localDW->Traj_loc[737].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[737];
  localDW->Traj_loc[738].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[738];
  localDW->Traj_loc[739].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[739];
  localDW->Traj_loc[740].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[740];
  localDW->Traj_loc[741].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[741];
  localDW->Traj_loc[742].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[742];
  localDW->Traj_loc[743].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[743];
  localDW->Traj_loc[744].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[744];
  localDW->Traj_loc[745].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[745];
  localDW->Traj_loc[746].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[746];
  localDW->Traj_loc[747].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[747];
  localDW->Traj_loc[748].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[748];
  localDW->Traj_loc[749].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[749];
  localDW->Traj_loc[750].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[750];
  localDW->Traj_loc[751].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[751];
  localDW->Traj_loc[752].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[752];
  localDW->Traj_loc[753].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[753];
  localDW->Traj_loc[754].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[754];
  localDW->Traj_loc[755].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[755];
  localDW->Traj_loc[756].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[756];
  localDW->Traj_loc[757].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[757];
  localDW->Traj_loc[758].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[758];
  localDW->Traj_loc[759].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[759];
  localDW->Traj_loc[760].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[760];
  localDW->Traj_loc[761].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[761];
  localDW->Traj_loc[762].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[762];
  localDW->Traj_loc[763].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[763];
  localDW->Traj_loc[764].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[764];
  localDW->Traj_loc[765].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[765];
  localDW->Traj_loc[766].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[766];
  localDW->Traj_loc[767].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[767];
  localDW->Traj_loc[768].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[768];
  localDW->Traj_loc[769].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[769];
  localDW->Traj_loc[770].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[770];
  localDW->Traj_loc[771].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[771];
  localDW->Traj_loc[772].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[772];
  localDW->Traj_loc[773].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[773];
  localDW->Traj_loc[774].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[774];
  localDW->Traj_loc[775].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[775];
  localDW->Traj_loc[776].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[776];
  localDW->Traj_loc[777].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[777];
  localDW->Traj_loc[778].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[778];
  localDW->Traj_loc[779].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[779];
  localDW->Traj_loc[780].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[780];
  localDW->Traj_loc[781].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[781];
  localDW->Traj_loc[782].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[782];
  localDW->Traj_loc[783].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[783];
  localDW->Traj_loc[784].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[784];
  localDW->Traj_loc[785].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[785];
  localDW->Traj_loc[786].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[786];
  localDW->Traj_loc[787].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[787];
  localDW->Traj_loc[788].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[788];
  localDW->Traj_loc[789].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[789];
  localDW->Traj_loc[790].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[790];
  localDW->Traj_loc[791].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[791];
  localDW->Traj_loc[792].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[792];
  localDW->Traj_loc[793].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[793];
  localDW->Traj_loc[794].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[794];
  localDW->Traj_loc[795].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[795];
  localDW->Traj_loc[796].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[796];
  localDW->Traj_loc[797].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[797];
  localDW->Traj_loc[798].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[798];
  localDW->Traj_loc[799].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[799];
  localDW->Traj_loc[800].f1[1] += -localDW->rtb_alpha_idx_3 * rtP.y[800];
  for (localDW->p1 = 0; localDW->p1 < 5; localDW->p1++) {
    localDW->x_ld[localDW->p1] = a[localDW->p1] * localDW->rtb_alpha_idx_3 +
      localDW->p;
    localDW->y_ld[localDW->p1] = a[localDW->p1] * localDW->Derivative2 +
      localDW->p_o;
  }

  localDW->b.f1[0] = localDW->x_ld[0];
  localDW->b.f1[1] = localDW->y_ld[0];
  localDW->c.f1[0] = localDW->x_ld[1];
  localDW->c.f1[1] = localDW->y_ld[1];
  localDW->d.f1[0] = localDW->x_ld[2];
  localDW->d.f1[1] = localDW->y_ld[2];
  localDW->e.f1[0] = localDW->x_ld[3];
  localDW->e.f1[1] = localDW->y_ld[3];
  localDW->f.f1[0] = localDW->x_ld[4];
  localDW->f.f1[1] = localDW->y_ld[4];
  localDW->PP[0] = localDW->b;
  localDW->PP[1] = localDW->c;
  localDW->PP[2] = localDW->d;
  localDW->PP[3] = localDW->e;
  localDW->PP[4] = localDW->f;
  for (localDW->p1 = 0; localDW->p1 < 5; localDW->p1++) {
    localDW->PP_loc[localDW->p1].f1[0] = 0.0;
    localDW->PP_loc[localDW->p1].f1[0] += localDW->rtb_alpha_idx_3 * localDW->
      PP[localDW->p1].f1[0];
    localDW->PP_loc[localDW->p1].f1[0] += localDW->Derivative2 * localDW->
      PP[localDW->p1].f1[1];
    localDW->PP_loc[localDW->p1].f1[1] = 0.0;
    localDW->PP_loc[localDW->p1].f1[1] += localDW->Derivative2 * localDW->
      PP[localDW->p1].f1[0];
    localDW->PP_loc[localDW->p1].f1[1] += -localDW->rtb_alpha_idx_3 *
      localDW->PP[localDW->p1].f1[1];
    localDW->x_ld[localDW->p1] = 0.0;
    localDW->rtb_alpha_idx_1 = 1000.0;
    for (localDW->p2 = 0; localDW->p2 < 801; localDW->p2++) {
      localDW->absx21 = fabs(localDW->Traj_loc[localDW->p2].f1[0] -
        localDW->PP_loc[localDW->p1].f1[0]);
      if (localDW->absx21 < localDW->rtb_alpha_idx_1) {
        localDW->rtb_alpha_idx_0 = localDW->Traj_loc[localDW->p2].f1[1] -
          localDW->PP_loc[localDW->p1].f1[1];
        if (fabs(localDW->rtb_alpha_idx_0) < localDW->rtb_alpha_idx_2) {
          localDW->rtb_alpha_idx_1 = localDW->absx21;
          localDW->x_ld[localDW->p1] = localDW->rtb_alpha_idx_0;
        }
      }
    }
  }

  localDW->rtb_alpha_idx_2 = 0.0;
  for (localDW->p1 = 0; localDW->p1 < 5; localDW->p1++) {
    localDW->rtb_alpha_idx_2 += localDW->x_ld[localDW->p1] * g[localDW->p1];
  }

  localDW->rtb_alpha_idx_0 = rtP.l1 + rtP.l2;
  localDW->rtb_alpha_idx_1 = localDW->vx * localDW->vx;
  localDW->rtb_alpha_idx_3 = localDW->rtb_alpha_idx_0 * localDW->BCD_idx_0;
  *rty_Steeringwheelangle = -((localDW->rtb_alpha_idx_0 - (rtP.l1 *
    localDW->BCD_idx_0 - rtP.l2 * localDW->BCD_idx_0) *
    (localDW->rtb_alpha_idx_1 * rtP.m) / (localDW->rtb_alpha_idx_3 *
    localDW->BCD_idx_0)) * 2.0 / (((rtP.l2 - rtP.l1 * rtP.m *
    localDW->rtb_alpha_idx_1 / localDW->rtb_alpha_idx_3) * 2.0 + 9.0) * 9.0)) *
    localDW->rtb_alpha_idx_2 * 16.0;

  /* End of MATLAB Function: '<S5>/FWS Controller' */

  /* Gain: '<S4>/1//mu2' incorporates:
   *  Constant: '<S4>/hzpc2'
   *  Sum: '<S4>/Sum4'
   */
  localDW->zu2_ddot = ((localDW->Ft[1] - localDW->rtb_Fs_idx_1) - rtP.mu2 *
                       rtP.g) * (rtP.sw_zu2 / rtP.mu2);

  /* Gain: '<S4>/1//mu3' incorporates:
   *  Constant: '<S4>/hzpc3'
   *  Sum: '<S4>/Sum5'
   */
  localDW->zu3_ddot = ((localDW->Ft[2] - localDW->rtb_Fs_idx_2) - rtP.mu3 *
                       rtP.g) * (rtP.sw_zu3 / rtP.mu3);

  /* Gain: '<S4>/1//mu4' incorporates:
   *  Constant: '<S4>/hzp4'
   *  Sum: '<S4>/Sum6'
   */
  localDW->zu4_ddot = ((localDW->Ft[3] - localDW->Derivative3) - rtP.mu4 * rtP.g)
    * (rtP.sw_zu4 / rtP.mu4);

  /* Gain: '<S4>/Gain10' incorporates:
   *  Constant: '<S4>/m*g'
   *  Sum: '<S4>/Sum of Elements2'
   *  Sum: '<S4>/Sum7'
   */
  localDW->Gain10 = ((((localDW->rtb_Fs_idx_0 + localDW->rtb_Fs_idx_1) +
                       localDW->rtb_Fs_idx_2) + localDW->Derivative3) - rtP.m *
                     rtP.g) * (rtP.sw_zs / rtP.m);

  /* Gain: '<S4>/`1//mu1' incorporates:
   *  Constant: '<S4>/hzpc1'
   *  Sum: '<S4>/Sum3'
   */
  localDW->zu1_ddot = ((localDW->Ft[0] - localDW->rtb_Fs_idx_0) - rtP.mu1 *
                       rtP.g) * (rtP.sw_zu1 / rtP.mu1);
  if (rtmIsMajorTimeStep(rtM)) {
    /* MATLAB Function: '<S113>/MATLAB Function' incorporates:
     *  Constant: '<S113>/Constant'
     *  Constant: '<S113>/Constant1'
     *  Constant: '<S113>/Constant2'
     *  Constant: '<S113>/Constant3'
     *  Constant: '<S113>/Constant4'
     */
    localDW->rtb_Fs_idx_0 = -(1.0 / rtP.Tc);
    memset(&localDW->B_c[0], 0, sizeof(real_T) << 4U);
    localDW->B_c[0] = -localDW->rtb_Fs_idx_0;
    localDW->B_c[1] = -0.0;
    localDW->B_c[4] = -0.0;
    localDW->B_c[5] = -0.0;
    localDW->Ft[0] = 1.0;
    localDW->Ft[1] = 0.0;
    localDW->Ft[2] = 0.0;
    localDW->Ft[3] = 1.0;
    localDW->rtb_Fs_idx_1 = rtP.B * rtP.B * 2.0 * 0.69314718055994529 /
      (0.59857471687193387 * rtP.Tc);
    localDW->rtb_Fs_idx_2 = rtP.K * rtP.K;
    for (localDW->p1 = 0; localDW->p1 < 2; localDW->p1++) {
      localDW->rtb_Ft_c[localDW->p1] = 0.0;
      localDW->rtb_Ft_c[localDW->p1] += localDW->Ft[localDW->p1] *
        localDW->rtb_Fs_idx_1;
      localDW->rtb_Ft_c[localDW->p1 + 2] = 0.0;
      localDW->rtb_Ft_c[localDW->p1 + 2] += localDW->Ft[localDW->p1 + 2] *
        localDW->rtb_Fs_idx_2;
      localDW->B_c[localDW->p1 + 8] = 0.0;
      localDW->B_c[localDW->p1 + 8] += localDW->rtb_Ft_c[localDW->p1];
      localDW->B_c[localDW->p1 + 8] += localDW->rtb_Ft_c[localDW->p1 + 2] * 0.0;
      localDW->B_c[localDW->p1 + 12] = 0.0;
      localDW->B_c[localDW->p1 + 12] += localDW->rtb_Ft_c[localDW->p1] * 0.0;
      localDW->B_c[localDW->p1 + 12] += localDW->rtb_Ft_c[localDW->p1 + 2];
    }

    localDW->B_c[10] = localDW->rtb_Fs_idx_0;
    localDW->B_c[11] = 0.0;
    localDW->B_c[14] = 0.0;
    localDW->B_c[15] = 0.0;
    for (localDW->p1 = 0; localDW->p1 < 16; localDW->p1++) {
      localDW->B_c_m[localDW->p1] = localDW->B_c[localDW->p1] / rtP.Fs;
    }

    expm(localDW->B_c_m, localDW->B_c, localDW);
    localDW->Phi[0] = localDW->B_c[10];
    localDW->Phi[1] = localDW->B_c[11];
    localDW->Phi[2] = localDW->B_c[14];
    localDW->Phi[3] = localDW->B_c[15];
    localDW->Qbk = localDW->B_c[8];
    localDW->Qkk = localDW->B_c[13];
    localDW->Qnk = rtP.N * rtP.N * rtP.Fs;

    /* End of MATLAB Function: '<S113>/MATLAB Function' */

    /* Memory: '<S113>/Memory' */
    localDW->Z_b1 = localDW->Memory_PreviousInput;

    /* Memory: '<S113>/Memory1' */
    localDW->Z_b2 = localDW->Memory1_PreviousInput;

    /* MATLAB Function: '<S113>/MATLAB Function1' incorporates:
     *  Constant: '<S113>/Constant'
     */
    localDW->rtb_Fs_idx_0 = 1.0 / rtP.Fs;
    localDW->rtb_alpha_idx_2 = localDW->absx11;
    if (localDW->rtb_Fs_idx_0 == 0.0) {
      if (localDW->absx11 == 0.0) {
        localDW->rtb_alpha_idx_2 = localDW->rtb_Fs_idx_0;
      }
    } else if (rtIsNaN(localDW->absx11)) {
      localDW->rtb_alpha_idx_2 = (rtNaN);
    } else if (rtIsNaN(localDW->rtb_Fs_idx_0)) {
      localDW->rtb_alpha_idx_2 = (rtNaN);
    } else if (rtIsInf(localDW->absx11)) {
      localDW->rtb_alpha_idx_2 = (rtNaN);
    } else if (localDW->absx11 == 0.0) {
      localDW->rtb_alpha_idx_2 = 0.0 / localDW->rtb_Fs_idx_0;
    } else if (rtIsInf(localDW->rtb_Fs_idx_0)) {
      if ((localDW->rtb_Fs_idx_0 < 0.0) != (localDW->absx11 < 0.0)) {
        localDW->rtb_alpha_idx_2 = localDW->rtb_Fs_idx_0;
      }
    } else {
      localDW->rtb_alpha_idx_2 = fmod(localDW->absx11, localDW->rtb_Fs_idx_0);
      rEQ0 = (localDW->rtb_alpha_idx_2 == 0.0);
      if ((!rEQ0) && (localDW->rtb_Fs_idx_0 > floor(localDW->rtb_Fs_idx_0))) {
        localDW->rtb_alpha_idx_1 = fabs(localDW->absx11 / localDW->rtb_Fs_idx_0);
        rEQ0 = !(fabs(localDW->rtb_alpha_idx_1 - floor(localDW->rtb_alpha_idx_1
                   + 0.5)) > 2.2204460492503131E-16 * localDW->rtb_alpha_idx_1);
      }

      if (rEQ0) {
        localDW->rtb_alpha_idx_2 = localDW->rtb_Fs_idx_0 * 0.0;
      } else if ((localDW->absx11 < 0.0) != (localDW->rtb_Fs_idx_0 < 0.0)) {
        localDW->rtb_alpha_idx_2 += localDW->rtb_Fs_idx_0;
      }
    }

    if (localDW->rtb_alpha_idx_2 == 0.0) {
      localDW->rtb_alpha_idx_2 = sqrt(localDW->Qbk);
      localDW->absx11 = sqrt(localDW->Qkk);
      localDW->rtb_Fs_idx_0 = sqrt(localDW->Qnk);
      localDW->rtb_alpha_idx_2 *= randn(localDW);
      localDW->absx11 *= randn(localDW);
      localDW->rtb_Fs_idx_0 *= randn(localDW);
      localDW->rtb_Fs_idx_1 = (localDW->Phi[0] * localDW->Z_b1 + localDW->Phi[2]
        * localDW->Z_b2) + localDW->rtb_alpha_idx_2;
      localDW->rtb_alpha_idx_0 = localDW->rtb_Fs_idx_1;
      localDW->rtb_Fs_idx_2 = localDW->rtb_Fs_idx_1;
      localDW->rtb_Fs_idx_1 = (localDW->Phi[1] * localDW->Z_b1 + localDW->Phi[3]
        * localDW->Z_b2) + localDW->absx11;
      localDW->Z_b1 = localDW->rtb_Fs_idx_2;
      localDW->Z_b2 = localDW->rtb_Fs_idx_1;
      localDW->rtb_Fs_idx_0 += localDW->rtb_alpha_idx_0 + localDW->rtb_Fs_idx_1;
    } else {
      localDW->rtb_Fs_idx_0 = 0.0;
    }

    /* Gain: '<S3>/Gain' incorporates:
     *  Gain: '<S113>/Gain'
     *  MATLAB Function: '<S113>/MATLAB Function1'
     */
    localDW->radnoise = localP->Gain_Gain * localDW->rtb_Fs_idx_0 *
      localP->Gain_Gain_h;
  }

  /* Sum: '<S3>/Sum' incorporates:
   *  Integrator: '<S4>/Integrator4'
   */
  localDW->Sum = localX->Integrator4_CSTATE + localDW->radnoise;
}

/* Update for atomic system: '<Root>/VandD' */
static void VandD_Update(RT_MODEL * const rtM, DW_VandD *localDW)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (localDW->TimeStampA == (rtInf)) {
    localDW->TimeStampA = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeA;
  } else if (localDW->TimeStampB == (rtInf)) {
    localDW->TimeStampB = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeB;
  } else if (localDW->TimeStampA < localDW->TimeStampB) {
    localDW->TimeStampA = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeA;
  } else {
    localDW->TimeStampB = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeB;
  }

  *lastU = localDW->Dz1;

  /* End of Update for Derivative: '<S4>/Derivative' */

  /* Update for Derivative: '<S4>/Derivative1' */
  if (localDW->TimeStampA_l == (rtInf)) {
    localDW->TimeStampA_l = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeA_k;
  } else if (localDW->TimeStampB_o == (rtInf)) {
    localDW->TimeStampB_o = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeB_a;
  } else if (localDW->TimeStampA_l < localDW->TimeStampB_o) {
    localDW->TimeStampA_l = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeA_k;
  } else {
    localDW->TimeStampB_o = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeB_a;
  }

  *lastU = localDW->Dz2;

  /* End of Update for Derivative: '<S4>/Derivative1' */

  /* Update for Derivative: '<S4>/Derivative2' */
  if (localDW->TimeStampA_p == (rtInf)) {
    localDW->TimeStampA_p = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeA_n;
  } else if (localDW->TimeStampB_n == (rtInf)) {
    localDW->TimeStampB_n = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeB_e;
  } else if (localDW->TimeStampA_p < localDW->TimeStampB_n) {
    localDW->TimeStampA_p = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeA_n;
  } else {
    localDW->TimeStampB_n = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeB_e;
  }

  *lastU = localDW->Dz3;

  /* End of Update for Derivative: '<S4>/Derivative2' */

  /* Update for Derivative: '<S4>/Derivative3' */
  if (localDW->TimeStampA_o == (rtInf)) {
    localDW->TimeStampA_o = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeA_e;
  } else if (localDW->TimeStampB_l == (rtInf)) {
    localDW->TimeStampB_l = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeB_h;
  } else if (localDW->TimeStampA_o < localDW->TimeStampB_l) {
    localDW->TimeStampA_o = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeA_e;
  } else {
    localDW->TimeStampB_l = rtM->Timing.t[0];
    lastU = &localDW->LastUAtTimeB_h;
  }

  *lastU = localDW->Dz4;

  /* End of Update for Derivative: '<S4>/Derivative3' */
  if (rtmIsMajorTimeStep(rtM)) {
    /* Update for Memory: '<S113>/Memory' */
    localDW->Memory_PreviousInput = localDW->Z_b1;

    /* Update for Memory: '<S113>/Memory1' */
    localDW->Memory1_PreviousInput = localDW->Z_b2;
  }
}

/* Derivatives for atomic system: '<Root>/VandD' */
static void VandD_Deriv(DW_VandD *localDW, P_VandD *localP, X_VandD *localX,
  XDot_VandD *localXdot)
{
  /* Derivatives for Integrator: '<S4>/Integrator20' */
  localXdot->Integrator20_CSTATE = localDW->euler_dot1;

  /* Derivatives for Integrator: '<S4>/Integrator21' */
  localXdot->Integrator21_CSTATE = localDW->euler_dot2;

  /* Derivatives for Integrator: '<S4>/Integrator22' */
  localXdot->Integrator22_CSTATE = localDW->euler_dot3;

  /* Derivatives for Integrator: '<S4>/Integrator' */
  localXdot->Integrator_CSTATE = localDW->Gain6;

  /* Derivatives for Integrator: '<S4>/Integrator2' */
  localXdot->Integrator2_CSTATE = localDW->Gain3;

  /* Derivatives for Integrator: '<S4>/Integrator1' */
  localXdot->Integrator1_CSTATE = localDW->Gain5;

  /* Derivatives for Integrator: '<S4>/Integrator15' */
  localXdot->Integrator15_CSTATE = localDW->zs_dot;

  /* Derivatives for Integrator: '<S4>/Integrator14' */
  localXdot->Integrator14_CSTATE = localDW->Gain10;

  /* Derivatives for Integrator: '<S4>/Integrator11' */
  localXdot->Integrator11_CSTATE = localDW->Integrator10;

  /* Derivatives for Integrator: '<S4>/Integrator10' */
  localXdot->Integrator10_CSTATE = localDW->zu3_ddot;

  /* Derivatives for Integrator: '<S4>/Integrator13' */
  localXdot->Integrator13_CSTATE = localDW->Integrator12;

  /* Derivatives for Integrator: '<S4>/Integrator12' */
  localXdot->Integrator12_CSTATE = localDW->zu4_ddot;

  /* Derivatives for Integrator: '<S4>/Integrator7' */
  localXdot->Integrator7_CSTATE = localDW->Integrator3;

  /* Derivatives for Integrator: '<S4>/Integrator3' */
  localXdot->Integrator3_CSTATE = localDW->zu1_ddot;

  /* Derivatives for Integrator: '<S4>/Integrator9' */
  localXdot->Integrator9_CSTATE = localDW->Integrator8;

  /* Derivatives for Integrator: '<S4>/Integrator8' */
  localXdot->Integrator8_CSTATE = localDW->zu2_ddot;

  /* Derivatives for Integrator: '<S48>/Integrator' */
  localXdot->Integrator_CSTATE_p = localDW->IProdOut;

  /* Derivatives for Integrator: '<S43>/Filter' */
  localXdot->Filter_CSTATE = localDW->NProdOut;

  /* Derivatives for Integrator: '<S96>/Integrator' */
  localXdot->Integrator_CSTATE_b = localDW->IntegralGain;

  /* Derivatives for Integrator: '<S91>/Filter' */
  localXdot->Filter_CSTATE_n = localDW->FilterCoefficient;

  /* Derivatives for Integrator: '<S4>/Integrator19' */
  localXdot->Integrator19_CSTATE = localDW->Gain12;

  /* Derivatives for Integrator: '<S4>/Integrator17' */
  localXdot->Integrator17_CSTATE = localDW->Integrator16;

  /* Derivatives for Integrator: '<S4>/Integrator16' */
  localXdot->Integrator16_CSTATE = localDW->Gain11;

  /* Derivatives for Integrator: '<S4>/Integrator18' */
  localXdot->Integrator18_CSTATE = localDW->Integrator19;

  /* Derivatives for Integrator: '<S4>/Integrator5' */
  localXdot->Integrator5_CSTATE = localDW->Product9[0];

  /* Derivatives for Integrator: '<S4>/Integrator6' */
  localXdot->Integrator6_CSTATE = localDW->Product9[1];

  /* Derivatives for TransferFcn: '<S3>/Transfer Fcn' */
  localXdot->TransferFcn_CSTATE = localP->TransferFcn_A *
    localX->TransferFcn_CSTATE;
  localXdot->TransferFcn_CSTATE += localDW->Sum;

  /* Derivatives for Integrator: '<S4>/Integrator4' */
  localXdot->Integrator4_CSTATE = localDW->psi_dot;

  /* Derivatives for Integrator: '<S4>/Integrator23' */
  localXdot->Integrator23_CSTATE = localDW->Product9[2];
}

/* Model step function */
void VandD_step(void)
{
  if (rtmIsMajorTimeStep(rtM)) {
    /* set solver stop time */
    rtsiSetSolverStopTime(&rtM->solverInfo,((rtM->Timing.clockTick0+1)*
      rtM->Timing.stepSize0));
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(rtM)) {
    rtM->Timing.t[0] = rtsiGetT(&rtM->solverInfo);
  }

  /* Outputs for Atomic SubSystem: '<Root>/VandD' */

  /* Outport: '<Root>/Steering wheel angle' incorporates:
   *  Inport: '<Root>/lamda'
   */
  VandD_g(rtM, rtU.lamda, rtDW.Fy, &rtY.Steeringwheelangle, &rtDW.VandD_gs,
          &rtP.VandD_gs, &rtX.VandD_gs);

  /* End of Outputs for SubSystem: '<Root>/VandD' */

  /* Outport: '<Root>/Fy1' */
  rtY.Fy1 = rtDW.Fy[0];

  /* Outport: '<Root>/Fy2' */
  rtY.Fy2 = rtDW.Fy[1];

  /* Outport: '<Root>/Fy3' */
  rtY.Fy3 = rtDW.Fy[2];

  /* Outport: '<Root>/Fy4' */
  rtY.Fy4 = rtDW.Fy[3];

  /* Outport: '<Root>/throttle force' */
  rtY.throttleforce[0] = rtDW.VandD_gs.ManualSwitch[0];
  rtY.throttleforce[1] = rtDW.VandD_gs.ManualSwitch[1];
  rtY.throttleforce[2] = rtDW.VandD_gs.ManualSwitch[2];
  rtY.throttleforce[3] = rtDW.VandD_gs.ManualSwitch[3];

  /* Outport: '<Root>/Actual Velocity' */
  rtY.ActualVelocity = rtDW.VandD_gs.vx;

  /* Outport: '<Root>/Lateral Acceleration' */
  rtY.LateralAcceleration = rtDW.VandD_gs.vx;

  /* Outport: '<Root>/Angular Velocity' */
  rtY.AngularVelocity = rtDW.VandD_gs.psi_dot;

  /* Outport: '<Root>/Lateral Velocity' */
  rtY.LateralVelocity = rtDW.VandD_gs.vy;

  /* Outport: '<Root>/x_position' */
  rtY.x_position = rtDW.VandD_gs.p;

  /* Outport: '<Root>/y_position' */
  rtY.y_position = rtDW.VandD_gs.p_o;

  /* Outport: '<Root>/orientation' */
  rtY.orientation = rtDW.VandD_gs.TransferFcn;
  if (rtmIsMajorTimeStep(rtM)) {
    /* Update for Atomic SubSystem: '<Root>/VandD' */
    VandD_Update(rtM, &rtDW.VandD_gs);

    /* End of Update for SubSystem: '<Root>/VandD' */
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(rtM)) {
    rt_ertODEUpdateContinuousStates(&rtM->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     */
    ++rtM->Timing.clockTick0;
    rtM->Timing.t[0] = rtsiGetSolverStopTime(&rtM->solverInfo);

    {
      /* Update absolute timer for sample time: [0.01s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.01, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       */
      rtM->Timing.clockTick1++;
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void VandD_derivatives(void)
{
  XDot *_rtXdot;
  _rtXdot = ((XDot *) rtM->derivs);

  /* Derivatives for Atomic SubSystem: '<Root>/VandD' */
  VandD_Deriv(&rtDW.VandD_gs, &rtP.VandD_gs, &rtX.VandD_gs, &_rtXdot->VandD_gs);

  /* End of Derivatives for SubSystem: '<Root>/VandD' */
}

/* Model initialize function */
void VandD_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&rtM->solverInfo, &rtM->Timing.simTimeStep);
    rtsiSetTPtr(&rtM->solverInfo, &rtmGetTPtr(rtM));
    rtsiSetStepSizePtr(&rtM->solverInfo, &rtM->Timing.stepSize0);
    rtsiSetdXPtr(&rtM->solverInfo, &rtM->derivs);
    rtsiSetContStatesPtr(&rtM->solverInfo, (real_T **) &rtM->contStates);
    rtsiSetNumContStatesPtr(&rtM->solverInfo, &rtM->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&rtM->solverInfo,
      &rtM->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&rtM->solverInfo,
      &rtM->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&rtM->solverInfo,
      &rtM->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&rtM->solverInfo, (&rtmGetErrorStatus(rtM)));
    rtsiSetRTModelPtr(&rtM->solverInfo, rtM);
  }

  rtsiSetSimTimeStep(&rtM->solverInfo, MAJOR_TIME_STEP);
  rtM->intgData.y = rtM->odeY;
  rtM->intgData.f[0] = rtM->odeF[0];
  rtM->intgData.f[1] = rtM->odeF[1];
  rtM->intgData.f[2] = rtM->odeF[2];
  rtM->intgData.f[3] = rtM->odeF[3];
  rtM->contStates = ((X *) &rtX);
  rtsiSetSolverData(&rtM->solverInfo, (void *)&rtM->intgData);
  rtsiSetIsMinorTimeStepWithModeChange(&rtM->solverInfo, false);
  rtsiSetSolverName(&rtM->solverInfo,"ode4");
  rtmSetTPtr(rtM, &rtM->Timing.tArray[0]);
  rtM->Timing.stepSize0 = 0.01;

  /* SystemInitialize for Atomic SubSystem: '<Root>/VandD' */
  VandD_Init(&rtDW.VandD_gs, &rtP.VandD_gs, &rtX.VandD_gs);

  /* End of SystemInitialize for SubSystem: '<Root>/VandD' */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
