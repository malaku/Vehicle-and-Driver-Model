/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: VandD_data.c
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

/* Block parameters (default storage) */
P rtP = {
  /* Variable: B
   * Referenced by: '<S113>/Constant3'
   */
  0.00416,

  /* Variable: Fs
   * Referenced by: '<S113>/Constant'
   */
  100.0,

  /* Variable: Ix
   * Referenced by: '<S4>/Gain12'
   */
  616.00993684887271,

  /* Variable: Ix_Iy
   * Referenced by: '<S4>/Ix-Iy'
   */
  -852.92722689800007,

  /* Variable: Iy
   * Referenced by: '<S4>/Gain11'
   */
  1468.9371637468728,

  /* Variable: Iy_Iz
   * Referenced by: '<S4>/Iy-Iz'
   */
  281.93716374687278,

  /* Variable: Iz
   * Referenced by: '<S4>/Gain3'
   */
  1187.0,

  /* Variable: Iz_Ix
   * Referenced by: '<S4>/Iz-Ix'
   */
  570.99006315112729,

  /* Variable: K
   * Referenced by: '<S113>/Constant2'
   */
  0.0014,

  /* Variable: Model
   * Referenced by: '<S4>/MATLAB Function2'
   */
  3.0,

  /* Variable: N
   * Referenced by: '<S113>/Constant1'
   */
  0.00166667,

  /* Variable: THETAo
   * Referenced by: '<S4>/Integrator4'
   */
  0.0,

  /* Variable: Tc
   * Referenced by: '<S113>/Constant4'
   */
  32.0,

  /* Variable: Xo
   * Referenced by: '<S4>/Integrator5'
   */
  0.0,

  /* Variable: Yo
   * Referenced by:
   *   '<S4>/Integrator23'
   *   '<S4>/Integrator6'
   */
  0.0,

  /* Variable: b1
   * Referenced by: '<S4>/MATLAB Function'
   */
  0.776,

  /* Variable: b2
   * Referenced by: '<S4>/MATLAB Function'
   */
  0.776,

  /* Variable: b3
   * Referenced by: '<S4>/MATLAB Function'
   */
  0.78,

  /* Variable: b4
   * Referenced by: '<S4>/MATLAB Function'
   */
  0.78,

  /* Variable: dpc
   * Referenced by:
   *   '<S4>/drc'
   *   '<S4>/m*g*dpc'
   */
  0.48781449179999997,

  /* Variable: drc
   * Referenced by:
   *   '<S4>/dpc'
   *   '<S4>/m*g*drc'
   */
  0.48781449179999997,

  /* Variable: g
   * Referenced by:
   *   '<S4>/hzp4'
   *   '<S4>/hzpc1'
   *   '<S4>/hzpc2'
   *   '<S4>/hzpc3'
   *   '<S4>/m*g'
   *   '<S4>/m*g*dpc'
   *   '<S4>/m*g*drc'
   */
  9.81,

  /* Variable: l1
   * Referenced by:
   *   '<S4>/MATLAB Function'
   *   '<S4>/MATLAB Function1'
   *   '<S5>/FWS Controller'
   */
  1.12695753,

  /* Variable: l2
   * Referenced by:
   *   '<S4>/MATLAB Function'
   *   '<S4>/MATLAB Function1'
   *   '<S5>/FWS Controller'
   */
  1.59304247,

  /* Variable: m
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
  1609.23,

  /* Variable: mu1
   * Referenced by:
   *   '<S4>/hzpc1'
   *   '<S4>/`1//mu1'
   */
  49.5,

  /* Variable: mu2
   * Referenced by:
   *   '<S4>/hzpc2'
   *   '<S4>/1//mu2'
   */
  49.5,

  /* Variable: mu3
   * Referenced by:
   *   '<S4>/hzpc3'
   *   '<S4>/1//mu3'
   */
  42.0,

  /* Variable: mu4
   * Referenced by:
   *   '<S4>/hzp4'
   *   '<S4>/1//mu4'
   */
  42.0,

  /* Variable: radius
   * Referenced by: '<S5>/FWS Controller'
   */
  60.0,

  /* Variable: ro
   * Referenced by: '<S4>/Integrator2'
   */
  0.0,

  /* Variable: sw_delta_f
   * Referenced by:
   *   '<S8>/Constant2'
   *   '<S9>/Constant2'
   */
  1.0,

  /* Variable: sw_phi
   * Referenced by: '<S4>/Gain12'
   */
  1.0,

  /* Variable: sw_psi
   * Referenced by: '<S4>/Gain3'
   */
  1.0,

  /* Variable: sw_theta
   * Referenced by: '<S4>/Gain11'
   */
  1.0,

  /* Variable: sw_x
   * Referenced by: '<S4>/Gain6'
   */
  1.0,

  /* Variable: sw_y
   * Referenced by: '<S4>/Gain5'
   */
  1.0,

  /* Variable: sw_zs
   * Referenced by: '<S4>/Gain10'
   */
  1.0,

  /* Variable: sw_zu1
   * Referenced by: '<S4>/`1//mu1'
   */
  1.0,

  /* Variable: sw_zu2
   * Referenced by: '<S4>/1//mu2'
   */
  1.0,

  /* Variable: sw_zu3
   * Referenced by: '<S4>/1//mu3'
   */
  1.0,

  /* Variable: sw_zu4
   * Referenced by: '<S4>/1//mu4'
   */
  1.0,

  /* Variable: vxo
   * Referenced by: '<S4>/Integrator'
   */
  10.0,

  /* Variable: vyo
   * Referenced by: '<S4>/Integrator1'
   */
  0.0,

  /* Variable: x
   * Referenced by: '<S5>/FWS Controller'
   */
  { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25,
    3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75,
    7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0,
    10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0,
    13.25, 13.5, 13.75, 14.0, 14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0,
    16.25, 16.5, 16.75, 17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0,
    19.25, 19.5, 19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0,
    22.25, 22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
    25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75, 28.0,
    28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5, 30.75, 31.0,
    31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25, 33.5, 33.75, 34.0,
    34.25, 34.5, 34.75, 35.0, 35.25, 35.5, 35.75, 36.0, 36.25, 36.5, 36.75, 37.0,
    37.25, 37.5, 37.75, 38.0, 38.25, 38.5, 38.75, 39.0, 39.25, 39.5, 39.75, 40.0,
    40.25, 40.5, 40.75, 41.0, 41.25, 41.5, 41.75, 42.0, 42.25, 42.5, 42.75, 43.0,
    43.25, 43.5, 43.75, 44.0, 44.25, 44.5, 44.75, 45.0, 45.25, 45.5, 45.75, 46.0,
    46.25, 46.5, 46.75, 47.0, 47.25, 47.5, 47.75, 48.0, 48.25, 48.5, 48.75, 49.0,
    49.25, 49.5, 49.75, 50.0, 50.25, 50.5, 50.75, 51.0, 51.25, 51.5, 51.75, 52.0,
    52.25, 52.5, 52.75, 53.0, 53.25, 53.5, 53.75, 54.0, 54.25, 54.5, 54.75, 55.0,
    55.25, 55.5, 55.75, 56.0, 56.25, 56.5, 56.75, 57.0, 57.25, 57.5, 57.75, 58.0,
    58.25, 58.5, 58.75, 59.0, 59.25, 59.5, 59.75, 60.0, 60.25, 60.5, 60.75, 61.0,
    61.25, 61.5, 61.75, 62.0, 62.25, 62.5, 62.75, 63.0, 63.25, 63.5, 63.75, 64.0,
    64.25, 64.5, 64.75, 65.0, 65.25, 65.5, 65.75, 66.0, 66.25, 66.5, 66.75, 67.0,
    67.25, 67.5, 67.75, 68.0, 68.25, 68.5, 68.75, 69.0, 69.25, 69.5, 69.75, 70.0,
    70.25, 70.5, 70.75, 71.0, 71.25, 71.5, 71.75, 72.0, 72.25, 72.5, 72.75, 73.0,
    73.25, 73.5, 73.75, 74.0, 74.25, 74.5, 74.75, 75.0, 75.25, 75.5, 75.75, 76.0,
    76.25, 76.5, 76.75, 77.0, 77.25, 77.5, 77.75, 78.0, 78.25, 78.5, 78.75, 79.0,
    79.25, 79.5, 79.75, 80.0, 80.25, 80.5, 80.75, 81.0, 81.25, 81.5, 81.75, 82.0,
    82.25, 82.5, 82.75, 83.0, 83.25, 83.5, 83.75, 84.0, 84.25, 84.5, 84.75, 85.0,
    85.25, 85.5, 85.75, 86.0, 86.25, 86.5, 86.75, 87.0, 87.25, 87.5, 87.75, 88.0,
    88.25, 88.5, 88.75, 89.0, 89.25, 89.5, 89.75, 90.0, 90.25, 90.5, 90.75, 91.0,
    91.25, 91.5, 91.75, 92.0, 92.25, 92.5, 92.75, 93.0, 93.25, 93.5, 93.75, 94.0,
    94.25, 94.5, 94.75, 95.0, 95.25, 95.5, 95.75, 96.0, 96.25, 96.5, 96.75, 97.0,
    97.25, 97.5, 97.75, 98.0, 98.25, 98.5, 98.75, 99.0, 99.25, 99.5, 99.75,
    100.0, 100.25, 100.5, 100.75, 101.0, 101.25, 101.5, 101.75, 102.0, 102.25,
    102.5, 102.75, 103.0, 103.25, 103.5, 103.75, 104.0, 104.25, 104.5, 104.75,
    105.0, 105.25, 105.5, 105.75, 106.0, 106.25, 106.5, 106.75, 107.0, 107.25,
    107.5, 107.75, 108.0, 108.25, 108.5, 108.75, 109.0, 109.25, 109.5, 109.75,
    110.0, 110.25, 110.5, 110.75, 111.0, 111.25, 111.5, 111.75, 112.0, 112.25,
    112.5, 112.75, 113.0, 113.25, 113.5, 113.75, 114.0, 114.25, 114.5, 114.75,
    115.0, 115.25, 115.5, 115.75, 116.0, 116.25, 116.5, 116.75, 117.0, 117.25,
    117.5, 117.75, 118.0, 118.25, 118.5, 118.75, 119.0, 119.25, 119.5, 119.75,
    120.0, 120.25, 120.5, 120.75, 121.0, 121.25, 121.5, 121.75, 122.0, 122.25,
    122.5, 122.75, 123.0, 123.25, 123.5, 123.75, 124.0, 124.25, 124.5, 124.75,
    125.0, 125.25, 125.5, 125.75, 126.0, 126.25, 126.5, 126.75, 127.0, 127.25,
    127.5, 127.75, 128.0, 128.25, 128.5, 128.75, 129.0, 129.25, 129.5, 129.75,
    130.0, 130.25, 130.5, 130.75, 131.0, 131.25, 131.5, 131.75, 132.0, 132.25,
    132.5, 132.75, 133.0, 133.25, 133.5, 133.75, 134.0, 134.25, 134.5, 134.75,
    135.0, 135.25, 135.5, 135.75, 136.0, 136.25, 136.5, 136.75, 137.0, 137.25,
    137.5, 137.75, 138.0, 138.25, 138.5, 138.75, 139.0, 139.25, 139.5, 139.75,
    140.0, 140.25, 140.5, 140.75, 141.0, 141.25, 141.5, 141.75, 142.0, 142.25,
    142.5, 142.75, 143.0, 143.25, 143.5, 143.75, 144.0, 144.25, 144.5, 144.75,
    145.0, 145.25, 145.5, 145.75, 146.0, 146.25, 146.5, 146.75, 147.0, 147.25,
    147.5, 147.75, 148.0, 148.25, 148.5, 148.75, 149.0, 149.25, 149.5, 149.75,
    150.0, 150.25, 150.5, 150.75, 151.0, 151.25, 151.5, 151.75, 152.0, 152.25,
    152.5, 152.75, 153.0, 153.25, 153.5, 153.75, 154.0, 154.25, 154.5, 154.75,
    155.0, 155.25, 155.5, 155.75, 156.0, 156.25, 156.5, 156.75, 157.0, 157.25,
    157.5, 157.75, 158.0, 158.25, 158.5, 158.75, 159.0, 159.25, 159.5, 159.75,
    160.0, 160.25, 160.5, 160.75, 161.0, 161.25, 161.5, 161.75, 162.0, 162.25,
    162.5, 162.75, 163.0, 163.25, 163.5, 163.75, 164.0, 164.25, 164.5, 164.75,
    165.0, 165.25, 165.5, 165.75, 166.0, 166.25, 166.5, 166.75, 167.0, 167.25,
    167.5, 167.75, 168.0, 168.25, 168.5, 168.75, 169.0, 169.25, 169.5, 169.75,
    170.0, 170.25, 170.5, 170.75, 171.0, 171.25, 171.5, 171.75, 172.0, 172.25,
    172.5, 172.75, 173.0, 173.25, 173.5, 173.75, 174.0, 174.25, 174.5, 174.75,
    175.0, 175.25, 175.5, 175.75, 176.0, 176.25, 176.5, 176.75, 177.0, 177.25,
    177.5, 177.75, 178.0, 178.25, 178.5, 178.75, 179.0, 179.25, 179.5, 179.75,
    180.0, 180.25, 180.5, 180.75, 181.0, 181.25, 181.5, 181.75, 182.0, 182.25,
    182.5, 182.75, 183.0, 183.25, 183.5, 183.75, 184.0, 184.25, 184.5, 184.75,
    185.0, 185.25, 185.5, 185.75, 186.0, 186.25, 186.5, 186.75, 187.0, 187.25,
    187.5, 187.75, 188.0, 188.25, 188.5, 188.75, 189.0, 189.25, 189.5, 189.75,
    190.0, 190.25, 190.5, 190.75, 191.0, 191.25, 191.5, 191.75, 192.0, 192.25,
    192.5, 192.75, 193.0, 193.25, 193.5, 193.75, 194.0, 194.25, 194.5, 194.75,
    195.0, 195.25, 195.5, 195.75, 196.0, 196.25, 196.5, 196.75, 197.0, 197.25,
    197.5, 197.75, 198.0, 198.25, 198.5, 198.75, 199.0, 199.25, 199.5, 199.75,
    200.0 },

  /* Variable: y
   * Referenced by: '<S5>/FWS Controller'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.00061304553327112821, 0.0024517621504112075, 0.0055148896868995045,
    0.0098003288281343082, 0.015305142548196926, 0.022025558122744426,
    0.029956969714644011, 0.039093941530586873, 0.049430211546505953,
    0.060958695799256279, 0.073671493241610708, 0.087559891157239722,
    0.10261437113197558, 0.1188246155772541, 0.13617951480127855,
    0.15466717462304458, 0.17427492452401738, 0.19498932633187582,
    0.21679618343036045, 0.23968055048892567, 0.26362674370552597,
    0.28861835155550897, 0.3146382460392485, 0.34166859442083086,
    0.36969087144971469, 0.398685872057011, 0.42863372451768056,
    0.45951390406962433, 0.49130524698032496, 0.52398596505141182,
    0.55753366055121067, 0.591925341565025, 0.62713743775264019,
    0.66314581650226989, 0.6999257994698358, 0.73745217949226927,
    0.77569923786324479, 0.81464076195950075, 0.85425006320565378,
    0.894499995365219, 0.93536297314528938, 0.976810991102108,
    1.0188156428346056, 1.061348140452731, 1.1043793343072306,
    1.1478797329673585, 1.1918195234328328, 1.2361685915661824,
    1.2808965427314614, 1.3259727226252214, 1.3713662382854415,
    1.4170459792640149, 1.4629806389482987, 1.5091387360171065,
    1.5554886360164291, 1.6019985730401076, 1.6486366715006002,
    1.6953709679749234, 1.7421694331107753, 1.7889999935778591,
    1.835830554049346, 1.8826290191983992, 1.9293633156947128,
    1.9760014141859732, 2.0225113512491739, 2.0688612512967461,
    2.1150193484224982, 2.1609540081723835, 2.2066337492251669,
    2.2520272649681554, 2.2971034449531893, 2.3418313962181814,
    2.3861804644596143, 2.4301202550414716, 2.4736206538262007,
    2.5166518478134341, 2.5591843455723349, 2.6011889974535549,
    2.6426370155669376, 2.6834999935113069, 2.7237499258427968,
    2.7633592272683796, 2.8023007515514462, 2.8405478101164885,
    2.8780741903401106, 2.9148541735158484, 2.9508625524804915,
    2.9860746488898147, 3.0204663301318764, 3.0540140258663095,
    3.0866947441782564, 3.1184860873358744, 3.1493662671406257,
    3.1793141198598218, 3.2083091207311822, 3.2363313980294897,
    3.2633617466856704, 3.2893816414489949, 3.3143732495833556,
    3.3383194430889338, 3.3612038104408786, 3.3830106678369423,
    3.4037250699463768, 3.4233328201527158, 3.4418204802834271,
    3.4591753798197669, 3.4753856245805155, 3.4904401048736586,
    3.5043285031104165, 3.5170413008763988, 3.5285697854550566,
    3.5389060557989369, 3.5480430279446713, 3.5559744398679674,
    3.5626948557752876, 3.5681996698292728, 3.5724851093053496,
    3.575548237177371, 3.577386954130505, 3.578, 3.578, 3.578, 3.578, 3.578,
    3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578,
    3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578,
    3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578,
    3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578,
    3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578,
    3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578,
    3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578,
    3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578,
    3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.578, 3.5746237326686008,
    3.5704065219289629, 3.5645082546081084, 3.5569347515837735,
    3.5476934869890773, 3.536793580836457, 3.524245790017317, 3.5100624976862731,
    3.4942577010404761, 3.4768469975060707, 3.4578475693454123,
    3.4372781667002532, 3.4151590890876169, 3.3915121653666227,
    3.3663607321960405, 3.3397296110038215, 3.3116450834913507,
    3.2821348656965781, 3.2512280806416327, 3.2189552295919248,
    3.1853481619550679, 3.1504400438493669, 3.1142653253728594, 3.07685970660523,
    3.0382601023761429, 2.9985046058347575, 2.9576324508563951,
    2.915683973323433, 2.8727005713186715, 2.8287246642704083,
    2.7837996510895944, 2.7379698673403463, 2.6912805414860954,
    2.6437777502545585, 2.5955083731655662, 2.54652004626664, 2.4968611151219604,
    2.4465805871011437, 2.3957280830148671, 2.3443537881451317,
    2.2925084027184455, 2.2402430918708189, 2.1876094351539566,
    2.1346593756324626, 2.0814451686223108, 2.0280193301211513,
    1.9744345849813751, 1.9207438148770304, 1.8670000061160081,
    1.8132561973489503, 1.7595654272265042, 1.705980682056579,
    1.6525548435132558, 1.5993406364489606, 1.5463905768613975,
    1.493756920066611, 1.4414916091292802, 1.3896462236011966,
    1.3382719286184721, 1.2874194244077255, 1.237138896251081,
    1.1874799649593526, 1.138491637902296, 1.0902222606442453,
    1.0427194692328996, 0.99603014318826222, 0.95020035923823642,
    0.90527534584645275, 0.86129943857723579, 0.818316036341754,
    0.77636755856853579, 0.73549540334061347, 0.69573990654061368,
    0.65714030204411167, 0.61973468300053247, 0.58355996423981116,
    0.54865184584191307, 0.5150447779051649, 0.48277192654816681,
    0.45186514117883714, 0.42235492306289324, 0.39427039522278318,
    0.36763927369678029, 0.34248784018659795, 0.318840916120524,
    0.29672183815766784, 0.276152435157495, 0.2571530066373785,
    0.23974230273942573, 0.22393750572635249, 0.20975421302466354,
    0.19720642183187484, 0.18630651530297446, 0.17706525032973572,
    0.16949174692497002, 0.16359347922217213, 0.15937626809945402,
    0.15684427543504106, 0.156, 0.15684427499225437, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156,
    0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156, 0.156 },

  /* Variable: zu1_0
   * Referenced by: '<S4>/Integrator7'
   */
  -0.025542532638346844,

  /* Variable: zu2_0
   * Referenced by: '<S4>/Integrator9'
   */
  -0.025542532638346844,

  /* Variable: zu3_0
   * Referenced by: '<S4>/Integrator11'
   */
  -0.01841190811165316,

  /* Variable: zu4_0
   * Referenced by: '<S4>/Integrator13'
   */
  -0.01841190811165316,

  /* Start of '<Root>/VandD' */
  {
    /* Mask Parameter: PIDController1_D
     * Referenced by: '<S90>/Derivative Gain'
     */
    9.58388407210146,

    /* Mask Parameter: PIDController1_I
     * Referenced by: '<S93>/Integral Gain'
     */
    316.552153422139,

    /* Mask Parameter: PIDController_InitialConditionForFilter
     * Referenced by: '<S43>/Filter'
     */
    0.0,

    /* Mask Parameter: PIDController1_InitialConditionForFilter
     * Referenced by: '<S91>/Filter'
     */
    0.0,

    /* Mask Parameter: PIDController_InitialConditionForIntegrator
     * Referenced by: '<S48>/Integrator'
     */
    0.0,

    /* Mask Parameter: PIDController1_InitialConditionForIntegrator
     * Referenced by: '<S96>/Integrator'
     */
    0.0,

    /* Mask Parameter: PIDController1_N
     * Referenced by: '<S99>/Filter Coefficient'
     */
    238.770606857548,

    /* Mask Parameter: PIDController1_P
     * Referenced by: '<S101>/Proportional Gain'
     */
    830.015869772746,

    /* Computed Parameter: Out1_Y0
     * Referenced by: '<S10>/Out1'
     */
    0.0,

    /* Computed Parameter: Out1_Y0_j
     * Referenced by: '<S11>/Out1'
     */
    0.0,

    /* Computed Parameter: Out1_Y0_c
     * Referenced by: '<S12>/Out1'
     */
    0.0,

    /* Computed Parameter: Out1_Y0_b
     * Referenced by: '<S13>/Out1'
     */
    0.0,

    /* Expression: 0
     * Referenced by: '<S4>/Integrator20'
     */
    0.0,

    /* Expression: 0
     * Referenced by: '<S4>/Integrator21'
     */
    0.0,

    /* Expression: 0
     * Referenced by: '<S4>/Integrator22'
     */
    0.0,

    /* Expression: [hzpc,hzpc,hzpc,hzpc]
     * Referenced by: '<S4>/[hzpc,hzpc,hzpc,hzpc]'
     */
    { 0.05, 0.05, 0.05, 0.05 },

    /* Expression: [+l1, +l1, -l2, -l2]
     * Referenced by: '<S4>/[l1,l1,-l2,-l2]'
     */
    { 1.12695753, 1.12695753, -1.59304247, -1.59304247 },

    /* Expression: [hzrc,hzrc,hzrc,hzrc]
     * Referenced by: '<S4>/hzrc'
     */
    { 0.05, 0.05, 0.05, 0.05 },

    /* Expression: [b1, -b2, b3, -b4]
     * Referenced by: '<S4>/b1, -b2, b3, -b4'
     */
    { 0.776, -0.776, 0.78, -0.78 },

    /* Expression: [+b1, -b2, +b3, -b4]
     * Referenced by: '<S4>/[b1, -b2, b3, -b4]'
     */
    { 0.776, -0.776, 0.78, -0.78 },

    /* Expression: [l1, l1, -l2, -l2]
     * Referenced by: '<S4>/[l1, l1, -l2, -l2]'
     */
    { 1.12695753, 1.12695753, -1.59304247, -1.59304247 },

    /* Expression: [-b1, +b2, -b3, +b4]
     * Referenced by: '<S4>/[-b1, b2, -b3, b4]'
     */
    { -0.776, 0.776, -0.78, 0.78 },

    /* Expression: -0.058
     * Referenced by: '<S4>/Integrator15'
     */
    -0.058,

    /* Expression: 0
     * Referenced by: '<S4>/Integrator14'
     */
    0.0,

    /* Expression: 0
     * Referenced by: '<S4>/Integrator10'
     */
    0.0,

    /* Expression: 0
     * Referenced by: '<S4>/Integrator12'
     */
    0.0,

    /* Expression: 0
     * Referenced by: '<S4>/Integrator3'
     */
    0.0,

    /* Expression: 0
     * Referenced by: '<S4>/Integrator8'
     */
    0.0,

    /* Expression: 10
     * Referenced by: '<S6>/Constant'
     */
    10.0,

    /* Expression: 10
     * Referenced by: '<S6>/Constant1'
     */
    10.0,

    /* Expression: 0
     * Referenced by: '<S4>/Integrator19'
     */
    0.0,

    /* Expression: 0
     * Referenced by: '<S4>/Integrator17'
     */
    0.0,

    /* Expression: 0
     * Referenced by: '<S4>/Integrator16'
     */
    0.0,

    /* Expression: 0
     * Referenced by: '<S4>/Integrator18'
     */
    0.0,

    /* Computed Parameter: TransferFcn_A
     * Referenced by: '<S3>/Transfer Fcn'
     */
    -30.0,

    /* Computed Parameter: TransferFcn_C
     * Referenced by: '<S3>/Transfer Fcn'
     */
    30.0,

    /* Expression: 0
     * Referenced by: '<S113>/Memory'
     */
    0.0,

    /* Expression: 0
     * Referenced by: '<S113>/Memory1'
     */
    0.0,

    /* Expression: pi/180
     * Referenced by: '<S113>/Gain'
     */
    0.017453292519943295,

    /* Expression: 1
     * Referenced by: '<S3>/Gain'
     */
    1.0,

    /* Expression: 1
     * Referenced by: '<S8>/Constant'
     */
    1.0,

    /* Expression: 1
     * Referenced by: '<S8>/Constant1'
     */
    1.0,

    /* Expression: 1
     * Referenced by: '<S9>/Constant'
     */
    1.0,

    /* Expression: 1
     * Referenced by: '<S9>/Constant1'
     */
    1.0,

    /* Computed Parameter: ManualSwitch_CurrentSetting
     * Referenced by: '<S6>/Manual Switch'
     */
    1U
  }
  /* End of '<Root>/VandD' */
};

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
