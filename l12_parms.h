
#ifndef _L12_PARMS_H
#define _L12_PARMS_H

#include <sensorDefs.h>

#define PROGRAM    "l2gen"

/* #define NBANDS        16 */
#define NBANDSIR       8
#define NQMIN          3
#define NQMAX        500
#define FILTMAX      200
#define MAXPIX     10000
#define MAX_OFILES    10
#define MAX_IFILES  1024
#define NFLAGS        32
#define NSSTFLAGS     16
#define NGIOPFLAGS    16
#define NINPRODS       3
#define NQSSTFLAGS     5

#define SUCCESS        0
#define FATAL_ERROR    1
#define LONLAT_ERROR 110
/*
#define SEAWIFS_GAC    0
#define SEAWIFS_LAC    1
#define MODISA_SUB     2
#define MODISA_FULL    3
#define MODISB_SUB     4
#define MODISB_FULL    5
 */
#define OFF            0
#define ON             1
#define NO             0
#define YES            1

/*
#define BAD_FLT  -32767.0
#define BAD_INT    -32767
#define BAD_UINT    65535
#define BAD_BYTE     -128 
#define BAD_UBYTE     255
 */

#define MAXAERMOD     100
#define AERWHITE        0
#define AERWANG        -1
#define AERRHNIR       -2
#define AERWANGNIR     -3
#define FIXMODPAIR     -4
#define FIXMODPAIRNIR  -5
#define FIXANGSTROM    -6
#define FIXANGSTROMNIR -7
#define FIXAOT         -8
#define AERWANGSWIR    -9
#define AERMUMM        -10
#define AERRHFRNIR     -13
#define AERRHSWIR      -14
#define AERRH          -15
#define AERRHMUMM      -16
#define AERRHMSEPS     -17
#define AERRHSM     -18
#define AERRHMSEPS_lin -19
#define AERSSP         -20  //Singh et al., 2019
#define AERSS14        -21  //Singh and Shanmugam, 2014
#define AERNULL        -99

#define XCALRVS          1
#define XCALPOL          2
#define XCALOLI          4 //Sudipta added for OLI SCA based XCAL

#define DEFAULT_CHL     0
#define CHL_MIN      0.00
#define CHL_MAX     100.0
#define AOT_MIN      0.00
#define AOT_MAX       1.0
#define BT_LO       -1000
#define BT_HI        1000
#define GLINT_MIN  0.0001

#define DEM_WIDTH  43200
#define DEM_HEIGHT 21600

#define BANDW 10

#define NOBRDF    0     /* brdf  */
#define FRESNSEN  1     /* bit 1 */
#define FRESNSOL  2     /* bit 2 */
#define FOQMOREL  4     /* bit 3 */
#define DTBRDF    8     /* bit 4 */
#define QMOREL   16     /* bit 5 */

#define O3_BIT     1
#define CO2_BIT    2
#define NO2_BIT    4
#define H2O_BIT    8
#define ATREM_BIT 16
#define GAS_TRANS_TBL_BIT 32

#define IOPNONE   0
#define IOPCARDER 1
#define IOPGSM    2
#define IOPQAA    3
#define IOPPML    4
#define IOPNIWA   5
#define IOPLAS    6
#define IOPGIOP   7
#define IOPSWIM   8
#define IOPDEFAULT IOPQAA

#define QAABLEND  0
#define QAA555    1
#define QAA640    2

#define NOMATCH_ERROR  110
#define FILESIZE_ERROR 111

#ifndef PI
#define PI  3.141592654
#endif
#define RADEG 57.29577951
#define STDPR 1013.25

#define STDPROC          0   /* evalmask bit definitions                 */
#define OLDAERMOD        1   /* init to old aerosol models               */
#define MODCLOUD         2   /* enables MODIS/MERIS cloud mask algorithm */
#define MODCIRRUS       16   /* enables MODIS cirrus mask                */
#define NEWSENSINFO     32   /* use test sensor info file                */
#define NEWRAYTAB       64   /* use test rayleigh tables                 */
#define NEWAERTAB      128   /* use test aerosol tables                  */
#define NEWPOLTAB      256   /* use test polarization tables             */
#define MSKMODMIR1    1024   /* mask modis mirror-side 1 (navfail)       */
#define MSKMODMIR2    2048   /* mask modis mirror-side 2 (navfail)       */
#define SSTMODS       4096   /* reserved for testing SST changes         */
#define ALTSENSORINFO 8192   /* use .alt sensor infor file in eval       */
#define TRANSSPHER   32768   /* enables spherical path geom for dtran    */

#define SOLZNIGHT      85.0
#define SOLZNIGHTA     80.0
#define DAYSCENE       0
#define NIGHTSCENE     1
#define DAYANDNIGHT    2
#define UNKNOWNSCENE   3

#define SWN  0
#define SWA  1
#define SWBB 2

#endif
