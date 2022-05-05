/* ======================================================================================== */
/* module aerosol.c  - functions to facilitate aerosol model selection and application      */
/*                                                                                          */
/* Description:                                                                             */
/*                                                                                          */
/* This code replaces the original set of fortran subroutines developed by M.Wang, H.Gordon,*/
/* and others (e.g., rho_a_sub_quad, linear_abc, funct_eps, load_aer, load_ss11) as well as */
/* the original get_aerosol() developed for MSl12.                                          */
/*                                                                                          */
/* The functions herein read and interpolate the aerosol model tables, which are now stored */
/* as individual HDF files per model.  Whfere sensor wavelengths differ from tabulated model */
/* wavelengths, interpolation is performed. Efficiencies are gained by remembering and      */
/* reusing computational results when applicable.                                           */
/*                                                                                          */
/* Primary Function:                                                                        */
/* ----------------                                                                         */
/* aerosol() - compute aerosol reflectance using specified algorithm (main interface)       */
/*                                                                                          */
/* Secondary Functions:                                                                     */
/* -------------------                                                                      */
/* wangaer() - compute aerosol reflectance using Gordon & Wang 1994 algorithm               */
/* fixedaot() - compute aerosol reflectance for fixed aot(lambda)                           */
/* fixedaer() - compute aerosol reflectance for fixed aerosol model                         */
/* get_angstrom() - compute angstrom coefficient (interface to l2_hdf_generic)              */
/* diff_tran() - compute Rayleigh-aerosol diffuse trans for selected model pair, both paths */
/*                                                                                          */
/* Supporting Functions:                                                                    */
/* --------------------                                                                     */
/* load_aermod() - loads the entire aerosol model table for the specified model list        */
/* ss_to_ms_coef() - for one model, return coefficients of function relating single         */
/*                   scattering to multiple scattering at the input geometry.               */
/* rhoas_to_rhoa() - SS aerosol reflectance to MS aerosol reflectance                       */
/* rhoa_to_rhoas() - MS aerosol reflectance to SS aerosol reflectance                       */
/* model_select_wang() - M. Wang aerosol model selection process                 .          */
/* model_select_angst() - Select model pair based on input Angstrom coefficient             */
/* model_phase() - return phase function at model wavelengths at the input geometry.        */
/* model_epsilon() - return model epsilon at input wavelengths for input geometry.          */
/* model_transmittance() - compute path Rayleigh-aerosol diffuse trans for specified model  */
/* model_taua() - compute AOT at input wavelengths using specified model                    */
/* aeroob - out-of-band water-vapor correction                                              */
/*                                                                                          */
/*                                                                                          */
/* Written By:  B. Franz, SAIC, NASA/GSFC Ocean Biology Processing Group, Summer 2004.      */
/* W. Robinson, SAIC, 24 MAr 2017, modularized and enhanced for band-dependent              */
/*  geometry                                                                                */
/*                                                                                          */
/* ======================================================================================== */

#include <float.h>

#include "l12_proto.h"

#define MAXMODEL MAXAERMOD
#define MAXSOLZ  33
#define MAXSENZ  35
#define MAXPHI   19
#define MAXSCATT 75
#define DTNTHETA 33

static float pi = PI;
static double radeg = RADEG;
static float p0 = STDPR;

static int have_ms = 0;
static int have_rh = 0;
static int have_sd = 0;
static int use_rh = 0;

static int32_t Nbands;
static int32_t Maxband; /* must be >= NBANDS */

typedef struct aermod_struct {
    char name[32];
    float rh;
    int sd;

    /* angstrom exponent (nbands+1)*/
    float *angstrom;

    /* single-scattering albedo(nbands+1), extinction coefficient(nbands+1), phase function */
    float *albedo;
    float *extc;
    float **phase;

    /* quadratic coefficients for SS to MS relationship */
    float *acost;
    float *bcost;
    float *ccost;

    /* cubic coefficients for ms_epsilon atmospheric correction ..ZA */
    float *ams_all;
    float *bms_all;
    float *cms_all;
    float *dms_all;
    float *ems_all;


    /* Rayleigh-aerosol diffuse transmittance coeffs */
    float **dtran_a;
    float **dtran_b;

    /* derived quantities */
    float **lnphase;
    float **d2phase;

} aermodstr;

aermodstr* alloc_aermodstr(int nbands, int nscatt, int nphi, int nsolz, int nsenz, int ntheta) {
    aermodstr *model;

    model = (aermodstr *) malloc(sizeof (aermodstr));
    model->angstrom = (float*) malloc(nbands * sizeof (float));
    model->albedo = (float*) malloc(nbands * sizeof (float));
    model->extc = (float*) malloc(nbands * sizeof (float));
    model->phase = allocate2d_float(nbands, nscatt);
    model->lnphase = allocate2d_float(nbands, nscatt);
    model->d2phase = allocate2d_float(nbands, nscatt);
    model->dtran_a = allocate2d_float(nbands, ntheta);
    model->dtran_b = allocate2d_float(nbands, ntheta);
    model->acost = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
    model->bcost = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
    model->ccost = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
    model->ams_all = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
    model->bms_all = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
    model->cms_all = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
    model->dms_all = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
    model->ems_all = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));

    return model;
}

typedef struct aermodtab_struct {
    int32_t sensorID;

    /* table dimensions */
    int32_t nwave;
    int32_t nmodel;
    int32_t nsolz;
    int32_t nsenz;
    int32_t nphi;
    int32_t nscatt;
    int32_t dtran_nwave;
    int32_t dtran_ntheta;

    /* table spectral bands and angles */
    float *wave;
    float *solz;
    float *senz;
    float *phi;
    float *scatt;

    /* diffuse transmittance spectral bands and angles */
    float *dtran_wave;
    float *dtran_theta;
    float *dtran_airmass;

    aermodstr **model;

} aermodtabstr;

typedef struct alphaT_struct {
    int32_t modnum;
    float angstrom;
} alphaTstr;

typedef struct epsilonT_struct {
    int32_t modnum;
    float eps_obs;
} epsilonTstr;

/* aerosol table */
static aermodtabstr *aertab = NULL;

/* structure for carrying the geometry information */
typedef struct geom_strdef {
    int gmult; /* band offset multiplier: 0 for nominal geometry
                 1 for band-dependent geometry  */
    float *senz;
    float *solz;
    float *phi;
    float *csolz; /* cosine of solz */
    float *csenz; /* cosine of senz */
    float *airmass;
    float *airmass_plp;
    float *airmass_sph;
} geom_str;

geom_str geom;

/* global variable declarations */
static int loaded = 0;
static int interpol = 0;
static int32_t *iwatab;
static int32_t *iwdtab;

static int32_t iwnir_s = -1;
static int32_t iwnir_l = -1;

static int imm50 = -1;
static int imc50 = -1;
static int imc70 = -1;
static int imt90 = -1;
static int imt99 = -1;
static int wang_modx = 0;

static float mu0;
static float mu;
static float airmass;

static int32_t evalmask = 0;
static int32_t aer_opt = 0;
static float airmass_plp;
static float airmass_sph;

int cmpfunc(const void * a, const void * b) {
    if (*(double*) a > *(double*) b) return 1;
    else if (*(double*) a < *(double*) b) return -1;
    else return 0;
}



/* ---------------------------------------------------------------------------------------- */
/* first_deriv() - returns first derivative (dy/dx) of 1st or last array indices using a    */
/*                 4-pt Lagrangian interpolation.  NOTE: It is assumed that 4 points exist. */

/* ---------------------------------------------------------------------------------------- */
float first_deriv(float x[], float y[], int n) {
    float a1, a2, a3, a4, a5, a6, d1;

    if (n == 0) {

        a1 = x[0] - x[1];
        a2 = x[0] - x[2];
        a3 = x[0] - x[3];
        a4 = x[1] - x[2];
        a5 = x[1] - x[3];
        a6 = x[2] - x[3];

        d1 = y[0]*(1.0 / a1 + 1.0 / a2 + 1.0 / a3)
                - a2 * a3 * y[1] / (a1 * a4 * a5)
                + a1 * a3 * y[2] / (a2 * a4 * a6)
                - a1 * a2 * y[3] / (a3 * a5 * a6);

    } else {

        a1 = x[n - 1] - x[n - 4];
        a2 = x[n - 1] - x[n - 3];
        a3 = x[n - 1] - x[n - 2];
        a4 = x[n - 2] - x[n - 4];
        a5 = x[n - 2] - x[n - 3];
        a6 = x[n - 3] - x[n - 4];

        d1 = -a2 * a3 * y[n - 4] / (a6 * a4 * a1)
                + a1 * a3 * y[n - 3] / (a6 * a5 * a2)
                - a1 * a2 * y[n - 2] / (a4 * a5 * a3)
                + y[n - 1]*(1.0 / a1 + 1.0 / a2 + 1.0 / a3);
    }

    return (d1);
}


/* ---------------------------------------------------------------------------------------- */
/* load_aermod() - loads the entire aerosol model table for the specified model list        */

/* ---------------------------------------------------------------------------------------- */
int load_aermod(int32_t sensorID, float wave[], int32_t nwave, char *aermodfile, char models[MAXAERMOD][32], int32_t nmodels) {
    int32 sd_id;
    int32 sds_id;
    int32 status;
    int32 rank;
    int32 nt;
    int32 dims[H4_MAX_VAR_DIMS];
    int32 nattrs;
    int32 start[4] = {0, 0, 0, 0};

    int32_t mwave, msolz, msenz, mphi, mscatt, dtran_mwave, dtran_mtheta;

    float d1phase1;
    float d1phaseN;
    float rh;
    int16 sd;

    char name [H4_MAX_NC_NAME] = "";
    char sdsname[H4_MAX_NC_NAME] = "";
    char file [FILENAME_MAX] = "";
    char path [FILENAME_MAX] = "";

    int iw, im, is, iwbase, i;
    static int firstCall = 1;

    if (firstCall == 1) {
        if ((iwatab = (int32_t *) calloc(nwave, sizeof (int32_t))) == NULL) {
            printf("Unable to allocate space for iwatab.\n");
            exit(1);
        }
        if ((iwdtab = (int32_t *) calloc(nwave, sizeof (int32_t))) == NULL) {
            printf("Unable to allocate space for iwdtab.\n");
            exit(1);
        }
        firstCall = 0;
    }

    printf("Loading aerosol models from %s\n", aermodfile);

    for (im = 0; im < nmodels + 1; im++) {

        if (im < nmodels) {

            strcpy(file, path);
            strcat(file, aermodfile);
            strcat(file, "_");
            strcat(file, models[im]);
            strcat(file, ".hdf");

            /* Open the file and initiate the SD interface */
            sd_id = SDstart(file, DFACC_RDONLY);
            if (sd_id == -1) {
                printf("-E- %s:  Error opening file %s.\n",
                        __FILE__, file);
                exit(1);
            }

        } else {

            strcpy(file, path);
            strcat(file, aermodfile);
            strcat(file, "_default.hdf");

            /* Open the file and initiate the SD interface */
            sd_id = SDstart(file, DFACC_RDONLY);
            if (sd_id == -1) {
                printf("-E- %s:  Error opening file %s.\n",
                        __FILE__, file);
                exit(1);
            }
        }

        /* read attributes which should be constant between models */

        status = SDreadattr(sd_id, SDfindattr(sd_id, "Number of Wavelengths"),
                &mwave);
        status = SDreadattr(sd_id, SDfindattr(sd_id,
                "Number of Solar Zenith Angles"), &msolz);
        status = SDreadattr(sd_id, SDfindattr(sd_id,
                "Number of View Zenith Angles"), &msenz);
        status = SDreadattr(sd_id, SDfindattr(sd_id,
                "Number of Relative Azimuth Angles"), &mphi);
        status = SDreadattr(sd_id, SDfindattr(sd_id,
                "Number of Scattering Angles"), &mscatt);
        status = SDreadattr(sd_id, SDfindattr(sd_id,
                "Number of Diffuse Transmittance Wavelengths"), &dtran_mwave);
        status = SDreadattr(sd_id, SDfindattr(sd_id,
                "Number of Diffuse Transmittance Zenith Angles"), &dtran_mtheta);

        if (im == 0) {

            int32_t nwave = mwave;
            int32_t nsolz = msolz;
            int32_t nsenz = msenz;
            int32_t nphi = mphi;
            int32_t nscatt = mscatt;
            int32_t dtran_nwave = dtran_mwave;
            int32_t dtran_ntheta = dtran_mtheta;

            printf("Number of Wavelengths                          %d\n", nwave);
            printf("Number of Solar Zenith Angles                  %d\n", nsolz);
            printf("Number of View Zenith Angles                   %d\n", nsenz);
            printf("Number of Relative Azimuth Angles              %d\n", nphi);
            printf("Number of Scattering Angles                    %d\n", nscatt);

            printf("Number of Diffuse Transmittance Wavelengths    %d\n", dtran_nwave);
            printf("Number of Diffuse Transmittance Zenith Angles  %d\n", dtran_ntheta);

            // allocate the aerosol table

            if ((aertab = (aermodtabstr *) calloc(1, sizeof (aermodtabstr))) == NULL) {
                printf("Unable to allocate space for aerosol table.\n");
                exit(1);
            }

            aertab->nmodel = nmodels;
            aertab->nwave = nwave;
            aertab->nsolz = nsolz;
            aertab->nsenz = nsenz;
            aertab->nphi = nphi;
            aertab->nscatt = nscatt;

            aertab->wave = (float *) malloc(nwave * sizeof (float));
            aertab->solz = (float *) malloc(nsolz * sizeof (float));
            aertab->senz = (float *) malloc(nsenz * sizeof (float));
            aertab->phi = (float *) malloc(nphi * sizeof (float));
            aertab->scatt = (float *) malloc(nscatt * sizeof (float));

            aertab->dtran_nwave = dtran_nwave;
            aertab->dtran_ntheta = dtran_ntheta;

            aertab->dtran_wave = (float *) malloc(dtran_nwave * sizeof (float));
            aertab->dtran_theta = (float *) malloc(dtran_ntheta * sizeof (float));
            aertab->dtran_airmass = (float *) malloc(dtran_ntheta * sizeof (float));

            // allocate the model tables

            if ((aertab->model = (aermodstr **) calloc(1, (nmodels + 1) * sizeof (aermodstr*))) == NULL) {
                printf("Unable to allocate space for %d aerosol models.\n", nmodels + 1);
                exit(1);
            }
            for (i = 0; i < nmodels + 1; i++) {
                if ((aertab->model[i] = alloc_aermodstr(nwave + 1, nscatt, MAXPHI, MAXSOLZ, MAXSENZ, dtran_ntheta)) == NULL) {
                    printf("Unable to allocate space for aerosol model %d.\n", im);
                    exit(1);
                }
            }

            /* read SDSes which are constant between models */

            strcpy(sdsname, "wave");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->wave);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname, "solz");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->solz);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname, "senz");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->senz);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname, "phi");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->phi);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname, "scatt");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->scatt);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname, "dtran_wave");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->dtran_wave);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname, "dtran_theta");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->dtran_theta);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

        } else {
            /*  check that all the aerosol files at least have the same
                main dimensions  */
            if ((aertab->nwave != mwave) || (aertab->nsolz != msolz) ||
                    (aertab->nsenz != msenz) || (aertab->nphi != mphi) ||
                    (aertab->nscatt != mscatt) ||
                    (aertab->dtran_nwave != dtran_mwave) ||
                    (aertab->dtran_ntheta != dtran_mtheta)) {
                printf("-E- %s, %d:  Error, Aerosol table %s\n",
                        __FILE__, __LINE__, file);
                printf("    has different dimensions from previous tables\n");
                exit(1);
            }
        }

        if (im < nmodels)
            strncpy(aertab->model[im]->name, models[im], 32);
        else
            strncpy(aertab->model[im]->name, "default", 32);

        status = SDreadattr(sd_id, SDfindattr(sd_id, "Relative Humidity"), &rh);
        if (status == 0) {
            aertab->model[im]->rh = rh;
            have_rh = 1;
        } else
            aertab->model[im]->rh = -1.0;

        status = SDreadattr(sd_id, SDfindattr(sd_id, "Size Distribution"), &sd);
        if (status == 0) {
            aertab->model[im]->sd = sd;
            have_sd = 1;
        } else
            aertab->model[im]->sd = -1;

        strcpy(sdsname, "albedo");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->albedo);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "extc");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->extc);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        /* now computed from model parameters
        strcpy(sdsname,"angstrom");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->angstrom);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }
         */

        strcpy(sdsname, "phase");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->phase[0]);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "acost");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->acost);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "bcost");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->bcost);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "ccost");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->ccost);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        // check for multi-scattering epsilon tables
        if (SDnametoindex(sd_id, "ams_all") != FAIL) {

            have_ms = 1;

            strcpy(sdsname, "ams_all");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->ams_all);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname, "bms_all");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->bms_all);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname, "cms_all");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->cms_all);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname, "dms_all");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->dms_all);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname, "ems_all");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->ems_all);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

        }

        strcpy(sdsname, "dtran_a");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->dtran_a[0]);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "dtran_b");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->dtran_b[0]);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        /* terminate access to the SD interface and close the file */
        status = SDend(sd_id);


        /* compute angstrom exponent for each model wavelength relative to max wavelength */
        iwbase = windex(865, aertab->wave, aertab->nwave);
        for (iw = 0; iw < aertab->nwave; iw++) {
            if (iw != iwbase)
                aertab->model[im]->angstrom[iw] = -log(aertab->model[im]->extc[iw] / aertab->model[im]->extc[iwbase]) /
                log(aertab->wave[iw] / aertab->wave[iwbase]);
        }
        aertab->model[im]->angstrom[iwbase] = aertab->model[im]->angstrom[iwbase - 1];

        /* precompute log of phase function and 2nd derivative (for cubic interp) */
        for (iw = 0; iw < aertab->nwave; iw++) {
            for (is = 0; is < aertab->nscatt; is++) {
                aertab->model[im]->lnphase[iw][is] = log(aertab->model[im]->phase[iw][is]);
            }
            d1phase1 = first_deriv(aertab->scatt, &aertab->model[im]->lnphase[iw][0], 0);
            d1phaseN = first_deriv(aertab->scatt, &aertab->model[im]->lnphase[iw][0], aertab->nscatt);
            spline(aertab->scatt,
                    &aertab->model[im]->lnphase[iw][0],
                    aertab->nscatt,
                    d1phase1,
                    d1phaseN,
                    &aertab->model[im]->d2phase[iw][0]);
        }

        /* precompute airmass for diffuse transmittance */
        for (iw = 0; iw < aertab->dtran_nwave; iw++) {
            for (is = 0; is < aertab->dtran_ntheta; is++) {
                aertab->dtran_airmass[is] = 1.0 / cos(aertab->dtran_theta[is] / radeg);
            }
        }
    }

    aertab->nmodel = nmodels;
    aertab->sensorID = sensorID;

    /* map sensor wavelengths to table wavelengths */
    for (iw = 0; iw < nwave; iw++) {
        iwatab[iw] = windex(wave[iw], aertab->wave, aertab->nwave);
        iwdtab[iw] = windex(wave[iw], aertab->dtran_wave, aertab->dtran_nwave);
        if (fabsf(wave[iw] - aertab->wave[iwatab[iw]]) > 0.5) {
            printf("Aerosol model coefficients will be interpolated for %5.1f nm band.\n", wave[iw]);
            interpol = 1;
        }
    }

    /* in case of more aertab wavelengths then sensor wavelengths */
    if (aertab->nwave != nwave)
        interpol = 1;

    loaded = 1;

    /* set-up for Wang model cross-over correction */
    for (im = 0; im < nmodels; im++) {
        if (strcmp(models[im], "m50") == 0) imm50 = im;
        else if (strcmp(models[im], "c50") == 0) imc50 = im;
        else if (strcmp(models[im], "c70") == 0) imc70 = im;
        else if (strcmp(models[im], "t90") == 0) imt90 = im;
        else if (strcmp(models[im], "t99") == 0) imt99 = im;
    }
    if (imm50 >= 0 && imc50 >= 0 && imc70 >= 0 && imt90 >= 0 && imt99 >= 0) {
        wang_modx = 1;
        printf("\nM. Wang model cross-over correction enabled.\n");
    }

    return (0);
}



#define INDEX(iw,isol,iphi,isen) (iw*aertab->nsolz*aertab->nphi*aertab->nsenz + isol*aertab->nphi*aertab->nsenz + iphi*aertab->nsenz + isen)

/* ---------------------------------------------------------------------------------------- */
/* ss_to_ms_coef() - for one model, return coefficients of function relating single         */
/*                   scattering to multiple scattering at the input geometry.               */
/*                                                                                          */
/* This is effectively a C version of M. Wangs linear_a_b_c.f.  The program optimizes for   */
/* multiple calls at the same geometry by computing for all models on the first call with a */
/* new geometry.  It returns pointers to the internal static arrays of coefficients.        */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */
/* W. Robinson, SAIC  adapt to band-ependent geometry                                       */

/* ---------------------------------------------------------------------------------------- */
void ss_to_ms_coef(int modnum, geom_str *geom, float **a, float **b, float **c) {
    static float lastsolz = -999.;
    static float lastsenz = -999.;
    static float lastphi = -999.;

    static int computed[MAXMODEL];

    static float *a_coef[MAXMODEL];
    static float *b_coef[MAXMODEL];
    static float *c_coef[MAXMODEL];

    static float *p, *q, *r;
    static float as000, as100, as010, as110, as001, as011, as101, as111;
    static float ai000, ai100, ai010, ai110, ai001, ai011, ai101, ai111;
    static float ac000, ac100, ac010, ac110, ac001, ac011, ac101, ac111;

    static int *isolz1, *isolz2;
    static int *isenz1, *isenz2;
    static int *iphi1, *iphi2;

    static float *p_ar, *q_ar, *r_ar;
    static float p_cnst, q_cnst, r_cnst;
    static int *isolz1_ar, *isolz2_ar, isolz1_cnst, isolz2_cnst;
    static int *isenz1_ar, *isenz2_ar, isenz1_cnst, isenz2_cnst;
    static int *iphi1_ar, *iphi2_ar, iphi1_cnst, iphi2_cnst;

    float aphi;
    float px, qx, rx;
    int im, iw, i, ig;
    static int firstCall = 1;
    static int gmult;

    if (firstCall == 1) {
        firstCall = 0;
        for (i = 0; i < MAXMODEL; i++) {
            if ((a_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for a_coef.\n");
                exit(1);
            }
            if ((b_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for rhoa.\n");
                exit(1);
            }
            if ((c_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for rhoa.\n");
                exit(1);
            }
        }
        /*  set up indicies, weights for band-dependent or nominal geometry */
        if ((geom->gmult == 0) || (interpol == 1)) {
            gmult = 0;
            p = &p_cnst;
            q = &q_cnst;
            r = &r_cnst;
            isolz1 = &isolz1_cnst;
            isolz2 = &isolz2_cnst;
            isenz1 = &isenz1_cnst;
            isenz2 = &isenz2_cnst;
            iphi1 = &iphi1_cnst;
            iphi2 = &iphi2_cnst;
        } else
            gmult = 1;
        {
            if (((p_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                    == NULL) ||
                    ((q_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                    == NULL) ||
                    ((r_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                    == NULL)) {
                printf("Unable to allocate space for p, q, r weights.\n");
                exit(1);
            }
            if (((isolz1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((isenz1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((iphi1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL)) {
                printf("Unable to allocate space for interp indicies 1.\n");
                exit(1);
            }
            if (((isolz2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((isenz2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((iphi2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL)) {
                printf("Unable to allocate space for interp indicies 2.\n");
                exit(1);
            }
            p = p_ar;
            q = q_ar;
            r = r_ar;
            isolz1 = isolz1_ar;
            isolz2 = isolz2_ar;
            isenz1 = isenz1_ar;
            isenz2 = isenz2_ar;
            iphi1 = iphi1_ar;
            iphi2 = iphi2_ar;
        }
    }

    if ((geom->solz[0] != lastsolz) || (geom->senz[0] != lastsenz) ||
            (geom->phi[0] != lastphi)) {
        for (im = 0; im < aertab->nmodel; im++)
            computed[im] = 0;

        for (iw = 0; iw < aertab->nwave; iw++) {
            ig = iw * gmult;
            /* find bracketing solar indices */
            for (i = 0; i < aertab->nsolz; i++) {
                if (geom->solz[ig] < aertab->solz[i])
                    break;
            }
            isolz1[iw] = MAX(i - 1, 0);
            isolz2[iw] = MIN(i, aertab->nsolz - 1);
            if (isolz2[iw] != isolz1[iw])
                r[iw] = (geom->solz[ig] - aertab->solz[isolz1[iw]]) /
                (aertab->solz[isolz2[iw]] - aertab->solz[isolz1[iw]]);
            else
                r[iw] = 0.0;

            /* find bracketing view indices */
            for (i = 0; i < aertab->nsenz; i++) {
                if (geom->senz[ig] < aertab->senz[i])
                    break;
            }
            isenz1[iw] = MAX(i - 1, 0);
            isenz2[iw] = MIN(i, aertab->nsenz - 1);
            if (isenz2[iw] != isenz1[iw])
                p[iw] = (geom->senz[ig] - aertab->senz[isenz1[iw]]) /
                (aertab->senz[isenz2[iw]] - aertab->senz[isenz1[iw]]);
            else
                p[iw] = 0.0;

            /* find bracketing azimuth indices */
            aphi = fabs(geom->phi[ig]);
            for (i = 0; i < aertab->nphi; i++) {
                if (aphi < aertab->phi[i])
                    break;
            }
            iphi1[iw] = MAX(i - 1, 0);
            iphi2[iw] = MIN(i, aertab->nphi - 1);
            if (iphi2[iw] != iphi1[iw])
                q[iw] = (aphi - aertab->phi[iphi1[iw]]) /
                (aertab->phi[iphi2[iw]] - aertab->phi[iphi1[iw]]);
            else
                q[iw] = 0.0;
            if (gmult == 0) break;
        }

        /* remember last geometry */
        lastsolz = geom->solz[0];
        lastsenz = geom->senz[0];
        lastphi = geom->phi[0];
    }

    if (!computed[modnum]) {
        im = modnum;
        computed[modnum] = 1;

        for (iw = 0; iw < aertab->nwave; iw++) {
            ig = iw * gmult;
            px = p[ig];
            qx = q[ig];
            rx = r[ig];
            if (isolz2[ig] == 0) {
                as000 = aertab->model[im]->acost[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                as100 = aertab->model[im]->acost[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                as001 = aertab->model[im]->acost[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                as101 = aertab->model[im]->acost[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                ai000 = aertab->model[im]->bcost[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                ai100 = aertab->model[im]->bcost[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                ai001 = aertab->model[im]->bcost[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                ai101 = aertab->model[im]->bcost[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                ac000 = aertab->model[im]->ccost[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                ac100 = aertab->model[im]->ccost[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                ac001 = aertab->model[im]->ccost[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                ac101 = aertab->model[im]->ccost[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                a_coef[im][iw] = (1. - px)*(1. - rx) * as000 + px * rx * as101
                        + (1. - px) * rx * as001 + px * (1. - rx) * as100;

                b_coef[im][iw] = (1. - px)*(1. - rx) * ai000 + px * rx * ai101
                        + (1. - px) * rx * ai001 + px * (1. - qx)*(1. - rx) * ai100;

                c_coef[im][iw] = (1. - px)*(1. - rx) * ac000 + px * rx * ac101
                        + (1. - px) * rx * ac001 + px * (1. - qx)*(1. - rx) * ac100;
            } else {
            	as000 = aertab->model[im]->acost[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
                as100 = aertab->model[im]->acost[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
                as010 = aertab->model[im]->acost[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
                as110 = aertab->model[im]->acost[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
                as001 = aertab->model[im]->acost[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
                as011 = aertab->model[im]->acost[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
                as101 = aertab->model[im]->acost[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
                as111 = aertab->model[im]->acost[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

                ai000 = aertab->model[im]->bcost[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
                ai100 = aertab->model[im]->bcost[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
                ai010 = aertab->model[im]->bcost[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
                ai110 = aertab->model[im]->bcost[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
                ai001 = aertab->model[im]->bcost[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
                ai011 = aertab->model[im]->bcost[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
                ai101 = aertab->model[im]->bcost[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
                ai111 = aertab->model[im]->bcost[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

                ac000 = aertab->model[im]->ccost[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
                ac100 = aertab->model[im]->ccost[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
                ac010 = aertab->model[im]->ccost[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
                ac110 = aertab->model[im]->ccost[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
                ac001 = aertab->model[im]->ccost[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
                ac011 = aertab->model[im]->ccost[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
                ac101 = aertab->model[im]->ccost[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
                ac111 = aertab->model[im]->ccost[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

                a_coef[im][iw] = (1. - px)*(1. - qx)*(1. - rx) * as000 + px * qx * rx * as111
                        + px * (1. - qx) * rx * as101 + (1. - px) * qx * (1. - rx) * as010
                        + px * qx * (1. - rx) * as110 + (1. - px)*(1. - qx) * rx * as001
                        + (1. - px) * qx * rx * as011 + px * (1. - qx)*(1. - rx) * as100;

                b_coef[im][iw] = (1. - px)*(1. - qx)*(1. - rx) * ai000 + px * qx * rx * ai111
                        + px * (1. - qx) * rx * ai101 + (1. - px) * qx * (1. - rx) * ai010
                        + px * qx * (1. - rx) * ai110 + (1. - px)*(1. - qx) * rx * ai001
                        + (1. - px) * qx * rx * ai011 + px * (1. - qx)*(1. - rx) * ai100;

                c_coef[im][iw] = (1. - px)*(1. - qx)*(1. - rx) * ac000 + px * qx * rx * ac111
                        + px * (1. - qx) * rx * ac101 + (1. - px) * qx * (1. - rx) * ac010
                        + px * qx * (1. - rx) * ac110 + (1. - px)*(1. - qx) * rx * ac001
                        + (1. - px) * qx * rx * ac011 + px * (1. - qx)*(1. - rx) * ac100;
            }
        }
    }

    /* return pointers to coeffs for this geometry */
    *a = &a_coef[modnum][0];
    *b = &b_coef[modnum][0];
    *c = &c_coef[modnum][0];

    return;
}

/* ---------------------------------------------------------------------------------------- */
/* fresnel_coef() - computes Fresnel reflectance coefficient for specified index of refr.   */

/* ---------------------------------------------------------------------------------------- */
float fresnel_coef(float mu, float index) {
    float sq, r2, q1;

    sq = sqrt(pow(index, 2.0) - 1.0 + pow(mu, 2.0));
    r2 = pow((mu - sq) / (mu + sq), 2.0);
    q1 = (1.0 - pow(mu, 2.0) - mu * sq) / (1.0 - pow(mu, 2.0) + mu * sq);

    return (r2 * (q1 * q1 + 1.0) / 2.0);
}



/* ---------------------------------------------------------------------------------------- */
/* ms_eps_coef() -   for a given model, returns  ams, bms, cms, dms and ems coefficients to */
/*                   compute ms_reflectance at the input geometry.                          */
/*                   Also, the program optimizes for multiple calls at the same geometry    */
/*                   by computing for all models on the first call with a new geometry.     */
/*                   It returns pointers to the internal static arrays of coefficients.     */
/*                                                                                          */
/* Z. Ahmad July 08, 2014                                                                   */
/* W. Robinson, SAIC  23 Mar 2017  add band-dependent geometry to code                      */
/* M. Zhang 06/19/2019   fixed the interpolation part                                      */

/* ---------------------------------------------------------------------------------------- */



void ms_eps_coef(int modnum, int32_t iwnir_l, float wave[], geom_str *geom,
        float **a, float **b, float **c, float **d, float **e)
 {
    static float lastsolz = -999.;
    static float lastsenz = -999.;
    static float lastphi = -999.;

    static int computed[MAXMODEL];

    static float *a_coef[MAXMODEL];
    static float *b_coef[MAXMODEL];
    static float *c_coef[MAXMODEL];
    static float *d_coef[MAXMODEL];
    static float *e_coef[MAXMODEL];

    static float *p, *q, *r;
    static float as000, as100, as010, as110, as001, as011, as101, as111;
    static float ai000, ai100, ai010, ai110, ai001, ai011, ai101, ai111;
    static float ac000, ac100, ac010, ac110, ac001, ac011, ac101, ac111;
    static float ad000, ad100, ad010, ad110, ad001, ad011, ad101, ad111;
    static float ae000, ae100, ae010, ae110, ae001, ae011, ae101, ae111;

    static int *isolz1, *isolz2;
    static int *isenz1, *isenz2;
    static int *iphi1, *iphi2;

    static float *p_ar, *q_ar, *r_ar;
    static float p_cnst, q_cnst, r_cnst;
    static int *isolz1_ar, *isolz2_ar, isolz1_cnst, isolz2_cnst;
    static int *isenz1_ar, *isenz2_ar, isenz1_cnst, isenz2_cnst;
    static int *iphi1_ar, *iphi2_ar, iphi1_cnst, iphi2_cnst;
    static int gmult;

    float aphi;
    float px, qx, rx;
    int im, iw, i, ig;
    static int firstCall = 1;

    if (firstCall == 1) {
        firstCall = 0;
        for (i = 0; i < MAXMODEL; i++) {
            if ((a_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for a_coef.\n");
                exit(1);
            }
            if ((b_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for b_coef.\n");
                exit(1);
            }
            if ((c_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for c_coef.\n");
                exit(1);
            }
            if ((d_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for d_coef.\n");
                exit(1);
            }

            if ((e_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for e_coef.\n");
                exit(1);
            }

        }
        /*  set up indicies, weights for band-dependent or nominal geometry */
        if ((geom->gmult == 0) || (interpol == 1)) {
            gmult = 0;
            p = &p_cnst;
            q = &q_cnst;
            r = &r_cnst;
            isolz1 = &isolz1_cnst;
            isolz2 = &isolz2_cnst;
            isenz1 = &isenz1_cnst;
            isenz2 = &isenz2_cnst;
            iphi1 = &iphi1_cnst;
            iphi2 = &iphi2_cnst;
        } else {
            gmult = 1;
            if (((p_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                    == NULL) ||
                    ((q_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                    == NULL) ||
                    ((r_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                    == NULL)) {
                printf("Unable to allocate space for p, q, r weights.\n");
                exit(1);
            }
            if (((isolz1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((isenz1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((iphi1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL)) {
                printf("Unable to allocate space for interp indicies 1.\n");
                exit(1);
            }
            if (((isolz2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((isenz2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((iphi2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL)) {
                printf("Unable to allocate space for interp indicies 2.\n");
                exit(1);
            }
            p = p_ar;
            q = q_ar;
            r = r_ar;
            isolz1 = isolz1_ar;
            isolz2 = isolz2_ar;
            isenz1 = isenz1_ar;
            isenz2 = isenz2_ar;
            iphi1 = iphi1_ar;
            iphi2 = iphi2_ar;
        }
    }

    if (geom->solz[0] != lastsolz || geom->senz[0] != lastsenz ||
            geom->phi[0] != lastphi) {
        for (im = 0; im < aertab->nmodel; im++)
            computed[im] = 0;

        for (iw = 0; iw < aertab->nwave; iw++) {
            ig = iw * gmult;
            /* find bracketing solar indices */
            for (i = 0; i < aertab->nsolz; i++) {
                if (geom->solz[ig] < aertab->solz[i])
                    break;
            }
            isolz1[iw] = MAX(i - 1, 0);
            isolz2[iw] = MIN(i, aertab->nsolz - 1);
            if (isolz2[iw] != isolz1[iw])
                r[iw] = (geom->solz[ig] - aertab->solz[isolz1[iw]]) /
                (aertab->solz[isolz2[iw]] - aertab->solz[isolz1[iw]]);
            else
                r[iw] = 0.0;

            /* find bracketing view indices */
            for (i = 0; i < aertab->nsenz; i++) {
                if (geom->senz[ig] < aertab->senz[i])
                    break;
            }
            isenz1[iw] = MAX(i - 1, 0);
            isenz2[iw] = MIN(i, aertab->nsenz - 1);
            if (isenz2[iw] != isenz1[iw])
                p[iw] = (geom->senz[ig] - aertab->senz[isenz1[iw]]) /
                (aertab->senz[isenz2[iw]] - aertab->senz[isenz1[iw]]);
            else
                p[iw] = 0.0;

            /* find bracketing azimuth indices */
            aphi = fabs(geom->phi[ig]);
            for (i = 0; i < aertab->nphi; i++) {
                if (aphi < aertab->phi[i])
                    break;
            }
            iphi1[iw] = MAX(i - 1, 0);
            iphi2[iw] = MIN(i, aertab->nphi - 1);
            if (iphi2[iw] != iphi1[iw])
                q[iw] = (aphi - aertab->phi[iphi1[iw]]) /
                (aertab->phi[iphi2[iw]] - aertab->phi[iphi1[iw]]);
            else
                q[iw] = 0.0;
            if (gmult == 0) break;
        }

        /* remember last geometry */
        lastsolz = geom->solz[0];
        lastsenz = geom->senz[0];
        lastphi = geom->phi[0];

    }

    im = modnum;

    if (!computed[modnum]) {
        im = modnum;
        computed[modnum] = 1;

        for (iw = 0; iw < aertab->nwave; iw++) {
            ig = iw * gmult;
            px = p[ig];
            qx = q[ig];
            rx = r[ig];
            if (isolz2[ig] == 0) {
                as000 = aertab->model[im]->ams_all[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                as100 = aertab->model[im]->ams_all[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                as001 = aertab->model[im]->ams_all[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                as101 = aertab->model[im]->ams_all[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                ai000 = aertab->model[im]->bms_all[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                ai100 = aertab->model[im]->bms_all[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                ai001 = aertab->model[im]->bms_all[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                ai101 = aertab->model[im]->bms_all[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                ac000 = aertab->model[im]->cms_all[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                ac100 = aertab->model[im]->cms_all[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                ac001 = aertab->model[im]->cms_all[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                ac101 = aertab->model[im]->cms_all[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                ad000 = aertab->model[im]->dms_all[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                ad100 = aertab->model[im]->dms_all[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                ad001 = aertab->model[im]->dms_all[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                ad101 = aertab->model[im]->dms_all[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                ae000 = aertab->model[im]->ems_all[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                ae100 = aertab->model[im]->ems_all[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                ae001 = aertab->model[im]->ems_all[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                ae101 = aertab->model[im]->ems_all[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                a_coef[im][iw] = (1. - px)*(1. - rx) * as000 + px * rx * as101
                        + (1. - px) * rx * as001 + px * (1. - rx) * as100;

                b_coef[im][iw] = (1. - px)*(1. - rx) * ai000 + px * rx * ai101
                        + (1. - px) * rx * ai001 + px * (1. - qx)*(1. - rx) * ai100;

                c_coef[im][iw] = (1. - px)*(1. - rx) * ac000 + px * rx * ac101
                        + (1. - px) * rx * ac001 + px * (1. - qx)*(1. - rx) * ac100;

                d_coef[im][iw] = (1. - px)*(1. - rx) * ad000 + px * rx * ad101
                        + (1. - px) * rx * ad001 + px * (1. - qx)*(1. - rx) * ad100;

                e_coef[im][iw] = (1. - px)*(1. - rx) * ae000 + px * rx * ae101
                        + (1. - px) * rx * ae001 + px * (1. - qx)*(1. - rx) * ae100;
            } else {
                /*       printf("coeffs: as000,ai000,ac000,ad000,ae000\n");   */

            as000 = aertab->model[im]->ams_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
            as100 = aertab->model[im]->ams_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
            as010 = aertab->model[im]->ams_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
            as110 = aertab->model[im]->ams_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
            as001 = aertab->model[im]->ams_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
            as011 = aertab->model[im]->ams_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
            as101 = aertab->model[im]->ams_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
            as111 = aertab->model[im]->ams_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

            ai000 = aertab->model[im]->bms_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
            ai100 = aertab->model[im]->bms_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
            ai010 = aertab->model[im]->bms_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
            ai110 = aertab->model[im]->bms_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
            ai001 = aertab->model[im]->bms_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
            ai011 = aertab->model[im]->bms_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
            ai101 = aertab->model[im]->bms_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
            ai111 = aertab->model[im]->bms_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

            ac000 = aertab->model[im]->cms_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
            ac100 = aertab->model[im]->cms_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
            ac010 = aertab->model[im]->cms_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
            ac110 = aertab->model[im]->cms_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
            ac001 = aertab->model[im]->cms_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
            ac011 = aertab->model[im]->cms_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
            ac101 = aertab->model[im]->cms_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
            ac111 = aertab->model[im]->cms_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

            ad000 = aertab->model[im]->dms_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
            ad100 = aertab->model[im]->dms_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
            ad010 = aertab->model[im]->dms_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
            ad110 = aertab->model[im]->dms_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
            ad001 = aertab->model[im]->dms_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
            ad011 = aertab->model[im]->dms_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
            ad101 = aertab->model[im]->dms_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
            ad111 = aertab->model[im]->dms_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

            ae000 = aertab->model[im]->ems_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
            ae100 = aertab->model[im]->ems_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
            ae010 = aertab->model[im]->ems_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
            ae110 = aertab->model[im]->ems_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
            ae001 = aertab->model[im]->ems_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
            ae011 = aertab->model[im]->ems_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
            ae101 = aertab->model[im]->ems_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
            ae111 = aertab->model[im]->ems_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];


           a_coef[im][iw] = (1. - rx)*(1. - qx)*(1. - px) * as000 + rx * qx * px * as111
                            + rx * (1. - qx) * px * as101 + (1. - rx) * qx * (1. - px) * as010
                            + rx * qx * (1. - px) * as110 + (1. - rx)*(1. - qx) * px * as001
                            + (1. - rx) * qx * px * as011 + rx * (1. - qx)*(1. - px) * as100;

           b_coef[im][iw] = (1. - rx)*(1. - qx)*(1. - px) * ai000 + rx * qx * px * ai111
                            + rx * (1. - qx) * px * ai101 + (1. - rx) * qx * (1. - px) * ai010
                            + rx * qx * (1. - px) * ai110 + (1. - rx)*(1. - qx) * px * ai001
                            + (1. - rx) * qx * px * ai011 + rx * (1. - qx)*(1. - px) * ai100;

           c_coef[im][iw] = (1. - rx)*(1. - qx)*(1. - px) * ac000 + rx * qx * px * ac111
                            + rx * (1. - qx) * px * ac101 + (1. - rx) * qx * (1. - px) * ac010
                            + rx * qx * (1. - px) * ac110 + (1. - rx)*(1. - qx) * px * ac001
                            + (1. - rx) * qx * px * ac011 + rx * (1. - qx)*(1. - px) * ac100;

           d_coef[im][iw] = (1. - rx)*(1. - qx)*(1. - px) * ad000 + rx * qx * px * ad111
                            + rx * (1. - qx) * px * ad101 + (1. - rx) * qx * (1. - px) * ad010
                            + rx * qx * (1. - px) * ad110 + (1. - rx)*(1. - qx) * px * ad001
                            + (1. - rx) * qx * px * ad011 + rx * (1. - qx)*(1. - px) * ad100;

           e_coef[im][iw] = (1. - rx)*(1. - qx)*(1. - px) * ae000 + rx * qx * px * ae111
                            + rx * (1. - qx) * px * ae101 + (1. - rx) * qx * (1. - px) * ae010
                            + rx * qx * (1. - px) * ae110 + (1. - rx)*(1. - qx) * px * ae001
                            + (1. - rx) * qx * px * ae011 + rx * (1. - qx)*(1. - px) * ae100;

            }
        }
    }

    // root finding is quadratic, but table includes cubic term, make sure it's zero
    if (fabs(d_coef[modnum][iwnir_l]) > 1e-9) {
        printf("non zero cubic term found in longest NIR wavelength of aerosol table. Zia!!\n");
        exit(1);
    }

    /* return pointers to coeffs for this geometry */
    *a = &a_coef[modnum][0];
    *b = &b_coef[modnum][0];
    *c = &c_coef[modnum][0];
    *d = &d_coef[modnum][0];
    *e = &e_coef[modnum][0];


    return;
}

/* ---------------------------------------------------------------------------------------- */
/* model_select_ahmad() - select two aerosol models whose epsilon values bracket the        */
/*                        the observed ms epsilon, eps_obs                                  */
/*                                                                                          */
/* Z Ahmad July 2014.                                                                       */

/* ---------------------------------------------------------------------------------------- */
int comp_epsilonT(epsilonTstr *x, epsilonTstr *y) {
    return (x->eps_obs < y->eps_obs ? -1 : 1);
}

void model_select_ahmad(int32_t nmodels, int32_t *mindx, float eps_pred[], float eps_obs, int32_t *modmin,
        int32_t *modmax, float *modrat) {
    static epsilonTstr epsilonT[MAXAERMOD];

    int im, im1, im2;

    // set-up table to keep model epsilon and model number pairs together

    for (im = 0; im < nmodels; im++) {
        epsilonT[im].eps_obs = eps_pred[im];
        epsilonT[im].modnum = im;
        /*     printf("%d %7.4f %7.4f\n",im,eps_pred[im],eps_obs);              */
    }

    // sort into ascending model epsilon order

    qsort(epsilonT, nmodels, sizeof (epsilonTstr),
            (int (*)(const void *, const void *)) comp_epsilonT);

    // find bounding epsilon indices in table

    for (im = 0; im < nmodels; im++) {
        if (eps_obs < epsilonT[im].eps_obs)
            break;
    }

    //if(im==0) //no lower bounding by M. Zhang
    //{
    //	*modmin=-1;
    //	return;
    //}


    im1 = MAX(MIN(im - 1, nmodels - 1), 0);
    im2 = MAX(MIN(im, nmodels - 1), 0);


    // convert table indices to model indices of the input order
    // compute model weighting

    *modmin = epsilonT[im1].modnum;
    *modmax = epsilonT[im2].modnum;
    *modrat = (eps_obs - epsilonT[im1].eps_obs) / (epsilonT[im2].eps_obs - epsilonT[im1].eps_obs);
    /*    printf("%d %d %7.4f %7.4f %7.4f\n",im1,im2,eps_obs,epsilonT[im1].eps_obs,epsilonT[im2].eps_obs); */

    // If eps_obs is higer or lower than epsilon from table, then set the weight to 1
    if (*modmin == *modmax)
        *modrat = 1;

    return;
}


/*------------------------------------------------------------------------------------------ */
/* comp_rhoa_ms_eps() -  for a given model, computes the AOT and aerosol reflectance         */
/*                                                                                           */
/* Z ahmad July2014                                                                          */

/*------------------------------------------------------------------------------------------ */
int comp_rhoa_ms_eps(int32_t nwave, float wave[], geom_str *geom,
        float tau_iwnir_l, int32_t modl, float tau_pred[], float rho_pred[])
 {
    float *ac, *bc, *cc, *dc, *ec;
    float ext_modl[nwave];
    float lg_tau_pred[nwave];
    float lg_rho_pred[nwave];
    int iw, iwtab;

    /* get the coefficients for lg_rho vs lg_aot  */

    // Zia's function ---> something is wrong
    ms_eps_coef(modl, iwnir_l, wave, geom, &ac, &bc, &cc, &dc, &ec);
    // ms_eps_coef_cal(modl,nwave,solz,senz,phi,&ac,&bc,&cc,&dc,&ec);


    /* get the extinction coefficients and compute AOT at all wavelengths */

    for (iw = 0; iw < nwave; iw++) {
        iwtab = iwatab[iw];
        ext_modl[iw] = aertab->model[modl]->extc[iwtab];
    }

    /*   printf("tau_pred[iw],tau_iwnir_l\n");    */
    for (iw = 0; iw < nwave; iw++) {
        tau_pred[iw] = (ext_modl[iw] / ext_modl[iwnir_l]) * tau_iwnir_l;
        lg_tau_pred[iw] = log(tau_pred[iw]);
    }

    /* compute rho_pred */

    for (iw = 0; iw < nwave; iw++) {
        iwtab = iwatab[iw];
        lg_rho_pred[iw] = ac[iwtab] +
                bc[iwtab] * lg_tau_pred[iw] +
                cc[iwtab] * pow(lg_tau_pred[iw], 2) +
                dc[iwtab] * pow(lg_tau_pred[iw], 3) +
                ec[iwtab] * pow(lg_tau_pred[iw], 4);
        rho_pred[iw] = exp(lg_rho_pred[iw]);
    }


    return (0);
}


/*------------------------------------------------------------------------------------------ */
/* comp_rhoa_ms_eps() -  for a given model, computes the AOT and aerosol reflectance         */
/*                                                                                           */
/* Z ahmad July2014                                                                          */
/* Modified by Amir Ibrahim January 2015 for linear space coefficients                       */

/*------------------------------------------------------------------------------------------ */
int comp_rhoa_ms_eps_lin(int32_t nwave, float wave[], geom_str *geom,
        float tau_iwnir_l, int32_t modl, float tau_pred[], float rho_pred[])
 {
    float *ac, *bc, *cc, *dc, *ec;
    float ext_modl[nwave];
    float lg_tau_pred[nwave];
    float lg_rho_pred[nwave];
    int iw, iwtab;

    /* get the coefficients for lg_rho vs lg_aot  */

    // Zia's function ---> something is wrong
    ms_eps_coef(modl, iwnir_l, wave, geom, &ac, &bc, &cc, &dc, &ec);
    // ms_eps_coef_cal(modl,nwave,geom,&ac,&bc,&cc,&dc,&ec);


    /* get the extinction coefficients and compute AOT at all wavelengths */

    for (iw = 0; iw <= nwave; iw++) {
        iwtab = iwatab[iw];
        ext_modl[iw] = aertab->model[modl]->extc[iwtab];
    }

    /*   printf("tau_pred[iw],tau_iwnir_l\n");    */
    for (iw = 0; iw <= nwave; iw++) {
        tau_pred[iw] = (ext_modl[iw] / ext_modl[iwnir_l]) * tau_iwnir_l;
        lg_tau_pred[iw] = (tau_pred[iw]);
    }

    /* compute rho_pred */

    for (iw = 0; iw <= nwave; iw++) {
        iwtab = iwatab[iw];
        lg_rho_pred[iw] = ac[iwtab] +
                bc[iwtab] * lg_tau_pred[iw] +
                cc[iwtab] * pow(lg_tau_pred[iw], 2) +
                dc[iwtab] * pow(lg_tau_pred[iw], 3) +
                ec[iwtab] * pow(lg_tau_pred[iw], 4);
        rho_pred[iw] = (lg_rho_pred[iw]);
    }


    return (0);
}


/*----------------------------------------------------------------------------------------------*/
/* spectral_matching() - calculates the best aerosol model that fits aerosol reflectance   from */
/* Red to SWIR channels. The spectral matching is done for 3 rh sets, and it assumes 1 band  in */
/* the NIR to be perfectly calibrated (i.e. 869 nm)                                             */
/* Amir Ibrahim - July 2016                                                                     */

/*----------------------------------------------------------------------------------------------*/
int spectral_matching(int32_t sensorID, float wave[], int32_t nwave,
        int32_t iwnir_s, int32_t iwnir_l, int32_t nmodels, int32_t mindx1[],
        int32_t mindx2[], int32_t mindx3[], geom_str *geom,
        float wv, float rhoa[], float rho_aer[], int32_t *mod1_indx,
        int32_t *mod2_indx, float *weight, float tau_pred_min[],
        float tau_pred_max[], float tg_sol_sm[], float tg_sen_sm[],
        float Lt_sm[], int32_t ip) {


    float *ac, *bc, *cc, *dc, *ec;
    float ext_iwnir_cal;
    double ax, bx, cx, fx;

    int iw, i;

    int status = 0;


    float tau_iwnir_cal;
    float *lg_tau_all_wav = (float*) malloc(nwave * sizeof (float));
    float *lg_rho_all_wav_pred = (float*) malloc(nwave * sizeof (float));
    float **tau_all_wav = (float**) malloc(nwave * sizeof (float*));
    float **rho_all_wav_pred = (float**) malloc(nwave * sizeof (float*));
    ;


    for (i = 0; i < nwave; i++) {
        tau_all_wav[i] = (float*) malloc((nmodels) * sizeof (float));
        rho_all_wav_pred[i] = (float*) malloc((nmodels) * sizeof (float));
    }

    float *ext_all_wav = (float*) malloc(nwave * sizeof (float));

    int iwtab, im, modl;
    float lg_tau_iwnir_cal;

    double *diff = malloc(nwave * sizeof (double));
    double *chi = malloc(nmodels * sizeof (double));
    double *chi_temp = malloc(nmodels * sizeof (double));

    double diff_2;
    double min1, min2;
    int min1_indx, min2_indx;

    //  float *SNR   =  malloc(nwave * sizeof(float));
    float *noise = malloc(nwave * sizeof (float));
    static float scaled_lt;
    int iwtab_l, iwtab_cal;


    if (sensorID == MODISA) {
        /* These are the noise coefficient for hmodisa from Erdem */

        float ltSnrFitCoef[16][2] = {
            {/*412:C0,C1*/0.05499859, 0.00008340},
            {/*443:*/0.02939470, 0.00009380},
            {/*469:*/0.11931482, 0.00008195},
            {/*488:*/0.01927545, 0.00009450},
            {/*531:*/0.01397522, 0.00010040},
            {/*547:*/0.01139088, 0.00016480},
            {/*555*/0.08769538, 0.00007000},
            {/*645:*/0.10406925, 0.00008533},
            {/*667L*/0.00496291, 0.00014050},
            {/*678L*/0.00427147, 0.00013160},
            {/*748*/0.00416994, 0.00021250},
            {/*859*/0.04055895, 0.000197550},
            {/*869*/0.00312263, 0.00018600},
            {/*1240*/0.07877732, 0.00049940},
            {/*1640*/100000.1, 100000000.1},
            {/*2130*/0.00628912, 0.00021160}
        };

        /* scaled noise based on MODISA spatial resolution */
        float snr_scale[16] = {
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 1, 2, 2, 2
        };
        /* Noise model */
        for (i = 0; i < nwave; i++) {
            scaled_lt = Lt_sm[ip * nwave + i];
            noise[i] = (ltSnrFitCoef[i][0] + ltSnrFitCoef[i][1] * scaled_lt * snr_scale[i]);
            //SNR[i] = scaled_lt * snr_mult / noise[i]; //noise model based on
        }
    } else {
        for (i = 0; i < nwave; i++)
            //SNR[i] = 1;
            noise[i] = 1;
    }



    int32_t *mindx = (int32_t*) malloc(nmodels * sizeof (int32_t));
    if (nmodels == 1) {
        *mindx = 0;
    } else {
        if (mindx1[0] == mindx2[0]) {
            nmodels = 10;
            for (i = 0; i < 10; i++)
                mindx[i] = mindx1[i];
        } else {
            if (nmodels < 30) {
                for (i = 0; i < 10; i++) {
                    mindx[i] = mindx1[i];
                    mindx[i + 10] = mindx2[i];
                }
            } else if (mindx3[0] >= mindx1[0]) {
                for (i = 0; i < 10; i++) {
                    mindx[i] = mindx1[i];
                    mindx[i + 10] = mindx2[i];
                    mindx[i + 20] = mindx3[i];
                }
            } else {
                for (i = 0; i < 10; i++) {
                    mindx[i] = mindx3[i];
                    mindx[i + 10] = mindx1[i];
                    mindx[i + 20] = mindx2[i];
                }
            }
        }
    }

    // compute MS epsilon for all nmodels.  note that nmodels may be only a subset of
    // the full model suite.  mindx maps from subset index to full index in suite
    /*
        float SNR[nwave];
        for (i=0;i<9;i++)
            SNR[i] = 0;
        SNR[10] = 1;
        SNR[11] = 0.3;
        SNR[12] = 1;
        SNR[13] = 0.5;
        SNR[14] = 0;
        SNR[15] = 0.2;
     */
    // printf("nmodels=%d\n",nmodels);
    for (im = 0; im < nmodels; im++) {

        modl = mindx[im];
        /* compute AOT at longest aerosol wavelength (iwnir_l) */
        int iwnir_cal;

        if (sensorID == AVIRIS) {
            iwnir_cal = 54; // this is hard coded (need to change) --> Amir
            iwtab_cal = iwatab[iwnir_cal];
        } else if (sensorID == MODISA) {
            iwnir_cal = 12;
            iwtab_cal = iwatab[iwnir_cal];

        } else {
            printf("Error: Sensor unkown %d\n", sensorID);
            return (1);
        }
        //printf("Using %f nm to estimate optical thickness\n",wave[iwnir_cal]);

        //ms_eps_coef_cal(modl,nwave,geom,&ac,&bc,&cc,&dc,&ec);
        ms_eps_coef(modl, iwnir_l, wave, geom, &ac, &bc, &cc, &dc, &ec);

        iwtab_l = iwatab[iwnir_l];
        //        iwtab_s = iwatab[iwnir_s];

        ax = (double) ac[iwtab_cal] - log((double) rhoa[iwnir_cal]);
        bx = (double) bc[iwtab_cal];
        cx = (double) cc[iwtab_cal];

        fx = bx * bx - 4.0 * ax*cx;
        if (fx > 0.0 && cx != 0.0) {
            lg_tau_iwnir_cal = 0.5 * (-bx + sqrt(fx)) / cx;
            tau_iwnir_cal = exp(lg_tau_iwnir_cal);
            status = 0;
        } else {
            status = 1;
            break;
        }

        /* compute AOT at shortest aerosol wavelength (iwnir_s) */

        /* compute reflectance at all wavelength */

        ext_iwnir_cal = aertab->model[modl]->extc[iwtab_cal];
        for (iw = 0; iw <= iwtab_l - 1; iw++) {
            iwtab = iwatab[iw];
            ext_all_wav[iw] = aertab->model[modl]->extc[iwtab];
            tau_all_wav[iw][im] = (ext_all_wav[iw] / ext_iwnir_cal) * tau_iwnir_cal;
            lg_tau_all_wav[iw] = log(tau_all_wav[iw][im]);
            lg_rho_all_wav_pred[iw] = ac[iwtab] + bc[iwtab] * lg_tau_all_wav[iw] +
                    cc[iwtab] * pow(lg_tau_all_wav[iw], 2) +
                    dc[iwtab] * pow(lg_tau_all_wav[iw], 3) +
                    ec[iwtab] * pow(lg_tau_all_wav[iw], 4);
            rho_all_wav_pred[iw][im] = exp(lg_rho_all_wav_pred[iw]);


        }
        for (iw = iwnir_s; iw <= iwnir_l; iw++) {
            // calculate the rms error
            diff[iw] = pow((rhoa[iw] - rho_all_wav_pred[iw][im]), 2);
            diff_2 += diff[iw]* (1 / noise[iw]) *(tg_sol_sm[ip * nwave + iw] * tg_sol_sm[ip * nwave + iw]);

        }
        chi[im] = sqrt(diff_2 / (iwnir_l - iwnir_s));
        diff_2 = 0.0;
        //  printf("Model number = %d",im);
    }

    min1 = 1000;
    min2 = 1000;
    min1_indx = 0;
    min2_indx = 0;

    for (i = 0; i < nmodels; i++)
        chi_temp[i] = chi[i];

    qsort(chi_temp, nmodels, sizeof (chi_temp[0]), cmpfunc);
    min1 = chi_temp[0];
    min2 = chi_temp[1];

    for (i = 0; i < nmodels; i++) {
        if (chi[i] == min1)
            min1_indx = i;
        else if (chi[i] == min2)
            min2_indx = i;
    }


    *mod1_indx = 0;
    *mod2_indx = 0;

    if (min2 == 0)
        min2 = 0.001;
    //    double weight1 = sqrt(pow(aertab->model[min2_indx]->extc[iwtab_cal] - aertab->model[min1_indx]->extc[iwtab_cal],2))/(aertab->model[min1_indx]->extc[iwtab_cal]);
    double weight1 = ((rhoa[iwnir_s] / rhoa[iwnir_l]) - (rho_all_wav_pred[iwnir_s][min1_indx] / rho_all_wav_pred[iwnir_l][min1_indx])) / ((rho_all_wav_pred[iwnir_s][min2_indx] / rho_all_wav_pred[iwnir_l][min2_indx])-(rho_all_wav_pred[iwnir_s][min1_indx] / rho_all_wav_pred[iwnir_l][min1_indx]));
    //double weight1 = (min2-min1)/min1 ;
    *weight = weight1;
    if (nmodels == 1) {
        *weight = 1;
        weight1 = 1;
    }

    *mod1_indx = mindx[min1_indx];
    *mod2_indx = mindx[min2_indx];
    for (iw = 0; iw <= iwnir_l; iw++) {
        rho_aer[iw] = (1 - weight1) * rho_all_wav_pred[iw][min1_indx] + (weight1) * rho_all_wav_pred[iw][min2_indx];
        tau_pred_min[iw] = tau_all_wav[iw][min1_indx];
        tau_pred_max[iw] = tau_all_wav[iw][min2_indx];
    }

    free(lg_tau_all_wav);
    free(lg_rho_all_wav_pred);
    for (i = 0; i < nwave; i++) {
        free(tau_all_wav[i]);
        free(rho_all_wav_pred[i]);
    }
    free(tau_all_wav);
    free(rho_all_wav_pred);
    free(ext_all_wav);
    free(chi_temp);
    free(chi);
    free(diff);
    free(mindx);
    //   free(SNR);
    free(noise);
    //free(tg_sol_sm);
    // free(tg_sen_sm);


    return (status);
}
/*------------------------------------------------------------------------------------------*/
/* ahmad_atm_corr() - compute aerosol reflectance at all wavelengths                        */
/*                                                                                          */
/* Z Ahmad. July 2014                                                                       */

/*----------------------------------------------------------------------------------------  */
int ahmad_atm_corr(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
        int32_t nmodels, int32_t mindx[],
        geom_str *geom, float wv, float rhoa[],
        int32_t *modmin, int32_t *modmax, float *modrat, float *epsnir,
        float tau_pred_max[], float tau_pred_min[], float rho_pred_max[], float rho_pred_min[],
        float tau_aer[], float rho_aer[]) {


    float *ac, *bc, *cc, *dc, *ec;
    float ext_iwnir_l, ext_iwnir_s;
    double ax, bx, cx, fx;
    float tau_iwnir_l[nmodels];
    float lg_tau_iwnir_s[nmodels], tau_iwnir_s[nmodels];
    float lg_rho_iwnir_s_pred[nmodels], rho_iwnir_s_pred[nmodels];
    float eps_pred[nmodels];
    int im, modl;
    int im1, im2;
    float mwt;
    float lg_tau_iwnir_l;
    int iwtab_l, iwtab_s;


    static float eps_obs;

    int iw;

    int status = 0.0;


    /* compute the observed epsilon */

    eps_obs = rhoa[iwnir_s] / rhoa[iwnir_l];

    /*             printf("rho_869,rho_748,eps_obs\n");                    */
    /*    printf("%10.5f %10.5f %10.5f\n",rhoa[iwnir_l],rhoa[iwnir_s],eps_obs);   */


    // compute MS epsilon for all nmodels.  note that nmodels may be only a subset of
    // the full model suite.  mindx maps from subset index to full index in suite


    for (im = 0; im < nmodels; im++) {

        modl = mindx[im]; // index in full model suite, needed for coefficient look-up

        /* compute AOT at longest aerosol wavelength (iwnir_l) */

        // Zia's function ---> something is wrong
        //ms_eps_coef(modl,iwnir_l,wave,solz,senz,phi,&ac,&bc,&cc,&dc,&ec);

        //        ms_eps_coef_cal(modl,nwave,solz,senz,phi,&ac,&bc,&cc,&dc,&ec);
        ms_eps_coef(modl, iwnir_l, wave, geom, &ac, &bc, &cc, &dc, &ec);
        iwtab_l = iwatab[iwnir_l];
        iwtab_s = iwatab[iwnir_s];

        ax = (double) ac[iwtab_l] - log((double) rhoa[iwnir_l]);
        bx = (double) bc[iwtab_l];
        cx = (double) cc[iwtab_l];

        fx = bx * bx - 4.0 * ax*cx;
        if (fx > 0.0 && cx != 0.0) {
            lg_tau_iwnir_l = 0.5 * (-bx + sqrt(fx)) / cx;
            tau_iwnir_l[im] = exp(lg_tau_iwnir_l);
        } else {
            status = 1;
            break;    // TODO: need to think about it
        }

        /* compute AOT at shortest aerosol wavelength (iwnir_s) */

        ext_iwnir_l = aertab->model[modl]->extc[iwtab_l];
        ext_iwnir_s = aertab->model[modl]->extc[iwtab_s];

        tau_iwnir_s[im] = (ext_iwnir_s / ext_iwnir_l) * tau_iwnir_l[im];
        lg_tau_iwnir_s[im] = log(tau_iwnir_s[im]);


        /* compute reflectance at (iwnir_s) */

        lg_rho_iwnir_s_pred[im] = ac[iwtab_s] + bc[iwtab_s] * lg_tau_iwnir_s[im] +
                cc[iwtab_s] * pow(lg_tau_iwnir_s[im], 2) +
                dc[iwtab_s] * pow(lg_tau_iwnir_s[im], 3) +
                ec[iwtab_s] * pow(lg_tau_iwnir_s[im], 4);

        rho_iwnir_s_pred[im] = exp(lg_rho_iwnir_s_pred[im]);

        /* compute model epsilon */

        eps_pred[im] = rho_iwnir_s_pred[im] / rhoa[iwnir_l];

    }


    *epsnir = eps_obs;


    /* now, find the bounding models, but skip this if we only have one */
    /* model (as when vicariously calibrating)                          */

    if (nmodels > 1) {

        /* locate two model_epsilons that bracket the observed epsilon */
        /* this will return the model numbers for the subset of models */

        model_select_ahmad(nmodels, mindx, eps_pred, eps_obs, &im1, &im2, &mwt);

    } else {

        /* single model case (allows fixed model by limiting model suite) */

        im1 = 0;
        im2 = 0;
        mwt = 0.0;
    }

    /* map selection to full suite for subsequent processing and return */

    //if(im1<0)  // no lower bounding, M. Zhang
    //{
    //	return 1;
    //	status=1;
    //}


    *modmin = mindx[im1];
    *modmax = mindx[im2];
    *modrat = mwt;

    /* compute tau_pred and rho_predicted */

    comp_rhoa_ms_eps(nwave, wave, geom, tau_iwnir_l[im1], *modmin, tau_pred_min, rho_pred_min);
    comp_rhoa_ms_eps(nwave, wave, geom, tau_iwnir_l[im2], *modmax, tau_pred_max, rho_pred_max);

    /* compute weighted tau_aer and rho_aer */

    for (iw = 0; iw < nwave; iw++) {
        tau_aer[iw] = (1.0 - mwt) * tau_pred_min[iw] + mwt * tau_pred_max[iw];
        rho_aer[iw] = (1.0 - mwt) * rho_pred_min[iw] + mwt * rho_pred_max[iw];
    }


    return (status);
}

/*------------------------------------------------------------------------------------------*/
/* ahmad_atm_corr_lin() - compute aerosol reflectance at all wavelengths in linear space    */
/*                                                                                          */
/* Z Ahmad. July 2014                                                                       */
/* Modified to linear spaced coefficient by Amir Ibrahim January 2017                       */

/*----------------------------------------------------------------------------------------  */
int ahmad_atm_corr_lin(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
        int32_t nmodels, int32_t mindx[],
        geom_str *geom, float wv, float rhoa[],
        int32_t *modmin, int32_t *modmax, float *modrat, float *epsnir,
        float tau_pred_max[], float tau_pred_min[], float rho_pred_max[], float rho_pred_min[],
        float tau_aer[], float rho_aer[]) {


    float *ac, *bc, *cc, *dc, *ec;
    float ext_iwnir_l, ext_iwnir_s;
    double ax, bx, cx, fx;
    float tau_iwnir_l[nmodels];
    float lg_tau_iwnir_s[nmodels], tau_iwnir_s[nmodels];
    float lg_rho_iwnir_s_pred[nmodels], rho_iwnir_s_pred[nmodels];
    float eps_pred[nmodels];
    int im, modl;
    int im1, im2;
    float mwt;
    float lg_tau_iwnir_l;
    int iwtab_l, iwtab_s;


    static float eps_obs;

    int iw;

    int status = 0.0;


    /* compute the observed epsilon */

    eps_obs = rhoa[iwnir_s] / rhoa[iwnir_l];

    /*             printf("rho_869,rho_748,eps_obs\n");                    */
    /*    printf("%10.5f %10.5f %10.5f\n",rhoa[iwnir_l],rhoa[iwnir_s],eps_obs);   */


    // compute MS epsilon for all nmodels.  note that nmodels may be only a subset of
    // the full model suite.  mindx maps from subset index to full index in suite


    for (im = 0; im < nmodels; im++) {

        modl = mindx[im]; // index in full model suite, needed for coefficient look-up

        /* compute AOT at longest aerosol wavelength (iwnir_l) */

        // Zia's function ---> something is wrong
        //ms_eps_coef(modl,iwnir_l,wave,solz,senz,phi,&ac,&bc,&cc,&dc,&ec);

        //        ms_eps_coef_cal(modl,nwave,solz,senz,phi,&ac,&bc,&cc,&dc,&ec);
        ms_eps_coef(modl, iwnir_l, wave, geom, &ac, &bc, &cc, &dc, &ec);
        iwtab_l = iwatab[iwnir_l];
        iwtab_s = iwatab[iwnir_s];

        ax = (double) ac[iwtab_l]-(double) rhoa[iwnir_l];
        bx = (double) bc[iwtab_l];
        cx = (double) cc[iwtab_l];

        fx = bx * bx - 4.0 * ax*cx;
        if (fx > 0.0 && cx != 0.0) {
            lg_tau_iwnir_l = 0.5 * (-bx + sqrt(fx)) / cx;
            tau_iwnir_l[im] = (lg_tau_iwnir_l);
        } else {
            status = 1;
            break;
        }

        /* compute AOT at shortest aerosol wavelength (iwnir_s) */

        ext_iwnir_l = aertab->model[modl]->extc[iwtab_l];
        ext_iwnir_s = aertab->model[modl]->extc[iwtab_s];

        tau_iwnir_s[im] = (ext_iwnir_s / ext_iwnir_l) * tau_iwnir_l[im];
        lg_tau_iwnir_s[im] = (tau_iwnir_s[im]);

        if (ax > 1e5) {
            printf("\nErroneous aerosol option, %d\n", aer_opt);
            printf("You are using a log-space LUT. The aerosol LUT coefficients need to be in linear-space. Use aer_opt=-18 instead. \n");
            exit(FATAL_ERROR);
        }


        /* compute reflectance at (iwnir_s) */

        lg_rho_iwnir_s_pred[im] = ac[iwtab_s] + bc[iwtab_s] * lg_tau_iwnir_s[im] +
                cc[iwtab_s] * pow(lg_tau_iwnir_s[im], 2) +
                dc[iwtab_s] * pow(lg_tau_iwnir_s[im], 3) +
                ec[iwtab_s] * pow(lg_tau_iwnir_s[im], 4);

        rho_iwnir_s_pred[im] = (lg_rho_iwnir_s_pred[im]);

        /* compute model epsilon */

        eps_pred[im] = rho_iwnir_s_pred[im] / rhoa[iwnir_l];

    }


    *epsnir = eps_obs;


    /* now, find the bounding models, but skip this if we only have one */
    /* model (as when vicariously calibrating)                          */

    if (nmodels > 1) {

        /* locate two model_epsilons that bracket the observed epsilon */
        /* this will return the model numbers for the subset of models */

        model_select_ahmad(nmodels, mindx, eps_pred, eps_obs, &im1, &im2, &mwt);

    } else {

        /* single model case (allows fixed model by limiting model suite) */

        im1 = 0;
        im2 = 0;
        mwt = 0.0;
    }

    /* map selection to full suite for subsequent processing and return */

    *modmin = mindx[im1];
    *modmax = mindx[im2];
    *modrat = mwt;

    /* compute tau_pred and rho_predicted */

    comp_rhoa_ms_eps_lin(nwave, wave, geom, tau_iwnir_l[im1], *modmin, tau_pred_min, rho_pred_min);
    comp_rhoa_ms_eps_lin(nwave, wave, geom, tau_iwnir_l[im2], *modmax, tau_pred_max, rho_pred_max);

    /* compute weighted tau_aer and rho_aer */

    for (iw = 0; iw < nwave; iw++) {
        tau_aer[iw] = (1.0 - mwt) * tau_pred_min[iw] + mwt * tau_pred_max[iw];
        rho_aer[iw] = (1.0 - mwt) * rho_pred_min[iw] + mwt * rho_pred_max[iw];
    }


    return (status);
}

/* ---------------------------------------------------------------------------------------- */
/* model_phase() - return phase function at model wavelengths at the input geometry.        */
/*                                                                                          */
/* This is effectively a C version of load_ss.f by M. Wang.  The program optimizes for      */
/* multiple calls at the same geometry by computing for all models on the first call with a */
/* new geometry.  It returns a pointer to the internal static array for the requested model.*/
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */
/*  W. Robinson, SAIC modified to use band-dependent viewing geometry, if avail             */

/* ---------------------------------------------------------------------------------------- */
float *model_phase(int modnum, geom_str *geom) {
    static float nw = 1.334;

    static int computed[MAXMODEL];
    static float lastsolz = -999.;
    static float lastsenz = -999.;
    static float lastphi = -999.;
    static float *phase[MAXMODEL];
    static float *fres1, *fres2;
    static float *scatt1, *scatt2;
    static float scatt1_cnst, scatt2_cnst, fres1_cnst, fres2_cnst;
    static float *scatt1_ar, *scatt2_ar, *fres1_ar, *fres2_ar;
    static int firstCall = 1, gmult = 1;

    float phase1, phase2;
    int im, iw, ig;

    if (firstCall == 1) {
        firstCall = 0;
        for (im = 0; im < MAXMODEL; im++) {
            if ((phase[im] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for phase.\n");
                exit(1);
            }
        }
        /* set up geometry based scattering and fresnel */
        if ((geom->gmult == 0) || (interpol == 1)) {
            gmult = 0;
            scatt1 = &scatt1_cnst;
            scatt2 = &scatt2_cnst;
            fres1 = &fres1_cnst;
            fres2 = &fres2_cnst;
        } else {
            gmult = 1;
            if ((scatt1_ar =
                    (float *) malloc(aertab->nwave * sizeof (float))) == NULL) {
                printf("Unable to allocate space for scatt1.\n");
                exit(1);
            }
            if ((scatt2_ar =
                    (float *) malloc(aertab->nwave * sizeof (float))) == NULL) {
                printf("Unable to allocate space for scatt2.\n");
                exit(1);
            }
            if ((fres1_ar =
                    (float *) malloc(aertab->nwave * sizeof (float))) == NULL) {
                printf("Unable to allocate space for fres1.\n");
                exit(1);
            }
            if ((fres2_ar =
                    (float *) malloc(aertab->nwave * sizeof (float))) == NULL) {
                printf("Unable to allocate space for fres2.\n");
                exit(1);
            }
            scatt1 = scatt1_ar;
            scatt2 = scatt2_ar;
            fres1 = fres1_ar;
            fres2 = fres2_ar;
        }
    }

    /* recalculate only if geometry changes */

    if ((geom->solz[0] != lastsolz) || (geom->senz[0] != lastsenz) ||
            (geom->phi[0] != lastphi)) {

        float temp;

        for (im = 0; im < aertab->nmodel; im++)
            computed[im] = 0;

        /* determine scattering angles (direct and surface reflected) */
        for (iw = 0; iw < aertab->nwave; iw++) {
            ig = gmult * iw;
            temp = sqrt((1.0 - geom->csenz[ig] * geom->csenz[ig]) *
                    (1.0 - geom->csolz[ig] * geom->csolz[ig])) *
                    cos(geom->phi[ig] / radeg);
            scatt1[iw] = acos(
                    MAX(-geom->csenz[ig] * geom->csolz[ig] + temp, -1.0)) * radeg;
            scatt2[iw] = acos(
                    MIN(geom->csenz[ig] * geom->csolz[ig] + temp, 1.0)) * radeg;

            /* compute Fresnel coefficients */
            fres1[iw] = fresnel_coef(geom->csenz[ig], nw);
            fres2[iw] = fresnel_coef(geom->csolz[ig], nw);
            if (gmult == 0) break;
        }

        lastsolz = geom->solz[0];
        lastsenz = geom->senz[0];
        lastphi = geom->phi[0];
    }

    if (!computed[modnum]) {

        im = modnum;
        computed[modnum] = 1;

        /* compute phase function for this geometry, all models */
        for (iw = 0; iw < aertab->nwave; iw++) {
            ig = gmult * iw;
            splint(aertab->scatt,
                    &aertab->model[im]->lnphase[iw][0],
                    &aertab->model[im]->d2phase[iw][0],
                    aertab->nscatt, scatt1[ig], &phase1);
            splint(aertab->scatt,
                    &aertab->model[im]->lnphase[iw][0],
                    &aertab->model[im]->d2phase[iw][0],
                    aertab->nscatt, scatt2[ig], &phase2);
            //          incident diffuse   reflected   diff  dir
            phase[im][iw] = exp(phase1) +
                    exp(phase2)*(fres1[ig] + fres2[ig]);
        }
    }

    return (&phase[modnum][0]);
}


/* ---------------------------------------------------------------------------------------- */
/* aeroob_cf() - out-of-band water-vapor scale factor                                       */

/* ---------------------------------------------------------------------------------------- */
float aeroob_cf(int modnum, geom_str *geom) {
    static int firstCall = 1;
    static int iw1;
    static int iw2;

    float *phase;
    float rhoas1, rhoas2;
    float eps;
    float cf;

    if (firstCall) {
        iw1 = windex(765, aertab->wave, aertab->nwave);
        iw2 = windex(865, aertab->wave, aertab->nwave);
        if (iw1 == iw2) iw1--;
        firstCall = 0;
    }

    phase = model_phase(modnum, geom);
    rhoas1 = aertab->model[modnum]->albedo[iw1] * phase[iw1] * aertab->model[modnum]->extc[iw1];
    rhoas2 = aertab->model[modnum]->albedo[iw2] * phase[iw2] * aertab->model[modnum]->extc[iw2];
    eps = rhoas1 / rhoas2;
    cf = log(eps) / (aertab->wave[iw2] - aertab->wave[iw1]);

    return (cf);
}


/* ---------------------------------------------------------------------------------------- */
/* aeroob - out-of-band water-vapor correction                                              */

/* ---------------------------------------------------------------------------------------- */
float aeroob(int32_t sensorID, int32_t iw, float airmass, float cf, float wv) {
    static float *a01;
    static float *a02;
    static float *a03;
    static float *a04;
    static float *a05;
    static float *a06;
    static float *a07;
    static float *a08;
    static float *a09;
    static float *a10;
    static float *a11;
    static float *a12;

    static int firstCall = 1;

    float f;

    if (firstCall) {
        firstCall = 0;
        rdsensorinfo(sensorID, evalmask, "oobwv01", (void **) &a01); /* coeff #1 per sensor wave */
        rdsensorinfo(sensorID, evalmask, "oobwv02", (void **) &a02);
        rdsensorinfo(sensorID, evalmask, "oobwv03", (void **) &a03);
        rdsensorinfo(sensorID, evalmask, "oobwv04", (void **) &a04);
        rdsensorinfo(sensorID, evalmask, "oobwv05", (void **) &a05);
        rdsensorinfo(sensorID, evalmask, "oobwv06", (void **) &a06);
        rdsensorinfo(sensorID, evalmask, "oobwv07", (void **) &a07);
        rdsensorinfo(sensorID, evalmask, "oobwv08", (void **) &a08);
        rdsensorinfo(sensorID, evalmask, "oobwv09", (void **) &a09);
        rdsensorinfo(sensorID, evalmask, "oobwv10", (void **) &a10);
        rdsensorinfo(sensorID, evalmask, "oobwv11", (void **) &a11);
        rdsensorinfo(sensorID, evalmask, "oobwv12", (void **) &a12);
        printf("\nLoading water-vapor correction coefficients.\n");
    }

    f = (a01[iw] + a02[iw] * airmass + cf * (a03[iw] + a04[iw] * airmass))
            + (a05[iw] + a06[iw] * airmass + cf * (a07[iw] + a08[iw] * airmass)) * wv
            + (a09[iw] + a10[iw] * airmass + cf * (a11[iw] + a12[iw] * airmass)) * wv*wv;

    return (f);
}


/* ---------------------------------------------------------------------------------------- */
/* rhoa_to_rhoas() - MS aerosol reflectance to SS aerosol reflectance                       */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
int rhoa_to_rhoas(int32_t sensorID, int modnum, geom_str *geom, float wv,
        float rhoa[], float wave[], int32_t nwave, int iw1, int iw2, float rhoas[]) {
    float *ac, *bc, *cc;
    double a, b, c;
    double f;
    int iw, ig;
    int iwtab;
    float cf;
    int status = 0;

    ss_to_ms_coef(modnum, geom, &ac, &bc, &cc);

    cf = aeroob_cf(modnum, geom);

    for (iw = iw1; iw <= iw2; iw++) {
        ig = iw * geom->gmult;
        if (rhoa[iw] < 1.e-20)
            rhoas[iw] = rhoa[iw];
        else {
            iwtab = iwatab[iw];
            a = (double) ac[iwtab];
            b = (double) bc[iwtab];
            c = (double) cc[iwtab];
            f = b * b - 4 * c * (a - log((double) rhoa[iw]));
            if (f > 0.00001) { // this was 0.0, but small values caused segfault (BAF, 9/2014)
                if (fabs(c) > 1.e-20) {
                    rhoas[iw] = exp(0.5 * (-b + sqrt(f)) / c);
                } else if (fabs(a) > 1.e-20 && fabs(b) > 1.e-20) {
                    rhoas[iw] = pow(rhoa[iw] / a, 1. / b);
                } else {
                    status = 1;
                    break;
                }
                rhoas[iw] = rhoas[iw] / aeroob(sensorID, iw, geom->airmass[ig], cf, wv);
                if (!isfinite(rhoas[iw]) || rhoas[iw] < 1.e-20) {
                    status = 1;
                    break;
                }
            } else {
                status = 1;
                break;
            }
        }
    }

    // return input values and failure status if any wavelengths failed

    if (status != 0) {
        for (iw = iw1; iw <= iw2; iw++)
            rhoas[iw] = rhoa[iw];
    }

    return (status);
}



/* ---------------------------------------------------------------------------------------- */
/* rhoas_to_rhoa() - SS aerosol reflectance to MS aerosol reflectance                       */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
void rhoas_to_rhoa(int32_t sensorID, int modnum, geom_str *geom, float wv,
        float rhoas[], float wave[], int32_t nwave, int iw1, int iw2, float rhoa[]) {
    float *ac, *bc, *cc;
    float a, b, c;
    float lnrhoas;
    int iw, ig;
    int iwtab;
    float cf;

    ss_to_ms_coef(modnum, geom, &ac, &bc, &cc);

    cf = aeroob_cf(modnum, geom);

    for (iw = iw1; iw <= iw2; iw++) {
        ig = iw * geom->gmult;

        /* these changes ensure that rhoa1 -> rhoas -> rhoa2 == rhoa1 */
        /* but, that changes everything, slightly (tau)               */

        if (rhoas[iw] < 1.e-20)
            rhoa[iw] = rhoas[iw];
        else {
            iwtab = iwatab[iw];
            a = ac[iwtab];
            b = bc[iwtab];
            c = cc[iwtab];
            lnrhoas = log(rhoas[iw] * aeroob(sensorID, iw, geom->airmass[ig], cf, wv));
            rhoa[iw] = exp(a + b * lnrhoas + c * lnrhoas * lnrhoas);
        }

        /*
          iwtab = iwatab[iw];
          a = ac[iwtab];
          b = bc[iwtab];
          c = cc[iwtab];
          lnrhoas = log(rhoas[iw]);
          rhoa[iw] = exp(a + b*lnrhoas + c*lnrhoas*lnrhoas)
         * aeroob(sensorID,iw, geom->airmass[ig],cf,wv);
         */
    }

    return;
}


/* ---------------------------------------------------------------------------------------- */
/* model_epsilon() - return model epsilon at input wavelengths for input geometry.          */
/*                                                                                          */
/* If the input wavelengths are not equal to the model wavelengths, the model epsilon will  */
/* be interpolated to the input wavelengths.  It is assumed that the longest (last) input   */
/* wavelength is equivalent to the longest wavelength of the model table.  Hence,           */
/* the function should always be called with the full array of input sensor wavelengths.    */
/*                                                                                          */
/* The program optimizes for multiple calls at the same geometry by computing for all       */
/* models on the first call with a new geometry.  It returns a pointer to the internal      */
/* static arrays of epsilon for the requested model.                               .        */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
float *model_epsilon(int modnum, int32_t iwnir_l, float wave[], int32_t nwave, geom_str *geom) {
    static float lastsolz = -999.;
    static float lastsenz = -999.;
    static float lastphi = -999.;
    static int32_t lastiwl = -999;
    static int lastmod = -999;
    static float *epsilon[MAXMODEL];
    static int firstCall = 1;
    int i;
    float maxwave;
    /* recalculate only if geometry changes */

    if (firstCall == 1) {
        firstCall = 0;
        maxwave = MAX(aertab->nwave, nwave);
        for (i = 0; i < MAXMODEL; i++) {
            if ((epsilon[i] = (float *) calloc(maxwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for epsilon.\n");
                exit(1);
            }

        }
    }

    /* recalculate only if geometry changes */

    if (modnum != lastmod || geom->solz[0] != lastsolz ||
            geom->senz[0] != lastsenz || geom->phi[0] != lastphi ||
            iwnir_l != lastiwl) {

        int iwnir = iwatab[iwnir_l];
        float *phase;
        float *lneps;
        float rhoas1, rhoas2;
        int im, iw, iwtab;

        if ((lneps = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
            printf("Unable to allocate space for lneps.\n");
            exit(1);
        }

        im = modnum;
        phase = model_phase(im, geom);
        for (iw = 0; iw < aertab->nwave; iw++) {
            rhoas1 = aertab->model[im]->albedo[iw] * phase[iw] * aertab->model[im]->extc[iw];
            rhoas2 = aertab->model[im]->albedo[iwnir] * phase[iwnir] * aertab->model[im]->extc[iwnir];
            epsilon[im][iw] = rhoas1 / rhoas2;
            if (interpol)
                lneps[iw] = log(epsilon[im][iw]);
        }
        if (interpol) {
            for (iw = 0; iw < nwave; iw++) {
                iwtab = iwatab[iw];
                if (aertab->wave[iwtab] != wave[iw] && wave[iw] > 0)
                    epsilon[im][iw] = exp(linterp(aertab->wave, lneps, aertab->nwave, wave[iw]));
                else
                    epsilon[im][iw] = exp(lneps[iwtab]);
            }
        }

        lastsolz = geom->solz[0];
        lastsenz = geom->senz[0];
        lastphi = geom->phi[0];
        lastiwl = iwnir_l;
        lastmod = modnum;
        free(lneps);

    }

    return (&epsilon[modnum][0]);
}


/* ---------------------------------------------------------------------------------------- */
/* model_select_wang() - M. Wang aerosol model selection process                 .          */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
int model_select_wang(int32_t sensorID, float wave[], int32_t nwave,
        int32_t nmodel, int32_t mindx[], geom_str *geom, float wv,
        float rhoa[], int32_t iwnir_s, int32_t iwnir_l, int32_t *modmin,
        int32_t *modmax, float *modrat, float *epsnir) {
    float *rhoas;
    float eps_ret [MAXMODEL];
    float eps_mod [MAXMODEL];
    float eps_err [MAXMODEL];
    int imod [MAXMODEL];
    int nmod = nmodel;
    float eps_ave;
    float *eps;
    float err_m;
    float err_p;
    int jm, im, iim;
    int eps_flg = 0;
    float wt;
    float tot_err;
    int itmp;

    *modmin = -1;
    *modmax = -1;
    *modrat = 0.0;
    *epsnir = 0.0;

    if ((rhoas = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas.\n");
        exit(1);
    }

    /* get theoretical and retrieved epsilon for each model, and */
    /* compute average retrieved epsilon                         */

    eps_ave = 0.0;
    for (jm = 0; jm < nmod; jm++) {

        im = mindx[jm];

        /* retrieve epsilon in NIR assuming each model */
        rhoa_to_rhoas(sensorID, im, geom, wv, rhoa, wave, nwave, iwnir_s, iwnir_l, rhoas);

        if (rhoas[iwnir_l] > 0.0000001)
            eps_ret[jm] = rhoas[iwnir_s] / rhoas[iwnir_l];
        else
            eps_ret[jm] = 0;

        /* get model epsilon for each model at this geometry */
        eps = model_epsilon(im, iwnir_l, wave, nwave, geom);
        eps_mod[jm] = eps[iwnir_s];

        eps_ave += eps_ret[jm];
    }
    if (isfinite(eps_ave))
        eps_ave /= nmod;
    else
        eps_ave = 1.0;


    /* determine average retrieved epsilon for the four retrievals which most */
    /* closely match the theoretical model epsilons. the model set is reduced */
    /* from the full suite to the final four by iterative outlier rejection.  */

    while (nmod > 4) {

        /* compute differences between retrieved and model epsilon */
        for (im = 0; im < nmodel; im++) {
            imod[im] = im;
            eps_err[im] = eps_ave - eps_mod[im];
        }

        /* sort model indices by smallest to largest absolute differences */
        for (im = 0; im < nmodel - 1; im++) {
            for (iim = im + 1; iim < nmodel; iim++)
                if (fabs(eps_err[imod[im]]) > fabs(eps_err[imod[iim]])) {
                    itmp = imod[im ];
                    imod[im ] = imod[iim];
                    imod[iim] = itmp;
                }
        }

        /* recompute average retrieved epsilon over the n-2 closest models  */
        /* averaging is done as a weighted mean with wt=1-|eps_err|/tot_err */
        /* note that the sum of the weights is equal to n-1                 */

        nmod = nmod - 2;

        tot_err = 0.0;
        for (iim = 0; iim < nmod; iim++) {
            im = imod[iim];
            tot_err += fabs(eps_err[im]);
        }

        eps_ave = 0.0;
        for (iim = 0; iim < nmod; iim++) {
            im = imod[iim];
            wt = 1.0 - fabs(eps_err[im]) / tot_err;
            eps_ave += eps_ret[im] * wt;
        }
        eps_ave /= (nmod - 1);
    }

    /* now select the two models which bracket eps_ave  */
    err_m = 0 - FLT_MAX;
    err_p = FLT_MAX;
    for (im = 0; im < nmodel; im++) {
        eps_err[im] = eps_ave - eps_mod[im];
        if (eps_err[im] >= 0.0) {
            if (eps_err[im] < err_p) {
                err_p = eps_err[im];
                *modmin = im;
            }
        } else {
            if (eps_err[im] > err_m) {
                err_m = eps_err[im];
                *modmax = im;
            }
        }
    }

    /* M. Wang's model cross-over correction */
    if (wang_modx && eps_mod[imc50] > eps_mod[imt99]) {
        if (*modmax == imt90)
            *modmin = imt99;
        else if (*modmax == imc50 && *modmin == imt99)
            *modmin = imc70;
        else if (*modmin == imm50)
            *modmax = imc50;
    }

    /* compute interpolation ratio */
    if (*modmin < 0) {
        /* no lower-bounding model */
        *modmin = *modmax;
        *modrat = 0.0;
        eps_flg = -1;
    } else if (*modmax < 0) {
        /* no upper-bounding model */
        *modmax = *modmin;
        *modrat = 0.0;
        eps_flg = 1;
    } else
        *modrat = (eps_ave - eps_mod[*modmin]) / (eps_mod[*modmax] - eps_mod[*modmin]);

    *modmin = mindx[*modmin];
    *modmax = mindx[*modmax];

    /* return retrieved epsilon */
    *epsnir = eps_ave;

    free(rhoas);

    return (eps_flg);
}


/* ---------------------------------------------------------------------------------------- */
/* model_select_angst() - Select model pair based on input Angstrom coefficient             */
/*                                                                                          */
/* B. Franz, 1 August 2004.                                                                 */

/* ---------------------------------------------------------------------------------------- */
int compalphaT(alphaTstr *x, alphaTstr *y) {
    return (x->angstrom < y->angstrom ? -1 : 1);
}

void model_select_angstrom(float angstrom, int32_t *modmin, int32_t *modmax, float *modrat) {
    static alphaTstr alphaT[MAXAERMOD];
    static int firstCall = 1;

    int im, im1, im2;

    if (firstCall) {

        int ib = windex(520, aertab->wave, aertab->nwave);

        /* grab angstrom coefficients and sort in ascending order */
        for (im = 0; im < aertab->nmodel; im++) {
            alphaT[im].angstrom = aertab->model[im]->angstrom[ib];
            alphaT[im].modnum = im;
        }
        qsort(alphaT, aertab->nmodel, sizeof (alphaTstr),
                (int (*)(const void *, const void *)) compalphaT);

        firstCall = 0;
    }

    for (im = 0; im < aertab->nmodel; im++) {
        if (angstrom < alphaT[im].angstrom)
            break;
    }
    im1 = MAX(MIN(im - 1, aertab->nmodel - 1), 0);
    im2 = MAX(MIN(im, aertab->nmodel - 1), 0);

    *modmin = alphaT[im1].modnum;
    *modmax = alphaT[im2].modnum;

    if (im1 == im2)
        *modrat = 1.0;
    else
        *modrat = (angstrom - alphaT[im1].angstrom) /
        (alphaT[im2].angstrom - alphaT[im1].angstrom);

    return;

}

/* ---------------------------------------------------------------------------------------- */
/* model_taua() - compute AOT at input wavelengths using specified model                    */
/*  Note by W. Robinson - the aot is determined for SENSOR wavelength nir_l
    using the mu, mu0 for the geometry for that band.  To get the aot at other
    TABLE bands, it starts with the aot(nir_l).  So, how would the geometry
    change get properly translated to the other table bands?  For now, I'm
    not accounting for that translation.  Possible method would be to (for
    interpol=0 only) compute the aot(nir_l, geom(sensor band)) and use
    that aot to get the aot(geom(sensor band))                                              */

/* ---------------------------------------------------------------------------------------- */
void model_taua(int32_t sensorID, int modnum, float wave[], int32_t nwave,
        int32_t iwnir_l, float rhoa[], geom_str *geom, float wv, float taua[]) {
    float *aot;
    float *lnaot;
    float *rhoas;

    int iwnir = iwatab[iwnir_l];
    float *phase = model_phase(modnum, geom);
    int iw, iwg, iwtab;
    float maxwave;

    iwg = iwnir * geom->gmult; /* index for geometry */

    maxwave = MAX(aertab->nwave, nwave);

    if ((rhoas = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas.\n");
        exit(1);
    }
    if ((aot = (float *) calloc(maxwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for aot.\n");
        exit(1);
    }
    if ((lnaot = (float *) calloc(maxwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for lnaot.\n");
        exit(1);
    }

    /* get SS aerosol reflectance at longest sensor wavelength */
    rhoa_to_rhoas(sensorID, modnum, geom, wv, rhoa, wave, nwave, iwnir_l, iwnir_l, rhoas);

    /* get aerosol optical thickness at longest sensor wavelength */
    aot[iwnir] = rhoas[iwnir_l]*(4.0 * geom->csolz[iwg] * geom->csenz[iwg])
            / (phase[iwnir] * aertab->model[modnum]->albedo[iwnir]);

    /* get aerosol optical thickness at all other table wavelengths */
    for (iw = 0; iw < aertab->nwave; iw++) {
        /* note to self: i actually have aot(865) for czcs at this point */
        aot[iw] = aot[iwnir] * aertab->model[modnum]->extc[iw] / aertab->model[modnum]->extc[iwnir];
        if (interpol)
            lnaot[iw] = log(aot[iw]);
    }

    /* interpolate aot from table to sensor wavelengths */
    if (interpol) {
        for (iw = 0; iw < nwave; iw++) {
            iwtab = iwatab[iw];
            if (aertab->wave[iwtab] != wave[iw] && wave[iw] > 0)
                taua[iw] = exp(linterp(aertab->wave, lnaot, aertab->nwave, wave[iw]));
            else
                taua[iw] = aot[iwtab];
        }
    } else
        for (iw = 0; iw < nwave; iw++)
            taua[iw] = aot[iw];

    free(rhoas);
    free(aot);
    free(lnaot);
    return;
}


/* ---------------------------------------------------------------------------------------- */
/* model_transmittance() - compute path Rayleigh-aerosol diffuse trans for specified model  */
/*                                                                                          */
/*  W. Robinson, SAIC, 24 Mar 2017, modify for band-dependent geometry                      */

/* ---------------------------------------------------------------------------------------- */
void model_transmittance(int modnum, float wave[], int32_t nwave,
        float *theta, int gmult, float taua[], float dtran[]) {
    static int firstCall = 1;
    static float *intexp;
    static int *inttst;

    int i1, i2, i, ig;
    int iw, iwtab;
    float a1, b1;
    float a2, b2;
    float x1, x2;
    float y1, y2;
    float xbar;
    float wt;

    if (firstCall) {
        float taur1, um1;
        float taur2, um2;
        firstCall = 0;

        intexp = (float *) malloc(nwave * sizeof (float));
        inttst = (int *) malloc(nwave * sizeof (int));

        for (iw = 0; iw < nwave; iw++) {
            intexp[iw] = 1.0;
            inttst[iw] = 0;
            iwtab = iwdtab[iw];
            if (fabs(wave[iw] - aertab->dtran_wave[iwtab]) > 0.51) {
                um1 = aertab->dtran_wave[iwtab] / 1000.0;
                um2 = wave[iw] / 1000.0;
                taur1 = 0.008569 * pow(um1, -4)*(1.0 + (0.0113 * pow(um1, -2))+(0.00013 * pow(um1, -4)));
                taur2 = 0.008569 * pow(um2, -4)*(1.0 + (0.0113 * pow(um2, -2))+(0.00013 * pow(um2, -4)));
                intexp[iw] = taur2 / taur1;
                inttst[iw] = 1;
                printf("Interpolating diffuse transmittance for %d from %f by %f\n",
                        (int) wave[iw], aertab->dtran_wave[iwtab], intexp[iw]);
            }
        }
    }

    /* use coefficients of nearest wavelength */
    for (iw = 0; iw < nwave; iw++) {
        if ((iw == 0) || (gmult != 0)) {
            ig = iw * gmult;
            /* find bracketing zenith angle indices */
            for (i = 0; i < aertab->dtran_ntheta; i++) {
                if (theta[ig] < aertab->dtran_theta[i])
                    break;
            }
            if (i == aertab->dtran_ntheta) {
                i1 = i - 1;
                i2 = i - 1;
                wt = 0.0;
            } else {
                i1 = MIN(MAX(i - 1, 0), aertab->dtran_ntheta - 2);
                i2 = i1 + 1;
                x1 = aertab->dtran_airmass[i1];
                x2 = aertab->dtran_airmass[i2];
                xbar = 1.0 / cos(theta[ig] / radeg);
                wt = (xbar - x1) / (x2 - x1);
            }
        }
        iwtab = iwdtab[iw];

        a1 = aertab->model[modnum]->dtran_a[iwtab][i1];
        b1 = aertab->model[modnum]->dtran_b[iwtab][i1];

        a2 = aertab->model[modnum]->dtran_a[iwtab][i2];
        b2 = aertab->model[modnum]->dtran_b[iwtab][i2];

        if (inttst[iw]) {
            a1 = pow(a1, intexp[iw]);
            a2 = pow(a2, intexp[iw]);
        }

        y1 = a1 * exp(-b1 * taua[iw]);
        y2 = a2 * exp(-b2 * taua[iw]);

        dtran[iw] = MAX(MIN((1.0 - wt) * y1 + wt*y2, 1.0), 1e-5);
    }

    return;
}


/* ---------------------------------------------------------------------------------------- */
/* diff_tran() - compute Rayleigh-aerosol diffuse trans for selected model pair, both paths */

/* ---------------------------------------------------------------------------------------- */
void diff_tran(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_l,
        geom_str *geom, float wv, float pr, float taur[], int32_t modmin,
        int32_t modmax, float modrat, float rhoa[], float taua[], float tsol[],
        float tsen[], float tauamin[], float tauamax[], int tauaopt) {
    int iw, gmult, ig;
    float *tsolmin;
    float *tsolmax;
    float *tsenmin;
    float *tsenmax;

    if ((tsolmin = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsolmin.\n");
        exit(1);
    }
    if ((tsolmax = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsolmax.\n");
        exit(1);
    }
    if ((tsenmin = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsenmin.\n");
        exit(1);
    }
    if ((tsenmax = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsenmax.\n");
        exit(1);
    }
    gmult = geom->gmult;
    if (interpol == 1) gmult = 0; /* geom used for model wavelengths */

    /* get AOT per band for each model */
    if (tauaopt != 1) {
        model_taua(sensorID, modmin, wave, nwave, iwnir_l, rhoa, geom, wv, tauamin);
        model_taua(sensorID, modmax, wave, nwave, iwnir_l, rhoa, geom, wv, tauamax);
    }

    /* get diff trans sun to ground, per band for each model */
    model_transmittance(modmin, wave, nwave, geom->solz, gmult, tauamin, tsolmin);
    model_transmittance(modmax, wave, nwave, geom->solz, gmult, tauamax, tsolmax);

    /* get diff trans ground to sensor, per band for each model */
    model_transmittance(modmin, wave, nwave, geom->senz, gmult, tauamin, tsenmin);
    model_transmittance(modmax, wave, nwave, geom->senz, gmult, tauamax, tsenmax);

    /* interpolate and pressure correct */
    for (iw = 0; iw < nwave; iw++) {
        ig = iw * geom->gmult; /* geom used for sensor wavelengths */
        taua[iw] = tauamin[iw]*(1.0 - modrat) + tauamax[iw] * modrat;
        tsol[iw] = tsolmin[iw]*(1.0 - modrat) + tsolmax[iw] * modrat;
        tsen[iw] = tsenmin[iw]*(1.0 - modrat) + tsenmax[iw] * modrat;

        /* correct for pressure difference from standard pressure */
        tsol[iw] = tsol[iw] * exp(-0.5 * taur[iw] / geom->csolz[ig]*(pr / p0 - 1));
        tsen[iw] = tsen[iw] * exp(-0.5 * taur[iw] / geom->csenz[ig] *(pr / p0 - 1));

        if ((evalmask & TRANSSPHER) != 0) {
            /* correct for airmass difference, plane-parallel to spherical atmosphere */
            tsol[iw] = pow(tsol[iw], geom->airmass_sph[ig] / geom->airmass_plp[ig]);
            tsen[iw] = pow(tsen[iw], geom->airmass_sph[ig] / geom->airmass_plp[ig]);
        }
    }

    free(tsolmin);
    free(tsolmax);
    free(tsenmin);
    free(tsenmax);

    return;
}



/* ---------------------------------------------------------------------------------------- */
/* smaer() - compute aerosol reflectance using MSEPS approach of Ahmad
 * perform spectral matching in reflectance space                                           */
/* output the reflectance of the two bracketing aerosol models                              */
/*
 * Amir Ibrahim, Jul 2016.                                                                  */

/* ---------------------------------------------------------------------------------------- */

int smaer(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s,
        int32_t iwnir_l, int32_t nmodels, int32_t mindx1[], int32_t mindx2[],
        int32_t mindx3[], geom_str *geom, float wv, float rhoa[],
        float rho_aer[], int32_t *modmin, int32_t *modmax, float *modrat,
        float tau_pred_min[], float tau_pred_max[], float tg_sol_sm[],
        float tg_sen_sm[], float Lt_sm[], int32_t ip) {
    int iw;

    if (!have_ms) {
        printf("\nThe multi-scattering spectral matching atmospheric correction method requires\n");
        printf("ams_all, bms_all, cms_all, dms_all, ems in the aerosol model tables.\n");
        exit(1);
    }

    // calculate the weighted aerosol reflectance from Red to SWIR
    if (spectral_matching(sensorID, wave, nwave, iwnir_s, iwnir_l, nmodels, mindx1, mindx2, mindx3, geom, wv, rhoa,
            rho_aer, modmin, modmax, modrat, tau_pred_min, tau_pred_max, tg_sol_sm, tg_sen_sm, Lt_sm, ip) == 0) {

        // extrapolate to VIS
        for (iw = 0; iw < nwave; iw++) {
            rhoa[iw] = rho_aer[iw];
        }
    } else
        return (1);

    return (0);
}

/* ---------------------------------------------------------------------------------------- */
/* ahmadaer() - compute aerosol reflectance using MSEPS approach of Ahmad                   */
/*                                                                                          */
/* Z. Ahmad, August 2014.                                                                   */

/* ---------------------------------------------------------------------------------------- */
int ahmadaer(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
        int32_t nmodels, int32_t mindx[],
        geom_str *geom, float wv, float rhoa[],
        int32_t *modmin, int32_t *modmax, float *modrat, float *epsnir,
        float tau_pred_min[], float tau_pred_max[]) {
    int iw;

    float rho_pred_min[nwave], rho_pred_max[nwave];
    float rho_aer[nwave], tau_aer[nwave];

    if (!have_ms) {
        printf("\nThe multi-scattering epsilon atmospheric correction method requires\n");
        printf("ams_all, bms_all, cms_all, dms_all, ems in the aerosol model tables.\n");
        exit(1);
    }
    if (aer_opt == AERRHMSEPS_lin) {
        if (ahmad_atm_corr_lin(sensorID, wave, nwave, iwnir_s, iwnir_l, nmodels, mindx, geom, wv,
                rhoa, modmin, modmax, modrat, epsnir,
                tau_pred_max, tau_pred_min, rho_pred_max, rho_pred_min, tau_aer, rho_aer) != 0)
            return (1);
    } else {
        /* use the ms_epsilon method to get rhoa */
        if (ahmad_atm_corr(sensorID, wave, nwave, iwnir_s, iwnir_l, nmodels, mindx, geom, wv,
                rhoa, modmin, modmax, modrat, epsnir,
                tau_pred_max, tau_pred_min, rho_pred_max, rho_pred_min, tau_aer, rho_aer) != 0)
            return (1);
    }

    for (iw = 0; iw < nwave; iw++) {
        rhoa[iw] = rho_aer[iw];
    }

    return (0);
}


/* ---------------------------------------------------------------------------------------- */
/* wangaer() - compute aerosol reflectance using Gordon & Wang 1994 algorithm               */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */

int wangaer(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s,
        int32_t iwnir_l, int32_t nmodels, int32_t mindx[], geom_str *geom,
        float wv, float rhoa[], int32_t *modmin, int32_t *modmax,
        float *modrat, float *epsnir, float tauamin[], float tauamax[]) {
    int modflg;
    float *epsmin;
    float *epsmax;

    float epsmin1;
    float epsmax1;

    float *rhoasmin;
    float *rhoasmax;
    float *rhoamin;
    float *rhoamax;

    float cc = 0.0;
    int iw;

    if ((rhoasmin = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoasmin.\n");
        exit(1);
    }
    if ((rhoasmax = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoasmax.\n");
        exit(1);
    }
    if ((rhoamin = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoamin.\n");
        exit(1);
    }
    if ((rhoamax = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoasmax.\n");
        exit(1);
    }

    /* find upper and lower-bounding models */
    modflg = model_select_wang(sensorID, wave, nwave, nmodels, mindx, geom, wv,
            rhoa, iwnir_s, iwnir_l, modmin, modmax, modrat, epsnir);

    /* if no lower-bounding aerosol model, set-up for extrapolation */
    if (modflg < 0)
        cc = log(*epsnir) / (wave[iwnir_l] - wave[iwnir_s]);

    /* get model epsilon for each bounding model, all wavelengths */
    epsmin = model_epsilon(*modmin, iwnir_l, wave, nwave, geom);
    epsmax = model_epsilon(*modmax, iwnir_l, wave, nwave, geom);

    /* get SS aerosol reflectance at longest wavelength for the two models */
    if (rhoa_to_rhoas(sensorID, *modmin, geom, wv, rhoa, wave, nwave, iwnir_l, iwnir_l, rhoasmin) != 0) {
        free(rhoamin);
        free(rhoasmin);
        free(rhoamax);
        free(rhoasmax);
        return (1);
    }
    if (rhoa_to_rhoas(sensorID, *modmax, geom, wv, rhoa, wave, nwave, iwnir_l, iwnir_l, rhoasmax) != 0) {
        free(rhoamin);
        free(rhoasmin);
        free(rhoamax);
        free(rhoasmax);
        return (1);
    }
    /* compute SS aerosol reflectance in all bands */
    for (iw = 0; iw < nwave; iw++) {

        epsmin1 = epsmin[iw];
        epsmax1 = epsmax[iw];

        if (modflg < 0) {
            epsmin1 = exp(cc * (wave[iwnir_l] - wave[iw]));
            epsmax1 = epsmin1;
        }

        rhoasmin[iw] = rhoasmin[iwnir_l] * epsmin1;
        rhoasmax[iw] = rhoasmax[iwnir_l] * epsmax1;
    }

    /* compute MS aerosol reflectance in visible bands */
    rhoas_to_rhoa(sensorID, *modmin, geom, wv, rhoasmin, wave, nwave, 0, nwave - 1, rhoamin);
    rhoas_to_rhoa(sensorID, *modmax, geom, wv, rhoasmax, wave, nwave, 0, nwave - 1, rhoamax);

    /* interpolate between upper and lower-bounding models */
    for (iw = 0; iw < nwave; iw++) {
        rhoa[iw] = rhoamin[iw]*(1.0 - (*modrat)) + rhoamax[iw]*(*modrat);
    }

    model_taua(sensorID, *modmin, wave, nwave, iwnir_l, rhoa, geom, wv, tauamin);
    model_taua(sensorID, *modmax, wave, nwave, iwnir_l, rhoa, geom, wv, tauamax);

    free(rhoamin);
    free(rhoasmin);
    free(rhoamax);
    free(rhoasmax);

    return (0);
}

/* ---------------------------------------------------------------------------------------- */
/* model_select_franz() - Franz aerosol model selection process.                            */
/*                                                                                          */
/* B. Franz, 1 February 2009.                                                               */

/* ---------------------------------------------------------------------------------------- */

typedef struct rhoaT_struct {
    int32_t modnum;
    float rhoa;
    float eps;
} rhoaTstr;

int comp_rhoaT(rhoaTstr *x, rhoaTstr *y) {
    return (x->rhoa < y->rhoa ? -1 : 1);
}

int model_select_franz(int32_t sensorID, float wave[], int32_t nwave,
        int32_t nmodel, int32_t mindx[], geom_str *geom, float wv,
        float rhoa[], int32_t iwnir_s, int32_t iwnir_l, int32_t *modmin,
        int32_t *modmax, float *modrat, float *epsnir) {

    float *rhoas;
    float *rhoa_tmp;
    rhoaTstr rhoa_tab[MAXMODEL];

    float *eps;
    int jm, im;
    int jm1, jm2;
    float wt;

    *modmin = -1;
    *modmax = -1;
    *modrat = 0.0;
    *epsnir = 0.0;

    if ((rhoas = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas.\n");
        exit(1);
    }
    if ((rhoa_tmp = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas_tmp.\n");
        exit(1);
    }

    // predict MS aerosol reflectance assuming each model

    for (jm = 0; jm < nmodel; jm++) {

        im = mindx[jm];

        // get model epsilon at this geometry
        eps = model_epsilon(im, iwnir_l, wave, nwave, geom);

        // get SS aerosol reflectance at iwnir_l
        if (rhoa_to_rhoas(sensorID, im, geom, wv, rhoa, wave, nwave, iwnir_l, iwnir_l, rhoas) != 0) {
            free(rhoas);
            free(rhoa_tmp);
            return (1);
        }
        // get SS aerosol reflectance at iwnir_s
        rhoas[iwnir_s] = rhoas[iwnir_l] * eps[iwnir_s];

        // get MS aerosol reflectance at iwnir_s and save
        rhoas_to_rhoa(sensorID, im, geom, wv, rhoas, wave, nwave, iwnir_s, iwnir_s, rhoa_tmp);

        rhoa_tab[jm].modnum = im;
        rhoa_tab[jm].rhoa = rhoa_tmp[iwnir_s];
        rhoa_tab[jm].eps = eps[iwnir_s];
    }

    // put results in ascending order of predicted rhoa[iwnir_s]
    qsort(rhoa_tab, nmodel, sizeof (rhoaTstr), (int (*)(const void *, const void *)) comp_rhoaT);

    // compare observed rhoa with model predictions at iwnir_s to select models
    for (jm = 0; jm < nmodel; jm++) {
        if (rhoa_tab[jm].rhoa > rhoa[iwnir_s])
            break;
    }
    if (jm == 0) {
        jm1 = 0;
        jm2 = 1;
    } else if (jm == nmodel) {
        jm1 = nmodel - 2;
        jm2 = nmodel - 1;
    } else {
        jm1 = jm - 1;
        jm2 = jm1 + 1;
    }
    wt = (rhoa[iwnir_s] - rhoa_tab[jm1].rhoa) / (rhoa_tab[jm2].rhoa - rhoa_tab[jm1].rhoa);

    *modmin = rhoa_tab[jm1].modnum;
    *modmax = rhoa_tab[jm2].modnum;
    *modrat = wt;
    *epsnir = rhoa_tab[jm1].eps * (1.0 - wt) + rhoa_tab[jm2].eps*wt;

    free(rhoas);
    free(rhoa_tmp);

    return (0);
}

/* ---------------------------------------------------------------------------------------- */
/* franzaer() - compute aerosol reflectance using modified G&W algorithm                    */
/*                                                                                          */
/* B. Franz, February 2009.                                                                 */

/* ---------------------------------------------------------------------------------------- */

int franzaer(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s,
        int32_t iwnir_l, int32_t nmodels, int32_t mindx[], geom_str *geom,
        float wv, float rhoa[], int32_t *modmin, int32_t *modmax,
        float *modrat, float *epsnir, float tauamin[], float tauamax[]) {
    float *epsmin;
    float *epsmax;

    float epsmin1;
    float epsmax1;

    float *rhoasmin;
    float *rhoasmax;
    float *rhoamin;
    float *rhoamax;

    int iw;

    if ((rhoasmin = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoasmin.\n");
        exit(1);
    }
    if ((rhoasmax = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoasmax.\n");
        exit(1);
    }
    if ((rhoamin = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoamin.\n");
        exit(1);
    }
    if ((rhoamax = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoamax.\n");
        exit(1);
    }


    /* find upper and lower-bounding models */
    if (model_select_franz(sensorID, wave, nwave, nmodels, mindx, geom, wv, rhoa, iwnir_s, iwnir_l,
            modmin, modmax, modrat, epsnir) != 0)
        return (1);

    /* get model epsilon for each bounding model, all wavelengths */
    epsmin = model_epsilon(*modmin, iwnir_l, wave, nwave, geom);
    epsmax = model_epsilon(*modmax, iwnir_l, wave, nwave, geom);

    /* get SS aerosol reflectance at longest wavelength for the two models */
    if (rhoa_to_rhoas(sensorID, *modmin, geom, wv, rhoa, wave, nwave, iwnir_l, iwnir_l, rhoasmin) != 0) {
        free(rhoasmin);
        free(rhoasmax);
        free(rhoamin);
        free(rhoamax);
        return (1);
    }
    if (rhoa_to_rhoas(sensorID, *modmax, geom, wv, rhoa, wave, nwave, iwnir_l, iwnir_l, rhoasmax) != 0) {
        free(rhoasmin);
        free(rhoasmax);
        free(rhoamin);
        free(rhoamax);
        return (1);
    }
    /* compute SS aerosol reflectance in all bands */
    for (iw = 0; iw < nwave; iw++) {

        epsmin1 = epsmin[iw];
        epsmax1 = epsmax[iw];

        rhoasmin[iw] = rhoasmin[iwnir_l] * epsmin1;
        rhoasmax[iw] = rhoasmax[iwnir_l] * epsmax1;
    }

    /* compute MS aerosol reflectance in visible bands */
    rhoas_to_rhoa(sensorID, *modmin, geom, wv, rhoasmin, wave, nwave, 0, nwave - 1, rhoamin);
    rhoas_to_rhoa(sensorID, *modmax, geom, wv, rhoasmax, wave, nwave, 0, nwave - 1, rhoamax);

    /* interpolate between upper and lower-bounding models */
    for (iw = 0; iw < nwave; iw++) {
        rhoa[iw] = rhoamin[iw]*(1.0 - (*modrat)) + rhoamax[iw]*(*modrat);
    }

    model_taua(sensorID, *modmin, wave, nwave, iwnir_l, rhoa, geom, wv, tauamin);
    model_taua(sensorID, *modmax, wave, nwave, iwnir_l, rhoa, geom, wv, tauamax);

    free(rhoasmin);
    free(rhoasmax);
    free(rhoamin);
    free(rhoamax);

    return (0);
}

/* ---------------------------------------------------------------------------------------- */
/* order_models is the sorting function used in rhaer                                       */

/* ---------------------------------------------------------------------------------------- */
static int order_models(const void *p1, const void *p2) {
    aermodstr *x = *(aermodstr **) p1;
    aermodstr *y = *(aermodstr **) p2;

    if (x->rh == y->rh) {
        if (x->sd > y->sd)
            return ( 1);
        else
            return (-1);
    } else {
        if (x->rh > y->rh)
            return ( 1);
        else
            return (-1);
    }
}


/* ---------------------------------------------------------------------------------------- */
/* rhaer() - compute aerosol reflectance using RH descrimination + desired selection scheme */

/* ---------------------------------------------------------------------------------------- */
int rhaer(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
        geom_str *geom, float wv, float rh, float pr, float taur[], float rhoa[],
        int32_t *modmin1, int32_t *modmax1, float *modrat1, int32_t *modmin2, int32_t *modmax2, float *modrat2,
        float *eps, float taua[], float tsol[], float tsen[], float tg_sol_sm[], float tg_sen_sm[], float Lt_sm[], int32_t ip) {
    static int firstCall = 1;
    static int nrh;
    static float rhtab[MAXAERMOD];
    static int nsd;
    static int sdtab[MAXAERMOD];

    float *rhoa1;
    float *rhoa2;
    float *taua1;
    float *taua2;
    float *tsol1;
    float *tsol2;
    float *tsen1;
    float *tsen2;
    float eps1;
    float eps2;
    //float modrat1;
    //float modrat2;

    float *tau_pred_min1;
    float *tau_pred_max1;
    float *tau_pred_min2;
    float *tau_pred_max2;

    //int32_t  modmin1;
    //int32_t  modmin2;
    //int32_t  modmax1;
    //int32_t  modmax2;
    int32_t mindx1[MAXAERMOD];
    int32_t mindx2[MAXAERMOD];
    int32_t mindx3[MAXAERMOD]; // Third model set defined --> Amir
    int irh1, irh2, irh;
    int irh3; // Third RH index --> Amir
    int isd;
    float wt;


    int iw, im;

    if (firstCall) {
        firstCall = 0;
        float lastrh = -1.0;
        int lastsd = -1;
        if (!have_rh || !have_sd) {
            printf("-E- %s line %d: This aerosol selection method requires models with a Relative Humidity attribute.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        // need in order of rh and sd within rh
        qsort(aertab->model, aertab->nmodel, sizeof (aermodstr*), (int (*)(const void *, const void *)) order_models);

        // count the number of model humidities and the number of model size distributions
        // note that use of a single model suite will yield nrh=1, which inherently avoids RH weighting that case

        nsd = 0;
        nrh = 0;

        for (im = 0; im < aertab->nmodel; im++) {
            if (aertab->model[im]->rh != lastrh) {
                rhtab[nrh] = aertab->model[im]->rh;
                lastrh = rhtab[nrh];
                nrh++;
            }
            if (nrh == 1 && aertab->model[im]->sd != lastsd) {
                sdtab[nsd] = aertab->model[im]->sd;
                lastsd = sdtab[nsd];
                nsd++;
            }
        }
        if (nrh * nsd != aertab->nmodel) {
            printf("-E- %s line %d: number of humidities (%d) x number of size distributions (%d) must equal number of models (%d).\n",
                    __FILE__, __LINE__, nrh, nsd, aertab->nmodel);
            exit(1);
        } else {
            printf("%d aerosol models: %d humidities x %d size fractions\n", aertab->nmodel, nrh, nsd);
            for (irh = 0; irh < nrh; irh++) {
                for (isd = 0; isd < nsd; isd++) {
                    im = irh * nsd + isd;
                    printf("model %d, rh=%f, sd=%d, alpha=%f, name=%s\n",
                            im, aertab->model[im]->rh, aertab->model[im]->sd, aertab->model[im]->angstrom[1], aertab->model[im]->name);
                }
            }
        }
    }

    // initialize
    if ((taua1 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for taua1.\n");
        exit(1);
    }
    if ((taua2 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for taua2.\n");
        exit(1);
    }
    if ((tsol1 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsol1.\n");
        exit(1);
    }
    if ((tsol2 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsol2.\n");
        exit(1);
    }
    if ((rhoa1 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa1.\n");
        exit(1);
    }
    if ((rhoa2 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa2.\n");
        exit(1);
    }
    if ((tsen1 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsen1.\n");
        exit(1);
    }
    if ((tsen2 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsen2.\n");
        exit(1);
    }
    if ((tau_pred_min1 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tau_pred_min1.\n");
        exit(1);
    }
    if ((tau_pred_min2 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tau_pred_min2.\n");
        exit(1);
    }
    if ((tau_pred_max1 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tau_pred_max1.\n");
        exit(1);
    }
    if ((tau_pred_max2 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tau_pred_max2.\n");
        exit(1);
    }

    //*modmin  = -1;
    //*modmax  = -1;
    //*modrat  =  0;
    //*modmin2 = -1;
    //*modmax2 = -1;
    //*modrat2 =  0;

    for (iw = 0; iw < nwave; iw++) {
        taua[iw] = -1.0;
        tsol[iw] = -1.0;
        tsen[iw] = -1.0;
        rhoa1[iw] = rhoa[iw];
        rhoa2[iw] = rhoa[iw];
        rhoa [iw] = BAD_FLT;
    }


    // adjust rh for spectral matching
    if (aer_opt == AERRHSM) {
        if (rh >= 95) {
            printf("Warning rh is greater than 95%%. Reset to 94%% rh=%f\n", rh);
            rh = 94;
        }
    }

    // find RH index and wts
    if (nrh == 1 || rhtab[0] > rh) { // actual RH < smallest model RH or only one model RH
        irh1 = 0;
        irh2 = 0;
        wt = 0.0;
    } else if (rhtab[nrh - 1] < rh) { // actual RH > smallestlargest model RH
        irh1 = nrh - 1;
        irh2 = nrh - 1;
        wt = 0.0;
    } else {
        for (irh = 0; irh < nrh; irh++) {
            if (rhtab[irh] > rh)
                break;
        }
        irh1 = MIN(MAX(0, irh - 1), nrh - 2);
        irh2 = irh1 + 1;
        wt = (rh - rhtab[irh1]) / (rhtab[irh2] - rhtab[irh1]);
    }

    // calculate the third closest RH index for the models suite -->Amir
    if (rh - rhtab[irh1] <= (rhtab[irh2] - rhtab[irh1]) / 2)
        if (irh1 > 1)
            irh3 = irh1 - 1;
        else
            irh3 = irh1;
    else
        if (irh2 < 10)
        irh3 = irh2 + 1;
    else
        irh3 = irh2;

    for (im = 0; im < nsd; im++) {
        mindx3[im] = irh3 * nsd + im;
    }

    // set indices of active model sets

    for (im = 0; im < nsd; im++) {
        mindx1[im] = irh1 * nsd + im;
        mindx2[im] = irh2 * nsd + im;
    }

    // compute aerosol reflectances, aot, diffuse trans, eps from first model set

    /* perform spectral matching from Red to SWIR in radiance space, based on Ahmad
     * Multi-scattering coeff method --> Amir*/
    if (aer_opt == AERRHSM) {

        int nmodels = 30;

        if (mindx3[0] == mindx1[0] || mindx3[0] == mindx2[0] || mindx3[0] > 79 || mindx3[0] < 0)
            nmodels = 20;

        else if (mindx1[0] == mindx2[0])
            nmodels = 10;

        int32_t vcal_flag = 1;
        if ((aertab->nmodel) == vcal_flag)
            nmodels = 1;

        int modmin, modmax;
        float modrat;

        float rho_aer[nwave];
        //      printf("%zu %zu %zu\n",mindx1[0],mindx2[0],mindx3[0]);

        if (smaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nmodels, mindx1, mindx1, mindx1,
                geom, wv, rhoa1, rho_aer, &modmin, &modmax, &modrat, tau_pred_min1, tau_pred_max1, tg_sol_sm, tg_sen_sm, Lt_sm, ip) == 0) {

            *modmin1 = modmin;
            *modmax1 = modmax;
            *modrat1 = modrat;
            *modmin2 = modmin;
            *modmax2 = modmax;
            *modrat2 = modrat;
            //printf("modmin=%f\n",modrat);
            //tau_pred_min2 = tau_pred_min1;
            //tau_pred_max2 = tau_pred_max1;
            diff_tran(sensorID, wave, nwave, iwnir_l, geom, wv, pr, taur,
                    *modmin1, *modmax1, *modrat1, rhoa1, taua1, tsol1, tsen1, tau_pred_min1, tau_pred_max1, 1);
            //return (0);
        } //second run
        if (smaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nmodels, mindx2, mindx2, mindx2,
                geom, wv, rhoa2, rho_aer, &modmin, &modmax, &modrat, tau_pred_min2, tau_pred_max2, tg_sol_sm, tg_sen_sm, Lt_sm, ip) == 0) {

            *modmin1 = modmin;
            *modmax1 = modmax;
            *modrat1 = modrat;
            *modmin2 = modmin;
            *modmax2 = modmax;
            *modrat2 = modrat;
            //printf("modmin=%f\n",modrat);
            //tau_pred_min2 = tau_pred_min1;
            //tau_pred_max2 = tau_pred_max1;
            diff_tran(sensorID, wave, nwave, iwnir_l, geom, wv, pr, taur,
                    *modmin2, *modmax2, *modrat2, rhoa2, taua2, tsol2, tsen2, tau_pred_min2, tau_pred_max2, 1);
            //return (0);
        } else {
            free(taua1);
            free(taua2);
            free(tsol1);
            free(tsol2);
            free(tsen1);
            free(tsen2);
            free(rhoa1);
            free(rhoa2);
            free(tau_pred_min1);
            free(tau_pred_max1);
            free(tau_pred_min2);
            free(tau_pred_max2);
            return (1);
        }
    } else if (aer_opt == AERRHMSEPS || aer_opt == AERRHMSEPS_lin) {
        if (ahmadaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nsd, mindx1,
                geom, wv, rhoa1, modmin1, modmax1, modrat1, &eps1, tau_pred_min1, tau_pred_max1) != 0) {
            free(taua1);
            free(taua2);
            free(tsol1);
            free(tsol2);
            free(tsen1);
            free(tsen2);
            free(rhoa1);
            free(rhoa2);
            free(tau_pred_min1);
            free(tau_pred_max1);
            free(tau_pred_min2);
            free(tau_pred_max2);
            return (1);

        }
    } else if (aer_opt == AERRHFRNIR) {
        if (franzaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nsd, mindx1,
                geom, wv, rhoa1, modmin1, modmax1, modrat1, &eps1, tau_pred_min1, tau_pred_max1) != 0) {
            free(taua1);
            free(taua2);
            free(tsol1);
            free(tsol2);
            free(tsen1);
            free(tsen2);
            free(rhoa1);
            free(rhoa2);
            free(tau_pred_min1);
            free(tau_pred_max1);
            free(tau_pred_min2);
            free(tau_pred_max2);
            return (1);
        }
    } else {
        if (wangaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nsd, mindx1,
                geom, wv, rhoa1, modmin1, modmax1, modrat1, &eps1, tau_pred_min1, tau_pred_max1) != 0) {
            free(taua1);
            free(taua2);
            free(tsol1);
            free(tsol2);
            free(tsen1);
            free(tsen2);
            free(rhoa1);
            free(rhoa2);
            free(tau_pred_min1);
            free(tau_pred_max1);
            free(tau_pred_min2);
            free(tau_pred_max2);
            return (1);
        }
    }

    if (aer_opt != AERRHSM)
        diff_tran(sensorID, wave, nwave, iwnir_l, geom, wv, pr, taur,
            *modmin1, *modmax1, *modrat1, rhoa1, taua1, tsol1, tsen1, tau_pred_min1, tau_pred_max1, 1);

    // compute aerosol reflectances, aot, diffuse trans, eps from second model set (if needed)

    if (irh2 != irh1) {

        if (aer_opt == AERRHSM) {
            for (iw = 0; iw < nwave; iw++) {
                rhoa[iw] = rhoa1[iw]*(1 - wt) + rhoa2[iw] * wt;
                taua[iw] = taua1[iw]*(1 - wt) + taua2[iw] * wt;
                tsol[iw] = tsol1[iw]*(1 - wt) + tsol2[iw] * wt;
                tsen[iw] = tsen1[iw]*(1 - wt) + tsen2[iw] * wt;
            }
            *eps = eps1;
            return (1); // A second model set interpolation is not needed if the spectral matching is active
        }

        if (aer_opt == AERRHMSEPS || aer_opt == AERRHMSEPS_lin) {
            if (ahmadaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nsd, mindx2,
                    geom, wv, rhoa2, modmin2, modmax2, modrat2, &eps2, tau_pred_min2, tau_pred_max2) != 0) {
                free(taua1);
                free(taua2);
                free(tsol1);
                free(tsol2);
                free(tsen1);
                free(tsen2);
                free(rhoa1);
                free(rhoa2);
                free(tau_pred_min1);
                free(tau_pred_max1);
                free(tau_pred_min2);
                free(tau_pred_max2);
                return (1);
            }
        } else if (aer_opt == AERRHFRNIR) {
            if (franzaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nsd, mindx2,
                    geom, wv, rhoa2, modmin2, modmax2, modrat2, &eps2, tau_pred_min2, tau_pred_max2) != 0) {
                free(taua1);
                free(taua2);
                free(tsol1);
                free(tsol2);
                free(tsen1);
                free(tsen2);
                free(rhoa1);
                free(rhoa2);
                free(tau_pred_min1);
                free(tau_pred_max1);
                free(tau_pred_min2);
                free(tau_pred_max2);
                return (1);
            }
        } else {
            if (wangaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nsd, mindx2,
                    geom, wv, rhoa2, modmin2, modmax2, modrat2, &eps2, tau_pred_min2, tau_pred_max2) != 0) {
                free(taua1);
                free(taua2);
                free(tsol1);
                free(tsol2);
                free(tsen1);
                free(tsen2);
                free(rhoa1);
                free(rhoa2);
                free(tau_pred_min1);
                free(tau_pred_max1);
                free(tau_pred_min2);
                free(tau_pred_max2);
                return (1);
            }
        }

        diff_tran(sensorID, wave, nwave, iwnir_l, geom, wv, pr, taur,
                *modmin2, *modmax2, *modrat2, rhoa2, taua2, tsol2, tsen2, tau_pred_min2, tau_pred_max2, 1);

        for (iw = 0; iw < nwave; iw++) {
            rhoa[iw] = rhoa1[iw]*(1 - wt) + rhoa2[iw] * wt;
            taua[iw] = taua1[iw]*(1 - wt) + taua2[iw] * wt;
            tsol[iw] = tsol1[iw]*(1 - wt) + tsol2[iw] * wt;
            tsen[iw] = tsen1[iw]*(1 - wt) + tsen2[iw] * wt;
        }
        *eps = eps1 * (1 - wt) + eps2*wt;

    } else {

        for (iw = 0; iw < nwave; iw++) {
            rhoa[iw] = rhoa1[iw];
            taua[iw] = taua1[iw];
            tsol[iw] = tsol1[iw];
            tsen[iw] = tsen1[iw];
        }
        *eps = eps1;
    }

    // store model info for the dominant RH only

    //if (wt < 0.5) {
    //    *modmin = modmin1;
    //    *modmax = modmax1;
    //    *modrat = modrat1;
    //} else {
    //    *modmin = modmin2;
    //    *modmax = modmax2;
    //    *modrat = modrat2;
    //}

    free(taua1);
    free(taua2);
    free(tsol1);
    free(tsol2);
    free(tsen1);
    free(tsen2);
    free(rhoa1);
    free(rhoa2);
    free(tau_pred_min1);
    free(tau_pred_max1);
    free(tau_pred_min2);
    free(tau_pred_max2);

    return (0);
}


/* ---------------------------------------------------------------------------------------- */
/* fixedaer() - compute aerosol reflectance for fixed aerosol model                         */
/*                                                                                          */
/* B. Franz, August 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
int fixedaer(int32_t sensorID, int32_t modnum, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
        char models[MAXAERMOD][32], int32_t nmodels,
        geom_str *geom, float wv, float rhoa[], float *epsnir) {
    float *eps;
    float *rhoas;
    int iw;

    if ((rhoas = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas.\n");
        exit(1);
    }

    if (rhoa[iwnir_l] < 0.0) {
        *epsnir = BAD_FLT;
        //        for (iw=0; iw<nwave; iw++)
        //            rhoas[iw] = BAD_FLT;
        free(rhoas);
        return (1);
    }

    /* get model epsilon for all wavelengths at this geometry */
    eps = model_epsilon(modnum, iwnir_l, wave, nwave, geom);

    /* get SS aerosol reflectance at longest wavelength */
    if (rhoa_to_rhoas(sensorID, modnum, geom, wv, rhoa, wave, nwave, iwnir_l, iwnir_l, rhoas) != 0) {
        printf("Error getting rhoas\n");
        free(rhoas);
        return (1);
    }

    /* compute SS aerosol reflectance in visible bands */
    for (iw = 0; iw < nwave; iw++) {
        rhoas[iw] = rhoas[iwnir_l] * eps[iw];
    }

    /* compute MS aerosol reflectance in visible bands */
    rhoas_to_rhoa(sensorID, modnum, geom, wv, rhoas, wave, nwave, 0, nwave - 1, rhoa);

    if (iwnir_s == iwnir_l)
        *epsnir = eps[iwnir_l - 1];
    else
        *epsnir = eps[iwnir_s];

    free(rhoas);

    return (0);
}


/* ---------------------------------------------------------------------------------------- */
/* fixedmodpair() - compute aerosol reflectance for fixed model pair                        */
/*                                                                                          */
/* B. Franz, August 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
int fixedmodpair(int32_t sensorID, float wave[], int32_t nwave,
        int32_t iwnir_s, int32_t iwnir_l, geom_str *geom, float wv,
        int32_t modmin, int32_t modmax, float modrat, float rhoa[], float *eps) {
    float *rhoa1;
    float *rhoa2;
    float eps1;
    float eps2;
    int iw;
    int status;

    if ((rhoa1 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa1.\n");
        exit(1);
    }
    if ((rhoa2 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa2.\n");
        exit(1);
    }

    if (modmin < 0 || modmin >= input->naermodels ||
            modmax < 0 || modmax >= input->naermodels ||
            modrat < 0.0 || modrat > 1.0) {
        printf("Bogus input for fixed model pair: %d %d %f\n", modmin + 1, modmax + 1, modrat);
        exit(1);
    }


    if (rhoa[iwnir_l] > input->rhoamin) {

        rhoa2[iwnir_l] = rhoa1[iwnir_l] = rhoa[iwnir_l];

        /* get aerosol reflectance in visible for fixed model, and epsilon */
        status = fixedaer(sensorID, modmin, wave, nwave, iwnir_s, iwnir_l,
                input->aermodels, input->naermodels, geom, wv, rhoa1, &eps1);
        status = fixedaer(sensorID, modmax, wave, nwave, iwnir_s, iwnir_l,
                input->aermodels, input->naermodels, geom, wv, rhoa2, &eps2);

        /* convert aerosol relectance to radiance */
        if (status == 0) {
            for (iw = 0; iw < nwave; iw++) {
                if (iw != iwnir_l) // without this check tLw-La may go slight negative
                    rhoa[iw] = MAX((1.0 - modrat) * rhoa1[iw] + modrat * rhoa2[iw], 0.0);
            }
            *eps = (1.0 - modrat) * eps1 + modrat*eps2;
        }

    } else if (rhoa[iwnir_l] > -(input->rhoamin)) {

        /* if input NIR is near zero, assume a small white aerosol signal */
        *eps = 1.0;
        for (iw = 0; iw < nwave; iw++) {
            rhoa[iw] = MAX(rhoa[iwnir_l], 1e-6);
        }

        status = 0;

    } else {

        /* if input NIR is very negative, fail the case */
        status = 1;
    }

    free(rhoa1);
    free(rhoa2);

    return (status);
}


/* ---------------------------------------------------------------------------------------- */
/* fixedaot() - compute aerosol reflectance for fixed aot(lambda)                           */
/*                                                                                          */
/* B. Franz, August 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */

int fixedaot(int32_t sensorID, float aot[], float wave[], int32_t nwave,
        int32_t iwnir_s, int32_t iwnir_l, geom_str *geom, float wv,
        int32_t *modmin, int32_t *modmax, float *modrat, float rhoa[],
        float *epsnir) {
    static int firstCall = 1;
    static int angst_band1 = -1;
    static int angst_band2 = -1;

    float *phase1;
    float *phase2;
    float *f1;
    float *f2;
    float *lnf1;
    float *lnf2;
    float *rhoas1;
    float *rhoas2;
    float *rhoa1;
    float *rhoa2;
    float eps1;
    float eps2;
    float angstrom;
    int ig, gmult, iw, iwtab;
    float maxwave;

    maxwave = MAX(aertab->nwave, nwave);

    if ((rhoa1 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa1.\n");
        exit(1);
    }
    if ((rhoa2 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa2.\n");
        exit(1);
    }
    if ((rhoas1 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas1.\n");
        exit(1);
    }
    if ((rhoas2 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas2.\n");
        exit(1);
    }
    if ((f1 = (float *) calloc(maxwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for f1.\n");
        exit(1);
    }
    if ((f2 = (float *) calloc(maxwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for f2.\n");
        exit(1);
    }
    if ((lnf1 = (float *) calloc(maxwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for lnf1.\n");
        exit(1);
    }
    if ((lnf2 = (float *) calloc(maxwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for lnf2.\n");
        exit(1);
    }

    if (firstCall) {
        angst_band1 = windex(520, wave, nwave);
        angst_band2 = windex(865, wave, nwave);
        firstCall = 0;
    }

    /* bail on negative input aot */
    for (iw = 0; iw < nwave; iw++) {
        if (aot[iw] < 0.0) {
            free(rhoa1);
            free(rhoa2);
            free(rhoas1);
            free(rhoas2);
            free(f1);
            free(f2);
            free(lnf1);
            free(lnf2);
            return (1);
        }
    }

    /* compute angstrom and use to select bounding models */
    if (aot[iwnir_l] > 0.0)
        angstrom = -log(aot[angst_band1] / aot[angst_band2]) /
        log(wave[angst_band1] / wave[angst_band2]);
    else
        angstrom = 0.0;

    model_select_angstrom(angstrom, modmin, modmax, modrat);


    /* get model phase function for all wavelengths at this geometry for the two models */
    phase1 = model_phase(*modmin, geom);
    phase2 = model_phase(*modmax, geom);

    gmult = (interpol == 1) ? 0 : geom->gmult;

    /* compute factor for SS approximation, set-up for interpolation */
    for (iw = 0; iw < aertab->nwave; iw++) {
        ig = gmult * iw;
        f1[iw] = aertab->model[*modmin]->albedo[iw] *
                phase1[iw] / 4.0 / geom->csolz[ig] / geom->csenz[ig];
        f2[iw] = aertab->model[*modmax]->albedo[iw] *
                phase2[iw] / 4.0 / geom->csolz[ig] / geom->csenz[ig];
        if (interpol) {
            lnf1[iw] = log(f1[iw]);
            lnf2[iw] = log(f2[iw]);
        }
    }

    /* compute SS aerosol reflectance */
    if (interpol) {
        for (iw = 0; iw < nwave; iw++) {
            iwtab = iwatab[iw];
            if (aertab->wave[iwtab] != wave[iw] && wave[iw] > 0) {
                rhoas1[iw] = aot[iw] * exp(linterp(aertab->wave, lnf1, aertab->nwave, wave[iw]));
                rhoas2[iw] = aot[iw] * exp(linterp(aertab->wave, lnf2, aertab->nwave, wave[iw]));
            } else {
                rhoas1[iw] = aot[iw] * f1[iwtab];
                rhoas2[iw] = aot[iw] * f2[iwtab];
            }
        }
    } else {
        for (iw = 0; iw < nwave; iw++) {
            iwtab = iwatab[iw];
            rhoas1[iw] = aot[iw] * f1[iwtab];
            rhoas2[iw] = aot[iw] * f2[iwtab];
        }
    }
    eps1 = rhoas1[iwnir_s] / rhoas1[iwnir_l];
    eps2 = rhoas2[iwnir_s] / rhoas2[iwnir_l];

    /* compute MS aerosol reflectance */
    rhoas_to_rhoa(sensorID, *modmin, geom, wv, rhoas1, wave, nwave, 0, nwave - 1, rhoa1);
    rhoas_to_rhoa(sensorID, *modmax, geom, wv, rhoas2, wave, nwave, 0, nwave - 1, rhoa2);

    for (iw = 0; iw < nwave; iw++) {
        rhoa[iw] = (1.0 - *modrat) * rhoa1[iw] + *modrat * rhoa2[iw];
    }
    *epsnir = (1.0 - *modrat) * eps1 + *modrat * eps2;

    free(rhoa1);
    free(rhoa2);
    free(rhoas1);
    free(rhoas2);
    free(f1);
    free(f2);
    free(lnf1);
    free(lnf2);

    return (0);
}


/* ---------------------------------------------------------------------------------------- */
/* aerosol() - compute aerosol reflectance using specified algorithm                        */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
int aerosol(l2str *l2rec, int32_t aer_opt_in, aestr *aerec, int32_t ip,
        float wave[], int32_t nwave, int32_t iwnir_s_in, int32_t iwnir_l_in,
        float F0_in[], float La1_in[], float La2_out[],
        float t_sol_out[], float t_sen_out[], float *eps, float taua_out[],
        int32_t *modmin, int32_t *modmax, float *modrat,
        int32_t *modmin2, int32_t *modmax2, float *modrat2) {
    static int firstCall = 1;
    static int32_t mindx[MAXAERMOD];

    int status = 1;
    float *rhoa;
    float *radref;
    float temp;
    float *aot;
    float angstrom;

    int iw, ipb;

    float *F0;
    float *taur;
    float *La1;
    float *La2;
    float *t_sol;
    float *t_sen;
    float *taua;
    float *taua_pred_min;
    float *taua_pred_max;

    l1str *l1rec = l2rec->l1rec;

    float *tg_sol_sm = l2rec->l1rec->tg_sol;
    float *tg_sen_sm = l2rec->l1rec->tg_sen;
    float *Lt_sm = l2rec->l1rec->Lt;

    int32_t sensorID = l1rec->l1file->sensorID;
    float solz = l1rec->solz [ip];
    float senz = l1rec->senz [ip];
    float wv = l1rec->wv [ip];
    float rh = l1rec->rh [ip];
    float pr = l1rec->pr [ip];

    if (firstCall == 1) {
        Nbands = nwave;
        Maxband = nwave + 1; /* Must be >= Nbands */
    }

    if ((rhoa = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa.\n");
        exit(1);
    }
    if ((radref = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for raderef.\n");
        exit(1);
    }
    if ((F0 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for F0.\n");
        exit(1);
    }
    if ((taur = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for taur.\n");
        exit(1);
    }
    if ((La1 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for La1.\n");
        exit(1);
    }
    if ((La2 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for La2.\n");
        exit(1);
    }
    if ((t_sol = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for t_sol.\n");
        exit(1);
    }
    if ((t_sen = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for t_sen.\n");
        exit(1);
    }
    if ((taua = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for taua.\n");
        exit(1);
    }
    if ((taua_pred_min = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for taua_pred_min.\n");
        exit(1);
    }
    if ((taua_pred_max = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for taua_pred_max.\n");
        exit(1);
    }

    /* set up the geometry structure to be passed  */
    if (l1rec->geom_per_band == NULL) {
        /* nominal geometry setup */
        geom.senz = &l1rec->senz[ip];
        geom.solz = &l1rec->solz[ip];
        geom.phi = &l1rec->delphi[ip];
        geom.csolz = &l1rec->csolz[ip];
        geom.csenz = &l1rec->csenz[ip];
        geom.gmult = 0;

        geom.airmass = (float *) malloc(sizeof (float));
        *geom.airmass = 1. / geom.csolz[0] + 1. / geom.csenz[0];
        if ((evalmask & TRANSSPHER) != 0) {
            geom.airmass_sph = (float *) malloc(sizeof (float));
            geom.airmass_plp = (float *) malloc(sizeof (float));
            *geom.airmass_sph = ky_airmass(geom.solz[0]) +
                    ky_airmass(geom.senz[0]);
            *geom.airmass_plp = pp_airmass(geom.solz[0]) +
                    pp_airmass(geom.senz[0]);
        }
    } else {
        ipb = ip * Nbands;
        geom.senz = &l1rec->geom_per_band->senz[ipb];
        geom.solz = &l1rec->geom_per_band->solz[ipb];
        geom.phi = &l1rec->geom_per_band->delphi[ipb];
        geom.csolz = &l1rec->geom_per_band->csolz[ipb];
        geom.csenz = &l1rec->geom_per_band->csenz[ipb];
        geom.gmult = 1;

        geom.airmass = (float *) malloc(Nbands * sizeof (float));
        for (iw = 0; iw < Nbands; iw++) {
            geom.airmass[iw] = 1. / geom.csolz[iw] + 1. / geom.csenz[iw];
        }
        if ((evalmask & TRANSSPHER) != 0) {
            geom.airmass_plp = (float *) malloc(Nbands * sizeof (float));
            geom.airmass_sph = (float *) malloc(Nbands * sizeof (float));
            for (iw = 0; iw < Nbands; iw++) {
                geom.airmass_plp[iw] = pp_airmass(geom.solz[iw]) +
                        pp_airmass(geom.senz[iw]);
                geom.airmass_sph[iw] = ky_airmass(geom.solz[iw]) +
                        ky_airmass(geom.senz[iw]);
            }
        }
    }

    /* set static global evaluation level */
    evalmask = input->evalmask;
    aer_opt = aer_opt_in;

    /* transfer inputs per band to inputs per wavelength */
    for (iw = 0; iw < nwave; iw++) {
        F0 [iw] = F0_in [iw];
        taur[iw] = l1rec->l1file->Tau_r[iw];
        La1 [iw] = La1_in[iw];
        if (iw == iwnir_s_in) iwnir_s = iw;
        if (iw == iwnir_l_in) iwnir_l = iw;
    }

    /* compute total airmass (static global) */
    mu0 = cos(solz / radeg);
    mu = cos(senz / radeg);
    airmass = 1.0 / mu0 + 1.0 / mu;
    if ((evalmask & TRANSSPHER) != 0) {
        airmass_plp = pp_airmass(solz) + pp_airmass(senz);
        airmass_sph = ky_airmass(solz) + ky_airmass(senz);
    }

    /* initialize epsilon and diffuse transmittance to defaults */
    *eps = 1.0;
    *modmin = BAD_INT;
    *modmax = BAD_INT;
    *modrat = BAD_FLT;
    *modmin2 = BAD_INT;
    *modmax2 = BAD_INT;
    *modrat2 = BAD_FLT;
    for (iw = 0; iw < nwave; iw++) {
        taua [iw] = 0.0;
        t_sol [iw] = 1.0;
        t_sen [iw] = 1.0;
        La2 [iw] = 0.0;
        radref[iw] = pi / F0[iw] / geom.csolz[iw * geom.gmult];
    }

    /* load aerosol model tables */
    if (!loaded) {
        int32_t im;
        load_aermod(sensorID, wave, nwave, input->aermodfile, input->aermodels, input->naermodels);
        for (im = 0; im < aertab->nmodel; im++) mindx[im] = im;
        if (have_rh && aertab->nmodel >= 30) {
            printf("\nLimiting aerosol models based on RH.\n");
            use_rh = 1;
        }
    }
    /* do not permit use of geom_per_band if interpolation is done */
    /* WDR note that it CAN work, but currently, we don't want users to
       use this combination */
    if ((interpol == 1) && (geom.gmult == 1)) {
        fprintf(stderr, "-E- %s line %d: Interpolated aerosol tables are\n",
                __FILE__, __LINE__);
        fprintf(stderr, "           not permitted for use with band-dependent geometry, set geom_per_band=0\n");
        exit(1);
    }

    /* Change the aerosol option if need be */
    if (use_rh)
        switch (aer_opt) {
        case AERWANG: aer_opt = AERRH;
            break;
        case AERWANGNIR: aer_opt = AERRHNIR;
            break;
        case AERWANGSWIR: aer_opt = AERRHSWIR;
            break;
        case AERMUMM: aer_opt = AERRHMUMM;
            break;
        }


    if (firstCall) {
        if (aer_opt > 0 && aer_opt <= MAXAERMOD) {
            printf("\nUsing fixed aerosol model #%d (%s)\n", aer_opt, input->aermodels[aer_opt - 1]);
            printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
        } else {
            switch (aer_opt) {
            case AERRHSM: // Spectral Matching --> Amir
                printf("\nUsing Spectral Matching of aerosols reflectance for\n");
                printf("wavelength from %4.1f nm to %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                break;
            case AERWHITE:
                printf("\nUsing white-aerosol approximation\n");
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case AERWANG:
            case AERRH:
                printf("\nUsing Gordon & Wang aerosol model selection\n");
                printf("Using bands at %4.1f and %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case AERWANGNIR:
            case AERRHNIR:
            case AERSS14:
                printf("\nUsing Gordon & Wang aerosol model selection\n");
                printf("  and NIR correction with up to %d iterations\n", input->aer_iter_max);
                printf("Using bands at %4.1f and %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case AERWANGSWIR:
            case AERRHSWIR:
                printf("\nUsing Gordon & Wang aerosol model selection with NIR/SWIR switching.\n");
                printf("NIR correction with up to %d iterations\n", input->aer_iter_max);
                printf("NIR bands at %d and %d nm\n", input->aer_wave_short, input->aer_wave_long);
                printf("SWIR bands at %d and %d nm\n\n", input->aer_swir_short, input->aer_swir_long);
                break;
            case AERMUMM:
            case AERRHMUMM:
                printf("\nUsing Gordon & Wang aerosol model selection\n");
                printf("  and MUMM correction\n");
                printf("Using bands at %4.1f and %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case AERRHFRNIR:
                printf("\nUsing Ahmad & Franz aerosol model selection.\n");
                printf("  and NIR correction with up to %d iterations\n", input->aer_iter_max);
                printf("Using bands at %4.1f and %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case AERRHMSEPS:
                printf("\nUsing multi-scattering aerosol model selection.\n");
                printf("  and NIR correction with up to %d iterations\n", input->aer_iter_max);
                printf("Using bands at %4.1f and %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case AERRHMSEPS_lin:
                printf("\nUsing multi-scattering aerosol model selection in linear space.\n");
                printf("  and NIR correction with up to %d iterations\n", input->aer_iter_max);
                printf("Using bands at %4.1f and %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case AERSSP:
                printf("\nUsing Gordon & Wang aerosol model selection\n");
                printf("  and Spectral shape parameter calculated using rhos\n");
                printf("Using bands at %4.1f and %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case FIXAOT:
                printf("\nUsing fixed, input aerosol optical thicknesses for aerosol selection.\n");
                break;
            case FIXMODPAIR:
                printf("\nUsing multiple scattering aerosols with fixed model pair\n");
                break;
            case FIXMODPAIRNIR:
                printf("\nUsing multiple scattering aerosols with fixed model pair\n");
                printf("  and NIR iteration with up to %d iterations\n", input->aer_iter_max);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case FIXANGSTROM:
                printf("\nUsing fixed aerosol model based on predefined Angstrom exponent\n");
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case FIXANGSTROMNIR:
                printf("\nUsing fixed aerosol model based on predefined Angstrom exponent\n");
                printf("  and NIR iteration with up to %d iterations\n", input->aer_iter_max);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            default:
                printf("\nErroneous aerosol option, %d\n", aer_opt);
                exit(FATAL_ERROR);
                break;
            }
        }
        firstCall = 0;
    }


    /* convert input aerosol radiances to relectance */
    for (iw = iwnir_s; iw <= iwnir_l; iw += MAX(iwnir_l - iwnir_s, 1))
        rhoa[iw] = La1[iw] * radref[iw];


    /* compute aerosol using requested method */
    /* -------------------------------------- */

    switch (aer_opt) {

    case AERWANG: case AERWANGNIR: case AERWANGSWIR: case AERMUMM:

        if (iwnir_l <= iwnir_s || wave[iwnir_s] < 600) {
            printf("Aerosol selection bands must be greater than 600nm with short wave less than long wave (%d,%d)\n", iwnir_l, iwnir_s);
            exit(1);
        }

        /* Multi-Scattering with Gordon & Wang Aerosol Selection */
        /* ----------------------------------------------------- */

        /* convert input NIR aerosol radiances to reflectance */
        for (iw = iwnir_s; iw <= iwnir_l; iw++) {
            rhoa[iw] = La1[iw] * radref[iw];
        }

        /* require sufficient signal in two NIR bands */
        if (rhoa[iwnir_s] > input->rhoamin && rhoa[iwnir_l] > input->rhoamin) {

            /* require MS epsilon to be reasonable */
            if (La1[iwnir_s] / La1[iwnir_l] > 0.1) {

                status = wangaer(sensorID, wave, nwave, iwnir_s, iwnir_l,
                        aertab->nmodel, mindx,
                        &geom, wv, rhoa, modmin, modmax, modrat, eps, taua, taua); // this taua is not used //

                if (status == 0)
                    for (iw = 0; iw < nwave; iw++)
                        La2[iw] = rhoa[iw] / radref[iw];
            }

        } else if (rhoa[iwnir_s] > -(input->rhoamin) && rhoa[iwnir_l] > -(input->rhoamin)) {

            /* if input NIR is near zero, assume a small white aerosol signal */
            *eps = 1.0;
            *modmin = aertab->nmodel;
            *modmax = aertab->nmodel;
            *modrat = 0.0;
            temp = MAX(rhoa[iwnir_l], 1e-6);
            for (iw = 0; iw < nwave; iw++) {
                rhoa[iw] = temp;
                La2 [iw] = rhoa[iw] / radref[iw];
            }

            status = 0;

        } else {

            /* if input NIR is very negative, fail the case */
            status = 1;
        }

        break;

    case AERRH:
    case AERRHNIR:
    case AERRHSWIR:
    case AERRHMUMM:
    case AERRHFRNIR:
    case AERRHMSEPS:
    case AERRHMSEPS_lin:
    case AERRHSM: // spectral matching --> Amir
    case AERSSP:
    case AERSS14:

        if (iwnir_l <= iwnir_s || wave[iwnir_s] < 600) {
            printf("Aerosol selection bands must be greater than 600nm with short wave less than long wave");
            exit(1);
        }

        /* Multi-Scattering with RH-based selection              */
        /* ----------------------------------------------------- */

        /* convert input NIR aerosol radiances to relectance */
        for (iw = iwnir_s; iw <= iwnir_l; iw++) {
            rhoa[iw] = La1[iw] * radref[iw];
        }

        /* require sufficient signal in two NIR bands */
        if (rhoa[iwnir_s] > input->rhoamin && rhoa[iwnir_l] > input->rhoamin) {

            /* require MS epsilon to be reasonable */
            if (rhoa[iwnir_s] / rhoa[iwnir_l] > 0.1 && rhoa[iwnir_s] / rhoa[iwnir_l] < 10.0) {

                status = rhaer(sensorID, wave, nwave, iwnir_s, iwnir_l,
                        &geom, wv, rh, pr, taur, rhoa,
                        modmin, modmax, modrat, modmin2, modmax2, modrat2, eps, taua, t_sol, t_sen, tg_sol_sm, tg_sen_sm, Lt_sm, ip);
                status = 0;
                if (status == 0)
                    for (iw = 0; iw < nwave; iw++)
                        La2[iw] = rhoa[iw] / radref[iw];
            }

        } else if (rhoa[iwnir_s] > -(input->rhoamin) && rhoa[iwnir_l] > -(input->rhoamin)) {

            /* if input NIR is near zero, assume a small white aerosol signal */
            *eps = 1.0;
            *modmin = aertab->nmodel;
            *modmax = aertab->nmodel;
            *modrat = 0.0;
            *modmin2 = aertab->nmodel;
            *modmax2 = aertab->nmodel;
            *modrat2 = 0.0;
            temp = MAX(rhoa[iwnir_l], 1e-6);
            for (iw = 0; iw < nwave; iw++) {
                rhoa[iw] = temp;
                La2 [iw] = rhoa[iw] / radref[iw];
            }

            diff_tran(sensorID, wave, nwave, iwnir_l, &geom, wv, pr, taur,
                    *modmin, *modmax, *modrat, rhoa, taua, t_sol, t_sen, taua_pred_min, taua_pred_max, 0);

            status = 0;

        } else {

            /* if input NIR is very negative, fail the case */
            status = 1;
        }

        break;

    case AERWHITE:

        /* White Aerosol */
        /* ------------- */

        if (La1[iwnir_l] > 0.0) {

            *eps = 1.0;
            *modmin = 0;
            *modmax = 0;
            *modrat = 0.0;

            for (iw = 0; iw < nwave; iw++) {
                La2 [iw] = *eps * F0[iw] / F0[iwnir_l] * La1[iwnir_l];
                rhoa[iw] = La2[iw] * radref[iw];
            }


            status = 0;
        }
        break;

    case FIXMODPAIR: case FIXMODPAIRNIR:

        /* Multi-Scattering with Fixed Model Pair */
        /* -------------------------------------- */

        if (aerec != NULL && aerec->mode == ON) {
            *modmin = aerec->mod_min[ip] - 1;
            *modmax = aerec->mod_max[ip] - 1;
            *modrat = aerec->mod_rat[ip];
        } else {
            *modmin = input->aermodmin - 1;
            *modmax = input->aermodmax - 1;
            *modrat = input->aermodrat;
        }

        status = fixedmodpair(sensorID, wave, nwave, iwnir_s, iwnir_l, &geom, wv,
                *modmin, *modmax, *modrat, rhoa, eps);
        if (status == 0) {
            for (iw = 0; iw < nwave; iw++)
                La2[iw] = rhoa[iw] / radref[iw];
        }

        break;

    case FIXANGSTROM: case FIXANGSTROMNIR:

        if (input->aer_angstrom > -2) {
            angstrom = input->aer_angstrom;
        } else {
            int16_t year, day;
            double sec;
            unix2yds(l2rec->l1rec->scantime, &year, &day, &sec);
            angstrom = bin_climatology(input->aerbinfile, day, l1rec->lon[ip], l1rec->lat[ip], "angstrom");
        }

        if (angstrom > -2) {

            model_select_angstrom(angstrom, modmin, modmax, modrat);

            status = fixedmodpair(sensorID, wave, nwave, iwnir_s, iwnir_l, &geom, wv,
                    *modmin, *modmax, *modrat, rhoa, eps);
        } else
            status = 1;

        if (status == 0) {
            for (iw = 0; iw < nwave; iw++)
                La2[iw] = rhoa[iw] / radref[iw];
        }

        break;

    case FIXAOT:

        /* Multi-Scattering with fixed AOT */
        /* ------------------------------- */
        if (aerec != NULL && aerec->mode == ON)
            aot = &aerec->taua[ip * Nbands];
        else
            aot = input->taua;

        status = fixedaot(sensorID, aot, wave, nwave, iwnir_s, iwnir_l, &geom, wv,
                modmin, modmax, modrat, rhoa, eps);
        if (status == 0) {
            for (iw = 0; iw < nwave; iw++)
                La2[iw] = rhoa[iw] / radref[iw];
        }

        break;

    default:

        /* Multi-Scattering with Fixed Model */
        /* --------------------------------- */

        *modmin = aer_opt - 1;
        *modmax = aer_opt - 1;
        *modrat = 0.0;

        if (*modmin < 0 || *modmin > input->naermodels - 1) {
            printf("Invalid aerosol option: %d\n", *modmin);
            exit(1);
        }

        /* convert input NIR aerosol radiance to relectance */
        rhoa[iwnir_l] = La1[iwnir_l] * radref[iwnir_l];

        /* get aerosol reflectance in visible for fixed model, and epsilon */
        status = fixedaer(sensorID, *modmin, wave, nwave, iwnir_s, iwnir_l,
                input->aermodels, input->naermodels,
                &geom, wv, rhoa, eps);

        /* convert aerosol relectance to radiance */
        if (status == 0)
            for (iw = 0; iw < nwave; iw++)
                La2[iw] = rhoa[iw] / radref[iw];
        break;
    }

    /* compute diffuse transmittance through aerosol and Rayleigh, if not yet computed */
    if (status == 0 && aer_opt != AERRHNIR && aer_opt != AERRHFRNIR && aer_opt != AERRHSWIR && aer_opt != AERRH && aer_opt != AERRHMSEPS && aer_opt != AERRHSM) {

        diff_tran(sensorID, wave, nwave, iwnir_l, &geom, wv, pr, taur,
                *modmin, *modmax, *modrat, rhoa, taua, t_sol, t_sen, taua_pred_min, taua_pred_max, 0);
    }

    /* transfer outputs per wavelength to outputs per band */
    for (iw = 0; iw < nwave; iw++) {
        La2_out [iw] = La2 [iw];
        taua_out [iw] = taua [iw];
        t_sol_out[iw] = t_sol[iw];
        t_sen_out[iw] = t_sen[iw];
    }

    /* model numbers are reported as 1-based */
    *modmin = *modmin + 1;
    *modmax = *modmax + 1;
    *modmin2 = *modmin2 + 1;
    *modmax2 = *modmax2 + 1;

    free(rhoa);
    free(radref);
    free(F0);
    free(taur);
    free(La1);
    free(La2);
    free(t_sol);
    free(t_sen);
    free(taua);
    free(taua_pred_min);
    free(taua_pred_max);
    free(geom.airmass);
    if ((evalmask & TRANSSPHER) != 0) {
        free(geom.airmass_plp);
        free(geom.airmass_sph);
    }

    return (status);
}



/* --------------------------------------------------------------- */
/* get_angstrom.c - compute angstrom coefficient                   */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/*     band  = band number 0-7                                     */
/* Outputs:                                                        */
/*     angst - angstrom coefficient                                */
/*                                                                 */
/* Written By: B. A. Franz, SAIC, August 2004                      */
/*                                                                 */

/* --------------------------------------------------------------- */
void get_angstrom(l2str *l2rec, int band, float angst[]) {
    static int firstCall = 1;
    static int32_t ib2;
    static float wave2;

    int32_t ip;
    int32_t ib1;
    float wave1;
    float aot1;
    float aot2;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;

    if (firstCall) {
        ib2 = windex(865.0, l1file->fwave, l1file->nbands);
        wave2 = l1file->fwave[ib2];
        firstCall = 0;
    }

    if (band < 0)
        ib1 = windex(443.0, l1file->fwave, l1file->nbands);
    else
        ib1 = band;
    wave1 = l1file->fwave[ib1];

    for (ip = 0; ip < l1rec->npix; ip++) {
        aot1 = l2rec->taua[ip * l1file->nbands + ib1];
        aot2 = l2rec->taua[ip * l1file->nbands + ib2];
        if (aot1 > 0.0 && aot2 > 0.0)
            angst[ip] = -log(aot1 / aot2) / log(wave1 / wave2);
        else if (aot1 == 0.0 && aot2 == 0.0)
            angst[ip] = 0.0;
        else {
            angst[ip] = BAD_FLT;
            l1rec->flags[ip] |= PRODFAIL;
        }
        //     printf("angst = %f aot1= %f aot2= %f ib1=%d ib2=%d\n",angst[ip],aot1,aot2,ib1,ib2);
    }

    return;
}


/* --------------------------------------------------------------- */
/* get_ms_epsilon.c -                                              */

/* --------------------------------------------------------------- */
void get_ms_epsilon(l2str *l2rec, float eps[]) {
    int32_t ip, ipb1, ipb2;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;

    for (ip = 0; ip < l1rec->npix; ip++) {
        ipb1 = ip * l1file->nbands + iwnir_s;
        ipb2 = ip * l1file->nbands + iwnir_l;
        if (l2rec->La[ipb2] > 0.0) {
            eps[ip] = l2rec->La[ipb1] / l2rec->La[ipb2]
                    * l1rec->Fo[iwnir_l] / l1rec->Fo[iwnir_s];
        } else {
            /* "should" indicate ATMFAIL, but just to be safe */
            eps[ip] = BAD_FLT;
            l1rec->flags[ip] |= PRODFAIL;
        }
    }

    return;
}
