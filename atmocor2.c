#include "l12_proto.h"
#include "atrem_corl1.h"

/* --------------------------------------------------------------------------------------- */
/* atmocor2() - atmospheric correction, converts Lt -> nLw                                 */
/*                                                                                         */
/* C-version, September 2004, B. Franz                                                     */
/* Added SS14 and SSP aerosol corredtion, October 2021, R. K. Singh                        */
/* --------------------------------------------------------------------------------------- */

int atmocor2(l2str *l2rec, aestr *aerec, int32_t ip) {
    static int firstCall = 1;
    static int want_nirRrs = 0;
    static int want_nirLw = 0;
    static int want_mumm = 0;
    static int want_ssp = 0;
    static int want_ramp = 1;
    static int32_t aer_iter_max = 1;
    static int32_t aer_iter_min = 1;

    static float pi = PI;
    static float radeg = RADEG;
    static float badval = BAD_FLT;
    static float badchl = BAD_FLT;
    static float p0 = STDPR;
    static float df = 0.33;
    static float cbot = 0.7;
    static float ctop = 1.3;
    static float seed_chl = 0.0;
    static float seed_green = 0.0;
    static float seed_red = 0.0;
    static float nir_chg = 0.02;
    static float glint_min = GLINT_MIN;
    static float cslp;
    static float cint;

    static int32_t green;
    static int32_t red;
    static int32_t nir_s;
    static int32_t nir_l;
    static int32_t swir_s;
    static int32_t swir_l;
    static int32_t aer_s;
    static int32_t aer_l;
    static int32_t aer_ref;   //Declaration for SSP
    static int32_t daer;
    static int32_t nwvis;
    static int32_t nwave;
    static float *wave;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;

    int32_t sensorID = l1file->sensorID;
    int32_t brdf_opt = input->brdf_opt;
    int32_t aer_opt = input->aer_opt;
    int32_t glint_opt = input->glint_opt;
    int32_t cirrus_opt = input->cirrus_opt;

    float *Fo = l1rec->Fo;
    float *Fobar = l1file->Fobar;
    float fsol = l1rec->fsol;
    float solz = l1rec->solz [ip];
    float senz = l1rec->senz [ip];
    float delphi = l1rec->delphi[ip];
    float ws = l1rec->ws [ip];
    float gc = l1rec->glint_coef[ip];

    int32_t *aermodmin = &l2rec->aermodmin[ip];
    int32_t *aermodmax = &l2rec->aermodmax[ip];
    float *aermodrat = &l2rec->aerratio [ip];
    int32_t *aermodmin2 = &l2rec->aermodmin2[ip];
    int32_t *aermodmax2 = &l2rec->aermodmax2[ip];
    float *aermodrat2 = &l2rec->aerratio2 [ip];
    float *eps = &l2rec->eps [ip];

    float *TLg = &l1rec->TLg [ip * l1file->nbands];
    float *Lw = &l2rec->Lw [ip * l1file->nbands];
    float *nLw = &l2rec->nLw [ip * l1file->nbands];
    float *La = &l2rec->La [ip * l1file->nbands];
    float *taua = &l2rec->taua [ip * l1file->nbands];
    float *t_sol = &l1rec->t_sol[ip * l1file->nbands];
    float *t_sen = &l1rec->t_sen[ip * l1file->nbands];
    float *brdf = &l2rec->brdf [ip * l1file->nbands];
    float *Rrs = &l2rec->Rrs [ip * l1file->nbands];

    float *nLw_unc = &l2rec->nLw_unc[ip * l1file->nbands];
    float *Rrs_unc = &l2rec->Rrs_unc[ip * l1file->nbands];

    float *tg_sol = &l1rec->tg_sol[ip * l1file->nbands];

    float *taur;
    float *tLw;
    float *rhown_nir;
    float *tLw_nir, *last_tLw_nir;
    int32_t ib, ipb;
    int32_t status;
    int32_t iter, iter_max, iter_min, last_iter, iter_reset;
    float mu, mu0;
    int32_t nneg;
    float airmass;
    float *Ltemp;
    float *Lunc;
    float chl;
    int want_glintcorr;
    float refl_nir = 0, last_refl_nir = 0;
    float tindx;
    float Ka = 0.8; /* 1.375 um channel transmittance for water vapor (<=1)  Gao et al. 1998 JGR */
    float kappa=0, Lref=0;	  /*Declaration for SSP & SS14*/
    float laer1, laer2, efact;	    /*Declaration for SSP & SS14*/

    if (firstCall == 1) {
        firstCall = 0;
        nwave = l1file->nbands;
        if ((wave = (float *) calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to wave\n");
            exit(FATAL_ERROR);
        }
        for (ib = 0; ib < nwave; ib++) {
            wave[ib] = l1file->fwave[ib];
        }
        if (sensorID == CZCS) {
            want_ramp = 0;
            seed_chl = 0.01;
            seed_green = 0.003;
            seed_red = 0.00125;
        }
        nwvis = rdsensorinfo(l1file->sensorID, input->evalmask, "NbandsVIS", NULL);
        nir_s = bindex_get(input->aer_wave_short);
        nir_l = bindex_get(input->aer_wave_long);
        if (nir_s < 0 || nir_l < 0) {
            printf("Aerosol selection bands %d and %d not available for this sensor\n",
                    input->aer_wave_short, input->aer_wave_long);
            exit(1);
        }
        if (nir_l < nir_s) {
            printf("Invalid aerosol selection bands: long (%d) must be greater than short (%d).\n",
                    input->aer_wave_long, input->aer_wave_short);
            exit(1);
        }
        if (wave[nir_s] < 600) {
            printf("Aerosol selection band(s) must be greater than 600nm");
            exit(1);
        }

        aer_s = nir_s;
        aer_l = nir_l;
        daer = MAX(nir_l - nir_s, 1);
        cslp = 1. / (ctop - cbot);
        cint = -cslp * cbot;

        if (aer_opt != AERSSP)
          printf("Aerosol selection bands %d and %d\n", l1file->iwave[aer_s], l1file->iwave[aer_l]);

        switch (aer_opt) {
        case AERWANGNIR:
        case AERRHNIR:
        case AERRHFRNIR:
        case FIXANGSTROMNIR:
        case FIXMODPAIRNIR:
            want_nirLw = 1;
            aer_iter_min = 1;
            aer_iter_max = input->aer_iter_max;
            printf("NIR correction enabled.\n");
            break;
        case AERMUMM:
        case AERRHMUMM:
            want_mumm = 1;
            aer_iter_min = 3;
            aer_iter_max = input->aer_iter_max;
            printf("MUMM correction enabled.\n");
            break;
        case AERWANGSWIR:
        case AERRHSWIR:
            want_nirLw = 1;
            aer_iter_min = 1;
            aer_iter_max = input->aer_iter_max;
            swir_s = bindex_get(input->aer_swir_short);
            swir_l = bindex_get(input->aer_swir_long);
            if (swir_s < 0 || swir_l < 0) {
                printf("Aerosol selection bands %d and %d not available for this sensor\n",
                        input->aer_swir_short, input->aer_swir_long);
                exit(1);
            }
            printf("NIR/SWIR switching correction enabled.\n");
            break;
        case AERRHMSEPS:
            want_nirLw = 1;
            aer_iter_min = 1;
            aer_iter_max = input->aer_iter_max;
            printf("NIR correction enabled --> for multi-scattering epsilon.\n");
            break;
        case AERRHSM:
            want_nirLw = 1; //This needs to be turned on, but a new SWIR water model is needed for it to work
            aer_iter_min = 1;
            aer_iter_max = input->aer_iter_max;
            printf("NIR correction enabled --> for spectral matching.\n");
            break;
        case AERSS14:
        case AERSSP:
            want_ssp = 1;
            //want_nirLw = 1; //This option is on but will not be used to estimate NIR Lw
            aer_iter_min = 1;
            aer_iter_max = 1;
            printf("SSP correction enabled --> Spectral shape parameter will be calculated using rhos.\n");
            break;
        default:
            if (input->aer_rrs_short >= 0.0 && input->aer_rrs_long >= 0.0) {
                want_nirRrs = 1;
                aer_iter_min = 3;
                aer_iter_max = input->aer_iter_max;
                printf("NIR correction via input Rrs enabled.\n");
            }
            break;
        }
        if (input->aer_iter_max < 1)
            want_nirLw = 0;
        if ( want_nirLw || want_nirRrs || want_ssp) {
            if ((red = bindex_get(670)) < 0) {
                if ((red = bindex_get(680)) < 0)
                    if ((red = bindex_get(620)) < 0)
                        if ((red = bindex_get(765)) < 0)
                            if ((red = bindex_get(655)) < 0)
                                if ((red = bindex_get(664)) < 0) /* added for MSI */ {
                                    printf("%s line %d: can't find red band\n", __FILE__, __LINE__);
                                    exit(1);
                                }
            }
            if ((green = bindex_get(550)) < 0) {
                if ((green = bindex_get(555)) < 0)
                    if ((green = bindex_get(560)) < 0)
                        if ((green = bindex_get(565)) < 0) {
                            printf("%s line %d: can't find green band\n", __FILE__, __LINE__);
                            exit(1);
                        }
            }
        }

        if (aer_opt == AERSSP)
        { aer_l = aer_s;
          aer_s = red;
          printf("Aerosol selection bands %d and %d\n", l1file->iwave[aer_s], l1file->iwave[aer_l]);
        }

        if (want_ssp)
        { aer_ref = red;
          printf("Aerosol reference band is %d\n\n",l1file->iwave[aer_ref]);
        }
    }

    if ((taur = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to taur\n");
        exit(FATAL_ERROR);
    }
    if ((tLw = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to tLw\n");
        exit(FATAL_ERROR);
    }
    if ((rhown_nir = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to rhown_nir\n");
        exit(FATAL_ERROR);
    }
    if ((tLw_nir = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to tLw_nir\n");
        exit(FATAL_ERROR);
    }
    if ((last_tLw_nir = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to last_tLw_nir\n");
        exit(FATAL_ERROR);
    }
    if ((Ltemp = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to Ltemp\n");
        exit(FATAL_ERROR);
    }
    if ((Lunc = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to Lunc\n");
        exit(FATAL_ERROR);
    }

    /* Initialize output values. If any input radiances are negative, just return. */

    status = 0;
    iter = 0;
    nneg = 0;

    for (ib = 0; ib < nwave; ib++) {

        ipb = ip * nwave + ib;

        //t_sol [ib]  = 1.0;    leave them as rayleigh only
        //t_sen [ib]  = 1.0;
        La [ib] = badval;
        tLw [ib] = badval;
        Lw [ib] = badval;
        nLw [ib] = badval;
        taua [ib] = badval;
        if (glint_opt != 2) TLg [ib] = 0.0;
        brdf [ib] = 1.0;
        Rrs [ib] = badval;

        l2rec->eps[ip] = badval;

        if (l2rec->l1rec->Lt[ipb] <= 0.0)
            if (wave[ib] < 1000) nneg++; /* don't test SWIR */
    }

    /* If any expected channels are negative */
    if (nneg > 0) {
        free(taur);
        free(tLw);
        free(rhown_nir);
        free(tLw_nir);
        free(last_tLw_nir);
        free(Ltemp);
        free(Lunc);

        status = 1;
        return (status);
    }

    mu0 = cos(solz / radeg);
    mu = cos(senz / radeg);
    airmass = 1.0 / mu0 + 1.0 / mu;

    /* Remove pre-computed atmospheric effects */
    /* --------------------------------------- */
    for (ib = 0; ib < nwave; ib++) {

        ipb = ip * nwave + ib;

        /* Pressure correct the Rayleigh optical thickness */
        taur[ib] = l1rec->pr[ip] / p0 * l1file->Tau_r[ib];

        /* Copy TOA radiance to temp var, eliminate missing bands */
        Ltemp[ib] = l2rec->l1rec->Lt[ipb];
        Lunc [ib] = l2rec->l1rec->Lt_unc[ipb];

        /* Correct for ozone absorption.  We correct for inbound and outbound here, */
        /* then we put the inbound back when computing Lw.                          */
        Ltemp[ib] = Ltemp[ib] / l1rec->tg_sol[ipb] / l1rec->tg_sen[ipb];
        Lunc [ib] = Lunc [ib] / l1rec->tg_sol[ipb] / l1rec->tg_sen[ipb];

        /* Do Cirrus correction - subtract off cirrus reflectance from Ltemp */
        /* Ka is the 1.375 um transmittance of water vapor above cirrus clouds */
        /* Add cirrus_opt to input */
        if (cirrus_opt) Ltemp[ib] -= l1rec->rho_cirrus[ip] / Ka * Fo[ib] * mu0 / pi;

        /*  Apply polarization correction */
        Ltemp[ib] /= l1rec->polcor[ipb];
        Lunc [ib] /= l1rec->polcor[ipb];

        /* Remove whitecap radiance */
        Ltemp[ib] -= l1rec->tLf[ipb];

        /* Subtract the Rayleigh contribution for this geometry. */
        Ltemp[ib] -= l1rec->Lr[ipb];

        /* Correct for O2 absorption (replace O2 losses) */
        if (input->oxaband_opt != 0 && !((input->gas_opt & ATREM_BIT) && (input->atrem_opt & ATREM_O2))) {
            Ltemp[ib] /= l1rec->t_o2[ipb];
            Lunc [ib] /= l1rec->t_o2[ipb];
        }
    }


    /* Compute BRDF correction, if not dependent on radiance */
    if (brdf_opt < FOQMOREL && brdf_opt > NOBRDF)
        ocbrdf(l2rec, ip, brdf_opt, wave, nwvis, solz, senz, delphi, ws, -1.0, NULL, NULL, brdf);

    /* Initialize iteration loop */
    chl = seed_chl;
    iter = 0;
    last_iter = 0;
    iter_max = aer_iter_max;
    iter_min = aer_iter_min;
    iter_reset = 0;
    last_refl_nir = 100.;
    want_glintcorr = 0;

    /*  new glint_opt usage - a 2 will use the simple TLg from atmocor1 */
    if (glint_opt == 1 && l1rec->glint_coef[ip] > glint_min) {
        iter_max = MAX(2, iter_max);
        iter_min = MAX(2, iter_min);
        want_glintcorr = 1;
    }
    if (glint_opt == 2) {
        want_glintcorr = 1;
    }

    if (aer_opt == AERWANGSWIR || aer_opt == AERRHSWIR) {
        tindx_shi(l2rec, ip, &tindx);
        if (tindx >= 1.05) {
            iter_max = 1;
            aer_s = swir_s;
            aer_l = swir_l;
            want_nirLw = 0;
        } else {
            aer_s = nir_s;
            aer_l = nir_l;
            want_nirLw = 1;
        }
        daer = MAX(aer_l - aer_s, 1);
    }
NIRSWIR:

    if (want_nirLw || want_nirRrs) {
        for (ib = 0; ib < nwave; ib++) {
            last_tLw_nir[ib] = 0.0;
            tLw_nir[ib] = 0.0;
            Rrs[ib] = 0.0;
        }
        Rrs[green] = seed_green;
        Rrs[red ] = seed_red;
    }

    if (want_ssp) {
        for (ib = 0; ib < nwave; ib++) {
            tLw[ib] = (Ltemp[ib]*pi)/Fo[ib];
        }
        kappa=calc_kappa(sensorID,tLw);
    }
    /* -------------------------------------------------------------------- */
    /* Begin interations for aerosol with corrections for non-zero nLw(NIR) */
    /* -------------------------------------------------------------------- */

    while (!last_iter) {

        iter++;
        status = 0;

        /* Initialize tLw as surface + aerosol radiance */
        for (ib = 0; ib < nwave; ib++) {
            tLw[ib] = Ltemp[ib];
        }

        /*  Compute and subtract glint radiance */
        if (want_glintcorr) {

            if (glint_opt == 1)
                glint_rad(iter, nwave, aer_s, aer_l, gc, airmass, mu0, Fo, taur, taua, tLw, TLg);


            for (ib = 0; ib < nwave; ib++) {
                tLw[ib] -= TLg[ib];
            }
        }

        /* Adjust for non-zero NIR water-leaving radiances using input Rrs */
        if (want_nirRrs) {

            rhown_nir[aer_s] = pi * input->aer_rrs_short;
            rhown_nir[aer_l] = pi * input->aer_rrs_long;

            for (ib = aer_s; ib <= aer_l; ib += daer) {

                /* Convert NIR reflectance to TOA W-L radiance */
                tLw_nir[ib] = rhown_nir[ib] / pi * Fo[ib] * mu0 * t_sol[ib] * t_sen[ib] / brdf[ib];

                /* Avoid over-subtraction */
                if (tLw_nir[ib] > tLw[ib] && tLw[ib] > 0.0)
                    tLw_nir[ib] = tLw[ib];

                /* Remove estimated NIR water-leaving radiance */
                tLw[ib] -= tLw_nir[ib];
            }
        }


        /* Adjust for non-zero NIR water-leaving radiances using MUMM */
        if (want_mumm) {

            get_rhown_mumm(l2rec, ip, aer_s, aer_l, rhown_nir);

            for (ib = aer_s; ib <= aer_l; ib += daer) {

                /* Convert NIR reflectance to TOA W-L radiance */
                tLw_nir[ib] = rhown_nir[ib] / pi * Fo[ib] * mu0 * t_sol[ib] * t_sen[ib] / brdf[ib];

                /* Remove estimated NIR water-leaving radiance */
                tLw[ib] -= tLw_nir[ib];
            }
        }

        /* Aerosol correction using SSP algorithm (Singh et al. 2019) */
        if (aer_opt == AERSSP) {

          last_iter = 1;

          laer1 = tLw[aer_s];
          laer2 = tLw[aer_l];

          efact=1.0;
          if (sensorID == L7ETMP || sensorID == L5TM || sensorID == OLI)
            efact=0.8; // was 0.7
          tLw[aer_l]/=efact;      // relating 86x nm band to 75x nm band

     			Lref = 3.585*tLw[aer_l]*exp(-531.0/547.0)/(kappa*kappa);//reference aerosol

          efact=Lref*exp(-pow((wave[aer_s]-531.0)/380.0,2.0));

            if (efact>tLw[aer_s])
     					Lref = Lref - 1.05*(efact-tLw[aer_s]);    // correcting for over-correction

          tLw[aer_s] = Lref*exp(-pow((wave[aer_s]-531.0)/380.0,2.0));
     			tLw[aer_l] = Lref*exp(-pow((wave[aer_l]-531.0)/380.0,2.0));

        }

        /* Aerosol correction using Gaussian extrapolation with spectral shape parameter (Singh and Shanmugam, 2014) */
        if (aer_opt == AERSS14) {

          last_iter = 1;

          efact=1.0;
          if (sensorID == L7ETMP || sensorID == L5TM || sensorID == OLI)
          { aer_l = aer_s;
            aer_s = red;
            efact=0.70;
          }

          tLw[aer_l]/=efact;      // relating 86x nm band to 75x nm band

			    Lref = 3.585*tLw[aer_s]*exp(-531.0/547.0)/(kappa*kappa);			//reference aerosol

			    if (((1.2/(kappa*kappa))*tLw[aer_s])>tLw[red])
				      Lref = Lref - (((1.3/(kappa*kappa))*tLw[aer_s])-tLw[red]); 			// correction at 667nm

          for(ib=0;ib<=nwave;ib++)
					     La[ib] = Lref*exp(-pow((wave[ib]-531.0)/380.0,2.0))+tLw[aer_l]*exp(-pow((wave[ib]-465.0)/(469.0-443.0),2.0))/kappa;		//final expression for La

        	La[aer_l] = tLw[aer_l]-0.77*tLw[aer_s]+0.77*La[aer_s];

          for(ib=0;ib<=nwave;ib++)
          {  Ltemp[ib] = tLw[ib];
             tLw[ib] = La[ib];
          }

        }

        /* Adjust for non-zero NIR water-leaving radiances using IOP model */
        if (want_nirLw) {

            ipb = ip*nwave;
            get_rhown_eval(input->fqfile, Rrs, wave, aer_s, aer_l, nwave, &l1rec->sw_a_avg[ipb], &l1rec->sw_bb_avg[ipb], chl, solz, senz, delphi, rhown_nir);

            for (ib = aer_s; ib <= aer_l; ib += daer) {

                /* Convert NIR reflectance to TOA W-L radiance */
                tLw_nir[ib] = rhown_nir[ib] / pi * Fo[ib] * mu0 * t_sol[ib] * t_sen[ib] / brdf[ib];

                /* Iteration damping */
                tLw_nir[ib] = (1.0 - df) * tLw_nir[ib] + df * last_tLw_nir[ib];

                /* Ramp-up ?*/
                if (want_ramp) {
                    if (chl > 0.0 && chl <= cbot)
                        tLw_nir[ib] = 0.0;
                    else if ((chl > cbot) && (chl < ctop))
                        tLw_nir[ib] *= (cslp * chl + cint);
                }

                /* Remove estimated NIR water-leaving radiance */
                tLw[ib] -= tLw_nir[ib];
            }
        }

        /*  Compute the aerosol contribution */
        if (status == 0) {
            if (aer_opt != AERNULL)
                status = aerosol(l2rec, aer_opt, aerec, ip, wave, nwave, aer_s, aer_l, Fo, tLw,
                    La, t_sol, t_sen, eps, taua, aermodmin, aermodmax, aermodrat,
                    aermodmin2, aermodmax2, aermodrat2);
            else {
                for (ib = 0; ib < nwave; ib++) {

                    ipb = ip * nwave + ib;

                    t_sol [ib] = 1.0;
                    t_sen [ib] = 1.0;
                    La [ib] = 0.0;
                    taua [ib] = 0.0;
                    *eps = 1.0;
                }
            }
        }


        /* Subtract aerosols and normalize */
        if (status == 0) {

            if (aer_opt == AERSSP) {
              tLw[aer_s] = laer1;
              tLw[aer_l] = laer2;
            }

            if (aer_opt == AERSS14)
              for (ib = 0; ib < nwave; ib++) {
                La[ib] = tLw[ib];
                tLw[ib] = Ltemp[ib];
              }


            for (ib = 0; ib < nwave; ib++) {

                /* subtract aerosol and normalize */
                tLw[ib] = tLw[ib] - La[ib];
                Lw [ib] = tLw[ib] / t_sen[ib] * tg_sol[ib];
                nLw[ib] = Lw [ib] / t_sol[ib] / tg_sol[ib] / mu0 / fsol * brdf[ib];
            }

            /* Compute new estimated chlorophyll */
            if (want_nirLw) {
                refl_nir = Rrs[red];
                for (ib = aer_s; ib <= aer_l; ib += daer)
                    last_tLw_nir[ib] = tLw_nir[ib];
                for (ib = 0; ib < nwvis; ib++) {
                    Rrs[ib] = nLw[ib] / Fobar[ib];
                }
                chl = get_default_chl(l2rec, Rrs);

                // if we passed atmospheric correction but the spectral distribution of
                // Rrs is bogus (chl failed), assume this is a turbid-water case and
                // reseed iteration as if all 670 reflectance is from water.

                if (chl == badchl && iter_reset == 0 && iter < iter_max) {
                    chl = 10.0;
                    Rrs[red] = 1.0 * (Ltemp[red] - TLg[red]) / t_sol[red] / tg_sol[red] / mu0 / fsol / Fobar[red];
                    iter_reset = 1;
                }

                // if we already tried a reset, and still no convergence, force one last
                // pass with an assumption that all red radiance is water component, and
                // force iteration to end.  this will be flagged as atmospheric correction
                // failure, but a qualitatively useful retrieval may still result.

                if (chl == badchl && iter_reset == 1 && iter < iter_max) {
                    chl = 10.0;
                    Rrs[red] = 1.0 * (Ltemp[red] - TLg[red]) / t_sol[red] / tg_sol[red] / mu0 / fsol / Fobar[red];
                    iter = iter_max; // so iter will trigger maxiter flag and ATMFAIL
                    iter_reset = 2;
                }
            }

        } else {

            /* Aerosol determination failure */
            for (ib = 0; ib < nwave; ib++) {
                La [ib] = badval;
                tLw[ib] = badval;
                Lw[ib] = badval;
                nLw[ib] = badval;
                Rrs[ib] = badval;
            }
        }

        /* If NIR/SWIR switching, do secondary test for turbidity and reset if needed */
        if (iter == 1 && (aer_opt == AERWANGSWIR || aer_opt == AERRHSWIR)) {
            if (tindx >= 1.05 && status == 0) {
                for (ib = 0; ib < nwvis; ib++) {
                    Rrs[ib] = nLw[ib] / Fobar[ib];
                }
                chl = get_default_chl(l2rec, Rrs);
                //printf("Checking turbidity %d %f %f %f\n",ip,tindx,chl,nLw[nir_l]);
                if (chl > 0 && (chl <= 1.0 || nLw[nir_l] < 0.08)) {
                    iter_max = aer_iter_max;
                    aer_s = nir_s;
                    aer_l = nir_l;
                    daer = MAX(aer_l - aer_s, 1);
                    want_nirLw = 1;
                    //printf("Reverting to NIR %d %f %f %f\n",ip,tindx,chl,nLw[nir_l]);
                    goto NIRSWIR;
                } else
                    l1rec->flags[ip] |= TURBIDW;
            }
        }


        /* Shall we continue iterating */
        if (status != 0) {
            last_iter = 1;
        } else if (iter < iter_min) {
            last_iter = 0;
        } else if (want_nirLw && (fabs(refl_nir - last_refl_nir) < fabs(nir_chg * refl_nir) || refl_nir < 0.0)) {
            last_iter = 1;
        } else if (want_mumm || want_nirRrs) {
            last_iter = 1;
        } else if (iter > iter_max) {
            last_iter = 1;
        }

        last_refl_nir = refl_nir;

    } /* end of iteration loop */

    l2rec->num_iter[ip] = iter;


    /* If the atmospheric correction failed, we don't need to do more. */
    if (status != 0) {
        free(taur);
        free(tLw);
        free(rhown_nir);
        free(tLw_nir);
        free(last_tLw_nir);
        free(Ltemp);
        free(Lunc);
        return (status);

    }

    /* If we used a NIR Lw correction, we record the tLw as it was predicted. */
    if (want_nirLw || want_nirRrs || want_mumm) {
        for (ib = aer_s; ib <= aer_l; ib += daer) {
            tLw[ib] = tLw_nir[ib];
            Lw [ib] = tLw[ib] / t_sen[ib] * tg_sol[ib];
            nLw[ib] = Lw [ib] / t_sol[ib] / tg_sol[ib] / mu0 / fsol * brdf[ib];
            Rrs[ib] = nLw[ib] / Fobar[ib];
        }
    }

    /* Convert water-leaving radiances from band-averaged to band-centered.  */
    /* Switch mean solar irradiance from band-averaged to band centered also.*/
    if (input->outband_opt >= 2) {
        nlw_outband(input->evalmask, sensorID, wave, nwave, Lw, nLw);
        Fobar = l1file->Fonom;
    }

    /* Compute f/Q correction and apply to nLw */
    if (brdf_opt >= FOQMOREL) {
        ocbrdf(l2rec, ip, brdf_opt, wave, nwvis, solz, senz, delphi, ws, -1.0, nLw, Fobar, brdf);
        for (ib = 0; ib < nwvis; ib++) {
            nLw[ib] = nLw[ib] * brdf[ib];
        }
    }

    /* Compute final Rrs */
    if (want_ssp) {
      for (ib = 0; ib < nwave; ib++) {
          Rrs[ib] = nLw[ib] / Fobar[ib];
          l2rec->Rrs[ip * nwave + ib] = Rrs[ib];

          nLw_unc[ib] = Lunc[ib] / t_sen[ib] / t_sol[ib] / mu0 / fsol * brdf[ib];
          Rrs_unc[ib] = nLw_unc[ib] / Fobar[ib];
      }
    }
    else {
      for (ib = 0; ib < nwave; ib++) {
          if (ib != aer_s && ib != aer_l) {
              Rrs[ib] = nLw[ib] / Fobar[ib];
              l2rec->Rrs[ip * nwave + ib] = Rrs[ib];

              nLw_unc[ib] = Lunc[ib] / t_sen[ib] / t_sol[ib] / mu0 / fsol * brdf[ib];
              Rrs_unc[ib] = nLw_unc[ib] / Fobar[ib];
            }
        }
  }

    /* Compute final chl from final nLw (needed for flagging) */
    l2rec->chl[ip] = get_default_chl(l2rec, Rrs);

    /*Determine Raman scattering contribution to Rrs*/
    run_raman_cor(l2rec, ip);

    free(taur);
    free(tLw);
    free(rhown_nir);
    free(tLw_nir);
    free(last_tLw_nir);
    free(Ltemp);
    free(Lunc);

    return (status);
}
