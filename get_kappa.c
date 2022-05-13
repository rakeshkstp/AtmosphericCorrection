#include "l12_proto.h"

float calc_kappa(int32_t sensorID, float rfr[]);

/* --------------------------------------------------------------------------------------- */
/* get_kappa() - get kappa values used in SSP algorithm                                    */
/* Rakesh Kumar Singh                                          		                         */
/* Version: 2020.02.13                                                                     */
/* --------------------------------------------------------------------------------------- */

void get_kappa(l2str *l2rec, float kappa[])
{ l1str *l1rec = l2rec->l1rec;
  filehandle *l1file = l1rec->l1file;
  int32_t sensorID = l1file->sensorID;
  int32_t ip, ipb;
  static float pi = PI;
  int ib;
  int nwave=l1file->nbands;
  float rfr[nwave];

  for (ip = 0; ip < l1rec->npix; ip++)  /*Converting and correcting rho_t to rho_s*/
  { for (ib=0; ib<nwave; ib++)
    {   ipb = ip * nwave + ib;
        rfr[ib] =l2rec->l1rec->Lt[ipb];
        rfr[ib] /= (l1rec->tg_sol[ipb] * l1rec->tg_sen[ipb]);
        rfr[ib] += l1rec->radcor[ipb];
        rfr[ib] /= l1rec->polcor[ipb];
        rfr[ib] -= l1rec->tLf[ipb];
        rfr[ib] -= l1rec->Lr[ipb];
        rfr[ib] *= (pi/l1file->Fobar[ib]);
    }
    kappa[ip] = calc_kappa(sensorID,rfr);
  }
}

/* --------------------------------------------------------------------------------------- */
/* calc_kappa() - calculate kappa values used in SSP algorithm using band ratios           */
/* Rakesh Kumar Singh                                          		                         */
/* Version: 2020.02.13                                                                     */
/* --------------------------------------------------------------------------------------- */

float calc_kappa(int32_t sensorID, float rfr[])
{  float k;

    if (sensorID == GOCI)
    {	k= 1.0;								//clear water
      if (k<rfr[3]/rfr[2])						//suspended sediment
        k = (rfr[3]/rfr[2]);
      if(k<rfr[5]/rfr[4])						//in-water bloom
        k = (rfr[5]/rfr[4]);
      if(k<rfr[6]/rfr[4])						//floating bloom
        k = (rfr[6]/rfr[4]);
      if(k<rfr[4]/rfr[0])		        //Extremely turbid waters
        k = (rfr[4]/rfr[0]);
    }

    if (sensorID == MODISA || sensorID == MODIST)
    {	k=1.0;								//clear water
      if (k<rfr[4]/rfr[3])						//suspended sediment
        k = (rfr[4]/rfr[3]);
      if(k<rfr[9]/rfr[8])						//in-water bloom
        k = (rfr[9]/rfr[8]);
      if(k<rfr[10]/rfr[8])						//floating bloom
        k = (rfr[10]/rfr[8]);
      if(k<rfr[7]/rfr[2] || k > 20)		        //Extremely turbid waters
        k = (rfr[7]/rfr[2]);
      if (k<1)
        k=1;
    }

    if (sensorID == MSIS2A || sensorID == MSIS2B)
    {	k=1.0;								//clear water
      if (k<rfr[2]/rfr[1])						//suspended sediment
        k = (rfr[2]/rfr[1]);
      if(k<rfr[4]/rfr[3])						//floating bloom
        k = (rfr[4]/rfr[3]);
      if(k<rfr[3]/rfr[0])		        //Extremely turbid waters
        k = (rfr[3]/rfr[0]);
    }

    if (sensorID == MERIS)
    {	k= 1.0;								//clear water
      if (k<rfr[3]/rfr[2])						//suspended sediment
        k = (rfr[3]/rfr[2]);
      if(k<rfr[7]/rfr[6])						//in-water bloom
        k = (rfr[7]/rfr[6]);
      if(k<rfr[9]/rfr[6])						//floating bloom
        k = (rfr[9]/rfr[6]);
      if(k<rfr[6]/rfr[0])		        //Extremely turbid waters
        k = (rfr[6]/rfr[0]);
    }

    if (sensorID == OCM2)
    {	k= 1.0;								//clear water
      if (k<rfr[3]/rfr[2])						//suspended sediment
        k = (rfr[3]/rfr[2]);
      if(k<rfr[6]/rfr[5])						//floating bloom
        k = (rfr[6]/rfr[5]);
      if(k<rfr[5]/rfr[0])		        //Extremely turbid waters
        k = (rfr[5]/rfr[0]);
    }

    if (sensorID == OLCIS3A || sensorID == OLCIS3B)
    {	k= 1.0;								//clear water
      if (k<rfr[4]/rfr[3])						//suspended sediment
        k = (rfr[4]/rfr[3]);
      if(k<rfr[9]/rfr[7])						//in-water bloom
        k = (rfr[9]/rfr[7]);
      if(k<rfr[11]/rfr[7])						//floating bloom
        k = (rfr[11]/rfr[7]);
      if(k<rfr[7]/rfr[1])		        //Extremely turbid waters
        k = (rfr[7]/rfr[1]);
    }

    if (sensorID == OLI)
    {	k= 1.0;								//clear water
      if (k<rfr[2]/rfr[1])						//suspended sediment
        k = (rfr[2]/rfr[1]);
      if(k<rfr[4]/rfr[3])						//floating blooms
        k = (rfr[4]/rfr[3]);
      if(k<rfr[3]/rfr[0])						//Extremely turbid waters
        k = (rfr[3]/rfr[0]);
      }

    if (sensorID == L7ETMP || sensorID == L5TM)
    {	k= 1.0;								//clear water
      if (k<rfr[1]/rfr[0])						//suspended sediment
        k = (rfr[1]/rfr[0]);
      if(k<rfr[3]/rfr[2])						//floating bloom
        k = (rfr[3]/rfr[2]);
      if(k<rfr[3]/rfr[0])						//Extremely turbid waters
        k = (rfr[3]/rfr[0]);
      }

    if (sensorID == VIIRSN || sensorID == VIIRSJ1)
    {	k= 1.0;								//clear water
      if (k<rfr[3]/rfr[2])						//suspended sediment
        k = (rfr[3]/rfr[2]);
      if(k<rfr[5]/rfr[4])						//floating bloom
        k = (rfr[5]/rfr[4]);
      if(k<rfr[4]/rfr[0])		        //Extremely turbid waters
        k = (rfr[4]/rfr[0]);
    }

    if (sensorID == HICO)
    {	k=1.0;								//clear water
      if (k<rfr[32]/rfr[24])						//suspended sediment
        k = (rfr[32]/rfr[24]);
      if(k<rfr[58]/rfr[55])						//in-water bloom
        k = (rfr[58]/rfr[55]);
      if(k<rfr[69]/rfr[55])						//floating bloom
        k = (rfr[69]/rfr[55]);
      if(k<rfr[55]/rfr[11])		        //Extremely turbid waters
        k = (rfr[55]/rfr[11]);
    }
  return(k);
}
