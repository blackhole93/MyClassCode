  /** @file cl.c Documented spectra module
 *
 * Julien Lesgourgues, 25.08.2010
 *
 * This module computes the anisotropy and Fourier power spectra
 * \f$ C_l^{X}, P(k), ... \f$'s given the transfer and Bessel functions
 * (for anisotropy spectra), the source functions (for Fourier spectra)
 * and the primordial spectra.
 *
 * The following functions can be called from other modules:
 *
 * -# spectra_init() at the beginning (but after transfer_init())
 * -# spectra_cl_at_l() at any time for computing C at any l
 * -# spectra_spectrum_at_z() at any time for computing P(k) at any z
 * -# spectra_spectrum_at_k_and z() at any time for computing P at any k and z
 * -# spectra_free() at the end
 */

#include "spectra.h"



int spectra_bandpower(struct spectra * psp,
                      int l1,
                      int l2,
                      double * TT_II,
                      double * TT_RI,
                      double * TT_RR
                      ) {

  int l;
  int index_md;
  double * cl_tot;
  double ** cl_md;
  double ** cl_md_ic;

  class_alloc(cl_tot,psp->ct_size*sizeof(double),psp->error_message);
  class_alloc(cl_md,psp->md_size*sizeof(double*),psp->error_message);
  class_alloc(cl_md_ic,psp->md_size*sizeof(double*),psp->error_message);
  for (index_md=0;index_md<psp->md_size; index_md++) {
    class_alloc(cl_md[index_md],psp->ct_size*sizeof(double),psp->error_message);
    class_alloc(cl_md_ic[index_md],psp->ct_size*psp->ic_ic_size[index_md]*sizeof(double),psp->error_message);
  }

  *TT_RR=0.;
  *TT_RI=0.;
  *TT_II=0.;

  for (l=l1; l<=l2; l++) {

    class_call(spectra_cl_at_l(psp,
                               (double)l,
                               cl_tot,
                               cl_md,
                               cl_md_ic),
               psp->error_message,
               psp->error_message);

    *TT_RR += (double)(2*l+1)*cl_md_ic[psp->index_md_scalars][index_symmetric_matrix(0,0,psp->ic_size[psp->index_md_scalars])*psp->ct_size+psp->index_ct_tt];
    *TT_RI += (double)(2*l+1)*cl_md_ic[psp->index_md_scalars][index_symmetric_matrix(0,1,psp->ic_size[psp->index_md_scalars])*psp->ct_size+psp->index_ct_tt]*2.;
    *TT_II += (double)(2*l+1)*cl_md_ic[psp->index_md_scalars][index_symmetric_matrix(1,1,psp->ic_size[psp->index_md_scalars])*psp->ct_size+psp->index_ct_tt];

  }

  for (index_md=0;index_md<psp->md_size; index_md++) {
    free(cl_md[index_md]);
    free(cl_md_ic[index_md]);
  }
  free(cl_tot);
  free(cl_md);
  free(cl_md_ic);

  return _SUCCESS_;

}

/**
 * Anisotropy power spectra C_l's for all types, modes and initial conditions.
 *
 * This routine evaluates all the C_l's at a given value of l by
 * interpolating in the pre-computed table. When relevant, it also
 * sums over all initial conditions for each mode, and over all modes.
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
 *
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param l          Input: multipole number
 * @param cl_tot     Ouput: total C_l's for all types (TT, TE, EE, etc..)
 * @param cl_md      Ouput: C_l's for all types (TT, TE, EE, etc..) decomposed mode by mode (scalar, tensor, ...) when relevant
 * @param cl_md_ic   Ouput: C_l's for all types (TT, TE, EE, etc..) decomposed by pairs of initial conditions (adiabatic, isocurvatures) for each mode (usually, only for the scalar mode) when relevant
 * @return the error status
 */

int spectra_cl_at_l(
                    struct spectra * psp,
                    double l,
                    double * cl_tot,    /* array with argument cl_tot[index_ct] (must be already allocated) */
                    double * * cl_md,   /* array with argument cl_md[index_md][index_ct] (must be already allocated only if several modes) */
                    double * * cl_md_ic /* array with argument cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct] (must be already allocated for a given mode only if several ic's) */
                    ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  int index_md;
  int index_ic1,index_ic2,index_ic1_ic2;
  int index_ct;

  /** A) treat case in which there is only one mode and one initial condition.
      Then, only cl_tot needs to be filled. */

  if ((psp->md_size == 1) && (psp->ic_size[0] == 1)) {
    index_md = 0;
    if ((int)l <= psp->l[psp->l_size[index_md]-1]) {

      /* interpolate at l */
      class_call(array_interpolate_spline(psp->l,
                                          psp->l_size[index_md],
                                          psp->cl[index_md],
                                          psp->ddcl[index_md],
                                          psp->ct_size,
                                          l,
                                          &last_index,
                                          cl_tot,
                                          psp->ct_size,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      /* set to zero for the types such that l<l_max */
      for (index_ct=0; index_ct<psp->ct_size; index_ct++)
        if ((int)l > psp->l_max_ct[index_md][index_ct])
          cl_tot[index_ct]=0.;
    }
    else {
      for (index_ct=0; index_ct<psp->ct_size; index_ct++)
        cl_tot[index_ct]=0.;
    }
  }

  /** B) treat case in which there is only one mode
      with several initial condition.
      Fill cl_md_ic[index_md=0] and sum it to get cl_tot. */

  if ((psp->md_size == 1) && (psp->ic_size[0] > 1)) {
    index_md = 0;
    for (index_ct=0; index_ct<psp->ct_size; index_ct++)
      cl_tot[index_ct]=0.;
    for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);
        if (((int)l <= psp->l[psp->l_size[index_md]-1]) &&
            (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_)) {

          class_call(array_interpolate_spline(psp->l,
                                              psp->l_size[index_md],
                                              psp->cl[index_md],
                                              psp->ddcl[index_md],
                                              psp->ic_ic_size[index_md]*psp->ct_size,
                                              l,
                                              &last_index,
                                              cl_md_ic[index_md],
                                              psp->ic_ic_size[index_md]*psp->ct_size,
                                              psp->error_message),
                     psp->error_message,
                     psp->error_message);

          for (index_ct=0; index_ct<psp->ct_size; index_ct++)
            if ((int)l > psp->l_max_ct[index_md][index_ct])
              cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct]=0.;
        }
        else {
          for (index_ct=0; index_ct<psp->ct_size; index_ct++)
            cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct]=0.;
        }

        /* compute cl_tot by summing over cl_md_ic */
        for (index_ct=0; index_ct<psp->ct_size; index_ct++) {
          if (index_ic1 == index_ic2)
            cl_tot[index_ct]+=cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct];
          else
            cl_tot[index_ct]+=2.*cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct];
        }
      }
    }
  }

  /** C) loop over modes */

  if (psp->md_size > 1) {

    for (index_ct=0; index_ct<psp->ct_size; index_ct++)
      cl_tot[index_ct]=0.;

    for (index_md = 0; index_md < psp->md_size; index_md++) {

      /** C.1) treat case in which the mode under consideration
          has only one initial condition.
          Fill cl_md[index_md]. */

      if (psp->ic_size[index_md] == 1) {
        if ((int)l <= psp->l[psp->l_size[index_md]-1]) {

          class_call(array_interpolate_spline(psp->l,
                                              psp->l_size[index_md],
                                              psp->cl[index_md],
                                              psp->ddcl[index_md],
                                              psp->ct_size,
                                              l,
                                              &last_index,
                                              cl_md[index_md],
                                              psp->ct_size,
                                              psp->error_message),
                     psp->error_message,
                     psp->error_message);

          for (index_ct=0; index_ct<psp->ct_size; index_ct++)
            if ((int)l > psp->l_max_ct[index_md][index_ct])
              cl_md[index_md][index_ct]=0.;
        }
        else {
          for (index_ct=0; index_ct<psp->ct_size; index_ct++)
            cl_md[index_md][index_ct]=0.;
        }
      }

      /** C.2) treat case in which the mode under consideration
          has several initial conditions.
          Fill cl_md_ic[index_md] and sum it to get cl_md[index_md] */

      if (psp->ic_size[index_md] > 1) {

        if ((int)l <= psp->l[psp->l_size[index_md]-1]) {

          /* interpolate all ic and ct */
          class_call(array_interpolate_spline(psp->l,
                                              psp->l_size[index_md],
                                              psp->cl[index_md],
                                              psp->ddcl[index_md],
                                              psp->ic_ic_size[index_md]*psp->ct_size,
                                              l,
                                              &last_index,
                                              cl_md_ic[index_md],
                                              psp->ic_ic_size[index_md]*psp->ct_size,
                                              psp->error_message),
                     psp->error_message,
                     psp->error_message);

          /* set to zero some of the components */
          for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
            for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
              index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);
              for (index_ct=0; index_ct<psp->ct_size; index_ct++) {

                if (((int)l > psp->l_max_ct[index_md][index_ct]) || (psp->is_non_zero[index_md][index_ic1_ic2] == _FALSE_))
                  cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct]=0.;
              }
            }
          }
        }
        /* if l was too big, set anyway all components to zero */
        else {
          for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
            for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
              index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);
              for (index_ct=0; index_ct<psp->ct_size; index_ct++) {
                cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct]=0.;
              }
            }
          }
        }

        /* sum up all ic for each mode */

        for (index_ct=0; index_ct<psp->ct_size; index_ct++) {

          cl_md[index_md][index_ct]=0.;

          for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
            for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
              index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);

              if (index_ic1 == index_ic2)
                cl_md[index_md][index_ct]+=cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct];
              else
                cl_md[index_md][index_ct]+=2.*cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct];
            }
          }
        }
      }

      /** C.3) add contribution of cl_md[index_md] to cl_tot */

      for (index_ct=0; index_ct<psp->ct_size; index_ct++)
        cl_tot[index_ct]+=cl_md[index_md][index_ct];
    }
  }

  return _SUCCESS_;

}

/**
 * Matter power spectrum for arbitrary redshift and for all initial conditions.
 *
 * This routine evaluates the matter power spectrum at a given value of z by
 * interpolating in the pre-computed table (if several values of z have been stored)
 * or by directly reading it (if it only contains values at z=0 and we want P(k,z=0))
 *
 *
 * Can be called in two modes: linear or logarithmic.
 *
 * - linear: returns P(k) (units: Mpc^3)
 *
 * - logarithmic: returns ln(P(k))
 *
 * One little subtlety: in case of several correlated initial conditions,
 * the cross-correlation spectrum can be negative. Then, in logarithmic mode,
 * the non-diagonal elements contain the cross-correlation angle P_12/sqrt(P_11 P_22)
 * (from -1 to 1) instead of ln(P_12)
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param mode       Input: linear or logarithmic
 * @param z          Input: redshift
 * @param output_tot Ouput: total matter power spectrum P(k) in Mpc**3 (linear mode), or its logarithms (logarithmic mode)
 * @param output_ic  Ouput: for each pair of initial conditions, matter power spectra P(k) in Mpc**3 (linear mode), or their logarithms and cross-correlation angles (logarithmic mode)
 * @return the error status
 */

int spectra_pk_at_z(
                    struct background * pba,
                    struct spectra * psp,
                    enum linear_or_logarithmic mode,
                    double z,
                    double * output_tot, /* array with argument output_tot[index_k] (must be already allocated) */
                    double * output_ic   /* array with argument output_tot[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] (must be already allocated only if more than one initial condition) */
                    ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int last_index;
  int index_k;
  double tau,ln_tau;
  int index_ic1,index_ic2,index_ic1_ic2;

  index_md = psp->index_md_scalars;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             psp->error_message);

  class_test(tau <= 0.,
             psp->error_message,
             "negative or null value of conformal time: cannot interpolate");

  ln_tau = log(tau);

  /** - second step: for both modes (linear or logarithmic), store the spectrum in logarithmic format in the output array(s) */

  /**   (a.) if only values at tau=tau_today are stored and we want P(k,z=0), no need to interpolate */

  if (psp->ln_tau_size == 1) {

    class_test(z != 0.,
               psp->error_message,
               "asked z=%e but only P(k,z=0) has been tabulated",z);

    for (index_k=0; index_k<psp->ln_k_size; index_k++)
      if (psp->ic_size[index_md] == 1) {
      	output_tot[index_k] = psp->ln_pk[index_k];
      }
      else {
        for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_md]; index_ic1_ic2++) {
          output_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] =
            psp->ln_pk[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2];
        }
      }
  }

  /**   (b.) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else {

    if (psp->ic_ic_size[index_md] == 1) {

      class_call(array_interpolate_spline(psp->ln_tau,
                                          psp->ln_tau_size,
                                          psp->ln_pk,
                                          psp->ddln_pk,
                                          psp->ln_k_size,
                                          ln_tau,
                                          &last_index,
                                          output_tot,
                                          psp->ln_k_size,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

    }
    else {

      class_call(array_interpolate_spline(psp->ln_tau,
                                          psp->ln_tau_size,
                                          psp->ln_pk,
                                          psp->ddln_pk,
                                          psp->ic_ic_size[index_md]*psp->ln_k_size,
                                          ln_tau,
                                          &last_index,
                                          output_ic,
                                          psp->ic_ic_size[index_md]*psp->ln_k_size,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);
    }
  }

  /** - third step: if there are several initial conditions, compute the total P(k) and set back all uncorrelated coefficients to exactly zero. Check positivity of total P(k). */

  if (psp->ic_size[index_md] > 1) {
    for (index_k=0; index_k<psp->ln_k_size; index_k++) {
      output_tot[index_k] = 0.;
      for (index_ic1=0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
        for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
          index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);
          if (index_ic1 == index_ic2) {
            output_tot[index_k] += exp(output_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2]);
          }
          else {
            if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {
              output_tot[index_k] +=
                2. * output_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] *
                sqrt(exp(output_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md])]) *
                     exp(output_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_md])]));
            }
            else
              output_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] = 0.;
          }
        }
      }

      class_test(output_tot[index_k] <= 0.,
                 psp->error_message,
                 "for k=%e, z=%e, the matrix of initial condition amplitudes was not positive definite, hence P(k)_total=%e results negative",
                 exp(psp->ln_k[index_k]),z,output_tot[index_k]);

    }
  }

  /** - fourth step: depending on requested mode (linear or logarithmic), apply necessary transformation to the output arrays */

  /**   (a.) linear mode: if only one initial condition, convert output_pk to linear format; if several initial conditions, convert output_ic to linear format, output_tot is already in this format */

  if (mode == linear) {

    if (psp->ic_size[index_md] == 1) {
      for (index_k=0; index_k<psp->ln_k_size; index_k++) {
        output_tot[index_k] = exp(output_tot[index_k]);
      }
    }

    else {
      for (index_k=0; index_k<psp->ln_k_size; index_k++) {
        for (index_ic1=0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
          index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md]);
          output_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2] = exp(output_ic[index_k * psp->ic_ic_size[index_md] + index_ic1_ic2]);
        }
        for (index_ic1=0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
          for (index_ic2 = index_ic1+1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {

            output_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md])] =
              output_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md])]
              *sqrt(output_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md])] *
                    output_ic[index_k * psp->ic_ic_size[index_md] + index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_md])]);
          }
        }
      }
    }
  }

  /**   (b.) logarithmic mode: if only one initial condition, nothing to be done; if several initial conditions, convert output_tot to logarithmic format, output_ic is already in this format */

  else {

    if (psp->ic_size[index_md] > 1) {
      for (index_k=0; index_k<psp->ln_k_size; index_k++) {
        /* we have already checked above that output_tot was positive */
        output_tot[index_k] = log(output_tot[index_k]);
      }
    }
  }

  return _SUCCESS_;

}

/**
 * Matter power spectrum for arbitrary wavenumber, redshift and initial condition.
 *
 * This routine evaluates the matter power spectrum at a given value of k and z by
 * interpolating in a table of all P(k)'s computed at this z by spectra_pk_at_z() (when kmin <= k <= kmax),
 * or eventually by using directly the primordial spectrum (when 0 <= k < kmin):
 * the latter case is an approximation, valid when kmin << comoving Hubble scale today.
 * Returns zero when k=0. Returns an error when k<0 or k > kmax.
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param ppm        Input: pointer to primordial structure (used only in the case 0 < k < kmin)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param k          Input: wavenumber in 1/Mpc
 * @param z          Input: redshift
 * @param pk_tot     Ouput: total matter power spectrum P(k) in Mpc**3
 * @param pk_ic      Ouput: for each pair of initial conditions, matter power spectra P(k) in Mpc**3
 * @return the error status
 */

int spectra_pk_at_k_and_z(
                          struct background * pba,
                          struct primordial * ppm,
                          struct spectra * psp,
                          double k,
                          double z,
                          double * pk_tot, /* pointer to a single number (must be already allocated) */
                          double * pk_ic   /* array of argument pk_ic[index_ic1_ic2] (must be already allocated only if several initial conditions) */
                          ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_k;
  int last_index;
  int index_ic1,index_ic2,index_ic1_ic2;

  double * spectrum_at_z = NULL;
  double * spectrum_at_z_ic = NULL;
  double * spline;
  double * pk_primordial_k = NULL;
  double kmin;
  double * pk_primordial_kmin = NULL;

  index_md = psp->index_md_scalars;

  /** - first step: check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_pk_at_z()) */

  class_test((k < 0.) || (k > exp(psp->ln_k[psp->ln_k_size-1])),
             psp->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,exp(psp->ln_k[psp->ln_k_size-1]));

  /** - deal with case 0 <= k < kmin */

  if (k < exp(psp->ln_k[0])) {

    /**   (a.) subcase k=0: then P(k)=0 */

    if (k == 0.) {
      if (psp->ic_size[index_md] == 1) {
        *pk_tot=0.;
      }
      else {
        for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_md]; index_ic1_ic2++) {
          pk_ic[index_ic1_ic2] = 0.;
        }
      }
    }

    /**    (b.) subcase 0<k<kmin: in this case we know that on super-Hubble scales:
     *          P(k) = [some number] * k  * P_primordial(k)
     *          so
     *          P(k) = P(kmin) * (k P_primordial(k)) / (kmin P_primordial(kmin))
     *          (note that the result is accurate only if kmin is such that [a0 kmin] << H0)
     */

    else {

      /* compute P(k,z) which contains P(kmin,z)*/
      class_alloc(spectrum_at_z,
                  psp->ln_k_size*sizeof(double),
                  psp->error_message);
      if (psp->ic_size[index_md] > 1) {
        class_alloc(spectrum_at_z_ic,
                    sizeof(double)*psp->ic_ic_size[index_md]*psp->ln_k_size,
                    psp->error_message);
      }
      class_call(spectra_pk_at_z(pba,
                                 psp,
                                 linear,
                                 z,
                                 spectrum_at_z,
                                 spectrum_at_z_ic),
                 psp->error_message,
                 psp->error_message);

      /* compute P_primordial(k) */
      class_alloc(pk_primordial_k,
                  sizeof(double)*psp->ic_ic_size[index_md],
                  psp->error_message);
      class_call(primordial_spectrum_at_k(ppm,
                                          index_md,
                                          linear,
                                          k,
                                          pk_primordial_k),
                 ppm->error_message,psp->error_message);

      /* compute P_primordial(kmin) */
      kmin = exp(psp->ln_k[0]);
      class_alloc(pk_primordial_kmin,
                  sizeof(double)*psp->ic_ic_size[index_md],
                  psp->error_message);
      class_call(primordial_spectrum_at_k(ppm,
                                          index_md,
                                          linear,
                                          kmin,
                                          pk_primordial_kmin),
                 ppm->error_message,
                 psp->error_message);

      /* apply above analytic approximation for P(k) */
      index_k=0;
      if (psp->ic_size[index_md] == 1) {
        index_ic1_ic2 = 0;
        *pk_tot = spectrum_at_z[index_k]
          *k*pk_primordial_k[index_ic1_ic2]
          /kmin/pk_primordial_kmin[index_ic1_ic2];
      }
      else {
      	for (index_ic1_ic2 = 0; index_ic1_ic2 < psp->ic_ic_size[index_md]; index_ic1_ic2++) {
          pk_ic[index_ic1_ic2] = spectrum_at_z_ic[index_ic1_ic2]
            *k*pk_primordial_k[index_ic1_ic2]
            /kmin/pk_primordial_kmin[index_ic1_ic2];
        }
      }

      free(spectrum_at_z);
      if (psp->ic_size[index_md] > 1)
        free(spectrum_at_z_ic);
      free(pk_primordial_k);
      free(pk_primordial_kmin);

    }
  }

  /** - deal with case kmin <= k <= kmax */

  else {

    /* compute P(k,z) (in logarithmic format for more accurate interpolation) */
    class_alloc(spectrum_at_z,
                psp->ln_k_size*sizeof(double),
                psp->error_message);
    if (psp->ic_size[index_md] > 1) {
      class_alloc(spectrum_at_z_ic,
                  sizeof(double)*psp->ic_ic_size[index_md]*psp->ln_k_size,
                  psp->error_message);
    }
    class_call(spectra_pk_at_z(pba,
                               psp,
                               logarithmic,
                               z,
                               spectrum_at_z,
                               spectrum_at_z_ic),
               psp->error_message,
               psp->error_message);

    /* get its second derivatives with spline, then interpolate, then convert to linear format */

    class_alloc(spline,
                sizeof(double)*psp->ic_ic_size[index_md]*psp->ln_k_size,
                psp->error_message);

    if (psp->ic_size[index_md] == 1) {

      class_call(array_spline_table_lines(psp->ln_k,
                                          psp->ln_k_size,
                                          spectrum_at_z,
                                          1,
                                          spline,
                                          _SPLINE_NATURAL_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      class_call(array_interpolate_spline(psp->ln_k,
                                          psp->ln_k_size,
                                          spectrum_at_z,
                                          spline,
                                          1,
                                          log(k),
                                          &last_index,
                                          pk_tot,
                                          1,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      *pk_tot = exp(*pk_tot);

    }
    else {

      class_call(array_spline_table_lines(psp->ln_k,
                                          psp->ln_k_size,
                                          spectrum_at_z_ic,
                                          psp->ic_ic_size[index_md],
                                          spline,
                                          _SPLINE_NATURAL_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      class_call(array_interpolate_spline(psp->ln_k,
                                          psp->ln_k_size,
                                          spectrum_at_z_ic,
                                          spline,
                                          psp->ic_ic_size[index_md],
                                          log(k),
                                          &last_index,
                                          pk_ic,
                                          psp->ic_ic_size[index_md],
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

      for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md]);
        pk_ic[index_ic1_ic2] = exp(pk_ic[index_ic1_ic2]);
      }
      for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
        for (index_ic2 = index_ic1+1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
          index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);
          if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {
            pk_ic[index_ic1_ic2] = pk_ic[index_ic1_ic2]*
              sqrt(pk_ic[index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md])]*
                   pk_ic[index_symmetric_matrix(index_ic2,index_ic2,psp->ic_size[index_md])]);
          }
          else {
            pk_ic[index_ic1_ic2] = 0.;
          }
        }
      }
      free(spectrum_at_z_ic);
    }

    free(spectrum_at_z);
    free(spline);
  }

  /** - last step: if more than one condition, sum over pk_ic to get pk_tot, and set back coefficients of non-correlated pairs to exactly zero. */

  if (psp->ic_size[index_md] > 1) {

    *pk_tot = 0.;

    for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);

        if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

          if (index_ic1 == index_ic2)
            *pk_tot += pk_ic[index_ic1_ic2];
          else
            *pk_tot += 2.*pk_ic[index_ic1_ic2];
        }
        else {
          pk_ic[index_ic1_ic2] = 0.;
        }
      }
    }

    class_test(*pk_tot <= 0.,
               psp->error_message,
               "for k=%e, the matrix of initial condition amplitudes was not positive definite, hence P(k)_total results negative",k);

  }

  return _SUCCESS_;

}

/**
 * Non-linear total matter power spectrum for arbitrary redshift.
 *
 * This routine evaluates the non-linear matter power spectrum at a given value of z by
 * interpolating in the pre-computed table (if several values of z have been stored)
 * or by directly reading it (if it only contains values at z=0 and we want P(k,z=0))
 *
 *
 * Can be called in two modes: linear or logarithmic.
 *
 * - linear: returns P(k) (units: Mpc^3)
 *
 * - logarithmic: returns ln(P(k))
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param mode       Input: linear or logarithmic
 * @param z          Input: redshift
 * @param output_tot Ouput: total matter power spectrum P(k) in Mpc**3 (linear mode), or its logarithms (logarithmic mode)
 * @return the error status
 */

int spectra_pk_nl_at_z(
                       struct background * pba,
                       struct spectra * psp,
                       enum linear_or_logarithmic mode,
                       double z,
                       double * output_tot /* array with argument output_tot[index_k] (must be already allocated) */
                       ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  int index_k;
  double tau,ln_tau;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             psp->error_message);

  class_test(tau <= 0.,
             psp->error_message,
             "negative or null value of conformal time: cannot interpolate");

  ln_tau = log(tau);

  /** - second step: for both modes (linear or logarithmic), store the spectrum in logarithmic format in the output array(s) */

  /**   (a.) if only values at tau=tau_today are stored and we want P(k,z=0), no need to interpolate */

  if (psp->ln_tau_size == 1) {

    class_test(z != 0.,
               psp->error_message,
               "asked z=%e but only P(k,z=0) has been tabulated",z);

    for (index_k=0; index_k<psp->ln_k_size; index_k++) {
      output_tot[index_k] = psp->ln_pk_nl[index_k];
    }
  }

  /**   (b.) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else {

    class_call(array_interpolate_spline(psp->ln_tau,
                                        psp->ln_tau_size,
                                        psp->ln_pk_nl,
                                        psp->ddln_pk_nl,
                                        psp->ln_k_size,
                                        ln_tau,
                                        &last_index,
                                        output_tot,
                                        psp->ln_k_size,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);
  }

  /** - fourth step: eventually convert to linear format */

  if (mode == linear) {
    for (index_k=0; index_k<psp->ln_k_size; index_k++) {
      output_tot[index_k] = exp(output_tot[index_k]);
    }
  }

  return _SUCCESS_;

}

/**
 * Non-linear total matter power spectrum for arbitrary wavenumber and redshift.
 *
 * This routine evaluates the matter power spectrum at a given value of k and z by
 * interpolating in a table of all P(k)'s computed at this z by spectra_pk_nl_at_z() (when kmin <= k <= kmax),
 * or eventually by using directly the primordial spectrum (when 0 <= k < kmin):
 * the latter case is an approximation, valid when kmin << comoving Hubble scale today.
 * Returns zero when k=0. Returns an error when k<0 or k > kmax.
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param ppm        Input: pointer to primordial structure (used only in the case 0 < k < kmin)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param k          Input: wavenumber in 1/Mpc
 * @param z          Input: redshift
 * @param pk_tot     Ouput: total matter power spectrum P(k) in Mpc**3
 * @return the error status
 */

int spectra_pk_nl_at_k_and_z(
                             struct background * pba,
                             struct primordial * ppm,
			     struct perturbs * ppt,
			     struct nonlinear * pnl,
                             struct spectra * psp,
                             double k,
                             double z,
                             double * pk_tot /* pointer to a single number (must be already allocated) */
                             ) {

  
  if ((pnl->method != nl_none)&&(pnl->method != nl_halomodel)){
  /** Summary: */
  /** - define local variables */

  int index_md;
  int last_index;

  double * spectrum_at_z = NULL;
  double * spline;

  index_md = psp->index_md_scalars;

  /** - first step: check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_pk_at_z()) */

  class_test((k < exp(psp->ln_k[0])) || (k > exp(psp->ln_k[psp->ln_k_size-1])),
             psp->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,exp(psp->ln_k[psp->ln_k_size-1]));

  /* compute P(k,z) (in logarithmic format for more accurate interpolation) */
  class_alloc(spectrum_at_z,
              psp->ln_k_size*sizeof(double),
              psp->error_message);

  class_call(spectra_pk_nl_at_z(pba,
                                psp,
                                logarithmic,
                                z,
                                spectrum_at_z),
             psp->error_message,
             psp->error_message);

  /* get its second derivatives with spline, then interpolate, then convert to linear format */

  class_alloc(spline,
              sizeof(double)*psp->ic_ic_size[index_md]*psp->ln_k_size,
              psp->error_message);

  class_call(array_spline_table_lines(psp->ln_k,
                                      psp->ln_k_size,
                                      spectrum_at_z,
                                      1,
                                      spline,
                                      _SPLINE_NATURAL_,
                                      psp->error_message),
             psp->error_message,
             psp->error_message);

  class_call(array_interpolate_spline(psp->ln_k,
                                      psp->ln_k_size,
                                      spectrum_at_z,
                                      spline,
                                      1,
                                      log(k),
                                      &last_index,
                                      pk_tot,
                                      1,
                                      psp->error_message),
             psp->error_message,
             psp->error_message);

  *pk_tot = exp(*pk_tot);

  free(spectrum_at_z);
  free(spline);
  }
  
  
  if (pnl->method == nl_halomodel){
    double pk;
	class_call(spectra_pk_hm_at_k_and_z(pba,ppm,psp,k,z,&pk),
               psp->error_message,
               psp->error_message);
	
	*pk_tot = pk;
  }
  return _SUCCESS_;

}


/**
 * Matter transfer functions T_i(k) for arbitrary redshift and for all
 * initial conditions.
 *
 * This routine evaluates the matter transfer functions at a given value of z by
 * interpolating in the pre-computed table (if several values of z have been stored)
 * or by directly reading it (if it only contains values at z=0 and we want T_i(k,z=0))
 *
 *
 * This function can be
 * called from whatever module at whatever time, provided that
 * spectra_init() has been called before, and spectra_free() has not
 * been called yet.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param z          Input: redshift
 * @param output     Ouput: matter transfer functions
 * @return the error status
 */

int spectra_tk_at_z(
                    struct background * pba,
                    struct spectra * psp,
                    double z,
                    double * output /* array with argument output[(index_k*psp->ic_size[index_md]+index_ic)*psp->tr_size+index_tr] (must be already allocated) */
                    ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int last_index;
  int index_k;
  int index_tr;
  double tau,ln_tau;
  int index_ic;

  index_md = psp->index_md_scalars;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             psp->error_message);

  class_test(tau <= 0.,
             psp->error_message,
             "negative or null value of conformal time: cannot interpolate");

  ln_tau = log(tau);

  /** - second step: store the matter transfer functions in the output array */

  /**   (a.) if only values at tau=tau_today are stored and we want T_i(k,z=0), no need to interpolate */

  if (psp->ln_tau_size == 1) {

    class_test(z != 0.,
               psp->error_message,
               "asked z=%e but only T_i(k,z=0) has been tabulated",z);

    for (index_k=0; index_k<psp->ln_k_size; index_k++)
      for (index_tr=0; index_tr<psp->tr_size; index_tr++)
        for (index_ic = 0; index_ic < psp->ic_size[index_md]; index_ic++)
          output[(index_k*psp->ic_size[index_md]+index_ic)*psp->tr_size+index_tr]
            = psp->matter_transfer[(index_k*psp->ic_size[index_md]+index_ic)*psp->tr_size+index_tr];

  }

  /**   (b.) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else {

    class_call(array_interpolate_spline(psp->ln_tau,
                                        psp->ln_tau_size,
                                        psp->matter_transfer,
                                        psp->ddmatter_transfer,
                                        psp->ic_size[index_md]*psp->tr_size*psp->ln_k_size,
                                        ln_tau,
                                        &last_index,
                                        output,
                                        psp->ic_size[index_md]*psp->tr_size*psp->ln_k_size,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);

  }

  return _SUCCESS_;

}

/**
 * Matter transfer functions T_i(k) for arbitrary wavenumber, redshift
 * and initial condition.
 *
 * This routine evaluates the matter transfer functions at a given
 * value of k and z by interpolating in a table of all T_i(k,z)'s
 * computed at this z by spectra_tk_at_z() (when kmin <= k <= kmax).
 * Returns an error when k<kmin or k > kmax.
 *
 * This function can be called from whatever module at whatever time,
 * provided that spectra_init() has been called before, and
 * spectra_free() has not been called yet.
 *
 * @param pba        Input: pointer to background structure (used for converting z into tau)
 * @param psp        Input: pointer to spectra structure (containing pre-computed table)
 * @param k          Input: wavenumber in 1/Mpc
 * @param z          Input: redshift
 * @param output     Ouput: matter transfer functions
 * @return the error status
 */

int spectra_tk_at_k_and_z(
                          struct background * pba,
                          struct spectra * psp,
                          double k,
                          double z,
                          double * output  /* array with argument output[index_ic*psp->tr_size+index_tr] (must be already allocated) */
                          ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int last_index;
  double * tks_at_z;
  double * ddtks_at_z;

  index_md = psp->index_md_scalars;

  /** - first step: check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_tk_at_z()) */

  class_test((k < 0.) || (k > exp(psp->ln_k[psp->ln_k_size-1])),
             psp->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,exp(psp->ln_k[psp->ln_k_size-1]));

  /* compute T_i(k,z) */

  class_alloc(tks_at_z,
              psp->ln_k_size*psp->tr_size*psp->ic_size[index_md]*sizeof(double),
              psp->error_message);

  class_call(spectra_tk_at_z(pba,
                             psp,
                             z,
                             tks_at_z),
             psp->error_message,
             psp->error_message);

  /* get its second derivatives w.r.t. k with spline, then interpolate */

  class_alloc(ddtks_at_z,
              psp->ln_k_size*psp->tr_size*psp->ic_size[index_md]*sizeof(double),
              psp->error_message);

  class_call(array_spline_table_lines(psp->ln_k,
                                      psp->ln_k_size,
                                      tks_at_z,
                                      psp->tr_size*psp->ic_size[index_md],
                                      ddtks_at_z,
                                      _SPLINE_NATURAL_,
                                      psp->error_message),
             psp->error_message,
             psp->error_message);

  class_call(array_interpolate_spline(psp->ln_k,
                                      psp->ln_k_size,
                                      tks_at_z,
                                      ddtks_at_z,
                                      psp->tr_size*psp->ic_size[index_md],
                                      log(k),
                                      &last_index,
                                      output,
                                      psp->tr_size*psp->ic_size[index_md],
                                      psp->error_message),
             psp->error_message,
             psp->error_message);

  free(tks_at_z);
  free(ddtks_at_z);

  return _SUCCESS_;

}

/**
 * This routine initializes the spectra structure (in particular,
 * computes table of anisotropy and Fourier spectra \f$ C_l^{X}, P(k), ... \f$)
 *
 * @param ppr Input : pointer to precision structure
 * @param pba Input : pointer to background structure (will provide H, Omega_m at redshift of interest)
 * @param ppt Input : pointer to perturbation structure
 * @param ptr Input : pointer to transfer structure
 * @param ppm Input : pointer to primordial structure
 * @param psp Output: pointer to initialized spectra structure
 * @return the error status
 */

int spectra_init(
                 struct precision * ppr,
                 struct background * pba,
                 struct perturbs * ppt,
                 struct primordial * ppm,
                 struct nonlinear *pnl,
                 struct transfers * ptr,
                 struct spectra * psp
                 ) {

  /** Summary: */

  double TT_II,TT_RI,TT_RR;
  int l1,l2;

  double k;

  /** - check that we really want to compute at least one spectrum */

  if ((ppt->has_cls == _FALSE_) &&
      (ppt->has_pk_matter == _FALSE_) &&
      (ppt->has_density_transfers == _FALSE_) &&
      (ppt->has_velocity_transfers == _FALSE_)) {
    psp->md_size = 0;
    if (psp->spectra_verbose > 0)
      printf("No spectra requested. Spectra module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (psp->spectra_verbose > 0)
      printf("Computing unlensed linear spectra\n");
  }

  /** - initialize indices and allocate some of the arrays in the
      spectra structure */

  class_call(spectra_indices(pba,ppt,ptr,ppm,psp),
             psp->error_message,
             psp->error_message);

  /** - deal with C_l's, if any */

  if (ppt->has_cls == _TRUE_) {

    class_call(spectra_cls(pba,ppt,ptr,ppm,psp),
               psp->error_message,
               psp->error_message);

  }
  else {
    psp->ct_size=0;
  }

  /** - deal with P(k,tau) and T_i(k,tau) */

  if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_)) {

    class_call(spectra_k_and_tau(pba,ppt,psp),
               psp->error_message,
               psp->error_message);

    if (ppt->has_pk_matter == _TRUE_) {

      class_call(spectra_pk(pba,ppt,ppm,pnl,psp),
                 psp->error_message,
                 psp->error_message);

    }
    else {
      psp->ln_pk=NULL;
    }

    if ((ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_)) {

      class_call(spectra_matter_transfers(pba,ppt,psp),
                 psp->error_message,
                 psp->error_message);
    }
    else {
      psp->matter_transfer=NULL;
    }

  }
  else {
    psp->ln_k_size=0;
  }

  /* if there is one isocurvature mode, compute and store in the psp
     structure the isocurvature contribution to some bandpowers in
     different ranges of l, and the contribution to the primordial
     spectrum at different wavenumbers (used in the Planck
     analysis) */

  if ((ppt->has_scalars == _TRUE_) && (ppt->has_cls == _TRUE_) && (ppt->ic_size[ppt->index_md_scalars] == 2)) {

    l1=2;
    l2=20;

    class_call(spectra_bandpower(psp,l1,l2,&TT_II,&TT_RI,&TT_RR),
               psp->error_message,
               psp->error_message);

    class_test(TT_II+TT_RI+TT_RR==0.,
               psp->error_message,
               "should never happen");

    psp->alpha_II_2_20=TT_II/(TT_II+TT_RI+TT_RR);
    psp->alpha_RI_2_20=TT_RI/(TT_II+TT_RI+TT_RR);
    psp->alpha_RR_2_20=TT_RR/(TT_II+TT_RI+TT_RR);

    l1=21;
    l2=200;

    class_call(spectra_bandpower(psp,l1,l2,&TT_II,&TT_RI,&TT_RR),
               psp->error_message,
               psp->error_message);

    class_test(TT_II+TT_RI+TT_RR==0.,
               psp->error_message,
               "should never happen");

    psp->alpha_II_21_200=TT_II/(TT_II+TT_RI+TT_RR);
    psp->alpha_RI_21_200=TT_RI/(TT_II+TT_RI+TT_RR);
    psp->alpha_RR_21_200=TT_RR/(TT_II+TT_RI+TT_RR);

    l1=201;
    l2=2500;

    class_call(spectra_bandpower(psp,l1,l2,&TT_II,&TT_RI,&TT_RR),
               psp->error_message,
               psp->error_message);

    class_test(TT_II+TT_RI+TT_RR==0.,
               psp->error_message,
               "should never happen");

    psp->alpha_II_201_2500=TT_II/(TT_II+TT_RI+TT_RR);
    psp->alpha_RI_201_2500=TT_RI/(TT_II+TT_RI+TT_RR);
    psp->alpha_RR_201_2500=TT_RR/(TT_II+TT_RI+TT_RR);

    l1=2;
    l2=2500;

    class_call(spectra_bandpower(psp,l1,l2,&TT_II,&TT_RI,&TT_RR),
               psp->error_message,
               psp->error_message);

    class_test(TT_II+TT_RI+TT_RR==0.,
               psp->error_message,
               "should never happen");

    psp->alpha_II_2_2500=TT_II/(TT_II+TT_RI+TT_RR);
    psp->alpha_RI_2_2500=TT_RI/(TT_II+TT_RI+TT_RR);
    psp->alpha_RR_2_2500=TT_RR/(TT_II+TT_RI+TT_RR);

    if (ppt->has_cdi==_TRUE_) {

      psp->alpha_kp=ppm->f_cdi*ppm->f_cdi
        /(1.+ppm->f_cdi*ppm->f_cdi);

      psp->alpha_k1=ppm->f_cdi*ppm->f_cdi*exp((ppm->n_cdi-ppm->n_s)*log(0.002/ppm->k_pivot))
        /(1.+ppm->f_cdi*ppm->f_cdi*exp((ppm->n_cdi-ppm->n_s)*log(0.002/ppm->k_pivot)));

      psp->alpha_k2=ppm->f_cdi*ppm->f_cdi*exp((ppm->n_cdi-ppm->n_s)*log(0.1/ppm->k_pivot))
        /(1.+ppm->f_cdi*ppm->f_cdi*exp((ppm->n_cdi-ppm->n_s)*log(0.1/ppm->k_pivot)));
    }

    if (ppt->has_nid==_TRUE_) {

      psp->alpha_kp=ppm->f_nid*ppm->f_nid
        /(1.+ppm->f_nid*ppm->f_nid);

      psp->alpha_k1=ppm->f_nid*ppm->f_nid*exp((ppm->n_nid-ppm->n_s)*log(0.002/ppm->k_pivot))
        /(1.+ppm->f_nid*ppm->f_nid*exp((ppm->n_nid-ppm->n_s)*log(0.002/ppm->k_pivot)));

      psp->alpha_k2=ppm->f_nid*ppm->f_nid*exp((ppm->n_nid-ppm->n_s)*log(0.1/ppm->k_pivot))
        /(1.+ppm->f_nid*ppm->f_nid*exp((ppm->n_nid-ppm->n_s)*log(0.1/ppm->k_pivot)));
    }

    if (ppt->has_niv==_TRUE_) {

      psp->alpha_kp=ppm->f_niv*ppm->f_niv
        /(1.+ppm->f_niv*ppm->f_niv);

      psp->alpha_k1=ppm->f_niv*ppm->f_niv*exp((ppm->n_niv-ppm->n_s)*log(0.002/ppm->k_pivot))
        /(1.+ppm->f_niv*ppm->f_niv*exp((ppm->n_niv-ppm->n_s)*log(0.002/ppm->k_pivot)));

      psp->alpha_k2=ppm->f_niv*ppm->f_niv*exp((ppm->n_niv-ppm->n_s)*log(0.1/ppm->k_pivot))
        /(1.+ppm->f_niv*ppm->f_niv*exp((ppm->n_niv-ppm->n_s)*log(0.1/ppm->k_pivot)));
    }
  }
  
  if (pnl->method == nl_halomodel){   
    /*printf("%d, %d, %d, %d",pnl->method, psp->hm_convention, psp->hm_convention_cm, psp->hm_highest_order);*/
    class_call(spectra_sigma_table(pba,ppm,psp),psp->error_message,psp->error_message);    
    class_call(halomodel(pba,psp,ppm,ppt),psp->error_message,psp->error_message);
    double pk, pk_ic, bk, bk_lin, q123, q123_lin;
    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,1,0,&pk,&pk_ic),psp->error_message,psp->error_message);
    if (psp->hm_highest_order >2){
    class_call(spectra_bk_lin_at_k_and_z(pba,ppm,psp,1,0,&bk_lin),psp->error_message,psp->error_message);
    class_call(spectra_bk_hm_at_k_and_z(pba,ppm,psp,1,0,&bk),psp->error_message,psp->error_message);
    class_call(spectra_Q123lin_at_k_and_z(pba,ppm,psp,1,0,&q123_lin),psp->error_message,psp->error_message);
    class_call(spectra_Q123_at_k_and_z(pba,ppm,psp,1,0,&q123),psp->error_message,psp->error_message);
    printf("PK(1,0)=%f\nBK_lin(1,0)=%f\nBK(1,0)=%f\nQ_lin(1,0)=%f\nQ(1,0)=%f\n",pk,bk_lin,bk,q123_lin,q123);
    }
    /*int i;
  double * bkarray;
  double * output;
  class_alloc(bkarray,sizeof(double)*psp->hm_k_slices,psp->error_message);
  class_call(spectra_bk_hm_at_z(pba,psp,1,&output),psp->error_message,psp->error_message)
  class_call(spectra_bk_hm_at_k_and_z(pba,ppm,psp,1,0,&bk),psp->error_message,psp->error_message);
  printf("bk1z0=%f\n",bk);
  for (i=0;i<psp->hm_k_slices;i++){
    bkarray[i] = output[i];
    printf("bz0=%f\n",bkarray[i]);
  }*/
  }
  
   return _SUCCESS_;
}



/**
 * This routine computes a table for the nonlinear power, bi and trispectra via the Halo model approach
 * 
 * In the input file one can choose 
 *  - the maximal order for the spectra & the geometric shape of the correlation fct
 *  - The convention for several ingredients of the spectrum )
 *     - multiplicity function
 *     - Bias function
 *     - Concentration mass relation
 * 
 * As the linear Pk is just computed up to k_max also the integral for the variance will be cut off
 * resulting in an unphysical behavior at small smoothing scales. We 'fix' this by introducing
 * a lower mass cutoff and using the renormalization procedure of Schmidt in order to tame it.
 *  
 *
 * @param psp Input: pointer to spectra structure (which fields must be freed)
 * @param pba Input: pointer to background structure
 * @param ppm Input: pointer to primordial structure
 * @param ppt Input: pointer to perturbation structure
 * @return the error status
 * 
 * */

int halomodel(struct background * pba, 
              struct spectra * psp,
	      struct primordial * ppm,
	      struct perturbs * ppt){
  
  if (psp->spectra_verbose > 0)
    printf("Computing non-linear matter power spectrum with Halomodel.\n");
  /* Check if input variables may cause problems in computation */
  /*class_test((ppt->k_max_for_pk < psp->hm_k_min) || (ppt->k_max_for_pk > ppt->k_max_for_pk ),
             psp->error_message,
             "ppt->k_max_for_pk =%e not in valid range [%e:%e]\n",ppt->k_max_for_pk,psp->hm_k_min,ppt->k_max_for_pk);*/
 
  class_test(psp->z_max_pk < 0. || (psp->z_max_pk > psp->z_max_pk ),
             psp->error_message,
             "psp->z_max_pk=%e not in valid range [%e:%e]\n",psp->z_max_pk,0.,psp->z_max_pk);             
  
  /* Make sure that highest order chosen is working in code. */
  int * highest_order_implementedarray;
  int number_conventions_n = 3; /*ST,T10,B15*/
  class_alloc(highest_order_implementedarray,sizeof(int)*number_conventions_n,psp->error_message);
  highest_order_implementedarray[0] = 4; /* ST */
  highest_order_implementedarray[1] = 2; /* T10 */
  highest_order_implementedarray[2] = 2; /* B15 */
  
  int highest_order_implemented = highest_order_implementedarray[psp->hm_convention - 1];
  
  free(highest_order_implementedarray);
  
  class_test((psp->hm_highest_order < 2),
             psp->error_message,
             "You should get at least want the power spectrum as output. Please change hm_highest_order (currently set to %d) to an int bigger than 1",psp->hm_highest_order);
  if (psp->hm_highest_order > highest_order_implemented){
    if (psp->hm_highest_order > highest_order_implemented)
    printf("For the chosen mass function this order (%d) has not been implemented yet.\n", psp->hm_highest_order);
    psp->hm_highest_order = highest_order_implemented;
  }
  if (psp->spectra_verbose > 0){
    if (psp->hm_highest_order == 2)
      printf("Calculations will be done up to the Power Spectrum.\n");
    if (psp->hm_highest_order == 3)
      printf("Calculations will be done up to the Bispectrum.\n");
    if (psp->hm_highest_order == 4)
      printf("Calculations will be done up to the Trispectrum.\n");
  }
  
  /* The parameters ofe ach convention for the mass function are fitted wrt a fixed halo definition has its own halo cut_off defined --> Check if choice of c(M) is ok.  */
  if (psp->hm_convention == 3)
    psp->hm_convention_cm == 3;
   
  if (psp->hm_convention_cm == 3)
    psp->hm_convention_overdensity = 200;
  
  if (psp->hm_convention != 3)
    psp->hm_convention_overdensity = 178;
  
  double per_z;
  if (psp->z_max_pk == 0)
    per_z = 1.;
  else
    per_z = ((double)(psp->hm_z_slices - 1) / (double)(psp->z_max_pk));
  
  
  
  /*printf("amount of m-slices: %d\n",psp->hm_m_slices);
  printf("amount of z-slices: %d\n",psp->hm_z_slices);
  printf("amount of k-slices: %d\n",psp->hm_k_slices);*/
  
  /* Iteration indices (i for masses, j for redshifts, l for wavenumbers) */
  int i,j,l;
  
 
/* Arrays of quantities that will be needed lateron */

  double * aarray;
  double * Deltaarray;
  double * rhobar_bgarray;
  double * Omega_marray;
  double * rhobar_marray;
  double * rhobar_cut_halo;/* Where do you cut off halo? --> Depends on chosen conventions! */
  double * D_plus_fitarray;
  double * D_plus_classarray;
  double * log10marray;
  double * m_to_Rarray;
  double * nuarray;
  double * nutolog10marray;
  double * nu1tolog10marray;
  double * carray;
  double * r_sarray;
  double * rhoarray;
  double * Siarray;
  double * Ciarray;
  double * ukmarray;
  
  class_alloc(aarray,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(Deltaarray,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(rhobar_bgarray,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(Omega_marray,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(rhobar_marray,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(rhobar_cut_halo,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(D_plus_fitarray,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(D_plus_classarray,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(log10marray,psp->hm_m_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(psp->log10karray,psp->hm_k_slices*sizeof(double),psp->error_message);
  class_alloc(m_to_Rarray,psp->hm_z_slices*psp->hm_m_slices*sizeof(double),psp->error_message);
  class_alloc(nuarray,psp->hm_z_slices*psp->hm_m_slices*sizeof(double),psp->error_message);
  class_alloc(nutolog10marray,psp->hm_z_slices*psp->hm_m_slices*sizeof(double),psp->error_message);
  class_alloc(nu1tolog10marray,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(carray,psp->hm_z_slices*psp->hm_m_slices*sizeof(double),psp->error_message);
  class_alloc(r_sarray,psp->hm_z_slices*psp->hm_m_slices*sizeof(double),psp->error_message);
  class_alloc(rhoarray,psp->hm_z_slices*psp->hm_m_slices*sizeof(double),psp->error_message);
  class_alloc(Siarray,psp->hm_z_slices*psp->hm_m_slices*sizeof(double),psp->error_message);
  class_alloc(Ciarray,psp->hm_z_slices*psp->hm_m_slices*sizeof(double),psp->error_message); 
  class_alloc(ukmarray,psp->hm_k_slices*psp->hm_z_slices*psp->hm_m_slices*psp->hm_highest_order*sizeof(double),psp->error_message);
  /* Todays matter energy density */
  double Omega0_m = pba->Omega0_cdm + pba->Omega0_b;
  
  
  /*int number_components = 2;
  double * hm_componentsz0;
  class_alloc(hm_componentsz0,number_components*sizeof(double),psp->error_message);
  double Omega0_stars = 0;
  double Omega0_gas = 0;*/
  
  /* Gives scale factor as a fct of redshift */
  for (j=0;j<psp->hm_z_slices;j++){
    aarray[j] = 1/(1 + (double)(j)/per_z);
    /*printf( "a (%f)= %f \n", (double)(j)/per_z,aarray[j]); */
  } 
  
  class_alloc(psp->hm_zarray,sizeof(double)*psp->hm_z_slices,psp->error_message);
  for(j=0;j<psp->hm_z_slices;j++){
    psp->hm_zarray[j] = (double)(j)/per_z;
    /*printf("z[%le] = %f",(double)(j),psp->hm_zarray[j]);*/
  }
  
  
  /* Gives Delta as defined in Massara (-> virial overdensity)*/  
  for (j=0;j<psp->hm_z_slices;j++){
    if (psp->hm_convention_overdensity == 178){
      double x = Omega0_m * pow(aarray[j],-3.) / 
              (Omega0_m * pow(aarray[j],-3.) + 
              (1. - Omega0_m)) - 1;
      Deltaarray[j] = (18. * pow(_PI_,2.) + 82.*x - 39.*pow(x,2.)) / (1. + x);
    }
    if (psp->hm_convention_overdensity == 200)
      Deltaarray[j] = 200.;
     /* printf( "rhobar_vir (%f)= %f \n", j/per_z,rhobar_cut_halo[j]); */
    /* printf( "Delta (%f)= %f \n", (double)(j)/per_z,Deltaarray[j]); */
  }    

   /* Gives matter energy density at z
      --> as background module works in tau at first get the tau values for desired points in chosen redshift slicing. 
          Afterwards insert them in H0 */
  int last_index_back;  
  double * pvecback_sp_long; 
  double * Harray;
  
  class_alloc(pvecback_sp_long,pba->bg_size*sizeof(double),psp->error_message);
  class_alloc(psp->ln_tau_hm,pba->bg_size*sizeof(double),psp->error_message);
  class_alloc(Harray,pba->bg_size*sizeof(double),psp->error_message);
  
  for (j=0;j<psp->hm_z_slices;j++){
    double tau;
    class_call(background_tau_of_z(pba,(double)(j)/per_z,&tau), pba->error_message, psp->error_message);
    psp->ln_tau_hm[j] = log(tau);
    /*printf("tau(%f) = %f \n",(double)(j)/per_z,psp->ln_tau_hm[j]); */
  }
 
  for (j=0 ; j < psp->hm_z_slices; j++) {     
    double tau;
    double Omega_m = Omega0_m * pow(aarray[j],-3) * pow(Harray[0]/Harray[j],2);
    class_call(background_tau_of_z(pba,(double)(j)/per_z,&tau), pba->error_message, psp->error_message);
    class_call(background_at_tau(pba,
                                 tau,
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback_sp_long),
               pba->error_message,
               psp->error_message);  
    Harray[j] = pvecback_sp_long[pba->index_bg_H];
    rhobar_marray[j]=Omega0_m * pow(1+(double)(j)/per_z,3) * pow(pvecback_sp_long[pba->index_bg_H],2)  / (8.0 / 3.0 *_PI_*_G_) * pow(_c_,2) / _M_SUN_ * _Mpc_over_m_ ;    
    /* printf("rhobar_m(%f) = %f\n",(double)(j)/per_z,rhobar_marray[j]); */
    /* printf("H(%f) = %f\n",(double)(j)/per_z,Harray[j]); */
  } 
  
  
  /* Gives background energy density inside cutoff radius */
  for (j=0;j<psp->hm_z_slices;j++){
    rhobar_cut_halo[j] = rhobar_marray[j] * (Deltaarray[j]-1);
  }
  
  
  /* LCDM growth fitting function from background module*/
  for (j=0;j<psp->hm_z_slices;j++){
    double tau, Norm;
    
    class_call(background_tau_of_z(pba,(double)(j)/per_z,&tau), pba->error_message, psp->error_message);
    class_call(background_at_tau(pba,
                                 tau,
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback_sp_long),
               pba->error_message,
               psp->error_message); 
    
    if (j == 0){
      Norm = pvecback_sp_long[pba->index_bg_D];
    }
          
    D_plus_classarray[j] =  pvecback_sp_long[pba->index_bg_D] / Norm;    
    /*printf( "D_plus_class (%f)= %f \n", (double)((j)/per_z),D_plus_classarray[j]);*/
  }
  
  /* Allocate maximal and minimal values (could be made z-dependent..) 
     Check if lowest mass cares about high-k-cutoff in the sigma integral */
  double * m_minarray;
  double * m_maxarray;
  class_alloc(m_minarray,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(m_maxarray,psp->hm_z_slices*sizeof(double),psp->error_message);
  for (j=0;j<psp->hm_z_slices;j++){
    m_minarray[j] = psp->hm_m_min /* In input file: 10**10 */ /** (1 + 10* psp->hm_zarray[j])*/;
    m_maxarray[j] = psp->hm_m_max /* In input file: 10**18 */ /** (1 + 10* psp->hm_zarray[j]) */;
    /*if (m_minarray[j] < pow(10,8.5)*pow(ppt->k_max_for_pk/30.,-3)){
      printf("Note: Low mass cutoff set from %f to %f in order to have a good table for variance. If you want a lower mass cutoff choose a higher P_k_max_h/Mpc.\n",m_minarray[j],pow(10,8.5)*pow(ppt->k_max_for_pk/30.,-3));
      m_minarray[j] = pow(10,8.5)*pow(ppt->k_max_for_pk/30.,-3);/* /* Set value to lowest possible mass that doesnt conclict with cutoff.
                                                             Just obtained from one cosmology but shouldnt change too much... */
    /*}*/
  }
   if (psp->spectra_verbose > 0)
      printf("Low mass cutoff at k=%f.\n",m_minarray[0]);
  
  
 
  
 /* Gives array of m values integrated over */
  for (i=0;i<psp->hm_m_slices;i++){
    for (j=0;j<psp->hm_z_slices;j++){
    log10marray[i*psp->hm_z_slices+j] = log10(m_minarray[j]) + log10(m_maxarray[j]/m_minarray[j]) * 
              (double)((double)(i)/(double)(psp->hm_m_slices - 1));
    /*printf( "log10marray (%f)= %f \n", (double)(i),log10marray[i*psp->hm_z_slices+j]); */
    }
  }
  
  /* Bi or Trispectrum depend on 3/4 k vectors, but everything can be reexpressed in terms of their norms (Cooray sheth 2002, eq (35)).
       --> To reduce computation time reexpress several combinations of vectors in terms of the longest vector */
    /* k_i is vector between point i and i+1 */
    int number_components = 3;
    int number_vectors;
    int index_order_max = 0; /* Just needed for Bi-/Trispectra */

    if (psp->hm_highest_order == 3)
      number_vectors = 3;
    if (psp->hm_highest_order == 4)
      number_vectors = 6;
    double * veccomponentarray;
    double * lengtharray;
    double * rel_lengtharray;
    double * scalar_prod;
    double * muarray;
    double * Fs_2array;
    
    /* If we choose pk as highest order just one absolute value of k interesting 
       --> Both k's have relative length 1. */
    if (psp->hm_highest_order == 2){
      number_vectors = 1;
      class_alloc(rel_lengtharray,number_vectors*sizeof(double),psp->error_message); 
      rel_lengtharray[0] = 1;
      /*rel_lengtharray[1] = 1;*/
    }
    
    
    if (psp->hm_highest_order > 2){    
      int index_order;
      int index_component, index_component_2;
      int vec_1, vec_2;   
      class_alloc(muarray,psp->hm_z_slices*sizeof(double),psp->error_message);
      class_alloc(Fs_2array,pow(number_vectors,2)*psp->hm_z_slices*sizeof(double),psp->error_message);
      class_alloc(veccomponentarray,number_vectors * number_components *sizeof(double),psp->error_message); 
      class_alloc(lengtharray,number_vectors*sizeof(double),psp->error_message); 

      class_alloc(rel_lengtharray,number_vectors*sizeof(double),psp->error_message); 
      class_alloc(scalar_prod,pow(number_vectors,2)*sizeof(double),psp->error_message);
  

      /* Allocate vectors and compute sums of vectors needed later on (See sheth 2002 p.41)
         --> For Bispectrum we just need vectors and consecutive sums (3+3 = 6 terms)
             --> 1st row: vectors (3 entries)
                 2nd row: conses sums (3 entries)
         --> For Trispectrum we need vectors, consecutive sums and all sums of two vectors
             (in total 4 + 4 + 4*4 = 24 terms) - Note that 8 components (0 + 1 + 7) are not needed but hard to implement this exclusion in loops... 
             --> 1st row: vectors (4 entries)
                 2nd row: consecutive sums (4 entries)
                 3rd+ith row: vec_i + {vec} (4 entries each) */
      veccomponentarray[0*number_components+0]=psp->hm_k_11;
      veccomponentarray[0*number_components+1]=psp->hm_k_12;
      veccomponentarray[0*number_components+2]=psp->hm_k_13;
      veccomponentarray[1*number_components+0]=psp->hm_k_21;
      veccomponentarray[1*number_components+1]=psp->hm_k_22;
      veccomponentarray[1*number_components+2]=psp->hm_k_23;
      veccomponentarray[2*number_components+0]=psp->hm_k_31;
      veccomponentarray[2*number_components+1]=psp->hm_k_32;
      veccomponentarray[2*number_components+2]=psp->hm_k_33;
      if (psp->hm_highest_order >3){
        veccomponentarray[3*number_components+0]=psp->hm_k_41;
        veccomponentarray[3*number_components+1]=psp->hm_k_42;
        veccomponentarray[3*number_components+2]=psp->hm_k_43;
	veccomponentarray[4*number_components+0]=psp->hm_k_11+psp->hm_k_21;
        veccomponentarray[4*number_components+1]=psp->hm_k_12+psp->hm_k_22;
        veccomponentarray[4*number_components+2]=psp->hm_k_13+psp->hm_k_23;
	veccomponentarray[5*number_components+0]=psp->hm_k_21+psp->hm_k_31;
        veccomponentarray[5*number_components+1]=psp->hm_k_22+psp->hm_k_32;
        veccomponentarray[5*number_components+2]=psp->hm_k_23+psp->hm_k_33;
	
      }
          
      
     /* Compute length of each vector, find maximal length of sides of n-eck and compute relative length to first vector (correct squeezing in bispectrum!!)*/
     double max = 0; 
     for (index_order = 0; index_order < number_vectors; index_order++){
        double  sum = 0;
        for (index_component = 0; index_component < number_components; index_component++)
          sum += pow(veccomponentarray[index_order*number_components+index_component],2);
        lengtharray[index_order] = pow(sum,.5);	
	if (lengtharray[index_order] > max){
	  max = lengtharray[index_order];
          index_order_max = index_order;
	}
        printf("length[%d] = %f\n",index_order,lengtharray[index_order]);
      }
      
      /* one could divide to max but then base is changed */
      for (index_order = 0; index_order < number_vectors; index_order++){
	rel_lengtharray[index_order] = lengtharray[index_order]/lengtharray[index_order_max];
        printf("rel_length[%d] = %f\n",index_order,rel_lengtharray[index_order]);
      }
      
      /* Compute all possible scalar products between vecs and store them in symmetric matrix */
      for (vec_1 = 0; vec_1 < number_vectors; vec_1++){
	for (vec_2 = 0; vec_2 < number_vectors; vec_2++){
	  scalar_prod[vec_1*number_vectors+vec_2] = 0;
	  for (index_component = 0;index_component<number_components;index_component++){
	    scalar_prod[vec_1*number_vectors+vec_2] += veccomponentarray[vec_1*number_components+index_component]*veccomponentarray[vec_2*number_components+index_component];	    
	  }
	  /*printf("vec_%d * vec_%d = %f\n",vec_1,vec_2,scalar_prod[vec_1*number_vectors+vec_2]);*/
	}
      }
      
      /* Allocate F^s_2 - kernel appearing in Bi-/Trispectrum calculation for all vec combinations (does not depend on magnitude, just orientation) */
      int slot_1, slot_2;
      for(j=0;j<psp->hm_z_slices;j++){
        muarray[j] = 3./7. * pow(Omega0_m*pow(D_plus_classarray[j],-3),-2./63.); 
	/*printf("muarray = %f\n",muarray[j]);*/
        for (slot_1=0;slot_1 < number_vectors;slot_1++){
	  for (slot_2=0;slot_2 < number_vectors;slot_2++){
	    if (lengtharray[slot_1] == 0|| lengtharray[slot_2] == 0)
	      Fs_2array[j*number_vectors*number_vectors+slot_1*number_vectors+slot_2] = 0.;
	    else {
              Fs_2array[j*number_vectors*number_vectors+slot_1*number_vectors+slot_2] 
              = 0.5*((1+muarray[j])+scalar_prod[slot_1*number_vectors+slot_2]/(lengtharray[slot_1]*lengtharray[slot_2])*(lengtharray[slot_1]/lengtharray[slot_2]+lengtharray[slot_2]/lengtharray[slot_1])
	            +(1-muarray[j])*pow(scalar_prod[slot_1*number_vectors+slot_2],2)/(pow(lengtharray[slot_1]*lengtharray[slot_2],2)));
	    /*printf("Fs_2 = %f\n",Fs_2array[j*psp->hm_highest_order*psp->hm_highest_order+slot_1*psp->hm_highest_order+slot_2]);*/
	    }
	  }
	}
      } 
  }
  
  double k_max_for_bk = ppt->k_max_for_pk;
  /* Gives array of k values integrated over */
  /* 1./rel_lengtharray[index_order_max] correction arises due to the fact that base might not be longest line and so we dont run into trouble calling pk. */
  for (l=0;l<psp->hm_k_slices;l++){
    if (psp->hm_highest_order < 3){
      psp->log10karray[l] = log10(psp->hm_k_min) + log10(k_max_for_bk/psp->hm_k_min) * 
              (double)((double)(l)/(double)(psp->hm_k_slices - 1));
    }
    else{
      if (rel_lengtharray[index_order_max] * k_max_for_bk > ppt->k_max_for_pk){
	printf("Warning: For this orientation we can just compute Bispectrum up to %f",1./rel_lengtharray[index_order_max] * ppt->k_max_for_pk);
	k_max_for_bk = 1./rel_lengtharray[index_order_max] * ppt->k_max_for_pk;
      }
      psp->log10karray[l] = log10(psp->hm_k_min) + log10(k_max_for_bk/psp->hm_k_min) * 
              (double)((double)(l)/(double)(psp->hm_k_slices - 1));
    }
	      
    /*printf( "log10karray (%f)= %f \n", (double)(l),psp->log10karray[l]); 
    printf ("k_max = %f\n",ppt->k_max_for_pk);*/
  }
  
  
  
  /* Gives array of R(m,z) in logspaced mass interval */  
  double c_1 =  log10(3. / 4. / _PI_);
  for (i=0;i<psp->hm_m_slices;i++){
    for (j=0;j<psp->hm_z_slices;j++){
      m_to_Rarray[i*psp->hm_z_slices+j] = pow(10,1./3. * 
                       (c_1 - log10(rhobar_marray[j]) + log10marray[i*psp->hm_z_slices+j]));
      /*printf( "R(10** %le, %f)=%f\n",log10marray[(int)(floor(i/psp->hm_z_slices))],j/per_z,m_to_Rarray[i*psp->hm_z_slices+j]);*/
      /* alternative way of calculating R_min & R_max (boundaries for splining of sigma) if weird cosmology
      if (log10(m_to_Rarray[i*psp->hm_z_slices+j]) > log10R_max)
	log10R_max = log10(m_to_Rarray[i*psp->hm_z_slices+j]);
      if (log10(m_to_Rarray[i*psp->hm_z_slices+j]) <= log10R_min)
	log10R_min = log10(m_to_Rarray[i*psp->hm_z_slices+j]); */    
    }
  }
  
  /* Calls sigma psp->hm_m_slices * psp->hm_z_slices times --> slow but correct! */
  double sigma1;
  for (i = 0;i<psp->hm_m_slices;i++){
    for (j = 0;j<psp->hm_z_slices;j++){
      class_call(spectra_sigma_at_R_and_z(pba,psp,
                              m_to_Rarray[i*psp->hm_z_slices+j],((double)(j))/per_z,
                              &sigma1),
                psp->error_message,
                psp->error_message);
    nuarray[i*psp->hm_z_slices+j]= (_DELTA_C_*(1+(0.353* pow(-1,4) + 1.044*pow(-1,3)+1.128*pow(-1,2)+0.555*(-1)+0.131  ) * log10(Omega0_m*D_plus_classarray[j]))) / sigma1 /*/ D_plus_classarray[j] / pow((1+psp->hm_zarray[j]),.5)*/;  
    /*printf( "nu1(1e%le,%f)=%f \n",log10marray[i*psp->hm_z_slices+j],(double)(j)/per_z,nuarray[i*psp->hm_z_slices+j]);*/
    } 
  }
  
 /* Gives log10m as a function of nu 
     Needs to be splined for 1 Halo term */
  for (i = 0;i<psp->hm_m_slices;i++){
    for (j=0;j<psp->hm_z_slices;j++){
      nutolog10marray[i*psp->hm_z_slices+j] = log10marray[i*psp->hm_z_slices+j];
    }
  }
  

  /* Find approximately log10(m(nu=1)) -> for each value of z gives such this value*/
  /* Assume linear function --> nu1tolog10m works good (better that 1%) */
  /* Underestimates vale (? -> as nu(m) is (mostly) convex function) */
  /* Check if sensitve under percent difference eg just take the one for z = 0 or spline  */
  
  /* First entry to check weather nu = 1 
   * gets lower for higher redshift 
   *-> this is safe for m_min ~ 5 & m_max ~17 */
  int nu1tom = 0;
  for (j = 0;j<psp->hm_z_slices;j++){
    if (nuarray[0*psp->hm_z_slices+j] > 1.0){
	nu1tolog10marray[j] = nutolog10marray[0*psp->hm_z_slices+j];
	printf("Nonlinear mass for z=%f < %f. First crossing chosen to be %f where nu(%f,%f) = %f.\n",psp->hm_zarray[j],m_minarray[j],m_minarray[j],m_minarray[j],psp->hm_zarray[j],nuarray[j]);
      }
    for (i=0;i<psp->hm_m_slices;i++){
      if (nuarray[i*psp->hm_z_slices+j] >= 1.0){      
        nu1tolog10marray[j] = nutolog10marray[i*psp->hm_z_slices+j] - (nutolog10marray[i*psp->hm_z_slices+j]-nutolog10marray[(i-1)*psp->hm_z_slices+j])*
        (nuarray[i*psp->hm_z_slices+j]-1.0)/(1.0-nuarray[(i-1)*psp->hm_z_slices+j] + nuarray[i*psp->hm_z_slices+j]-1.0);
        /*printf("log10m(nu=1,z= %f) ~ %le \n",(double)(j)/per_z, nu1tolog10marray[j]); */
        nu1tom = i;
	break;
      }
    }
  }
  
 
 /*int hm_cm = 1;*/
  /* Gives concentration mass relation of Ludlow et al (2016) (Fitted to Planck) (https://arxiv.org/pdf/1601.02624v2.pdf)  */
  /* Too low but right asymptotics for higher z! */
  if (psp->hm_convention_cm == 2){
    if (psp->spectra_verbose > 0)
      printf("Using concentration mass relation fitted to Planck data by Ludlow et al. 2016.\n");
    double * c0;
    double * beta;
    double * gamma_1;
    double * gamma_2;
    double * nu0;
    class_alloc(c0,psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(beta,psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(gamma_1,psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(gamma_2,psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(nu0,psp->hm_z_slices*sizeof(double),psp->error_message);
    
    
    for (j = 0;j<psp->hm_z_slices;j++){
      c0[j] = 3.395*pow(1.+ psp->hm_zarray[j],-.215);
      beta[j] = .307*pow(1.+ psp->hm_zarray[j],.540);
      gamma_1[j] = .628*pow(1.+ psp->hm_zarray[j],-.047);
      gamma_2[j] = .317*pow(1.+ psp->hm_zarray[j],-.893);
      nu0[j] = (4.135-.564*pow(1.+ psp->hm_zarray[j],1)-.210*pow(1.+ psp->hm_zarray[j],2)
	       +.0557*pow(1.+ psp->hm_zarray[j],3)-.00348*pow(1.+ psp->hm_zarray[j],4))/D_plus_classarray[j];
      for (i=0;i<psp->hm_m_slices;i++){
        carray[i*psp->hm_z_slices+j] = 
          c0[j]*pow((nuarray[i*psp->hm_z_slices+j]/nu0[j]),-1.*gamma_1[j]) * pow((
	  1.+pow((nuarray[i*psp->hm_z_slices+j]/nu0[j]),1./beta[j])),-beta[j]*(gamma_2[j]-gamma_1[j]));  
          
        /*printf("c(1e%f, %f) ~ %f \n", log10marray[i*psp->hm_z_slices+j],(double)(j)/per_z,carray[i*psp->hm_z_slices+j]); */
      }
    }
    free(c0);
    free(beta);
    free(gamma_1);
    free(gamma_2);
    free(nu0);
  }
  
  /* Gives concentration mass relatinon of cooray sheth review */
  if (psp->hm_convention_cm == 1){
    if (psp->spectra_verbose > 0)
      printf("Using concentration mass relation of Cooray & Sheth 2002.\n");
    for (j = 0;j<psp->hm_z_slices;j++){
      for (i=0;i<psp->hm_m_slices;i++){	
        carray[i*psp->hm_z_slices+j] = 9*pow(pow(10,nutolog10marray[i*psp->hm_z_slices+j] - nu1tolog10marray[j]),-.13) * aarray[j];
        /*printf("c(1e%f, %f) ~ %f \n", log10marray[i*psp->hm_z_slices+j],(double)(j)/per_z,carray[i*psp->hm_z_slices+j]); */
      }
    }
  }
  
  /* Gives concentration mass relation of Duffy 2008 with unrelaxed parameteres */
  if (psp->hm_convention_cm == 3){
    if (psp->spectra_verbose > 0)
      printf("Using concentration mass relation of Duffy 2008.\n");
    double A = 10.14;
    double B = -.081;
    double C = -1.01;
    double M_piv = 2*pow(10,12);
    for (j = 0;j<psp->hm_z_slices;j++){
      for (i=0;i<psp->hm_m_slices;i++){	
        carray[i*psp->hm_z_slices+j] = A*pow(pow(10,nutolog10marray[i*psp->hm_z_slices+j])/M_piv,B) * pow(1+psp->hm_zarray[j],C);
        /*printf("c(1e%f, %f) ~ %f \n", log10marray[i*psp->hm_z_slices+j],(double)(j)/per_z,carray[i*psp->hm_z_slices+j]); */
      }
    }
  }
  
  /* Gives core radius */
  for (i = 0;i<psp->hm_m_slices;i++){
    for (j=0;j<psp->hm_z_slices;j++){
      r_sarray[i*psp->hm_z_slices+j] = pow((3./(4. * _PI_ * rhobar_cut_halo[j]) * pow(10,log10marray[i*psp->hm_z_slices+j])/ 
      pow(carray[i*psp->hm_z_slices+j],3.)),1./3.);
      /* printf("r_s(1e%f, %f) ~ %f \n", log10marray[i*psp->hm_z_slices+j],(double)(j)/per_z,r_sarray[i*psp->hm_z_slices+j]); */
    }
  }
  
  /* Gives density in core radius of halo */
  for (i = 0;i<psp->hm_m_slices;i++){
    for (j = 0;j<psp->hm_z_slices;j++){
      rhoarray[i*psp->hm_z_slices+j] =  pow(10,log10marray[i*psp->hm_z_slices+j]) / (4. * _PI_ * pow(r_sarray[i*psp->hm_z_slices+j],3) * (
      log(1. + carray[i*psp->hm_z_slices+j]) - carray[i*psp->hm_z_slices+j]/(1 + carray[i*psp->hm_z_slices+j])));
       /*printf("rho_s(1e%f, %f) ~ %f \n", log10marray[i*psp->hm_z_slices+j],(double)(j)/per_z,log10(rhoarray[i*psp->hm_z_slices+j])); */
    }
  }
  
  
  /* Sici from scipy routine. Found on
  /* http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/mcisia_cpp.txt */ 

      double CI,SI,CI1,SI1,X,Y;

      void cisi(double X, double *CI, double *SI) {

        double EL, EPS, P2, X2, XA, XA0, XA1, XCS, XF, XG, XG1, XG2, XR, XS, XSS;
        double BJ[101];
	int K, M;
        P2=1.570796326794897;
        EL=.5772156649015329;
        EPS=1.0e-15;
        X2=X*X;
        if (X == 0.0) { 
           *CI=-1.0e+100;
           *SI=0.0;
        }
        else if (X <= 16.0) {
           XR=-0.25*X2;
           *CI=EL+log(X)+XR;
           for (K=2; K<41; K++) {
              XR=-0.5*XR*(K-1)/(K*K*(2*K-1))*X2;
              *CI=*CI+XR;
              if (fabs(XR) < fabs(*CI)*EPS) goto e15;
	   }
  e15:       XR=X;
             *SI=X;
             for (K=1; K<41; K++) { 
                XR=-0.5*XR*(2*K-1)/K/(4*K*K+4*K+1)*X2;
                *SI=*SI+XR;
                if (fabs(XR) < fabs(*SI)*EPS) return;
	     }
          }
        else if (X <= 32.0) {
           M=(int)(47.2+.82*X);
           XA1=0.0;
           XA0=1.0e-100;
	   for (K=M; K>0; K--) {
              XA=4.0*K*XA0/X-XA1;
              BJ[K]=XA;
              XA1=XA0;
              XA0=XA;
           }
           XS=BJ[1];
           for (K=3; K<=M; K+=2)  XS=XS+2.0*BJ[K];
           BJ[1]=BJ[1]/XS;
           for (K=2; K<=M; K++)   BJ[K]=BJ[K]/XS;
           XR=1.0;
           XG1=BJ[1];
	   for (K=2; K<=M; K++) {
              XR=0.25*XR*pow(2.0*K-3.0,2.0)/((K-1.0)*pow(2.0*K-1.0,2.0))*X;
              XG1=XG1+BJ[K]*XR;
           }
           XR=1.0;
           XG2=BJ[1];
           for (K=2; K<=M; K++) { 
              XR=0.25*XR*pow(2.0*K-5.0,2.0)/((K-1.0)*pow(2.0*K-3.0,2.0))*X;
              XG2=XG2+BJ[K]*XR;
	   }
           XCS=cos(X/2.0);
           XSS=sin(X/2.0);
           *CI=EL+log(X)-X*XSS*XG1+2*XCS*XG2-2*XCS*XCS;
           *SI=X*XCS*XG1+2*XSS*XG2-sin(X);
        }
        else {
           XR=1.0;
           XF=1.0;
	   for (K=1; K<10; K++) {
              XR=-2.0*XR*K*(2*K-1)/X2;
              XF=XF+XR;
           }
           XR=1.0/X;
           XG=XR;
           for (K=1; K<9; K++) {
              XR=-2.0*XR*(2*K+1)*K/X2;
              XG=XG+XR;
           }
           *CI=XF*sin(X)/X-XG*cos(X)/X;
           *SI=P2-XF*cos(X)/X-XG*sin(X)/X;
        }
      }
  
  
  /* Gives Fourier Transform of NFW profile for each vector needed lateron */
  /* Here we need to recall that we work relative to the longest vec */
  int order;
  for (i=0;i<psp->hm_m_slices;i++){
    for (j=0;j<psp->hm_z_slices;j++){
      for (l=0;l<psp->hm_k_slices;l++){
	for (order = 0; order < number_vectors; order++){	
       	  double x,y,ukm_k, ci, si, ci1, si1;
	  ukm_k = pow(10,psp->log10karray[l])*rel_lengtharray[order];
	  y = ukm_k * r_sarray[i*psp->hm_z_slices+j];
          x = (1+carray[i*psp->hm_z_slices+j]) * y;
          cisi(x, &ci, &si);         
          cisi(y, &ci1, &si1);	    
	  ukmarray[order*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] =
	       4 * _PI_ * rhoarray[i*psp->hm_z_slices + j] * pow(r_sarray[i*psp->hm_z_slices + j],3) / pow(10,log10marray[i*psp->hm_z_slices+j]) * (
	       sin(y) * (si - si1) - sin(carray[i*psp->hm_z_slices+j] * y) / x + cos(y) * (ci - ci1) );			 
	  /*printf("ukm(%f,%f,%f) = %f  \n",log10marray[i*psp->hm_z_slices+j],ukm_k,(double)(j)/per_z,ukmarray[order*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l]);*/ 
        }
      }
    }
  }
  
  
  
  /* Arrays of multiplicity and bias functions up to fourth order that will be allocated and filled after convention has been chosen. */
  double * multiplicityarray;
  double * biasarray;
  double * normmultiplicityarray;
  double * normfb0array;
  double * normfb1array;
  double * normfb2array;
  double * multipl_renormarray; /* This is the constant n^bar_s (eq. A2 from Schmidt 2016) */
  double * b0_renormarray; /* This is the constant b_1(M) (eq. A4 from Schmidt 2016) */
  double * b1_renormarray;
  double * b2_renormarray;
  double * normfb3array;
  double * b3_renormarray;   
  double * normfb4array; 
  double * b4_renormarray;
  
  double * bias_pert_a;
  double * bias_pert_eps;
  double * bias_pert_E;
  double * bias_pert;
  double * M_00;
  double * M_10;
  double * M_20;
  double * M_30;
  double * M_40;
  double * M_01;
  double * M_11;
  double * M_21;
  double * M_31;
  double * M_02;
  double * M_12;
  double * M_22;
  double * M_03;
  double * M_13;
  double * M_04;
    
  double * P_1harray;
  double * P_2harray;
  
  double * B_1harray;
  double * B_2harray;
  double * B_3harray;
  
  double * T_1harray;
  double * T_2harray;
  double * T_3harray;
  double * T_4harray;
  
  /* Use Sheth Tormen multiplicity and bias */
  if (psp->hm_convention == 1){
    if (psp->spectra_verbose > 0)
      printf("Using mass function and bias of Sheth Tormen (universal parameters)\n");
    
    class_alloc(multiplicityarray,psp->hm_z_slices * psp->hm_m_slices*sizeof(double),psp->error_message);
    class_alloc(biasarray,psp->hm_m_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
    
    /* Unnormed bias function of ST */
    /* Only the values covered by particular (nu,z) pairs are taken as arguments */
    for (i = 0;i<psp->hm_m_slices;i++){
      for (j = 0;j<psp->hm_z_slices;j++){
	double a = 0.7689;
	double b = 0.2536;
	biasarray[i*psp->hm_z_slices+j] = 
	  1. + (a * pow(nuarray[i*psp->hm_z_slices+j],2) -1) / (_DELTA_C_)  + 2.*b / (_DELTA_C_)  / (1 + pow(a * pow(nuarray[i*psp->hm_z_slices+j],2),b));
	/* printf("bst_unnormed(%f,%f) = %f  \n",pow(nuarray[i*psp->hm_z_slices+j],2),(double)(j)/per_z,biasarray[i*psp->hm_z_slices + j]); */
      }
    }
    
    
    /* Unnormed Multiplicity function of ST */
    for (i = 0;i<psp->hm_m_slices;i++){
      for (j = 0;j<psp->hm_z_slices;j++){
	double A = 0.3295; /*0.322*/
	double a = 0.7689; /*.707*/
	double b = 0.2536; /*.3*/
	multiplicityarray[i * psp->hm_z_slices + j] = 
	  A * (1. + pow((a * pow(nuarray[i*psp->hm_z_slices+j],2)),-.3) ) * pow(a * pow(nuarray[i*psp->hm_z_slices+j],2) / 2. / _PI_ , .5) * 
	  exp(- a * pow(nuarray[i*psp->hm_z_slices+j],2) / 2.) / pow(nuarray[i*psp->hm_z_slices+j],2);
	/*printf("fst_unnormed(%f,%f) = %f  \n",pow(nuarray[i*psp->hm_z_slices+j],2),(double)(j)/per_z,multiplicityarray[i*psp->hm_z_slices + j]);*/
      }
    }
    
    
    /* Expand the bias of chosen mass function according to Cooray review p. 25f  */
    /* Note different convention: nu_sheth = pow(nu_code,2) */ 
    /* TO CHECK: LOW MASS CUTOFF ENFORCED BY CLASS --> RENORMALIZATION NEEDED PER SE (SEE SCHMIDT 2015 FOR POWER & BISPECTRUM)
    /* CHECK QUICK, DIRTY AND EVENTUALLY INCORRECT SOLUTION FROM PYTHON CODE WHERE bias[i] -> bias[i]/bias[0] for each z slice  */
    class_alloc(bias_pert_a,psp->hm_highest_order*sizeof(double),psp->error_message);
    class_alloc(bias_pert_eps,psp->hm_m_slices*psp->hm_z_slices*psp->hm_highest_order*sizeof(double),psp->error_message);
    class_alloc(bias_pert_E,psp->hm_m_slices*psp->hm_z_slices*psp->hm_highest_order*sizeof(double),psp->error_message);
    class_alloc(bias_pert,psp->hm_m_slices*psp->hm_z_slices*(psp->hm_highest_order+1)*sizeof(double),psp->error_message); /* Zeroth compoent := 1 globally (see sheth review) */
    
    /* Found in appendix of http://arxiv.org/pdf/0906.1314.pdf */
    bias_pert_a[0] = 1;
    bias_pert_a[1] = -17./24.;
    if (psp->hm_highest_order > 2)
      bias_pert_a[2] = 341./567.;
    if (psp->hm_highest_order > 3)
      bias_pert_a[3] = -55805./130977;
    
    /* To Check: Here 4th order  divided by delta_c**4 - as in Appendix of Manera et al 2009 (* delta_c**3) or Sheth review (( / delta_c**3)) */
    for (i = 0;i<psp->hm_m_slices;i++){
      for (j = 0;j<psp->hm_z_slices;j++){
	double q = 0.707;	
	bias_pert_eps[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	 = (q*pow(nuarray[i*psp->hm_z_slices+j],2)-1)/(_DELTA_C_) ;
	bias_pert_eps[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	 = (q*pow(nuarray[i*psp->hm_z_slices+j],2) * (q*pow(nuarray[i*psp->hm_z_slices+j],2)-3)) / ((_DELTA_C_)  *(_DELTA_C_) );
	if (psp->hm_highest_order > 2){
          bias_pert_eps[2*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	   = (q*pow(nuarray[i*psp->hm_z_slices+j],2)*(pow(q*pow(nuarray[i*psp->hm_z_slices+j],2),2)-6*q*pow(nuarray[i*psp->hm_z_slices+j],2)+3)) / (pow((_DELTA_C_) ,3 )); 	
	}
	if (psp->hm_highest_order > 3){
	  bias_pert_eps[3*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	   = pow(q*pow(nuarray[i*psp->hm_z_slices+j],2),2) * (pow(q*pow(nuarray[i*psp->hm_z_slices+j],2),2) - 10*q*pow(nuarray[i*psp->hm_z_slices+j],2)+15) / (pow((_DELTA_C_) ,3 ));      
	}	
      }
    }
    

    for (i = 0;i<psp->hm_m_slices;i++){
      for (j = 0;j<psp->hm_z_slices;j++){
	double p = 0.3;
	double q = 0.707;
	bias_pert_E[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	 = (2*p/(_DELTA_C_) )/(1+pow(q*pow(nuarray[i*psp->hm_z_slices+j],2),p));
	bias_pert_E[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	 = ((1+2*p)/(_DELTA_C_)  + 2*bias_pert_eps[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j])
	  *((2*p/(_DELTA_C_) )/(1+pow(q*pow(nuarray[i*psp->hm_z_slices+j],2),p)));
	if (psp->hm_highest_order > 2){
	  bias_pert_E[2*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	   = ((4*(p*p-1)+6*p*q*pow(nuarray[i*psp->hm_z_slices+j],2)) / pow((_DELTA_C_) ,2) + 3*pow(bias_pert_eps[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j],2))
	    *((2*p/(_DELTA_C_) )/(1+pow(q*pow(nuarray[i*psp->hm_z_slices+j],2),p)));
	}
	if (psp->hm_highest_order > 3){
	  bias_pert_E[3*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	 = (2*q*pow(nuarray[i*psp->hm_z_slices+j],2)/(pow((_DELTA_C_) ,2)) * ((2*pow(q*pow(nuarray[i*psp->hm_z_slices+j],2),2)/(_DELTA_C_) ) - 15*bias_pert_eps[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j])
	  + 2*(1+p)/pow((_DELTA_C_) ,2) * ((4*(p*p-1)+8*(p-1)*q*pow(nuarray[i*psp->hm_z_slices+j],2)+3) / (_DELTA_C_)  + 6*q*pow(nuarray[i*psp->hm_z_slices+j],2)*bias_pert_eps[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]))
	  *((2*p/(_DELTA_C_) )/(1+pow(q*pow(nuarray[i*psp->hm_z_slices+j],2),p)));
	}
      }
    }
     
    
    /* Fill expanded bias vector */
    /* WATCH OUT! This vector has (by construction of the matrix M_ij filled lateron) one entry more than the bias_pert_E & bias_pert_eps arrays */
    for (i = 0;i<psp->hm_m_slices;i++){
      for (j = 0;j<psp->hm_z_slices;j++){
	bias_pert[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	          = 1.;
        bias_pert[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	          = 1. + bias_pert_eps[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] + bias_pert_E[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j];
	bias_pert[2*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	          = 2.*(1+bias_pert_a[1])*(bias_pert_eps[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] + bias_pert_E[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j])
		    + bias_pert_eps[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] + bias_pert_E[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j];
	if (psp->hm_highest_order > 2){
	  bias_pert[3*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	          = 6.* (bias_pert_a[1] + bias_pert_a[2]) * (bias_pert_eps[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] + bias_pert_E[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]) 
                    + 3.*(1.+2.*bias_pert_a[1])*((bias_pert_eps[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] + bias_pert_E[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]))
                    + bias_pert_eps[2*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] + bias_pert_E[2*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j];	
	}
	if (psp->hm_highest_order > 3){
          bias_pert[4*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	          = 24.* (bias_pert_a[2] + bias_pert_a[3]) * (bias_pert_eps[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] + bias_pert_E[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j])
		    + 12.* (pow(bias_pert_a[1],2) + 2*(bias_pert_a[1] + bias_pert_a[1])) * (bias_pert_eps[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] + bias_pert_E[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j])
		    + 4.*(1.+3.*bias_pert_a[1]) * (bias_pert_eps[2*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] + bias_pert_E[2*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j])
		    + (bias_pert_eps[3*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] + bias_pert_E[3*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]);
	}
	/*printf("bias_1[%d,%d] = %f\n",i,j,bias_pert[1*psp->hm_m_slices*psp->hm_z_slices+i*psp->hm_z_slices+j]);
	printf("bias_2[%d,%d] = %f\n",i,j,bias_pert[2*psp->hm_m_slices*psp->hm_z_slices+i*psp->hm_z_slices+j]);
	printf("bias_3[%d,%d] = %f\n",i,j,bias_pert[3*psp->hm_m_slices*psp->hm_z_slices+i*psp->hm_z_slices+j]);
	printf("bias_4[%d,%d] = %f\n",i,j,bias_pert[4*psp->hm_m_slices*psp->hm_z_slices+i*psp->hm_z_slices+j]);*/
      }
    }
    free(bias_pert_a);
    free(bias_pert_E);
    free(bias_pert_eps);
  }
  
  /* Multiplicity and bias in Tinker2010 convention */
  if (psp->hm_convention == 2){
    if (psp->spectra_verbose > 0)
      printf("Using mass function and bias of Tinker 2010.\n");
    
    class_alloc(multiplicityarray,psp->hm_z_slices * psp->hm_m_slices*sizeof(double),psp->error_message);
    class_alloc(biasarray,psp->hm_m_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
    
    /* Unnormed bias function of T10 */
    /* Only the values covered by particular (nu,z) pairs are taken as arguments */
    for (i = 0;i<psp->hm_m_slices;i++){
      for (j = 0;j<psp->hm_z_slices;j++){
        double y = log10(Deltaarray[j]);
        double A = 1.0 + 0.24*y*exp(-pow((4./y),4.));
        double a = 0.44*y - 0.88;
        double B = 0.183;
        double b = 1.5;
        double C = 0.019 + 0.107*y + 0.19*exp(-pow((4./y),4.));
        double c = 2.4;
        biasarray[i*psp->hm_z_slices+j] = 
          1. - A*pow(nuarray[i*psp->hm_z_slices+j],a)/(pow(nuarray[i*psp->hm_z_slices+j],a) + pow((_DELTA_C_) ,a)) + 
	       B*pow(nuarray[i*psp->hm_z_slices+j],b) + C*pow(nuarray[i*psp->hm_z_slices+j],c);	    
	/*if (j == 0)
          printf("%f, %f\n",biasarray[i*psp->hm_z_slices+j],y);*/
      }
    }
    
    
    /* For Tinker 2010 multiplicity interpolate interpolation table of Tinker 2010 paper
       and evaluate it for each redshift at the corresponding virial overdensities (0th column) values to get the appropriate parameters. */
    int n_lines = 9;
    int n_columns = 6; 
    double params_tinker[54] = {200,0.368,0.589,0.864,-0.729,-0.243,
                                300,0.363,0.585,0.922,-0.789,-0.261,
				400,0.386,0.544,0.987,-0.910,-0.261,
				600,0.389,0.543,1.090,-1.050,-0.273,
				800,0.393,0.564,1.200,-1.200,-0.278,
				1200,0.365,0.623,1.34,-1.260,-0.301,
				1600,0.379,0.637,1.50,-1.450,-0.301,
				2400,0.355,0.673,1.68,-1.500,-0.319,
				3200,0.327,0.702,1.81,-1.490,-0.336};

    double * tinker_params_vir;
    class_alloc(tinker_params_vir,psp->hm_z_slices*n_columns*sizeof(double),psp->error_message);
    
    int last_index;
    for (j=0;j<psp->hm_z_slices;j++){
      double * tinker_params;
      double tinker_params1, Delta;
      class_alloc(tinker_params,sizeof(double)*n_columns*psp->hm_z_slices,psp->error_message);
      Delta = Deltaarray[j];
      
      if (Deltaarray[j] < 200){
	printf("Warning: For z > %f virial overdensity is at %f < 200 = min. of interpolation table in Tinker 2010. Parameters for 200 are used in following computation for z > %f. No guarantee of correctness of Pk for z > %f.\n",psp->hm_zarray[j],Deltaarray[j],psp->hm_zarray[j],psp->hm_zarray[j]);
	Delta = 200;
	int c,index;
	for (c=j;c<psp->hm_z_slices;c++){
	  for (index=0;index<n_columns;index++)
	    tinker_params_vir[index*psp->hm_z_slices+c] = params_tinker[index];	
	}
       goto end; 	
      }
	    
      class_call(array_interpolate(
		   params_tinker, 
		   n_columns,
		   n_lines,
		   0,   
		   Delta,
		   &last_index,
		   tinker_params,
		   n_columns,  
		   psp->error_message),psp->error_message,psp->error_message);
      
      int index;
      for (index=0;index<n_columns;index++){
        tinker_params_vir[index*psp->hm_z_slices+j] = tinker_params[index];	
	/* printf("%f\n",tinker_params[index]);*/
      }   
    }    
    end:
    
    for (j=0;j<psp->hm_z_slices;j++){
      for (i=0;i<psp->hm_m_slices;i++){
	multiplicityarray[i*psp->hm_z_slices+j] = 
	  tinker_params_vir[1*psp->hm_z_slices+j] * 
	  (1.+ pow(tinker_params_vir[2*psp->hm_z_slices+0] * pow(aarray[j],-0.2) * nuarray[i*psp->hm_z_slices+j],-2.*tinker_params_vir[4*psp->hm_z_slices+0] * pow(aarray[j],0.08))) * 
	  pow(nuarray[i*psp->hm_z_slices+j],2.*tinker_params_vir[5*psp->hm_z_slices+0] * pow(aarray[j],-0.27)) *
	  exp(-.5*tinker_params_vir[3*psp->hm_z_slices+0] * pow(aarray[j],0.01) * pow(nuarray[i*psp->hm_z_slices+j],2));
	  /*if (j == 0)
            printf("%f, \n",multiplicityarray[i*psp->hm_z_slices+j]);*/
      }
    }
    
    free(tinker_params_vir);
    
    /* Bias up to 2nd order. 
       Note that 1st order equals 2nd order. This is exactly the approximation of Cooray/Sheth 2002 that one can do if one is solely interested in Pk
       When higher orders are being implemented this needs to be changed! */
    class_alloc(bias_pert,psp->hm_m_slices*psp->hm_z_slices*(psp->hm_highest_order+1)*sizeof(double),psp->error_message);
    for (i = 0;i<psp->hm_m_slices;i++){
      for (j = 0;j<psp->hm_z_slices;j++){
	bias_pert[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	          = 1.;
        bias_pert[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	          = biasarray[i*psp->hm_z_slices+j];
	bias_pert[2*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	          = biasarray[i*psp->hm_z_slices+j];	
      }
    }
  }
  
  /* Mass function of Boquet et al 2015 (arxiv.org/pdf/1502.07357.pdf) fitted for baryons & its bias via the beak background split */
  if (psp->hm_convention == 3 ){
    if (psp->spectra_verbose > 0)
      printf("Using mass function of Boquet 2015 and its associated pbs bias.\n");
    
    double * A_boc;
    double * a_boc;
    double * b_boc;
    double * c_boc;
    
    class_alloc(multiplicityarray,psp->hm_z_slices * psp->hm_m_slices*sizeof(double),psp->error_message);
    class_alloc(biasarray,psp->hm_m_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(A_boc,psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(a_boc,psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(b_boc,psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(c_boc,psp->hm_z_slices*sizeof(double),psp->error_message);
    
    int baryons = 0;
    
    for (i = 0;i<psp->hm_m_slices;i++){
      for (j = 0;j<psp->hm_z_slices;j++){
	if (baryons == 0){
	  /* Params for DMonly */
          A_boc[j] = 0.175 * pow(1.+ psp->hm_zarray[j],-.012);
	  a_boc[j] = 1.530 * pow(1.+ psp->hm_zarray[j],-.040);
	  b_boc[j] = 2.550 * pow(1.+ psp->hm_zarray[j],-.194);
	  c_boc[j] = 1.190 * pow(1.+ psp->hm_zarray[j],-.021);
	}
	
	else {
	  /* Params for baryon inclusion */
          A_boc[j] = 0.228 * pow(1.+ psp->hm_zarray[j], .286);
	  a_boc[j] = 2.150 * pow(1.+ psp->hm_zarray[j],-.058);
	  b_boc[j] = 1.690 * pow(1.+ psp->hm_zarray[j],-.366);
	  c_boc[j] = 1.300 * pow(1.+ psp->hm_zarray[j],-.045);
	}
        multiplicityarray[i*psp->hm_z_slices+j] = A_boc[j]*(pow(b_boc[j]*nuarray[i*psp->hm_z_slices+j]/(_DELTA_C_) ,a_boc[j])+1.) * exp(-c_boc[j]*pow(nuarray[i*psp->hm_z_slices+j]/(_DELTA_C_) ,2));
        biasarray[i*psp->hm_z_slices+j] = 1 + 2*c_boc[j]*pow(nuarray[i*psp->hm_z_slices+j],2)/pow((_DELTA_C_) ,3) - a_boc[j]/(_DELTA_C_) /(1 + pow((_DELTA_C_) /(b_boc[j]*nuarray[i*psp->hm_z_slices+j]),a_boc[j]));
        /*printf("f(%f,%f)=%f\n",log10marray[i*psp->hm_z_slices+j],psp->hm_zarray[j],multiplicityarray[i*psp->hm_z_slices+j]);
	printf("b(%f,%f)=%f\n",log10marray[i*psp->hm_z_slices+j],psp->hm_zarray[j],biasarray[i*psp->hm_z_slices+j]);*/
      }     
    }
    
    /* Bias up to 2nd order. 
       Note that 1st order equals 2nd order. This is exactly the approximation of Cooray/Sheth 2002 that one can do if one is solely interested in Pk
       When higher orders are being implemented this needs to be changed! */
    class_alloc(bias_pert,psp->hm_m_slices*psp->hm_z_slices*(psp->hm_highest_order+1)*sizeof(double),psp->error_message);
    for (i = 0;i<psp->hm_m_slices;i++){
      for (j = 0;j<psp->hm_z_slices;j++){
	bias_pert[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	          = 1.;
        bias_pert[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	          = biasarray[i*psp->hm_z_slices+j];
	bias_pert[2*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j]
	          = biasarray[i*psp->hm_z_slices+j];
      }
    }
    
    free(A_boc);
    free(a_boc);
    free(b_boc);
    free(c_boc);
  }
  
  
    
  /*  Calculate renormalization constants for multiplicity and bias factors according to Schmidt (2016) */ 
  
  class_alloc(normmultiplicityarray, psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(normfb0array,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(normfb1array,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(normfb2array,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(multipl_renormarray, psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(b0_renormarray, psp->hm_z_slices*sizeof(double),psp->error_message); 
  class_alloc(b1_renormarray, psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(b2_renormarray, psp->hm_z_slices*sizeof(double),psp->error_message);
    
  if (psp->hm_highest_order > 2){
    class_alloc(normfb3array,psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(b3_renormarray, psp->hm_z_slices*sizeof(double),psp->error_message);
  }
    
  if (psp->hm_highest_order > 3){
    class_alloc(normfb4array,psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(b4_renormarray, psp->hm_z_slices*sizeof(double),psp->error_message);
  }
    
    /* Stuff one needs for class version of spline & spline integration */
   int index_num;
   int index_k;
   int index_y;
   int index_ddy;
   int i2;
    
    /* For each redshift value obtain the fst & fbst constants */
    /* The procedure inside the loop is similar to the one in the spectra_sigma function */
    for (j=0;j<psp->hm_z_slices;j++){

     double * array_for_fst;
     double * array_for_fb0st;
     double * array_for_fb1st;
     double * array_for_fb2st;
     double * array_for_fb3st;
     double * array_for_fb4st;     
     double bst, nu, fb0st,fb1st, fb2st, fb3st, fb4st, normfst, normfb0st, normfb1st, normfb2st, normfb3st, normfb4st;

     i2=0;
     index_k=i2;
     i2++;
     index_y=i2;
     i2++;
     index_ddy=i2;
     i2++;
     index_num=i2;

     class_alloc(array_for_fst,
                 psp->hm_m_slices*index_num*sizeof(double),
                 psp->error_message);
     
     class_alloc(array_for_fb0st,
                 psp->hm_m_slices*index_num*sizeof(double),
                 psp->error_message);
     
     class_alloc(array_for_fb1st,
                 psp->hm_m_slices*index_num*sizeof(double),
                 psp->error_message);
     
     class_alloc(array_for_fb2st,
                 psp->hm_m_slices*index_num*sizeof(double),
                 psp->error_message);
     if (psp->hm_highest_order > 2)
       class_alloc(array_for_fb3st,
                 psp->hm_m_slices*index_num*sizeof(double),
                 psp->error_message);
     
     if (psp->hm_highest_order > 3)
       class_alloc(array_for_fb4st,
                 psp->hm_m_slices*index_num*sizeof(double),
                 psp->error_message);

     for (i=0;i<psp->hm_m_slices;i++) {
       nu=nuarray[i*psp->hm_z_slices+j];
       fb0st = bias_pert[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * multiplicityarray[i*psp->hm_z_slices+j];
       fb1st = bias_pert[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * multiplicityarray[i*psp->hm_z_slices+j];
       fb2st = bias_pert[2*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * multiplicityarray[i*psp->hm_z_slices+j];
       if (psp->hm_highest_order > 2)
         fb3st = bias_pert[3*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * multiplicityarray[i*psp->hm_z_slices+j];
       if (psp->hm_highest_order > 3)
         fb4st = bias_pert[4*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * multiplicityarray[i*psp->hm_z_slices+j];
      
       if (psp->hm_convention == 1){
       array_for_fst[i*index_num+index_k]=nu*nu;
       array_for_fst[i*index_num+index_y]=multiplicityarray[i*psp->hm_z_slices+j];
       array_for_fb0st[i*index_num+index_k]=nu*nu;
       array_for_fb0st[i*index_num+index_y]=fb0st;
       array_for_fb1st[i*index_num+index_k]=nu*nu;
       array_for_fb1st[i*index_num+index_y]=fb1st;
       array_for_fb2st[i*index_num+index_k]=nu*nu;
       array_for_fb2st[i*index_num+index_y]=fb2st;
       if (psp->hm_highest_order > 2){
         array_for_fb3st[i*index_num+index_k]=nu*nu;
         array_for_fb3st[i*index_num+index_y]=fb3st;
       }
       if (psp->hm_highest_order > 3){
         array_for_fb4st[i*index_num+index_k]=nu*nu;
         array_for_fb4st[i*index_num+index_y]=fb4st;
       }
       }
       
       if (psp->hm_convention == 2|| psp->hm_convention == 3){
       array_for_fst[i*index_num+index_k]=nu;
       array_for_fst[i*index_num+index_y]=multiplicityarray[i*psp->hm_z_slices+j];
       array_for_fb0st[i*index_num+index_k]=nu;
       array_for_fb0st[i*index_num+index_y]=fb0st;
       array_for_fb1st[i*index_num+index_k]=nu;
       array_for_fb1st[i*index_num+index_y]=fb1st;
       array_for_fb2st[i*index_num+index_k]=nu;
       array_for_fb2st[i*index_num+index_y]=fb2st;
       if (psp->hm_highest_order > 2){
         array_for_fb3st[i*index_num+index_k]=nu;
         array_for_fb3st[i*index_num+index_y]=fb3st;
       }
       if (psp->hm_highest_order > 3){
         array_for_fb4st[i*index_num+index_k]=nu;
         array_for_fb4st[i*index_num+index_y]=fb4st;
       }
       }
       /*printf("array_for_fst(%f, %f) = %f\n",(double)(i*index_num+index_k),(double)(j)/per_z,array_for_fst[i*index_num+index_k]);
       printf("array_for_fst(%f, %f) = %f\n",(double)(i*index_num+index_y),(double)(j)/per_z,array_for_fst[i*index_num+index_y]);*/
     }
     

     class_call(array_spline(array_for_fst,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
     
     class_call(array_integrate_all_spline(array_for_fst,
                                           index_num,
                                           psp->hm_m_slices,
                                           index_k,
                                           index_y,
                                           index_ddy,
                                           & normfst,
                                           psp->error_message),
                psp->error_message,
                psp->error_message);
     
     
      class_call(array_spline(array_for_fb0st,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
      
      class_call(array_integrate_all_spline(array_for_fb0st,
                                           index_num,
                                           psp->hm_m_slices,
                                           index_k,
                                           index_y,
                                           index_ddy,
                                           & normfb0st,
                                           psp->error_message),
                psp->error_message,
                psp->error_message);
      
      
      class_call(array_spline(array_for_fb1st,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
      
      class_call(array_integrate_all_spline(array_for_fb1st,
                                           index_num,
                                           psp->hm_m_slices,
                                           index_k,
                                           index_y,
                                           index_ddy,
                                           & normfb1st,
                                           psp->error_message),
                psp->error_message,
                psp->error_message);
      
      
      class_call(array_spline(array_for_fb2st,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
      
      class_call(array_integrate_all_spline(array_for_fb2st,
                                           index_num,
                                           psp->hm_m_slices,
                                           index_k,
                                           index_y,
                                           index_ddy,
                                           & normfb2st,
                                           psp->error_message),
                psp->error_message,
                psp->error_message);
      
      if (psp->hm_highest_order > 2){
        class_call(array_spline(array_for_fb3st,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	
	 class_call(array_integrate_all_spline(array_for_fb3st,
                                           index_num,
                                           psp->hm_m_slices,
                                           index_k,
                                           index_y,
                                           index_ddy,
                                           & normfb3st,
                                           psp->error_message),
                psp->error_message,
                psp->error_message);
      }
      
      
      if (psp->hm_highest_order > 3){
        class_call(array_spline(array_for_fb4st,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	
	class_call(array_integrate_all_spline(array_for_fb4st,
                                           index_num,
                                           psp->hm_m_slices,
                                           index_k,
                                           index_y,
                                           index_ddy,
                                           & normfb4st,
                                           psp->error_message),
                psp->error_message,
                psp->error_message);
      }


     normmultiplicityarray[j] = normfst;
     normfb0array[j] = normfb0st;
     normfb1array[j] = normfb1st;
     normfb2array[j] = normfb2st;
     multipl_renormarray[j] =  1. - normfst;
     b0_renormarray[j] = (1. - normfb0st) / (1. - normfst);
     b1_renormarray[j] = (1. - normfb1st) / (1. - normfst);
     if (highest_order_implemented == 2)
       b2_renormarray[j] = (1. - normfb2st) / (1. - normfst);
     else
       b2_renormarray[j] = (0. - normfb2st) / (1. - normfst);
       
     /*printf("normf(z = %f) = %f\n",(double)(j)/per_z,normfst);
     printf("normfb1(z = %f) = %f\n",(double)(j)/per_z,normfb1st);
     printf("normfb2(z = %f) = %f\n",(double)(j)/per_z,normfb2st);
     printf("b2_renormarray(z = %f) = %f\n",(double)(j)/per_z,b2_renormarray[j]);
     printf("multipl_renormarray(z = %f) = %f\n",(double)(j)/per_z,multipl_renormarray[j]);*/
     free(array_for_fst);
     free(array_for_fb0st);
     free(array_for_fb1st);
     free(array_for_fb2st);
     
     if (psp->hm_highest_order > 2){
     normfb3array[j] = normfb3st;
     b3_renormarray[j] = (0. - normfb3st) / (1. - normfst);
     /*printf("normfbst(z = %f) = %f\n",(double)(j)/per_z,normfb3st);*/
     free(array_for_fb3st);
     }
     if (psp->hm_highest_order > 3){
     normfb4array[j] = normfb4st;   
     b4_renormarray[j] = (0. - normfb4st) / (1. - normfst);    
     /*printf("normfbst(z = %f) = %f\n",(double)(j)/per_z,normfb4st); */
     free(array_for_fb4st);
     }
    } 
    
    free(normfb0array);
    free(normfb1array);
    free(normfb2array);
    if (psp->hm_highest_order > 2)
      free(normfb3array);
    if (psp->hm_highest_order > 3)
      free(normfb4array);
    
  int include_fideli = 0;
  
  if (include_fideli ==1){
  /* Exteded Halo model from Fideli */
  /* Fractions of masses of several components of extended Halo model */
  double rhobar_s0 = 7*pow(10,8); /* Estimate for average density of stars */
  double Omega0_s = rhobar_s0 / rhobar_marray[0];
  double rhobar_g0 = rhobar_bgarray[0]*pba->Omega0_b-rhobar_s0;
  double Omega0_g = pba->Omega0_b - Omega0_s;
  double * f_g;
  double * f_s_unnormed;
  double * f_s;
  double * m_0_g;
  double * sharpness_g;
  double * m_0_s;
  double * sharpness_s;
  class_alloc(f_g,psp->hm_m_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(f_s_unnormed,psp->hm_m_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(f_s,psp->hm_m_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(m_0_g,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(m_0_s,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(sharpness_g,psp->hm_z_slices*sizeof(double),psp->error_message);
  class_alloc(sharpness_s,psp->hm_z_slices*sizeof(double),psp->error_message);

  
  /* Fraction of dark matter wrt total matter */
  double f_DM = 1 - pba->Omega0_b/Omega0_m; 
  
  /* Fraction of gas and stars wrt total matter */
  for (i=0;i<psp->hm_m_slices;i++){
    for (j=0;j<psp->hm_z_slices;j++){
      m_0_g[j] = pow(10,12.0); /* Value where gas fraction drops */
      m_0_s[j] = 5*pow(10,12.0); /* Value where star fraction peaks */
      sharpness_g[j] = 3.0; /* Sharpness of drop */
      sharpness_s[j] = 1.2; /* Sharpness of peak */
      if (pow(10,log10marray[i*psp->hm_z_slices+j])>m_0_g[j])
        f_g[i*psp->hm_z_slices+j] = (1 - f_DM)*erf(log10(pow(10,log10marray[i*psp->hm_z_slices+j])/m_0_g[j])/sharpness_g[j]);
      else
	f_g[i*psp->hm_z_slices+j] = 0;
      f_s_unnormed[i*psp->hm_z_slices+j] = exp(-pow(log10(pow(10,log10marray[i*psp->hm_z_slices+j])/m_0_s[j]),2)/(2.*pow(sharpness_s[j],2)));
      /*printf("f_g(%f) = %f\n",log10marray[i*psp->hm_z_slices+j],f_g[i*psp->hm_z_slices+j]);   
      printf("f_s(%f) = %f\n",log10marray[i*psp->hm_z_slices+j],f_s_unnormed[i*psp->hm_z_slices+j]);  */    
    }
  }
  
  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;
  
  double nu1, int_f_s;
  double * array_for_f_s;
  class_alloc(array_for_f_s,
                 psp->hm_m_slices*index_num*sizeof(double),
                 psp->error_message);
  
  for (i=0;i<psp->hm_m_slices;i++){
    nu1 = nuarray[i*psp->hm_z_slices+0];
    if (psp->hm_convention == 1)
      array_for_f_s[i*index_num+index_k]=nu1*nu1;
    if (psp->hm_convention == 2 || psp->hm_convention == 3)
      array_for_f_s[i*index_num+index_k]=nu1;
    array_for_f_s[i*index_num+index_y]=multiplicityarray[i*psp->hm_z_slices+0] * f_s_unnormed[i*psp->hm_z_slices+0];
    /*printf("base=%f\n",array_for_f_s[i*index_num+index_k]);
    printf("int=%f\n",array_for_f_s[i*index_num+index_y]);*/
    
  }
  
  /* Find normalization constant A for f_s (see Fideli 2014) */
  class_call(array_spline(array_for_f_s,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
     
  class_call(array_integrate_all_spline(array_for_f_s,
                                           index_num,
                                           psp->hm_m_slices,
                                           index_k,
                                           index_y,
                                           index_ddy,
                                           & int_f_s,
                                           psp->error_message),
                psp->error_message,
                psp->error_message);
  
  double Amplitude_s = rhobar_s0/rhobar_marray[0] / int_f_s;
  /*printf("A=%f,int_f_s = %f",Amplitude_s,int_f_s);*/
  
   for (i=0;i<psp->hm_m_slices;i++){
    for (j=0;j<psp->hm_z_slices;j++){
      f_s[i*psp->hm_z_slices+j] = Amplitude_s * f_s_unnormed[i*psp->hm_z_slices+j];
      /*printf("star fraction (%f,%f) = %f\n",log10marray[i*psp->hm_z_slices+j],psp->hm_zarray[j],f_s[i*psp->hm_z_slices+j]);
      printf("gas fraction (%f,%f) = %f\n",log10marray[i*psp->hm_z_slices+j],psp->hm_zarray[j],f_g[i*psp->hm_z_slices+j]);*/
    }
   }
    
  /* Fourier transformed gas profile (need to normalize!)*/
  double * profile_gas;
  class_alloc(profile_gas,psp->hm_k_slices*psp->hm_z_slices*psp->hm_m_slices*sizeof(double),psp->error_message);
  /* Fourier transform of gas profile (for simplicity we set \beta = 3/2) */
  for (i=0;i<psp->hm_m_slices;i++){
    for (j=0;j<psp->hm_z_slices;j++){
      for (l=0;l<psp->hm_k_slices;l++){
	profile_gas[i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] = pow(_PI_/2.,.5)*exp(-fabs(pow(10,psp->log10karray[l])));
        /*printf("profile(%f,%f,%f)=%f",log10marray[i*psp->hm_z_slices+j],psp->hm_zarray[j],psp->log10karray[l],profile_gas[i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l]);*/	
      }
    }
  }
  }
  

    /* Define integrands for 'matrix' M_ij from which components the Fourier transforms of the n point cfs are built 
       --> As no sums of vectors are in these we only need to distribute the leading highest_order - elements of 
           rel_lengtharray to j free spots of the M_ij --> higest_order^j possibilities*/
    /* In the last step construct the renomalized matrix element according to Schmidt 2016 with help of the constants defined before */
    int l_2, l_3, l_4; /* Subindex gives order of Spectrum where first appears. Eg l_2 appears in Pk,l_3 in bispectrum, ... */
    double nu;
    double I_M_00, I_M_10, I_M_20, I_M_30, I_M_40, I_M_01, I_M_11, I_M_21, I_M_31, I_M_02, I_M_12, I_M_22, I_M_03, I_M_13, I_M_04;
    double M00, M10, M20, M30, M40;
    double M01, M11, M21, M31;
    double M02, M12, M22;
    double M03, M13;
    double M04;
               
    class_alloc(M_00,psp->hm_z_slices*sizeof(double),psp->error_message); 
    class_alloc(M_10,psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(M_20,psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(M_30,psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(M_40,psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(M_01,psp->hm_z_slices*psp->hm_k_slices*number_vectors*sizeof(double),psp->error_message); 
    class_alloc(M_11,psp->hm_z_slices*psp->hm_k_slices*number_vectors*sizeof(double),psp->error_message); 
    class_alloc(M_21,psp->hm_z_slices*psp->hm_k_slices*number_vectors*sizeof(double),psp->error_message); 
    class_alloc(M_31,psp->hm_z_slices*psp->hm_k_slices*number_vectors*sizeof(double),psp->error_message); 
    class_alloc(M_02,psp->hm_z_slices*psp->hm_k_slices*pow(number_vectors,2)*sizeof(double),psp->error_message);
    class_alloc(M_12,psp->hm_z_slices*psp->hm_k_slices*pow(number_vectors,2)*sizeof(double),psp->error_message); 
    class_alloc(M_22,psp->hm_z_slices*psp->hm_k_slices*pow(number_vectors,2)*sizeof(double),psp->error_message); 
    
    if (psp->hm_highest_order > 2){
      class_alloc(M_03,psp->hm_z_slices*psp->hm_k_slices*pow(number_vectors,3)*sizeof(double),psp->error_message); 
      class_alloc(M_13,psp->hm_z_slices*psp->hm_k_slices*pow(number_vectors,3)*sizeof(double),psp->error_message); 
    }
    
    if (psp->hm_highest_order > 3){
      class_alloc(M_04,psp->hm_z_slices*psp->hm_k_slices*pow(number_vectors,4)*sizeof(double),psp->error_message); 
    }
    
    for (j=0;j<psp->hm_z_slices;j++){
      int order_1, order_2, order_3, order_4;
      
      int index_num;
      int index_k;
      int index_y;
      int index_ddy;
      int i1;
	      
      i1=0;
      index_k=i1;
      i1++;
      index_y=i1;
      i1++;
      index_ddy=i1;
      i1++;
      index_num=i1;
      
      double * array_for_M_00;
      double * array_for_M_10;
      double * array_for_M_20;
      double * array_for_M_30;
      double * array_for_M_40;
      
      class_alloc(array_for_M_00,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
	      
      class_alloc(array_for_M_10,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
	      
      class_alloc(array_for_M_20,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
      
      for (i=0;i<psp->hm_m_slices;i++){
	nu = nuarray[i*psp->hm_z_slices+j];
	I_M_00 = 0.;
	I_M_10 = 0.;
	I_M_20 = multiplicityarray[i*psp->hm_z_slices+j] * pow( pow(10,log10marray[i*psp->hm_z_slices+j])/rhobar_marray[j],-1) * bias_pert[2*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] ;
	/* printf("I_M_20(%f,%f) = %f\n",nuarray[i],(double)(j)/per_z,I_M_20);*/
      
	if (psp->hm_convention == 1){
          array_for_M_00[i*index_num+index_k] = nu*nu;
          array_for_M_00[i*index_num+index_y] = I_M_00;
          array_for_M_10[i*index_num+index_k] = nu*nu;
          array_for_M_10[i*index_num+index_y] = I_M_10;
          array_for_M_20[i*index_num+index_k] = nu*nu;
          array_for_M_20[i*index_num+index_y] = I_M_20;
	}
	
	if (psp->hm_convention == 2 || psp->hm_convention == 3){
          array_for_M_00[i*index_num+index_k] = nu;
          array_for_M_00[i*index_num+index_y] = I_M_00;
          array_for_M_10[i*index_num+index_k] = nu;
          array_for_M_10[i*index_num+index_y] = I_M_10;
          array_for_M_20[i*index_num+index_k] = nu;
          array_for_M_20[i*index_num+index_y] = I_M_20;
	}
	
	
      }
      
      
      class_call(array_spline(array_for_M_00,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
      class_call(array_integrate_all_spline(array_for_M_00,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M00,
                              psp->error_message),
                psp->error_message,
                psp->error_message);  
      
      
      class_call(array_spline(array_for_M_10,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
      class_call(array_integrate_all_spline(array_for_M_10,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M10,
                              psp->error_message),
                psp->error_message,
                psp->error_message);  
      
      
      class_call(array_spline(array_for_M_20,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
      class_call(array_integrate_all_spline(array_for_M_20,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M20,
                              psp->error_message),
                psp->error_message,
                psp->error_message); 
      
      M_00[j] = M00 +  multipl_renormarray[j] * b0_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],0-1);
      M_10[j] = M10 + multipl_renormarray[j] * b1_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],0-1);
      M_20[j] = M20 + multipl_renormarray[j] * b2_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],0-1);
      
      free(array_for_M_00);
      free(array_for_M_10);
      free(array_for_M_20);

      if (psp->hm_highest_order > 2){
        class_alloc(array_for_M_30,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
	
	for (i=0;i<psp->hm_m_slices;i++){
		nu = nuarray[i*psp->hm_z_slices+j];
		I_M_30 = multiplicityarray[i*psp->hm_z_slices+j] * pow( pow(10,log10marray[i*psp->hm_z_slices+j])/rhobar_marray[j],-1) * bias_pert[3*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] ; 
	        if (psp->hm_convention == 1){
		  array_for_M_30[i*index_num+index_k] = nu*nu;
                  array_for_M_30[i*index_num+index_y] = I_M_30;
		}
		if (psp->hm_convention == 2|| psp->hm_convention == 3){
		  array_for_M_30[i*index_num+index_k] = nu;
                  array_for_M_30[i*index_num+index_y] = I_M_30;
		}
	}
	
	class_call(array_spline(array_for_M_30,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
        class_call(array_integrate_all_spline(array_for_M_30,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M30,
                              psp->error_message),
                psp->error_message,
                psp->error_message);  
	
	M_30[j] = 4*M30 + multipl_renormarray[j] * b3_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],0-1);
        
	free(array_for_M_30);
      }
      
      
      if (psp->hm_highest_order > 3){
        class_alloc(array_for_M_40,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
        for (i=0;i<psp->hm_m_slices;i++){
	  nu = nuarray[i*psp->hm_z_slices+j];
          I_M_40 = multiplicityarray[i*psp->hm_z_slices+j] * pow( pow(10,log10marray[i*psp->hm_z_slices+j])/rhobar_marray[j],-1) * bias_pert[4*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] ;
          if (psp->hm_convention == 1){
	  array_for_M_40[i*index_num+index_k] = nu*nu;
          array_for_M_40[i*index_num+index_y] = I_M_40;
	  }
	  if (psp->hm_convention == 2 || psp->hm_convention == 3){
	  array_for_M_40[i*index_num+index_k] = nu;
          array_for_M_40[i*index_num+index_y] = I_M_40;
	  }
	}
	
	class_call(array_spline(array_for_M_40,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
        class_call(array_integrate_all_spline(array_for_M_40,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M40,
                              psp->error_message),
                psp->error_message,
                psp->error_message);  
	
	 M_40[j] = M40 + 4*multipl_renormarray[j] * b4_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],0-1);
	
	 free(array_for_M_40);	    
    }
      for (order_1 =0; order_1<number_vectors;order_1++){
      for (l=0;l<psp->hm_k_slices;l++){
	
	double * array_for_M_01;
	double * array_for_M_11;
	double * array_for_M_21;
	double * array_for_M_31;
	
	
        class_alloc(array_for_M_01,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
	      
	class_alloc(array_for_M_11,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
	      
        class_alloc(array_for_M_21,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
	
	
	
	for (i=0;i<psp->hm_m_slices;i++){
	  nu = nuarray[i*psp->hm_z_slices+j];
	  I_M_01 = 0.; 
	  I_M_11 = multiplicityarray[i*psp->hm_z_slices+j] * pow( pow(10,log10marray[i*psp->hm_z_slices+j])/rhobar_marray[j],0) * bias_pert[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l]; 
	  I_M_21 = multiplicityarray[i*psp->hm_z_slices+j] * pow( pow(10,log10marray[i*psp->hm_z_slices+j])/rhobar_marray[j],0) * bias_pert[2*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];        
	  /*printf("I_M_11(%f,%f,%f) = %f\n",nuarray[i*psp->hm_z_slices+j],pow(10,psp->log10karray[l]),(double)(j)/per_z,I_M_11);*/
	  if (psp->hm_convention == 1){
	    array_for_M_01[i*index_num+index_k] = nu*nu;
	    array_for_M_01[i*index_num+index_y] = I_M_01;
       	    array_for_M_11[i*index_num+index_k] = nu*nu;
	    array_for_M_11[i*index_num+index_y] = I_M_11;
	    array_for_M_21[i*index_num+index_k] = nu*nu;
	    array_for_M_21[i*index_num+index_y] = I_M_21;
	  }
	  if (psp->hm_convention == 2 || psp->hm_convention == 3){
	    array_for_M_01[i*index_num+index_k] = nu;
	    array_for_M_01[i*index_num+index_y] = I_M_01;
       	    array_for_M_11[i*index_num+index_k] = nu;
	    array_for_M_11[i*index_num+index_y] = I_M_11;
	    array_for_M_21[i*index_num+index_k] = nu;
	    array_for_M_21[i*index_num+index_y] = I_M_21;
	  }
	}
	
	class_call(array_spline(array_for_M_01,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
        class_call(array_integrate_all_spline(array_for_M_01,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M01,
                              psp->error_message),
                psp->error_message,
                psp->error_message);  
      
      
        class_call(array_spline(array_for_M_11,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
        class_call(array_integrate_all_spline(array_for_M_11,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M11,
                              psp->error_message),
                psp->error_message,
                psp->error_message);  
      
      
        class_call(array_spline(array_for_M_21,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
        class_call(array_integrate_all_spline(array_for_M_21,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M21,
                              psp->error_message),
                psp->error_message,
                psp->error_message); 
      
	M_01[order_1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] = M01
	                               + multipl_renormarray[j] * b0_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],1-1) * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];
        M_11[order_1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] = M11
                                       + multipl_renormarray[j] * b1_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],1-1) * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];
        M_21[order_1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] = M21
	                               + multipl_renormarray[j] * b2_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],1-1) * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];   
        /*printf("M_21[%d,%d] = %f\n" ,order_1,l, M_21[order_1*0*psp->hm_k_slices+0*psp->hm_k_slices + l]);
	printf("M_11[%d,%d] = %f\n" ,order_1,l, M_11[order_1*0*psp->hm_k_slices+0*psp->hm_k_slices + l]);*/
	         
        free(array_for_M_01);
        free(array_for_M_11);
        free(array_for_M_21);
        

	if (psp->hm_highest_order > 2){
	  
          class_alloc(array_for_M_31,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
	
	  for (i=0;i<psp->hm_m_slices;i++){
	    nu = nuarray[i*psp->hm_z_slices+j];
	    I_M_31 = multiplicityarray[i*psp->hm_z_slices+j] * pow( pow(10,log10marray[i*psp->hm_z_slices+j])/rhobar_marray[j],0) * bias_pert[3*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];	  
	    if (psp->hm_convention == 1){
	    array_for_M_31[i*index_num+index_k] = nu*nu;
	    array_for_M_31[i*index_num+index_y] = I_M_31;	
	    }
	    if (psp->hm_convention == 2 || psp->hm_convention == 3){
	    array_for_M_31[i*index_num+index_k] = nu;
	    array_for_M_31[i*index_num+index_y] = I_M_31;	
	    }
	  }
	  
	  class_call(array_spline(array_for_M_31,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
          class_call(array_integrate_all_spline(array_for_M_31,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M31,
                              psp->error_message),
                psp->error_message,
                psp->error_message); 
	
	   M_31[order_1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] = M31;
	                               + multipl_renormarray[j] * b3_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],0-1) * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];
     
	   free(array_for_M_31);
	}
	
	
	for (order_2=0;order_2<number_vectors;order_2++){
	  
	  double * array_for_M_02;
	  double * array_for_M_12;
	  double * array_for_M_22;
	  
	  
	  class_alloc(array_for_M_02,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
	      
	  class_alloc(array_for_M_12,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
	      
	  class_alloc(array_for_M_22,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
	  
	  
	  for(i=0;i<psp->hm_m_slices;i++){
	    nu = nuarray[i*psp->hm_z_slices+j];
	    I_M_02 = multiplicityarray[i*psp->hm_z_slices+j] * pow( pow(10,log10marray[i*psp->hm_z_slices+j])/rhobar_marray[j],1) * bias_pert[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_2*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] /*- pow(1.+pow(pow(10,psp->log10karray[l]+psp->log10karray[l_2])*m_to_Rarray[i*psp->hm_z_slices+j]/2.,4),-1)*pow(10,log10marray[i*psp->hm_z_slices+j])/rhobar_marray[j]*/;
            I_M_12 = multiplicityarray[i*psp->hm_z_slices+j] * pow( pow(10,log10marray[i*psp->hm_z_slices+j])/rhobar_marray[j],1) * bias_pert[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_2*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l]; 
	    I_M_22 = multiplicityarray[i*psp->hm_z_slices+j] * pow( pow(10,log10marray[i*psp->hm_z_slices+j])/rhobar_marray[j],1) * bias_pert[2*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_2*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l]; 
	  
	    if (psp->hm_convention == 1){
	    array_for_M_02[i*index_num+index_k] = nu*nu;
	    array_for_M_02[i*index_num+index_y] = I_M_02;
       	    array_for_M_12[i*index_num+index_k] = nu*nu;
	    array_for_M_12[i*index_num+index_y] = I_M_12;
	    array_for_M_22[i*index_num+index_k] = nu*nu;
	    array_for_M_22[i*index_num+index_y] = I_M_22;
	    }
	    if (psp->hm_convention == 2 || psp->hm_convention == 3){
	    array_for_M_02[i*index_num+index_k] = nu;
	    array_for_M_02[i*index_num+index_y] = I_M_02;
       	    array_for_M_12[i*index_num+index_k] = nu;
	    array_for_M_12[i*index_num+index_y] = I_M_12;
	    array_for_M_22[i*index_num+index_k] = nu;
	    array_for_M_22[i*index_num+index_y] = I_M_22;
	    }
	    
	  }
	    
	  class_call(array_spline(array_for_M_02,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
          class_call(array_integrate_all_spline(array_for_M_02,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M02,
                              psp->error_message),
                psp->error_message,
                psp->error_message);  
      
      
          class_call(array_spline(array_for_M_12,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
          class_call(array_integrate_all_spline(array_for_M_12,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M12,
                              psp->error_message),
                psp->error_message,
                psp->error_message); 
      
          class_call(array_spline(array_for_M_22,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
          class_call(array_integrate_all_spline(array_for_M_22,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M22,
                              psp->error_message),
                psp->error_message,
                psp->error_message); 
	
	  
	M_02[order_1*number_vectors*psp->hm_z_slices*psp->hm_k_slices+order_2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] = M02
	                                                                     + 4*multipl_renormarray[j] * b0_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],1) * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_2*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];
        M_12[order_1*number_vectors*psp->hm_z_slices*psp->hm_k_slices+order_2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] = M12
                                                                             + 4*multipl_renormarray[j] * b1_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],1) * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_2*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];
        M_22[order_1*number_vectors*psp->hm_z_slices*psp->hm_k_slices+order_2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] = M22
                                                                             + 4*multipl_renormarray[j] * b2_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],1) * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_2*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];  
	/*printf("M_12[%d] = %f\n" ,0*number_vectors*psp->hm_z_slices*psp->hm_k_slices+1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l, M_12[0*number_vectors*psp->hm_z_slices*psp->hm_k_slices+1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l]); */ 
	    
	free(array_for_M_02);
	free(array_for_M_12);
	free(array_for_M_22);
	
	
	  
	  if (psp->hm_highest_order > 2){
	  for (order_3=0;order_3<number_vectors;order_3++){
	    
	  double * array_for_M_03;
	  double * array_for_M_13;
	  
	  
	  class_alloc(array_for_M_03,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
	      
	  class_alloc(array_for_M_13,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);	     
	  
	  
	  for(i=0;i<psp->hm_m_slices;i++){
	    nu = nuarray[i*psp->hm_z_slices+j];
	    I_M_03 = multiplicityarray[i*psp->hm_z_slices+j] * pow( pow(10,log10marray[i*psp->hm_z_slices+j])/rhobar_marray[j],2) * bias_pert[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_2*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_3*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];
	    I_M_13 = multiplicityarray[i*psp->hm_z_slices+j] * pow( pow(10,log10marray[i*psp->hm_z_slices+j])/rhobar_marray[j],2) * bias_pert[1*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_2*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_3*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];
	  
	    if (psp->hm_convention == 1){
	      array_for_M_03[i*index_num+index_k] = nu*nu;
	      array_for_M_03[i*index_num+index_y] = I_M_03;
       	      array_for_M_13[i*index_num+index_k] = nu*nu;
	      array_for_M_13[i*index_num+index_y] = I_M_13;
	    }
	    if (psp->hm_convention == 2 || psp->hm_convention == 3){
	      array_for_M_03[i*index_num+index_k] = nu;
	      array_for_M_03[i*index_num+index_y] = I_M_03;
       	      array_for_M_13[i*index_num+index_k] = nu;
	      array_for_M_13[i*index_num+index_y] = I_M_13;
	    }
	  }
	    
	  class_call(array_spline(array_for_M_03,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
          class_call(array_integrate_all_spline(array_for_M_03,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M03,
                              psp->error_message),
                psp->error_message,
                psp->error_message);  
      
      
          class_call(array_spline(array_for_M_13,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
          class_call(array_integrate_all_spline(array_for_M_13,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M13,
                              psp->error_message),
                psp->error_message,
                psp->error_message); 
	  
          
	M_03[(order_1*number_vectors+order_2)*number_vectors*psp->hm_z_slices*psp->hm_k_slices+order_3*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] = M03
	     + 4*multipl_renormarray[j] * b0_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],2) * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_2*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_3*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];
        M_13[(order_1*number_vectors+order_2)*number_vectors*psp->hm_z_slices*psp->hm_k_slices+order_3*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] = M13
             + 4*multipl_renormarray[j] * b1_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],2) * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_2*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_3*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];
	/*printf("M_13[%d] = %f\n",(order_1*number_vectors+order_2)*number_vectors*psp->hm_z_slices*psp->hm_k_slices+order_3*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l,M_13[(order_1*number_vectors+order_2)*number_vectors*psp->hm_z_slices*psp->hm_k_slices+order_3*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l]);
	printf("M_03[%d] = %f\n",(order_1*number_vectors+order_2)*number_vectors*psp->hm_z_slices*psp->hm_k_slices+order_3*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l,M_03[(order_1*number_vectors+order_2)*number_vectors*psp->hm_z_slices*psp->hm_k_slices+order_3*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l]);*/
	        
	free(array_for_M_03);
	free(array_for_M_13);
	    
	    if (psp->hm_highest_order > 3){
	    for (order_4=0;order_4<number_vectors;order_4++){

	      double * array_for_M_04;;
     	      	     
	      class_alloc(array_for_M_04,
                    psp->hm_m_slices*index_num*sizeof(double),
                    psp->error_message);
	      
	      for (i=0;i<psp->hm_m_slices;i++){
		nu = nuarray[i*psp->hm_z_slices+j];		
		I_M_04 = multiplicityarray[i*psp->hm_z_slices+j] * pow( pow(10,log10marray[i*psp->hm_z_slices+j])/rhobar_marray[j],3) * bias_pert[0*psp->hm_m_slices*psp->hm_z_slices + i*psp->hm_z_slices + j] * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_2*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_3*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_4*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+i*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];
		if (psp->hm_convention == 1){
		  array_for_M_04[i*index_num+index_k] = nu*nu;
		  array_for_M_04[i*index_num+index_y] = I_M_04;
		}
		if (psp->hm_convention == 2 || psp->hm_convention == 3){
		  array_for_M_04[i*index_num+index_k] = nu;
		  array_for_M_04[i*index_num+index_y] = I_M_04;
		}
	    }

              class_call(array_spline(array_for_M_04,
                             index_num,
                             psp->hm_m_slices,
                             index_k,
                             index_y,
                             index_ddy,
                             _SPLINE_EST_DERIV_,
                             psp->error_message),
                psp->error_message,
                psp->error_message);
	  
              class_call(array_integrate_all_spline(array_for_M_04,
                              index_num,
                              psp->hm_m_slices,
                              index_k,
                              index_y,
                              index_ddy,
                              &M04,
                              psp->error_message),
                psp->error_message,
                psp->error_message); 
	  
          
	      M_04[(order_1*number_vectors+order_2)*number_vectors*number_vectors*psp->hm_z_slices*psp->hm_k_slices+(order_3*number_vectors+order_4)*psp->hm_z_slices*psp->hm_k_slices+j + l] = M04
                  + multipl_renormarray[j] * b0_renormarray[j] * pow(pow(10,log10marray[j])/rhobar_marray[j],3) * ukmarray[order_1*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_2*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_3*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l] * ukmarray[order_4*psp->hm_m_slices*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices + j*psp->hm_k_slices+l];
	     free(array_for_M_04);
	     
	    }
	  }
	}	   
        }
      }
    } 
  }
    }
  
  free(b0_renormarray);
  free(b1_renormarray);
  free(b2_renormarray);
  if (psp->hm_highest_order > 2)
    free(b3_renormarray);
  if (psp->hm_highest_order > 3)
    free(b4_renormarray);
  
  /* Check:Find maximum of arrays */
  /*double max_M_04 = 0;
  for (j=0;j<psp->hm_z_slices;j++){
    for (l=0;j<psp->hm_k_slices;l++){
      for (l_2=0;j<psp->hm_k_slices;l_2++){
	for (l_3=0;l_3<psp->hm_k_slices;l_3++){
	  for (l_4=0;l_3<psp->hm_k_slices;l_4++){
	    if (M_04[j*psp->hm_k_slices*psp->hm_k_slices*psp->hm_k_slices*psp->hm_k_slices + l*psp->hm_k_slices*psp->hm_k_slices*psp->hm_k_slices+l_2*psp->hm_k_slices*psp->hm_k_slices+l_3*psp->hm_k_slices+l_4] > max_M_04)
	      max_M_04 = M_04[j*psp->hm_k_slices*psp->hm_k_slices*psp->hm_k_slices*psp->hm_k_slices + l*psp->hm_k_slices*psp->hm_k_slices*psp->hm_k_slices+l_2*psp->hm_k_slices*psp->hm_k_slices+l_3*psp->hm_k_slices+l_4];
	  }
	}
      }
    }
  }
  printf("%f\n",max_M_04);*/
  
  /* Fill arrays of Power - Bi - and Trispectrum as defined in Sheth review */
  /* If necessary, pick just relevant elements of M_ij arrays (order: j.l.l_2,l_3,l_4) eg if two k's are same in M_02 set l_2 := l */
  
  class_alloc(psp->pk_linarray,psp->hm_k_slices*psp->hm_z_slices*sizeof(double),psp->error_message); 
  class_alloc(P_1harray,psp->hm_z_slices*pow(psp->hm_k_slices,1)*sizeof(double),psp->error_message);
  class_alloc(P_2harray,psp->hm_z_slices*pow(psp->hm_k_slices,1)*sizeof(double),psp->error_message);
  
  if (psp->hm_highest_order > 2){
    class_alloc(psp->bk_linarray,pow(number_vectors,psp->hm_highest_order)*psp->hm_k_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(psp->bk_fixlinarray,pow(number_vectors,psp->hm_highest_order)*psp->hm_k_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(psp->Q_123array,pow(number_vectors,psp->hm_highest_order)*psp->hm_k_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(psp->Q_123fixarray,psp->hm_k_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(psp->Q_123fixlinarray,psp->hm_k_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
    class_alloc(B_1harray,psp->hm_z_slices*pow(psp->hm_k_slices,1)*sizeof(double),psp->error_message);
    class_alloc(B_2harray,psp->hm_z_slices*pow(psp->hm_k_slices,1)*sizeof(double),psp->error_message);
    class_alloc(B_3harray,psp->hm_z_slices*pow(psp->hm_k_slices,1)*sizeof(double),psp->error_message);	
  }
  
  if (psp->hm_highest_order > 3){
    class_alloc(psp->tk_linarray,psp->hm_k_slices*psp->hm_z_slices*sizeof(double),psp->error_message); /* Note: This is just ok for squares! */
    class_alloc(T_1harray,psp->hm_z_slices*pow(psp->hm_k_slices,1)*sizeof(double),psp->error_message);
    class_alloc(T_2harray,psp->hm_z_slices*pow(psp->hm_k_slices,1)*sizeof(double),psp->error_message);
    class_alloc(T_3harray,psp->hm_z_slices*pow(psp->hm_k_slices,1)*sizeof(double),psp->error_message);
    class_alloc(T_4harray,psp->hm_z_slices*pow(psp->hm_k_slices,1)*sizeof(double),psp->error_message);
  }
  
  
  /* Fill entries of Power, Bi and Trispectra */
  /* Note that in deafault mode we just look at equilateral triangles for Bispectrum & Squares for the Trispectrum */
  /* If in input specific vectors were chosen we proceed as follows:   
   * i)   Get norms of vectors and find the longest vector needed for Pk calculation of all spectra calculated
       --> log10k array allocated st no problems appear in calling pk
     ii)  calculated all the M_ij's with its norm wrt longest vector 
       --> just factor pow(number_vectors,j) times longer as if just one config had been hard coded
     iii) Caluclated the higher order Spectra in terms of the 3 vectors
       --> Here some more memory will be needed --> Factor: pow(number_possibilities (<=2*psp->hm_highest_order),number_spots (<=psp->hm_highest_order)) 
       --> Reflects possibilities of ordering the number_possibilities vectors in positions of quanitity with number_spots free spots */
  
  int R_of_rhobar_m = 0;
  for(j=0;j<psp->hm_z_slices;j++){
    
    /* Find element of our massarray nearest to background matter density  */
    for (i=0;i<psp->hm_m_slices;i++){
      if (pow(10,log10marray[i]) >= rhobar_marray[j]){      
        R_of_rhobar_m = i;
	break;
      }
    }
    
    for(l=0;l<psp->hm_k_slices;l++){
      double k = pow(10,psp->log10karray[l]);
      double z = (double)(j)/per_z;
      /* To cut off 1 halo terms at large scales. R_lin value related to eq (32) from Schmidt 2016 */
      double cut_1h = pow(tanh(Deltaarray[j]*m_to_Rarray[R_of_rhobar_m*psp->hm_z_slices+j]* pow(10,psp->log10karray[l])),2); 
      double * pk_ic = NULL;
      
      /* For linear power spectrum at k */
      double pk;
      class_call(spectra_pk_at_k_and_z(pba,ppm,psp,k,z,&pk,pk_ic),
               psp->error_message,
               psp->error_message);

      
      P_1harray[j*psp->hm_k_slices+l] = (M_02[index_order_max*number_vectors*psp->hm_z_slices*psp->hm_k_slices+index_order_max*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l]) * cut_1h;
      P_2harray[j*psp->hm_k_slices+l] = pow(M_11[index_order_max*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l],2) * pk  ;
      
      /*printf("%f\n",P_1harray[j*psp->hm_k_slices+l]);*/
      
      if (psp->hm_highest_order > 2){
        /* Linear & Reduced Bispectrum: Note that we allocate in each (k,z)-element a (symmetric) number_vectors**highest_order matrix reflecting the choice of ordering the different k-values */
	double bk;
	int slot_1, slot_2, slot_3;
	double pk_1,pk_2,pk_3;
	
	for (slot_1 = 0; slot_1<number_vectors; slot_1++){
	  class_call(spectra_pk_at_k_and_z(pba,ppm,psp,rel_lengtharray[slot_1]*k,z,&pk_1,pk_ic),psp->error_message,psp->error_message);
	  for (slot_2 = 0; slot_2<number_vectors; slot_2++){
	    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,rel_lengtharray[slot_2]*k,z,&pk_2,pk_ic),psp->error_message,psp->error_message);
	    for (slot_3 = 0; slot_3 < number_vectors; slot_3++){
              class_call(spectra_pk_at_k_and_z(pba,ppm,psp,rel_lengtharray[slot_3]*k,z,&pk_3,pk_ic),psp->error_message,psp->error_message);	      
	      psp->bk_linarray[(j*psp->hm_k_slices+l)*number_vectors*number_vectors*number_vectors+(slot_1*number_vectors+slot_2)*number_vectors+slot_3]
	        = 2. * ( Fs_2array[j*number_vectors*number_vectors+slot_1*number_vectors+slot_2] * pk_1 * pk_2
	               + Fs_2array[j*number_vectors*number_vectors+slot_1*number_vectors+slot_3] * pk_1 * pk_3
	               + Fs_2array[j*number_vectors*number_vectors+slot_2*number_vectors+slot_3] * pk_2 * pk_3);	      
	      
	      psp->Q_123array[(j*psp->hm_k_slices+l)*number_vectors*number_vectors*number_vectors+(slot_1*number_vectors+slot_2)*number_vectors+slot_3]
	        = psp->bk_linarray[(j*psp->hm_k_slices+l)*number_vectors*number_vectors*number_vectors+(slot_1*number_vectors+slot_2)*number_vectors+slot_3] 
	        / (pk_1 * pk_2 + pk_1 * pk_3 + pk_2 * pk_3);
	      if (j == 0 && slot_1 == 0 && slot_2 == 1 && slot_3 == 2){
	        /*printf("bk_linarray[k,0][1,2,3] = %f\n",psp->bk_linarray[(0*psp->hm_k_slices+l)*number_vectors*number_vectors*number_vectors+(0*number_vectors+1)*number_vectors+2] );*/
		/*printf("Q_123array[k,0][1,2,3] = %f\n",psp->Q_123array[(0*psp->hm_k_slices+l)*number_vectors*number_vectors*number_vectors+(0*number_vectors+1)*number_vectors+2] );*/
	      }
	    }
	  }
	}
	
	/* For output: Allocate bk_lin and Q_123 / Q_123 for fixed vectors 1,2,3 */
	double pk_11, pk_12, pk_13;
	class_call(spectra_pk_at_k_and_z(pba,ppm,psp,rel_lengtharray[0]*k,z,&pk_11,pk_ic),psp->error_message,psp->error_message);
	class_call(spectra_pk_at_k_and_z(pba,ppm,psp,rel_lengtharray[1]*k,z,&pk_12,pk_ic),psp->error_message,psp->error_message);
	class_call(spectra_pk_at_k_and_z(pba,ppm,psp,rel_lengtharray[2]*k,z,&pk_13,pk_ic),psp->error_message,psp->error_message);
		
	psp->bk_fixlinarray[j*psp->hm_k_slices+l] = psp->bk_linarray[(j*psp->hm_k_slices+l)*number_vectors*number_vectors*number_vectors+(0*number_vectors+1)*number_vectors+2];
	psp->Q_123fixlinarray[j*psp->hm_k_slices+l] = psp->bk_fixlinarray[j*psp->hm_k_slices+l]/(pk_11 * pk_12 + pk_11 * pk_13 + pk_12*pk_13);
	
	B_1harray[j*psp->hm_k_slices+l] = M_03[(0*number_vectors+1)*number_vectors*psp->hm_z_slices*psp->hm_k_slices+2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * pow(cut_1h,2) ;
        B_2harray[j*psp->hm_k_slices+l] = M_11[0*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * M_12[1*number_vectors*psp->hm_z_slices*psp->hm_k_slices+2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * cut_1h * pk_11
                                        + M_11[1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * M_12[2*number_vectors*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * cut_1h * pk_12
                                        + M_11[2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * M_12[0*number_vectors*psp->hm_z_slices*psp->hm_k_slices+1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * cut_1h * pk_13;
        B_3harray[j*psp->hm_k_slices+l] = M_11[0*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * M_11[1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * M_11[2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] 
                                          * psp->bk_linarray[(j*psp->hm_k_slices+l)*number_vectors*number_vectors*number_vectors+(0*number_vectors+1)*number_vectors+2]
                                        + M_11[0*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * M_11[1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * M_21[2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * pk_11 * pk_12
                                        + M_11[1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * M_11[2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * M_21[0*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * pk_12 * pk_13
                                        + M_11[2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * M_11[0*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * M_21[1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l] * pk_13 * pk_11;
        
	/* Reduced Bispectrum where we use Halo Model Bi - and Powerspectrum */
	/* Where does the 3 come from? */
	psp->Q_123fixarray[j*psp->hm_k_slices+l] = (B_1harray[j*psp->hm_k_slices+l] + B_2harray[j*psp->hm_k_slices+l] + B_3harray[j*psp->hm_k_slices+l])
	                                     / (
					       ( (M_02[0*number_vectors*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l]) * cut_1h + pow(M_11[0*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l],2) * pk_11 )
					       *((M_02[1*number_vectors*psp->hm_z_slices*psp->hm_k_slices+1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l]) * cut_1h + pow(M_11[1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l],2) * pk_12 )
					      +( (M_02[1*number_vectors*psp->hm_z_slices*psp->hm_k_slices+1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l]) * cut_1h + pow(M_11[1*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l],2) * pk_12 )
					       *((M_02[2*number_vectors*psp->hm_z_slices*psp->hm_k_slices+2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l]) * cut_1h + pow(M_11[2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l],2) * pk_13 ) 
					      +( (M_02[2*number_vectors*psp->hm_z_slices*psp->hm_k_slices+2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l]) * cut_1h + pow(M_11[2*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l],2) * pk_13 )
					       *((M_02[0*number_vectors*psp->hm_z_slices*psp->hm_k_slices+0*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l]) * cut_1h + pow(M_11[0*psp->hm_z_slices*psp->hm_k_slices+j*psp->hm_k_slices + l],2) * pk_11 )
					     );
        /*printf("B_1h_new(%f,%f)=%f\n",k,z,B_1harray[j*psp->hm_k_slices+l]);
        printf("B_2h_new(%f,%f)=%f\n",k,z,B_2harray[j*psp->hm_k_slices+l]);
        printf("B_3h_new(%f,%f)=%f\n",k,z,B_3harray[j*psp->hm_k_slices+l]);*/
	if (j==0){	  
	  printf("Q_lin(%d) = %f\n",l,psp->Q_123fixlinarray[j*psp->hm_k_slices+l]);
	  printf("Q_nl = %f\n",psp->Q_123fixarray[j*psp->hm_k_slices+l]);
	}
      				
      }
      
      if (psp->hm_highest_order > 3){
        /* For linear power spectrum at abs(vec(k_1) + vec(k_2)) (Needed for T_2h) */
        /* We want Pk(|vec(k_i)+vec(k_(i+1))|) --> assuming just square configurations for trispectra the absolute value is always sqrt(2) * |k|  */
        double pk_for_T_2h;
        class_call(spectra_pk_at_k_and_z(pba,ppm,psp,k*pow(2,0.5),z,&pk_for_T_2h,pk_ic),
               psp->error_message,
               psp->error_message);
	
	T_1harray[j*psp->hm_k_slices+l] = M_04[j*psp->hm_k_slices/* *psp->hm_k_slices*psp->hm_k_slices * psp->hm_k_slices + l*psp->hm_k_slices*psp->hm_k_slices*psp->hm_k_slices+l*psp->hm_k_slices*psp->hm_k_slices+l*psp->hm_k_slices*/+l] ;				
        T_2harray[j*psp->hm_k_slices+l] = 4. * ( M_11[j*psp->hm_k_slices + l] * M_13[j*psp->hm_k_slices/* *psp->hm_k_slices*psp->hm_k_slices + l*psp->hm_k_slices*psp->hm_k_slices+l*psp->hm_k_slices+l] * psp->pk_linarray[j*psp->hm_k_slices*/+l] 			
                                     +     pow(M_12[j*psp->hm_k_slices/* *psp->hm_k_slices + l*psp->hm_k_slices */ + l],2) * pk_for_T_2h);
        T_3harray[j*psp->hm_k_slices+l] = ( pow(M_11[j*psp->hm_k_slices + l],2) * M_12[j*psp->hm_k_slices/* *psp->hm_k_slices + l*psp->hm_k_slices */+ l] * 1. /* To be done */
                                       + pow(M_11[j*psp->hm_k_slices + l],2) * M_22[j*psp->hm_k_slices/* *psp->hm_k_slices + l*psp->hm_k_slices */ + l] * pow(psp->pk_linarray[j*psp->hm_k_slices+l],2));
        T_4harray[j*psp->hm_k_slices+l] = ( pow(M_11[j*psp->hm_k_slices + l],4) * psp->tk_linarray[j*psp->hm_k_slices+l]
                                        + 4. * pow(M_11[j*psp->hm_k_slices + l],3) * M_21[j*psp->hm_k_slices + l] * pow(psp->pk_linarray[j*psp->hm_k_slices+l],3));
	
      }
    }
  }
    
    
    /* Free non-psp & non output - arrays that were allocated in main part of Pk convention */
    free(multiplicityarray);
    free(biasarray);
    free(bias_pert);
    free(normmultiplicityarray);
    free(multipl_renormarray);
    if (psp->hm_highest_order > 2){
      free(muarray);
      free(Fs_2array);
      free(veccomponentarray);
      free(lengtharray);
      free(rel_lengtharray);
      free(scalar_prod);
    }
    
    free(M_00);
    free(M_10);
    free(M_20);    
    free(M_01);
    free(M_11);
    free(M_21);       
    free(M_02);
    free(M_12);
    free(M_22);
    if (psp->hm_highest_order > 2){
      free(M_03);
      free(M_30);
      free(M_31);
      free(M_13);
    }
    if (psp->hm_highest_order > 3)
      free(M_04);
      free(M_40);
          

  /* Write all necessary data in output array + allocate ln of spectra to be splined lateron */
  /* All indices are enumerated in the spectra_indices function */
  /* printf("Number_output=%f\n",(double)(psp->number_output_hm));*/
  class_alloc(psp->halomodel_pk,psp->hm_k_slices*psp->hm_z_slices*psp->number_output_hm*sizeof(double),psp->error_message);
  class_alloc(psp->ln_halomodel_Pk,psp->hm_k_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
  if (psp->hm_highest_order > 2)
    class_alloc(psp->ln_halomodel_Bk,psp->hm_k_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
  if (psp->hm_highest_order > 3)
    class_alloc(psp->ln_halomodel_Tk,psp->hm_k_slices*psp->hm_z_slices*sizeof(double),psp->error_message);
  
  for (j=0;j<psp->hm_z_slices;j++){
    for (l=0;l<psp->hm_k_slices;l++){
	psp->halomodel_pk[l*psp->hm_z_slices*psp->number_output_hm + j*psp->number_output_hm+psp->index_z_out]
	                = psp->hm_zarray[j];
	psp->halomodel_pk[l*psp->hm_z_slices*psp->number_output_hm + j*psp->number_output_hm+psp->index_k_out]
	                = pow(10,psp->log10karray[l]);
	psp->halomodel_pk[l*psp->hm_z_slices*psp->number_output_hm + j*psp->number_output_hm+psp->index_p1h_out]
	                = P_1harray[j*psp->hm_k_slices+l];
	psp->halomodel_pk[l*psp->hm_z_slices*psp->number_output_hm + j*psp->number_output_hm+psp->index_p2h_out]
	                =P_2harray[j*psp->hm_k_slices+l];
	psp->ln_halomodel_Pk[j*psp->hm_k_slices+l] = log(P_1harray[j*psp->hm_k_slices+l]
                                                   + P_2harray[j*psp->hm_k_slices+l]);
	/*printf("ln_pk_hm(%f,%f)=%f\n",psp->log10karray[l],psp->hm_zarray[j],psp->ln_halomodel_Pk[j*psp->hm_k_slices+l]);*/
	
	if (psp->hm_highest_order > 2){
	psp->halomodel_pk[l*psp->hm_z_slices*psp->number_output_hm + j*psp->number_output_hm+psp->index_b1h_out]
	                = B_1harray[j*psp->hm_k_slices+l];
	psp->halomodel_pk[l*psp->hm_z_slices*psp->number_output_hm + j*psp->number_output_hm+psp->index_b2h_out]
	                = B_2harray[j*psp->hm_k_slices+l];
	psp->halomodel_pk[l*psp->hm_z_slices*psp->number_output_hm + j*psp->number_output_hm+psp->index_b3h_out]
	                = B_3harray[j*psp->hm_k_slices+l];
	psp->ln_halomodel_Bk[j*psp->hm_k_slices+l] = log(B_1harray[j*psp->hm_k_slices+l]
                                                   + B_2harray[j*psp->hm_k_slices+l]+B_3harray[j*psp->hm_k_slices+l]);		
	
	/*printf("ln_bk_hm(%f,%f)=%f\n",psp->log10karray[l],psp->hm_zarray[j],psp->ln_halomodel_Bk[j*psp->hm_k_slices+l]);*/
	}
	if (psp->hm_highest_order > 3){
	psp->halomodel_pk[l*psp->hm_z_slices*psp->number_output_hm + j*psp->number_output_hm+psp->index_t1h_out]
	                = T_1harray[j*psp->hm_k_slices+l];
	psp->halomodel_pk[l*psp->hm_z_slices*psp->number_output_hm + j*psp->number_output_hm+psp->index_t2h_out]
	                = T_2harray[j*psp->hm_k_slices+l];
	psp->halomodel_pk[l*psp->hm_z_slices*psp->number_output_hm + j*psp->number_output_hm+psp->index_t3h_out]
	                = T_3harray[j*psp->hm_k_slices+l];
	psp->halomodel_pk[l*psp->hm_z_slices*psp->number_output_hm + j*psp->number_output_hm+psp->index_t4h_out]
	                = T_4harray[j*psp->hm_k_slices+l];
	psp->ln_halomodel_Tk[j*psp->hm_k_slices+l] = log(T_1harray[j*psp->hm_k_slices+l]
                                                   + T_2harray[j*psp->hm_k_slices+l]  + T_3harray[j*psp->hm_k_slices+l]  + T_4harray[j*psp->hm_k_slices+l]);
	}
      }
    }
    
    
  /* Compute array of second derivatives needed lateron for spline interpolation */  
  int index_md;
  class_alloc(psp->ddln_halomodel_Pk,sizeof(double)*psp->hm_z_slices*psp->hm_k_slices,psp->error_message);
  class_call(array_spline_table_lines(psp->ln_tau_hm,
                                        psp->hm_z_slices,
                                        psp->ln_halomodel_Pk,
                                        psp->hm_k_slices,
                                        psp->ddln_halomodel_Pk,
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);
  if (psp->hm_highest_order > 2){
    class_alloc(psp->ddln_halomodel_Bk,sizeof(double)*psp->hm_z_slices*psp->hm_k_slices,psp->error_message);
    class_call(array_spline_table_lines(psp->ln_tau_hm,
                                        psp->hm_z_slices,
                                        psp->ln_halomodel_Bk,
                                        psp->hm_k_slices,
                                        psp->ddln_halomodel_Bk,
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);
    
    class_alloc(psp->ddln_Q123fixarray,sizeof(double)*psp->hm_z_slices*psp->hm_k_slices,psp->error_message);
    class_call(array_spline_table_lines(psp->ln_tau_hm,
                                        psp->hm_z_slices,
                                        psp->Q_123fixarray,
                                        psp->hm_k_slices,
                                        psp->ddln_Q123fixarray,
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);
    
    class_alloc(psp->ddln_Q123fixlinarray,sizeof(double)*psp->hm_z_slices*psp->hm_k_slices,psp->error_message);
    class_call(array_spline_table_lines(psp->ln_tau_hm,
                                        psp->hm_z_slices,
                                        psp->Q_123fixlinarray,
                                        psp->hm_k_slices,
                                        psp->ddln_Q123fixlinarray,
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);
    
    class_alloc(psp->ddln_halomodel_bk_lin,sizeof(double)*psp->hm_z_slices*psp->hm_k_slices,psp->error_message);
    class_call(array_spline_table_lines(psp->ln_tau_hm,
                                        psp->hm_z_slices,
                                        psp->bk_fixlinarray,
                                        psp->hm_k_slices,
                                        psp->ddln_halomodel_Bk,
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);
  }

  if (psp->hm_highest_order > 3){
    class_alloc(psp->ddln_halomodel_Tk,sizeof(double)*psp->hm_z_slices*psp->hm_k_slices,psp->error_message);
    class_call(array_spline_table_lines(psp->ln_tau_hm,
                                        psp->hm_z_slices,
                                        psp->ln_halomodel_Tk,
                                        psp->hm_k_slices,
                                        psp->ddln_halomodel_Tk,
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);
  }
  
  /*class_alloc(psp->ddlnsigmatable,sizeof(double)*psp->sigma_z_slices*psp->sigma_R_slices,psp->error_message);
  class_call(array_spline_table_lines(psp->sigma_ln_tauarray,
                                        psp->sigma_z_slices,
                                        psp->lnsigmatable,
                                        psp->sigma_R_slices,
                                        psp->ddlnsigmatable,
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);*/
  
  /*for(j=0;j<psp->hm_z_slices;j++){
      for(l=0;l<psp->hm_k_slices;l++){
        printf("l=%f,j=%f\n",(double)(l),(double)(j));
        printf("ddlnBk[k=%f,z=%f] = %f\n",pow(10,psp->log10karray[l]),psp->hm_zarray[j],psp->ddln_halomodel_Bk[j*psp->hm_k_slices+l]);	
      }
    }*/

  /* Free all 'short' arrays */
  free(aarray);
  free(Deltaarray);
  free(rhobar_bgarray);
  free(Omega_marray);
  free(rhobar_marray);
  free(Harray);
  free(pvecback_sp_long);
  free(rhobar_cut_halo);
  free(D_plus_fitarray);
  free(D_plus_classarray);
  
  /* Free all 'long' arrays */
  free(log10marray);  
  free(m_to_Rarray);
  free(nuarray);
  free(nutolog10marray);
  free(nu1tolog10marray);
  free(carray);
  free(r_sarray);
  free(rhoarray);
  free(Siarray);
  free(Ciarray);
  free(ukmarray);
  free(P_1harray);
  free(P_2harray);
  if (psp->hm_highest_order > 2){
    free(B_1harray);
    free(B_2harray);
    free(B_3harray);
  }
  if (psp->hm_highest_order > 3){
    free(T_1harray);
    free(T_2harray);
    free(T_3harray);
    free(T_4harray);
  }
  
  return _SUCCESS_;
  
}




int spectra_pk_hm_at_z(
                    struct background * pba,
                    struct spectra * psp,
                    double z,
                    double * output_tot /* array with argument output_tot[index_k] (must be already allocated) */
                    ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  int index_k;
  double tau1,ln_tau1;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau1),
             pba->error_message,
             psp->error_message);

  class_test(tau1 <= 0.,
             psp->error_message,
             "negative or null value of conformal time: cannot interpolate");

  ln_tau1 = log(tau1);
  
  
  /** - second step: for both modes (linear or logarithmic), store the spectrum in logarithmic format in the output array(s) */

  /**   (a.) if only values at tau=tau_today are stored and we want P(k,z=0), no need to interpolate */

  if (psp->hm_z_slices == 1) {
    class_test(z != 0.,
               psp->error_message,
               "asked z=%e but only P(k,z=0) has been tabulated",z);
    for (index_k=0; index_k<psp->hm_k_slices; index_k++)
      	output_tot[index_k] = psp->ln_halomodel_Pk[index_k];
           
  }

  /**   (b.) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else { 
    
    int j;
    for (j=0;j<psp->hm_z_slices;j++)
      /*printf("tau(%f) = %f \n",psp->hm_zarray[j],psp->ln_tau_hm[j]);*/
    
    class_call(array_interpolate_spline(psp->ln_tau_hm,
                                          psp->hm_z_slices,
                                          psp->ln_halomodel_Pk,
                                          psp->ddln_halomodel_Pk,
                                          psp->hm_k_slices,
                                          ln_tau1,
                                          &last_index,
                                          output_tot,
                                          psp->hm_k_slices,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);   
    
  }
  
     
  return _SUCCESS_;

}


int spectra_bk_hm_at_z(
                    struct background * pba,
                    struct spectra * psp,
                    double z,
                    double * output_tot /* array with argument output_tot[index_k] (must be already allocated) */
                    ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  int index_k;
  double tau1,ln_tau1;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau1),
             pba->error_message,
             psp->error_message);

  class_test(tau1 <= 0.,
             psp->error_message,
             "negative or null value of conformal time: cannot interpolate");

  ln_tau1 = log(tau1);
  
  
  /** - second step: for both modes (linear or logarithmic), store the spectrum in logarithmic format in the output array(s) */

  /**   (a.) if only values at tau=tau_today are stored and we want P(k,z=0), no need to interpolate */

  if (psp->hm_z_slices == 1) {
    class_test(z != 0.,
               psp->error_message,
               "asked z=%e but only P(k,z=0) has been tabulated",z);
    for (index_k=0; index_k<psp->hm_k_slices; index_k++)
      	output_tot[index_k] = psp->ln_halomodel_Bk[index_k];
           
  }

  /**   (b.) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else { 
    
    int j;
    for (j=0;j<psp->hm_z_slices;j++)
      /*printf("tau(%f) = %f \n",psp->hm_zarray[j],psp->ln_tau_hm[j]);*/
    
    class_call(array_interpolate_spline(psp->ln_tau_hm,
                                          psp->hm_z_slices,
                                          psp->ln_halomodel_Bk,
                                          psp->ddln_halomodel_Bk,
                                          psp->hm_k_slices,
                                          ln_tau1,
                                          &last_index,
                                          output_tot,
                                          psp->hm_k_slices,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);   
    
    
  }
     
  return _SUCCESS_;

}

int spectra_bk_lin_at_z(
                    struct background * pba,
                    struct spectra * psp,
                    double z,
                    double * output_tot /* array with argument output_tot[index_k] (must be already allocated) */
                    ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  int index_k;
  double tau1,ln_tau1;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau1),
             pba->error_message,
             psp->error_message);

  class_test(tau1 <= 0.,
             psp->error_message,
             "negative or null value of conformal time: cannot interpolate");

  ln_tau1 = log(tau1);
  
  
  /** - second step: for both modes (linear or logarithmic), store the spectrum in logarithmic format in the output array(s) */

  /**   (a.) if only values at tau=tau_today are stored and we want P(k,z=0), no need to interpolate */

  if (psp->hm_z_slices == 1) {
    class_test(z != 0.,
               psp->error_message,
               "asked z=%e but only P(k,z=0) has been tabulated",z);
    for (index_k=0; index_k<psp->hm_k_slices; index_k++)
      	output_tot[index_k] = psp->bk_linarray[index_k];
           
  }

  /**   (b.) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else { 
    
    int j;
    for (j=0;j<psp->hm_z_slices;j++)
      /*printf("tau(%f) = %f \n",psp->hm_zarray[j],psp->ln_tau_hm[j]);*/
    
    class_call(array_interpolate_spline(psp->ln_tau_hm,
                                          psp->hm_z_slices,
                                          psp->bk_fixlinarray,
                                          psp->ddln_halomodel_bk_lin,
                                          psp->hm_k_slices,
                                          ln_tau1,
                                          &last_index,
                                          output_tot,
                                          psp->hm_k_slices,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);   
    
    
  }
     
  return _SUCCESS_;

}

int spectra_Q123lin_at_z(
                    struct background * pba,
                    struct spectra * psp,
                    double z,
                    double * output_tot /* array with argument output_tot[index_k] (must be already allocated) */
                    ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  int index_k;
  double tau1,ln_tau1;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau1),
             pba->error_message,
             psp->error_message);

  class_test(tau1 <= 0.,
             psp->error_message,
             "negative or null value of conformal time: cannot interpolate");

  ln_tau1 = log(tau1);
  
  
  /** - second step: for both modes (linear or logarithmic), store the spectrum in logarithmic format in the output array(s) */

  /**   (a.) if only values at tau=tau_today are stored and we want P(k,z=0), no need to interpolate */

  if (psp->hm_z_slices == 1) {
    class_test(z != 0.,
               psp->error_message,
               "asked z=%e but only P(k,z=0) has been tabulated",z);
    for (index_k=0; index_k<psp->hm_k_slices; index_k++)
      	output_tot[index_k] = psp->Q_123fixlinarray[index_k];
           
  }

  /**   (b.) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else { 
    
    int j;
    for (j=0;j<psp->hm_z_slices;j++)
      /*printf("tau(%f) = %f \n",psp->hm_zarray[j],psp->ln_tau_hm[j]);*/
    
    class_call(array_interpolate_spline(psp->ln_tau_hm,
                                          psp->hm_z_slices,
                                          psp->Q_123fixlinarray,
                                          psp->ddln_Q123fixlinarray,
                                          psp->hm_k_slices,
                                          ln_tau1,
                                          &last_index,
                                          output_tot,
                                          psp->hm_k_slices,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);   
    
    
  }
     
  return _SUCCESS_;

}


int spectra_Q123_at_z(
                    struct background * pba,
                    struct spectra * psp,
                    double z,
                    double * output_tot /* array with argument output_tot[index_k] (must be already allocated) */
                    ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  int index_k;
  double tau1,ln_tau1;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau1),
             pba->error_message,
             psp->error_message);

  class_test(tau1 <= 0.,
             psp->error_message,
             "negative or null value of conformal time: cannot interpolate");

  ln_tau1 = log(tau1);
  
  
  /** - second step: for both modes (linear or logarithmic), store the spectrum in logarithmic format in the output array(s) */

  /**   (a.) if only values at tau=tau_today are stored and we want P(k,z=0), no need to interpolate */

  if (psp->hm_z_slices == 1) {
    class_test(z != 0.,
               psp->error_message,
               "asked z=%e but only P(k,z=0) has been tabulated",z);
    for (index_k=0; index_k<psp->hm_k_slices; index_k++)
      	output_tot[index_k] = psp->Q_123fixarray[index_k];
           
  }

  /**   (b.) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else { 
    
    int j;
    for (j=0;j<psp->hm_z_slices;j++)
      /*printf("tau(%f) = %f \n",psp->hm_zarray[j],psp->ln_tau_hm[j]);*/
    
    class_call(array_interpolate_spline(psp->ln_tau_hm,
                                          psp->hm_z_slices,
                                          psp->Q_123fixarray,
                                          psp->ddln_Q123fixarray,
                                          psp->hm_k_slices,
                                          ln_tau1,
                                          &last_index,
                                          output_tot,
                                          psp->hm_k_slices,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);   
    
    
  }
     
  return _SUCCESS_;

}

int spectra_tk_hm_at_z(
                    struct background * pba,
                    struct spectra * psp,
                    double z,
                    double * output_tot /* array with argument output_tot[index_k] (must be already allocated) */
                    ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  int index_k;
  double tau1,ln_tau1;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau1),
             pba->error_message,
             psp->error_message);

  class_test(tau1 <= 0.,
             psp->error_message,
             "negative or null value of conformal time: cannot interpolate");

  ln_tau1 = log(tau1);
  
  
  /** - second step: for both modes (linear or logarithmic), store the spectrum in logarithmic format in the output array(s) */

  /**   (a.) if only values at tau=tau_today are stored and we want P(k,z=0), no need to interpolate */

  if (psp->hm_z_slices == 1) {
    class_test(z != 0.,
               psp->error_message,
               "asked z=%e but only P(k,z=0) has been tabulated",z);
    for (index_k=0; index_k<psp->hm_k_slices; index_k++)
      	output_tot[index_k] = psp->ln_halomodel_Tk[index_k];
           
  }

  /**   (b.) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else { 
    
    int j;
    for (j=0;j<psp->hm_z_slices;j++)
      /*printf("tau(%f) = %f \n",psp->hm_zarray[j],psp->ln_tau_hm[j]);*/
    
    class_call(array_interpolate_spline(psp->ln_tau_hm,
                                          psp->hm_z_slices,
                                          psp->ln_halomodel_Tk,
                                          psp->ddln_halomodel_Tk,
                                          psp->hm_k_slices,
                                          ln_tau1,
                                          &last_index,
                                          output_tot,
                                          psp->hm_k_slices,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);   
    
  }
  
     
  return _SUCCESS_;

}


int spectra_pk_hm_at_k_and_z(
                          struct background * pba,
                          struct primordial * ppm,
                          struct spectra * psp,
                          double k,
                          double z,
                          double * pk_hm_tot /* pointer to a single number (must be already allocated) */
                          ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_k;
  int last_index;
  int index_ic1,index_ic2,index_ic1_ic2;

  double * spectrum_hm_at_z = NULL;
  double * spline;
  double * pk_primordial_k = NULL;
  double kmin;
  double * pk_primordial_kmin = NULL;

  index_md = psp->index_md_scalars;

  /** - first step: check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_pk_at_z()) */

  class_test((k < 0.) || (k > pow(10,(psp->log10karray[psp->hm_k_slices-1]))),
             psp->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,pow(10,(psp->log10karray[psp->hm_k_slices-1])));

  /** - deal with case 0 <= k < kmin */

  if (k < pow(10,(psp->log10karray[0]))) {

    /**   (a.) subcase k=0: then P(k)=0 */

    if (k == 0.) {
        *pk_hm_tot=0.;
    }
 

    /**    (b.) subcase 0<k<kmin: in this case we know that on super-Hubble scales:
     *          P(k) = [some number] * k  * P_primordial(k)
     *          so
     *          P(k) = P(kmin) * (k P_primordial(k)) / (kmin P_primordial(kmin))
     *          (note that the result is accurate only if kmin is such that [a0 kmin] << H0)
     */

    else {

      /* compute P(k,z) which contains P(kmin,z)*/
      class_alloc(spectrum_hm_at_z,
                  psp->hm_k_slices*sizeof(double),
                  psp->error_message);
      class_call(spectra_pk_hm_at_z(pba,
                                 psp,
                                 z,
                                 spectrum_hm_at_z),
                 psp->error_message,
                 psp->error_message);

      /* compute P_primordial(k) */
      class_alloc(pk_primordial_k,
                  sizeof(double)*psp->ic_ic_size[index_md],
                  psp->error_message);
      class_call(primordial_spectrum_at_k(ppm,
                                          index_md,
                                          linear,
                                          k,
                                          pk_primordial_k),
                 ppm->error_message,psp->error_message);

      /* compute P_primordial(kmin) */
      kmin = pow(10,(psp->log10karray[psp->hm_k_slices-1]));
      class_alloc(pk_primordial_kmin,
                  sizeof(double)*psp->ic_ic_size[index_md],
                  psp->error_message);
      class_call(primordial_spectrum_at_k(ppm,
                                          index_md,
                                          linear,
                                          kmin,
                                          pk_primordial_kmin),
                 ppm->error_message,
                 psp->error_message);

      /* apply above analytic approximation for P(k) */
      index_k=0;
      index_ic1_ic2 = 0;
      *pk_hm_tot = spectrum_hm_at_z[index_k]
        *k*pk_primordial_k[index_ic1_ic2]
        /kmin/pk_primordial_kmin[index_ic1_ic2];
      
  
      free(spectrum_hm_at_z);
      free(pk_primordial_k);
      free(pk_primordial_kmin);

    }
  }

  /** - deal with case kmin <= k <= kmax */

  else {

    /* compute P(k,z) (in logarithmic format for more accurate interpolation) */
    class_alloc(spectrum_hm_at_z,
                psp->hm_k_slices*sizeof(double),
                psp->error_message);
    class_call(spectra_pk_hm_at_z(pba,
                               psp,
                               z,
                               spectrum_hm_at_z),
               psp->error_message,
               psp->error_message);

    /* get its second derivatives with spline, then interpolate, then convert to linear format */

    class_alloc(spline,
                sizeof(double)*psp->hm_k_slices,
                psp->error_message);

    class_call(array_spline_table_lines(psp->log10karray,
                                          psp->hm_k_slices,
                                          spectrum_hm_at_z,
                                          1,
                                          spline,
                                          _SPLINE_NATURAL_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

    class_call(array_interpolate_spline(psp->log10karray,
                                          psp->hm_k_slices,
                                          spectrum_hm_at_z,
                                          spline,
                                          1,
                                          log10(k),
                                          &last_index,
                                          pk_hm_tot,
                                          1,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

    *pk_hm_tot = exp(*pk_hm_tot);
    }

    free(spectrum_hm_at_z);
    free(spline);


  return _SUCCESS_;

}


int spectra_bk_hm_at_k_and_z(
                          struct background * pba,
                          struct primordial * ppm,
                          struct spectra * psp,
                          double k,
                          double z,
                          double * bk_hm_tot /* pointer to a single number (must be already allocated) */
                          ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  double * spectrum_hm_at_z = NULL;
  double * spline;
 
  /** - first step: check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_pk_at_z()) */

  class_test((k < 0.) || (k > pow(10,(psp->log10karray[psp->hm_k_slices-1]))),
             psp->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,pow(10,(psp->log10karray[psp->hm_k_slices-1])));

  /** - deal with case 0 <= k < kmin */

  if (k < pow(10,(psp->log10karray[0]))) {
      printf("k = %f < %f = Minimum of interpolation table for Halomodel spectra. Bispectrum set to zero for k<%f",k,pow(10,(psp->log10karray[0])),pow(10,(psp->log10karray[0])));
      *bk_hm_tot=0.;
  }

  else {

    /* compute P(k,z) (in logarithmic format for more accurate interpolation) */
    class_alloc(spectrum_hm_at_z,
                psp->hm_k_slices*sizeof(double),
                psp->error_message);
    class_call(spectra_bk_hm_at_z(pba,
                               psp,
                               z,
                               spectrum_hm_at_z),
               psp->error_message,
               psp->error_message);

    /* get its second derivatives with spline, then interpolate, then convert to linear format */

    class_alloc(spline,
                sizeof(double)*psp->hm_k_slices,
                psp->error_message);

    class_call(array_spline_table_lines(psp->log10karray,
                                          psp->hm_k_slices,
                                          spectrum_hm_at_z,
                                          1,
                                          spline,
                                          _SPLINE_NATURAL_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

    class_call(array_interpolate_spline(psp->log10karray,
                                          psp->hm_k_slices,
                                          spectrum_hm_at_z,
                                          spline,
                                          1,
                                          log10(k),
                                          &last_index,
                                          bk_hm_tot,
                                          1,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

    *bk_hm_tot = exp(*bk_hm_tot);

    }

    free(spectrum_hm_at_z);
    free(spline);


  return _SUCCESS_;

}

int spectra_bk_lin_at_k_and_z(
                          struct background * pba,
                          struct primordial * ppm,
                          struct spectra * psp,
                          double k,
                          double z,
                          double * bk_lin_tot /* pointer to a single number (must be already allocated) */
                          ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  double * spectrum_hm_at_z = NULL;
  double * spline;
 
  /** - first step: check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_pk_at_z()) */

  class_test((k < 0.) || (k > pow(10,(psp->log10karray[psp->hm_k_slices-1]))),
             psp->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,pow(10,(psp->log10karray[psp->hm_k_slices-1])));

  /** - deal with case 0 <= k < kmin */

  if (k < pow(10,(psp->log10karray[0]))) {
      printf("k = %f < %f = Minimum of interpolation table for Halomodel spectra. Bispectrum set to zero for k<%f",k,pow(10,(psp->log10karray[0])),pow(10,(psp->log10karray[0])));
      *bk_lin_tot=0.;
  }

  else {

    /* compute P(k,z) (in logarithmic format for more accurate interpolation) */
    class_alloc(spectrum_hm_at_z,
                psp->hm_k_slices*sizeof(double),
                psp->error_message);
    class_call(spectra_bk_lin_at_z(pba,
                               psp,
                               z,
                               spectrum_hm_at_z),
               psp->error_message,
               psp->error_message);

    /* get its second derivatives with spline, then interpolate, then convert to linear format */

    class_alloc(spline,
                sizeof(double)*psp->hm_k_slices,
                psp->error_message);

    class_call(array_spline_table_lines(psp->log10karray,
                                          psp->hm_k_slices,
                                          spectrum_hm_at_z,
                                          1,
                                          spline,
                                          _SPLINE_NATURAL_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

    class_call(array_interpolate_spline(psp->log10karray,
                                          psp->hm_k_slices,
                                          spectrum_hm_at_z,
                                          spline,
                                          1,
                                          log10(k),
                                          &last_index,
                                          bk_lin_tot,
                                          1,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

    }

    free(spectrum_hm_at_z);
    free(spline);


  return _SUCCESS_;

}


int spectra_Q123lin_at_k_and_z(
                          struct background * pba,
                          struct primordial * ppm,
                          struct spectra * psp,
                          double k,
                          double z,
                          double * Q123lin_tot /* pointer to a single number (must be already allocated) */
                          ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  double * spectrum_hm_at_z = NULL;
  double * spline;
 
  /** - first step: check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_pk_at_z()) */

  class_test((k < 0.) || (k > pow(10,(psp->log10karray[psp->hm_k_slices-1]))),
             psp->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,pow(10,(psp->log10karray[psp->hm_k_slices-1])));

  /** - deal with case 0 <= k < kmin */

  if (k < pow(10,(psp->log10karray[0]))) {
      printf("k = %f < %f = Minimum of interpolation table for Halomodel spectra. Bispectrum set to zero for k<%f",k,pow(10,(psp->log10karray[0])),pow(10,(psp->log10karray[0])));
      *Q123lin_tot=0.;
  }

  else {

    /* compute P(k,z) (in logarithmic format for more accurate interpolation) */
    class_alloc(spectrum_hm_at_z,
                psp->hm_k_slices*sizeof(double),
                psp->error_message);
    class_call(spectra_Q123lin_at_z(pba,
                               psp,
                               z,
                               spectrum_hm_at_z),
               psp->error_message,
               psp->error_message);

    /* get its second derivatives with spline, then interpolate, then convert to linear format */

    class_alloc(spline,
                sizeof(double)*psp->hm_k_slices,
                psp->error_message);

    class_call(array_spline_table_lines(psp->log10karray,
                                          psp->hm_k_slices,
                                          spectrum_hm_at_z,
                                          1,
                                          spline,
                                          _SPLINE_NATURAL_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

    class_call(array_interpolate_spline(psp->log10karray,
                                          psp->hm_k_slices,
                                          spectrum_hm_at_z,
                                          spline,
                                          1,
                                          log10(k),
                                          &last_index,
                                          Q123lin_tot,
                                          1,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);
    }

    free(spectrum_hm_at_z);
    free(spline);


  return _SUCCESS_;

}

int spectra_Q123_at_k_and_z(
                          struct background * pba,
                          struct primordial * ppm,
                          struct spectra * psp,
                          double k,
                          double z,
                          double * Q123_tot /* pointer to a single number (must be already allocated) */
                          ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  double * spectrum_hm_at_z = NULL;
  double * spline;
 
  /** - first step: check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_pk_at_z()) */

  class_test((k < 0.) || (k > pow(10,(psp->log10karray[psp->hm_k_slices-1]))),
             psp->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,pow(10,(psp->log10karray[psp->hm_k_slices-1])));

  /** - deal with case 0 <= k < kmin */

  if (k < pow(10,(psp->log10karray[0]))) {
      printf("k = %f < %f = Minimum of interpolation table for Halomodel spectra. Bispectrum set to zero for k<%f",k,pow(10,(psp->log10karray[0])),pow(10,(psp->log10karray[0])));
      *Q123_tot=0.;
  }

  else {

    /* compute P(k,z) (in logarithmic format for more accurate interpolation) */
    class_alloc(spectrum_hm_at_z,
                psp->hm_k_slices*sizeof(double),
                psp->error_message);
    class_call(spectra_Q123_at_z(pba,
                               psp,
                               z,
                               spectrum_hm_at_z),
               psp->error_message,
               psp->error_message);

    /* get its second derivatives with spline, then interpolate, then convert to linear format */

    class_alloc(spline,
                sizeof(double)*psp->hm_k_slices,
                psp->error_message);

    class_call(array_spline_table_lines(psp->log10karray,
                                          psp->hm_k_slices,
                                          spectrum_hm_at_z,
                                          1,
                                          spline,
                                          _SPLINE_NATURAL_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

    class_call(array_interpolate_spline(psp->log10karray,
                                          psp->hm_k_slices,
                                          spectrum_hm_at_z,
                                          spline,
                                          1,
                                          log10(k),
                                          &last_index,
                                          Q123_tot,
                                          1,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);
    }

    free(spectrum_hm_at_z);
    free(spline);


  return _SUCCESS_;

}

int spectra_tk_hm_at_k_and_z(
                          struct background * pba,
                          struct primordial * ppm,
                          struct spectra * psp,
                          double k,
                          double z,
                          double * tk_hm_tot /* pointer to a single number (must be already allocated) */
                          ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  double * spectrum_hm_at_z = NULL;
  double * spline;


  /** - first step: check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_pk_at_z()) */

  class_test((k < 0.) || (k > pow(10,(psp->log10karray[psp->hm_k_slices-1]))),
             psp->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,pow(10,(psp->log10karray[psp->hm_k_slices-1])));

  /** - deal with case 0 <= k < kmin */

  if (k < pow(10,(psp->log10karray[0]))) {
      printf("k = %f < %f = Minimum of interpolation table for Halomodel spectra. Trispectrum set to zero for k<%f",k,pow(10,(psp->log10karray[0])),pow(10,(psp->log10karray[0])));
      *tk_hm_tot=0.;
  }

  else {

    /* compute P(k,z) (in logarithmic format for more accurate interpolation) */
    class_alloc(spectrum_hm_at_z,
                psp->hm_k_slices*sizeof(double),
                psp->error_message);
    class_call(spectra_tk_hm_at_z(pba,
                               psp,
                               z,
                               spectrum_hm_at_z),
               psp->error_message,
               psp->error_message);

    /* get its second derivatives with spline, then interpolate, then convert to linear format */

    class_alloc(spline,
                sizeof(double)*psp->hm_k_slices,
                psp->error_message);

    class_call(array_spline_table_lines(psp->log10karray,
                                          psp->hm_k_slices,
                                          spectrum_hm_at_z,
                                          1,
                                          spline,
                                          _SPLINE_NATURAL_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

    class_call(array_interpolate_spline(psp->log10karray,
                                          psp->hm_k_slices,
                                          spectrum_hm_at_z,
                                          spline,
                                          1,
                                          log10(k),
                                          &last_index,
                                          tk_hm_tot,
                                          1,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

    *tk_hm_tot = exp(*tk_hm_tot);

    }

    free(spectrum_hm_at_z);
    free(spline);


  return _SUCCESS_;	

}


/**
 * This routine frees all the memory space allocated by spectra_init().
 *
 * To be called at the end of each run, only when no further calls to
 * spectra_cls_at_l(), spectra_pk_at_z(), spectra_pk_at_k_and_z() are needed.
 *
 * @param psp Input: pointer to spectra structure (which fields must be freed)
 * @return the error status
 */

int spectra_free(
                 struct spectra * psp
                 ) {

  int index_md;

  if (psp->md_size > 0) {

    if (psp->ct_size > 0) {

      for (index_md = 0; index_md < psp->md_size; index_md++) {
        free(psp->l_max_ct[index_md]);
        free(psp->cl[index_md]);
        free(psp->ddcl[index_md]);
      }
      free(psp->l);
      free(psp->l_size);
      free(psp->l_max_ct);
      free(psp->l_max);
      free(psp->cl);
      free(psp->ddcl);
    }

    if (psp->ln_k_size > 0) {

      free(psp->ln_tau);
      free(psp->ln_k);

      if (psp->ln_pk != NULL) {

        free(psp->ln_pk);

        if (psp->ln_tau_size > 1) {
          free(psp->ddln_pk);
        }

        if (psp->ln_pk_nl != NULL) {

          free(psp->ln_pk_nl);

          if (psp->ln_tau_size > 1) {
            free(psp->ddln_pk_nl);
          }
        }
      }

      if (psp->matter_transfer != NULL) {

        free(psp->matter_transfer);
        if (psp->ln_tau_size > 1) {
          free(psp->ddmatter_transfer);
        }
      }
    }
  }

  for (index_md=0; index_md < psp->md_size; index_md++)
    free(psp->is_non_zero[index_md]);
  free(psp->is_non_zero);
  free(psp->ic_size);
  free(psp->ic_ic_size);
  
  /* Free Halomodel quantities */
  
  /* From halomodel table */
  if (psp->log10karray != NULL){
  free(psp->log10karray);
  free(psp->hm_zarray);  
  free(psp->ln_tau_hm);  
  free(psp->pk_linarray);
  if (psp->hm_highest_order > 2)
    free(psp->bk_linarray);
    free(psp->bk_fixlinarray);
  if (psp->hm_highest_order > 3)
    free(psp->tk_linarray);  
  free(psp->halomodel_pk);
  free(psp->ln_halomodel_Pk);
  if (psp->hm_highest_order > 2)
    free(psp->ln_halomodel_Bk);
  if (psp->hm_highest_order > 3)
    free(psp->ln_halomodel_Tk); 
  free(psp->ddln_halomodel_Pk);
  if (psp->hm_highest_order > 2){
    free(psp->ddln_halomodel_Bk);
    free(psp->ddln_halomodel_bk_lin);
    free(psp->Q_123array);
    free(psp->Q_123fixarray);
    free(psp->ddln_Q123fixarray);   
  }
  if (psp->hm_highest_order > 3)
    free(psp->ddln_halomodel_Tk);
  
  /* From sigmatable */
  free(psp->sigma_zarray);
  free(psp->sigma_ln_tauarray);
  free(psp->sigma_Rarray);
  free(psp->lnsigmatable);
  free(psp->ddlnsigmatable);
  }

  return _SUCCESS_;

}

/**
 * This routine defines indices and allocates tables in the spectra structure
 *
 * @param ppt  Input : pointer to perturbation structure
 * @param ptr  Input : pointer to transfers structure
 * @param ppm  Input : pointer to primordial structure
 * @param psp  Input/output: pointer to spectra structure
 * @return the error status
 */

int spectra_indices(
                    struct background * pba,
                    struct perturbs * ppt,
                    struct transfers * ptr,
                    struct primordial * ppm,
                    struct spectra * psp
                    ){

  int index_ct;
  int index_md;
  int index_ic1_ic2;
  int index_tr;

  psp->md_size = ppt->md_size;
  if (ppt->has_scalars == _TRUE_)
    psp->index_md_scalars = ppt->index_md_scalars;

  class_alloc(psp->ic_size,
              sizeof(int)*psp->md_size,
              psp->error_message);

  class_alloc(psp->ic_ic_size,
              sizeof(int)*psp->md_size,
              psp->error_message);

  class_alloc(psp->is_non_zero,
              sizeof(short *)*psp->md_size,
              psp->error_message);

  for (index_md=0; index_md < psp->md_size; index_md++) {
    psp->ic_size[index_md] = ppm->ic_size[index_md];
    psp->ic_ic_size[index_md] = ppm->ic_ic_size[index_md];
    class_alloc(psp->is_non_zero[index_md],
                sizeof(short)*psp->ic_ic_size[index_md],
                psp->error_message);
    for (index_ic1_ic2=0; index_ic1_ic2 < psp->ic_ic_size[index_md]; index_ic1_ic2++)
      psp->is_non_zero[index_md][index_ic1_ic2] = ppm->is_non_zero[index_md][index_ic1_ic2];
  }

  if (ppt->has_cls == _TRUE_) {

    /* types of C_l's relevant for both scalars and tensors: TT, EE, TE */

    index_ct=0;

    if (ppt->has_cl_cmb_temperature == _TRUE_) {
      psp->has_tt = _TRUE_;
      psp->index_ct_tt=index_ct;
      index_ct++;
    }
    else {
      psp->has_tt = _FALSE_;
    }

    if (ppt->has_cl_cmb_polarization == _TRUE_) {
      psp->has_ee = _TRUE_;
      psp->index_ct_ee=index_ct;
      index_ct++;
    }
    else {
      psp->has_ee = _FALSE_;
    }

    if ((ppt->has_cl_cmb_temperature == _TRUE_) &&
        (ppt->has_cl_cmb_polarization == _TRUE_)) {
      psp->has_te = _TRUE_;
      psp->index_ct_te=index_ct;
      index_ct++;
    }
    else {
      psp->has_te = _FALSE_;
    }

    if (ppt->has_cl_cmb_polarization == _TRUE_) {
      psp->has_bb = _TRUE_;
      psp->index_ct_bb=index_ct;
      index_ct++;
    }
    else {
      psp->has_bb = _FALSE_;
    }

    /* types of C_l's relevant only for scalars: phi-phi, T-phi, E-phi, d-d, T-d */

    if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_pp = _TRUE_;
      psp->index_ct_pp=index_ct;
      index_ct++;
    }
    else {
      psp->has_pp = _FALSE_;
    }

    if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_tp = _TRUE_;
      psp->index_ct_tp=index_ct;
      index_ct++;
    }
    else {
      psp->has_tp = _FALSE_;
    }

    psp->ct_size = index_ct;

    if ((ppt->has_cl_cmb_polarization == _TRUE_) && (ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_ep = _TRUE_;
      psp->index_ct_ep=index_ct;
      index_ct++;
    }
    else {
      psp->has_ep = _FALSE_;
    }

    if ((ppt->has_scalars == _TRUE_) &&
        ((ppt->has_cl_number_count == _TRUE_) || (ppt->has_cl_lensing_potential == _TRUE_)))
      psp->d_size=ppt->selection_num;
    else
      psp->d_size=0;

    if ((ppt->has_cl_number_count == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_dd = _TRUE_;
      psp->index_ct_dd=index_ct;
      index_ct+=(psp->d_size*(psp->d_size+1)-(psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag))/2;
    }
    else {
      psp->has_dd = _FALSE_;
    }

    /* the computation of C_l^Td would require a very good sampling of
       transfer functions over a wide range, and a huge computation
       time. In the current version, we prefer to switch it off, rather
       than either slowing down the code considerably, or producing
       very inaccurate spectra.

       if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_number_count == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
       psp->has_td = _TRUE_;
       psp->index_ct_td=index_ct;
       index_ct+=psp->d_size;
       }
       else {
       psp->has_td = _FALSE_;
       }
    */
    psp->has_td = _FALSE_;

    if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (ppt->has_cl_number_count == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_pd = _TRUE_;
      psp->index_ct_pd=index_ct;
      index_ct+=psp->d_size;
    }
    else {
      psp->has_pd = _FALSE_;
    }

    psp->has_td = _FALSE_;

    if ((ppt->has_cl_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_ll = _TRUE_;
      psp->index_ct_ll=index_ct;
      index_ct+=(psp->d_size*(psp->d_size+1)-(psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag))/2;
    }
    else {
      psp->has_ll = _FALSE_;
    }

    /* the computation of C_l^Tl would require a very good sampling of
       transfer functions over a wide range, and a huge computation
       time. In the current version, we prefer to switch it off, rather
       than either slowing down the code considerably, or producing
       very inaccurate spectra.

       if ((ppt->has_cl_cmb_temperature == _TRUE_) && (ppt->has_cl_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
       psp->has_tl = _TRUE_;
       psp->index_ct_tl=index_ct;
       index_ct+=psp->d_size;
       }
       else {
       psp->has_tl = _FALSE_;
       }
    */
    psp->has_tl = _FALSE_;

    if ((ppt->has_cl_number_count == _TRUE_) && (ppt->has_cl_lensing_potential == _TRUE_) && (ppt->has_scalars == _TRUE_)) {
      psp->has_dl = _TRUE_;
      psp->index_ct_dl=index_ct;
      index_ct += psp->d_size*psp->d_size - (psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag);
    }
    else {
      psp->has_dl = _FALSE_;
    }

    psp->ct_size = index_ct;

    /* infer from input quantities the l_max for each mode and type,
       l_max_ct[index_md][index_type].  Maximize it over index_ct, and
       then over index_md. */

    class_alloc(psp->l_max,sizeof(int*)*psp->md_size,psp->error_message);
    class_alloc(psp->l_max_ct,sizeof(int*)*psp->md_size,psp->error_message);
    for (index_md=0; index_md<psp->md_size; index_md++) {
      class_calloc(psp->l_max_ct[index_md],psp->ct_size,sizeof(int),psp->error_message);
    }

    if (ppt->has_scalars == _TRUE_) {

      /* spectra computed up to l_scalar_max */

      if (psp->has_tt == _TRUE_) psp->l_max_ct[ppt->index_md_scalars][psp->index_ct_tt] = ppt->l_scalar_max;
      if (psp->has_ee == _TRUE_) psp->l_max_ct[ppt->index_md_scalars][psp->index_ct_ee] = ppt->l_scalar_max;
      if (psp->has_te == _TRUE_) psp->l_max_ct[ppt->index_md_scalars][psp->index_ct_te] = ppt->l_scalar_max;
      if (psp->has_pp == _TRUE_) psp->l_max_ct[ppt->index_md_scalars][psp->index_ct_pp] = ppt->l_scalar_max;
      if (psp->has_tp == _TRUE_) psp->l_max_ct[ppt->index_md_scalars][psp->index_ct_tp] = ppt->l_scalar_max;
      if (psp->has_ep == _TRUE_) psp->l_max_ct[ppt->index_md_scalars][psp->index_ct_ep] = ppt->l_scalar_max;

      /* spectra computed up to l_lss_max */

      if (psp->has_dd == _TRUE_)
        for (index_ct=psp->index_ct_dd;
             index_ct<psp->index_ct_dd+(psp->d_size*(psp->d_size+1)-(psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag))/2;
             index_ct++)
          psp->l_max_ct[ppt->index_md_scalars][index_ct] = ppt->l_lss_max;

      if (psp->has_td == _TRUE_)
        for (index_ct=psp->index_ct_td;
             index_ct<psp->index_ct_td+psp->d_size;
             index_ct++)
          psp->l_max_ct[ppt->index_md_scalars][index_ct] = MIN(ppt->l_scalar_max,ppt->l_lss_max);

      if (psp->has_pd == _TRUE_)
        for (index_ct=psp->index_ct_pd;
             index_ct<psp->index_ct_pd+psp->d_size;
             index_ct++)
          psp->l_max_ct[ppt->index_md_scalars][index_ct] = MIN(ppt->l_scalar_max,ppt->l_lss_max);

      if (psp->has_ll == _TRUE_)
        for (index_ct=psp->index_ct_ll;
             index_ct<psp->index_ct_ll+(psp->d_size*(psp->d_size+1)-(psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag))/2;
             index_ct++)
          psp->l_max_ct[ppt->index_md_scalars][index_ct] = ppt->l_lss_max;

      if (psp->has_tl == _TRUE_)
        for (index_ct=psp->index_ct_tl;
             index_ct<psp->index_ct_tl+psp->d_size;
             index_ct++)
          psp->l_max_ct[ppt->index_md_scalars][index_ct] = MIN(ppt->l_scalar_max,ppt->l_lss_max);

      if (psp->has_dl == _TRUE_)
        for (index_ct=psp->index_ct_dl;
             index_ct < psp->index_ct_dl+(psp->d_size*psp->d_size - (psp->d_size-psp->non_diag)*(psp->d_size-1-psp->non_diag));
             index_ct++)
          psp->l_max_ct[ppt->index_md_scalars][index_ct] = ppt->l_lss_max;

    }
    if (ppt->has_tensors == _TRUE_) {

      /* spectra computed up to l_tensor_max */

      if (psp->has_tt == _TRUE_) psp->l_max_ct[ppt->index_md_tensors][psp->index_ct_tt] = ppt->l_tensor_max;
      if (psp->has_ee == _TRUE_) psp->l_max_ct[ppt->index_md_tensors][psp->index_ct_ee] = ppt->l_tensor_max;
      if (psp->has_te == _TRUE_) psp->l_max_ct[ppt->index_md_tensors][psp->index_ct_te] = ppt->l_tensor_max;
      if (psp->has_bb == _TRUE_) psp->l_max_ct[ppt->index_md_tensors][psp->index_ct_bb] = ppt->l_tensor_max;
    }

    /* maximizations */
    psp->l_max_tot = 0.;
    for (index_md=0; index_md < psp->md_size; index_md++) {
      psp->l_max[index_md] = 0.;
      for (index_ct=0.; index_ct<psp->ct_size; index_ct++)
        psp->l_max[index_md] = MAX(psp->l_max[index_md],psp->l_max_ct[index_md][index_ct]);
      psp->l_max_tot = MAX(psp->l_max_tot,psp->l_max[index_md]);
    }
  }

  /* indices for species associated with a matter transfer function in Fourier space */

  index_tr=0;
  class_define_index(psp->index_tr_delta_g,ppt->has_source_delta_g,index_tr,1);
  class_define_index(psp->index_tr_delta_b,ppt->has_source_delta_b,index_tr,1);
  class_define_index(psp->index_tr_delta_cdm,ppt->has_source_delta_cdm,index_tr,1);
  class_define_index(psp->index_tr_delta_dcdm,ppt->has_source_delta_dcdm,index_tr,1);
  class_define_index(psp->index_tr_delta_scf,ppt->has_source_delta_scf,index_tr,1);
  class_define_index(psp->index_tr_delta_fld,ppt->has_source_delta_fld,index_tr,1);
  class_define_index(psp->index_tr_delta_ur,ppt->has_source_delta_ur,index_tr,1);
  class_define_index(psp->index_tr_delta_dr,ppt->has_source_delta_dr,index_tr,1);
  class_define_index(psp->index_tr_delta_ncdm1,ppt->has_source_delta_ncdm,index_tr,pba->N_ncdm);
  class_define_index(psp->index_tr_delta_tot,ppt->has_density_transfers,index_tr,1);

  /* indices for species associated with a velocity transfer function in Fourier space */

  class_define_index(psp->index_tr_theta_g,ppt->has_source_theta_g,index_tr,1);
  class_define_index(psp->index_tr_theta_b,ppt->has_source_theta_b,index_tr,1);
  class_define_index(psp->index_tr_theta_cdm,ppt->has_source_theta_cdm,index_tr,1);
  class_define_index(psp->index_tr_theta_dcdm,ppt->has_source_theta_dcdm,index_tr,1);
  class_define_index(psp->index_tr_theta_scf,ppt->has_source_theta_scf,index_tr,1);
  class_define_index(psp->index_tr_theta_fld,ppt->has_source_theta_fld,index_tr,1);
  class_define_index(psp->index_tr_theta_ur,ppt->has_source_theta_ur,index_tr,1);
  class_define_index(psp->index_tr_theta_dr,ppt->has_source_theta_ur,index_tr,1);
  class_define_index(psp->index_tr_theta_ncdm1,ppt->has_source_theta_ncdm,index_tr,pba->N_ncdm);
  class_define_index(psp->index_tr_theta_tot,ppt->has_velocity_transfers,index_tr,1);

  psp->tr_size = index_tr;
  
  /* Indices for halomodel */
  /* Depend on the maximal order of n-point cf one wants to have as output */
  int i2; 
  i2=0;
  psp->index_z_out=i2;
  i2++;
  psp->index_k_out=i2;
  i2++;
  psp->index_p1h_out=i2;
  i2++;
  psp->index_p2h_out=i2;
  if (psp->hm_highest_order > 2){
    i2++;
    psp->index_b1h_out=i2;
    i2++;
    psp->index_b2h_out=i2;
    i2++;
    psp->index_b3h_out=i2;
  }
  if (psp->hm_highest_order > 3){ 
    i2++;
    psp->index_t1h_out=i2;
    i2++;
    psp->index_t2h_out=i2;
    i2++;
    psp->index_t3h_out=i2;
    i2++;
    psp->index_t4h_out=i2;
  }
  i2++;
  psp->number_output_hm =i2; /*Amount of output values (currently z,k,p1h,p2h,b1h,b2h,b3h,t1h,t2h,t3h,t4h) */
  
  
  int index_component;
  index_component = 0;
  psp->index_comp_cdm = index_component;
  index_component++;
  psp->index_comp_stars = index_component;
  index_component++;
  psp->index_comp_gas = index_component;
  index_component++;
  psp->index_comp_nu = index_component;
  index_component++;
  psp->index_comp_tot = index_component;
  

  return _SUCCESS_;
  
}

/**
 * This routine computes a table of values for all harmonic spectra C_l's,
 * given the transfer functions and primordial spectra.
 *
 * @param ppt Input : pointer to perturbation structure
 * @param ptr Input : pointer to transfers structure
 * @param ppm Input : pointer to primordial structure
 * @param psp Input/Output: pointer to spectra structure
 * @return the error status
 */

int spectra_cls(
                struct background * pba,
                struct perturbs * ppt,
                struct transfers * ptr,
                struct primordial * ppm,
                struct spectra * psp
                ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_ic1,index_ic2,index_ic1_ic2;
  int index_l;
  int index_ct;
  int cl_integrand_num_columns;

  double * cl_integrand; /* array with argument cl_integrand[index_k*cl_integrand_num_columns+1+psp->index_ct] */
  double * transfer_ic1; /* array with argument transfer_ic1[index_tt] */
  double * transfer_ic2; /* idem */
  double * primordial_pk;  /* array with argument primordial_pk[index_ic_ic]*/

  /* This code can be optionally compiled with the openmp option for parallel computation.
     Inside parallel regions, the use of the command "return" is forbidden.
     For error management, instead of "return _FAILURE_", we will set the variable below
     to "abort = _TRUE_". This will lead to a "return _FAILURE_" jus after leaving the
     parallel region. */
  int abort;

#ifdef _OPENMP
  /* instrumentation times */
  double tstart, tstop;
#endif

  /** - allocate pointers to arrays where results will be stored */

  class_alloc(psp->l_size,sizeof(int)*psp->md_size,psp->error_message);
  class_alloc(psp->cl,sizeof(double *)*psp->md_size,psp->error_message);
  class_alloc(psp->ddcl,sizeof(double *)*psp->md_size,psp->error_message);

  psp->l_size_max = ptr->l_size_max;
  class_alloc(psp->l,sizeof(double)*psp->l_size_max,psp->error_message);

  /** - store values of l */
  for (index_l=0; index_l < psp->l_size_max; index_l++) {
    psp->l[index_l] = (double)ptr->l[index_l];
  }

  /** - loop over modes (scalar, tensors, etc). For each mode: */

  for (index_md = 0; index_md < psp->md_size; index_md++) {

    /** - a) store number of l values for this mode */

    psp->l_size[index_md] = ptr->l_size[index_md];

    /** - b) allocate arrays where results will be stored */

    class_alloc(psp->cl[index_md],sizeof(double)*psp->l_size[index_md]*psp->ct_size*psp->ic_ic_size[index_md],psp->error_message);
    class_alloc(psp->ddcl[index_md],sizeof(double)*psp->l_size[index_md]*psp->ct_size*psp->ic_ic_size[index_md],psp->error_message);
    cl_integrand_num_columns = 1+psp->ct_size*2; /* one for k, ct_size for each type, ct_size for each second derivative of each type */

    /** d) loop over initial conditions */

    for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
      for (index_ic2 = index_ic1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {
        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);

        /* non-diagonal coefficients should be computed only if non-zero correlation */
        if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

          /* initialize error management flag */
          abort = _FALSE_;

          /* beginning of parallel region */

#pragma omp parallel                                                    \
  shared(ptr,ppm,index_md,psp,ppt,cl_integrand_num_columns,index_ic1,index_ic2,abort) \
  private(tstart,cl_integrand,primordial_pk,transfer_ic1,transfer_ic2,index_l,tstop)

          {

#ifdef _OPENMP
            tstart = omp_get_wtime();
#endif

            class_alloc_parallel(cl_integrand,
                                 ptr->q_size*cl_integrand_num_columns*sizeof(double),
                                 psp->error_message);

            class_alloc_parallel(primordial_pk,
                                 psp->ic_ic_size[index_md]*sizeof(double),
                                 psp->error_message);

            class_alloc_parallel(transfer_ic1,
                                 ptr->tt_size[index_md]*sizeof(double),
                                 psp->error_message);

            class_alloc_parallel(transfer_ic2,
                                 ptr->tt_size[index_md]*sizeof(double),
                                 psp->error_message);

#pragma omp for schedule (dynamic)

            /** - loop over l values defined in the transfer module.
                For each l, compute the C_l's for all types (TT, TE, ...)
                by convolving primordial spectra with transfer  functions.
                This elementary task is assigned to spectra_compute_cl() */

            for (index_l=0; index_l < ptr->l_size[index_md]; index_l++) {

#pragma omp flush(abort)

              class_call_parallel(spectra_compute_cl(pba,
                                                     ppt,
                                                     ptr,
                                                     ppm,
                                                     psp,
                                                     index_md,
                                                     index_ic1,
                                                     index_ic2,
                                                     index_l,
                                                     cl_integrand_num_columns,
                                                     cl_integrand,
                                                     primordial_pk,
                                                     transfer_ic1,
                                                     transfer_ic2),
                                  psp->error_message,
                                  psp->error_message);

            } /* end of loop over l */

#ifdef _OPENMP
            tstop = omp_get_wtime();
            if (psp->spectra_verbose > 1)
              printf("In %s: time spent in parallel region (loop over l's) = %e s for thread %d\n",
                     __func__,tstop-tstart,omp_get_thread_num());
#endif
            free(cl_integrand);

            free(primordial_pk);

            free(transfer_ic1);

            free(transfer_ic2);

          } /* end of parallel region */

          if (abort == _TRUE_) return _FAILURE_;

        }
        else {

          /* set non-diagonal coefficients to zero if pair of ic's uncorrelated */

          for (index_l=0; index_l < ptr->l_size[index_md]; index_l++) {
            for (index_ct=0; index_ct<psp->ct_size; index_ct++) {
              psp->cl[index_md]
                [(index_l * psp->ic_ic_size[index_md] + index_ic1_ic2) * psp->ct_size + index_ct]
                = 0.;
            }
          }
        }
      }
    }

    /** - e) now that for a given mode, all possible C_l's have been computed,
        compute second derivative of the array in which they are stored,
        in view of spline interpolation. */

    class_call(array_spline_table_lines(psp->l,
                                        psp->l_size[index_md],
                                        psp->cl[index_md],
                                        psp->ic_ic_size[index_md]*psp->ct_size,
                                        psp->ddcl[index_md],
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);
  }

  return _SUCCESS_;

}

/**
 * This routine computes the C_l's for a given mode, pair of initial conditions
 * and multipole, but for all types (TT, TE...), by convolving the
 * transfer functions with the primordial spectra.
 *
 * @param ppt           Input : pointer to perturbation structure
 * @param ptr           Input : pointer to transfers structure
 * @param ppm           Input : pointer to primordial structure
 * @param psp           Input/Output: pointer to spectra structure (result stored here)
 * @param index_md    Input : index of mode under consideration
 * @param index_ic1     Input : index of first initial condition in the correlator
 * @param index_ic2     Input : index of second initial condition in the correlato
 * @param index_l       Input : index of multipole under consideration
 * @param cl_integrand_num_column Input : number of columns in cl_integrand
 * @param cl_integrand  Input : an allocated workspace
 * @param primordial_pk Input : table of primordial spectrum values
 * @param transfer_ic1  Input : table of transfer function values for first initial condition
 * @param transfer_ic2  Input : table of transfer function values for second initial condition
 * @return the error status
 */

int spectra_compute_cl(
                       struct background * pba,
                       struct perturbs * ppt,
                       struct transfers * ptr,
                       struct primordial * ppm,
                       struct spectra * psp,
                       int index_md,
                       int index_ic1,
                       int index_ic2,
                       int index_l,
                       int cl_integrand_num_columns,
                       double * cl_integrand,
                       double * primordial_pk,
                       double * transfer_ic1,
                       double * transfer_ic2
                       ) {

  int index_q;
  int index_tt;
  int index_ct;
  int index_d1,index_d2;
  double k;
  double clvalue;
  int index_ic1_ic2;
  double transfer_ic1_temp=0.;
  double transfer_ic2_temp=0.;
  double * transfer_ic1_nc=NULL;
  double * transfer_ic2_nc=NULL;
  double factor;
  int index_q_spline=0;

  index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);

  if (ppt->has_cl_number_count == _TRUE_) {
    class_alloc(transfer_ic1_nc,psp->d_size*sizeof(double),psp->error_message);
    class_alloc(transfer_ic2_nc,psp->d_size*sizeof(double),psp->error_message);
  }

  for (index_q=0; index_q < ptr->q_size; index_q++) {

    //q = ptr->q[index_q];
    k = ptr->k[index_md][index_q];

    cl_integrand[index_q*cl_integrand_num_columns+0] = k;

    class_call(primordial_spectrum_at_k(ppm,index_md,linear,k,primordial_pk),
               ppm->error_message,
               psp->error_message);

    /* above routine checks that k>0: no possible division by zero below */

    for (index_tt=0; index_tt < ptr->tt_size[index_md]; index_tt++) {

      transfer_ic1[index_tt] =
        ptr->transfer[index_md]
        [((index_ic1 * ptr->tt_size[index_md] + index_tt)
          * ptr->l_size[index_md] + index_l)
         * ptr->q_size + index_q];

      if (index_ic1 == index_ic2) {
        transfer_ic2[index_tt] = transfer_ic1[index_tt];
      }
      else {
        transfer_ic2[index_tt] = ptr->transfer[index_md]
          [((index_ic2 * ptr->tt_size[index_md] + index_tt)
            * ptr->l_size[index_md] + index_l)
           * ptr->q_size + index_q];
      }
    }

    /* define combinations of transfer functions */

    if (ppt->has_cl_cmb_temperature == _TRUE_) {

      if (_scalars_) {

        transfer_ic1_temp = transfer_ic1[ptr->index_tt_t0] + transfer_ic1[ptr->index_tt_t1] + transfer_ic1[ptr->index_tt_t2];
        transfer_ic2_temp = transfer_ic2[ptr->index_tt_t0] + transfer_ic2[ptr->index_tt_t1] + transfer_ic2[ptr->index_tt_t2];

      }

      if (_vectors_) {

        transfer_ic1_temp = transfer_ic1[ptr->index_tt_t1] + transfer_ic1[ptr->index_tt_t2];
        transfer_ic2_temp = transfer_ic2[ptr->index_tt_t1] + transfer_ic2[ptr->index_tt_t2];

      }

      if (_tensors_) {

        transfer_ic1_temp = transfer_ic1[ptr->index_tt_t2];
        transfer_ic2_temp = transfer_ic2[ptr->index_tt_t2];

      }
    }

    if (ppt->has_cl_number_count == _TRUE_) {

      for (index_d1=0; index_d1<psp->d_size; index_d1++) {

        transfer_ic1_nc[index_d1] = 0.;
        transfer_ic2_nc[index_d1] = 0.;

        if (ppt->has_nc_density == _TRUE_) {
          transfer_ic1_nc[index_d1] += transfer_ic1[ptr->index_tt_density+index_d1];
          transfer_ic2_nc[index_d1] += transfer_ic2[ptr->index_tt_density+index_d1];
        }

        if (ppt->has_nc_rsd     == _TRUE_) {
          transfer_ic1_nc[index_d1]
            += transfer_ic1[ptr->index_tt_rsd+index_d1]
            + transfer_ic1[ptr->index_tt_d0+index_d1]
            + transfer_ic1[ptr->index_tt_d1+index_d1];
          transfer_ic2_nc[index_d1]
            += transfer_ic2[ptr->index_tt_rsd+index_d1]
            + transfer_ic2[ptr->index_tt_d0+index_d1]
            + transfer_ic2[ptr->index_tt_d1+index_d1];
        }

        if (ppt->has_nc_lens == _TRUE_) {
          transfer_ic1_nc[index_d1] +=
            psp->l[index_l]*(psp->l[index_l]+1.)*transfer_ic1[ptr->index_tt_nc_lens+index_d1];
          transfer_ic2_nc[index_d1] +=
            psp->l[index_l]*(psp->l[index_l]+1.)*transfer_ic2[ptr->index_tt_nc_lens+index_d1];
        }

        if (ppt->has_nc_gr == _TRUE_) {
          transfer_ic1_nc[index_d1]
            += transfer_ic1[ptr->index_tt_nc_g1+index_d1]
            + transfer_ic1[ptr->index_tt_nc_g2+index_d1]
            + transfer_ic1[ptr->index_tt_nc_g3+index_d1]
            + transfer_ic1[ptr->index_tt_nc_g4+index_d1]
            + transfer_ic1[ptr->index_tt_nc_g5+index_d1];
          transfer_ic2_nc[index_d1]
            += transfer_ic2[ptr->index_tt_nc_g1+index_d1]
            + transfer_ic2[ptr->index_tt_nc_g2+index_d1]
            + transfer_ic2[ptr->index_tt_nc_g3+index_d1]
            + transfer_ic2[ptr->index_tt_nc_g4+index_d1]
            + transfer_ic2[ptr->index_tt_nc_g5+index_d1];
        }

      }
    }

    /* integrand of Cl's */

    /* note: we must integrate

       C_l = int [4 pi dk/k calP(k) Delta1_l(q) Delta2_l(q)]

       where calP(k) is the dimensionless
       power spectrum equal to a constant in the scale-invariant case,
       and to P(k) = A_s k^(ns-1) otherwise and q=sqrt(k2+K) (scalars)
       or sqrt(k2+2K) (vectors) or sqrt(k2+3K) (tensors)

       In the literature, people often rewrite the integral in terms
       of q and absorb the Jacobian of the change of variables in a redefinition of the primodial
       spectrum. Let us illustrate this for scalars:

       dk/k = kdk/k2 = qdq/k2 = dq/q * (q/k)^2 = dq/q * [q2/(q2-K)] = q2dq * 1/[q(q2-K)]

       This factor 1/[q(q2-K)] is commonly absorbed in the definition of calP. Then one would have

       C_l = int [4 pi q2 dq {A_s k^(ns-1)/[q(q2-K)]} Delta1_l(q) Delta2_l(q)]

       Sometimes in the literature, the factor (k2-3K)=(q2-4K) present
       in the initial conditions of scalar transfer functions (if
       normalized to curvature R=1) is also absorbed in the definition
       of the power spectrum. Then the curvature power spectrum reads

       calP = (q2-4K)/[q(q2-K)] * (k/k)^ns

       In CLASS we prefer to define calP = (k/k)^ns like in the flat
       case, to have the factor (q2-4K) in the initialk conditions,
       and the factor 1/[q(q2-K)] doesn't need to be there since we
       integrate over dk/k.

       For tensors, the change of variable described above gives a slightly different result:

       dk/k = kdk/k2 = qdq/k2 = dq/q * (q/k)^2 = dq/q * [q2/(q2-3K)] = q2dq * 1/[q(q2-3K)]

       But for tensors there are extra curvature-related correction factors to
       take into account. See the comments in the perturbation module,
       related to initial conditions for tensors.

    */

    factor = 4. * _PI_ / k;

    if (psp->has_tt == _TRUE_)
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_tt]=
        primordial_pk[index_ic1_ic2]
        * transfer_ic1_temp
        * transfer_ic2_temp
        * factor;

    if (psp->has_ee == _TRUE_)
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_ee]=
        primordial_pk[index_ic1_ic2]
        * transfer_ic1[ptr->index_tt_e]
        * transfer_ic2[ptr->index_tt_e]
        * factor;

    if (psp->has_te == _TRUE_)
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_te]=
        primordial_pk[index_ic1_ic2]
        * 0.5*(transfer_ic1_temp * transfer_ic2[ptr->index_tt_e] +
               transfer_ic1[ptr->index_tt_e] * transfer_ic2_temp)
        * factor;

    if (_tensors_ && (psp->has_bb == _TRUE_))
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_bb]=
        primordial_pk[index_ic1_ic2]
        * transfer_ic1[ptr->index_tt_b]
        * transfer_ic2[ptr->index_tt_b]
        * factor;

    if (_scalars_ && (psp->has_pp == _TRUE_))
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_pp]=
        primordial_pk[index_ic1_ic2]
        * transfer_ic1[ptr->index_tt_lcmb]
        * transfer_ic2[ptr->index_tt_lcmb]
        * factor;

    if (_scalars_ && (psp->has_tp == _TRUE_))
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_tp]=
        primordial_pk[index_ic1_ic2]
        * 0.5*(transfer_ic1_temp * transfer_ic2[ptr->index_tt_lcmb] +
               transfer_ic1[ptr->index_tt_lcmb] * transfer_ic2_temp)
        * factor;

    if (_scalars_ && (psp->has_ep == _TRUE_))
      cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_ep]=
        primordial_pk[index_ic1_ic2]
        * 0.5*(transfer_ic1[ptr->index_tt_e] * transfer_ic2[ptr->index_tt_lcmb] +
               transfer_ic1[ptr->index_tt_lcmb] * transfer_ic2[ptr->index_tt_e])
        * factor;

    if (_scalars_ && (psp->has_dd == _TRUE_)) {
      index_ct=0;
      for (index_d1=0; index_d1<psp->d_size; index_d1++) {
        for (index_d2=index_d1; index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++) {
          cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_dd+index_ct]=
            primordial_pk[index_ic1_ic2]
            * transfer_ic1_nc[index_d1]
            * transfer_ic2_nc[index_d2]
            * factor;
          index_ct++;
        }
      }
    }

    if (_scalars_ && (psp->has_td == _TRUE_)) {
      for (index_d1=0; index_d1<psp->d_size; index_d1++) {
        cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_td+index_d1]=
          primordial_pk[index_ic1_ic2]
          * 0.5*(transfer_ic1_temp * transfer_ic2_nc[index_d1] +
                 transfer_ic1_nc[index_d1] * transfer_ic2_temp)
          * factor;
      }
    }

    if (_scalars_ && (psp->has_pd == _TRUE_)) {
      for (index_d1=0; index_d1<psp->d_size; index_d1++) {
        cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_pd+index_d1]=
          primordial_pk[index_ic1_ic2]
          * 0.5*(transfer_ic1[ptr->index_tt_lcmb] * transfer_ic2_nc[index_d1] +
                 transfer_ic1_nc[index_d1] * transfer_ic2[ptr->index_tt_lcmb])
          * factor;
      }
    }

    if (_scalars_ && (psp->has_ll == _TRUE_)) {
      index_ct=0;
      for (index_d1=0; index_d1<psp->d_size; index_d1++) {
        for (index_d2=index_d1; index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++) {
          cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_ll+index_ct]=
            primordial_pk[index_ic1_ic2]
            * transfer_ic1[ptr->index_tt_lensing+index_d1]
            * transfer_ic2[ptr->index_tt_lensing+index_d2]
            * factor;
          index_ct++;
        }
      }
    }

    if (_scalars_ && (psp->has_tl == _TRUE_)) {
      for (index_d1=0; index_d1<psp->d_size; index_d1++) {
        cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_tl+index_d1]=
          primordial_pk[index_ic1_ic2]
          * 0.5*(transfer_ic1_temp * transfer_ic2[ptr->index_tt_lensing+index_d1] +
                 transfer_ic1[ptr->index_tt_lensing+index_d1] * transfer_ic2_temp)
          * factor;
      }
    }

    if (_scalars_ && (psp->has_dl == _TRUE_)) {
      index_ct=0;
      for (index_d1=0; index_d1<psp->d_size; index_d1++) {
        for (index_d2=MAX(index_d1-psp->non_diag,0); index_d2<=MIN(index_d1+psp->non_diag,psp->d_size-1); index_d2++) {
          cl_integrand[index_q*cl_integrand_num_columns+1+psp->index_ct_dl+index_ct]=
            primordial_pk[index_ic1_ic2]
            * transfer_ic1_nc[index_d1] * transfer_ic2[ptr->index_tt_lensing+index_d2]
            * factor;
          index_ct++;
        }
      }
    }
  }

  for (index_ct=0; index_ct<psp->ct_size; index_ct++) {

    /* treat null spectra (C_l^BB of scalars, C_l^pp of tensors, etc. */

    if ((_scalars_ && (psp->has_bb == _TRUE_) && (index_ct == psp->index_ct_bb)) ||
        (_tensors_ && (psp->has_pp == _TRUE_) && (index_ct == psp->index_ct_pp)) ||
        (_tensors_ && (psp->has_tp == _TRUE_) && (index_ct == psp->index_ct_tp)) ||
        (_tensors_ && (psp->has_ep == _TRUE_) && (index_ct == psp->index_ct_ep)) ||
        (_tensors_ && (psp->has_dd == _TRUE_) && (index_ct == psp->index_ct_dd)) ||
        (_tensors_ && (psp->has_td == _TRUE_) && (index_ct == psp->index_ct_td)) ||
        (_tensors_ && (psp->has_pd == _TRUE_) && (index_ct == psp->index_ct_pd)) ||
        (_tensors_ && (psp->has_ll == _TRUE_) && (index_ct == psp->index_ct_ll)) ||
        (_tensors_ && (psp->has_tl == _TRUE_) && (index_ct == psp->index_ct_tl)) ||
        (_tensors_ && (psp->has_dl == _TRUE_) && (index_ct == psp->index_ct_dl))
        ) {

      psp->cl[index_md]
        [(index_l * psp->ic_ic_size[index_md] + index_ic1_ic2) * psp->ct_size + index_ct] = 0.;

    }
    /* for non-zero spectra, integrate over q */
    else {

      /* spline the integrand over the whole range of k's */

      class_call(array_spline(cl_integrand,
                              cl_integrand_num_columns,
                              ptr->q_size,
                              0,
                              1+index_ct,
                              1+psp->ct_size+index_ct,
                              _SPLINE_EST_DERIV_,
                              psp->error_message),
                 psp->error_message,
                 psp->error_message);

      /* Technical point: we will now do a spline integral over the
         whole range of k's, excepted in the closed (K>0) case. In
         that case, it is a bad idea to spline over the values of k
         corresponding to nu<nu_flat_approximation. In this region, nu
         values are integer values, so the steps dq and dk have some
         discrete jumps. This makes the spline routine less accurate
         than a trapezoidal integral with finer sampling. So, in the
         closed case, we set index_q_spline to
         ptr->index_q_flat_approximation, to tell the integration
         routine that below this index, it should treat the integral
         as a trapezoidal one. For testing, one is free to set
         index_q_spline to 0, to enforce spline integration
         everywhere, or to (ptr->q_size-1), to enforce trapezoidal
         integration everywhere. */

      if (pba->sgnK == 1) {
        index_q_spline = ptr->index_q_flat_approximation;
      }

      class_call(array_integrate_all_trapzd_or_spline(cl_integrand,
                                                      cl_integrand_num_columns,
                                                      ptr->q_size,
                                                      index_q_spline,
                                                      0,
                                                      1+index_ct,
                                                      1+psp->ct_size+index_ct,
                                                      &clvalue,
                                                      psp->error_message),
                 psp->error_message,
                 psp->error_message);

      /* in the closed case, instead of an integral, we have a
         discrete sum. In practise, this does not matter: the previous
         routine does give a correct approximation of the discrete
         sum, both in the trapezoidal and spline regions. The only
         error comes from the first point: the previous routine
         assumes a weight for the first point which is too small
         compared to what it would be in the an actual discrete
         sum. The line below correct this problem in an exact way.
      */

      if (pba->sgnK == 1) {
        clvalue += cl_integrand[1+index_ct] * ptr->q[0]/ptr->k[0][0]*sqrt(pba->K)/2.;
      }

      /* we have the correct C_l now. We can store it in the transfer structure. */

      psp->cl[index_md]
        [(index_l * psp->ic_ic_size[index_md] + index_ic1_ic2) * psp->ct_size + index_ct]
        = clvalue;

    }
  }

  if (ppt->has_cl_number_count == _TRUE_) {
    free(transfer_ic1_nc);
    free(transfer_ic2_nc);
  }

  return _SUCCESS_;

}

/**
 * This routine computes the values of k and tau at which the matter
 * power spectra P(k,tau) and the matter transfer functions T_i(k,tau)
 * will be stored.
 *
 * @param pba Input : pointer to background structure (for z to tau conversion)
 * @param ppt Input : pointer to perturbation structure (contain source functions)
 * @param psp Input/Output: pointer to spectra structure
 * @return the error status
 */

int spectra_k_and_tau(
                      struct background * pba,
                      struct perturbs * ppt,
                      struct spectra * psp
                      ) {

  /** Summary: */

  /** - define local variables */

  int index_k;
  int index_tau;
  double tau_min;

  /** - check the presence of scalar modes */

  class_test((ppt->has_scalars == _FALSE_),
             psp->error_message,
             "you cannot ask for matter power spectrum since you turned off scalar modes");

  /** - check the maximum redshift z_max_pk at which P(k,z) and T_i(k,z) should be
      computable by interpolation. If it is equal to zero, only P(k,z=0)
      needs to be computed. If it is higher, we will store in a table
      various P(k,tau) at several values of tau generously encompassing
      the range 0<z<z_max_pk */

  /* if z_max_pk<0, return error */
  class_test((psp->z_max_pk < 0),
             psp->error_message,
             "asked for negative redshift z=%e",psp->z_max_pk);

  /* if z_max_pk=0, there is just one value to store */
  if (psp->z_max_pk == 0.) {
    psp->ln_tau_size=1;
  }

  /* if z_max_pk>0, store several values (with a confortable margin above z_max_pk) in view of interpolation */
  else{

    /* find the first relevant value of tau (last value in the table tau_ampling before tau(z_max)) and infer the number of values of tau at which P(k) must be stored */

    class_call(background_tau_of_z(pba,psp->z_max_pk,&tau_min),
               pba->error_message,
               psp->error_message);

    index_tau=0;
    class_test((tau_min < ppt->tau_sampling[index_tau]),
               psp->error_message,
               "you asked for zmax=%e, i.e. taumin=%e, smaller than first possible value =%e",psp->z_max_pk,tau_min,ppt->tau_sampling[0]);

    while (ppt->tau_sampling[index_tau] < tau_min){
      index_tau++;
    }
    index_tau --;
    /* whenever possible, take a few more values in to avoid boundary effects in the interpolation */
    if (index_tau>0) index_tau--;
    if (index_tau>0) index_tau--;
    if (index_tau>0) index_tau--;
    if (index_tau>0) index_tau--;
    psp->ln_tau_size=ppt->tau_size-index_tau;

  }

  /** - allocate and fill table of tau values at which P(k,tau) and T_i(k,tau) are stored */

  class_alloc(psp->ln_tau,sizeof(double)*psp->ln_tau_size,psp->error_message);

  for (index_tau=0; index_tau<psp->ln_tau_size; index_tau++) {
    psp->ln_tau[index_tau]=log(ppt->tau_sampling[index_tau-psp->ln_tau_size+ppt->tau_size]);
  }

  /** - allocate and fill table of k values at which P(k,tau) is stored */

  psp->ln_k_size = ppt->k_size[ppt->index_md_scalars];
  class_alloc(psp->ln_k,sizeof(double)*psp->ln_k_size,psp->error_message);

  for (index_k=0; index_k<psp->ln_k_size; index_k++) {
    class_test(ppt->k[ppt->index_md_scalars][index_k] <= 0.,
               psp->error_message,
               "stop to avoid segmentation fault");
    psp->ln_k[index_k]=log(ppt->k[ppt->index_md_scalars][index_k]);
  }

  return _SUCCESS_;
}

/**
 * This routine computes a table of values for all matter power spectra P(k),
 * given the source functions and primordial spectra.
 *
 * @param pba Input : pointer to background structure (will provide H, Omega_m at redshift of interest)
 * @param ppt Input : pointer to perturbation structure (contain source functions)
 * @param ppm Input : pointer to primordial structure
 * @param psp Input/Output: pointer to spectra structure
 * @return the error status
 */

int spectra_pk(
               struct background * pba,
               struct perturbs * ppt,
               struct primordial * ppm,
               struct nonlinear * pnl,
               struct spectra * psp
               ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_ic1,index_ic2,index_ic1_ic2;
  int index_k;
  int index_tau;
  int last_index;
  double * primordial_pk; /* array with argument primordial_pk[index_ic_ic] */
  double source_ic1;
  double source_ic2;
  double ln_pk_tot;

  /** - check the presence of scalar modes */

  class_test((ppt->has_scalars == _FALSE_),
             psp->error_message,
             "you cannot ask for matter power spectrum since you turned off scalar modes");

  index_md = psp->index_md_scalars;

  /** - allocate temporary vectors where the primordial spectrum and the background quantitites will be stored */

  class_alloc(primordial_pk,psp->ic_ic_size[index_md]*sizeof(double),psp->error_message);

  /** - allocate and fill array of P(k,tau) values */

  class_alloc(psp->ln_pk,
              sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_ic_size[index_md],
              psp->error_message);

  if (pnl->method != nl_none) {
    class_alloc(psp->ln_pk_nl,
                sizeof(double)*psp->ln_tau_size*psp->ln_k_size,
                psp->error_message);
  }
  else {
    psp->ln_pk_nl = NULL;
  }

  for (index_tau=0 ; index_tau < psp->ln_tau_size; index_tau++) {
    for (index_k=0; index_k<psp->ln_k_size; index_k++) {

      class_call(primordial_spectrum_at_k(ppm,index_md,logarithmic,psp->ln_k[index_k],primordial_pk),
                 ppm->error_message,
                 psp->error_message);

      ln_pk_tot =0;

      /* curvature primordial spectrum:
         P_R(k) = 1/(2pi^2) k^3 <R R>
         so, primordial curvature correlator:
         <R R> = (2pi^2) k^-3 P_R(k)
         so, delta_m correlator:
         P(k) = <delta_m delta_m> = (2pi^2) k^-3 (source_m)^2 P_R(k)

         For isocurvature or cross adiabatic-isocurvature parts,
         replace one or two 'R' by 'S_i's */

      /* part diagonal in initial conditions */
      for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {

        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,psp->ic_size[index_md]);

        source_ic1 = ppt->sources[index_md]
          [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
          [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

        psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2] =
          log(2.*_PI_*_PI_/exp(3.*psp->ln_k[index_k])
              *source_ic1*source_ic1
              *exp(primordial_pk[index_ic1_ic2]));

        ln_pk_tot += psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2];
	

      }

      /* part non-diagonal in initial conditions */
      for (index_ic1 = 0; index_ic1 < psp->ic_size[index_md]; index_ic1++) {
        for (index_ic2 = index_ic1+1; index_ic2 < psp->ic_size[index_md]; index_ic2++) {

          index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,psp->ic_size[index_md]);

          if (psp->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

            source_ic1 = ppt->sources[index_md]
              [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            source_ic2 = ppt->sources[index_md]
              [index_ic2 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2] =
              primordial_pk[index_ic1_ic2]*SIGN(source_ic1)*SIGN(source_ic2);

            ln_pk_tot += psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2];
          }
          else {
            psp->ln_pk[(index_tau * psp->ln_k_size + index_k)* psp->ic_ic_size[index_md] + index_ic1_ic2] = 0.;
          }
        }
      }

      /* if non-linear corrections required, compute the total non-linear matter power spectrum */

      if (pnl->method == nl_halofit) {

        psp->ln_pk_nl[index_tau * psp->ln_k_size + index_k] =
          ln_pk_tot
          + 2.*log(pnl->nl_corr_density[(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k]);

      }
      /*
      if (pnl->method == nl_halomodel) {
	
	for (index_tau=0;index_tau<psp->ln_tau_size;index_tau++){
	  double * pvecback;
	  double * pvecback_sp_long;
	  double pk, pk_ic;
	  class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);
	  class_call(background_at_tau(pba,exp(psp->ln_tau[index_tau]),pba->long_info,pba->inter_normal,&last_index,pvecback_sp_long),
             pba->error_message,
             pnl->error_message);
          double z = pba->a_today/pvecback[pba->index_bg_a]-1.;
	  for (index_k=0;index_k<psp->ln_k_size;index_k++){    
	    double k = exp(psp->ln_k[index_k]);
            class_call(spectra_pk_hm_at_k_and_z(pba,ppm,psp,k,z,&pk),
               psp->error_message,
               psp->error_message);
            psp->ln_pk_hm[index_tau * psp->ln_k_size + index_k] = log(pk);
          }
        }
      }*/
    }
  }

  /**- if interpolation of P(k,tau) will be needed (as a function of tau),
     compute array of second derivatives in view of spline interpolation */

  if (psp->ln_tau_size > 1) {

    class_alloc(psp->ddln_pk,sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_ic_size[index_md],psp->error_message);

    class_call(array_spline_table_lines(psp->ln_tau,
                                        psp->ln_tau_size,
                                        psp->ln_pk,
                                        psp->ic_ic_size[index_md]*psp->ln_k_size,
                                        psp->ddln_pk,
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);

  }
  
  /*int l,j;
      for(j=0;j<psp->ln_tau_size;j++){  
        for(l=0;l<psp->ln_k_size;l++){
          printf("l=%f,j=%f\n",(double)(l),(double)(j));
          printf("ddlnPklin[k=%f,z=%f] = %f\n",exp(psp->ln_k[psp->ln_k_size-1])/pba->h,j,psp->ddln_pk[j*psp->ln_k_size+l]);
        } 
      }*/

  /* compute sigma8 (mean variance today in sphere of radius 8/h Mpc */
  class_call(spectra_sigma(pba,ppm,psp,8./pba->h,0.,&(psp->sigma8)),
             psp->error_message,
             psp->error_message);

  if (psp->spectra_verbose>0)
    fprintf(stdout," -> sigma8=%g (computed till k = %g h/Mpc)\n",
            psp->sigma8,
            exp(psp->ln_k[psp->ln_k_size-1])/pba->h);
  

  /**- if interpolation of P_NL(k,tau) will be needed (as a function of tau),
     compute array of second derivatives in view of spline interpolation */

  if ((pnl->method != nl_none)&&(pnl->method != nl_halomodel)) {
    if (psp->ln_tau_size > 1) {

      class_alloc(psp->ddln_pk_nl,sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_ic_size[index_md],psp->error_message);

      class_call(array_spline_table_lines(psp->ln_tau,
                                          psp->ln_tau_size,
                                          psp->ln_pk_nl,
                                          psp->ln_k_size,
                                          psp->ddln_pk_nl,
                                          _SPLINE_EST_DERIV_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);
      
      
    } 
  }

  free (primordial_pk);

  return _SUCCESS_;
}

/**
 * This routine computes sigma(R) given P(k) (does not check that k_max is large
 * enough)
 *
 * @param pba   Input: pointer to background structure
 * @param ppm   Input: pointer to primordial structure
 * @param psp   Input: pointer to spectra structure
 * @param z     Input: redhsift
 * @param R     Input: radius in Mpc
 * @param sigma Output: variance in a sphere of radius R (dimensionless)
 */

int spectra_sigma(
                  struct background * pba,
                  struct primordial * ppm,
                  struct spectra * psp,
                  double R,
                  double z,
                  double * sigma
                  ) {

  double pk;
  double * pk_ic = NULL;

  double * array_for_sigma;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i;

  double k,W,x;

  if (psp->ic_ic_size[psp->index_md_scalars]>1)
    class_alloc(pk_ic,
                psp->ic_ic_size[psp->index_md_scalars]*sizeof(double),
                psp->error_message);

  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              psp->ln_k_size*index_num*sizeof(double),
              psp->error_message);

  for (i=0;i<psp->ln_k_size;i++) {
    k=exp(psp->ln_k[i]);
    if (i == (psp->ln_k_size-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
    W=3./x/x/x*(sin(x)-x*cos(x));
    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,k,z,&pk,pk_ic),
               psp->error_message,
               psp->error_message);
    array_for_sigma[i*index_num+index_k]=k;
    array_for_sigma[i*index_num+index_y]=k*k*pk*W*W;
  }

  class_call(array_spline(array_for_sigma,
                          index_num,
                          psp->ln_k_size,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          psp->error_message),
             psp->error_message,
             psp->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        psp->ln_k_size,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma,
                                        psp->error_message),
             psp->error_message,
             psp->error_message);

  free(array_for_sigma);

  if (psp->ic_ic_size[psp->index_md_scalars]>1)
    free(pk_ic);

  *sigma = sqrt(*sigma/(2.*_PI_*_PI_));

  return _SUCCESS_;

}

/* Gives interpolation table for spectra_sigma; called only if halomodel chosen */ 
int spectra_sigma_table(struct background * pba,
	         	struct primordial * ppm,
			struct spectra * psp				
                        ) {
     
  int index_R, index_z, last_index_back;
  double Norm, sigma1, tau; 
  double * pvecback_sp_long;
  double * D_plus;
  double * sigmaarrayz0;
  class_alloc(D_plus,sizeof(double)*psp->sigma_z_slices,psp->error_message);
  class_alloc(psp->sigma_zarray,sizeof(double)*psp->sigma_z_slices,psp->error_message);
  class_alloc(psp->sigma_ln_tauarray,sizeof(double)*psp->sigma_z_slices,psp->error_message);
  class_alloc(pvecback_sp_long,pba->bg_size*sizeof(double),psp->error_message);
  
  /* Calculate maximal and minimal R_to_m */
  /* m_min/m_max must be superset of m_min(m_max of hm in order to work lateron (equal by default) */
  for (index_z=0;index_z<psp->sigma_z_slices;index_z++){
    double z = (double)(index_z)/((double)(psp->sigma_z_slices - 1) / (double)(psp->z_max_pk));
    psp->sigma_zarray[index_z] = z;
    
    class_call(background_tau_of_z(pba,z,&tau), pba->error_message, psp->error_message);  
    psp->sigma_ln_tauarray[index_z] = log(tau);
    
    class_call(background_at_tau(pba,
                                 tau,
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback_sp_long),
               pba->error_message,
               psp->error_message); 
    
    if (index_z == 0){
      Norm = pvecback_sp_long[pba->index_bg_D];
      psp->sigma_R_max = pow(10,1./3. * (log10(3. / 4. / _PI_) - log10((pba->Omega0_cdm + pba->Omega0_b) * pow(1+0,3) * pow(pvecback_sp_long[pba->index_bg_H],2)  / (8.0 / 3.0 *_PI_*_G_) * pow(_c_,2) / _M_SUN_ * _Mpc_over_m_ ) + log10(psp->sigma_m_max)));
      /*printf("r_max=%f\n",psp->sigma_R_max);*/
      
    }
    
    if (index_z == psp->sigma_z_slices - 1){
      psp->sigma_R_min = pow(10,1./3. * (log10(3. / 4. / _PI_) - log10((pba->Omega0_cdm + pba->Omega0_b) * pow(1+psp->z_max_pk,3) * pow(pvecback_sp_long[pba->index_bg_H],2)  / (8.0 / 3.0 *_PI_*_G_) * pow(_c_,2) / _M_SUN_ * _Mpc_over_m_ ) + log10(psp->sigma_m_min)));
     /* printf("r_min=%f\n",psp->sigma_R_min);*/ 
    }
    
    D_plus[index_z] =  pvecback_sp_long[pba->index_bg_D] / Norm;
    /*printf("D+=%f\n",D_plus[index_z]);*/
  }
  
  /* Allocate table of sigma for different R & z */
  class_alloc(psp->sigma_Rarray,psp->sigma_R_slices*sizeof(double),psp->error_message);
  class_alloc(sigmaarrayz0,psp->sigma_R_slices*sizeof(double),psp->error_message);
  class_alloc(psp->lnsigmatable,psp->sigma_R_slices*psp->sigma_z_slices*sizeof(double),psp->error_message);
  
  for (index_R=0;index_R<psp->sigma_R_slices;index_R++){
    psp->sigma_Rarray[index_R] = pow(10,log10(psp->sigma_R_min) + log10(psp->sigma_R_max/psp->sigma_R_min) * (double)((double)(index_R)/(double)(psp->sigma_R_slices - 1)));
    /*printf("R=%f\n", pow(10,log10(psp->sigma_R_min) + log10(psp->sigma_R_max/psp->sigma_R_min) * (double)((double)(index_R)/(double)(psp->sigma_R_slices - 1))));*/
  }	      
  for (index_R=0;index_R<psp->sigma_R_slices;index_R++){
    class_call(spectra_sigma(pba,ppm,psp,
                              psp->sigma_Rarray[index_R],0,
                              &sigma1),
                psp->error_message,
                psp->error_message);
       
    sigmaarrayz0[index_R] = sigma1;
    /*printf("sigma = %f\n",sigma1);*/
  }
       
  for (index_z=0;index_z<psp->sigma_z_slices;index_z++){     
    for (index_R=0;index_R<psp->sigma_R_slices;index_R++){
	 psp->lnsigmatable[index_z*psp->sigma_R_slices+index_R] = log(sigmaarrayz0[index_R]*D_plus[index_z]);   
	 /*printf("sigma(z=%f, R=%f) = %f\n",psp->sigma_zarray[index_z],psp->sigma_Rarray[index_R],psp->lnsigmatable[index_z*psp->sigma_R_slices+index_R]);*/
    }
  }
  
  /* Calculate 2nd derivatives needed for splining afterwards */
  class_alloc(psp->ddlnsigmatable,sizeof(double)*psp->sigma_z_slices*psp->sigma_R_slices,psp->error_message);
  class_call(array_spline_table_lines(psp->sigma_ln_tauarray,
                                        psp->sigma_z_slices,
                                        psp->lnsigmatable,
                                        psp->sigma_R_slices,
                                        psp->ddlnsigmatable,
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);
  
  /*for(index_z=0;index_z<psp->sigma_z_slices;index_z++){
    for(index_R = 0;index_R<psp->sigma_R_slices;index_R++){
      printf("dd(z=%f, R=%f) = %f\n",psp->sigma_zarray[index_z],psp->sigma_Rarray[index_R],psp->ddlnsigmatable[index_z*psp->sigma_R_slices+index_R]);
    }
  }*/
  
  free(pvecback_sp_long);
  free(D_plus);
  free(sigmaarrayz0);
  return _SUCCESS_;
}

int spectra_sigma_at_z(
                    struct background * pba,
                    struct spectra * psp,
                    double z,
                    double * output_tot /* array with argument output_tot[index_k] (must be already allocated) */
                    ) {

  /** Summary: */

  /** - define local variables */

  int last_index;
  int index_R;
  double tau1,ln_tau1;

  /** - first step: convert z into ln(tau) */

  class_call(background_tau_of_z(pba,z,&tau1),
             pba->error_message,
             psp->error_message);

  class_test(tau1 <= 0.,
             psp->error_message,
             "negative or null value of conformal time: cannot interpolate");

  ln_tau1 = log(tau1);
  
  
  /** - second step: for both modes (linear or logarithmic), store the spectrum in logarithmic format in the output array(s) */

  /**   (a.) if only values at tau=tau_today are stored and we want P(k,z=0), no need to interpolate */

  if (psp->sigma_z_slices == 1) {
    class_test(z != 0.,
               psp->error_message,
               "asked z=%e but only P(k,z=0) has been tabulated",z);
    for (index_R=0; index_R<psp->sigma_R_slices; index_R++)
      	output_tot[index_R] = psp->lnsigmatable[index_R];
           
  }

  /**   (b.) if several values of tau have been stored, use interpolation routine to get spectra at correct redshift */

  else { 
    
    int index_z;
    for (index_z=0;index_z<psp->sigma_z_slices;index_z++)
      /*printf("tau(%f) = %f \n",psp->hm_zarray[j],psp->ln_tau_hm[j]);*/
    
    class_call(array_interpolate_spline(psp->sigma_ln_tauarray,
                                          psp->sigma_z_slices,
                                          psp->lnsigmatable,
                                          psp->ddlnsigmatable,
                                          psp->sigma_R_slices,
                                          ln_tau1,
                                          &last_index,
                                          output_tot,
                                          psp->sigma_R_slices,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);   
    
  }
  
     
  return _SUCCESS_;

}


int spectra_sigma_at_R_and_z(
                          struct background * pba,
                          struct spectra * psp,
                          double R,
                          double z,
                          double * sigma /* pointer to a single number (must be already allocated) */
                          ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_k;
  int last_index;
  int index_ic1,index_ic2,index_ic1_ic2;

  double * sigmaarray_at_z = NULL;
  double * spline;

  /** - first step: check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_pk_at_z()) */

  class_test((R < psp->sigma_R_min ) || (R > psp->sigma_R_max ),
             psp->error_message,
             "R=%e out of bounds [%e:%e]",R,psp->sigma_R_min,psp->sigma_R_max);

  class_test((z < 0 ) || (z > psp->z_max_pk ),
             psp->error_message,
             "R=%e out of bounds [%e:%e]",z,0.,psp->z_max_pk);


    /* compute P(k,z) (in logarithmic format for more accurate interpolation) */
    class_alloc(sigmaarray_at_z,
                psp->sigma_R_slices*sizeof(double),
                psp->error_message);
    
    class_call(spectra_sigma_at_z(pba,
                           psp,
                           z,
                           sigmaarray_at_z),
               psp->error_message,
               psp->error_message);

    /* get its second derivatives with spline, then interpolate, then convert to linear format */

    class_alloc(spline,
                sizeof(double)*psp->sigma_R_slices,
                psp->error_message);

    class_call(array_spline_table_lines(psp->sigma_Rarray,
                                          psp->sigma_R_slices,
                                          sigmaarray_at_z,
                                          1,
                                          spline,
                                          _SPLINE_NATURAL_,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

    class_call(array_interpolate_spline(psp->sigma_Rarray,
                                          psp->sigma_R_slices,
                                          sigmaarray_at_z,
                                          spline,
                                          1,
                                          R,
                                          &last_index,
                                          sigma,
                                          1,
                                          psp->error_message),
                 psp->error_message,
                 psp->error_message);

    *sigma = exp(*sigma);

    free(sigmaarray_at_z);
    free(spline);


  return _SUCCESS_;

}

            

/**
 * This routine computes a table of values for all matter power spectra P(k),
 * given the source functions and primordial spectra.
 *
 * @param pba Input : pointer to background structure (will provide density of each species)
 * @param ppt Input : pointer to perturbation structure (contain source functions)
 * @param psp Input/Output: pointer to spectra structure
 * @return the error status
 */

int spectra_matter_transfers(
                             struct background * pba,
                             struct perturbs * ppt,
                             struct spectra * psp
                             ) {

  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_ic;
  int index_k;
  int index_tau;
  int last_index_back;
  double * pvecback_sp_long; /* array with argument pvecback_sp_long[pba->index_bg] */
  double delta_i,theta_i,rho_i;
  double delta_rho_tot,rho_tot;
  double rho_plus_p_theta_tot,rho_plus_p_tot;
  int n_ncdm;

  /** - check the presence of scalar modes */

  class_test((ppt->has_scalars == _FALSE_),
             psp->error_message,
             "you cannot ask for matter power spectrum since you turned off scalar modes");

  index_md = psp->index_md_scalars;

  /** - allocate and fill array of T_i(k,tau) values */

  class_alloc(psp->matter_transfer,sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_size[index_md]*psp->tr_size,psp->error_message);

  /** - allocate temporary vectors where the background quantitites will be stored */

  class_alloc(pvecback_sp_long,pba->bg_size*sizeof(double),psp->error_message);

  for (index_tau=0 ; index_tau < psp->ln_tau_size; index_tau++) {

    class_call(background_at_tau(pba,
                                 ppt->tau_sampling[index_tau-psp->ln_tau_size+ppt->tau_size],
                                 /* for this last argument we could have passed
                                    exp(psp->ln_tau[index_tau]) but we would then loose
                                    precision in the exp(log(x)) operation) */
                                 pba->long_info,
                                 pba->inter_normal,
                                 &last_index_back,
                                 pvecback_sp_long),
               pba->error_message,
               psp->error_message);

    for (index_k=0; index_k<psp->ln_k_size; index_k++) {

      for (index_ic = 0; index_ic < psp->ic_size[index_md]; index_ic++) {

        delta_rho_tot=0.;
        rho_tot=0.;
        rho_plus_p_theta_tot=0.;
        rho_plus_p_tot=0.;

        /* T_g(k,tau) */

        rho_i = pvecback_sp_long[pba->index_bg_rho_g];

        if (ppt->has_source_delta_g == _TRUE_) {

          delta_i = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_g]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_g] = delta_i;

          delta_rho_tot += rho_i * delta_i;

          rho_tot += rho_i;

        }

        if (ppt->has_source_theta_g == _TRUE_) {

          theta_i = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_g]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_g] = theta_i;

          rho_plus_p_theta_tot += 4./3. * rho_i * theta_i;

          rho_plus_p_tot += 4./3. * rho_i;

        }

        /* T_b(k,tau) */

        rho_i = pvecback_sp_long[pba->index_bg_rho_b];

        if (ppt->has_source_delta_b == _TRUE_) {

          delta_i = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_b]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_b] = delta_i;

          delta_rho_tot += rho_i * delta_i;

          rho_tot += rho_i;

        }

        if (ppt->has_source_theta_b == _TRUE_) {

          theta_i = ppt->sources[index_md]
            [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_b]
            [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_b] = theta_i;

          rho_plus_p_theta_tot += rho_i * theta_i;

          rho_plus_p_tot += rho_i;

        }

        /* T_cdm(k,tau) */

        if (pba->has_cdm == _TRUE_) {

          rho_i = pvecback_sp_long[pba->index_bg_rho_cdm];

          if (ppt->has_source_delta_cdm == _TRUE_) {

            delta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_cdm]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_cdm] = delta_i;

            delta_rho_tot += rho_i * delta_i;

            rho_tot += rho_i;

          }

          if (ppt->has_source_theta_cdm == _TRUE_) {

            theta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_cdm]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_cdm] = theta_i;

            rho_plus_p_theta_tot += rho_i * theta_i;

            rho_plus_p_tot += rho_i;

          }

        }

        /* T_dcdm(k,tau) */

        if (pba->has_dcdm == _TRUE_) {

          rho_i = pvecback_sp_long[pba->index_bg_rho_dcdm];

          if (ppt->has_source_delta_dcdm == _TRUE_) {

            delta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_dcdm]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_dcdm] = delta_i;

            delta_rho_tot += rho_i * delta_i;

            rho_tot += rho_i;

          }

          if (ppt->has_source_theta_dcdm == _TRUE_) {

            theta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_dcdm]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_dcdm] = theta_i;

            rho_plus_p_theta_tot += rho_i * theta_i;

            rho_plus_p_tot += rho_i;

          }

        }

        /* T_scf(k,tau) */

        if (pba->has_scf == _TRUE_) {

          rho_i = pvecback_sp_long[pba->index_bg_rho_scf];

          if (ppt->has_source_delta_scf == _TRUE_) {

            delta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_scf]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_scf] = delta_i;

            delta_rho_tot += rho_i * delta_i;

            rho_tot += rho_i;
          }

          if (ppt->has_source_theta_scf == _TRUE_) {

            theta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_scf]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_scf] = theta_i;

            rho_plus_p_theta_tot += (rho_i + pvecback_sp_long[pba->index_bg_p_scf]) * theta_i;

            rho_plus_p_tot += (rho_i + pvecback_sp_long[pba->index_bg_p_scf]);
          }

        }


        /* T_fld(k,tau) */

        if (pba->has_fld == _TRUE_) {

          rho_i = pvecback_sp_long[pba->index_bg_rho_fld];

          if (ppt->has_source_delta_fld == _TRUE_) {

            delta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_fld]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_fld] = delta_i;

            delta_rho_tot += rho_i * delta_i;

            rho_tot += rho_i;

          }

          if (ppt->has_source_theta_fld == _TRUE_) {

            theta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_fld]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_fld] = theta_i;

            rho_plus_p_theta_tot += (1. + pba->w0_fld + pba->wa_fld * (1. - pvecback_sp_long[pba->index_bg_a] / pba->a_today)) * rho_i * theta_i;

            rho_plus_p_tot += (1. + pba->w0_fld + pba->wa_fld * (1. - pvecback_sp_long[pba->index_bg_a] / pba->a_today)) * rho_i;

          }

        }

        /* T_ur(k,tau) */

        if (pba->has_ur == _TRUE_) {

          rho_i = pvecback_sp_long[pba->index_bg_rho_ur];

          if (ppt->has_source_delta_ur == _TRUE_) {

            delta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_ur]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_ur] = delta_i;

            delta_rho_tot += rho_i * delta_i;

            rho_tot += rho_i;

          }

          if (ppt->has_source_theta_ur == _TRUE_) {

            theta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_ur]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_ur] = theta_i;

            rho_plus_p_theta_tot += 4./3. * rho_i * theta_i;

            rho_plus_p_tot += 4./3. * rho_i;

          }

        }

        /* T_dr(k,tau) */

        if (pba->has_dr == _TRUE_) {

          rho_i = pvecback_sp_long[pba->index_bg_rho_dr];

          if (ppt->has_source_delta_dr == _TRUE_) {

            delta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_dr]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_dr] = delta_i;

            delta_rho_tot += rho_i * delta_i;

            rho_tot += rho_i;

          }

          if (ppt->has_source_theta_dr == _TRUE_) {

            theta_i = ppt->sources[index_md]
              [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_dr]
              [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

            psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_dr] = theta_i;

            rho_plus_p_theta_tot += 4./3. * rho_i * theta_i;

            rho_plus_p_tot += 4./3. * rho_i;

          }

        }

        /* T_ncdm_i(k,tau) */

        if (pba->has_ncdm == _TRUE_) {

          for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {

            rho_i = pvecback_sp_long[pba->index_bg_rho_ncdm1+n_ncdm];

            if (ppt->has_source_delta_ncdm == _TRUE_) {

              delta_i = ppt->sources[index_md]
                [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_ncdm1+n_ncdm]
                [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

              psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_ncdm1+n_ncdm] = delta_i;

              delta_rho_tot += rho_i * delta_i;

              rho_tot += rho_i;
            }

            if (ppt->has_source_theta_ncdm == _TRUE_) {

              theta_i = ppt->sources[index_md]
                [index_ic * ppt->tp_size[index_md] + ppt->index_tp_theta_ncdm1+n_ncdm]
                [(index_tau-psp->ln_tau_size+ppt->tau_size) * ppt->k_size[index_md] + index_k];

              psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_ncdm1+n_ncdm] = theta_i;

              rho_plus_p_theta_tot += (rho_i + pvecback_sp_long[pba->index_bg_p_ncdm1+n_ncdm]) * theta_i;

              rho_plus_p_tot += (rho_i + pvecback_sp_long[pba->index_bg_p_ncdm1+n_ncdm]);
            }

          }

        }

        /* could include homogeneous component in rho_tot if uncommented (leave commented to match CMBFAST/CAMB definition) */

        /* 	if (pba->has_lambda == _TRUE_) { */

        /* 	  rho_i = pvecback_sp_long[pba->index_bg_rho_lambda]; */

        /* 	  rho_tot += rho_i; */
        /* 	} */

        /* T_tot(k,tau) */

        if (ppt->has_density_transfers == _TRUE_) {

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_delta_tot] = delta_rho_tot/rho_tot;

        }

        if (ppt->has_velocity_transfers == _TRUE_) {

          psp->matter_transfer[((index_tau*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + psp->index_tr_theta_tot] = rho_plus_p_theta_tot/rho_plus_p_tot;

        }

      }
    }
  }

  /**- if interpolation of P(k,tau) will be needed (as a function of tau),
     compute array of second derivatives in view of spline interpolation */

  if (psp->ln_tau_size > 1) {

    class_alloc(psp->ddmatter_transfer,sizeof(double)*psp->ln_tau_size*psp->ln_k_size*psp->ic_size[index_md]*psp->tr_size,psp->error_message);

    class_call(array_spline_table_lines(psp->ln_tau,
                                        psp->ln_tau_size,
                                        psp->matter_transfer,
                                        psp->ic_size[index_md]*psp->ln_k_size*psp->tr_size,
                                        psp->ddmatter_transfer,
                                        _SPLINE_EST_DERIV_,
                                        psp->error_message),
               psp->error_message,
               psp->error_message);

  }

  free (pvecback_sp_long);

  return _SUCCESS_;
}

int spectra_output_tk_titles(struct background *pba,
                             struct perturbs *ppt,
                             enum file_format output_format,
                             char titles[_MAXTITLESTRINGLENGTH_]
                             ){
  int n_ncdm;
  char tmp[40];

  if (output_format == class_format) {
    class_store_columntitle(titles,"k (h/Mpc)",_TRUE_);
    if (ppt->has_density_transfers == _TRUE_) {
      class_store_columntitle(titles,"d_g",_TRUE_);
      class_store_columntitle(titles,"d_b",_TRUE_);
      class_store_columntitle(titles,"d_cdm",pba->has_cdm);
      class_store_columntitle(titles,"d_fld",pba->has_fld);
      class_store_columntitle(titles,"d_ur",pba->has_ur);
      if (pba->has_ncdm == _TRUE_) {
        for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
          sprintf(tmp,"d_ncdm[%d]",n_ncdm);
          class_store_columntitle(titles,tmp,_TRUE_);
        }
      }
      class_store_columntitle(titles,"d_dcdm",pba->has_dcdm);
      class_store_columntitle(titles,"d_dr",pba->has_dr);
      class_store_columntitle(titles,"d_scf",pba->has_scf);
      class_store_columntitle(titles,"d_tot",_TRUE_);
    }
    if (ppt->has_velocity_transfers == _TRUE_) {
      class_store_columntitle(titles,"t_g",_TRUE_);
      class_store_columntitle(titles,"t_b",_TRUE_);
      class_store_columntitle(titles,"t_cdm",((pba->has_cdm == _TRUE_) && (ppt->gauge != synchronous)));
      class_store_columntitle(titles,"t_fld",pba->has_fld);
      class_store_columntitle(titles,"t_ur",pba->has_ur);
      if (pba->has_ncdm == _TRUE_) {
        for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
          sprintf(tmp,"t_ncdm[%d]",n_ncdm);
          class_store_columntitle(titles,tmp,_TRUE_);
        }
      }
      class_store_columntitle(titles,"t_dcdm",pba->has_dcdm);
      class_store_columntitle(titles,"t_dr",pba->has_dr);
      class_store_columntitle(titles,"t__scf",pba->has_scf);
      class_store_columntitle(titles,"t_tot",_TRUE_);
    }
  }

  else if (output_format == camb_format) {

    class_store_columntitle(titles,"k (h/Mpc)",_TRUE_);
    class_store_columntitle(titles,"-T_cdm/k2",_TRUE_);
    class_store_columntitle(titles,"-T_b/k2",_TRUE_);
    class_store_columntitle(titles,"-T_g/k2",_TRUE_);
    class_store_columntitle(titles,"-T_ur/k2",_TRUE_);
    class_store_columntitle(titles,"-T_ncdm/k2",_TRUE_);
    class_store_columntitle(titles,"-T_tot/k2",_TRUE_);

  }

  return _SUCCESS_;

}

int spectra_output_tk_data(
                          struct background * pba,
                          struct perturbs * ppt,
                          struct spectra * psp,
                          enum file_format output_format,
                          double z,
                          int number_of_titles,
                          double *data
                          ) {

  int n_ncdm;
  double k, k_over_h, k2;
  double * tkfull=NULL;  /* array with argument
                   pk_ic[(index_k * psp->ic_size[index_md] + index_ic)*psp->tr_size+index_tr] */
  double *tk;
  double *dataptr;

  int index_md=0;
  int index_ic;
  int index_k;
  int index_tr;
  int storeidx;

  if (psp->ln_k_size*psp->ic_size[index_md]*psp->tr_size > 0){
  class_alloc(tkfull,
              psp->ln_k_size*psp->ic_size[index_md]*psp->tr_size*sizeof(double),
              psp->error_message);
  }

    /** - compute T_i(k) for each k (if several ic's, compute it for each ic; if z_pk = 0, this is done by directly reading inside the pre-computed table; if not, this is done by interpolating the table at the correct value of tau. */

    /* if z_pk = 0, no interpolation needed */

    if (z == 0.) {

      for (index_k=0; index_k<psp->ln_k_size; index_k++) {
        for (index_tr=0; index_tr<psp->tr_size; index_tr++) {
          for (index_ic=0; index_ic<psp->ic_size[index_md]; index_ic++) {
            tkfull[(index_k * psp->ic_size[index_md] + index_ic) * psp->tr_size + index_tr] = psp->matter_transfer[(((psp->ln_tau_size-1)*psp->ln_k_size + index_k) * psp->ic_size[index_md] + index_ic) * psp->tr_size + index_tr];
          }
        }
      }
    }

    /* if 0 <= z_pk <= z_max_pk, interpolation needed, */
    else {

      class_call(spectra_tk_at_z(pba,
                                 psp,
                                 z,
                                 tkfull),
                 psp->error_message,
                 psp->error_message);
    }

    /** - store data */

    for (index_ic = 0; index_ic < psp->ic_size[index_md]; index_ic++) {

      for (index_k=0; index_k<psp->ln_k_size; index_k++) {

        storeidx = 0;
        dataptr = data+index_ic*(psp->ln_k_size*number_of_titles)+index_k*number_of_titles;
        tk = &(tkfull[(index_k * psp->ic_size[index_md] + index_ic) * psp->tr_size]);
        k = exp(psp->ln_k[index_k]);
        k2 = k*k;
        k_over_h = k/pba->h;

        class_store_double(dataptr, k_over_h, _TRUE_,storeidx);

        /* indices for species associated with a velocity transfer function in Fourier space */

        if (output_format == class_format) {

          if (ppt->has_density_transfers == _TRUE_) {

            class_store_double(dataptr,tk[psp->index_tr_delta_g],ppt->has_source_delta_g,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_delta_b],ppt->has_source_delta_b,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_delta_cdm],ppt->has_source_delta_cdm,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_delta_fld],ppt->has_source_delta_fld,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_delta_ur],ppt->has_source_delta_ur,storeidx);
            if (pba->has_ncdm == _TRUE_){
              for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++){
                class_store_double(dataptr,tk[psp->index_tr_delta_ncdm1+n_ncdm],ppt->has_source_delta_ncdm,storeidx);
              }
            }
            class_store_double(dataptr,tk[psp->index_tr_delta_dcdm],ppt->has_source_delta_dcdm,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_delta_dr],ppt->has_source_delta_dr,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_delta_scf],ppt->has_source_delta_scf,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_delta_tot],_TRUE_,storeidx);

          }
          if (ppt->has_velocity_transfers == _TRUE_) {

            class_store_double(dataptr,tk[psp->index_tr_theta_g],ppt->has_source_theta_g,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_theta_b],ppt->has_source_theta_b,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_theta_cdm],ppt->has_source_theta_cdm,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_theta_fld],ppt->has_source_theta_fld,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_theta_ur],ppt->has_source_theta_ur,storeidx);
            if (pba->has_ncdm == _TRUE_){
              for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++){
                class_store_double(dataptr,tk[psp->index_tr_theta_ncdm1+n_ncdm],ppt->has_source_theta_ncdm,storeidx);
              }
            }
            class_store_double(dataptr,tk[psp->index_tr_theta_dcdm],ppt->has_source_theta_dcdm,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_theta_dr],ppt->has_source_theta_dr,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_theta_scf],ppt->has_source_theta_scf,storeidx);
            class_store_double(dataptr,tk[psp->index_tr_theta_tot],_TRUE_,storeidx);

          }

        }
        else if (output_format == camb_format) {

          /* rescale and reorder the matter transfer functions following the CMBFAST/CAMB convention */
          class_store_double_or_default(dataptr,-tk[psp->index_tr_delta_cdm]/k2,ppt->has_source_delta_cdm,storeidx,0.0);
          class_store_double_or_default(dataptr,-tk[psp->index_tr_delta_b]/k2,ppt->has_source_delta_b,storeidx,0.0);
          class_store_double_or_default(dataptr,-tk[psp->index_tr_delta_g]/k2,ppt->has_source_delta_g,storeidx,0.0);
          class_store_double_or_default(dataptr,-tk[psp->index_tr_delta_ur]/k2,ppt->has_source_delta_ur,storeidx,0.0);
          class_store_double_or_default(dataptr,-tk[psp->index_tr_delta_ncdm1]/k2,ppt->has_source_delta_ncdm,storeidx,0.0);
          class_store_double_or_default(dataptr,-tk[psp->index_tr_delta_tot]/k2,_TRUE_,storeidx,0.0);

        }
      }
    }

    //Neccessary because the size could be zero (if psp->tr_size is zero)
    if (tkfull != NULL)
      free(tkfull);

    return _SUCCESS_;
}

 int spectra_firstline_and_ic_suffix(struct perturbs *ppt,
                                    int index_ic,
                                    char first_line[_LINE_LENGTH_MAX_],
                                    FileName ic_suffix){

  first_line[0]='\0';
  ic_suffix[0]='\0';


  if ((ppt->has_ad == _TRUE_) && (index_ic == ppt->index_ic_ad)) {
    strcpy(ic_suffix,"ad");
    strcpy(first_line,"for adiabatic (AD) mode (normalized to initial curvature=1) ");
  }

  if ((ppt->has_bi == _TRUE_) && (index_ic == ppt->index_ic_bi)) {
    strcpy(ic_suffix,"bi");
    strcpy(first_line,"for baryon isocurvature (BI) mode (normalized to initial entropy=1)");
  }

  if ((ppt->has_cdi == _TRUE_) && (index_ic == ppt->index_ic_cdi)) {
    strcpy(ic_suffix,"cdi");
    strcpy(first_line,"for CDM isocurvature (CDI) mode (normalized to initial entropy=1)");
  }

  if ((ppt->has_nid == _TRUE_) && (index_ic == ppt->index_ic_nid)) {
    strcpy(ic_suffix,"nid");
    strcpy(first_line,"for neutrino density isocurvature (NID) mode (normalized to initial entropy=1)");
  }

  if ((ppt->has_niv == _TRUE_) && (index_ic == ppt->index_ic_niv)) {
    strcpy(ic_suffix,"niv");
    strcpy(first_line,"for neutrino velocity isocurvature (NIV) mode (normalized to initial entropy=1)");
  }
  return _SUCCESS_;
}
