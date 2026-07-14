/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Nickolas B. Moccetti, Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
/**
Function to call solution phase Minimization        
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#ifdef USE_MPI
	#include "mpi.h"
#endif

#include "nlopt.h"
#include "MAGEMin.h"
#include "gem_function.h"
#include "dump_function.h"
#include "toolkit.h"
#include "GH_database/GH_gem_function.h"
#include "phase_update_function.h"
#include "all_solution_phases.h"
/** 
Function to update xi and sum_xi during local minimization.
*/
SS_ref SS_UPDATE_function(		global_variable 	 gv,
								SS_ref 				 SS_ref_db, 
								bulk_info 	 		 z_b,
								char    			*name){

	/* sf_ok?*/
	if (strcmp(gv.research_group, "tc") 	== 0 ){
		SS_ref_db.sf_ok = 1;
		for (int i = 0; i < SS_ref_db.n_sf; i++){
			if (SS_ref_db.sf[i] < 0.0 || isnan(SS_ref_db.sf[i]) == 1|| isinf(SS_ref_db.sf[i]) == 1){
				SS_ref_db.sf_ok = 0;
				SS_ref_db.sf_id = i;
				break;
			}
		}

		/* xi calculation (phase fraction expression for PGE) */
		SS_ref_db.sum_xi 	= 0.0;	
		for (int i = 0; i < SS_ref_db.n_em; i++){ 
			SS_ref_db.xi_em[i] = exp(-SS_ref_db.mu[i]/(SS_ref_db.R*SS_ref_db.T));
			SS_ref_db.sum_xi  += SS_ref_db.xi_em[i]*SS_ref_db.p[i]*SS_ref_db.z_em[i];
		}

		/* get composition of solution phase */
		for (int j = 0; j < gv.len_ox; j++){
			SS_ref_db.ss_comp[j] = 0.0;
			for (int i = 0; i < SS_ref_db.n_em; i++){
				SS_ref_db.ss_comp[j] += SS_ref_db.Comp[i][j]*SS_ref_db.p[i]*SS_ref_db.z_em[i];
			} 
		}
	}
	else if (strcmp(gv.research_group, "sb") 	== 0 ){
		SS_ref_db.sf_ok = 1;
	}
	else if (strcmp(gv.research_group, "gh") 	== 0 ){
		/* "gh" is LP-only for now (no PGE/mu[] support yet, same reason "sb"
		   forces LP-only - see SetupDatabase), so sf_ok is unused; keep it
		   satisfied like "sb" rather than the "tc" PGE-oriented check.     */
		SS_ref_db.sf_ok = 1;
	}

	return SS_ref_db;
};


/** 
Function to update xi and sum_xi for the considered phases list (during the inner loop of the PGE stage).
NOTE: When the phase is "liq", the normalization factor is also updated as it depends on the endmember fractions
*/
csd_phase_set CP_UPDATE_function(		global_variable 	gv,
										SS_ref 				SS_ref_db,
										csd_phase_set  		cp, 
										bulk_info 	z_b			){

	/* sf_ok?*/
	cp.sf_ok = 1;
	for (int i = 0; i < cp.n_sf; i++){
		if (cp.sf[i] < 0.0 || isnan(cp.sf[i]) == 1|| isinf(cp.sf[i]) == 1){
			cp.sf_ok = 0;	
			break;
		}
	}
	cp.sum_xi 	= 0.0;	
	for (int i = 0; i < cp.n_em; i++){ 
		cp.xi_em[i] = exp(-cp.mu[i]/(SS_ref_db.R*SS_ref_db.T));
		cp.sum_xi  += cp.xi_em[i]*cp.p_em[i]*SS_ref_db.z_em[i];
	}

	/* get composition of solution phase */
	for (int j = 0; j < gv.len_ox; j++){
		cp.ss_comp[j] = 0.0;
		for (int i = 0; i < cp.n_em; i++){
		   cp.ss_comp[j] += SS_ref_db.Comp[i][j]*cp.p_em[i]*SS_ref_db.z_em[i];
	   } 
	}

	return cp;
};

/** 
	This function ensures that if we drift away from the set of x-eos obtained during levelling, a copy will of the phase will be added to the considered set of phase
	- Drifting occurs when tilting of the hyperplane moves the x-eos far away from their initial guess
	- Note that each instance of the phase, initialized during levelling can only be split once
*/
global_variable split_cp(		global_variable 	 gv,
								SS_ref 			    *SS_ref_db,
								csd_phase_set  		*cp
){
	int id_cp;
	int ph_id;
	double distance;

	for (int i = 0; i < gv.len_cp; i++){ 
		if (cp[i].ss_flags[0] == 1){

			ph_id= cp[i].id;
			
			distance 	= euclidean_distance( cp[i].xeos, cp[i].dguess, SS_ref_db[ph_id].n_xeos);

			if (distance > 2.0*gv.SS_PC_stp[ph_id]*pow((double)SS_ref_db[ph_id].n_xeos,0.5) && cp[i].split == 0){
				id_cp 					= gv.len_cp;
						
				cp[id_cp].split 		= 1;							/* set split number to one */
				cp[i].split 			= 1;							/* set split number to one */

				strcpy(cp[id_cp].name,gv.SS_list[ph_id]);				/* get phase name */	
				
				cp[id_cp].id 			= ph_id;						/* get phaseid */
				cp[id_cp].n_xeos		= SS_ref_db[ph_id].n_xeos;		/* get number of compositional variables */
				cp[id_cp].n_em			= SS_ref_db[ph_id].n_em;		/* get number of endmembers */
				cp[id_cp].n_sf			= SS_ref_db[ph_id].n_sf;		/* get number of site fractions */
				
				cp[id_cp].df			= 0.0;
				cp[id_cp].factor		= 0.0;	
				
				cp[id_cp].ss_flags[0] 	= 1;							/* set flags */
				cp[id_cp].ss_flags[1] 	= 0;
				cp[id_cp].ss_flags[2] 	= 1;
				
				cp[id_cp].ss_n          = 0.0;							/* get initial phase fraction */
				
				for (int ii = 0; ii < SS_ref_db[ph_id].n_em; ii++){
					cp[id_cp].p_em[ii]      = 0.0;
				}
				for (int ii = 0; ii < SS_ref_db[ph_id].n_xeos; ii++){
					cp[id_cp].dguess[ii]    = cp[i].dguess[ii];
					cp[id_cp].xeos[ii]      = cp[i].dguess[ii];
					cp[i].dguess[ii]    	= cp[i].xeos[ii];
				}

				gv.n_solvi[ph_id] 	+= 1;
				gv.len_cp 			+= 1;
				
				if (gv.verbose == 1){
					printf("\n  {FYI} %4s cp#%d is grazing away from its field, a copy has been added (xeos = dguess)\n",gv.SS_list[ph_id],i);
				}
				
				if (gv.len_cp == gv.max_n_cp){
					printf(" !! Maxmimum number of allowed phases under consideration reached !!\n    -> check your problem and potentially increase gv.max_n_cp\n");
				}
				
			}
		}	
	}
	return gv;


};

/**
	copy the minimized phase informations to cp structure, if the site fractions are respected
*/
void copy_to_cp(		int 				 i, 
						int 				 ph_id,
						global_variable 	 gv,
						SS_ref 			    *SS_ref_db,
						csd_phase_set  		*cp					){

	cp[i].min_time			= SS_ref_db[ph_id].LM_time;
	cp[i].df				= SS_ref_db[ph_id].df_raw;
	cp[i].factor			= SS_ref_db[ph_id].factor;
	cp[i].sum_xi			= SS_ref_db[ph_id].sum_xi;

	for (int ii = 0; ii < cp[i].n_xeos; ii++){
		cp[i].xeos_0[ii]	= cp[i].xeos[ii]; 
		cp[i].xeos[ii]		= SS_ref_db[ph_id].iguess[ii]; 
		cp[i].xeos_1[ii]	= SS_ref_db[ph_id].iguess[ii]; 
		cp[i].dfx[ii]		= SS_ref_db[ph_id].dfx[ii]; 
	}
	
	for (int ii = 0; ii < cp[i].n_em; ii++){
		cp[i].p_em[ii]		= SS_ref_db[ph_id].p[ii];
		cp[i].xi_em[ii]		= SS_ref_db[ph_id].xi_em[ii];
		cp[i].mu[ii]		= SS_ref_db[ph_id].mu[ii];
	}
	for (int ii = 0; ii < gv.len_ox; ii++){
		cp[i].ss_comp[ii]	= SS_ref_db[ph_id].ss_comp[ii];
	}
	
	for (int ii = 0; ii < cp[i].n_sf; ii++){
		cp[i].sf[ii]		= SS_ref_db[ph_id].sf[ii];
	}
}



/**
	add minimized phase to LP PGE pseudocompound list 
*/
void copy_to_Ppc(		int 				 pc_check,
						int 				 add,
						int 				 ph_id,
						global_variable 	 gv,

						obj_type 			*SS_objective,
						SS_ref 			    *SS_ref_db					){

		double G;
		int    m_Ppc;
		int    save_mSS;


		if (pc_check != 2 || SS_ref_db[ph_id].df < gv.mSS_df_min_add || SS_ref_db[ph_id].df > gv.mSS_df_max_add){
			save_mSS = 0;
		}
		else{
			save_mSS = 1;
		}

		/* get unrotated gbase */
		SS_ref_db[ph_id] = non_rot_hyperplane(	gv, 
												SS_ref_db[ph_id]			);

		/* get unrotated minimized point informations */
		G 	=		 (*SS_objective[ph_id])(	SS_ref_db[ph_id].n_xeos,
												SS_ref_db[ph_id].iguess,
												NULL,
												&SS_ref_db[ph_id]			);

		/* check where to add the new phase PC */
		if (SS_ref_db[ph_id].id_Ppc >= SS_ref_db[ph_id].n_Ppc){ SS_ref_db[ph_id].id_Ppc = 0; printf("SS_LP, MAXIMUM STORAGE SPACE FOR PC IS REACHED for %4s, INCREASED #PC_MAX\n",gv.SS_list[ph_id]);}
		m_Ppc = SS_ref_db[ph_id].id_Ppc;

		if (save_mSS == 1){
			SS_ref_db[ph_id].info_Ppc[m_Ppc]   = 9;
		}
		else{
			SS_ref_db[ph_id].info_Ppc[m_Ppc]   = 0;
		}

		SS_ref_db[ph_id].DF_Ppc[m_Ppc]     = G;
		
		/* get pseudocompound composition */
		for (int j = 0; j < gv.len_ox; j++){				
			SS_ref_db[ph_id].comp_Ppc[m_Ppc][j] = SS_ref_db[ph_id].ss_comp[j]*SS_ref_db[ph_id].factor;	/** composition */
		}
		for (int j = 0; j < SS_ref_db[ph_id].n_em; j++){												/** save coordinates */
			SS_ref_db[ph_id].p_Ppc[m_Ppc][j]  = SS_ref_db[ph_id].p[j];												
			SS_ref_db[ph_id].mu_Ppc[m_Ppc][j] = SS_ref_db[ph_id].mu[j]*SS_ref_db[ph_id].z_em[j];										
		}
		/* save xeos */
		for (int j = 0; j < SS_ref_db[ph_id].n_xeos; j++){		
			SS_ref_db[ph_id].xeos_Ppc[m_Ppc][j] = SS_ref_db[ph_id].iguess[j];							/** compositional variables */
		}	
		SS_ref_db[ph_id].G_Ppc[m_Ppc] = G;
		
		/* add increment to the number of considered phases */
		if (SS_ref_db[ph_id].tot_Ppc < SS_ref_db[ph_id].n_Ppc){
			SS_ref_db[ph_id].tot_Ppc += 1;
		}
		else{
			SS_ref_db[ph_id].tot_Ppc = SS_ref_db[ph_id].n_Ppc;
		}
		SS_ref_db[ph_id].id_Ppc  += 1;
}

/** 
	Minimization function for PGE 
*/
void ss_min_PGE(		global_variable 	 gv,
						PC_type             *PC_read,

						obj_type 			*SS_objective,
						NLopt_type			*NLopt_opt,
						bulk_info 	 		 z_b,
						SS_ref 			    *SS_ref_db,
						csd_phase_set  		*cp
){
	int 	ph_id;
	int 	pc_check;
	clock_t u;

	for (int i = 0; i < gv.len_cp; i++){ 
		if (cp[i].ss_flags[0] == 1){
			pc_check = gv.PC_checked;
			ph_id = cp[i].id;
			cp[i].min_time		  		= 0.0;								/** reset local minimization time to 0.0 */
			u = clock(); 
			/**
				set the iguess of the solution phase to the one of the considered phase 
			*/
			for (int k = 0; k < cp[i].n_xeos; k++) {
				SS_ref_db[ph_id].iguess[k] = cp[i].xeos[k];
			}

			/**
				Rotate G-base hyperplane
			*/
			SS_ref_db[ph_id] = rotate_hyperplane(	gv, 
													SS_ref_db[ph_id]			);

			/**
				Define a sub-hypervolume for the solution phases bounds
			*/
			SS_ref_db[ph_id] = restrict_SS_HyperVolume(	gv, 
														SS_ref_db[ph_id],
														gv.box_size_mode_PGE	);
			
			/**
				call to NLopt for non-linear + inequality constraints optimization
			*/
			SS_ref_db[ph_id] = (*NLopt_opt[ph_id])(		gv,
														SS_ref_db[ph_id]		);										
			
			/**
				establish a set of conditions to update initial guess for next round of local minimization 
			*/
			for (int k = 0; k < cp[i].n_xeos; k++) {
				SS_ref_db[ph_id].iguess[k]   =  SS_ref_db[ph_id].xeos[k];
			}
			
			
			SS_ref_db[ph_id] = PC_function(				gv,
														PC_read,
														SS_ref_db[ph_id], 
														z_b,
														ph_id 					);
													
			SS_ref_db[ph_id] = SS_UPDATE_function(		gv, 
														SS_ref_db[ph_id], 
														z_b, 
														gv.SS_list[ph_id]		);

			u = clock() - u;
			SS_ref_db[ph_id].LM_time = ((double)u)/CLOCKS_PER_SEC*1000.0; 

			;
			/** 
				print solution phase informations (print has to occur before saving PC)
			*/
			if (gv.verbose == 1){
				print_SS_informations(  				gv,
														SS_ref_db[ph_id],
														ph_id					);
			}


			/* if site fractions are respected then save the minimized point */
			if (SS_ref_db[ph_id].sf_ok == 1){
				/**
					copy the minimized phase informations to cp structure
				*/
				copy_to_cp(								i, 
														ph_id,
														gv,
														SS_ref_db,
														cp						);	

				// here we need to save the pseudocompound to have an estimate of the LP Matrix
				if (pc_check == 1){
					copy_to_Ppc(							pc_check,
															0,
															ph_id,
															gv,

															SS_objective,
															SS_ref_db						);
				}					
			}
			else{
				if (gv.verbose == 1){
					printf(" !> SF [:%d] not respected for %4s (SS not updated)\n",SS_ref_db[ph_id].sf_id,gv.SS_list[ph_id]);
				}											
			}
		}
	}

};


/** 
	Minimization function for PGE 
*/
void init_PGE_from_LP(	global_variable 	 gv,
						PC_type				*PC_read,

						obj_type 			*SS_objective,
						bulk_info 	 		 z_b,
						SS_ref 			    *SS_ref_db,
						csd_phase_set  		*cp
){
	int 	ph_id;

	for (int i = 0; i < gv.len_cp; i++){ 
		if (cp[i].ss_flags[0] == 1){
			ph_id = cp[i].id;

			/**
				set the iguess of the solution phase to the one of the considered phase 
			*/
			for (int k = 0; k < cp[i].n_xeos; k++) {
				SS_ref_db[ph_id].iguess[k] = cp[i].xeos[k];
			}

			/**
				Rotate G-base hyperplane
			*/
			SS_ref_db[ph_id] = rotate_hyperplane(	gv, 
													SS_ref_db[ph_id]			);

			
			SS_ref_db[ph_id] = PC_function(				gv,
														PC_read,
														SS_ref_db[ph_id], 
														z_b,
														ph_id 		);
													
			SS_ref_db[ph_id] = SS_UPDATE_function(		gv, 
														SS_ref_db[ph_id], 
														z_b, 
														gv.SS_list[ph_id]		);

			copy_to_cp(									i, 
														ph_id,
														gv,
														SS_ref_db,
														cp						);	

		}
	}

};

/* ====================================================================
   "liq" redundant-occurrence pseudocompound synthesis (gh and tc)
   (see docs/liq_pseudocompound_shortcut.md for the full derivation)

   When many simultaneous liq occurrences (N = count of active cp[] entries
   with id == liq's phase index) share the same rotated hyperplane, they
   routinely converge to the same minimum - every NLopt solve after the
   first is redundant. liq is minimized once in ss_min_LP's main loop
   (gated below), a refined Gamma is fit by weighted least squares across
   liq's own endmembers, and 2*N synthetic pseudocompounds are added to
   the Ppc candidate pool directly on (a linearization of) that refined
   Gamma hyperplane, evenly spread (antipodal pairs) around the single
   real minimum - all allocation-free (fixed-size stack scratch only).

   Two research groups are supported, each via its own eval/synth pair,
   because they differ in ways that go beyond a simple dispatch:
   - GH_liq_eval_raw / GH_liq_pc_synth: gh's direct p=xeos phases (raw mu
     hand-reconstructed from mu_Gex+gb_lvl+entropy term, since gh never
     populates .mu[]; g_eff lives directly in xeos-space; a Sigma(x)=1
     tangent-projection is required, since gh has that equality constraint).
   - TC_eval_raw / TC_liq_pc_synth: tc's reduced-basis phases (.mu[]/
     .df_raw/.p[] are already the raw quantities, no reconstruction needed;
     but n_xeos != n_em with p(x) genuinely nonlinear, so the endmember-
     space residual has to be chain-ruled into xeos-space via .dp_dx[][],
     populated by the objective function itself when called with grad!=NULL;
     no tangent-projection, since tc uses box bounds only, no equality
     constraint).

   gv.liq_pc_synth_active is a global on/off switch (default 1) falling
   back to the legacy per-occurrence NLopt path when set to 0.
   ==================================================================== */
#define LIQ_PC_SYNTH_MAX_DIM 16

/**
    Target xeos step size for the synthetic pseudocompound spread, scaled by
    how far the PGE/LP cycle currently is from convergence: sqrt of the most
    recently recorded Gamma-update norm, gv.gamma_norm[gv.global_ite-1] (the
    same per-iteration convergence diagnostic already used elsewhere, e.g.
    PGE_function.c's own solver-switching logic) - large when Gamma is still
    moving a lot (early iterations, coarser exploration is fine/desirable),
    shrinking automatically as the cycle approaches convergence (finer,
    more localized synthetic points once Gamma itself is barely changing).
    Clamped to [1e-6, 1e-2] so it never vanishes or explodes regardless of
    gv.gamma_norm's own scale. gv.global_ite==0 (no recorded norm yet, e.g.
    the very first LP cycle of a point) falls back to gamma_norm=0, i.e.
    the clamped floor 1e-6.
*/
static double GH_liq_pc_synth_step(	global_variable 	 gv			){
	int    gi    = gv.global_ite - 1;
	double gnorm = (gi >= 0) ? gv.gamma_norm[gi] : 0.0;
	double h     = gv.gh_liq_pc_synth_h * pow(gnorm,1.0/2.0);
	if (h > 1e-2){ h = 1e-2; }
	if (h < 1e-6){ h = 1e-6; }
	return h;
}

/**
    Evaluate a phase's real Gibbs energy at a given point (unrotated: Gamma=0),
    with no NLopt search, and reconstruct its RAW (pre-"boil-down") xeos-
    gradient mu_raw_out[j] = mu_Gex[j] + gb_lvl[j] + R*T*(log(p[j]+d_em[j])+1),
    matching obj_gh_liq's own dSi[] term exactly. Returns df_raw, not df.

    obj_gh_liq (like every gh phase) returns df = df_raw(x) * factor(x), with
    factor(x) = fbc/sum_apep(x) an "oxide boil-down" renormalization that is
    itself x-dependent - so NLopt's own returned gradient (d(df)/dx, via the
    grad[] parameter) is contaminated by factor's own derivative and is NOT
    simply the endmember chemical potential entering the Gamma-hyperplane
    condition. Since df(x) = df_raw(x)*factor(x) and factor(x) > 0 always,
    df(x)=0 and df_raw(x)=0 have the IDENTICAL zero-crossing set, so working
    with the raw, un-renormalized df_raw/mu_raw throughout (as this function
    does) sidesteps the factor(x) chain-rule complexity entirely rather than
    having to differentiate through it.

    LIQ-SPECIFIC: the dSi formula reconstructed here is obj_gh_liq's own
    (single-site, multiplicity-1) configurational entropy term - other gh
    phases with different entropy formulas (e.g. obj_gh_ol's factor-of-2
    two-site term) would need their own version of this reconstruction.
*/
static double GH_liq_eval_raw(		global_variable 	 gv,
									obj_type 			*SS_objective,
									SS_ref 			    *SS_ref_db,
									int 				 ph_id,
									double 				*xeos_in,
									double 				*mu_raw_out		){

	SS_ref_db[ph_id] = non_rot_hyperplane(gv, SS_ref_db[ph_id]);
	for (int k = 0; k < SS_ref_db[ph_id].n_xeos; k++){
		SS_ref_db[ph_id].iguess[k] = xeos_in[k];
	}
	(*SS_objective[ph_id])(	SS_ref_db[ph_id].n_xeos,
							SS_ref_db[ph_id].iguess,
							NULL,
							&SS_ref_db[ph_id]			);

	double R = SS_ref_db[ph_id].R;
	double T = SS_ref_db[ph_id].T;
	for (int j = 0; j < SS_ref_db[ph_id].n_em; j++){
		double dSi_j = R * T * (log(SS_ref_db[ph_id].p[j] + SS_ref_db[ph_id].d_em[j]) + 1.0);
		mu_raw_out[j] = SS_ref_db[ph_id].mu_Gex[j] + SS_ref_db[ph_id].gb_lvl[j] + dSi_j;
	}
	return SS_ref_db[ph_id].df_raw;
}

/**
    Accumulate one phase's weighted contribution to the stage-B normal
    equations MtM*Gamma = Mtmu (n_ox x n_ox, allocation-free, M itself
    never formed) using its real minimum's already-reconstructed raw
    endmember chemical potentials g_star[] (see GH_liq_eval_raw - NOT
    NLopt's own factor-scaled gradient) and its fixed Comp[][] table.
    Endmember rows are weighted by molar abundance p_j*z_em[j] so vanishing
    endmembers (mu dominated by a blowing-up log term) don't destabilize
    the fit, and z_em[j]==0 (deactivated) endmembers drop out entirely.

    Restricted to z_b.nzEl_array[] (the oxides actually present in the
    bulk, z_b.nzEl_val of them) rather than all gv.len_ox - matching the
    same restriction the rest of the codebase's own levelling/simplex
    machinery already applies (gv.gam_tot itself is only ever written at
    nzEl_array positions elsewhere, e.g. simplex_levelling.c). An oxide
    absent from the bulk has no meaningful Gamma to fit in the first place;
    excluding it here (rather than solving for it and discovering the
    system is singular) is both more principled and cheaper.
*/
static void GH_liq_gamma_ls_accumulate(	bulk_info 			 z_b,
											SS_ref 				 SS_ref_db_ph,
											double 				*g_star,
											double 				 MtM[n_ox_all][n_ox_all],
											double 				 Mtmu[n_ox_all]			){

	for (int j = 0; j < SS_ref_db_ph.n_em; j++){
		double w = SS_ref_db_ph.p[j] * SS_ref_db_ph.z_em[j];
		// if (w <= 0.0){ continue; }
		for (int ka = 0; ka < z_b.nzEl_val; ka++){
			int a = z_b.nzEl_array[ka];
			Mtmu[a] += w * SS_ref_db_ph.Comp[j][a] * g_star[j];
			for (int kb = 0; kb < z_b.nzEl_val; kb++){
				int b = z_b.nzEl_array[kb];
				MtM[a][b] += w * SS_ref_db_ph.Comp[j][a] * SS_ref_db_ph.Comp[j][b];
			}
		}
	}
}

/**
    Solve the weighted normal equations for Gamma_new by in-place Gauss
    elimination with partial pivoting - fixed-size, no allocation.

    Restricted to z_b.nzEl_array[] (see GH_liq_gamma_ls_accumulate): builds
    and solves only the nzEl_val x nzEl_val reduced system over oxides
    actually present in the bulk, using compact local indices 0..nzEl_val-1
    internally, then scatters the solved values back into the full-size
    gamma_new[] via nzEl_array. Every oxide NOT in the bulk is left exactly
    at gamma_prior (gamma_new is pre-filled with it) - the rest of the
    codebase never touches gv.gam_tot at those positions either.

    Tikhonov (ridge) regularized on top of that restriction, anchored at
    gamma_prior rather than at zero - a light numerical safety net for
    near-degeneracies *within* the present-oxide subspace (e.g. two oxides
    that happen to be almost perfectly correlated across every one of
    liq's endmembers), not the primary mechanism for handling absent
    oxides (the nzEl_array restriction already does that exactly). lambda
    is scaled to the reduced MtM's own average diagonal magnitude. Still
    returns 0 (should essentially never happen) if the caller should fall
    back to gamma_prior entirely.
*/
static int GH_liq_gamma_solve(		bulk_info 			 z_b,
									double 				 MtM[n_ox_all][n_ox_all],
									double 				*Mtmu,
									double 				*gamma_prior,
									double 				*gamma_new				){

	int    n = z_b.nzEl_val;
	double A[n_ox_all][n_ox_all];
	double b[n_ox_all];
	double x[n_ox_all];

	double trace = 0.0;
	for (int ka = 0; ka < n; ka++){ trace += MtM[z_b.nzEl_array[ka]][z_b.nzEl_array[ka]]; }
	double lambda = 1e-8 * (trace / (double)n);
	if (lambda <= 0.0){ lambda = 1e-8; }	/* MtM all-zero (degenerate) safety floor */

	for (int ka = 0; ka < n; ka++){
		int a = z_b.nzEl_array[ka];
		b[ka] = Mtmu[a] + lambda * gamma_prior[a];
		for (int kc = 0; kc < n; kc++){
			int c = z_b.nzEl_array[kc];
			A[ka][kc] = MtM[a][c] + ((ka == kc) ? lambda : 0.0);
		}
	}

	for (int col = 0; col < n; col++){
		int    piv  = col;
		double best = fabs(A[col][col]);
		for (int r = col+1; r < n; r++){
			if (fabs(A[r][col]) > best){ best = fabs(A[r][col]); piv = r; }
		}
		if (best < 1e-10){ return 0; }
		if (piv != col){
			for (int c = 0; c < n; c++){ double t = A[col][c]; A[col][c] = A[piv][c]; A[piv][c] = t; }
			double t = b[col]; b[col] = b[piv]; b[piv] = t;
		}
		for (int r = col+1; r < n; r++){
			double f = A[r][col] / A[col][col];
			for (int c = col; c < n; c++){ A[r][c] -= f * A[col][c]; }
			b[r] -= f * b[col];
		}
	}
	for (int r = n-1; r >= 0; r--){
		double s = b[r];
		for (int c = r+1; c < n; c++){ s -= A[r][c] * x[c]; }
		x[r] = s / A[r][r];
	}
	for (int r = 0; r < n; r++){ gamma_new[z_b.nzEl_array[r]] = x[r]; }
	return 1;
}

/**
    Stage C: for "liq" only, when its active occurrence count N (cp[] entries
    with id == ph_id) reaches gv.gh_liq_pc_synth_threshold, generate 2*N
    synthetic pseudocompounds lying on (a linearization of) the refined
    Gamma_new hyperplane, spread as antipodal pairs around the single real
    minimum x* found earlier this cycle - added to the Ppc pool via the
    existing PC_function + SS_UPDATE_function + copy_to_Ppc sequence (the
    same 3-call pattern this file's own gv.n_ss_ph[ph_id]==1 block uses,
    generalized here from a fixed 2-point shift-blend to 2*N hyperplane
    points). See docs/liq_pseudocompound_shortcut.md for the full derivation.
*/
static void GH_liq_pc_synth(		global_variable 	 gv,
									PC_type				*PC_read,
									obj_type 			*SS_objective,
									bulk_info 	 		 z_b,
									SS_ref 			    *SS_ref_db,
									csd_phase_set  		*cp,
									int 				 ph_id,
									double 				*gamma_new			){

	int n_xeos = SS_ref_db[ph_id].n_xeos;
	if (n_xeos > LIQ_PC_SYNTH_MAX_DIM){ return; }	/* safety: never index past the fixed-size scratch below */

	double x_star[LIQ_PC_SYNTH_MAX_DIM];
	double g_star[LIQ_PC_SYNTH_MAX_DIM];
	double g_eff[LIQ_PC_SYNTH_MAX_DIM];
	double x_mean[LIQ_PC_SYNTH_MAX_DIM];
	double x_base[LIQ_PC_SYNTH_MAX_DIM];
	double e_k[LIQ_PC_SYNTH_MAX_DIM];
	double v_k[LIQ_PC_SYNTH_MAX_DIM];
	double x_plus[LIQ_PC_SYNTH_MAX_DIM];
	double x_minus[LIQ_PC_SYNTH_MAX_DIM];

	for (int k = 0; k < n_xeos; k++){ x_star[k] = SS_ref_db[ph_id].xeos[k]; }

	/* G*, g* at x* (single direct, no-NLopt evaluation, raw/un-"boiled-down"
	   quantities - see GH_liq_eval_raw for why raw rather than df/grad[]) */
	double G_star = GH_liq_eval_raw(gv, SS_objective, SS_ref_db, ph_id, x_star, g_star);

	/* g_eff = g* - C^T.Gamma_new (raw, z_em-weighted Comp),  dG0 = G* - Gamma_new.Comp_raw(x*).
	   Restricted to z_b.nzEl_array[] (oxides present in the bulk) - see
	   GH_liq_gamma_ls_accumulate's own header comment for why: an endmember
	   can carry a nonzero Comp[][] entry for an oxide the CURRENT bulk
	   doesn't have (e.g. a water-bearing endmember's table entry, even
	   when z_em deactivates it for this bulk), so summing gamma_new over
	   every oxide rather than just the present ones could pull in a stale/
	   prior value through a nonzero Comp[][] entry that should not matter. */
	double dG0 = G_star;
	for (int k = 0; k < n_xeos; k++){ g_eff[k] = g_star[k]; }
	for (int ja = 0; ja < z_b.nzEl_val; ja++){
		int j = z_b.nzEl_array[ja];
		double gam = gamma_new[j];
		for (int k = 0; k < n_xeos; k++){
			double c   = SS_ref_db[ph_id].Comp[k][j] * SS_ref_db[ph_id].z_em[k];
			dG0       -= c * gam * x_star[k];
			g_eff[k]  -= c * gam;
		}
	}

	/* project g_eff onto the simplex's tangent space (vectors summing to 0):
	   every xeos in cp[] satisfies sum_k x_k = 1, so x_star + delta stays a
	   valid composition only if sum_k delta_k = 0 - this holds for e_k
	   (difference of two sum-to-1 vectors) but is NOT automatic for g_eff,
	   so it must be enforced here before g_eff drives delta_0/P below. */
	{
		double mean_g = 0.0;
		for (int k = 0; k < n_xeos; k++){ mean_g += g_eff[k]; }
		mean_g /= (double)n_xeos;
		for (int k = 0; k < n_xeos; k++){ g_eff[k] -= mean_g; }
	}

	double g_eff_norm2 = 0.0;
	for (int k = 0; k < n_xeos; k++){ g_eff_norm2 += g_eff[k]*g_eff[k]; }

	if (g_eff_norm2 > 1e-24){
		double lambda = -dG0 / g_eff_norm2;
		for (int k = 0; k < n_xeos; k++){ x_base[k] = x_star[k] + lambda * g_eff[k]; }
	}
	else{
		for (int k = 0; k < n_xeos; k++){ x_base[k] = x_star[k]; }
	}

	/* pass 1: mean of the N previous occurrences' xeos (cp[i].xeos is never
	   mutated inside ss_min_LP, so it is safe to read here unchanged) */
	int N = 0;
	for (int k = 0; k < n_xeos; k++){ x_mean[k] = 0.0; }
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[0] == 1 && cp[i].id == ph_id){
			for (int k = 0; k < n_xeos; k++){ x_mean[k] += cp[i].xeos[k]; }
			N += 1;
		}
	}
	if (N < 2){ return; }
	for (int k = 0; k < n_xeos; k++){ x_mean[k] /= (double)N; }

	/* pass 2: max norm of the projected deviations P*e_k, to normalize the step size */
	double max_norm = 0.0;
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[0] == 1 && cp[i].id == ph_id){
			double dot = 0.0;
			for (int k = 0; k < n_xeos; k++){
				e_k[k] = cp[i].xeos[k] - x_mean[k];
				dot   += e_k[k] * g_eff[k];
			}
			double f  = (g_eff_norm2 > 1e-24) ? (dot / g_eff_norm2) : 0.0;
			double n2 = 0.0;
			for (int k = 0; k < n_xeos; k++){
				double v = e_k[k] - f * g_eff[k];
				n2 += v*v;
			}
			double nrm = sqrt(n2);
			if (nrm > max_norm){ max_norm = nrm; }
		}
	}
	if (max_norm < 1e-12){ return; }	/* degenerate previous spread - nothing to reuse this cycle */

	double s = GH_liq_pc_synth_step(gv) / max_norm;

	/* pass 3: generate + add the 2*N antipodal points, one pair per real
	   occurrence's own deviation direction from the group mean */
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[0] == 1 && cp[i].id == ph_id){
			double dot = 0.0;
			for (int k = 0; k < n_xeos; k++){
				e_k[k] = cp[i].xeos[k] - x_mean[k];
				dot   += e_k[k] * g_eff[k];
			}
			double f = (g_eff_norm2 > 1e-24) ? (dot / g_eff_norm2) : 0.0;
			int ok_plus = 1, ok_minus = 1;
			for (int k = 0; k < n_xeos; k++){
				v_k[k]     = e_k[k] - f * g_eff[k];
				x_plus[k]  = x_base[k] + s * v_k[k];
				x_minus[k] = x_base[k] - s * v_k[k];
				if (x_plus[k]  < SS_ref_db[ph_id].bounds[k][0] || x_plus[k]  > SS_ref_db[ph_id].bounds[k][1]){ ok_plus  = 0; }
				if (x_minus[k] < SS_ref_db[ph_id].bounds[k][0] || x_minus[k] > SS_ref_db[ph_id].bounds[k][1]){ ok_minus = 0; }
			}
			if (!ok_plus || !ok_minus){ continue; }	/* keep pairs symmetric: drop both rather than clamp asymmetrically */

			for (int pass = 0; pass < 2; pass++){
				double *x_syn = (pass == 0) ? x_plus : x_minus;
				for (int k = 0; k < n_xeos; k++){ SS_ref_db[ph_id].iguess[k] = x_syn[k]; }

				SS_ref_db[ph_id] = PC_function(			gv,
															PC_read,
															SS_ref_db[ph_id],
															z_b,
															ph_id					);

				SS_ref_db[ph_id] = SS_UPDATE_function(		gv,
															SS_ref_db[ph_id],
															z_b,
															gv.SS_list[ph_id]		);

				if (SS_ref_db[ph_id].sf_ok == 1){
					copy_to_Ppc(							0,
															1,
															ph_id,
															gv,

															SS_objective,
															SS_ref_db				);
				}
			}
		}
	}
}

/**
    tc analogue of GH_liq_eval_raw - genuinely phase-generic (not just
    liq): unlike gh, EVERY tc phase's objective function already populates
    d->mu[j] (raw per-endmember chemical potential) and d->df_raw directly
    as a side effect of every call - no hand-reconstruction of a phase-
    specific entropy term needed (see obj_ig_liq: mu[j] = R*T*log(...) +
    gb[j] + mu_Gex[j], already stored; SS_UPDATE_function's own tc branch
    already relies on this generically for xi_em, not liq-specifically).
    What tc needs instead is d->dp_dx[j][i] = dp_j/dx_i (tc's xeos are a
    reduced, generally nonlinear basis, n_xeos != n_em) - this is only
    populated when the objective is called with grad != NULL (e.g.
    obj_ig_liq calls dpdx_ig_liq(...) inside "if (grad){...}"), so grad
    must be requested here even though its own contents are unused (the
    factor-scaled dfx[] it produces is not what this needs - see
    TC_liq_pc_synth for how the raw mu[]/dp_dx[][] are combined instead).
    Callable for any already-minimized tc phase, not just liq - see the
    "other minimized phases" block in ss_min_LP.
*/
static double TC_eval_raw(		global_variable 	 gv,
								obj_type 			*SS_objective,
								SS_ref 			    *SS_ref_db,
								int 				 ph_id,
								double 				*xeos_in,
								double 				*mu_raw_out		){

	SS_ref_db[ph_id] = non_rot_hyperplane(gv, SS_ref_db[ph_id]);
	for (int k = 0; k < SS_ref_db[ph_id].n_xeos; k++){
		SS_ref_db[ph_id].iguess[k] = xeos_in[k];
	}
	double dummy_grad[LIQ_PC_SYNTH_MAX_DIM];
	(*SS_objective[ph_id])(	SS_ref_db[ph_id].n_xeos,
							SS_ref_db[ph_id].iguess,
							dummy_grad,
							&SS_ref_db[ph_id]			);

	for (int j = 0; j < SS_ref_db[ph_id].n_em; j++){
		mu_raw_out[j] = SS_ref_db[ph_id].mu[j];
	}
	return SS_ref_db[ph_id].df_raw;
}

/**
    tc analogue of GH_liq_pc_synth - same occurrence-direction, antipodal-
    pair construction (see that function's own comments for the shared
    parts), but adapted for tc's reduced-basis phases:
    - the driving-force residual is built in ENDMEMBER space first
      (r[j] = mu_raw[j] - z_em[j]*(Comp[j].Gamma_new), j=1..n_em), then
      chain-ruled into XEOS space via dp_dx: g_eff[i] = sum_j dp_dx[j][i]*r[j]
      (dp_dx populated by TC_eval_raw's grad-requesting call);
    - no Sigma(x)=1 tangent-projection: tc's reduced-basis xeos have only
      box bounds, no such equality constraint (unlike gh).
*/
static void TC_liq_pc_synth(		global_variable 	 gv,
									PC_type				*PC_read,
									obj_type 			*SS_objective,
									bulk_info 	 		 z_b,
									SS_ref 			    *SS_ref_db,
									csd_phase_set  		*cp,
									int 				 ph_id,
									double 				*gamma_new			){

	int n_xeos = SS_ref_db[ph_id].n_xeos;
	int n_em   = SS_ref_db[ph_id].n_em;
	if (n_xeos > LIQ_PC_SYNTH_MAX_DIM || n_em > LIQ_PC_SYNTH_MAX_DIM){ return; }	/* safety: never index past the fixed-size scratch below */

	double x_star[LIQ_PC_SYNTH_MAX_DIM];
	double mu_raw[LIQ_PC_SYNTH_MAX_DIM];		/* endmember-space */
	double r_em[LIQ_PC_SYNTH_MAX_DIM];			/* endmember-space residual */
	double g_eff[LIQ_PC_SYNTH_MAX_DIM];		/* xeos-space, via dp_dx chain rule */
	double x_mean[LIQ_PC_SYNTH_MAX_DIM];
	double x_base[LIQ_PC_SYNTH_MAX_DIM];
	double e_k[LIQ_PC_SYNTH_MAX_DIM];
	double v_k[LIQ_PC_SYNTH_MAX_DIM];
	double x_plus[LIQ_PC_SYNTH_MAX_DIM];
	double x_minus[LIQ_PC_SYNTH_MAX_DIM];

	for (int k = 0; k < n_xeos; k++){ x_star[k] = SS_ref_db[ph_id].xeos[k]; }

	/* G*, raw endmember mu* at x* (grad-requesting call so dp_dx[][] gets populated too) */
	double G_star = TC_eval_raw(gv, SS_objective, SS_ref_db, ph_id, x_star, mu_raw);

	/* r[j] = mu_raw[j] - z_em[j]*(Comp[j].Gamma_new),  dG0 = G* - sum_j p*[j]*z_em[j]*(Comp[j].Gamma_new).
	   Restricted to z_b.nzEl_array[] (oxides present in the bulk) - see
	   GH_liq_gamma_ls_accumulate's header comment for why. */
	double dG0 = G_star;
	for (int j = 0; j < n_em; j++){
		double cg = 0.0;
		for (int ka = 0; ka < z_b.nzEl_val; ka++){
			int a = z_b.nzEl_array[ka];
			cg += SS_ref_db[ph_id].Comp[j][a] * gamma_new[a];
		}
		cg *= SS_ref_db[ph_id].z_em[j];
		r_em[j] = mu_raw[j] - cg;
		dG0    -= SS_ref_db[ph_id].p[j] * cg;
	}

	/* chain rule into xeos-space: g_eff[i] = sum_j dp_dx[j][i]*r_em[j] (no
	   tangent-projection - see function header comment) */
	for (int k = 0; k < n_xeos; k++){
		double sum = 0.0;
		for (int j = 0; j < n_em; j++){ sum += SS_ref_db[ph_id].dp_dx[j][k] * r_em[j]; }
		g_eff[k] = sum;
	}

	double g_eff_norm2 = 0.0;
	for (int k = 0; k < n_xeos; k++){ g_eff_norm2 += g_eff[k]*g_eff[k]; }

	if (g_eff_norm2 > 1e-24){
		double lambda = -dG0 / g_eff_norm2;
		for (int k = 0; k < n_xeos; k++){ x_base[k] = x_star[k] + lambda * g_eff[k]; }
	}
	else{
		for (int k = 0; k < n_xeos; k++){ x_base[k] = x_star[k]; }
	}

	/* pass 1: mean of the N previous occurrences' xeos */
	int N = 0;
	for (int k = 0; k < n_xeos; k++){ x_mean[k] = 0.0; }
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[0] == 1 && cp[i].id == ph_id){
			for (int k = 0; k < n_xeos; k++){ x_mean[k] += cp[i].xeos[k]; }
			N += 1;
		}
	}
	if (N < 2){ return; }
	for (int k = 0; k < n_xeos; k++){ x_mean[k] /= (double)N; }

	/* pass 2: max norm of the projected deviations P*e_k, to normalize the step size */
	double max_norm = 0.0;
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[0] == 1 && cp[i].id == ph_id){
			double dot = 0.0;
			for (int k = 0; k < n_xeos; k++){
				e_k[k] = cp[i].xeos[k] - x_mean[k];
				dot   += e_k[k] * g_eff[k];
			}
			double f  = (g_eff_norm2 > 1e-24) ? (dot / g_eff_norm2) : 0.0;
			double n2 = 0.0;
			for (int k = 0; k < n_xeos; k++){
				double v = e_k[k] - f * g_eff[k];
				n2 += v*v;
			}
			double nrm = sqrt(n2);
			if (nrm > max_norm){ max_norm = nrm; }
		}
	}
	if (max_norm < 1e-12){ return; }	/* degenerate previous spread - nothing to reuse this cycle */

	double s = GH_liq_pc_synth_step(gv) / max_norm;

	/* pass 3: generate + add the 2*N antipodal points */
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[0] == 1 && cp[i].id == ph_id){
			double dot = 0.0;
			for (int k = 0; k < n_xeos; k++){
				e_k[k] = cp[i].xeos[k] - x_mean[k];
				dot   += e_k[k] * g_eff[k];
			}
			double f = (g_eff_norm2 > 1e-24) ? (dot / g_eff_norm2) : 0.0;
			int ok_plus = 1, ok_minus = 1;
			for (int k = 0; k < n_xeos; k++){
				v_k[k]     = e_k[k] - f * g_eff[k];
				x_plus[k]  = x_base[k] + s * v_k[k];
				x_minus[k] = x_base[k] - s * v_k[k];
				if (x_plus[k]  < SS_ref_db[ph_id].bounds[k][0] || x_plus[k]  > SS_ref_db[ph_id].bounds[k][1]){ ok_plus  = 0; }
				if (x_minus[k] < SS_ref_db[ph_id].bounds[k][0] || x_minus[k] > SS_ref_db[ph_id].bounds[k][1]){ ok_minus = 0; }
			}
			if (!ok_plus || !ok_minus){ continue; }	/* keep pairs symmetric: drop both rather than clamp asymmetrically */

			for (int pass = 0; pass < 2; pass++){
				double *x_syn = (pass == 0) ? x_plus : x_minus;
				for (int k = 0; k < n_xeos; k++){ SS_ref_db[ph_id].iguess[k] = x_syn[k]; }

				SS_ref_db[ph_id] = PC_function(				gv,
															PC_read,
															SS_ref_db[ph_id],
															z_b,
															ph_id					);

				SS_ref_db[ph_id] = SS_UPDATE_function(		gv,
															SS_ref_db[ph_id],
															z_b,
															gv.SS_list[ph_id]		);

				if (SS_ref_db[ph_id].sf_ok == 1){
					copy_to_Ppc(							0,
															1,
															ph_id,
															gv,

															SS_objective,
															SS_ref_db				);
				}
			}
		}
	}
}

/**
	Minimization function for PGE
*/
void ss_min_LP(			global_variable 	 gv,
						PC_type				*PC_read,

						obj_type 			*SS_objective,
						NLopt_type 			*NLopt_opt,
						bulk_info 	 		 z_b,
						SS_ref 			    *SS_ref_db,
						csd_phase_set  		*cp
){

	double r;
	int 	ph_id;
	int     pc_check;
	int 	act;
	clock_t u;
	for (int i = 0; i < gv.len_ss; i++){
		gv.n_min[i] = 0;
	}

	/* "liq" redundant-occurrence pseudocompound synthesis, gh and tc (see
	   docs/liq_pseudocompound_shortcut.md): locate liq's phase index and
	   decide whether it fires this cycle, and zero the allocation-free
	   accumulators for the weighted-least-squares Gamma refinement. Gated
	   on gv.liq_pc_synth_active so it can be fully disabled (legacy
	   per-occurrence NLopt path) without touching any other logic. */
	int    is_gh            = (strcmp(gv.research_group, "gh") == 0);
	int    is_tc             = (strcmp(gv.research_group, "tc") == 0);
	int    ph_id_liq        = -1;
	if (gv.liq_pc_synth_active && (is_gh || is_tc)){
		for (int iss = 0; iss < gv.len_ss; iss++){
			if (strcmp(gv.SS_list[iss], "liq") == 0){ ph_id_liq = iss; break; }
		}
	}
	/* gv.n_ss_ph[] is only ever written by LP_pc_composite (PGE_function.c),
	   which is not called anywhere in the ss_min_LP/run_LP path this
	   database uses - it is always 0 here (confirmed empirically), so N is
	   counted directly from cp[] instead, which is what's actually live. */
	int    N_liq = 0;
	if (ph_id_liq >= 0){
		for (int i = 0; i < gv.len_cp; i++){
			if (cp[i].ss_flags[0] == 1 && cp[i].id == ph_id_liq){ N_liq += 1; }
		}
	}
	int    liq_synth_active = (ph_id_liq >= 0) && (N_liq >= gv.gh_liq_pc_synth_threshold);
	int    liq_real_min_found = 0;
	double MtM[n_ox_all][n_ox_all];
	double Mtmu[n_ox_all];
	if (liq_synth_active){
		for (int a = 0; a < n_ox_all; a++){
			Mtmu[a] = 0.0;
			for (int b = 0; b < n_ox_all; b++){ MtM[a][b] = 0.0; }
		}
	}

	pc_check = gv.PC_checked;
	for (int i = 0; i < gv.len_cp; i++){

		if (cp[i].ss_flags[0] == 1){
			ph_id = cp[i].id;

			// deactivating the next part helps for IGAD database at VHT
			int is_liq_synth_candidate = (liq_synth_active && ph_id == ph_id_liq && !liq_real_min_found);
			if (liq_synth_active && ph_id == ph_id_liq){
				act = is_liq_synth_candidate ? 1 : 0;
			}
			else if ( strcmp( gv.SS_list[ph_id], "liq") == 0 && gv.n_min[ph_id] > gv.n_max_val){
				act = 0;
			}
			else{
				act = 1;
			}
			gv.n_min[ph_id] += 1;

			if (act == 1){
				cp[i].min_time		  		= 0.0;								/** reset local minimization time to 0.0 */
				u = clock(); 
				/**
					set the iguess of the solution phase to the one of the considered phase 
				*/
				for (int k = 0; k < cp[i].n_xeos; k++) {
					SS_ref_db[ph_id].iguess[k] 	= cp[i].xeos[k];
					cp[i].xeos_0[k] 			= cp[i].xeos[k];
					// SS_ref_db[ph_id].dguess[k] = cp[i].xeos[k];			//dguess can be used of LP, it is used for PGE to check for drifting
				}

				/**
					Rotate G-base hyperplane
				*/
				SS_ref_db[ph_id] = rotate_hyperplane(		gv, 
															SS_ref_db[ph_id]		);

				/**
					Define a sub-hypervolume for the solution phases bounds
				*/
				SS_ref_db[ph_id] = restrict_SS_HyperVolume(	gv, 
															SS_ref_db[ph_id],
															gv.box_size_mode_LP		);

				/**
					call to NLopt for non-linear + inequality constraints optimization
				*/
				SS_ref_db[ph_id] = (*NLopt_opt[ph_id])(		gv,
															SS_ref_db[ph_id]		);

				/**
					gh only: gh phases use n_xeos==n_em direct p[]=x[] with
					Sigma(x)=1 enforced only as a soft NLopt equality
					constraint (no dependent-variable basis reduction like
					tc/sb use, so they can't need this). If NLopt returns an
					iterate that doesn't satisfy the constraint well (can
					happen on abnormal termination, e.g. NLOPT_ROUNDOFF_LIMITED
					- see [[gh-pseudocompound-validity-bugs]]), project the
					point back onto the simplex (renormalize) and retry once
					from there, rather than accepting/discarding a poorly-
					constrained result outright. Added 2026-07-14 per user
					request ("force a better satisfaction of the equality
					constraint") as a proactive complement to the reactive
					Sigma(p)=1 sf_ok validity check.
				*/
				if (strcmp(gv.research_group, "gh") == 0){
					double sum_x = 0.0;
					for (int k = 0; k < SS_ref_db[ph_id].n_xeos; k++){
						sum_x += SS_ref_db[ph_id].xeos[k];
					}
					if (fabs(sum_x - 1.0) > 1.0e-6 && sum_x > 0.0){
						for (int k = 0; k < SS_ref_db[ph_id].n_xeos; k++){
							SS_ref_db[ph_id].iguess[k] = SS_ref_db[ph_id].xeos[k] / sum_x;
						}
						SS_ref_db[ph_id] = (*NLopt_opt[ph_id])(		gv,
																	SS_ref_db[ph_id]		);
					}
				}
				/**
					print solution phase informations (print has to occur before saving PC)
				*/
			
				u = clock() - u;
				SS_ref_db[ph_id].LM_time = ((double)u)/CLOCKS_PER_SEC*1000.0; 

				if (gv.verbose == 1){
					SS_ref_db[ph_id] = SS_UPDATE_function(		gv, 
																SS_ref_db[ph_id], 
																z_b, 
																gv.SS_list[ph_id]		);

					print_SS_informations(  				gv,
															SS_ref_db[ph_id],
															ph_id					);
				}
				if (SS_ref_db[ph_id].status == -4){
					/* NLOPT_ROUNDOFF_LIMITED: NLopt still writes its best feasible
					   iterate back into .xeos[] on return regardless of status, so
					   that remains the right array to read here. Fixed 2026-07-14:
					   this branch used to substitute SS_ref_db[ph_id].p[] instead,
					   which is NOT the final optimum - each phase's objective
					   function overwrites .p[] on every internal call (rejected
					   line-search steps, gradient finite-difference probes, etc.),
					   so by the time NLopt returns it can hold an arbitrary, often
					   physically invalid (sum(p) not 1, pinned at boiled-out
					   epsilon bounds) intermediate trial point instead of the
					   accepted solution. Confirmed via a real reproduction: rhm
					   (whose Gex evaluation is numerically stiffer due to its own
					   internal 3-parameter order Newton solve, so it hits
					   ROUNDOFF_LIMITED often - 11 times in one single-point solve,
					   liq 24 times) was retaining these corrupted compositions as
					   legitimate pseudocompounds, producing spurious near-identical
					   "solvus" instances with endmember fractions that didn't even
					   sum to 1. */
					if (gv.verbose == 1){
						printf(" Round-off error in the minimization of %4s\n",gv.SS_list[ph_id]);
					}
				}
				for (int k = 0; k < cp[i].n_xeos; k++) {
					cp[i].xeos_1[k] 			 =  SS_ref_db[ph_id].xeos[k];
				}
				/** 
					Here if the number of phase occurence in the LP matrix is equal to we add 2 pseudocompounds
				*/
				// /*
				// if (gv.n_ss_ph[ph_id] == 1){

				// 	double shift = 0.0;
				// 	double sh_array[] = {-0.01,0.001,0.001,0.01};
				// 	for (int add = 0; add < 2; add++){
				// 		shift = sh_array[add];
				// 		for (int k = 0; k < cp[i].n_xeos; k++) {
				// 			SS_ref_db[ph_id].iguess[k]   =  cp[i].xeos_1[k] * (1.0-shift) + cp[i].xeos_0[k] * (shift);
				// 		}

				// 		SS_ref_db[ph_id] = PC_function(				gv,
				// 													PC_read,
				// 													SS_ref_db[ph_id], 
				// 													z_b,
				// 													ph_id 		);
																
				// 		SS_ref_db[ph_id] = SS_UPDATE_function(		gv, 
				// 													SS_ref_db[ph_id], 
				// 													z_b, 
				// 													gv.SS_list[ph_id]		);

				// 		if (SS_ref_db[ph_id].sf_ok == 1){
				// 			copy_to_Ppc(							0,
				// 													1,
				// 													ph_id,
				// 													gv,

				// 													SS_objective,
				// 													SS_ref_db					);	
				// 		}
				// 	}
				// }
				// */

				for (int k = 0; k < cp[i].n_xeos; k++) {
					SS_ref_db[ph_id].iguess[k]   =  cp[i].xeos_1[k];
				}
				SS_ref_db[ph_id] = PC_function(				gv,
															PC_read,
															SS_ref_db[ph_id], 
															z_b,
															ph_id 		);
														
				SS_ref_db[ph_id] = SS_UPDATE_function(		gv,
															SS_ref_db[ph_id],
															z_b,
															gv.SS_list[ph_id]		);

				int liq_candidate_sf_ok = SS_ref_db[ph_id].sf_ok;

				/**
					add minimized phase to LP PGE pseudocompound list
				*/
				if (SS_ref_db[ph_id].sf_ok == 1){
					copy_to_Ppc(							pc_check,
															1,
															ph_id,
															gv,

															SS_objective,
															SS_ref_db					);	
				}
				else{
					for (int k = 0; k < cp[i].n_xeos; k++) {
						SS_ref_db[ph_id].iguess[k]   =  cp[i].xeos_0[k];
					}
					
					SS_ref_db[ph_id] = PC_function(				gv,
																PC_read,
																SS_ref_db[ph_id], 
																z_b,
																ph_id 		);
															
					SS_ref_db[ph_id] = SS_UPDATE_function(		gv, 
																SS_ref_db[ph_id], 
																z_b, 
																gv.SS_list[ph_id]		);

					copy_to_Ppc(								0,
																1,
																ph_id,
																gv,

																SS_objective,
																SS_ref_db					);	
			
				}


				if (is_liq_synth_candidate && liq_candidate_sf_ok == 1){
					liq_real_min_found = 1;
				}

				if (is_liq_synth_candidate && liq_real_min_found && ph_id == ph_id_liq){
					double g_row[LIQ_PC_SYNTH_MAX_DIM];
					if (is_gh && SS_ref_db[ph_id].n_xeos <= LIQ_PC_SYNTH_MAX_DIM){
						GH_liq_eval_raw(			gv, SS_objective, SS_ref_db, ph_id,
													SS_ref_db[ph_id].iguess, g_row		);
						GH_liq_gamma_ls_accumulate(z_b, SS_ref_db[ph_id], g_row, MtM, Mtmu);
					}
					else if (is_tc && SS_ref_db[ph_id].n_em <= LIQ_PC_SYNTH_MAX_DIM){
						TC_eval_raw(			gv, SS_objective, SS_ref_db, ph_id,
													SS_ref_db[ph_id].iguess, g_row		);
						GH_liq_gamma_ls_accumulate(z_b, SS_ref_db[ph_id], g_row, MtM, Mtmu);
					}
				}


			}

		}	
	}

	/* stage C: refine Gamma and, for liq only, replace the redundant
	   per-occurrence NLopt solves with 2*N analytic pseudocompounds on
	   the refined hyperplane (see docs/liq_pseudocompound_shortcut.md).
	   Requires liq_real_min_found too: if every occurrence's real
	   minimization turned out invalid (see the per-candidate validity
	   check above), there is no valid x* to build any of this from, and
	   every occurrence already got a full legacy minimization instead. */
	if (liq_synth_active && liq_real_min_found){
		double gamma_new[n_ox_all];
		for (int a = 0; a < n_ox_all; a++){ gamma_new[a] = gv.gam_tot[a]; }	/* oxides absent from the bulk keep this prior untouched */
		if (GH_liq_gamma_solve(z_b, MtM, Mtmu, gv.gam_tot, gamma_new) == 1){
			if (gv.verbose == 1){
				printf(" [liq pc synth] Gamma_new = [");
				for (int a = 0; a < gv.len_ox; a++){ printf(" %+10f", gamma_new[a]); }
				printf(" ]\n");
			}
			if 		(is_gh){ GH_liq_pc_synth(gv, PC_read, SS_objective, z_b, SS_ref_db, cp, ph_id_liq, gamma_new); }
			else if (is_tc){ TC_liq_pc_synth(gv, PC_read, SS_objective, z_b, SS_ref_db, cp, ph_id_liq, gamma_new); }
		}
		else{
			if (gv.verbose == 1){
				printf(" [liq pc synth] Gamma_new solve failed (singular), falling back to gv.gam_tot = [");
				for (int a = 0; a < gv.len_ox; a++){ printf(" %+10f", gv.gam_tot[a]); }
				printf(" ]\n");
			}
			if 		(is_gh){ GH_liq_pc_synth(gv, PC_read, SS_objective, z_b, SS_ref_db, cp, ph_id_liq, gv.gam_tot); }
			else if (is_tc){ TC_liq_pc_synth(gv, PC_read, SS_objective, z_b, SS_ref_db, cp, ph_id_liq, gv.gam_tot); }
		}
	}

};

/**
  initialize solution phase database
**/
global_variable init_ss_db(		int 				 EM_database,
								bulk_info 	 		 z_b,
								global_variable 	 gv,
								SS_ref 				*SS_ref_db
){


	if (EM_database == 0){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_mp_EM_function(	gv, 
													SS_ref_db[i], 
													gv.EM_dataset, 
													z_b, 
													gv.SS_list[i]		);
		}
	}
	else if (EM_database == 1){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_mb_EM_function(	gv, 
													SS_ref_db[i], 
													gv.EM_dataset, 
													z_b, 
													gv.SS_list[i]		);
		}
	}
	else if (EM_database == 11){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_mb_ext_EM_function(	gv, 
														SS_ref_db[i], 
														gv.EM_dataset, 
														z_b, 
														gv.SS_list[i]		);
		}
	}
	else if (EM_database == 2){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_ig_EM_function(	gv, 
													SS_ref_db[i], 
													gv.EM_dataset, 
													z_b, 
													gv.SS_list[i]		);
		}
	}
	else if (EM_database == 22){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_igd_EM_function(	gv, 
													SS_ref_db[i], 
													gv.EM_dataset, 
													z_b, 
													gv.SS_list[i]		);
		}
	}
	else if (EM_database == 3){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_igad_EM_function(	gv, 
														SS_ref_db[i], 
														gv.EM_dataset, 
														z_b, 
														gv.SS_list[i]		);
		}
	}
	else if (EM_database == 4 ){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_um_EM_function(	gv, 
													SS_ref_db[i], 
													gv.EM_dataset, 
													z_b, 
													gv.SS_list[i]		);
		}
	}
	else if (EM_database == 5 ){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_um_ext_EM_function(	gv, 
														SS_ref_db[i], 
														gv.EM_dataset, 
														z_b, 
														gv.SS_list[i]		);
		}
	}
	else if (EM_database == 6 ){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_mtl_EM_function(		gv, 
														SS_ref_db[i], 
														gv.EM_dataset, 
														z_b, 
														gv.SS_list[i]		);
		}
	}
	else if (EM_database == 7 ){
		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_mpe_EM_function(		gv, 
														SS_ref_db[i], 
														gv.EM_dataset, 
														z_b, 
														gv.SS_list[i]		);
		}
	}

	return gv;
};



/**
  initialize solution phase database
**/
global_variable init_ss_db_sb(	int 				 EM_database,
								bulk_info 	 		 z_b,
								global_variable 	 gv,
								SS_ref 				*SS_ref_db
){

	if (EM_database == 0){

		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_sb11_EM_function(	gv, 
														SS_ref_db[i], 
														gv.EM_dataset, 
														z_b, 
														gv.SS_list[i]		);
											
										/** can become a global variable instead */
		}
	}
	else if (EM_database == 1){

		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_sb21_EM_function(	gv, 
														SS_ref_db[i], 
														gv.EM_dataset, 
														z_b, 
														gv.SS_list[i]		);
											
										/** can become a global variable instead */
		}
	}
	else if (EM_database == 2){

		for (int i = 0; i < gv.len_ss; i++){
			SS_ref_db[i].P  = z_b.P;									/** needed to pass to local minimizer, allows for P variation for liq/sol */
			SS_ref_db[i].T  = z_b.T;		
			SS_ref_db[i].R  = 0.0083144;

			// if (SS_ref_db[i].is_liq == 1){
			// 	SS_ref_db[i].P  = z_b.P + gv.melt_pressure;
			// }

			SS_ref_db[i]    = G_SS_sb24_EM_function(	gv, 
														SS_ref_db[i], 
														gv.EM_dataset,
														z_b,
														gv.SS_list[i]		);

										/** can become a global variable instead */
		}
	}
	return gv;
};

/**
  initialize solution phase database for the "gh" (Ghiorso/MELTS) research group
**/
global_variable init_ss_db_gh(	int 				 EM_database,
								bulk_info 	 		 z_b,
								global_variable 	 gv,
								SS_ref 				*SS_ref_db
){
	/* see init_em_db_gh's own comment on why this is also set here, not
	   just in GH_SS_objective_init_function. */
	GH_actual_EM_database = gv.EM_database;
	for (int i = 0; i < gv.len_ss; i++){
		SS_ref_db[i].P  = z_b.P;
		SS_ref_db[i].T  = z_b.T;
		SS_ref_db[i].R  = 0.0083144;

		SS_ref_db[i]    = G_SS_gh_EM_function(	gv,
												SS_ref_db[i],
												gv.EM_dataset,
												z_b,
												gv.SS_list[i]		);
	}
	return gv;
};
