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
#ifndef __GH_FLUID_EOS_H_
#define __GH_FLUID_EOS_H_

double GH_pitzer_sterner_G(int is_H2O, double T, double Pbar);
double GH_pitzer_sterner_H2O_hTs_G(double T, double Pbar);
double GH_pitzer_sterner_mix_G(double x_h2o, double T, double Pbar, double *dGdx_h2o);
double GH_haar_H2O_G(double T, double P);
double GH_wdh78_G(double T, double P);
double GH_duan_CO2_G(double T);

#endif
