********************************************************************************
********                	     STRATIFY - CVD risk            	  	********
********************************************************************************
** Software: STATA/MP 16.1									   	   			  **
** This programme is used to conduct calibration & validation of the CVD model**

* calibration of CVD model *
stset TimetoEvent, failure(CVD_comp= 1) scale(365) 

********** calibration at 10 years ********************
*******************************************************

* create evenly sized groups
egen groups = cut(combPI), group(50)
replace groups = groups + 1

local k oxford
global grps 50
forvalues i = 1/$grps {
stpci if groups == `i', at(10) gen(pajcif`i'_`k'_d)  competingvalues(2)                
}

* create the pseudo observation into one column

gen pseudo = .
global grps 50 
local k oxford
forvalues i = 1/$grps {
    replace pseudo = pajcif`i'_`k'_d if groups == `i'
}

* calibration curve

running pseudo risk10y , ci leg(off) xtitle(Predicted probabilities) ytitle(Observed probabilities)  xlabel(0(.1)1) nopts addplot(function y = x, ra(risk10y) clpat(solid) clcol(green)) ylabel(0(0.2)1) xlabel(0(0.2)1) aspect(1) yscale(range(0 1)) xscale(range(0 1)) title("") graphregion(color(white)) lineopts(lcolor(maroon))

snapshot save

drop groups
egen groups = cut(combPI), group(10)
replace groups = groups + 1
gen cal_obs_1y=.
forvalues i=1/10 {
    stcompet np_cif_1y`i'=ci if groups==`i', compet1(2) 
       sum np_cif_1y`i' if _d==1 & groups==`i' & _t<=10
    replace cal_obs_1y=r(max) if _d==1 & groups==`i'
}

* calculate average predicted in each group
gen cal_pred_1y = .
forvalues i=1/10 {
    sum risk10y if groups==`i' 
    replace cal_pred_1y=r(mean) if groups==`i'
}
				   
running pseudo risk10y , ci leg(off) xtitle("Predicted probabilities") ytitle("Observed probabilities") nopts lcolor(maroon)  addplot(scatter cal_obs_1y cal_pred_1y, msize(small) msymbol(circle_hollow) color(orange) || function y = x, clpat(solid) lcolor(green)) title("")  graphregion(color(white)) aspect(1) saving(cal10, replace)
						
hist risk10y, xtitle(Distribution of predictions) ytitle("") ylabel(0 (40) 40) col(maroon) xlabel(0 (0.1)0.4) yscale(range(-1 35)) xscale(range(0 0.4)) graphregion(color(white)) fysize(18) fxsize(62.2) saving(hist10ys, replace) 

graph combine cal10.gph hist10ys.gph , col(1) graphregion(color(white)) plotregion(margin(zero)) 
							  
							
* Recalibration of 10 year model
						
* linear recalibrate 
snapshot restore 2

glm pseudo combPI , link(logit)  vce(robust) 
predict logit_P_new , xb
gen P_modelnew = exp(logit_P_new)/(1 + exp(logit_P_new))
gen check = invlogit(logit_P_new)
 
running pseudo P_modelnew , ci leg(off) xtitle(Predicted probabilities) ytitle(Observed probabilities)  xlabel(0(.1)1) nopts addplot(function y = x, ra(P_modelnew) clpat(solid) clcol(black)) ylabel(0(0.1)1) xlabel(0(0.1)1) aspect(1) yscale(range(-0.05 1)) xscale(range(-0.05 1)) title("") graphregion(color(white)) lineopts(lcolor(gs8))

* nonlinear recalibrate 
fracpoly, degree(2) center(no) noscaling powers(1 1) :   glm pseudo combPI , link(logit)  vce(robust) 
predict logit_P_new2, xb
gen P_modelnew2 = exp(logit_P_new2)/(1 + exp(logit_P_new2))
running pseudo P_modelnew2 , ci leg(off) xtitle(Predicted probabilities) ytitle(Observed probabilities)  xlabel(0(.1)1) nopts addplot(function y = x, ra(P_modelnew2) clpat(solid) clcol(green)) ylabel(0(0.1)0.6) xlabel(0(0.1)0.6) aspect(1) yscale(range(0 0.6)) xscale(range(0 0.6)) title("") graphregion(color(white))  lineopts(lcolor(maroon))

gen logitp_check =  -7.228996 + (3.281498 * (combPI+.1708773532176397)) + ( -1.094527  *(combPI+.1708773532176397)*ln(combPI+.1708773532176397))
gen p_check = exp(logitp_check)/(1 + exp(logitp_check))	

gen P_modelnew2_v2 = exp(logit_P_new2)/(1 + exp(logit_P_new2))
running pseudo p_check , ci leg(off) xtitle(Predicted probabilities) ytitle(Observed probabilities)  xlabel(0(.1)1) nopts addplot(function y = x, ra(p_check) clpat(solid) clcol(green)) ylabel(0(0.1)0.6) xlabel(0(0.1)0.6) aspect(1) yscale(range(0 0.4)) xscale(range(0 0.6)) title("") graphregion(color(white))  lineopts(lcolor(maroon))					
						
drop groups 
egen groups = cut(logitp_check), group(10)
replace groups = groups + 1
gen cal_obs_1y=.
forvalues i=1/10 {
    stcompet np_cif_1y`i'=ci if groups==`i', compet1(2) 
       sum np_cif_1y`i' if _d==1 & groups==`i' & _t<=10
    replace cal_obs_1y=r(max) if _d==1 & groups==`i'
}

* calculate average predicted in each group
gen cal_pred_1y = .
forvalues i=1/10 {
    sum P_modelnew2 if groups==`i' 
    replace cal_pred_1y=r(mean) if groups==`i'
}

running  pseudo p_check, ci leg(off) xtitle("Predicted probabilities") ytitle("Observed probabilities") nopts lcolor(maroon)  xlabel(0(.1)0.6) ylabel(0(0.1)0.6) addplot(scatter cal_obs_1y cal_pred_1y, msize(small) msymbol(circle_hollow) color(orange)  xlabel(0(.1)0.6) ylabel(0(0.1)0.6)|| function y = x, clpat(solid) ra(p_check) lcolor(green)) title("")  graphregion(color(white)) aspect(1)  xlabel(0(.1)0.6) ylabel(0(0.1)0.6)  saving(recal10, replace)
						
hist P_modelnew2, xtitle(Distribution of predictions) ytitle("") ylabel(0 (40) 40) col(maroon) xlabel(0 (0.1)0.4) yscale(range(-1 35)) xscale(range(0 0.4)) graphregion(color(white)) fysize(18) fxsize(62.2) saving(hist10ys, replace) 

graph combine recal10.gph hist10ys.gph , col(1) graphregion(color(white)) plotregion(margin(zero)) 
							  
	
*External validation	
*generate STRATIFY-CVD score
**Standardise varibles across datasets***

*** Qrisk2 variables ****

generate Qrisk2age=age
replace Qrisk2age=84 if Qrisk2age>84  /* maximum value for age in Qrisk2 is 84. We assume maximal risk for people over the age of 84 */

*generate STRATIFY-CVD score
**Standardise varibles across datasets, use imputed dataset***

generate Black_ethnicity=1 if IMP_Ethnicity==1
replace Black_ethnicity=0 if Black_ethnicity==.

generate SouthAsian_ethnicity=1 if IMP_Ethnicity==2
replace SouthAsian_ethnicity=0 if SouthAsian_ethnicity==.

generate Other_ethnicity=1 if IMP_Ethnicity==3 
replace Other_ethnicity=0 if Other_ethnicity==.

generate IMD2=1 if IMP_IMD==2
replace IMD2=0 if IMD2==.

generate IMD3=1 if IMP_IMD==3
replace IMD3=0 if IMD3==.

generate IMD4=1 if IMP_IMD==4
replace IMD4=0 if IMD4==.

generate IMD5=1 if IMP_IMD==5
replace IMD5=0 if IMD5==.

generate Ex_smoker=1 if IMP_Smoking==1
replace Ex_smoker=0 if Ex_smoker==.

generate Current_smoker=1 if IMP_Smoking==2
replace Current_smoker=0 if Current_smoker==.

generate underweight=1 if IMP_BMI==0
replace underweight=0 if underweight==.

generate overweight=1 if IMP_BMI==2  
replace overweight=0 if overweight==.

generate obese=1 if IMP_BMI==3
replace obese=0 if obese==.

generate Mobobese=1 if IMP_BMI==4
replace Mobobese=0 if Mobobese==.

gen TC_HDL= IMP_TC_HDL 

****Estimate CVD Risk*****
gen Age=((age/100)^3)-0.2419

gen combPI=((Age * 4.1188) + (Gender*(-0.3629)) + (TC_HDL * 0.0709) + (sysbp_inclusion * 0.0062) + (IMD2 * 0.0824) + (IMD3 * 0.1367) + (IMD4 * 0.2218) + (IMD5 * 0.3609) + (underweight * 0.1129) + (overweight * (-0.0142)) + (obese * (-0.0042)) + (Mobobese * 0.0878) + (Black_ethnicity * (-0.0392)) + (SouthAsian_ethnicity * 0.1656) +  (Other_ethnicity * (-0.0159)) + (Ex_smoker * 0.1073) + (Current_smoker * 0.4761)  + (TreatedHypertension * 0.4310) + (RheumatoidArthritis * 0.3795) + (Diabetes * 0.4312) + (Qrisk2FH * 0.0173) + (CKD * 0.0089) + (Atrialfib * 0.4017)  + (prevStroke * 0.9702) + (prevMI * 0.6906) + (HF * 0.5893) + (statins * 0.1730))

gen risk10y= 1-(((1- 0.0146165))^exp(combPI))
gen risk5y= 1-(((1- 0.00699506))^exp(combPI))
gen risk1y= 1-(((1- 0.001500368))^exp(combPI))



gen recal_10yr =  -7.228996 + (3.281498 * (combPI+.1708773532176397)) + ( -1.094527  *(combPI+.1708773532176397)*ln(combPI+.1708773532176397))
gen recal_risk10y = exp(recal_10yr)/(1 + exp(recal_10yr))	


gen recal_5yr =  -8.06039  + (3.165212 * (combPI+.1708773532176397)) + (-.9804859 *(combPI+.1708773532176397)*ln(combPI+.1708773532176397))
gen recal_risk5y = exp(recal_5yr)/(1 + exp(recal_5yr))




********** calibration at 10 years ********************
*******************************************************
stset TimetoEvent, failure(CVD_comp= 1) scale(365) id(patid)

* create evenly sized groups
egen groups = cut(combPI), group(50)
replace groups = groups + 1

local k oxford
global grps 50
forvalues i = 1/$grps {
stpci if groups == `i', at(10) gen(pajcif`i'_`k'_d)  competingvalues(2)                
}

* create the pseudo observation into one column

gen pseudo = .
global grps 50 
local k oxford
forvalues i = 1/$grps {
    replace pseudo = pajcif`i'_`k'_d if groups == `i'
}

* calibration curve

running pseudo risk10y , ci leg(off) xtitle(Predicted probabilities) ytitle(Observed probabilities)  xlabel(0(.1)1) nopts addplot(function y = x, ra(risk10y) clpat(solid) clcol(green)) ylabel(0(0.2)1) xlabel(0(0.2)1) aspect(1) yscale(range(0 1)) xscale(range(0 1)) title("") graphregion(color(white)) lineopts(lcolor(maroon))

snapshot save

drop groups
egen groups = cut(combPI), group(10)
replace groups = groups + 1
gen cal_obs_1y=.
forvalues i=1/10 {
    stcompet np_cif_1y`i'=ci if groups==`i', compet1(2) 
       sum np_cif_1y`i' if _d==1 & groups==`i' & _t<=10
    replace cal_obs_1y=r(max) if _d==1 & groups==`i'
}

* calculate average predicted in each group
gen cal_pred_1y = .
forvalues i=1/10 {
    sum risk10y if groups==`i' 
    replace cal_pred_1y=r(mean) if groups==`i'
}
				   
running pseudo risk10y , ci leg(off) xtitle("Predicted probabilities") ytitle("Observed probabilities") nopts lcolor(maroon)  addplot(scatter cal_obs_1y cal_pred_1y, msize(small) msymbol(circle_hollow) color(orange) || function y = x, clpat(solid) lcolor(green)) title("") graphregion(color(white)) aspect(1) saving(cal10, replace)
						
hist risk10y, xtitle(Distribution of predictions) ytitle("") ylabel(0 (40) 40) col(maroon) xlabel(0 (0.1)0.4) yscale(range(-1 35)) xscale(range(0 0.4)) graphregion(color(white)) fysize(18) fxsize(62.2) saving(hist10ys, replace) 

graph combine cal10.gph hist10ys.gph , col(1) graphregion(color(white)) plotregion(margin(zero)) 


********** Recalibration at 10 years ******************
*******************************************************

* calculate average predicted in each group - recalibrate
gen recal_pred_1y = .
forvalues i=1/10 {
    sum recal_risk10y if groups==`i' 
    replace recal_pred_1y=r(mean) if groups==`i'
}
				   
running pseudo recal_risk10y , ci leg(off) xtitle("Predicted probabilities") ytitle("Observed probabilities") nopts lcolor(maroon) addplot(scatter cal_obs_1y recal_pred_1y, msize(small) msymbol(circle_hollow) color(orange) || function y = x, clpat(solid) lcolor(green)) title("")  graphregion(color(white)) aspect(1) saving(recal10, replace)
						
hist recal_risk10y, xtitle(Distribution of predictions) ytitle("") ylabel(0 (40) 40) col(maroon) xlabel(0 (0.1)0.4) yscale(range(-1 35)) xscale(range(0 0.4)) graphregion(color(white)) fysize(18) fxsize(62.2) saving(recalhist10ys, replace) 

graph combine recal10.gph recalhist10ys.gph , col(1) graphregion(color(white)) plotregion(margin(zero)) 


	