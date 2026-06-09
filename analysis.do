/*===============================================
  Replication File
  Empowering Women Through Radio: 
  Evidence from Occupied Japan
  
  Author: Yoko Okuyama (Uppsala University)
  Journal of Development Economics, 2026
  
  This file replicates all main tables 
  Table 1, 2, 3, and 4
  ===============================================*/

clear all

*==============================================================================*
* USER: CHANGE THIS PATH TO YOUR REPLICATION PACKAGE LOCATION
*==============================================================================*
global replication_path "C:/Users/yokok634/Dropbox/00_Japan_population/replication_package/"
* Example: "C:/Users/YourName/Desktop/okuyama_replication"
* Use forward slashes (/) not backslashes (\)

*==============================================================================*
* The following paths are set automatically - DO NOT CHANGE
*==============================================================================*
global data    "${replication_path}/data"
global output  "${replication_path}/output"

* Create output subdirectories if they don't exist
capture mkdir "${output}"
capture mkdir "${output}/tables"

  
set more off, permanently
set scheme s1mono, permanently
graph set window fontface "Times New Roman"
set linesize 150

*
* Set the basic set of covariates
*
global var_base "i.NEAR_DIST_DECILEBIN i.radio_station i.prefecture LOG_HH_10000 HEAVY_DAMAGE"
global setstarlevels "* 0.10 ** 0.05 *** 0.01"

/*----- Table 2: Political effects  -----*/

clear all
use "${data}/analysis_data_clean.dta",clear
keep if sample_baseline == 1 
keep if WOMEN_TURNOUT != .

sum HH_SUBSCRBE_TOTAL
local total_hh_freq = r(sum)
gen popu_weight     = HH_SUBSCRBE_TOTAL / `total_hh_freq'


* First Stage
reg RADIO_SUBSCRIPTION1946_SDUNIT wavg_signal_sdunit  ${var_base}, cl(radio_station) 

est store fstage_base
	local r2c: display %4.2f `e(r2)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	local maincoef      =  _b[wavg_signal_sdunit]

***** WOMEN *****
* Reduced form
reg WOMEN_TURNOUT wavg_signal_sdunit ${var_base}, cl(radio_station) 

est store rf_women_base
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ WOMEN_TURNOUT
	local mu : di %4.2f r(mean)
	estadd local mu `mu'

* OLS
reg WOMEN_TURNOUT RADIO_SUBSCRIPTION1946_SDUNIT ${var_base}, cl(radio_station) 
est store ols_women_base
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ WOMEN_TURNOUT
	local mu : di %4.2f r(mean)
	estadd local mu `mu'

*tsls
ivreg2  WOMEN_TURNOUT ${var_base} (RADIO_SUBSCRIPTION1946_SDUNIT = wavg_signal_sdunit), cl(radio_station) 

est store tsls_women_base
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"	
	estadd local CL "Yes"
	
	qui summ WOMEN_TURNOUT
	local mu_turnout = r(mean)
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
    local sddev =  r(sd)
	
***** MEN *****
* Reduced form
drop if missing(MEN_TURNOUT)
reg MEN_TURNOUT wavg_signal_sdunit ${var_base}, cl(radio_station)
est store rf_men_base
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ MEN_TURNOUT
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
	
* OLS
reg MEN_TURNOUT RADIO_SUBSCRIPTION1946_SDUNIT ${var_base}, cl(radio_station) 
est store ols_men_base
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ MEN_TURNOUT
	local mu : di %4.2f r(mean)
	estadd local mu `mu'

ivreg2  MEN_TURNOUT ${var_base} (RADIO_SUBSCRIPTION1946_SDUNIT = wavg_signal_sdunit), cl(radio_station) 
est store tsls_men_base
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ MEN_TURNOUT
	local mu : di %4.2f r(mean)
	estadd local mu `mu'

*
* Export
*
label var RADIO_SUBSCRIPTION1946_SDUNIT "\shortstack[l]{Radio subscription \\ in std.dev. unit}" 
label var wavg_signal_sdunit "\shortstack[l]{Field strength \\ in std.dev. unit}"  
label var LOG_HH_10000 "log N.of HH"
label var WOMEN_TURNOUT "Women's turnout"
label var MEN_TURNOUT "Men's turnout"
local labwomen: variable label `var' WOMEN_TURNOUT
local labmen: variable label `var' MEN_TURNOUT
local labradiosub: variable label `var' RADIO_SUBSCRIPTION1946_SDUNIT

#delimit;
esttab  fstage_base rf_women_base ols_women_base tsls_women_base rf_men_base ols_men_base tsls_men_base 
        using "${output}/tables/TURNOUT-OLS-IV-result-for-replication.tex",
        replace
		compress 
		b(a3) 
		se(a3) 
		star(${setstarlevels})
		noconstant 
		sfmt(%4.2f)
		mtitle("FS" "RF" "OLS" "TSLS" "RF" "OLS" "TSLS" )
		mgroups("`labradiosub'" "`labwomen'" "`labmen'", pattern(1 1 0 0  1 0 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}))
		keep(wavg_signal_sdunit RADIO_SUBSCRIPTION1946_SDUNIT)
		order(wavg_signal_sdunit RADIO_SUBSCRIPTION1946_SDUNIT LOG_HH_10000)
		label
		booktabs
		substitute(\_ _)
	    scalar("FSTAT First stage F-stat"   
			"PFE Prefecture FE" 
			"SFE Transmitter FE" 
			"DT Distance control"  
			"CL Log N of hh and WWII damage"
			"mu Mean outcome");
#delimit cr  


/*---------- Table 4: Birth, marriage, employment ---------*/
discard
clear all
use "${data}/analysis_data_clean.dta", clear
keep if sample_baseline == 1 

***** first stage ***** 
reg RADIO_SUBSCRIPTION1946_SDUNIT wavg_signal_sdunit  ${var_base}, cl(radio_station) 
est store fstage_base
	local r2c: display %4.2f `e(r2)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE  "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"

***** reduced form *****   
reg CBR1947 wavg_signal_sdunit ${var_base}, cl(radio_station) 
est store rf47_base
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CBR1947
	local mu : di %4.2f r(mean)
	estadd local mu `mu'

reg CBR1947 wavg_signal_sdunit ${var_base} MFRATIO1947, cl(radio_station) 
est store rf47_add
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CBR1947
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
	
reg CBR1950 wavg_signal_sdunit ${var_base}, cl(radio_station) 
est store rf50_base
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CBR1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
	
reg CBR1950 wavg_signal_sdunit ${var_base} MFRATIO1950, cl(radio_station) 
est store rf50_add
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CBR1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'		
	
** CSBR1950	
reg CSBR1950 wavg_signal_sdunit ${var_base}, cl(radio_station) 
est store rf_csbr50_base
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CSBR1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'	

** Employment Rate
reg LFP1950_NOFAMILYEMP_SEX2 wavg_signal_sdunit ${var_base}, cl(radio_station) 
est store rf_lfp_base
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ LFP1950_NOFAMILYEMP_SEX2
	local mu : di %4.2f r(mean)
	estadd local mu `mu' 	
	
** CRMARRIAGE1950
reg CRMARRIAGE1950 wavg_signal_sdunit ${var_base}, cl(radio_station) 
est store rf_marriage50_base
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CRMARRIAGE1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'

** 	CRDIVORCE1950
reg CRDIVORCE1950 wavg_signal_sdunit ${var_base}, cl(radio_station) 
est store rf_divorce50_base
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CRDIVORCE1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
	 
	
***** OLS *****
reg CBR1947 RADIO_SUBSCRIPTION1946_SDUNIT ${var_base}, cl(radio_station) 
est store ols_cbr1947_base
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CBR1947
	local mu : di %4.2f r(mean)
	estadd local mu `mu'

reg CBR1947 RADIO_SUBSCRIPTION1946_SDUNIT ${var_base} MFRATIO1947, cl(radio_station) 
est store ols_cbr1947_add
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CBR1947
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
	
	
reg CBR1950 RADIO_SUBSCRIPTION1946_SDUNIT ${var_base}, cl(radio_station) 
est store ols_cbr1950_base
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CBR1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
	
reg CBR1950 RADIO_SUBSCRIPTION1946_SDUNIT ${var_base} MFRATIO1950, cl(radio_station) 
est store ols_cbr1950_add
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CBR1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'

* labor force participation
reg  LFP1950_NOFAMILYEMP_SEX2 RADIO_SUBSCRIPTION1946_SDUNIT ${var_base}, cl(radio_station)
est store ols_lfp_base
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ LFP1950_NOFAMILYEMP_SEX2
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
	
reg CRMARRIAGE1950 RADIO_SUBSCRIPTION1946_SDUNIT ${var_base}, cl(radio_station) 
est store ols_marriage1950_base
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CRMARRIAGE1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
	
reg CRMARRIAGE1950 RADIO_SUBSCRIPTION1946_SDUNIT ${var_base} MFRATIO1950, cl(radio_station) 
est store ols_marriage1950_add
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CRMARRIAGE1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'		
	
reg CRDIVORCE1950 RADIO_SUBSCRIPTION1946_SDUNIT ${var_base}, cl(radio_station) 
est store ols_divorce1950_base
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CRDIVORCE1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
	
reg CRDIVORCE1950 RADIO_SUBSCRIPTION1946_SDUNIT ${var_base} MFRATIO1950, cl(radio_station) 
est store ols_divorce1950_add
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CRDIVORCE1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
	
***** TSLS *****
ivreg2  CBR1947 ${var_base} (RADIO_SUBSCRIPTION1946_SDUNIT = wavg_signal_sdunit), cl(radio_station) 
est store tsls_cbr1947_base
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CBR1947
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
	
ivreg2  CBR1950 ${var_base} (RADIO_SUBSCRIPTION1946_SDUNIT = wavg_signal_sdunit), cl(radio_station) 
est store tsls_cbr1950_base
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CBR1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
	
ivreg2  CSBR1950 ${var_base} (RADIO_SUBSCRIPTION1946_SDUNIT = wavg_signal_sdunit), cl(radio_station)
est store tsls_csbr1950_base
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CSBR1950 
	local mu : di %4.2f r(mean)
	estadd local mu `mu'
	
ivreg2  CRMARRIAGE1950 ${var_base} (RADIO_SUBSCRIPTION1946_SDUNIT = wavg_signal_sdunit), cl(radio_station)
est store tsls_marriage1950_base
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CRMARRIAGE1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'	
	
ivreg2  CRDIVORCE1950 ${var_base} (RADIO_SUBSCRIPTION1946_SDUNIT = wavg_signal_sdunit), cl(radio_station)
est store tsls_divorce1950_base
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ CRDIVORCE1950
	local mu : di %4.2f r(mean)
	estadd local mu `mu'		

* Employment rate (excluding family employment)
ivreg2  LFP1950_NOFAMILYEMP_SEX2 ${var_base} (RADIO_SUBSCRIPTION1946_SDUNIT =wavg_signal_sdunit), cl(radio_station)
est store tsls_lfp_base
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CL "Yes"
	qui summ LFP1950_NOFAMILYEMP_SEX2
	local mu : di %4.2f r(mean)
	estadd local mu `mu'	

***** var label *****
label var LFP1950_NOFAMILYEMP_SEX2 "\shortstack[c]{Women's employment rate \\ excld. family emp}"
label var FEMALE_LABOR_RATE1940 "Women's LFP in 1940"
label var wavg_signal_sdunit "\shortstack[l]{Field strength \\ in std.dev. unit}" 
label var RADIO_SUBSCRIPTION1946_SDUNIT "\shortstack[l]{Radio subscription \\ in std.dev. unit}" 
label var MFRATIO1950 "\shortstack[l]{Male to Female Ratio \\ in 1950}" 

label var RADIO_SUBSCRIPTION1946_SDUNIT "\shortstack[l]{Radio subscription \\ in std.dev. unit}" 
label var wavg_signal_sdunit "\shortstack[l]{Field strength \\ in std.dev. unit}" 
label var MFRATIO1947 "\shortstack[l]{Male to Female Ratio \\ in 1947}" 
label var MFRATIO1950 "\shortstack[l]{Male to Female Ratio \\ in 1950}" 
label var CBR1935 "\shortstack[l]{No. live births \\ per 1000 population in 1935}"
label var CBR1947 "\shortstack[l]{No. live births \\ per 1000 population in 1947}"
label var CBR1950 "\shortstack[l]{No. live births \\ per 1000 population in 1950}"
label var CRMARRIAGE1935 "\shortstack[l]{No. marriages \\ per 1000 population in 1935}"
label var CRMARRIAGE1950 "\shortstack[l]{No. marriages \\ per 1000 population in 1950}"
label var CRDIVORCE1935 "\shortstack[l]{No. divorces \\ per 1000 population in 1935}"
label var CRDIVORCE1950 "\shortstack[l]{No. divorces \\ per 1000 population in 1950}"
label var CSBR1950 "\shortstack[l]{No. still births \\ per 1000 population in 1950}"


local labradiosub: variable label `var' RADIO_SUBSCRIPTION1946_SDUNIT
local labbirth35: variable label `var' CBR1935
local labbirth50: variable label `var' CBR1950
local labbirth47: variable label `var' CBR1947
local labbirth50: variable label `var' CBR1950
local labsb50: variable label `var' CSBR1950
local labmar35: variable label `var' CRMARRIAGE1935
local labmar50: variable label `var' CRMARRIAGE1950
local labdiv35: variable label `var' CRDIVORCE1935
local labdiv50: variable label `var' CRDIVORCE1950
local lablabor2: variable label `var' LFP1950_NOFAMILYEMP_SEX2


*
* Export
*	

* Employment Rate, Marriages and divorces			
#delimit;
		esttab  ols_lfp_base tsls_lfp_base  rf_lfp_base ols_marriage1950_base tsls_marriage1950_base rf_marriage50_base  ols_divorce1950_base tsls_divorce1950_base rf_divorce50_base
		using "${output}/tables/MARRIAGE-OLS-IV-result-for-revsion.tex",  
				replace
				compress 
				booktabs
				b(a3) 
				se(a3) 
				star(${setstarlevels})
				noconstant 
				sfmt(%4.2f)
				mtitle("OLS"  "TSLS" "RF" "OLS"  "TSLS" "RF" "OLS"  "TSLS"  "RF")
				mgroups("`lablabor2'" "`labmar50'" "`labdiv50'", pattern(1 0 0 1 0 0 1 0 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}))
				keep(wavg_signal_sdunit RADIO_SUBSCRIPTION1946_SDUNIT)
				order(RADIO_SUBSCRIPTION1946_SDUNIT  wavg_signal_sdunit)
				label
				substitute(\_ _)
	    scalar("FSTAT First stage F-stat"   
			"PFE Prefecture FE" 
			"SFE Transmitter FE" 
			"DT Distance control"  
			"CL Log N of hh and WWII damage"
			"mu Mean outcome");
#delimit cr	

* Fertility main table
#delimit;
		esttab  ols_cbr1947_base tsls_cbr1947_base rf47_base ols_cbr1950_base tsls_cbr1950_base rf50_base tsls_csbr1950_base
		using "${output}/tables/CBR-OLS-IV-result-for-revsion.tex",  
				replace
				compress 
				booktabs
				b(a3) 
				se(a3) 
				star(${setstarlevels})
				noconstant 
				sfmt(%4.2f)
				mtitle("OLS" "TSLS" "RF" "OLS" "TSLS" "RF" "TSLS")
				mgroups("`labbirth47'" "`labbirth50'" "`labsb50'", pattern(1 0 0 1 0 0 1) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}))
				keep(wavg_signal_sdunit RADIO_SUBSCRIPTION1946_SDUNIT)
				label
				substitute(\_ _)
	    scalar("FSTAT First stage F-stat"   
			"PFE Prefecture FE" 
			"SFE Transmitter FE" 
			"DT Distance control"  
			"CL Log N of hh and WWII damage"
			"mu Mean outcome");
#delimit cr


/*----- Table 3 election outcome -----*/
 
* A unit of obs = district * candidate pair

use "${data}/analysis_data2_clean.dta", clear
gen weight     = 1 / n_district_per_candidate

keep if wavg_signal_sdunit != .
keep if WOMEN_TURNOUT != .

* Reduced form
reg VOTE_SHARE wavg_signal_sdunit ${var_base} [aweight = weight], cl(radio_station)
est store rf_voteshare_base
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CFE "Yes"
	estadd local CL "Yes"
	qui summ VOTE_SHARE
	local mu : di %4.2f r(mean)
	estadd local mu `mu'


* OLS
reg VOTE_SHARE RADIO_SUBSCRIPTION1946_SDUNIT ${var_base} [aweight = weight],  cl(radio_station) 
est store ols_voteshare_base
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"
	estadd local CFE "Yes"
	estadd local CL "Yes"
	qui summ VOTE_SHARE
	local mu : di %4.2f r(mean)
	estadd local mu `mu'

	
*tsls
ivreg2 VOTE_SHARE ${var_base} (RADIO_SUBSCRIPTION1946_SDUNIT = wavg_signal_sdunit) [aweight = weight], cl(radio_station) 
est store tsls_voteshare_base
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"	
	estadd local CFE "Yes"
	estadd local CL "Yes"
	qui summ VOTE_SHARE
	local mu : di %4.2f r(mean)
	estadd local mu `mu'

local b_pct = 100*_b[RADIO_SUBSCRIPTION1946_SDUNIT]

*back of the envelop calculation
gen VOTE_SHARE_CF = 	VOTE_SHARE - _b[RADIO_SUBSCRIPTION1946_SDUNIT]*RADIO_SUBSCRIPTION1946_SDUNIT
summ VOTE_SHARE_CF
local mu_share_cf =  r(mean)

	
* Calculate Persuasion Rate
ivreg2 VOTE_SHARE ${var_base} (RADIO_SUBSCRIPTION1946_PCT = wavg_signal_sdunit) [aweight = weight], cl(radio_station) 
summ RADIO_SUBSCRIPTION1946_PCT
local mean_exposure = r(mean)
display `mean_exposure' 

* mean turnout
summ WOMEN_TURNOUT
local mean_turnout = r(mean)
display `mean_turnout' 

* the  number of eligible voters who changed their voting behavior as a result of the introduction of the  newspaper, as a fraction of all those who could have changed their behavior.

local persuasion_rate = 100*`mean_exposure'*_b[RADIO_SUBSCRIPTION1946_PCT]*`mean_turnout'/((1 -  `mu_share_cf'))
display `persuasion_rate' 

ivreg2 VOTE_SHARE ${var_base} (WOMEN_TURNOUT_PCT = wavg_signal_sdunit) [aweight = weight], cl(radio_station) 
ereturn list
	est store tsls_voteshare_2
	local fstat: display %4.2f `e(widstat)'
	local r2c: display %4.2f `e(r2c)'
	estadd local R2C `r2c'
	estadd local FSTAT `fstat'
	estadd local DT "decile bins"
	estadd local PFE "Yes"
	estadd local SFE "Yes"	
	estadd local CFE "Yes"
	estadd local CL "Yes"
	qui summ VOTE_SHARE
	local mu : di %4.2f r(mean)
	estadd local mu `mu'

	
label var RADIO_SUBSCRIPTION1946_SDUNIT "\shortstack[l]{Radio subscription \\ in std.dev. unit}" 
label var wavg_signal_sdunit "\shortstack[l]{Field strength \\ in std.dev. unit}" 	
label var WOMEN_TURNOUT_PCT "\shortstack[l]{Women's turnout \%}" 
	
#delimit;
esttab  rf_voteshare_base ols_voteshare_base tsls_voteshare_base tsls_voteshare_2
        using "${output}/tables/FEMALE-VOTE-SHARE-result-for-revision.tex",
        replace
		compress 
		b(a3) 
		se(a3) 
		star(${setstarlevels})
		noconstant 
		sfmt(%4.2f)
		mtitle("RF" "OLS" "TSLS" "TSLS")
		mgroups("Female Candidate's Vote Share", pattern(1 0 0 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}))
		keep(wavg_signal_sdunit RADIO_SUBSCRIPTION1946_SDUNIT  WOMEN_TURNOUT_PCT)
		order(wavg_signal_sdunit RADIO_SUBSCRIPTION1946_SDUNIT WOMEN_TURNOUT_PCT)
		label
		booktabs
		substitute(\_ _)
	    scalar(
		    "FSTAT First stage F-stat"  
			"CFE Candidate FE"		     
			"PFE Prefecture FE" 
			"SFE Transmitter FE" 
			"DT Distance control"  
			"CL Log N of hh and WWII damage"
			"mu Mean outcome");
#delimit cr  


/*----- Table 1 Balance test -----*/

global var_base "i.NEAR_DIST_DECILEBIN i.radio_station i.prefecture LOG_HH_10000  HEAVY_DAMAGE"

use "${data}/analysis_data_clean.dta", clear
keep if sample_baseline == 1 
sum HH_SUBSCRBE_TOTAL
local total_hh_freq = r(sum)
gen popu_weight     = HH_SUBSCRBE_TOTAL / `total_hh_freq'

label var wavg_signal_sdunit "\shortstack[l]{Residual Radio Signal}" 
reg wavg_signal_sdunit ${var_base}
	predict predict_wavg_signal
	predict double resid_signal, residuals

/* Political variables*/
	
	* Men's turnout in the last election pre WWII
    reg TURNOUT1942  resid_signal  [aweight = popu_weight], cl(radio_station)
	est store balturnout1942
	sum TURNOUT1942 
	local meandep: display %4.2f  `r(mean)'
	estadd local meandep  `meandep'
		
	* The share of votes cast on candidates from the Imperial Rule Assistance Association (Taiseki Yokusan Kai)
	reg VOTE1942_SHARE_YOKUSAN resid_signal   [aweight = popu_weight], cl(radio_station)
	est store balIRAA
	sum VOTE1942_SHARE_YOKUSAN
	local meandep: display %4.2f  `r(mean)'
	estadd local meandep  `meandep'
	
	* Whether or not WSL candidate existed in the district
	reg   D1928_WSLSUPPORT resid_signal  [aweight = popu_weight], cl(radio_station)
	est store balWSL
	sum D1928_WSLSUPPORT 
	local meandep: display %4.2f  `r(mean)'
	estadd local meandep  `meandep'
	
	* The share of votes cast on candidates that Women's Suffrage League supported in 1928 election	
	reg  VOTE1928_SHARE_WSLSUPPORT resid_signal  [aweight = popu_weight] if D1928_WSLSUPPORT == 1, cl(radio_station)
	est store balWSLVOTE
	sum VOTE1928_SHARE_WSLSUPPORT
	local meandep: display %4.2f  `r(mean)'
	estadd local meandep  `meandep'
	
	reg TURNOUT1928 resid_signal   [aweight = popu_weight], cl(radio_station)
	est store balturnout1928
	sum TURNOUT1928
	local meandep: display %4.2f  `r(mean)'
	estadd local meandep  `meandep'


	reg TURNOUT1937  resid_signal   [aweight = popu_weight], cl(radio_station)
	est store balturnout1937
	sum TURNOUT1937
	local meandep: display %4.2f  `r(mean)'
	estadd local meandep  `meandep'

local mymodels balturnout1928 balWSL balWSLVOTE balturnout1942 balIRAA balall
local nmodels : word count `mymodels'
local colcount = `=`nmodels' + 1'

#delimit;
esttab  balturnout1928 balWSL balWSLVOTE balturnout1937 balturnout1942 balIRAA 
    using "${output}/tables/balancing-test-for-revision.tex",
    replace
    booktabs
    b(a3) 
    se(a3) 
    star(${setstarlevels})
	mgroups("Political variables pre WWII" "Political variables during WWII", pattern(1 0 0 0 1 0)
	prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}))
    noconstant 
    sfmt(%4.2f)
    mtitle("\shortstack[c]{Turnout \\ (1928)}" "\shortstack[c]{WSL district \\ (1928)}" "\shortstack[c]{WSL vote share \\ (1928)}" "\shortstack[c]{Turnout \\ (1937)}" "\shortstack[c]{Turnout \\ (1942)}" "\shortstack[c]{IRAA vote share \\ (1942)}" )
    keep(resid_signal)
    label
    nonotes
    substitute(\_ _)
    obslast
    scalar("meandep Mean Dep. Var.")
    prehead("{"
	"\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}"
"\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}lrrrrrr}"
"\multicolumn{7}{c}{\textbf{Panel A}} \\[0.5em]"
"\toprule")
	postfoot("\bottomrule" "\end{tabular*}" "}");
#delimit cr  


/* Demographic variables */

* CBR1935
reg CBR1935 resid_signal  [aweight = popu_weight], cl(radio_station)
ereturn list
est store balCBR1935
sum CBR1935
local meandep: display %4.2f `r(mean)'
estadd local meandep `meandep'

* MFRATIO at birth (the number of boys per girl)
reg  MFRATIO_BIRTH1935 resid_signal [aweight = popu_weight], cl(radio_station)
est store balMFBIRTH1935
sum MFRATIO_BIRTH1935
local meandep: display %4.2f `r(mean)'
estadd local meandep `meandep'

* CRDIVORCE1935
reg CRDIVORCE1935 resid_signal  [aweight = popu_weight], cl(radio_station)
est store balCRDIVORCE1935
sum CRDIVORCE1935
local meandep: display %4.2f `r(mean)'
estadd local meandep `meandep'

* CRMARRIAGE1935
reg CRMARRIAGE1935 resid_signal   [aweight = popu_weight], cl(radio_station)
est store balCRMARRIAGE1935
sum CRMARRIAGE1935
local meandep: display %4.2f `r(mean)'
estadd local meandep `meandep'

* MFRATIO_1940
reg MFRATIO_1940 resid_signal  [aweight = popu_weight], cl(radio_station)
est store balMFRATIO1940
sum MFRATIO_1940
local meandep: display %4.2f `r(mean)'
estadd local meandep `meandep'

* Agriculture labor share
reg LABORSHARE1940_PRIMARY resid_signal  [aweight = popu_weight], cl(radio_station)
est store balPRIM1940
sum LABORSHARE1940_PRIMARY
local meandep: display %4.2f `r(mean)'
estadd local meandep `meandep'


* FEMALE_LABOR_RATE1940
reg FEMALE_LABOR_RATE1940 resid_signal  [aweight = popu_weight], cl(radio_station)
est store balFLP1940
sum FEMALE_LABOR_RATE1940
local meandep: display %4.2f `r(mean)'
estadd local meandep `meandep'

* Sex industrial segregation
reg SEGREGATION1940_H resid_signal [aweight = popu_weight], cl(radio_station)
est store balGEG1940
sum SEGREGATION1940_H
local meandep: display %4.2f `r(mean)'
estadd local meandep `meandep'

local mymodels balMFBIRTH1935 balMFRATIO1940 balPRIM1940  balGEG1940 balFLP1940  balCRMARRIAGE1935 balCBR1935 
local nmodels : word count `mymodels'
local colcount = `=`nmodels' + 1'
display  `colcount'

#delimit ;
esttab balMFBIRTH1935 balMFRATIO1940  balGEG1940 balPRIM1940  balFLP1940  balCRMARRIAGE1935 balCBR1935 
    using "${output}/tables/balancing-test-for-revision-2.tex",
    replace
    booktabs
    b(a3)
    se(a3)
    star(${setstarlevels})
    noconstant
    sfmt(%4.2f)
	mgroups("Demographic and Economic Variables pre/during WWII" "Pre US Occupation Outcomes", pattern(1 0 0 0 1 0 0)
	prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}))
    mtitle("\shortstack[c]{Male to female \\ ratio at birth \\ (1935)}" 
	"\shortstack[c]{Male to female \\ ratio \\ (1940)}" 
	"\shortstack[c]{Gender industrial \\ segregation \\ (1940)}"
	"\shortstack[c]{Primary sector \\ labor share \\  (1940)}" 
	"\shortstack[c]{Female \\ Employment Rate \\ (1940)}"
	"\shortstack[c]{Marriage \\ rate  \\ (1935)}"
	"\shortstack[c]{Birth \\ rate \\ (1935)}")
    keep(resid_signal)
    label
    nonotes
    obslast
    substitute(\_ _)
    scalar("meandep Mean Dep. Var.")
	prehead("{"
	"\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}"
"\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}lrrrrrrrr}"
"\multicolumn{8}{c}{\textbf{Panel B}} \\[0.5em]"
"\toprule")
postfoot("\bottomrule" "\end{tabular*}" "}");
#delimit cr

