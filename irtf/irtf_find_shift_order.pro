;+
; NAME:
;  irtf_find_shift_order
;
; PURPOSE:
;
;  Find the RV shift for a single continuous region of IRTF spectrum
;
; CATEGORY:
;
;  IRTF Analysis
;
; CALLING SEQUENCE:
;
;  irtf_find_shift_order( wl, fl, er, wl_temp, fl_temp, er_temp)
;
; INPUTS:
;
;	wl,fl,er: The input spectrum, with no gaps or NaNs
;
;	wl,fl,er_temp: The template spectrum, with no gaps or NaNs
;	
; OUTPUTS:
;	
;	A structure with the following:
;		rv: The final measured RV shift in m/s
;		rv1: The initial center RV shift in m/s
;		ccf1_even/odd: The initial fit parameters for the odd/even# of CCF points
;		ccf2_even/odd: The final fit parameters for the odd/even# of CCF points
;
; KEYWORD PARAMETERS:
;	
;	np1: Number of parameters in the first RV fit (4)
;
;	np2: Number of parameters in the second RV fit (4)
;
;	range1: Range of RV shifts for first RV fit (-500 to 500km/s) in m/s
;
;	range2: Range of RV shifts for second RV fit (-200 to 200km/s) in m/s
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-11-2014
;-


function irtf_find_shift_order, wl, fl, er, wl_temp, fl_temp, er_temp, np1 = np1, np2 = np2, range1 = range1, range2 = range2, diag = diag

	;Set defaults and intialize the output
	
	cm = 299792458d ;meter/sec
	
	if n_elements(np1) eq 0 then np1 = 5
	if n_elements(np2) eq 0 then np2 = 5
	if n_elements(range1) eq 0 then range1 = [-5d5,5d5]
	if n_elements(Range2) eq 0 then range2 = [-2d5,2d5]

	output = {rv:!values.f_nan, rv1:!values.f_nan, ccf1_odd:ptr_new(), ccf1_even:ptr_new(), ccf2_odd:ptr_new(),ccf2_even:ptr_new()}
	
	nn = n_elements(wl)
	
	output.ccf1_odd = ptr_new(/allocate_heap)
	output.ccf2_odd = ptr_new(/allocate_heap)
	output.ccf1_even = ptr_new(/allocate_heap)
	output.ccf2_even = ptr_new(/allocate_heap)
	
	wl_spec_1 = double(wl)
	fl_spec_1 = double(fl)
	er_spec_1 = double(er)
	wl_temp_1 = double(wl_temp)
	fl_temp_1 = double(fl_temp)
	er_temp_1 = double(er_temp)
	
	;Find overlap between input spectrum and template
	lims = [min(wl_temp_1) > min(wl_spec_1), max(wl_temp_1) < max(wl_spec_1)]
	
	;Construct a WL array in logspace that has even bins
	loglims = alog(lims)
	ran = max(loglims) - min(loglims)
	npix = floor(ran / alog(1.d + 1.d4/cm))
	pix = alog(1.d + 1.d4/cm) ; Define a bin to be 10km/s
	reg_lwl = dindgen(npix)*pix + loglims[0]
	
	;Select appropriate regions of the spectrum and template
	;If there are NaNs in the spectrum, get rid of them before resampling, otherwise CCF
	;is entirely NaNs
	i_spec = where(wl_spec_1 ge lims[0] and wl_spec_1 le lims[1] and finite(fl_spec_1))
	i_temp = where(wl_temp_1 ge lims[0] and wl_temp_1 le lims[1])
	
	wl_spec_2 = wl_spec_1[i_spec]
	fl_spec_2 = fl_spec_1[i_spec]
	er_spec_2 = er_spec_1[i_spec]
	wl_temp_2 = wl_temp_1[i_temp]
	fl_temp_2 = fl_temp_1[i_temp]
	er_temp_2 = er_temp_1[i_temp]
	
	;Do the log-resample
	fl_log = irtf_resample_log(wl_spec_2,fl_spec_2,reg_lwl)
	fl_temp_log = irtf_resample_log(wl_temp_2,fl_temp_2,reg_lwl)
	
	;Normalize the log spectra
	fl_log_norm = irtf_legendre_norm(reg_lwl,fl_log,er_spec_2,20,10,buffer_ends=50,/filter_spikes, diag = ndiag)
	fl_temp_log_norm = irtf_legendre_norm(reg_lwl,fl_temp_log,er_temp_2,20,10,buffer_ends=50,/filter_spikes)
	
	;trim the end pixels
	fl_log_norm = fl_log_norm[1:-2]
	fl_temp_log_norm = fl_temp_log_norm[1:-2]

	;CCF round 1, use parameters defined in input
	
	nlags_1 = (range1[1] - range1[0]) / (pix * 3d8)
	if nlags_1 mod 2 eq 0 then nlags_1 += 1 ;Make sure lags are odd
	lags_odd1 = indgen(nlags_1) - floor(nlags_1/2.)
	xc_result_odd1 = c_correlate(fl_log_norm,fl_temp_log_norm,lags_odd1,/double)
	;now do the even, figure out which side to add the point to and cut a pixel from the odd one
	b11 = xc_result_odd1[0]
	b21 = xc_result_odd1[-1]
	if b11 lt b21 then begin
		xc_result_even1 = xc_result_odd1[1:*] 
		lags_even1 = lags_odd1[1:*]
	endif else begin
		xc_result_even1 = xc_result_odd1[0:-2]
		lags_even1 = lags_odd1[0:-2]
	endelse

	;fit both peaks with gaussians
	fit_odd1 = mpfitpeak(lags_odd1,xc_result_odd1,fit_result_odd1,nterms=np1)
	fit_even1 = mpfitpeak(lags_even1,xc_result_even1,fit_result_even1,nterms=np1)

	shift_odd1 = fit_result_odd1[1]
	shift_even1 = fit_result_even1[1]

	shift_avg1 = shift_odd1 / (1. + 2.*(abs(shift_odd1) - abs(shift_even1)))
	
	;Save results
	*output.ccf1_odd = fit_result_odd1
	*output.ccf1_even = fit_result_even1
	dlogwl1 = shift_avg1 * pix
	dv1 = -1.d * (exp(dlogwl1) - 1.d) * cm
	output.rv1 = dv1

	;CCF round 2, use parameters defined in input
	nlags_2 = (range2[1] - range2[0]) / (pix * 3d8)
	if nlags_2 mod 2 eq 0 then nlags_2 += 1
	lags_odd = indgen(nlags_2) - floor(nlags_2/2.) + shift_avg1

	xc_result_odd = c_correlate(fl_log_norm,fl_temp_log_norm,lags_odd,/double)
	;now do the even, figure out which side to add the point to and cut a pixel from the odd one
	b1 = xc_result_odd[0]
	b2 = xc_result_odd[-1]
	if b1 lt b2 then begin
		xc_result_even = xc_result_odd[1:*] 
		lags_even = lags_odd[1:*]
	endif else begin
		xc_result_even = xc_result_odd[0:-2]
		lags_even = lags_odd[0:-2]
	endelse

	;fit both peaks with gaussians
	fit_odd = mpfitpeak(lags_odd,xc_result_odd,fit_result_odd,nterms=np2)
	fit_even = mpfitpeak(lags_even,xc_result_even,fit_result_even,nterms=np2)

	shift_odd = fit_result_odd[1]
	shift_even = fit_result_even[1]

	shift_avg = shift_odd / (1. + 2.*(abs(shift_odd) - abs(shift_even)))
	
	;store results
	*output.ccf2_odd = fit_result_odd
	*output.ccf2_even = fit_result_even

	;convert pixel shift to wavelength shift
	dlogwl = shift_avg * pix

	;convert to velocity
	dv = -1.d * (exp(dlogwl) - 1.d) * cm
	
	output.rv = dv

	;Populate diag	
	diag = {wl:wl, fl:fl, lwl:reg_lwl[1:-2], lfl_norm:fl_log_norm[1:-2], lags1:lags_odd1, ccf1:xc_result_odd1, fit1:fit_odd1, lags2:lags_odd, ccf2:xc_result_odd, fit2:fit_odd, norm_wl:ndiag.wl, norm_fl:ndiag.fl, norm_fl_nvals:ndiag.fl_nvals, norm_ms:ndiag.ms}
	
	
	return,output
		
end