;+
; NAME:
;  irtf_legendre_norm
;
; PURPOSE:
;
;  Normalize a spectrum using legendre polynomials
;
; CATEGORY:
;
;  IRTF Analysis
;
; CALLING SEQUENCE:
;
;  out = irtf_legendre_norm(wls, fls, ers, boxcar_param, legendre_order)
;
; INPUTS:
;
;	wls,fls,ers: The input spectrum
;
;	boxcar_param: The size of the boxcar to use for max/medians
;
;	legendre_order: The order of the Legendre polynomial
;	
; OUTPUTS:
;	
;	A normalized spectrum
;
; KEYWORD PARAMETERS:
;
;	buffer_ends: The number of pixels on each end of the spectrum to replace with a median
;
;	filter_spikes: Do a rough exclusion of anomalously high-flux pixels
;
;	diag: A structure containing parameters of the normalization for plotting/etc
;	
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-11-2014
;-


function irtf_legendre_norm, wls, fls, ers, boxcar_param, legendre_order,buffer_ends = buffer_ends, filter_spikes = filter_spikes, diag = diag

	;If the boxcar median is needed:
	if n_elements(buffer_ends) ne 0 or n_elements(filter_spikes) ne 0 then $
		irtf_boxcar_med,wls,fls,boxcar_param,wlm0,ms0,/keep_nan
		
	;Filter spikes by masking anything > 2 * [sliding STDEV] + [sliding MED]
	if n_elements(filter_spikes) ne 0 then begin
		fls_orig = fls
		for i=0, n_elements(fls)-1 do begin
			i1 = 0. > (i - 20)
			i2 = (n_elements(fls)-1) < (i + 20)
			sd = stddev(fls_orig[i1:i2],/nan)
			if fls_orig[i] gt (2.*sd + ms0[i]) then fls[i] = ms0[i]
		endfor
	endif
	
	;Buffer the ends with the sliding median
	if n_elements(buffer_ends) ne 0 then begin
		fls1 = fls
		fls1[-1*buffer_ends:*] = ms0[-1*buffer_ends:*]
		fls1[0:buffer_ends] = ms0[0:buffer_ends]
	endif else begin
		fls1 = fls
	endelse
	
	;Find the sliding max
	irtf_boxcar_max,wls,fls1,boxcar_param,wlm1,ms1,/keep_nan	
	
	;Select out the finite pixels
	good = where(finite(fls))
	ms = double(ms1[good])
	wlm = double(wlm1[good])
	ers = double(ers[good])
	
	;A small number of spectra have NaN errs but finite fluxes. The svdfit does not
	;converge if fed NaNs, so replace these with dummy errors for the RV purpose
	naners = where(~finite(ers),nnan)
	if nnan ne 0 then ers[naners] = median(ers,/double)
	
	;Set up and perform the Legendre fit
	wl_fit = ((wlm - min(wlm)) / (max(wlm) - min(wlm)))*2d - 1d
	fl_fit = svdfit(wl_fit,ms,legendre_order,measure_errors = ers,/legendre,/double,yfit=fl_nvals_nonan)
	wl_all = ((wls - min(wls)) / (max(wls) - min(wls)))*2d - 1d
	fl_nvals = irtf_rleg(wl_all,fl_fit)
	fln = fls / fl_nvals
	
	;Populate diag
	diag = {wl:wls, fl:fls, ms:ms, fl_nvals:fl_nvals}
	
	return,fln
end