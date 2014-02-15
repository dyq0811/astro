;+
; NAME:
;  irtf_boxcar_max
;
; PURPOSE:
;
;  Perform a sliding boxcar maximum
;
; CATEGORY:
;
;  IRTF Analysis
;
; CALLING SEQUENCE:
;
;  irtf_boxcar_max, wl, fl, pix, wl_out, ms_out
;
; INPUTS:
;
;	wl, fl: Input spectrum 
;
;	pix: The size of the boxcar
;
;	wl_out,ms_out: The output max-spectrum
;	
; OUTPUTS:
;	
;	An array with the same size as the number of finite input elements
;
; KEYWORD PARAMETERS:
;
;	keep_nan: Do not delete the NaN elements, leave the array as its fully size
;	
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-11-2014
;-


pro irtf_boxcar_max, wl, fl, pix, wl_out, ms_out, keep_nan = keep_nan
	;boxcar smooth a function
	nel = n_elements(fl)
	ms = dblarr(nel)
	for i=0, nel-1 do begin
		if ~finite(fl[i]) then ms[i] = !values.f_nan else begin
			temp = max(fl[0 > (i - pix) : (nel-1) < (i + pix)],/nan)
			ms[i] = temp
		endelse
	endfor
	if keyword_set(keep_nan) then begin
		wl_out = wl
		ms_out = ms
		return
	endif
	g = where(finite(ms))
	wl_out = wl[g]
	ms_out = ms[g]
end
