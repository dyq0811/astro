;+
; NAME:
;  irtf_resample_log
;
; PURPOSE:
;
;  Resample a spectrum onto log(wl), conserving flux
;
; CATEGORY:
;
;  IRTF Analysis
;
; CALLING SEQUENCE:
;
;  out = irtf_resample_log(orig_wl, orig_fl, new_wl)
;
; INPUTS:
;
;	orig_wl, orig_fl: The original spectrum. Don't feed it NaNs. No req. for even WL input.
;
;	new_wl: The new wl array (regularly spaced in log(wl))
;	
; OUTPUTS:
;	
;	An array of energies. The output quantities are no longer /wavelength
;
; KEYWORD PARAMETERS:
;	
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-11-2014
;-


function irtf_resample_log, orig_wl, orig_fl, new_wl

fl_out = dblarr(n_elements(new_wl))
sides = irtf_mkbins(new_wl) ;in LOG
left = exp(sides[0:-2]) ;NOTLOG
right = exp(sides[1:*]) ;NOTLOG

orig_sides = irtf_mkbins(orig_wl)
orig_left = orig_sides[0:-2]
orig_right = orig_sides[1:*]

for i=0, n_elements(new_wl)-1 do begin
	;for each pixel in the new log wl array
	;go to linear space and figure out which pixels are btwn left and right sides
	;pixels entirely within limits
	inds = where(orig_left ge left[i] and orig_right le right[i],ni)
	lind = where(orig_right gt left[i] and orig_left lt left[i],nlind)
	rind = where(orig_right gt right[i] and orig_left lt right[i],nrind)
	if nlind gt 1 or nrind gt 1 then message,'this cant happen'
	;if there are multiple pixels
	if ni ge 1 then begin
		tmid = orig_wl[inds] * orig_fl[inds]
		if nlind gt 0 then tl = (orig_right[lind] - left[i])*orig_fl[lind] else tl = 0.d
		if nrind gt 0 then tr = (right[i] - orig_left[rind])*orig_fl[rind] else tr = 0.d
		tot = tl + total(tmid) + tr
	endif
	;if they are together in one pixel
	if ni eq 0 and lind eq rind then begin
		tot = (right[i] - left[i])*orig_fl[lind]
	endif
	;if there is a pixel boundary between them
	if ni eq 0 and lind ne rind then begin
		if nlind gt 0 then tl = (orig_right[lind] - left[i])*orig_fl[lind] else tl = 0.d
		if nrind gt 0 then tr = (right[i] - orig_left[rind])*orig_fl[rind] else tr = 0.d
		tot = tl + tr
	endif
	fl_out[i] = tot
endfor

return, fl_out

end