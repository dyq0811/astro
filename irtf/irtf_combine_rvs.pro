;+
; NAME:
;  irtf_combine_rvs
;
; PURPOSE:
;
;  Combine RV measurements from all orders of an IRTF spectrum
;
; CATEGORY:
;
;  IRTF Analysis
;
; CALLING SEQUENCE:
;
;  irtf_combine_rvs( spec, rv_struct)
;
; INPUTS:
;
;	spec: A (1024,3,n_orders) array from IRTF
;
;	rv_struct: A output RV structure with n elements, each having an .rv tag
;	
; OUTPUTS:
;	
;	A two-element array with the combined RV (weighted by SN) and unweighted scatter
;
; KEYWORD PARAMETERS:
;	
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-11-2014
;-

function irtf_combine_rvs, spec, rv_struct

	;initialize the variables
	n_orders = (size(spec))[3]
	sns = dblarr(n_orders)
	rvs = dblarr(n_orders)
	
	; pull out the rvs and s/n values
	for i=0, n_orders-1 do begin
		sns[i] = median(spec[*,1,i]/spec[*,2,i])
		rvs[i] = rv_struct[i].rv
	endfor
	
	; calculate the weighted mean
	sntot = total(sns)
	
	weights = sns / sntot
	
	outmean = total(weights * rvs,/double)
	
	;for now just report the standard deviation
	outscatter = stddev(rvs)
	
	out = [outmean, outscatter]
	
	return, out
end