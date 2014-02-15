;+
; NAME:
;  irtf_rleg
;
; PURPOSE:
;
;  A convenience function to apply a legendre polynomial
;
; CATEGORY:
;
;  IRTF Analysis
;
; CALLING SEQUENCE:
;
;  out = irtf_rleg(x,c)
;
; INPUTS:
;
;	x: The independant variable
;
;	c: The coefficients for each Legendre term
;	
; OUTPUTS:
;	
;	A legendre polynomial
;
; KEYWORD PARAMETERS:
;	
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-11-2014
;-


function irtf_rleg, x, c
	;return legendre polynomials (m=0)
	;x is indep variable
	;c is coeffs
	nx = n_elements(x)
	nc = n_elements(c)
	res = dblarr(nx)
	for i=0, nc-1 do begin
		res += c[i] * legendre(x,i)
	endfor
	return,res
end
