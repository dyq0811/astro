;+
; NAME:
;  irtf_mkbins
;
; PURPOSE:
;
;  Find bin edge values as midpoints between pixels
;
; CATEGORY:
;
;  IRTF Analysis
;
; CALLING SEQUENCE:
;
;  out = irtf_mkbins(inp)
;
; INPUTS:
;
;	inp: Array to interpret, monotonic increasing or decreasing
;	
; OUTPUTS:
;	
;	Array of bin edge estimates. One more element than input (1|2|3)
;
; KEYWORD PARAMETERS:
;	
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-11-2014
;-


function irtf_mkbins, inp

arr = double(inp)

shift = arr[1:*] - arr
bin1 = arr + shift/2.d
last = n_elements(arr) - 1
bina = arr[0] - shift[0]/2.d
binz = arr[last] + shift[last-1]/2.d

bin = [bina,bin1,binz]

return,bin

end

