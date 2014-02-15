function irtf_rv_shift, spec, rv

	cm = 299792458d ;meter/sec

	;rv in m/s

	out = spec
	beta = rv/cm
	out[*,0,*] = out[*,0,*] * sqrt( (1d + beta) / (1d - beta) ) 
	
	return, out
end