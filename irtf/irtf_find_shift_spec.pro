;+
; NAME:
;  irtf_find_shift_spec
;
; PURPOSE:
;
;  Find the RV shift for an IRTF spectrum based on all the orders
;
; CATEGORY:
;
;  IRTF Analysis
;
; CALLING SEQUENCE:
;
;  irtf_find_shift_spec( spec, lims = lims)
;
; INPUTS:
;
;	spec: A (1024,3,n_orders) array from IRTF
;	
; OUTPUTS:
;	
;	A structure array with the following:
;		rv: The final measured RV shift in m/s
;		rv1: The initial center RV shift in m/s
;		ccf1_even/odd: The initial fit parameters for the odd/even# of CCF points
;		ccf2_even/odd: The final fit parameters for the odd/even# of CCF points
;
; KEYWORD PARAMETERS:
;	
;	lims = A (2, n_orders) array giving the limits of good spectrum for each
;
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-11-2014
;-

function irtf_find_shift_spec, spec, lims = lims, diagplot_file = diagplot_file
	
	;Set up defaults, initialize outputs
	
	n_orders = (size(spec))[3]
	
	if n_elements(lims) eq 0 then begin
		lims = dblarr(2,n_orders)
		lims[1,*] = 10d
	endif

	output_1 = {rv:!values.f_nan, rv1:!values.f_nan, ccf1_odd:ptr_new(), ccf1_even:ptr_new(), ccf2_odd:ptr_new(),ccf2_even:ptr_new()}
	
	output = replicate(output_1,n_orders)
	
	diags = ptrarr(n_orders,/allocate_heap)
		
	;get template array. Only using one for now.
	info = get_login_info()
	hostname = info.machine_name
	case hostname of
		'mlb.astro.psu.edu': template_file = 'M1.5V_HD36395.fits'
		else: template_file = '/Users/ryan/research/metals/irtf_sl/M_fits_091201/M1.5V_HD36395.fits'
	endcase
	template = mrdfits(template_file,0,template_head,/silent)
	
	;Template is not separated by order so it can be outside the order loop
	wl_temp1 = reform(double(template[*,0]))
	fl_temp1 = reform(double(template[*,1]))
	er_temp1 = reform(double(template[*,2]))
	
	;For each order, find the RV measurement parameters
	for i=0, n_orders-1 do begin
		;Extract the spectra
		wl1 = reform(double(spec[*,0,i]))
		fl1 = reform(double(spec[*,1,i]))
		er1 = reform(double(spec[*,2,i]))
		
		;Only take parts of the spectra with finite values
		wlreal = where(finite(fl1))
		wl1r = wl1[wlreal]
		
		;Select which limit is appropriate to use (a small number of spectra have >1
		;aperture and you can't count on [*,0,i] = order i
		limsy = lims[*,0] lt median(wl1r) and lims[*,1] gt median(wl1r)
		limsind = where(limsy,nli)
		if nli eq 0 then stop
		limsind = limsind[0] ; The last two orders overlap, so it doesn't matter which
		
		;Find where spectra are inside manual limits
		xx = where(wl1 ge lims[limsind,0] and wl1 le lims[limsind,1],nxx)
		xx_temp = where(wl_temp1 ge lims[limsind,0] and wl_temp1 le lims[limsind,1],nxx_temp)
		if nxx eq 0 or nxx_temp eq 0 then stop
		
		wl2 = wl1[xx]
		fl2 = fl1[xx]
		er2 = er1[xx]
		wl_temp2 = wl_temp1[xx_temp]
		fl_temp2 = fl_temp1[xx_temp]
		er_temp2 = er_temp1[xx_temp]
		
		;Find the RV shift info for this order
		out = irtf_find_shift_order(wl2, fl2, er2, wl_temp2, fl_temp2, er_temp2, diag = diag)
		
		;Store the result
		output[i] = out
		
		*diags[i] = diag
			
	endfor
	
	;If diagnostic plots are being made then construct them in the appropriate place
	if n_elements(diagplot_file) ne 0 then begin
		specfile = file_dirname(diagplot_file)+'/'+file_basename(diagplot_file,'.fits')+'_rvspec'
		ccffile = file_dirname(diagplot_file)+'/'+file_basename(diagplot_file,'.fits')+'_rvccf'

		; Show the spectrum and the normalization
		openpps,specfile
		!p.multi = [0,3,2]
		for i=0, n_orders-1 do begin
			cgplot,(*diags[i]).norm_wl,(*diags[i]).norm_fl,chars=1,thick=1.,/xs,/ys
			cgplot,(*diags[i]).norm_wl,(*diags[i]).norm_ms,color='blue',thick=1,/overplot
			cgplot,(*diags[i]).norm_wl,(*diags[i]).norm_fl_nvals,thick=1,color='red',/overplot
		endfor
		closepps
		
		; Show the CCF fits for each order
		openpps,ccffile
		!p.multi = [0,3,2]
		for i=0, n_orders-1 do begin
			cgplot,(*diags[i]).lags1,(*diags[i]).ccf1,ps=1.,chars=1,thick=1,/xs
			cgplot,(*diags[i]).lags1,(*diags[i]).fit1,color='red',thick=1,/overplot
			cgplot,(*diags[i]).lags2,(*diags[i]).ccf2,ps=6,thick=1,/overplot,color='blue'
			vline,(*output[i].ccf2_odd)[1],lines=2
		endfor
		closepps
	endif
		
	
	return, output
	
end