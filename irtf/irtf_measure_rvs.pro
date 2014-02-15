;+
; NAME:
;  irtf_measure_rvs
;
; PURPOSE:
;
;  Top level script to chew through all IRTF spectra and append RV measurements to the headers
;
; CATEGORY:
;
;  IRTF Analysis
;
; CALLING SEQUENCE:
;
;  irtf_measure_rvs
;
; INPUTS:
;
; OUTPUTS:
;
; KEYWORD PARAMETERS:
;	
; MODIFICATION HISTORY:
;
;  Written by: Ryan Terrien 02-11-2014
;-


pro irtf_measure_rvs

	;file locations on ryan's laptop
	
	files = file_search('~/research/metals/f_analysis/data_annotated/*.fits',count=nf)
	rv_time = systime()
	rv_dataloc = '~/research/metals/f_analysis/rv_data/'+rv_time
	diagloc = '~/research/metals/f_analysis/rv_diag/'+rv_time
	output_loc = '~/research/metals/f_analysis/data_annotated_analyze/'
	
	file_mkdir,rv_dataloc
	file_mkdir,diagloc
	
	;manual limits, to exclude telluricy regions
	lims = transpose([[1.95,2.37],[1.5,1.75],[1.13,1.33],[.96,1.2],[.81,1.02],[.81,.88d]])

	
	for i=0, nf-1 do begin
		diagplot_file = diagloc + '/'+file_basename(files[i])
		
		;Load data
		a = mrdfits(files[i],0,hdr)
		print,files[i]
		
		;Perform the measurements and find the average
		res = irtf_find_shift_spec(a,lims=lims,diagplot_file = diagplot_file)
		rv_final = irtf_combine_rvs(a,res)
		
		;Store the results in the FITS headers
		sxaddpar,hdr,'RV',rv_final[0],'Measured RV (m/s)',before='HISTORY'
		sxaddpar,hdr,'RV_sct',rv_final[1],'RV Scatter (m/s)',before='HISTORY'
		
		mwrfits,res,rv_dataloc+'/'+rv_filename,rv_hdr,/create
		aout = output_loc + file_basename(files[i])
		mwrfits,a,aout,hdr,/create

		;Store the CCF parameters
		rv_filename = file_basename(files[i],'.fits')+'_RV.fits'
		fxbhmake,rv_hdr,n_elements(res)
		sxaddpar,rv_hdr,'RV_time',rv_time,'Time of RV Measurement'
		
	endfor
	
end