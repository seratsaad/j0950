
;-------------------------------------------------------------
;+
; NAME:
;	EXTRACT_UVSPEC 
;
; PURPOSE:                                                 
;	Extract UV spectra from COS and STIS from HST Cycle 18
;	observations and generate ASCII spectra and headers. 
;
; CALLING SEQUENCE:
;       EXTRACT_UVSPEC 
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; NOTES:
;        This code is under construction.   
;
;
; EXAMPLE:
; 
;
; MODIFICATION HISTORY:
;        JCR, 25 July 2016: VERSION 1.00
;           -  Started writing the first version of the code.
;
;-
;-------------------------------------------------------------


function ang_sep,ra1,dec1,ra2,dec2
	RAD_PER_DEG = !pi / 180.            ; radians per degree

    	sra1  = sin(ra1 * RAD_PER_DEG)		
    	cra1  = cos(ra1 * RAD_PER_DEG)		
    	sdec1 = sin(dec1 * RAD_PER_DEG)	
    	cdec1 = cos(dec1 * RAD_PER_DEG)	

	sra2  = sin(ra2 * RAD_PER_DEG)
	cra2  = cos(ra2 * RAD_PER_DEG)
	sdec2 = sin(dec2 * RAD_PER_DEG)
	cdec2 = cos(dec2 * RAD_PER_DEG)


    	csep = cdec1*cdec2*(cra1*cra2 + sra1*sra2) + sdec1*sdec2
    	degsep = acos(csep) / RAD_PER_DEG

	return,degsep
end


pro extract_uvspec,HARDCOPY=HARDCOPY 

	; grab all the COS spectra
	spawn,'ls ./fits/*x1dsum.fits',raw_cos_files

	; get the list of objects that were observed
	readcol,'./hst_targets.dat',format='I,F,A,A',num,z,ra_str,dec_str,/silent	
	bbh_ra  = hms2dec(ra_str)*15.	
	bbh_dec = hms2dec(dec_str)		
	radius  = 0.5/3600.			; spherematch radius in degrees

	; loop over the COS spectra
	for i=0,n_elements(raw_cos_files)-1 do begin
		; read in the fits file
		dat = mrdfits(raw_cos_files[i],0,hdr_info,/silent)
		dat = mrdfits(raw_cos_files[i],1,header,/silent)

		; split them into NUVA and NUVB segments
		wave_long     = (dat.wavelength)[0:(dat.nelem)[0]-1.]	
		flux_long     = (dat.flux)[0:(dat.nelem)[0]-1.]	
		fluxerr_long  = (dat.error)[0:(dat.nelem)[0]-1.]	
		dq_long       = (dat.dq)[0:(dat.nelem)[0]-1.]

		wave_short    = (dat.wavelength)[(dat.nelem)[0]:n_elements(dat.wavelength)-1]	
		flux_short    = (dat.flux)[(dat.nelem)[0]:n_elements(dat.wavelength)-1]	
		fluxerr_short = (dat.error)[(dat.nelem)[0]:n_elements(dat.wavelength)-1]		
		dq_short      = (dat.dq)[(dat.nelem)[0]:n_elements(dat.wavelength)-1]	

		; throw out out-of-bounds pixels
		; these have dq flag values of 128 
		; (http://www.stsci.edu/hst/cos/documents/handbooks/datahandbook/ch2_cos_data8.html#464798)
		;wave_long     = wave_long[where(dq_long ne 128.)]
		;flux_long     = flux_long[where(dq_long ne 128.)] 
		;fluxerr_long  = fluxerr_long[where(dq_long ne 128.)]
		;dq_long       = dq_long[where(dq_long ne 128.)]

		;wave_short    = wave_short[where(dq_short ne 128.)]
		;flux_short    = flux_short[where(dq_short ne 128.)] 
		;fluxerr_short = fluxerr_short[where(dq_short ne 128.)]
		;dq_short      = dq_short[where(dq_short ne 128.)]

		; get the observation date
		regex = 'DATE-OBS'
		date_pos = STREGEX(header,regex)
		datestr = header[where(date_pos eq 0.)]
		obsdate = strmid(datestr,11,4)+strmid(datestr,16,2)+strmid(datestr,19,2)
		utstr   = 'UT '+strmid(datestr,11,4)+'/'+strmid(datestr,16,2)+'/'+strmid(datestr,19,2)

		; get the coordinates 
		regex = 'RA_TARG'
		ra_pos = STREGEX(hdr_info,regex)
		rastr = hdr_info[where(ra_pos eq 0.)]
		ra    = double(strmid(rastr,12,18))

		regex = 'DEC_TARG'
		dec_pos = STREGEX(hdr_info,regex)
		decstr = hdr_info[where(dec_pos eq 0.)]
		dec    = double(strmid(decstr,12,18))

		; get the target name
		regex   = 'TARGNAME'
		obj_pos = STREGEX(hdr_info,regex)
		objstr  = hdr_info[where(obj_pos eq 0.)]
		target  = strmid(objstr,11,23) ; JCR 04/13/22 edited for SDSSJ instead of BBH num name
		m2      = 0 ; JCR 04/13/22 hardcoded for one target with SDSSJ and no BBH match 
		bbhnum  = num[m2]	
		bbhz    = (z[m2])[0]

		; get exposure times
		regex = 'EXPTIME ='
		expt  = float(strmid(header(where(stregex(header,regex) eq 0.)),20,10))

		; match to BBH objects
		; JCR 04/13/22 only one object, no matching needed
		; spherematch doesn't run on Ilarion, so use
		; a klugey workaround that John Ruan wrote
		;degsep  = ang_sep(replicate(ra,n_elements(bbh_ra)),replicate(dec,n_elements(bbh_dec)),bbh_ra,bbh_dec)
		;m2      = where(degsep eq min(degsep))
		;if degsep[m2] ge radius then print, 'Heads up, yo!  The closest positional match is larger than the search radius.'
		;bbhnum  = num[m2]	
		;bbhz    = (z[m2])[0]

		; print an error if the target and object
		; matching produced different results
		if string(bbhnum,format='(I03)') ne target then begin
			print, 'The coordinate and target matching produced different results.'
		endif

		; convert the flux units to 1d-17 erg/s/cm^2/A
		flux_long     = (1d17)*flux_long
		fluxerr_long  = (1d17)*fluxerr_long
		flux_short    = (1d17)*flux_short
		fluxerr_short = (1d17)*fluxerr_short

		; convert to rest wavelengths
		wave_long     = wave_long/(1.+bbhz)
		wave_short    = wave_short/(1.+bbhz)

		; generate the output filenames
		file_hdr0     = './ascii/bbh'+strtrim(string(bbhnum,format='(I03)'),2)+'_'+obsdate+'_COS_info.hdr'
		file_hdr1     = './ascii/bbh'+strtrim(string(bbhnum,format='(I03)'),2)+'_'+obsdate+'_COS_dat.hdr'
		file_coss     = './ascii/bbh'+strtrim(string(bbhnum,format='(I03)'),2)+'_'+obsdate+'_COSs.dat'
		file_coss_2ca = './2ca/bbh'+strtrim(string(bbhnum,format='(I03)'),2)+'_'+obsdate+'_COSs.2ca'
		file_coss_2cu = './2cu/bbh'+strtrim(string(bbhnum,format='(I03)'),2)+'_'+obsdate+'_COSs.2cu'
		file_cosl     = './ascii/bbh'+strtrim(string(bbhnum,format='(I03)'),2)+'_'+obsdate+'_COSl.dat'
		file_cosl_2ca = './2ca/bbh'+strtrim(string(bbhnum,format='(I03)'),2)+'_'+obsdate+'_COSl.2ca'
		file_cosl_2cu = './2cu/bbh'+strtrim(string(bbhnum,format='(I03)'),2)+'_'+obsdate+'_COSl.2cu'
		plotfile_coss = './plots/bbh'+strtrim(string(bbhnum,format='(I03)'),2)+'_'+obsdate+'_COSs.eps'
		plotfile_cosl = './plots/bbh'+strtrim(string(bbhnum,format='(I03)'),2)+'_'+obsdate+'_COSl.eps'

		; set the filter
		filter = bbhz le 0.4 ? 'G140L' : 'G230L'

		; write out all the files
		openw,1,file_hdr0,width=1600
			printf,1,hdr_info
		close,1
		openw,1,file_hdr1,width=1600
			printf,1,header
		close,1
		openw,1,file_coss,width=800
			printf,1,'BBH '+strtrim(string(bbhnum,format='(I03)'),2)               
			printf,1,strtrim(utstr,2)+' COS/FUV '+strtrim(filter,2)                                                            
			printf,1,' '+strtrim(string(expt),2)+' second exposure'
			printf,1,strtrim(string(n_elements(wave_short),format='(I)'),2)+' points'
			printf,1,''
			printf,1,'     Wavelength (A)            Flux (erg/s/cm^2/A)            Flux Err erg/s/cm^2/A()'
			printf,1,'     --------------            ------------------             -----------------------'
			for j=0,n_elements(wave_short)-1 do printf,1,wave_short[j],'        ',flux_short[j],'               ',fluxerr_short[j]
		close,1
		openw,1,file_coss_2ca,width=800
			printf,1,'BBH '+strtrim(string(bbhnum,format='(I03)'),2)               
			printf,1,strtrim(utstr,2)+' COS/FUV '+strtrim(filter,2)                                                            
			printf,1,' '+strtrim(string(expt),2)+' second exposure'
			printf,1,strtrim(string(n_elements(wave_short),format='(I)'),2)+' points'
			printf,1,''
			printf,1,'     Wavelength (A)            Flux (erg/s/cm^2/A)'
			printf,1,'     --------------            ------------------ '
			for j=0,n_elements(wave_short)-1 do printf,1,wave_short[j],'        ',flux_short[j]
		close,1
		openw,1,file_cosl,width=800
			printf,1,'BBH '+strtrim(string(bbhnum,format='(I03)'),2)               
			printf,1,strtrim(utstr,2)+' COS/NUV '+strtrim(filter,2)                                                            
			printf,1,' '+strtrim(string(expt),2)+' second exposure'
			printf,1,strtrim(string(n_elements(wave_long),format='(I)'),2)+' points'
			printf,1,''
			printf,1,'     Wavelength (A)            Flux (erg/s/cm^2/A)            Flux Err erg/s/cm^2/A()'
			printf,1,'     --------------            ------------------             -----------------------'
			for j=0,n_elements(wave_long)-1 do printf,1,wave_long[j],'        ',flux_long[j],'               ',fluxerr_long[j]
		close,1
		openw,1,file_cosl_2ca,width=800
			printf,1,'BBH '+strtrim(string(bbhnum,format='(I03)'),2)               
			printf,1,strtrim(utstr,2)+' COS/NUV '+strtrim(filter,2)                                                            
			printf,1,' '+strtrim(string(expt),2)+' second exposure'
			printf,1,strtrim(string(n_elements(wave_long),format='(I)'),2)+' points'
			printf,1,''
			printf,1,'     Wavelength (A)            Flux (erg/s/cm^2/A)'
			printf,1,'     --------------            ------------------ '
			for j=0,n_elements(wave_long)-1 do printf,1,wave_long[j],'        ',flux_long[j]
		close,1

		; make a temporary command file for 
		; the 2ca-->2cu conversion
		; JCR 04/13/2022 this isn't working because 2ca_2cu.x throws a bus error
;		openw,1,'scratch.com',width=400
;			printf,1,'/Users/runnojc1/Software/mce/spec/2ca_2cu.x << inputs'
;			printf,1,file_coss_2ca
;			printf,1,file_coss_2cu
;			printf,1,'inputs'
;			printf,1,'/Users/runnojc1/Software/mce/spec/2ca_2cu.x << inputs'
;			printf,1,file_cosl_2ca
;			printf,1,file_cosl_2cu
;			printf,1,'inputs'
;		close,1
;		spawn,'chmod a+x scratch.com'
;		spawn,'./scratch.com'



		; make a figure
		;xmin = min(wave_short) le min(wave_long) ? min(wave_short) : min(wave_long)
		;xmax = max(wave_short) ge max(wave_long) ? max(wave_short) : max(wave_long)
		;ymin = 0.4*min(flux_short[where(dq_short eq 0)]) le 0.4*min(flux_long[where(dq_long eq 0)]) ? 0.4*min(flux_short[where(dq_short eq 0)]) : 0.4*min(flux_long[where(dq_long eq 0)]) 
		;ymax = 1.2*max(flux_short[where(dq_short eq 0)]) ge 1.2*max(flux_long[where(dq_long eq 0)]) ? 1.2*max(flux_short[where(dq_short eq 0)]) : 1.2*max(flux_long[where(dq_long eq 0)]) 

		;if keyword_set(hardcopy) then begin
		;	cgps_open,plotfile_stis,/encapsulated,xsize=7.5,ysize=6
		;	; set the font to make Mike happy
		;	DEVICE, SET_FONT = 'Times-Roman'
		;endif
		;	ang  = cgSymbol('Angstrom')
		;	lam  = cgSymbol('lambda')
		;	plot,wave_short,flux_short,/nodata,psym=10,ystyle=1,xstyle=1,ytit='f!d'+lam+'!n [10!u-17!n erg s!u-1!n cm!u-2!n '+ang+'!u-1!n]',xtit='Rest Wavelength ['+ang+']',charsize=2,xra=[xmin,xmax],yra=[ymin,ymax]
		;	oplot,wave_short,flux_short,psym=10,color=cgcolor('red')
		;	oplot,wave_long,flux_long,psym=10
		;if keyword_set(hardcopy) then begin
		;	cgps_close
		;	!p.font=-1
		;endif

		; make a figure
		xmin = min(wave_long[where(dq_long ne 128.)])-50.
		xmax = max(wave_long[where(dq_long ne 128.)])+50.
		ymin = 0.4*min(flux_long[where(dq_long eq 0)])
		ymax = 1.2*max(flux_long[where(dq_long eq 0)])

		if keyword_set(hardcopy) then begin
			cgps_open,plotfile_cosl,/encapsulated,xsize=7.5,ysize=6
			; set the font to make Mike happy
			DEVICE, SET_FONT = 'Times-Roman'
		endif
			ang  = cgSymbol('Angstrom')
			lam  = cgSymbol('lambda')
			plot,wave_long,flux_long,psym=10,ystyle=1,xstyle=1,ytit='f!d'+lam+'!n [10!u-17!n erg s!u-1!n cm!u-2!n '+ang+'!u-1!n]',xtit='Rest Wavelength ['+ang+']',charsize=2,xra=[xmin,xmax],yra=[ymin,ymax]
			oploterror,wave_long,flux_long,fluxerr_long,psym=3,color=cgcolor('gray'),errcolor=cgcolor('gray'),/nohat
			oplot,wave_long,flux_long,psym=10
			plot,wave_long,flux_long,psym=10,ystyle=1,xstyle=1,ytit='f!d'+lam+'!n [10!u-17!n erg s!u-1!n cm!u-2!n '+ang+'!u-1!n]',xtit='Rest Wavelength ['+ang+']',charsize=2,xra=[xmin,xmax],yra=[ymin,ymax],/noerase
		if keyword_set(hardcopy) then begin
			cgps_close
			!p.font=-1
		endif	
		xmin = min(wave_short[where(dq_short eq 0)]) gt 800. ? min(wave_short[where(dq_short eq 0)]) : 800.
		xmax = max(wave_short[where(dq_short eq 0)])+50.
		ymin = 0.4*min(flux_short[where(dq_short eq 0)])
		ymax = 1.0*max(flux_short[where(dq_short eq 0)])

		if keyword_set(hardcopy) then begin
			cgps_open,plotfile_coss,/encapsulated,xsize=7.5,ysize=6
			; set the font to make Mike happy
			DEVICE, SET_FONT = 'Times-Roman'
		endif
			ang  = cgSymbol('Angstrom')
			lam  = cgSymbol('lambda')
			plot,wave_short,flux_short,psym=10,ystyle=1,xstyle=1,ytit='f!d'+lam+'!n [10!u-17!n erg s!u-1!n cm!u-2!n '+ang+'!u-1!n]',xtit='Rest Wavelength ['+ang+']',charsize=2,xra=[800,xmax],yra=[ymin,ymax]
			oploterror,wave_short,flux_short,fluxerr_short,psym=3,color=cgcolor('gray'),errcolor=cgcolor('gray'),/nohat
			oplot,wave_short,flux_short,psym=10
			plot,wave_short,flux_short,psym=10,ystyle=1,xstyle=1,ytit='f!d'+lam+'!n [10!u-17!n erg s!u-1!n cm!u-2!n '+ang+'!u-1!n]',xtit='Rest Wavelength ['+ang+']',charsize=2,xra=[800,xmax],yra=[ymin,ymax],/noerase
		if keyword_set(hardcopy) then begin
			cgps_close
			!p.font=-1
		endif	
	
	endfor
	
		
	

stop,'.cont to continue'
end
