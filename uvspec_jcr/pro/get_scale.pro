
;-------------------------------------------------------------
;+
; NAME:
;       GET_SCALE 
;
; PURPOSE:                                                 
;       Explore methodology for comparing spectra from different
;       emission lines.
;
; CALLING SEQUENCE:
;       get_scale,[/hardcopy,/dolines/slow] 
;
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; NOTES:
;
;
; EXAMPLE:
; 
;
; MODIFICATION HISTORY:
;        JCR, 21 May 2021: VERSION 1.00
;           -  Started writing the first version of the code.
;
;-
;-------------------------------------------------------------


function get_fwzm,wave,flux,pcent
;-------------------------------------------------------------
;+
; PURPOSE:                                                 
;	Numerically calculate the width of a line profile at
;       the percent of peak given as an input.
;
;	fwzm     = get_fwzm(wave,flux,pcent)	
;
; MODIFICATION HISTORY:
;	 JCR, 2 June 2021: VERSION 1.00
;           -  Started writing the first version of the code.
;
;-
;-------------------------------------------------------------
	c        = 2.99792458d5	; km/s
	pflux    = max(flux)
	cwave    = wave[where(flux eq pflux)]

	use      = where(flux ge pcent*pflux)
	wavelo   = min(wave[use]) 
	wavehi   = max(wave[use]) 
	fwzm = c * (wavehi-wavelo)/cwave

	return, fwzm
end

function line, x, p
;-------------------------------------------------------------
;+
; PURPOSE:                                                 
;	Line function for MPFIT.
;
;	y = b+m*x	
;
; MODIFICATION HISTORY:
;	 JCR, 2 June 2021: VERSION 1.00
;           -  Started writing the first version of the code.
;
;-
;-------------------------------------------------------------
  b = p[0]   ; Intercept
  m = p[1]   ; Slope
  f = b + m*x
  return, f
end

pro get_scale,HARDCOPY=HARDCOPY,DOLINES=DOLINES,SLOW=SLOW 

	; set some constants
	c       = 2.99792458d5	; km/s
	w0_lya  = 1215.67
	w0_civ  = 1549.48
	w0_mgii = 2799.94
	w0_hb   = 4862.721
	w0_ha   = 6564.614

	; read in the list of spectra for each object
	readcol,'spec_obs.lis',format='A,A,A,A,A,A,A,A',skipline=1,num,sdssj,ha_file,hb_file,lya_file,civ_file,mgii_file,hbuv_file,/silent

	; generate the sdss j names
	jname = 'J'+strtrim(strmid(sdssj,0,6),2)

	; set some parameters for each object
	; these are for smoothing the UV spectra
	readcol,'bins.dat',dum,smth_width_lya,smth_width_civ,smth_width_mgii,smth_width_ha,/silent

	; save the scale factors
	scales = 0.*num-9999.
	

	; loop over the objects	
	for i=1,n_elements(num)-1 do begin
		;if i eq 11 then goto,skip ; BBH080 has no usable data

		; read in the spectra and specfit models
		readcol,'./specfit/Lya/plot-'+lya_file[i],wave_lya_hires,flux_lya,err_lya,cont_lya,/silent
		if civ_file[i] ne 'nodata' then begin
			readcol,'./specfit/CIV/plot-'+civ_file[i],wave_civ_hires,flux_civ,err_civ,cont_civ,/silent
		endif else begin
			wave_civ_hires = wave_lya_hires 
			flux_civ = flux_lya
			cont_civ = cont_lya
		endelse
		if mgii_file[i] ne 'nodata' then begin
			readcol,'./specfit/MgII/plot-'+mgii_file[i],wave_mgii_hires,flux_mgii,err_mgii,cont_mgii,/silent
		endif else begin
			wave_mgii_hires = wave_lya_hires
			flux_mgii = flux_lya
			cont_mgii = cont_lya
		endelse
		if ha_file[i] ne 'nodata' then begin
			readcol,'./specfit/Ha/plot-'+ha_file[i],wave_ha_hires,flux_ha,err_ha,cont_ha,/silent
		endif else begin
			wave_ha_hires = wave_lya_hires
			flux_ha = flux_lya
			cont_ha = cont_lya
		endelse
		readcol,'./specfit/Hb/'+hb_file[i]+'_mod.txt',wave_hb,flux_hb,model_hb,pl_hb,fe2_hb,dum_hb,hbb_hb,hbn_hb,o3s_hb,o3l_hb,/silent

		; generate continuum-subtracted spectra
		lya_hires  = flux_lya-cont_lya
		civ_hires  = flux_civ-cont_civ
		mgii_hires = flux_mgii-cont_mgii
		ha_hires   = flux_ha-cont_ha 
		hb         = flux_hb-model_hb+hbb_hb
		hb         = flux_hb-pl_hb-fe2_hb-o3l_hb-o3s_hb

		; rebin some of the UV spectra
		boxsmooth_s,wave_lya_hires,lya_hires,smth_width_lya[i],wave_lya,lya
		boxsmooth_s,wave_civ_hires,civ_hires,smth_width_civ[i],wave_civ,civ
		boxsmooth_s,wave_mgii_hires,mgii_hires,smth_width_mgii[i],wave_mgii,mgii
		boxsmooth_s,wave_ha_hires,ha_hires,smth_width_ha[i],wave_ha,ha

		; convert the wavelength to velocity
		vel_lya = c*(((wave_lya/w0_lya)^2.)-1.)/(((wave_lya/w0_lya)^2.)+1.)
		vel_civ = c*(((wave_civ/w0_civ)^2.)-1.)/(((wave_civ/w0_civ)^2.)+1.)
		vel_mgii = c*(((wave_mgii/w0_mgii)^2.)-1.)/(((wave_mgii/w0_mgii)^2.)+1.)
		vel_hb = c*(((wave_hb/w0_hb)^2.)-1.)/(((wave_hb/w0_hb)^2.)+1.)
		vel_ha = c*(((wave_ha/w0_ha)^2.)-1.)/(((wave_ha/w0_ha)^2.)+1.)
		vel_nv = c*(((1240./w0_lya)^2.)-1.)/(((1240./w0_lya)^2.)+1.)	

		; reset some junk
		if civ_file[i] eq 'nodata' then civ = 0.*civ-1000.
		if mgii_file[i] eq 'nodata' then mgii = 0.*mgii-1000.
		if ha_file[i] eq 'nodata' then ha = 0.*ha-1000.

		; set the two lines to use
		thisline = 'lya'
		v1 = vel_hb
		v2 = vel_lya
		f1 = hb
		f2 = lya
	        epsfile = './plots/bbh'+num[i]+'_'+thisline+'scale_spec.eps'
		fitsfile = './'+thisline+'_scales.fits'	
		
		; skip if there's no data
		if thisline eq 'civ' and civ_file[i] eq 'nodata' then goto,skip
		if thisline eq 'mgii' and mgii_file[i] eq 'nodata' then goto,skip
		if thisline eq 'ha' and ha_file[i] eq 'nodata' then goto,skip

		; grab only the line velocities
		; this is necessary because the resampling
		; breaks if the v arrays are too different
		if i le 8 then begin
			vlim = 25000.
			vuse = where((v1 ge -vlim) and (v1 le vlim))
			v1   = v1[vuse]
			f1   = f1[vuse]
			vuse = where((v2 ge -vlim) and (v2 le vlim))
			v2   = v2[vuse]
			f2   = f2[vuse]	
		endif

		if num[i] eq 027 then begin
			vlim = 15000.
			vuse = where((v1 ge -vlim) and (v1 le vlim))
			v1   = v1[vuse]
			f1   = f1[vuse]
			vuse = where((v2 ge -vlim) and (v2 le vlim))
			v2   = v2[vuse]
			f2   = f2[vuse]	
		endif

		if num[i] eq 062 then begin
			vlim = 20000.
			vuse = where((v1 ge -vlim) and (v1 le vlim))
			v1   = v1[vuse]
			f1   = f1[vuse]
			vuse = where((v2 ge min(v2)) and (v2 le vlim))
			if thisline eq 'mgii' then vuse = where((v2 ge -vlim) and (v2 le vlim))
			if thisline eq 'ha' then begin
		       		vuse = where((v1 ge -vlim) and (v1 le max(v2)))
				v1   = v1[vuse]
				f1   = f1[vuse]	
				vuse = where((v2 ge -vlim) and (v2 le max(v2)))
			endif
			v2   = v2[vuse]
			f2   = f2[vuse]	
		endif

		if num[i] eq 072 then begin
			vlim = 20000.
			if thisline eq 'ha' then vlim = 10000
			vuse = where((v1 ge min(v2)) and (v1 le max(v1)))
			v1   = v1[vuse]
			f1   = f1[vuse]
			vuse = where((v2 ge min(v2)) and (v2 le max(v1)))
			v2   = v2[vuse]
			f2   = f2[vuse]	
		endif

		if num[i] eq 080 then begin
			vlim = 5000.
			vuse = where((v1 ge -vlim) and (v1 le vlim))
			v1   = v1[vuse]
			f1   = f1[vuse]
			vuse = where((v2 ge -vlim) and (v2 le vlim))
			v2   = v2[vuse]
			f2   = f2[vuse]	
		endif

		if num[i] eq 085 then begin
			vlim = 20000.
			vuse = where((v1 ge -vlim) and (v1 le vlim))
			v1   = v1[vuse]
			f1   = f1[vuse]
			vuse = where((v2 ge -vlim) and (v2 le vlim))
			v2   = v2[vuse]
			f2   = f2[vuse]	
		endif

		; resample to same velocity grid
		f1_old = f1
		f2_old = f2
		if v1[1]-v1[0] ge v2[1]-v2[0] then begin
			v = v1
		       f2_new = 0.*v
		       rebin_linear_flam,n_elements(v2),v2,f2,n_elements(v),v,f2_new 
		       f2     = f2_new

			use = where(f2 ne 0.)
		       f1  = f1[use]
		       f2  = f2[use]
		       v   = v[use] 
		endif else begin
			v = v2
		       f1_new = 0.*v
		       rebin_linear_flam,n_elements(v1),v1,f1,n_elements(v),v,f1_new
		       f1     = f1_new

		       use = where(f1 ne 0.)
		       f1  = f1[use]
		       f2  = f2[use]
		       v   = v[use]
		endelse	

		; choose wavelength windows for comparing the profiles
		; first, choose a window that goes down to 1% of the broad Hβ peak
		; then, excise a region around narrow Hβ with width of [O III] λ5007 at 10%
		; then run a 2/3σ outlier rejection algorithm to convergence to remove core emission
		; or intervening absorption from Lyα or NV.
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		; pick a wavelength window to avoid the narrow lines
		; do this by finding the [O III] λ5007 width
		; get the FWZM (actually, do 10% of the peak)
		; this is pretty narrow, but if it should be wider those 
		; points will be rejected later by the outlier algorithm

		; only excise the narrow Hβ region, not [O III], becuase 
		; otherwise we lose the whole red wing of the line profile
		fwzm_o3l = get_fwzm(wave_hb,o3l_hb,0.1)
		vlo      = -fwzm_o3l[0]
		vhi      = fwzm_o3l[0]

		; pick a window from the broad-line peak
		cwave_hbb = wave_hb[where(hbb_hb eq max(hbb_hb))]
		cdv_hbb   = c*(cwave_hbb-w0_hb)/w0_hb

		; do it by finding where the broad line is within 1% of max on each side
		; this will allow for asymmetric line profiles (which we have)
		; it's also very inclusive (and sometimes even overlaps C III*, but this 
		; will be handled by the outlier rejection algorithm)
		lo = c*(min(wave_hb[where(hbb_hb ge 0.01*max(hbb_hb))])-w0_hb)/w0_hb
		hi = c*(max(wave_hb[where(hbb_hb ge 0.01*max(hbb_hb))])-w0_hb)/w0_hb
		;if (thisline eq 'lya') and (num[i] eq 013) then hi = c*(max(wave_hb[where(hbb_hb ge 0.07*max(hbb_hb))])-w0_hb)/w0_hb

		vlo = cdv_hbb[0]+lo[0]
		vhi = cdv_hbb[0]+hi[0]

		hbb_lo = cdv_hbb[0]+lo 
		hbn_lo = -fwzm_o3l[0]/2. 
		hbn_hi = fwzm_o3l[0]/2.
		o3n_lo = (c*(4959.-w0_hb)/w0_hb)-fwzm_o3l[0] 
		o3n_hi = (c*(5007.-w0_hb)/w0_hb)+fwzm_o3l[0] 
		hbb_hi = cdv_hbb[0]+hi 
		;reg = where( ((v ge hbb_lo) and (v le hbn_lo)) or ((v ge hbn_hi) and (v le o3n_lo)) or ((v ge o3n_hi) and (v le hbb_hi)))
	        if num[i] ne 080 then begin	
			reg = where( ((v ge hbb_lo) and (v le hbn_lo)) or ((v ge hbn_hi) and (v le hbb_hi))) 
		endif else begin
			fwzm_o3l = get_fwzm(wave_hb,hbn_hb,0.1)
			hbn_lo      = -fwzm_o3l[0]
			hbn_hi      = fwzm_o3l[0]

			reg = where( ((v ge min(v)) and (v le hbn_lo)) or ((v ge hbn_hi) and (v le max(v))))

		endelse

		; do a preliminary fit to determine the scale factor
	        ; since these are continuum subtracted already, there should be no 
		; flux offset thus, only allow a scale factor between the line profiles
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		; set up the fit	
		x         = f1[reg]
	  	y         = f2[reg]
	  	xerr      = x/x 
	  	yerr      = y/y
	  	start     = [0.,50.]
		;if num[i] eq 027 then start = [0,1.]
    	  	pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},2)	
	  	pi(0).fixed=1 ; fit only a slope

		; perform fit
	  	result = mpfitfun('line',x,y,xerr,start,PARINFO=pi,BESTNORM=chi,PERROR=err,dof=dof)
		xfit   = dindgen(101)*(max(x)-min(x))/100+min(x)
	  	fit    = result[0] + result[1]*xfit

		; dο the first round of outlier rejection
		siglo    = 2. 
		sighi    = 3.

		if thisline eq 'lya' then begin
			if (num[i] eq 085) then begin
				siglo = 1
				sighi = 10
			endif
			if (num[i] eq 004) or (num[i] eq 027) then begin
				siglo = 1.
				sighi = 5.
			endif
			if (num[i] eq 013) then begin
				siglo = 0.2 
				sighi = 5. 
			endif
			if (num[i] eq 072) then begin
				siglo = 2
				sighi = 3. 
			endif
			if (num[i] eq 080) then begin
				siglo = 3.
				sighi = 2 
			endif
		endif
		if thisline eq 'civ' then begin
			if (num[i] eq 003) then begin
				siglo = 3.
				sighi = 3.
			endif
			if (num[i] eq 004) then begin
				siglo = 1.
				sighi = 3.
			endif
			if (num[i] eq 008) then begin
				siglo = 1.5
				sighi = 2.
			endif
			if (num[i] eq 013) then begin
				siglo = 0.2
				sighi = 5.
			endif

			if (num[i] eq 027) then begin
				siglo = 1.
				sighi = 4.
			endif
		endif
		if thisline eq 'mgii' then begin
			if (num[i] eq 004) then begin
				siglo = 1.
				sighi = 5.
			endif
			if (num[i] eq 008) then begin
				siglo = 2.
				sighi = 1.5
			endif
			if (num[i] eq 013) then begin
				siglo = 1. 
				sighi = 5. 
			endif
			if (num[i] eq 026) then begin
				siglo = 1.
				sighi = 3.
			endif
			if (num[i] eq 062) then begin
				siglo = 1.
				sighi = 1.9
			endif
			if (num[i] eq 085) then begin
				siglo = 0.8
				sighi = 3.
			endif
		endif


		resid  = y-(result[0]+result[1]*x)
		threshlo = siglo*stddev(resid)
		threshhi = sighi*stddev(resid)
		good = where((resid le threshhi) and (resid ge -threshlo))
		bad  = where((resid gt threshhi) or (resid lt -threshlo))
		
		; plot it
		ylo0 = -threshlo le min(resid) ? -threshlo : min(resid)
		yhi0 = threshhi ge max(resid) ? threshhi : max(resid)
		xtit = '!8f!X!dopt!n [10!u-17!n erg s!u-1!n cm!u-2!n '+cgsymbol('Angstrom')+'!u-1!n]' 
		ytit = '!8f!X!dUV-fit!n [10!u-17!n erg s!u-1!n cm!u-2!n '+cgsymbol('Angstrom')+'!u-1!n]' 
		ytit = 'Residual' 
		csize = 2 
		erase
		multiplot,[1,2],gap=0.01
		plot,x,y,/nodata,yra=[ylo0,yhi0],xra=[min(x),max(x)],xtit='',ytit='!8f!X!dUV!n [10!u-17!n erg s!u-1!n cm!u-2!n '+cgsymbol('Angstrom')+'!u-1!n]',charsize=csize
		oplot,xfit,fit,color='0000FF'x
		oplot,x,y,psym=1	
		al_legend,'BBH'+num[i],charsize=csize,/right,box=0
		multiplot
		plot,x,y-fit,/nodata,yra=[ylo0,yhi0],xra=[min(x),max(x)],xtit=xtit,ytit=ytit,charsize=csize
		oplot,[-1e5,1e5],[threshhi,threshhi],linestyle=1,color='00FF00'x
		oplot,[-1e5,1e5],[-threshlo,-threshlo],linestyle=1,color='00FF00'x
		oplot,[-1e5,1e5],[0,0],linestyle=2,color='FFFFFF'x
		oplot,x[good],resid[good],psym=1
		oplot,x[bad],resid[bad],psym=1,color='0000FF'x
	
		; update rejected pixels	
		rej   = bad
		keep  = good
		ind   = indgen(n_elements(y))
		niter = 10
		j     = 0
		if keyword_set(slow) then stop,'.cont to continue' else wait,1

		; now iterate the outlier rejection
		while ((j le niter-1) and (bad[0] ne -1)) do begin
			; redo the fit without the rejected pixels
	  		result = mpfitfun('line',x[keep],y[keep],xerr[keep],start,PARINFO=pi,BESTNORM=chi,PERROR=err,dof=dof)
	  		fit    = result[0] + result[1]*xfit

			; find the outliers
			resid  = y[keep]-(result[0]+result[1]*x[keep])
			threshlo = siglo*stddev(resid)
			threshhi = sighi*stddev(resid)
			bad  = where((resid gt threshhi) or (resid lt -threshlo))

			; plot it
			ylo = -threshlo le min(resid) ? -threshlo : min(resid)
			yhi = threshhi ge max(resid) ? threshhi : max(resid)
		 	erase	
			multiplot,[1,2],ygap=0.01
			plot,x,y,/nodata,yra=[ylo0,yhi0],xra=[min(x),max(x)],xtit='',ytit='!8f!X!dUV!n [10!u-17!n erg s!u-1!n cm!u-2!n '+cgsymbol('Angstrom')+'!u-1!n]',charsize=csize
			oplot,xfit,fit,color='0000FF'x
			oplot,x[rej],y[rej],psym=7,color='0000FF'x
			oplot,x[keep],y[keep],psym=1	
			al_legend,'BBH'+num[i],charsize=csize,/right,box=0
			multiplot
			plot,x[keep],resid,/nodata,yra=[ylo,yhi],xra=[min(x),max(x)],xtit=xtit,ytit=ytit,charsize=csize
			oplot,[-1e5,1e5],[threshhi,threshhi],linestyle=1,color='00FF00'x
			oplot,[-1e5,1e5],[-threshlo,-threshlo],linestyle=1,color='00FF00'x
			oplot,[-1e5,1e5],[0,0],linestyle=2,color='FFFFFF'x
			oplot,x[keep],resid,psym=1
			oplot,x[rej],y[rej]-(result[0]+result[1]*x[rej]),psym=7,color='0000FF'x

			; update rejected pixels
			if bad[0] ne -1 then begin
				oplot,x[keep[bad]],resid[bad],psym=1,color='0000FF'x
				rej  = [rej,(ind[keep])[bad]]
				remove,bad,keep 
			endif
			j++
			;if keyword_set(slow) then stop,'.cont to continue' else wait,1

		endwhile
		if keyword_set(slow) then stop,'.cont to continue' else wait,1

		erase
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
			clr = cgcolor('white')
		if keyword_set(hardcopy) then begin
			cgps_open,epsfile,/encapsulated,/nomatch,xsize=9,ysize=6.5;,xoffset=1.8,yoffset=-5
 	    		DEVICE, SET_FONT = 'Times-Roman'
			clr = cgcolor('black')
		endif
			
			; plot the scaled spectra
			fscale = 10
			ytit = '!C!8f!X!d'+cgsymbol('lambda')+'!n [10!u-16!n erg s!u-1!n cm!u-2!n '+cgsymbol('Angstrom')+'!u-1!n]' 
			xtit = 'Velocity [10!u3!n km s!u-1!n]' 	
			multiplot,[1,1]
			dx  = (hbb_hi-hbb_lo)
			xlo = abs(hbb_hi) ge abs(hbb_lo) ? -ceil(abs(hbb_hi))-0.1*dx : -ceil(abs(hbb_lo))-0.1*dx
			xhi = abs(hbb_hi) ge abs(hbb_lo) ? ceil(abs(hbb_hi))+0.1*dx : ceil(abs(hbb_lo))+0.1*dx

			fuse = max(f2[where(v ge -5000.)]) ge max(result[0]+result[1]*f1[where(v ge -5000.)]) ? f2 : result[0]+result[1]*f1 
			dy  = max(fuse[where(v ge -5000.)]) - min(fuse[where(v ge -5000.)]) 
			ylo = min(fuse[where(v ge -5000.)]) - 0.1*dy
			yhi = max(fuse[where(v ge -5000.)]) + 0.1*dy 

			if (((num[i] eq 023) or (num[i] eq 026)) and ((thisline eq 'civ') or (thisline eq 'mgii'))) then begin
				fuse = max(f2[where(v le 2000.)]) ge max(result[0]+result[1]*f1[where(v le 2000.)]) ? f2 : result[0]+result[1]*f1 
				dy  = max(fuse[where(v le 2000.)]) - min(fuse[where(v le 2000.)]) 
				ylo = min(fuse[where(v le 2000.)]) - 0.1*dy
				yhi = max(fuse[where(v le 2000.)]) + 0.1*dy 
			endif

			if ((num[i] eq 026) and (thisline eq 'mgii')) then begin
				fuse = max(f2[where(v le 0.)]) ge max(result[0]+result[1]*f1[where(v le 0.)]) ? f2 : result[0]+result[1]*f1 
				dy  = max(fuse[where(v le 0.)]) - min(fuse[where(v le 2000.)]) 
				ylo = min(fuse[where(v le 0.)]) - 0.1*dy
				yhi = max(fuse[where(v le 0.)]) + 0.1*dy 
			endif

			if ((num[i] eq 027) and (thisline eq 'civ')) then begin
				ylo = min(fuse[where(v ge -5000.)]) - 0.22*dy
			endif
			if ((num[i] eq 062) and (thisline eq 'ha')) then begin
				ylo = min(fuse[where(v le 5000.)]) - 0.1*dy
			endif
			if (num[i] eq 080) then begin
				xlo = -6d3
				xhi = 6d3	
			endif

			plot,v/1d3,f2/fscale,/nodata,psym=10,xra=[xlo/1d3,xhi/1d3],yra=[ylo/fscale,yhi/fscale],charsize=csize,xtit=xtit,ytit=ytit,xstyle=1,ystyle=1;,xtitoffset=1.8,ytitoffset=-2
	    		if num[i] ne 080 then begin
				polyfill,[xlo,hbb_lo,hbb_lo,xlo]/1d3,[ylo,ylo,yhi,yhi]/fscale,color=cgcolor('light gray') 		; horizontal zero line, left
	    			polyfill,[hbb_hi,xhi,xhi,hbb_hi]/1d3,[ylo,ylo,yhi,yhi]/fscale,color=cgcolor('light gray') 		; right
			endif
			polyfill,[hbn_lo,hbn_hi,hbn_hi,hbn_lo]/1d3,[ylo,ylo,yhi,yhi]/fscale,color=cgcolor('light gray') 	; narrow Hb
			oplot,[xlo/1d3,xhi/1d3],[0,0],color=clr,linestyle=2					; horizontal zero line
			oplot,[0,0],[ylo,yhi]/fscale,color=clr,linestyle=2					; vertical zero line 
			;oplot,(v[reg])[rej]/1d3,(f2[reg])[rej],psym=7,color=cgcolor('red'),symsize=0.5			; plot spectra and rejected points
			;oplot,(v[reg])[rej]/1d3,result[0]+(f1[reg])[rej]*result[1],psym=7,color=cgcolor('red'),symsize=0.5
			oplot,v/1d3,f2/fscale,psym=10,color=clr,thick=5							; plot the spectra
			oplot,v/1d3,(result[0]+f1*result[1])/fscale,psym=10,color=cgcolor('red'),thick=5

			if keyword_set(dolines) then begin								; mark some UV lines
				f_mark = (result[0]+result[1]*mean(f1[where((v gt vel_nv-150.) and (v lt vel_nv+150.))]))[0]								 
				oplot,[vel_nv,vel_nv]/1d3,[f_mark+0.1*dy,f_mark+0.2*dy]/fscale,color=cgcolor('black')
				xyouts,vel_nv/1d3-0.87,(f_mark+0.21*dy)/fscale,'N V',color=cgcolor('black'),charsize=1.5
				f_mark = yhi-0.06*dy								 
				oplot,[0.,0.]/1d3,[f_mark-0.11*dy,f_mark-0.01*dy]/fscale,color=cgcolor('black')
				xyouts,-0.87,f_mark,'Ly'+cgsymbol('alpha'),color=cgcolor('black'),charsize=1.5
				if num[i] eq 003 then begin								; mark CIII* for BBH003
					vel_c3 = v[where(f2 eq max(f2))] 
					oplot,[vel_c3,vel_c3]/1d3,[f_mark-0.11*dy,f_mark-0.01*dy]/fscale,color=cgcolor('black')
					xyouts,vel_c3/1d3-1.15,f_mark/fscale,'C III!u*!n',color=cgcolor('black'),charsize=1.5
				endif
			endif
			al_legend,jname[i],box=0,margin=-0.25,/right,charsize=2
			plot,v/1d3,f2/fscale,/nodata,psym=10,xra=[xlo/1d3,xhi/1d3],yra=[ylo/fscale,yhi/fscale],charsize=csize,xtit='',ytit='',ystyle=1,xstyle=1,/noerase
			if keyword_set(hardcopy) then begin
				cgps_close
				!p.font=-1
			endif

		; save the scale factor
		scales[i] = result[1]
		
		skip:
		if keyword_set(slow) then stop,'.cont to continue' else wait,1
		;stop,'.cont to continue'
	endfor

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	print, 'Writing properties to a fits file...'

	; add the results to the data structure
 	out_dat = CREATE_STRUCT('num',num,$
				'sdssj',sdssj,$
				'hb_file',hb_file,$
				'lya_file',lya_file,$
				'hbuv_file',hbuv_file,$	
				'lya_scales',scales)
	
       	; write everything to a fits file
	mwrfits,out_dat,fitsfile,/create

	; print progress message
	print,'...properties table saved.'



stop,'.cont to continue'
end
