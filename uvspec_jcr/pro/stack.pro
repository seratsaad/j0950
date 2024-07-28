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
pro stack,HARDCOPY=HARDCOPY,DOMOD=DOMOD
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

	; set filenames
	if keyword_set(domod) then begin
		plotfile = './plots/bbh'+num+'_vel_wmod.eps'
	endif else begin
		plotfile = './plots/bbh'+num+'_vel.eps'
	endelse

	; get the scale factors
	scales = mrdfits('./lya_scales.fits',1,hdr)

	; set some parameters for each object
	readcol,'bins.dat',dum,smth_width_lya,smth_width_civ,smth_width_mgii,smth_width_ha,/silent
	
	yscale          = 10. ; changes the figure y axis scale, make sure to change ytit accordingly

	;ymin_ha         = -ymax_ha*0.0 
	;ymin_hb         = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
	;ymin_lya        = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
	;ymin_civ        = [+50,+29,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]-49.
	;ymin_mgii       = [0.,-6.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]-2.
	ymax_ha         = [160.,100.,250,250,200.,100.,175.,300.,175.,120.,250.,200.,99.]
	ymax_hb         = [49.,19.,110,23,59.,39.,79.,119.,199.,24.,89.,3.6,36.]
	ymax_lya        = [170.,799.,3999,2099,1510.,2099.,399.,1650.,999.,299.,2200.,359.,499.]
	ymax_civ	= [99.,299.,449,999,609.,699.,199.,399.,399.,99.,99.,900.,99.]
	ymax_mgii       = [36.,119.,199,199,119.,99.,99.,99.,119.,79.,99.,11.,79.]

	scale_lya       = [2.95,6.,24.,28.,55.,26.,5.,10.,4.0,9.,25.,30.,22.]
	scale_civ       = [-10.,6.,5.,15.,40.,19.,2.5,4.,2.,-10.,-10.,-10.,-10.]
	scale_mgii      = [0.6,2.5,2.,4.,7.,3.,1.7,1.1,1.,3.,-10.,2.4,2.6]
	scale_ha        = [3.,2.,0.6,4.,4.,3.,1.5,2.5,0.7,3.7,1.9,-10.,-10.]

	ymin_ha         = -ymax_ha*0.05
	ymin_hb         = -ymax_hb*0.05 
	ymin_lya        = -ymax_lya*0.05
	ymin_civ        = -ymax_civ*0.05
	ymin_mgii       = -ymax_mgii*0.05

	; 085 mgii could also be 12. to scale on blue side

	; loop over the objects	
	for i=0,n_elements(num)-1 do begin
		;if i eq 11 then goto,skip
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
		hb         = flux_hb-pl_hb-fe2_hb

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

		; reset some junk
		if civ_file[i] eq 'nodata' then civ = 0.*civ-1000.
		if mgii_file[i] eq 'nodata' then mgii = 0.*mgii-1000.
		if ha_file[i] eq 'nodata' then ha = 0.*ha-1000.


;		plot,vel_civ,civ,xra=[-10000,10000],yra=[ymin_civ[i],ymax_civ[i]],xstyle=1,ystyle=1,psym=10
;		oplot,vel_hb,hbb_hb*10.,thick=1,color=cgcolor('red')                                       
;		stop,'.cont to continue'

;		plot,vel_mgii,mgii,xra=[-10000,10000],yra=[ymin_mgii[i],ymax_mgii[i]],xstyle=1,ystyle=1,psym=10
;		oplot,vel_hb,hbb_hb*10.,thick=1,color=cgcolor('red')                                       
;		stop,'.cont to continue'

;		plot,vel_ha,ha,xra=[-10000,10000],yra=[ymin_ha[i],ymax_ha[i]],xstyle=1,ystyle=1,psym=10
;		oplot,vel_hb,hbb_hb*10.,thick=1,color=cgcolor('red')                                       
;		stop,'.cont to continue'


		;if i eq 11 then stop,'.cont to continue'


		; make the velocity stack figure
		; plot the total raw light curve with a running median 
		cgps_open,plotfile[i],/encapsulated,xsize=7.5,ysize=11.
			; set the font to make Mike happy
			DEVICE, SET_FONT = 'Times-Roman'
		
	  		ang = cgsymbol("angstrom")
			l   = textoidl('\lambda')
			a   = textoidl('\alpha')
			b   = textoidl('\beta')

			ytit = '!8f!X!d'+l+'!n [10!u-16!n erg s!u-1!n cm!u-2!n '+ang+'!u-1!n]' ; input units are 1d-17, but I multiplited by yscale=10
			xtit = 'Velocity [km s!u-1!n]'

			multiplot,[1,5],mXtitle=xtit,mYtitle=ytit,mxtitsize=2,mytitsize=2,mxtitoffset=1.8,mytitoffset=-2.
	 		p1=!p.position 	; at each plot I save the position keyword for later use
			p1[0] = 1.5*p1[0]
			p1[2] = 0.9*p1[2] 
			multiplot
	 		p2=!p.position
			p2[0] = 1.5*p2[0]
			p2[2] = 0.9*p2[2] 
	 		multiplot
	 		p3=!p.position
			p3[0] = 1.5*p3[0]
			p3[2] = 0.9*p3[2] 
	 		multiplot
	 		p4=!p.position
			p4[0] = 1.5*p4[0]
			p4[2] = 0.9*p4[2] 
	 		multiplot
	 		p5=!p.position
			p5[0] = 1.5*p5[0]
			p5[2] = 0.9*p5[2] 
	 		multiplot
 			!p.position=p1
			!p.position[1]=p3[3]; here I edit the position keyword so it spans 2 plot boxes

			; make the plot box
			multiplot
			!p.position=p1
			plot,vel_ha,ha,/nodata,xra=[-10000,10000],yra=[ymin_ha[i]/yscale,ymax_ha[i]/yscale],xstyle=1,ystyle=1,charsize=1.9,xtickname=replicate(' ',9)
			oplot,[0,0],[ymin_ha[i],ymax_ha[i]]/yscale,linestyle=2
			oplot,[-10000,10000],[0,0],linestyle=2
			oplot,vel_ha,ha/yscale,psym=10,thick=4
			if keyword_set(domod) then oplot,vel_hb,hbb_hb*scale_ha[i]/yscale,thick=4,color=cgcolor('red')
			al_legend,'H'+a,charsize=2,/left,box=0,margin=-0.25
			;al_legend,'BBH '+num[i],box=0,margin=-0.25
			al_legend,jname[i],box=0,margin=-0.25,/right,charsize=2

			multiplot
			!p.position=p2
			plot,vel_hb,hb,/nodata,xra=[-10000,10000],yra=[ymin_hb[i]/yscale,ymax_hb[i]/yscale],xstyle=1,ystyle=1,charsize=1.9,xtickname=replicate(' ',9)
			oplot,[0,0],[ymin_hb[i],ymax_hb[i]]/yscale,linestyle=2
			oplot,[-10000,10000],[0,0],linestyle=2
			if keyword_set(domod) then oplot,vel_hb,hbb_hb/yscale,thick=4,color=cgcolor('red')
			oplot,vel_hb,hb/yscale,psym=10,thick=4
			if keyword_set(domod) then begin
				al_legend,'H'+b,charsize=2,/left,box=0,margin=-0.25,textcolor=cgcolor('red')
			endif else begin 
				al_legend,'H'+b,charsize=2,/left,box=0,margin=-0.25
			endelse

			multiplot
			!p.position=p3
		        if total(where(mgii ge 0.)) ne -1 then ymin_mgii[i] = min(mgii[where((vel_mgii ge -1d4) and (vel_mgii le 1d4))]) le ymin_mgii[i] ? min(mgii[where((vel_mgii ge -1d4) and (vel_mgii le 1d4))]) : ymin_mgii[i]	
			plot,vel_mgii,mgii,/nodata,xra=[-10000,10000],yra=[ymin_mgii[i]/yscale,ymax_mgii[i]/yscale],xstyle=1,ystyle=1,charsize=1.9,xtickname=replicate(' ',9)
			oplot,[0,0],[ymin_mgii[i],ymax_mgii[i]]/yscale,linestyle=2
			oplot,[-10000,10000],[0,0],linestyle=2
			oplot,vel_mgii,mgii/yscale,psym=10,thick=4
			if keyword_set(domod) then oplot,vel_hb,hbb_hb*scale_mgii[i]/yscale,thick=4,color=cgcolor('red')
			al_legend,'Mg II',charsize=2,/left,box=0,margin=-0.25	

			multiplot
			!p.position=p4
		        if total(where(civ ge 0.)) ne -1 then ymin_civ[i] = min(civ[where((vel_civ ge -1d4) and (vel_civ le 1d4))]) le ymin_civ[i] ? min(civ[where((vel_civ ge -1d4) and (vel_civ le 1d4))]) : ymin_civ[i]	
			plot,vel_civ,civ,/nodata,xra=[-10000,10000],yra=[ymin_civ[i]/yscale,ymax_civ[i]/yscale],xstyle=1,ystyle=1,charsize=1.9,xtickname=replicate(' ',9)
			oplot,[0,0],[ymin_civ[i],ymax_civ[i]]/yscale,linestyle=2
			oplot,[-10000,10000],[0,0],linestyle=2
			oplot,vel_civ,civ/yscale,psym=10,thick=4
			if keyword_set(domod) then oplot,vel_hb,hbb_hb*scale_civ[i]/yscale,thick=4,color=cgcolor('red')
			al_legend,'C IV',charsize=2,/left,box=0,margin=-0.25

			multiplot
			!p.position=p5
			ymin_lya[i] = min(lya[where((vel_lya ge -1d4) and (vel_lya le 1d4))]) le ymin_lya[i] ? min(lya[where((vel_lya ge -1d4) and (vel_lya le 1d4))]) : ymin_lya[i]	
			plot,vel_lya,lya,/nodata,xra=[-10000,10000],yra=[ymin_lya[i]/yscale,ymax_lya[i]/yscale],xstyle=1,ystyle=1,charsize=1.9;,xtickname=replicate(' ',9)
			oplot,[0,0],[ymin_lya[i],ymax_lya[i]]/yscale,linestyle=2
			oplot,[-10000,10000],[0,0],linestyle=2
			oplot,vel_lya,lya/yscale,psym=10,thick=4
			if keyword_set(domod) then oplot,vel_hb,hbb_hb*scale_lya[i]/yscale,thick=4,color=cgcolor('red')
			al_legend,'Ly'+a,charsize=2,/left,box=0,margin=-0.25

			if keyword_set(domod) then begin
				oplot,vel_hb,hbb_hb*(scales.lya_scales)[i]/yscale,thick=4,color=cgcolor('red')
				al_legend,'   , H'+b,charsize=2,/right,box=0,margin=-0.25,corners=c1,textcolor=cgcolor('red')
       	  			al_legend,'Ly'+a,linestyle=c2_line,charsize=2,corners=c2,box=0,margin=-0.25,pos=[0.755*c1[0],c1[1]],/norm 
			endif; else begin
			;	al_legend,'Ly'+a,charsize=2,/right,box=0,margin=-0.25
			;endelse



			; clean up
			;;;;;;;;;;;;;;
 			multiplot
 			ERASE
			;;;;;;;;;;;;;;


	        cgps_close
		!p.font=-1
	
	










		skip:


	


		;stop,'.cont to continue'
	endfor






















;		if keyword_set(hardcopy) then begin
;			cgps_open,plotfile_coss,/encapsulated,xsize=7.5,ysize=6
;			; set the font to make Mike happy
;			DEVICE, SET_FONT = 'Times-Roman'
;		endif
;			ang  = cgSymbol('Angstrom')
;			lam  = cgSymbol('lambda')
;			plot,wave_short,flux_short,psym=10,ystyle=1,xstyle=1,ytit='f!d'+lam+'!n [10!u-17!n erg s!u-1!n cm!u-2!n '+ang+'!u-1!n]',xtit='Rest Wavelength ['+ang+']',charsize=2,xra=[800,xmax],yra=[ymin,ymax]
;			oploterror,wave_short,flux_short,fluxerr_short,psym=3,color=cgcolor('gray'),errcolor=cgcolor('gray'),/nohat
;			oplot,wave_short,flux_short,psym=10
;			plot,wave_short,flux_short,psym=10,ystyle=1,xstyle=1,ytit='f!d'+lam+'!n [10!u-17!n erg s!u-1!n cm!u-2!n '+ang+'!u-1!n]',xtit='Rest Wavelength ['+ang+']',charsize=2,xra=[800,xmax],yra=[ymin,ymax],/noerase
;		if keyword_set(hardcopy) then begin
;			cgps_close
;			!p.font=-1
;		endif	
	
	
		
	

stop,'.cont to continue'
end
