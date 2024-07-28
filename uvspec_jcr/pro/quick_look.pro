; J. Runnoe 
; 04/13/2022
; Take a quick look at the data after downloading them

pro go
	; object data
	z = 0.2144
	num = 019
	smth_width_lya = 10
	smth_width_civ = 10
	jname = 'SDSS J095036.75+512838.1'

	; set some constants
	c       = 2.99792458d5	; km/s
	w0_lya  = 1215.67
	w0_civ  = 1549.48
	w0_mgii = 2799.94
	w0_hb   = 4862.721
	w0_ha   = 6564.614

	; read in the optical spectra 
	ha = mrdfits('/Users/runnojc1/Dropbox/Shared/SMBBH/HETspec/PSU22-2-010/spectrum_20220325_0000009_exp01_red.fits',0)
	wave_ha = ha[*,0]/(1.+z)
	ha      = ha[*,1]*1d17 ; 1d-17 erg/s/cm2/A

	hb = mrdfits('/Users/runnojc1/Dropbox/Shared/SMBBH/HETspec/PSU22-2-010/spectrum_20220325_0000008_exp01_orange.fits',0)
	wave_hb = hb[*,0]/(1.+z)
	hb      = hb[*,1]*1d17 ; 1d-17 erg/s/cm2/A

	; read in the UV spectra
	readcol,'./2ca/bbh019_20220412_COSl.2ca',wave_lya_hires,lya_hires

	; smooth the UV spectra
	boxsmooth_s,wave_lya_hires,lya_hires,smth_width_lya,wave_lya,lya
	boxsmooth_s,wave_lya_hires,lya_hires,smth_width_civ,wave_civ,civ
	
	; convert the wavelength to velocity
	vel_lya = c*(((wave_lya/w0_lya)^2.)-1.)/(((wave_lya/w0_lya)^2.)+1.)
	vel_civ = c*(((wave_civ/w0_civ)^2.)-1.)/(((wave_civ/w0_civ)^2.)+1.)
	vel_hb = c*(((wave_hb/w0_hb)^2.)-1.)/(((wave_hb/w0_hb)^2.)+1.)
	vel_ha = c*(((wave_ha/w0_ha)^2.)-1.)/(((wave_ha/w0_ha)^2.)+1.)

	; set the two lines to use
	thisline = 'lya'
	v1 = vel_hb
	v2 = vel_lya
	f1 = hb
	f2 = lya
	
	plotfile = './plots/bbh019_quicklook_stack.eps'
	yscale          = 1. ; changes the figure y axis scale, make sure to change ytit accordingly

	ymax_ha         = 45.
	ymax_hb         = 21.
	ymax_lya        = 50.
	ymax_civ	= 50.
	ymax_mgii       = 0. 

	scale_lya       = 1.
	scale_civ       = 1.
	scale_mgii      = 1.
	scale_ha        = 1.

	ymin_ha         = -ymax_ha*0.05
	ymin_hb         = 11.
	ymin_lya        = -ymax_lya*0.05
	ymin_civ        = -ymax_civ*0.05
	ymin_mgii       = -ymax_mgii*0.05


	; make the velocity stack figure
	; plot the total raw light curve with a running median 
	cgps_open,plotfile,/encapsulated,xsize=7.5,ysize=11.
		; set the font to make Mike happy
		DEVICE, SET_FONT = 'Times-Roman'
	
		ang = cgsymbol("angstrom")
		l   = textoidl('\lambda')
		a   = textoidl('\alpha')
		b   = textoidl('\beta')

		ytit = '!8f!X!d'+l+'!n [10!u-17!n erg s!u-1!n cm!u-2!n '+ang+'!u-1!n]' ; input units are 1d-17, but I multiplited by yscale=1
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
		plot,vel_ha,ha,/nodata,xra=[-10000,10000],yra=[ymin_ha/yscale,ymax_ha/yscale],xstyle=1,ystyle=1,charsize=1.9,xtickname=replicate(' ',9)
		oplot,[0,0],[ymin_ha,ymax_ha]/yscale,linestyle=2
		oplot,[-10000,10000],[0,0],linestyle=2
		oplot,vel_ha,ha/yscale,psym=10,thick=4
		if keyword_set(domod) then oplot,vel_hb,hbb_hb*scale_ha/yscale,thick=4,color=cgcolor('red')
		al_legend,'H'+a,charsize=2,/left,box=0,margin=-0.25
		;al_legend,'BBH '+num,box=0,margin=-0.25
		xyouts,-10000,1.02*ymax_ha/yscale,jname;,box=0,margin=-0.25,/right,charsize=2

		multiplot
		!p.position=p2
		plot,vel_hb,hb,/nodata,xra=[-10000,10000],yra=[ymin_hb/yscale,ymax_hb/yscale],xstyle=1,ystyle=1,charsize=1.9,xtickname=replicate(' ',9)
		oplot,[0,0],[ymin_hb,ymax_hb]/yscale,linestyle=2
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
		plot,vel_hb,hb,/nodata,xra=[-10000,10000],yra=[-0.1,1],xstyle=1,ystyle=1,charsize=1.9,xtickname=replicate(' ',9)
		oplot,[0,0],[-0.1,1],linestyle=2
		oplot,[-10000,10000],[0,0],linestyle=2
		al_legend,'Mg II',charsize=2,/left,box=0,margin=-0.25	

		multiplot
		!p.position=p4
	        if total(where(civ ge 0.)) ne -1 then ymin_civ = min(civ[where((vel_civ ge -1d4) and (vel_civ le 1d4))]) le ymin_civ ? min(civ[where((vel_civ ge -1d4) and (vel_civ le 1d4))]) : ymin_civ	
		plot,vel_civ,civ,/nodata,xra=[-10000,10000],yra=[ymin_civ/yscale,ymax_civ/yscale],xstyle=1,ystyle=1,charsize=1.9,xtickname=replicate(' ',9)
		oplot,[0,0],[ymin_civ,ymax_civ]/yscale,linestyle=2
		oplot,[-10000,10000],[0,0],linestyle=2
		oplot,vel_civ,civ/yscale,psym=10,thick=4
		if keyword_set(domod) then oplot,vel_hb,hbb_hb*scale_civ/yscale,thick=4,color=cgcolor('red')
		al_legend,'C IV',charsize=2,/left,box=0,margin=-0.25

		multiplot
		!p.position=p5
		ymin_lya = min(lya[where((vel_lya ge -1d4) and (vel_lya le 1d4))]) le ymin_lya ? min(lya[where((vel_lya ge -1d4) and (vel_lya le 1d4))]) : ymin_lya	
		plot,vel_lya,lya,/nodata,xra=[-10000,10000],yra=[ymin_lya/yscale,ymax_lya/yscale],xstyle=1,ystyle=1,charsize=1.9;,xtickname=replicate(' ',9)
		oplot,[0,0],[ymin_lya,ymax_lya]/yscale,linestyle=2
		oplot,[-10000,10000],[0,0],linestyle=2
		oplot,vel_lya,lya/yscale,psym=10,thick=4
		if keyword_set(domod) then oplot,vel_hb,hbb_hb*scale_lya/yscale,thick=4,color=cgcolor('red')
		al_legend,'Ly'+a,charsize=2,/left,box=0,margin=-0.25

		;if keyword_set(domod) then begin
		;	oplot,vel_hb,hbb_hb*(scales.lya_scales)/yscale,thick=4,color=cgcolor('red')
		;	al_legend,'   , H'+b,charsize=2,/right,box=0,margin=-0.25,corners=c1,textcolor=cgcolor('red')
       		;	al_legend,'Ly'+a,linestyle=c2_line,charsize=2,corners=c2,box=0,margin=-0.25,pos=[0.755*c1[0],c1[1]],/norm 
		;endif; else begin
		;	al_legend,'Ly'+a,charsize=2,/right,box=0,margin=-0.25
		;endelse



		; clean up
		;;;;;;;;;;;;;;
 		multiplot
 		ERASE
		;;;;;;;;;;;;;;


	cgps_close
	!p.font=-1




stop,'.cont to continue'
end
