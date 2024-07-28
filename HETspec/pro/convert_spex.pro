; takes an input directory
; finds fits files in it
; and prints out a ASCII file
; with the same filename
;
; June 30, 2021
; J. Runnoe
pro go, dir, HARDCOPY=HARDCOPY

	spawn,'ls ./'+dir+'/*fits',files
	filename = strarr(n_elements(files))
	for i=0,n_elements(files)-1 do begin
		filename[i]=(StrSplit(files[i],'.fits', /Regex, /Extract))[0]

		dat = mrdfits(files[i],0,hdr)

		file_hdr0 = 'col1: wavelength (air)'
		file_hdr1 = 'col2: extracted object spectrum (f_lambda: ergs/s/cm^2/A)'
		file_hdr2 = 'col3: extracted sky spectrum from same aperture and weighting as object (s_lambda: ergs/s/cm^2/A)'
		file_hdr3 = 'col4: error for extracted object spectrum (e_f_lambda: ergs/s/cm^2/A)'
		file_hdr4 = 'col5: error for extracted sky spectrum (e_s_lambda: ergs/s/cm^2/A)'
		file_hdr5 = 'col6: response function (ergs / e-)'

		; write out all the files
		if keyword_set(hardcopy) then begin
			openw,1,strtrim(filename[i],2)+'.txt',width=1600
				printf,1,file_hdr0
				printf,1,file_hdr1
				printf,1,file_hdr2
				printf,1,file_hdr3
				printf,1,file_hdr4
				printf,1,file_hdr5
				printf,1,'================================================================================================================================='
				for j=0,n_elements(dat[*,0])-1 do printf,1,strtrim(string(dat[j,0]),2),'        ',dat[j,1],'        ',dat[j,2],'        ',dat[j,3],'        ',dat[j,4],'        ',dat[j,5]
			close,1
			openw,1,strtrim(filename[i],2)+'.hdr',width=1600
				for j=0,n_elements(hdr)-1 do printf,1,hdr[j]
			close,1
		endif
	endfor



stop,'.cont to continue'
end
