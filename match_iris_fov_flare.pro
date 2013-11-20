;search for AIA flares and match them to IRIS observing sequences.
;Written: P.Higgins - 20131119
;Requires the GITHUB repository: https://github.com/pohuigin/gen_library
;Example: IDL> match_iris_fov_flare,/dowrite,outstruct=outstruct
;Notes: All outputs will be to './'

pro match_iris_fov_flare,dowrite=dowrite,doprint=doprint, outstruct=outstruct, res1tore=res1tore

outpath='./'

if keyword_set(doprint) then doprint=1 else doprint=0
if keyword_set(dowrite) then dowrite=1 else dowrite=0

;initialise output structure
strblank={foundmatch:0,obsid:'',obstart:'',obsend:'',obsx1:0d,obsx2:0d,obsy1:0d,obsy2:0d,flareclass:'',flarepeak:'',flarestart:'',flareend:'',flarex:0d,flarey:0d}

;Initialise output file
foutput=outpath+'iris_allobs_hek_flares.txt'

if dowrite then write_data_file,strblank,filecsv=foutput,/nodata,header=['#This data file was created by MATCH_IRIS_FOV_FLARE.PRO - P.Higgins','#The SSW LastEvents flares corresponding to each IRIS observation sequence.', '#The OBSID corresponds to that in the IRIS_ALLOBS_HEK data file.']

;Read in the last events list

if not keyword_set(res1tore) then begin

	montharr=anytim(timegrid('1-jan-2013',anytim(anytim(systim(/utc))+62.*24.*3600.,/vms),/day),/vms)
	nmonth=n_elements(montharr)
	flrsttim=0.
	flrentim=0.
	flrpktim=0.
	flrxcen=0.
	flrycen=0.
	flrcls=''
	
	for j=0,nmonth-2 do begin
;		sswfl=ssw_her_query(ssw_her_make_query(montharr[j]+' 00:00',montharr[j+1],/fl,search="FRM_NAME=SSW Latest Events",result_limit=500),/padref,/struct4)

		sswfl=ssw_her_query(ssw_her_make_query(montharr[j]+' 00:00',montharr[j+1],/fl,search="FRM_NAME=Flare Detective - Trigger Module",result_limit=500),/padref,/struct4)

		if data_type(sswfl) ne 8 then continue
help,sswfl.fl

;Filter flares by peak flux level and the wavelength
;plot_hist,sswfl.fl.optional.FL_PEAKFLUX,bin=100
;		wgood=where(sswfl.fl.optional.FL_PEAKFLUX gt 300 and sswfl.fl.required.OBS_CHANNELID eq '131')
;help,wgood
;if wgood[0] eq -1 then continue

		flrentim=[flrentim,anytim(sswfl.fl.required.EVENT_ENDTIME)]
	
		flrsttim=[flrsttim,anytim(sswfl.fl.required.EVENT_STARTTIME)]
		flrpktim=[flrpktim,anytim(sswfl.fl.required.EVENT_PEAKTIME)]
		
;		thispos=hel2arcmin(sswfl.fl.required.EVENT_COORD1,sswfl.fl.required.EVENT_COORD2,date=anytim((anytim(montharr[j])+anytim(montharr[j+1]))/2.))
		
;		flrxcen=[flrxcen,reform(thispos[1,*]/60.)]
;		flrycen=[flrycen,reform(thispos[0,*]/60.)]
		
		flrxcen=[flrxcen,float(sswfl.fl.required.EVENT_COORD1)]
		flrycen=[flrycen,float(sswfl.fl.required.EVENT_COORD2)]
		
;		flrcls=[flrcls,sswfl.fl.optional.FL_GOESCLS]
		flrcls=[flrcls,strtrim(fix(sswfl.fl.optional.FL_PEAKFLUX),2)]
		
		help,flrpktim
		
	endfor
	
	flrsttim=flrsttim[1:*]
	flrentim=flrentim[1:*]
	flrpktim=flrpktim[1:*]
	flrxcen=flrxcen[1:*]
	flrycen=flrycen[1:*]
	flrcls=flrcls[1:*]

	save,flrsttim,flrentim,flrpktim,flrxcen,flrycen,flrcls,sswfl,file=outpath+'match_iris_fov_flare.res1.sav'
	
endif else restore,outpath+'match_iris_fov_flare.res1.sav'


;top='/archive/sswdb/latest_events_genxcat'
;read_genxcat,'1-jan-2013','1-jan-2014',lescat,top=top

;flrsttim=anytim(lescat.FSTART)
;flrentim=anytim(lescat.FSTOP)
;flrpktim=anytim(lescat.FPEAK)

;flrxcen=lescat.XCEN
;flrycen=lescat.YCEN

;flrcls=lescat.class

;Read in the flare detective list

;qstring=ssw_hcr_make_query('2013-07-18 0:00',reltime(/now),instrument='IRIS')
;irisobs=ssw_hcr_query(qstring,/ssw,count=count)
;nlines = n_elements(irisobs.starttime)
;obsid = strmid(irisobs.eventid,57,10)
;x1 = irisobs.xcen - irisobs.xfov/2.
;x2 = irisobs.xcen + irisobs.xfov/2.
;y1 = irisobs.ycen - irisobs.yfov/2.
;y2 = irisobs.ycen + irisobs.yfov/2.

readcol,outpath+'iris_allobs_hek.txt',STARTTIME,STARTTIME2,STOPTIME,STOPTIME2,OBSID,X1,X2,Y1,Y2,comment='#',delim=' ',form='(A,A,A,A,A,A,A)' ;format='(A21,A21,A10,A7,A7,A7,A7)'

OBSID=strtrim(OBSID,2)

STARTTIME=STARTTIME+' '+STARTTIME2
sttim=anytim(STARTTIME)
STOPTIME=STOPTIME+' '+STOPTIME2
entim=anytim(STOPTIME)
xran=float([transpose(X1),transpose(X2)])
yran=float([transpose(Y1),transpose(Y2)])

stop

nobs=n_elements(OBSID)
outstruct=replicate(strblank,nobs)

for i=0,nobs-1 do begin

;Load the obs info into the structure	
	thisstr=strblank
	thisstr.obsid=OBSID[i]
	thisstr.obstart=anytim(sttim[i],/vms) ;STARTTIME[i]
	thisstr.obsend=anytim(entim[i],/vms) ;STOPTIME[i]
	thisstr.obsx1=xran[0,i] ;X1[i]
	thisstr.obsx2=xran[1,i] ;X2[i]
	thisstr.obsy1=yran[0,i] ;y1[i]
	thisstr.obsy2=yran[1,i] ;y2[i]
	outstruct[i]=thisstr

;Check flare time match	
	wflaretime=where(flrentim ge sttim[i] and flrsttim le entim[i])
	
	if wflaretime[0] eq -1 then begin
	
;if no match was found then enter the info of the closest matched flare
;		wbest=where(abs(flrpktim-(entim[i]+sttim[i])/2.) eq min(abs(flrpktim-(entim[i]+sttim[i])/2.)))
;		thisstr.flareclass=strjoin(flrcls[wbest],',')
;		thisstr.flarestart=strjoin(anytim(flrsttim[wbest],/vms),',')	
;		thisstr.flareend=strjoin(anytim(flrentim[wbest],/vms),',')	
;		thisstr.flarepeak=strjoin(anytim(flrpktim[wbest],/vms),',')	
;		thisstr.flarex=strjoin(strtrim(flrxcen[wbest],2),',')
;		thisstr.flarey=strjoin(strtrim(flrycen[wbest],2),',')
;		outstruct[i]=thisstr		
	
		if doprint then print,'No Flare time match for: OBSID=',OBSID[i]
		if dowrite then write_data_file,thisstr,filecsv=foutput,/append
		continue
	endif


;Check flare location match	

buffpos=50. ;buffer for position matching in arcseconds

	wflareloc=where(flrxcen[wflaretime] ge xran[0,i]-buffpos and flrxcen[wflaretime] le xran[1,i]+buffpos $
		and flrycen[wflaretime] ge yran[0,i]-buffpos and flrycen[wflaretime] le yran[1,i]+buffpos)

	if wflareloc[0] eq -1 then begin
		if doprint then print,'No Flare location match for: OBSID=',OBSID[i]
		if dowrite then write_data_file,thisstr,filecsv=foutput,/append
		continue
	endif
	
;Output what was found
	wmatch=wflaretime[wflareloc]

	thisstr.foundmatch=1

help,wmatch
	thisstr.flareclass=strjoin(flrcls[wmatch],',')
	thisstr.flarestart=strjoin(anytim(flrsttim[wmatch],/vms),',')	
	thisstr.flareend=strjoin(anytim(flrentim[wmatch],/vms),',')	
	thisstr.flarepeak=strjoin(anytim(flrpktim[wmatch],/vms),',')	
	thisstr.flarex=strjoin(strtrim(flrxcen[wmatch],2),',')
	thisstr.flarey=strjoin(strtrim(flrycen[wmatch],2),',')
	outstruct[i]=thisstr
	
	if doprint then print,'Flare match for OBSID=',OBSID[i],'Class=',strjoin(flrcls[wmatch],',')
	if dowrite then write_data_file,thisstr,filecsv=foutput,/append

endfor


stop

end