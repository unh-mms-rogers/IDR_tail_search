;Created by Anthony J. Rogers (Univ. of New Hampshire)
;Change History:
;   2018-03-31 (ajr) -- added efil for low-pass filter of e-field data
;   2018-04-15 (ajr) -- Corrected velocity ofset reporting; add reporting for S1; added checks and reporting for S3
;   2018-04-17 (ajr) -- corrcted Hall E-field detection to include sign of E_z; Modified Hall B-field detection to detecting each quadrant and require ge 2 quads for pass
;   2018-06-26 (ajr) -- copied to _3ndel variant - removed delta requirement from S2  
;   2018-06-27 (ajr) -- increase needed Bz bounds for S1 detection to [-2,2] nT  
;   2018-07-10 (ajr) -- changed the ephemeris MEC datatype used to T89d as it is kept far more current.                    
;   2019-03-14 (ajr) -- cleaned comments and allowed stage 2 and stage 3 to be reversed in order
;   2019-03-22 (ajr) -- requre Hall E to search for Hall B; added requirement for correct polarity of By to be some ratio of total By in quadrant for Hall B

pro idr_tail_search_3ndel_clean, trange=trange, probe=probe, data_rate=data_rate, level=level, workdir=workdir, suffix=suffix, $
    coord=coord, efil=efil, swap=swap, s3bratio=s3bratio
  if undefined(probe) then probe = '1'
  if undefined(data_rate) then data_rate = 'fast'
  if undefined(level) then level = 'l2'
  if data_rate eq 'fast' or data_rate eq 'slow' or data_rate eq 'srvy' then bdatarate = 'srvy'
  if data_rate eq 'brst' then bdatarate = 'brst'
  if undefined(load) then load = 1
  if undefined(makeplot) then makeplot = 0
  if undefined(suffix) then suffix = ''
  if undefined(coord) then coord = 'gsm'
  if undefined(hpca) then hpca = 0
  if undefined(workdir) then workdir = '/home/buck/Work/MMS-IDR-2017'
  if undefined(efil) then efil = 1
  if undefined(swap) then swap = 0	; Switch to change order of checking for Hall fields and |E|
  if undefined(s3bratio) then s3bratio = 0.1	; ratio for of correct polarity By to total By in a quadrant to be Hall B
  Re = 6371.2                                 ; Earth radius in km for calculating
  
  CD, workdir   ; this should change to desired working directory
  mms_load_mec, trange = trange, probes = probe, level = level, data_rate = 'srvy', datatype='epht89d', /time_clip
  
  if hpca eq 0 then begin     ; Chooses FPI data
    mms_load_fpi, trange=trange, probes=[probe], level=level, data_rate=data_rate, datatype='dis-moms'
    mms_qcotrans, 'mms'+probe+'_dis_bulkv_gse_'+data_rate, 'mms'+probe+'_dis_bulkv_'+coord+'_'+data_rate, out_coord=coord
    get_data, 'mms'+probe+'_dis_bulkv_'+coord+'_'+data_rate, data = vdata
    get_data, 'mms'+probe+'_dis_tempperp_'+data_rate, data = tperpdata             ; perpendicular ion temp data product
  endif

  if size(vdata, /type) ne 8 or hpca eq 1 then begin  ; checks for successful structure in 'vdata', grabs ion data from HPCA if FPI data missing (Phase 1x)
    hpca = 1
    mms_load_hpca, trange=trange, probes=[probe], level=level, data_rate=bdatarate
    mms_qcotrans, 'mms'+probe+'_hpca_hplus_ion_bulk_velocity', 'mms'+probe+'_hpca_hplus_ion_bulk_velocity_'+coord, out_coord=coord
    options, 'mms'+probe+'_hpca_hplus_ion_bulk_velocity_'+coord, 'labels', ['Vx (H!U+!N) '+coord, 'Vy (H!U+!N) '+coord, 'Vz (H!U+!N) '+coord]
    get_data, 'mms'+probe+'_hpca_hplus_ion_bulk_velocity_'+coord, data=vdata                 ; HPCA ion velocity data - assumed majority protons
    get_data, 'mms'+probe+'_hpca_hplus_tperp', data=tperpdata                             ; perpendicular ion temp data product from HPCA protons
  endif

  ;  Load B-field and E-field data and get data for working variables
  mms_load_fgm, trange=trange, probes=probe, data_rate=bdatarate, instrument='fgm'
  get_data, 'mms'+probe+'_fgm_b_'+coord+'_'+bdatarate+'_'+level, data = bdata          ; B-field data -- Assumes desired coord exists in FGM data.  True for GSE, GSM, etc
  mms_load_edp, trange=trange, probes=[probe], level=level, data_rate=data_rate
  mms_qcotrans, 'mms'+probe+'_edp_dce_gse_'+data_rate+'_'+level, 'mms'+probe+'_edp_dce_'+coord+'_'+data_rate+'_'+level, out_coord=coord
  get_data, 'mms'+probe+'_edp_dce_'+coord+'_'+data_rate+'_'+level, data = edata         ; E-field data - using coordinate transformation

  ;  Apply low-pass filter to E-Field data for determination of Hall E fields and the adiabatic expansion coefficient
  if efil eq 1 then begin
    efil_coeff = digital_filter(0.0, 0.03125, 50, 16)  ;  0.5Hz LP filter; Assumes fast survey data; more subtle application needed to be robust for other data rates
    efil_data_x = convol(edata.y[*,0], efil_coeff)
    efil_data_y = convol(edata.y[*,1], efil_coeff)
    efil_data_z = convol(edata.y[*,2], efil_coeff)
    edata = {x:edata.x, y:[[efil_data_x], [efil_data_y], [efil_data_z]]}
  endif

  ; Need to use vdata and bdata which have the same time steps to check for correlated reversals, 
  ; so we'll interpolate both of them here as we would if we were already processing stage 2 soas
  ; not to do the work twice if we get a hit.  
  bdata_interp = {x:edata.x, y:[[interpol(bdata.y[*,0], bdata.x, edata.x)], [interpol(bdata.y[*,1], bdata.x, edata.x)], [interpol(bdata.y[*,2], bdata.x, edata.x)], [interpol(bdata.y[*,3], bdata.x, edata.x)]]}
  vdata_interp = {x:edata.x, y:[[interpol(vdata.y[*,0], vdata.x, edata.x)], [interpol(vdata.y[*,1], vdata.x, edata.x)], [interpol(vdata.y[*,2], vdata.x, edata.x)]]}
  
  start_time = time_double(trange[0])     ; converts the string time from the trange into unix time in seconds as is used in the data structures
  stop_time = start_time + 180
  
  s2flag = 0  ; Flag to detrmine if stage 2 parameters need to be calculated
    
  openw, lun, trange[0]+'_candidates'+suffix+'.txt', /get_lun ; Opens text file in working dir to record candidate intervals
  
  s1ctr = 0L  ; counter to track number of correlated Vx, Bz reversals in a time span
  esize = size(edata.x, /n_elements)-1  ; # of time steps in the trange.  Useful to have pre-calculated
  start_index = 0L  ; Must be a long interger to avoid overrun errors when the trange is more than a few hours
  
  ; ***
  ; Start search loop
  ; ***
  
  while stop_time le time_double(trange[1]) do begin  ; Control loop to run within the trange specified
    tctr = start_index  ; start looking roughly where the previous block was to avoid looking at the entire trange each loop
    while (edata.x[tctr] lt start_time) and (tctr lt esize) do tctr++  ; Loops to find indicies for the search period
    start_index = tctr
    while (edata.x[tctr] lt stop_time) and (tctr lt esize) do tctr++
    stop_index = tctr
    
    corpos = 0L ; # of time steps with both positive Bz and Vx
    corneg = 0L ; # of time steps with both negetive Bz and Vx
    
    for j=start_index, stop_index do begin  ; check at each time step in rage for S1 conditions
      if (vdata_interp.y[j,0] gt 100.) and (bdata_interp.y[j,2] gt 2.) then corpos++  ; Check for positive correlation, includes minimum values for both Bz and Vx, 
      										      ; assumed in nT and km/s
      if (vdata_interp.y[j,0] lt -100.) and (bdata_interp.y[j,2] lt -2.) then corneg++ ; check for negetive correlation
    endfor

    ; **************************************
    ; Conditional for Stage 2
    ; **************************************
    
    off_str = ''    ; string for later application of velocity offsets
    if (corpos ne 0) and (corneg ne 0) then off_str = off_str + 'VxBz,'  ; S1 pass for zero-velocity crossing
    
    if off_str ne '' then begin ; Requires a correctly correlated reversal to be true
    ;;;;;;;;;;;;;;;;;;;;;;;; Begin S2 & 3 checks ;;;;;;;;;;;;;;;;;;;;;;  
      
      s1ctr++ ; increment correlated Vx Bz reversal counter for later reporting
      
      s2ctr = 0L  ; counts number of time steps where S2 conditions are met

      ; ************************************
      ; Stage 2 Calculations and processing
      ; ************************************
      ;;;; --------NOTE: All S2 calculations done in spacecraft frame = NO adjustment for x-line motion

      if swap eq 0 then begin
        if s2flag eq 0 then begin     ; calculates stage 2 parameters if they have not been calculated for this trange yet

          mag_e = dblarr(esize+1)
          for m=0, esize do begin
            mag_e[m] = norm(transpose(edata.y[m,*]))
          endfor
          s2flag = 1  ; Sets flag so S2 parameters aren't re-calculated next time
        endif   ; Done calculating stage 2 parameters
      
        for m=start_index, stop_index do begin  ; Checks for existance of sufficient value of |E|
          if mag_e[m] ge 10. then begin
            s2ctr++
            ; print, mag_e[m]  ; Disgnostic line
          endif
        endfor
        
        if s2ctr gt 0 then s2tag = ' - |E|' else s2tag = ''  ; creates label for reporting if S2 satisfied for time period
      endif
      ; ***********************************
      ; End of Stage 2
      ; ***********************************
      
      ; ***********************************
      ; Beginning of Stage 3
      ; ***********************************
      
      s3tag = ''
      if s2ctr gt 0 or swap eq 1 then begin  ; Begin checks for stage 3
        s3e = 0   ; initialize counters for Hall electric
        s3b = 0   ; and Magnetic fields
        s3q1flag = 0    ; initialize flags for detecting Hall magnetic fields 
        s3q2flag = 0    ; in each of the four canonical 2D quadrants
        s3q3flag = 0
        s3q4flag = 0
	q1t = 0
	q2t = 0
	q3t = 0
	q4t = 0
        
      	bfil_coeff = digital_filter(0.0, 0.0625, 50, 16)  ;  0.5Hz LP filter; Assumes fast survey data; more subtle application needed to be robust for other data rates
      	bfil_data_x = convol(bdata_interp.y[*,0], bfil_coeff)
      	bfil_data_y = convol(bdata_interp.y[*,1], bfil_coeff)
      	bfil_data_z = convol(bdata_interp.y[*,2], bfil_coeff)
      	bdata_fil = {x:edata.x, y:[[bfil_data_x], [bfil_data_y], [bfil_data_z]]}
      
        for m=start_index, stop_index do begin  ; step through for hall E field
          if (edata.y[m,2] * signum(bdata_fil.y[m,0]) le -2.) then begin ; check for Hall E-field while ions demagnitized 
	    s3e++
	    ;;;;; Checking for all four Hall-B quadrants. Check requires that there exist a measurement that could satisfy Hall B_y requirements in at least two 
	    ;		quadrants.  DOES NOT check for sustained or average B_y in correct polarity across time window
            if abs(bdata_fil.y[m,1]) ge 1. then begin
              if (signum(vdata_interp.y[m,0]) gt 0.) and (signum(bdata_fil.y[m,0]) gt 0.) then begin; Q1
		      q1t++							; Total number of data points in the quadrant
		      if (signum(bdata_fil.y[m,1]) gt 0.) then s3q1flag++	; Number of data points with correct By polarity for Hall B
	      endif
              if (signum(vdata_interp.y[m,0]) lt 0.) and (signum(bdata_fil.y[m,0]) gt 0.) then begin; Q4
		      q4t++
		      if (signum(bdata_fil.y[m,1]) lt 0.) then s3q4flag++
	      endif
              if (signum(vdata_interp.y[m,0]) lt 0.) and (signum(bdata_fil.y[m,0]) lt 0.) then begin; Q2 
		      q2t++
		      if (signum(bdata_fil.y[m,1]) gt 0.) then s3q2flag++
	      endif
              if (signum(vdata_interp.y[m,0]) gt 0.) and (signum(bdata_fil.y[m,0]) lt 0.) then begin; Q3
		      q3t++
		      if (signum(bdata_fil.y[m,1]) lt 0.) then s3q3flag++
	      endif

            endif
	  endif
        endfor
        
	      if q1t ne 0 then q1ratio = float(s3q1flag)/q1t else q1ratio = 0 
	      if q2t ne 0 then q2ratio = float(s3q2flag)/q2t else q2ratio = 0
	      if q3t ne 0 then q3ratio = float(s3q3flag)/q3t else q3ratio = 0
	      if q4t ne 0 then q4ratio = float(s3q4flag)/q4t else q4ratio = 0

	      if q1ratio gt s3bratio then s3b++
	      if q2ratio gt s3bratio then s3b++
	      if q3ratio gt s3bratio then s3b++
	      if q4ratio gt s3bratio then s3b++

        if (s3e ge 2) and (s3b ge 2) then s3tag = ' - Hall('+string(s3b)+'---'+string(q1ratio)+','+string(q2ratio)+','+string(q3ratio)+','+string(q4ratio)+')' else s3tag = ''   ; creates label for reporting if S3 satisfied for time period
      endif

      if swap eq 1 then begin
      if s2flag eq 0 then begin     ; calculates stage 2 parameters if they have not been calculated for this trange yet

          mag_e = dblarr(esize+1)
          for m=0, esize do begin
            mag_e[m] = norm(transpose(edata.y[m,*]))
          endfor
          s2flag = 1  ; Sets flag so S2 parameters aren't re-calculated next time
        endif   ; Done calculating stage 2 parameters

        s2ctr = 0L  ; counts number of time steps where S2 conditions are met
      
        for m=start_index, stop_index do begin  ; Checks for existance of sufficient value of |E|
          if mag_e[m] ge 10. then begin
            s2ctr++
            ; print, mag_e[m]  ; Disgnostic line
          endif
        endfor
        
        if s2ctr gt 0 then s2tag = ' - |E|' else s2tag = ''  ; creates label for reporting if S2 satisfied for time period
      endif
    
      ; ***********************************
      ; End of Stage 3
      ; ***********************************
      
      ; Reporting of results for each segment.  time segment must pass at least S1 to be reported
      print, "MMS"+probe+", "+time_string(start_time) + ", " + time_string(stop_time) + ', ' + off_str + s2tag + s3tag + "; "       ; prints to console each candidate segment
      printf, lun, "MMS"+probe+", "+time_string(start_time) + ", " + time_string(stop_time) + ', ' + off_str + s2tag + s3tag + "; " ; prints to file each candidate segment
     
    endif   
    ;;;;;;;;;;;;;;;; End S2 & 3 checks ;;;;;;;;;;;;;;;;;;;;;;;;

    
    start_time = start_time + 60    ;  increments to next search period (1 min)
    stop_time = stop_time + 60
    
  endwhile
  
  ; ***********************************
  ; Done with search loop for given trange
  ; ***********************************
  
  del_data, '*'         ; tplot clean-up
  
  printf, lun, "Total correlated reversals:  "+string(s1ctr)  ; reports # of reversals on last line in file for future statistics
  free_lun, lun         ; Close candidate file
  
  print, ''
  print, "***************************************"
  print, "Done with "+trange[0]+' - '+trange[1]   ; Report done.
  print, "***************************************"
  print, ''
    
end
  
  
