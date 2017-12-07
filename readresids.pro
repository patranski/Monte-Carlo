;+
; NAME:
;       ReadResids
;
; PURPOSE:
;       Reads a TOA residual file output by TEMPO (the resid2.tmp file,
;       typically).  Returns a structure containing the file's data.
;
; CALLING SEQUENCE:
;       resids = ReadResids(files, count =, period =, /word32)
;
; INPUTS:
;       files   The name of the residual file to read.  Optionally, an array
;               of file names can be given, and the returned structure will
;               include the data from all of them.
;
; OPTIONAL OUTPUTS:
;       count   The number of TOA residuals returned.
;
;       period  The period of the pulsar.  Principally useful for converting
;               the TOA uncertainties stored in TimingUncertainty from
;               microseconds to cycles.
;
; KEYWORD PARAMETERS:
;       /word32 Indicates that the residual file was created on a 32-bit
;               machine.  The format of the TEMPO data files uses the machine
;               word length.  By default, a 64-bit word length is assumed.
;
; RETURNS:
;       An array of structures containing the data from the given file(s).
;       The notable structure entries are TOA (UTC MJD), PhaseResid (cycles),
;       TimeResid (s), OrbitalPhase (cycles), ObsFreq (unused for X-ray
;       timing), FitWeight (???), TimingUncertainty (us), and PrefitTimeResid
;       (cycles).  All entries are doubles, and hopefully the names are pretty
;       self-explanatory.
;
; MODIFICATION HISTORY:
;       Jake Hartman, written ... sometime in 2006, I guess
;-


function ReadResids, Files, COUNT = nRecords, WORD32 = Word32, PERIOD = Period
  
  ; Setup the structure as it is on the disk.  Note that it differs depending
  ; on whether IDL was run on a 32- or 64-bit machine.
  if keyword_set(Word32) then begin
    Dummy = { ResidRecord32, $
              Length1           : 0ul, $
              TOA               : 0d, $
              PhaseResid        : 0d, $
              TimeResid         : 0d, $
              OrbitalPhase      : 0d, $
              ObsFreq           : 0d, $
              FitWeight         : 0d, $
              TimingUncertainty : 0d, $
              PrefitTimeResid   : 0d, $
              Unused            : 0d, $
              Length2           : 0ul }
    nRecordLen = 80l
  endif else begin
    Dummy = { ResidRecord64, $
              Length1           : 0ull, $
              TOA               : 0d, $
              PhaseResid        : 0d, $
              TimeResid         : 0d, $
              OrbitalPhase      : 0d, $
              ObsFreq           : 0d, $
              FitWeight         : 0d, $
              TimingUncertainty : 0d, $
              PrefitTimeResid   : 0d, $
              Unused            : 0d, $
              Length2           : 0ull }
    nRecordLen = 88l
  endelse
  
  ; Setup for the read.
  nFiles = n_elements(Files)
  FileLens = (file_info(Files)).size / nRecordLen
  nRecords = total(FileLens)
  Data = replicate(Dummy, nRecords)
  i = 0
  
  ; Read in each file's data.  Make sure the formatting is as expected.
  for iFile = 0, nFiles-1 do begin
    NewData = replicate(Dummy, FileLens[iFile])
    openr, LUN, Files[iFile], /get_lun
    readu, LUN, NewData
    close, LUN
    free_lun, LUN
    
    iBadRecords = where( (NewData.Length1 NE 72) OR $
                         (NewData.Length2 NE 72) OR $
                         (NewData.Unused NE 0), nBadRecords )
    if nBadRecords GT 0 then $
      message, 'Found improperly formatted records in "' + Files[iFile] + $
               '".  Are you using the right machine word length?'
    
    if NewData[0].TOA LT 30000d then NewData.TOA += 39126
    
    Data[i : i + FileLens[iFile] - 1] = NewData
    i += FileLens[iFile]
  endfor
  
  ; If the period was requested, figure it out.
  if arg_present(Period) then begin
    iFirstGood = min(where(Data.PhaseResid NE 0))
    Period = Data[iFirstGood].TimeResid / $
             Data[iFirstGood].PhaseResid
  endif
  
  return, Data
end
