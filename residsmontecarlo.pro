; INPUT:
;   Resids  --  The array of timing residuals that we're analyzing, as output
;       by ReadResids().
;   Period  --  The period of the source, in seconds.
;   IntLen  --  The integration time for each TOA measurement, in seconds.
;   HalfMaxPeriod  --  The timescale at which the low-pass filter is halfway
;       between the 100% pass level (as it is for long timescales) and the
;       MinFilterLevel level (as it is for short timescales).
;   MinFilterLevel  -- The minimum fraction of signal to let through the low-
;       pass filter at short timescales.  Should be greater than zero to allow
;       through at least some short-timescale noise.
;   nTrials  --  Number of MC trials to run.
;   YRange  --  Range of the phase residuals to display.  Not used if
;       ShowSpectra = 1.
;   nDegree  --  The degree of the polynomial to fit.  Should be one more than
;       the degree used in the TEMPO fit.
;
; OPTIONAL INPUT:
;   ShowSpectra  --  Set to 1 to plot the spectra of the residuals.
;   TimeRange  --  Specify the time range to include in the analysis.


; Handle the optional arguments.
if n_elements(ShowSpectra) EQ 0 then ShowSpectra = 0

if n_elements(TimeRange) EQ 2 then begin
  i = where((TimeRange[0] LE Resids.TOA) AND (Resids.TOA LE TimeRange[1]))
  Resids = Resids[i]
endif

; Create an array of equally spaced phase residuals by rounding off the TOAs
; to the nearest IntLen.  If there are any points that get rounded off to the
; same time bin, increment the subsequent points until there's a gap.
nPts = n_elements(Resids.TOA)
iGoodPts = round((Resids.TOA - min(Resids.TOA)) / IntLen)

for i = 1, nPts-1 do begin
  if iGoodPts[i] EQ iGoodPts[i-1] then iGoodPts[i]++
endfor

; Copy the phase residuals into an FFT'able array, X.
nLen = ceil((max(iGoodPts) + 1) / 2) * 2
T = dindgen(nLen) * IntLen + min(Resids.TOA)
X = dblarr(nLen)
X[iGoodPts] = Resids.PhaseResid
XErrs = Resids.TimingUncertainty / (1d6 * Period)

; Get the magnitude and phases of the original phase residual array.
OrigStdDev = stdev(X[iGoodPts])
OrigSpec   = fft(X, +1)
Freqs      = dindgen(nLen / 2) / nLen / IntLen
Magnitudes = sqrt(double(OrigSpec)^2 + imaginary(OrigSpec)^2)
OrigPhases = atan(imaginary(OrigSpec), double(OrigSpec))

; Create the frequency filter and apply it to the magnitudes.
FreqSigma = 1 / (HalfMaxPeriod * sqrt(alog(2d)))
Filter = (1 - MinFilterLevel) * exp( -(Freqs / FreqSigma)^2 ) + MinFilterLevel
Unused = check_math(mask = 32)   ; Clear the underflow error from exp(-big).
Magnitudes[0:nLen/2-1] *= Filter
T0 = mean(Resids.TOA)

FitFreqMult = 1d6 / 86400d      ; Convert from cycles per day to uHz.
FitFDotMult = 1d14 / 86400d^2   ; Convert from cycles / day^2 to 10^-14 Hz/s.
FitFDDotMult = 1d20 / 86400d^3  ; Convert from cycles / day^3 to 10^-20 Hz/s^2.
; Setup the plots.  We just plot the original data and five trials for
; comparison.
!p.multi = [0, 1, 3]  &  !p.charsize = 2
!y.margin = [0,0]  &  !y.omargin = [4,2]
AllFitParams = dblarr(nTrials, nDegree)
XRange = [min(Resids.TOA) - 2, max(Resids.TOA) + 2]

; Run through the trials!  Note that trial 0 is actually the original data, so
; we can compare it with the MC trials.
for iTrial = 0, nTrials do begin
  
  ; For each trial, calculate new random phases.  Create a complex spectrum
  ; with them, then run an inverse FFT to get 
  NewPhases  = 2*!dpi * randomu(Seed, nLen/2-1, /double)
  NewSpec    = dcomplexarr(nLen)
  NewSpec[0] = Magnitudes[0]
  NewSpec[1:nLen/2-1] = Magnitudes[1:nLen/2-1] * exp(dcomplex(0, NewPhases))
  NewSpec[nLen/2+1:nLen-1] = reverse(conj(NewSpec[1:nLen/2-1]))
  NewX       = double((fft(NewSpec, -1))[iGoodPts])
  NewX      *= OrigStdDev / stdev(NewX)
  
  ; Throw out all the points not present in the original data.
  FitT = T[iGoodPts]
  FitX = (iTrial EQ 0) ? X[iGoodPts] : NewX
  
  ; If this is one of the first trials, plot it.  Plot either the power
  ; spectrum or the phase residuals, depending on whether ShowSpectra is set.
  if iTrial LT !p.multi[2] then begin
    if ShowSpectra then begin
      if iTrial EQ 0 then begin
        Spec = OrigSpec
      endif else begin
        X = dblarr(nLen)
        X[iGoodPts] = NewX
        Spec= fft(X, +1)
    endelse
;    Powers= sqrt(double(Spec)^2 + imaginary(Spec)^2)
    Powers=dblarr(n_elements(Freqs))
    for counter = 0, nLen/2-1 do Powers[counter] = sqrt(double(Spec[counter])^2 + imaginary(Spec[counter])^2)
;      Powers=Powers2
      if iTrial EQ 0 then MaxActualPower = max(Powers[1:*])
      plot, Freqs[1:*], Powers[1:*], /xlog, /ylog, xrange= [0.001, 100.0], yrange = [0.001, MaxActualPower], $
        xtickname = (iTrial EQ !p.multi[2]-1) ? '' : replicate(' ', 60), $
        ytickname = (iTrial EQ !p.multi[2]-1) ? '' : ' ' 
  endif else begin
      plot, FitT, FitX, psym = -1, $
        xrange = XRange, xstyle = 1, yrange = YRange, ystyle = 1, $
        xtickname = (iTrial EQ !p.multi[2]-1) ? '' : replicate(' ', 60), $
        xtickformat = (iTrial EQ !p.multi[2]-1) ? '(I0)' : '', $
        ytickname = (iTrial EQ !p.multi[2]-1) ? '' : ' '
    endelse
    xyouts, norm(!X, 0.9), norm(!Y, 0.8), string(format = '(I0)', iTrial)
  endif
  
  ; OK, here's the important bit: fit a degree nDegree polynomial to the fake
  ; residuals.  If we're past trial 0 (the original data), then save the
  ; results.  Print out the fit parameters and their errors.
  FitParams = poly_fit(FitT - T0, FitX, nDegree, measure_errors = XErrs, $
                       /double, sigma = FitParamErrs)
  Unused = check_math(mask = 32)  ; Another underflow from zero chisqr prob.
  if iTrial GT 0 then $
    AllFitParams[iTrial-1,*] = FitParams[1:*]
  
  case nDegree of
    1:  print, format = '(I4, T7, 2(D10.6))', $
          iTrial, FitParams[1] * FitFreqMult, FitParamErrs[1] * FitFreqMult
    2:  print, format = '(I4, T7, 2(D10.6), "  ", 2(D10.6))', $
          iTrial, FitParams[1] * FitFreqMult, FitParamErrs[1] * FitFreqMult, $
          FitParams[2] * FitFDotMult, FitParamErrs[1] * FitFDotMult
    3:  print, format = '(I4, T7, 2(D10.6), "  ", 2(D10.6), "  ", 2(D10.6))', $
          iTrial, FitParams[1] * FitFreqMult, FitParamErrs[1] * FitFreqMult, $
          FitParams[2] * FitFDotMult, FitParamErrs[1] * FitFDotMult, $
          FitParams[3] * FitFDDotMult, FitParamErrs[1] * FitFDDotMult
    else:  message, 'Unsupported value of nDegree.'
  endcase
endfor

; Calculate the standard deviation of the fit parameters.  These are our
; uncertainties.  We also calculate the means, but they should all be
; consistent with zero.  (If not, something is wrong.)
case nDegree of
  1:  print, format = '(" avg", T7, 2(D10.6))', $
        mean(AllFitParams) * FitFreqMult, stdev(AllFitParams) * FitFreqMult
  2:  print, format = '(" avg", T7, 2(D10.6), "  ", 2(D10.6))', $
        mean(AllFitParams[*,0]) * FitFreqMult, $
        stdev(AllFitParams[*,0]) * FitFreqMult, $
        mean(AllFitParams[*,1]) * FitFDotMult, $
        stdev(AllFitParams[*,1]) * FitFDotMult
  3:  print, format = '(" avg", T7, 2(D10.6), "  ", 2(D10.6),"  ", 2(D10.6))', $
        mean(AllFitParams[*,0]) * FitFreqMult, $
        stdev(AllFitParams[*,0]) * FitFreqMult, $
        mean(AllFitParams[*,1]) * FitFDotMult, $
        stdev(AllFitParams[*,1]) * FitFDotMult, $
        mean(AllFitParams[*,2]) * FitFDDotMult, $
        stdev(AllFitParams[*,2]) * FitFDDotMult 
  else:  message, 'Unsupported value of nDegree.'
endcase
end
