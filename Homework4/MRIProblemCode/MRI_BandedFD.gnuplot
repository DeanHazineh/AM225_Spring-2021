## Plot the saved source
dataSourceReal = 'Out/BandedSource_real'
dataSourceImag = 'Out/BandedSource_imag'

set output 'Out/BandedSource_real.png'
set terminal png size 1250,1000 font 'Helvetica,30"
set size square
set pm3d map
set title "Real Source"
splot dataSourceReal matrix binary

set output 'Out/BandedSource_imag.png'
set terminal png size 1250,1000 font 'Helvetica,30"
set title "Imaginary Source"
splot dataSourceImag matrix binary


## Plot the saved solution
dataSolReal = 'Out/BandedSol_real'
dataSolImag = 'Out/BandedSol_imag'

set output 'Out/BandedSol_real.png'
set title "FD Real Solution"
splot dataSolReal matrix binary

set output 'Out/BandedSol_imag.png'
set title "FD Imaginary Solution"
splot dataSolImag matrix binary


## When MMS is used as a check, plot the true solution and the error
#dataSolMMSReal = 'Out/mmsSolReal'
#dataSolMMSError = 'Out/BandedMMSError_real'
#
#set output 'Out/mmsSolReal.png'
#set title "Real MMS Solution"
#splot dataSolMMSReal matrix binary
#
#set output 'Out/BandedMMSError.png'
#set title "Difference MMS Analytic vs FD Approximated Solve"
#splot dataSolMMSError matrix binary

