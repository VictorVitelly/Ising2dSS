set title 'Peak location of C_v' font ',24'
set xlabel '1/L' font ',22' offset 0,-1
set ylabel 'T_{max}' font ',22' offset -6,0 rotate by 0
set lmargin 20
set bmargin 5
set xtics font ',16'
set ytics font ',16'
set yrange [2:3]
set grid x,y

set key left 
set key font ',16'
nn=10
#set logscale x
#set logscale y

array column1[nn]
array column2[nn]
array ta[nn]
array taerr[nn]
array tb[nn]
array tberr[nn]
array q[nn]
array chi2[nn]

f(x,a,b)=a*(1+b*x)

a=2.
b=0.3

do for [i=1:nn] {
set style line i pt 1
column1[i]=2*i
column2[i]=1+2*i
}

do for [i=1:5] {
q[i]=0.88+0.02*i
}
do for [i=6:10] {
q[i]=1.0+0.02*(i-5)
}

do for [i=1:nn] {
fit f(x,a,b) '../data/tccvtot.dat' u 1:column1[i]:column2[i] via a,b
ta[i]=a 
taerr[i]=a_err
tb[i]=b 
tberr[i]=b_err
chi2[i]=(FIT_STDFIT*FIT_STDFIT)
}

set print "results.dat"
do for [i=1:nn] {
print q[i], ta[i], taerr[i]
}
set print

plot for [i=1:nn] '../data/tccvtot.dat' u 1:column1[i]:column2[i] w errorbars notitle linestyle i lw 2, for [i=1:nn] f(x,ta[i],tb[i]) title sprintf('q=%.2f, T_c=%.4f±%.4f, a=%.4f±%.4f, χ^2/dof=%.2f',q[i],ta[i],taerr[i],tb[i],tberr[i],chi2[i]) linestyle i

pause -1
