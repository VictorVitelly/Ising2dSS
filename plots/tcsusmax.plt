set title 'Peak location of χ_m' font ',24'
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
array tc[nn]
array tcerr[nn]
array q[nn]
array chi2[nn]

f(x,a,b,c)=a*(x**(1/b))+c

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
fit f(x,a,b,c) '../data/tcsustot.dat' u 1:column1[i]:column2[i] via a,b,c
ta[i]=a 
taerr[i]=a_err
tb[i]=b 
tberr[i]=b_err
tc[i]=c 
tcerr[i]=a_err
chi2[i]=(FIT_STDFIT*FIT_STDFIT)
}

plot for [i=1:nn] '../data/tcsustot.dat' u 1:column1[i]:column2[i] w errorbars notitle linestyle i lw 2, for [i=1:nn] f(x,ta[i],tb[i],tc[i]) title sprintf('q=%.2f, T_c=%.2f±%.2f, ν=%.4f±%.4f, χ^2/dof=%.2f',q[i],tc[i],tcerr[i],tb[i],tberr[i],chi2[i]) linestyle i

pause -1
