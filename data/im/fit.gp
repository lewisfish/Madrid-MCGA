set terminal pngcairo size 1920,1080 font ",20"
set output 'testy.png'
#set fit quiet
#set fit logfile '/dev/null'
f(x) = a * exp(-(x-mu)**2 / (2.*sigma)**2)
a = 1500
mu = 100
sigma = 25
fit f(x) 'test.dat' via a, mu, sigma

p  f(x) lw 5,'test.dat' lw 5 w p
