function Jz = Jz_BHW(t,T)
a0 = 0.35322222;
a1 = -0.488;
a2 = 0.145;
a3 = -0.01022222;

Jz = a0 + a1*cos(2*pi()*(t/T)) + a2*cos(4*pi()*(t/T)) + a3*cos(6*pi()*(t/T)); 








