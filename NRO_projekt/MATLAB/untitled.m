syms u(t) v(t)

de1=diff(u)==-3*u+2*v;
de2=diff(v)==4*u-v;

de=[de1;de2];

S=dsolve(de);
