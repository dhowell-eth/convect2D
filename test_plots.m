i = load('T_initial.dat');
f = load('T_final.dat');
full = load('T_field.dat')';
i = i';
f = f';
phi = load('phi.dat');
phi = phi';

figure(1)
h1 = subplot(2,1,1)
contourf(phi)
h2 = subplot(2,1,2)
contourf(f)
%title(h1,'Final Stream Function (Ra=1e5; k=1; total_time=0.1; err=1e-3)')
%title(h2,'T Final (Ra=1e5; k=1; total_time=0.1; err=1e-3)')

%quiver(X,Y,vx',vy')

