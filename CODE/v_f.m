function [v_f] = v_f(u,v,nx,ny,dx,dy)
% Calculate the Nonlinear step
% for Dvdt=Fv
v_f=v;
cx=1/2/dx;
cy=1/2/dy;
for j=2:ny-1 %rows
    for k=2:nx-1 %cols
        a=u(j,k)*(v(j,k+1)-v(j,k-1))*cx;
        b=v(j,k)*(v(j+1,k)-v(j-1,k))*cy;
        v_f(j,k)=-a-b;        
    end
end