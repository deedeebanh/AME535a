function [u,v,fu,fv] = CONVEC(u,v,nx,ny,dx,dy,re,dt,fu1,fv1,uin)
% Second-step and so on use 
% Adam Bashford to update 
% velocity u and v
fu = u_f(u,v,nx,ny,dx,dy,re);
fv = v_f(u,v,nx,ny,dx,dy);
u = u + dt*(3/2*fu-1/2*fu1); %Ustar for nonlinear step
v = v + dt*(3/2*fv-1/2*fv1); %Vstar for nonlinear step
u(2:ny-1,1)=uin;
end