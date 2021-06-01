function [uex qex]=exactvelocity(nx,ny)
dy=2/(ny-1);
y=-1:dy:1;
for j=1:ny
    for k=1:nx
        uex(j,k)=3/2*(1-y(j)^2);        
    end
end
sum=0;
for j=1:ny
    sum=sum+uex(j,1);
    for k=1:nx
        qex(k)=sum*dy;        
    end
end