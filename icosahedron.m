function [vertice]=icosahedron(iter_times)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    t=(1+sqrt(5))/2;
    e=5^(1/4);
    a=sqrt(t)/5^(1/4);
    b=1/(sqrt(t)*5^(1/4));
    c=a+2*b;
    d=a+b;

    vertice = [
      0,  a,  b;
      0, -a,  b;
      0,  a, -b;
      0, -a, -b;
      b,  0,  a;
      b,  0, -a;
     -b,  0,  a;
     -b,  0, -a;
      a,  b,  0;
      a, -b,  0;
     -a,  b,  0;
     -a, -b,  0];

    surface_midpoints = [
      d,  d,  d;
      d,  d, -d;
      d, -d,  d;
      d, -d, -d;
     -d,  d,  d;
     -d,  d, -d;
     -d, -d,  d;
     -d, -d, -d;
      0,  a,  c;
      0,  a, -c;
      0, -a,  c;
      0, -a, -c;
      c,  0,  a;
      c,  0, -a;
     -c,  0,  a;
     -c,  0, -a;
      a,  c,  0;
      a, -c,  0;
     -a,  c,  0;
     -a, -c,  0
   ]/3;
    
    for t=1:iter_times
        surf_size = size(surface_midpoints,1);
        new_surf_midp = rand(0,3);
        new_vert = rand(0,3);
        for i=1:surf_size

            d=(vertice(:,1)-surface_midpoints(i,1)).^2+(vertice(:,2)-surface_midpoints(i,2)).^2+(vertice(:,3)-surface_midpoints(i,3)).^2;
            [~,order]=sort(d);
            surf_vert=vertice(order(1:3),:);
            
            p3 = (surf_vert(1,:)+surf_vert(2,:))./2;
            p1 = (surf_vert(2,:)+surf_vert(3,:))./2;
            p2 = (surf_vert(3,:)+surf_vert(1,:))./2;
            p1 = p1./norm(p1);
            p2 = p2./norm(p2);
            p3 = p3./norm(p3);
            
            
            new_vert = unique([new_vert; p1; p2; p3],'rows');
            m1 = (surf_vert(1,:)+p2+p3)/3;
            m2 = (surf_vert(2,:)+p3+p1)/3;
            m3 = (surf_vert(3,:)+p1+p2)/3;
            m4 = (p1+p2+p3)/3;
            new_surf_midp = unique([new_surf_midp; m1; m2; m3; m4],'rows'); 
        end
        surface_midpoints = new_surf_midp;
        vertice = unique([vertice;new_vert],'rows');
    end

    disp(size(vertice,1));    
    vertice = vertice(vertice(:,3)>=0,:);
    disp(size(vertice,1));
    
end

