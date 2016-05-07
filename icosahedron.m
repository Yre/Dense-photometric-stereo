function icosahedron()
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
global vertice;
global surface_midpoints;
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
    
    for t=1:4
        vert_size = size(vertice,1);
        surf_size = size(surface_midpoints,1);
%         max_dot = 0;
%         for i=1:vert_size
%             if (max_dot < dot(vertice(i,:),surface_midpoints(1,:)))
%                 max_dot = dot(vertice(i,:),surface_midpoints(1,:));
%             end
%         end
%         min_dis = 1;
%         for i=1:vert_size
%             if (min_dis>norm(vertice(i,:)-surface_midpoints(1,:)))
%                 min_dis=norm(vertice(i,:)-surface_midpoints(1,:));
%             end
%         end
%         fprintf('min_dis = %.3f\n',min_dis);
        new_surf_midp = rand(0,3);
        for i=1:surf_size
            surf_vert=rand(3,3);
%             disp('the surface midpoint is');
%             disp(surface_midpoints(i,:));
%             for j=1:vert_size
%                 fprintf('vertice[%d] %.3f %.3f %.3f %.4f\n',j,vertice(j,1:3), dot(vertice(j,:),surface_midpoints(i,:)));
%                 if (dot(vertice(j,:),surface_midpoints(i,:))==max_dot)
%                     surf_vert = [surf_vert; vertice(j,:)];
%                 end
%                 fprintf('vertice[%d] %.4f\n',j, norm(vertice(j,:)-surface_midpoints(i,:)));
%                 disp(vertice(j,:));
%                 if (norm(vertice(j,:)-surface_midpoints(i,:))==min_dis)
%                     surf_vert = [surf_vert; vertice(j,:)];
%                 end
%                 if (size(surf_vert,1)==3)
%                     break;
%                 end
%             end

            d=(vertice(:,1)-surface_midpoints(i,1)).^2+(vertice(:,2)-surface_midpoints(i,2)).^2+(vertice(:,3)-surface_midpoints(i,3)).^2;
            [~,order]=sort(d);
            surf_vert=vertice(order(1:3),:);
            
            p3 = (surf_vert(1,:)+surf_vert(2,:))./2;
            p1 = (surf_vert(2,:)+surf_vert(3,:))./2;
            p2 = (surf_vert(3,:)+surf_vert(1,:))./2;
            p1 = p1./norm(p1);
            p2 = p2./norm(p1);
            p3 = p3./norm(p1);
            
            
            vertice = unique([vertice; p1; p2; p3],'rows');
%             fprintf('now vertice size: %d\n',size(vertice,1));
            m1 = (surf_vert(1)+p2+p3)/3;
            m2 = (surf_vert(2)+p3+p1)/3;
            m3 = (surf_vert(3)+p1+p2)/3;
            m4 = (p1+p2+p3)/3;
            new_surf_midp = unique([new_surf_midp; m1; m2; m3; m4],'rows'); 
        end
        surface_midpoints = new_surf_midp;
        disp('size of surface_midpoints:');
        disp(size(surface_midpoints));
    end

    disp(size(vertice,1));    
    vertice = vertice(vertice(:,3)>=0,:);
    disp(size(vertice,1));
    
end

