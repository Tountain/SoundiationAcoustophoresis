function cone3(X1, X2, r, color)
%%
% plot the cone part of the arrow
% arrow's start point 'X1';
% arrow's ending point 'X2';
% the radius of the cone 'r';
% color of the cylinder 'color'.
%%

% height of the cone
length_cyl = norm(X2-X1);
[x,y,z] = cylinder(linspace(r,0,50),100);
z = z * length_cyl;
% plot the bottom of the cone
hold on;
EndPlate1 = fill3(x(1,:),y(1,:),z(1,:),'r');
Cylinder = mesh(x,y,z);
% the rotational angle of the cone
unit_V=[0 0 1];
angle_X1X2=acos(dot( unit_V,(X2-X1) )/( norm(unit_V)*norm(X2-X1)) )*180/pi;
% rotational axis
axis_rot=cross(unit_V,(X2-X1));
% rotate the cylinder to desirable direction
if angle_X1X2~=0 % Rotation is not needed if required direction is along X
    rotate(Cylinder,axis_rot,angle_X1X2,[0 0 0])
    rotate(EndPlate1,axis_rot,angle_X1X2,[0 0 0])
end
% translate the cylinder to desirable position
set(EndPlate1,'XData',get(EndPlate1,'XData')+X1(1))
set(EndPlate1,'YData',get(EndPlate1,'YData')+X1(2))
set(EndPlate1,'ZData',get(EndPlate1,'ZData')+X1(3))
set(Cylinder,'XData',get(Cylinder,'XData')+X1(1))
set(Cylinder,'YData',get(Cylinder,'YData')+X1(2))
set(Cylinder,'ZData',get(Cylinder,'ZData')+X1(3))
% the color of the cylinder
set(Cylinder,'FaceColor',color)
set(EndPlate1,'FaceColor',color)
set(Cylinder,'EdgeAlpha',0)
set(EndPlate1,'EdgeAlpha',0)
% axis equal;
% view(3)

