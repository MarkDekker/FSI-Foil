%  Quick script to test FEM approach
clc; clear;

%% Define Parameters and Constants for Computation
% Geometrical Parameters (Just to avoid errors - not used):
nose_angle  = 0.0;
chord       = 0.36;
top_layer   = 1.047;
bot_layer   = 1.0075;
Shift_SP(1) = 0;  %how many nodes should top separation point be shifted toward nose  
Shift_SP(2) = 0; %how many nodes should bottom separation point be shifted away from nose 

% Structural Constants:
E  = 1*10^9;% Pa
v  = 0.3;
A  = 0.00008;% m^2  (Taken as thickness in m by 1m width - Scaled by chord)
                         
%% Instantiate objects
FEM         = CALFEMNonL;

Couple      = Coupling;

Geom        = Geometry(nose_angle,chord,top_layer,bot_layer,Shift_SP);
%% Choose appropriate geometry
p1 = [-0.5 0];
p2 = [ 0.5 0];
r  = 0.5; 

Nodes = 100;
Shape = Geom.Catenary(p1,p2,r,Nodes);
Shape = Shape';

Lengths = ((Shape(1:end-1,1)-Shape(2:end,1)).^2+...
           (Shape(1:end-1,2)-Shape(2:end,2)).^2).^0.5;
f_rel   =  1/2*(Lengths(2:end,1)+Lengths(1:end-1,1));

% Define forces
theta = 270; %Force direction in degrees
F  = 1/mean(f_rel)*[f_rel*cosd(theta), f_rel*sind(theta)];
F  = [0 0; F; 0 0];
                      
%% Computation
FEM.CatenaryRun(Shape,F,[A E]);       

% Find Catenary of length pi*r
error = 1; error_1 = error;
a = 1; c=0; s=pi*r;
a_1 = a;

while abs(error)>0.0000001
    a_2 = a_1;
    a_1 = a;

    c = c+1;
    v = a*cosh(p2(1)/a)-a*cosh(p1(1)/a); % we know this is zero now
    PART_1 = sqrt(s^2-v^2);
    PART_2 = 2*a*sinh((p2(1)-p1(1))/(2*a));
    error  = PART_2-PART_1;
    
    error_2 = error_1;
    error_1 = error;
    
    if c<=2
        a = a-0.3;
    else
        a = a_1 - error_1*(a_1-a_2)/(error_1-error_2);
    end
end

%% Plot Results
Final_Shape = [[FEM.Ex(:,1);FEM.Ex(end,2)], ...
               [FEM.Ey(:,1);FEM.Ey(end,2)]];
s
Length = sum(((Final_Shape(1:end-1,1)-Final_Shape(2:end,1)).^2+...
              (Final_Shape(1:end-1,2)-Final_Shape(2:end,2)).^2).^0.5)

figure(2)
hold on
x = p1(1):0.01:p2(1);
plot(x,a*cosh(x./a)-a*cosh(p1(1)/a),':k')
plot(Shape(:,1),Shape(:,2),'k')
%plot(Final_Shape(:,1),Final_Shape(:,2),'b')
plot([FEM.Ex(:,1);FEM.Ex(end,2)], ...
     [FEM.Ey(:,1);FEM.Ey(end,2)],'g')
quiver([FEM.Ex(:,1);FEM.Ex(end,2)], ...
       [FEM.Ey(:,1);FEM.Ey(end,2)], ...
       FEM.r(1:3:end,1), FEM.r(2:3:end,1),0.25,'r','filled')
hold off
set(gca,'DataAspectRatio',[1 1 1])