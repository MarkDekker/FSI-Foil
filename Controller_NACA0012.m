%----------------------------------------------------------------------------------------%
%                                                                                        %
%               Master Program for a double-layered sail FSI Simulation                  %
%                            by Mark Dekker (January 2013)                               %
%                              written in MATLAB R2010a                                  %
%                                                                                        %
%----------------------------------------------------------------------------------------%

clc; clear;

% Fluid Constants:
Re     =  6000000;              % Reynolds number
Ma     =  0.15;                 % Mach Number
X_Tr   =  1;                 % Transition location (top and bottom)
N      =  9;          % Log of amplification factor in XFOIL: affects transition & turbulence

%% Instantiate objects
Geom        = Geometry(0,1,1,1,[0 0]);

Solution    = Results;
Solution.TitlePath;
Solution.readConstants(0,0,0,Re,N,0,0,0,0);
 
CFD           = XFOIL;
CFD.KeepFiles = false;  % Set it to false to delete all intermediate files created by XFOIL
CFD.Visible   = false;  % Set it to false to hide XFOIL plotting window

%% Create a NACA 4-series airfoil
Nodes = 280; 
NACA  = 0012;
Uni   = 2;
Geom.NACA0012Gen(NACA,Nodes,Uni);
Solution.conform{1,1}   = Geom.conform;
Solution.Airfoil{1,1}   = Geom.Airfoil;
Solution.Normals{1,1}   = Geom.Normals;
Solution.Tangents{1,1}  = Geom.Tangents;

%% Define Parameters for Computation
L80    = importdata('./ValidationData/Ladson80.mat');  %Angles of attack from validation case
Alpha  = L80(:,1)';

%% Computation
h = waitbar(0,'Running Computation...');

tStart=tic;
for i = 1:length(Alpha)
    CFD.addActionSet(0,Alpha(i),Re,Ma,N,0.01,X_Tr,X_Tr)
    Progress = i/length(Alpha);
    CFD.run                      % Now we're ready to run XFOIL
    waitbar(Progress)
    finished = CFD.wait(100);    %Wait up to 100 seconds for it to finish...
    if finished
    else
        CFD.kill;
    end
    Solution.readXFOIL(i,1,Alpha(i));
    Solution.Iter(i) = 1;
    Solution.Airfoil{i+1,1} = Solution.Airfoil{i,1};
    Solution.Normals{i+1,1} = Solution.Normals{i,1};
end
Solution.tElapsed=toc(tStart);
close(h)

disp('');
disp(repmat('-',1,70));
disp(['Analysis Completed in ',num2str(Solution.tElapsed),' seconds.'])
disp(repmat('-',1,70));
disp('');


%% Save Results File and Delete Working files
save(Solution.path,'Solution')

delete('Polar.txt');
delete('Dump.txt');
delete('Pressure.txt');
%delete('Current_Airfoil.dat');
%delete('Actions.txt');

%% Plot Results
Post = PlotTool;

figure(1)
Post.NACA0012Results(Solution);
%%
alp =5;
figure(2)
hold on
Post.PlotBaseProfile(Solution);
Post.PlotBL(Solution,alp,1,[1 1 1]*0.5,2,2);
hold off

figure(3)
plot(Solution.Airfoil{alp,1}(1:end/2,1),Solution.Cf{alp,1}(1:length(Solution.Airfoil{alp,1})/2),'b',...
     Solution.Airfoil{alp,1}(end/2:end,1),Solution.Cf{alp,1}(length(Solution.Airfoil{alp,1})/2:length(Solution.Airfoil{alp,1})),'g');

