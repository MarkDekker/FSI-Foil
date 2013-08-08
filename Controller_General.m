%----------------------------------------------------------------------------------------%
%                                                                                        %
%               Master Program for a double-layered sail FSI Simulation                  %
%                            by Mark Dekker (January 2013)                               %
%                              written in MATLAB R2010a                                  %
%                                                                                        %
%----------------------------------------------------------------------------------------%

clc; clear;

%% Define Parameters and Constants for Computation
% Geometrical Parameters:
nose_angle  = 12.5;
chord       = 0.36;
top_layer   = 1.052;
bot_layer   = 1.008;
Shift_SP(1) = 0; %how many nodes should top separation point be shifted toward nose  
Shift_SP(2) = 0; %how many nodes should bottom separation point be shifted away from nose 
X_Tr        = 1; %Set transition location for top layer

Alpha  =   [1 3 5 7 9];        % Angles of Attack

Cp_int =   -[0.19 0.225 0.28 0.33 0.38];     % Pressure Coefficient within sail


% Fluid Constants:
V      =  24.78 %m/s  		   Wind Speed - Drives Reynolds Number
mu     =  1.82128*10^-5;%kg/m*s  (At 17C and Sea Level - No humidity)
rho    =  1.225;%kg/m^3          (At 15C and Sea Level - No humidity)
Re     =  V*chord*rho/mu ;       % Reynolds number
Ma     =  0.1;                   % Mach Number
N      =  9;      		   % Log of amplification factor in XFOIL: affects transition & turbulence

% Structural Constants:
E      = 1*10^9;%Pa
v      = 0.3;
Area   = 0.000082/chord;%m^2  (Taken as thickness in m by 1m width - Scaled by chord)
                         
%% Instantiate objects
FEM         = CALFEMNonL;

Couple      = Coupling;

Geom        = Geometry(nose_angle,chord,top_layer,bot_layer,Shift_SP);

Solution    = Results;
Solution.TitlePath;
Solution.readConstants(V,rho,mu,Re,N,Cp_int,v,Area,E,chord);
 
CFD           = XFOIL;
CFD.KeepFiles = false;  % Set it to false to delete all intermediate files created by XFOIL
CFD.Visible   = false;  % Set it to false to hide XFOIL plotting window

%% Choose appropriate airfoil geometry
MastType   = 'Conformed NACA';
Nodes      = 180;  
Scale      = 1;
Geom.Generate(MastType,Nodes,0012,Scale);
Solution.readGeometry(1,1,Geom);
                      
%% Computation
h = waitbar(0,'Running Computation...');    last = false;     Progress = 0;
SimStart = tic;         % Start measuring total execution time
AirfoilLE = find(Solution.Airfoil{1,1}(:,1)==0);

for i = 1:length(Alpha)
    AlphaStart = tic;   % Start measuring execution time per Alpha
    j=0; uncon = 0;
    
    while j<100
        j = j+1;
        waitbar(Progress)
        % Run CFD    
        CFD.addActionSet(0,Alpha(i),Re,Ma,N,0.01,X_Tr,1,200); % Run Standard Action Set in XFOIL
        % (see class file for other options)
        CFD.run        % Now we're ready to run XFOIL
        
        
        
        % Wait up to 100 seconds for it to finish, otherwise kill process and log.
        finished = CFD.wait(200);
        if ~finished;
            CFD.kill;
            Solution.Unconv{i,j} = 'The CFD solver did not converge';
            disp(['The CFD solver did not converge at i=',num2str(i),...
                ' and j=',num2str(j),'.'])
        end
        
        Solution.readXFOIL(i,j,Alpha(i),Couple);
        
        if Solution.Unconv{i,j}==0
            % Update internal Cp
            Cp_i = (mean(Solution.Cp{i,j}(1:102,1))+mean(Solution.Cp{i,j}(102:end,1)))/2 ...
                - 0.00;
            Solution.Cp_int(1,i) = Cp_i;
        end
        
        if Solution.Unconv{i,j} == 0 % Check if XFOIL converged
            if j>1
                Geom.Separation(Solution,i,j) 
            end
            % Run FEM
            FEM.run(Solution,i,j);

            %Prepare Geometry for next iteration
            %Solution.U{i,j}  = Solution.U_fem{i,j};
            Couple.relaxU(Solution,i,j);            % Aitken relaxation            
            Couple.MeasureLyrs(Solution,i,j);       % Measure layer lengths after def.
            Couple.updateAirfoil(Solution,Geom,i,j);% Update airfoil definition
            Couple.convCheck(Solution,i,j);         % Check whether FSI sim. has converged
            
            % Record latest deformation in case of divergence in the next iteration. This 
            % can then be used to keep the simulation running.
            Solution.U_spare  =  Solution.U{i,j};   

            % Interrupt loop if Converged
            if j>2
                conv = sum(Solution.meanDev(i,j-2:j));
                Solution.Iter(i) = j;
                if conv<10^-4; Couple.newAlpha(Solution,Geom,i,j); break; end
            end
            
        else    % Alter geometry to try again (Adds robustness)
            uncon = uncon+1;    % Counter for unconverged solutions at this oper. point
            
            if (uncon>=2 && j<=3) || uncon>=6  % Trigger if XFOIL repeatedly diverges
                Solution.Unconv{i,j} = ['The CFD solver repeatedly failed to ',...
                                        'converge -> trying new value for alpha.'];
                
                % Use different alpha - obtained from random number generator.  
                % (In the vicinity of the original alpha of interest)
                Alpha(i) = Alpha(i)+random('Normal',0,0.2);
                
                disp(['The CFD solver repeatedly failed converge at i=',num2str(i),...
                      ' and j=',num2str(j),', hence trying an alternative alpha (',...
                        num2str(Alpha(i)),').'])
                    
                j=0;  uncon = 0;
                
            else % Simply alter geometry slightly to try again 
                Couple.reEval(Geom,Solution,i,j);
                Couple.convCheck(Solution,i,j);
            end
        end
        Progress                = (i-1+j/10)/length(Alpha);
    end

    Solution.tElapsed(i)    = toc(AlphaStart);
    save(Solution.path, 'Solution')  % Save throughout in case of program failure
end

Solution.tElapsed(i+1)=toc(SimStart);
close(h);

disp('');
disp(repmat('-',1,70));
disp(['Analysis Completed in ',num2str(Solution.tElapsed(i+1)),' seconds.'])
disp(repmat('-',1,70));
disp('');

%% Save Results File
save(Solution.path,'Solution')

disp('');
disp(repmat('-',1,70));
disp(['Results file saved as ',Solution.descriptor{1,1},'.m'])
disp(repmat('-',1,70));
disp('');

disp('');
disp(repmat('-',1,70));
disp(['Top layer length ',num2str(Solution.Length{1,1}(end,1)),...
      ' and engineering strain ',num2str((Solution.Length{1,1}(end,1) ...
      - Solution.Length{1,1}(1,1)) / Solution.Length{1,1}(1,1)),'.'])
disp(['Bottom layer length ',num2str(Solution.Length{1,1}(end,2)),...
      ' and engineering strain ',num2str((Solution.Length{1,1}(end,2) ...
      - Solution.Length{1,1}(1,2)) / Solution.Length{1,1}(1,2)),'.'])
disp(repmat('-',1,70));
disp('');

%delete('Polar.txt');
%delete('Dump.txt');
%delete('Pressure.txt');
%delete('Current_Airfoil.dat');
delete('Actions.txt');

%% Plot Results

Post = PlotTool;

% for i = 1:8
%     if i==2 || i==4 || i>=6
%      Post.vdBProfileCompare(Solution,i,i)
%     end
% end

figure(2)
Post.vdBForceCompare(Solution)

figure(1)
Post.ProfileOverlay(Solution);

