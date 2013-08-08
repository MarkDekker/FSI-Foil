 classdef Results < handle
    properties
        
        %-----------------------------------------------------------------%
        % Properties are all structured in a similar manner.              %
        % However some properties need to be recorded more often than     %
        % others.  Therefore the following system is adhered to:          %
        % 'Property'{i,j}(k,l) or {i,j} *prop(k,l)                        %
        % Where:
        %           i - angle of attack index                             %
        %           j - iteration index                                   %
        %           top - top or bottom half of profile                   %
        %           k - index for recorded property (eg pressure at point)%
        %           l - index for new set of properties (eg upper/lower)  %
        %           *prop - property name cell (eg alpha or itera)        %
        %-----------------------------------------------------------------%

        % ----------------   Naming and Current Simulation Description   --------------- %
        
        descriptor  % Contains inputs from user about Title, Subfolder and Comments (Saved
                    % as a data structure)
        path        % Contains original path where file was saved
        
        % --------------------   Geometric/Operational Properties   -------------------- %
        
        Alphas      % List angles of attack                     (i,j)
                    % (include alpha adjustments due to geometry changes)
        Iter        % No of iterations at for each alpha        {i}(k)
        Mast        % X- & Y-Coordinates of mast at alpha       {i}{1 or 2}(k,l)
        Layers      % X- & Y-Coordinates of sail layers         {i,j}{1 or 2}(k,l)
        Airfoil     % X- & Y-Coordinates of effective profile   {i,j}(k,l)
        AirfoilXF   % Coordinates produced by XFOIL
        CurAirfoil  % Current Airfoil for FEM
        Sep_Point   % Point where sails separate from mast      {i,j}{top}
        s           % Secant arc length parameters              {i,j}(k)
        Unconv      % Unconverged solutions                     *prop
        conform     % Normalising factor and rotation angle     {i,j}(k)
        Areas       % Area corresponding to each node (forces)  {i,j}(k,l)
        Tangents    % Tangent Vectors at each Node              {i,j}(k,l)
        Normals     % Normal Vectors at each Node               {i,j}(k,l)
        tElapsed    % Time elapsed                              {i}
        Length      % Length of top and bottom layer            {i,j}
        chord       % Initial chord length
        
        % ------------------------   Aerodynamic Properties -----------------------------%
        
        V           % Freestream Velocity [m/s]
        rho         % Density             [kg/m^3]
        mu          % Dynamic Viscosity   [kg/m*s]
        Ncrit       % Amplification factor for turbulent transition
        Re          % Reynolds Number (Based on expected chord length)
        CL          % Lift coefficient                          {i,j}
        CD          % Drag coefficient                          {i,j}
        CDp         % Drag coefficient due to pressure          {i,j}
        CM          % Moment coefficient                        {i,j}
        Top_Xtr     % Top transition point                      {i,j}
        Bot_Xtr     % Bottom transition point                   {i,j}
        Vrat        % Ratio between Ue (BL edge velocity) and   {i,j}(k)
                    % Vinf (freestream velocity)
        Theta       % Momentum thickness                        {i,j}(k)
        Cf          % Friction coefficient defined as           {i,j}(k)
                    % Cf = tau / 0.5 rho Qinf^2
        H           % Shape parameter                           {i,j}(k)
        Cp          % Pressure coefficient                      {i,j}(k)
        Cp_int      % Internal pressure coefficient
        F           % Forces exerted at every node              {i,j}(k,l)
        %Dstar       
        
        % -------------------------   Structural Properties -----------------------------%
        
        Area_Cross  % Crosssectional Area of Sail (per meter)   
        E           % Young's modulus of sail material
        v           % Poisson's Ratio
        U           % Movement of each node in x and y          {i,j}(k,l)
        U_fem       % Displacements caculated in Struct. Sol.   {i,j}(k,l)
        U_spare     % Last U that worked before unconv.         (k,l)
        Strain      % Deformation of each element               {i,j}(k,l)
        Stress      % Stress present in each element            {i,j}(k,l)
        N           % Internal axial load (positive in tension) {i,j}
        TESpring    % Trailing edge spring properties           {i} *prop
        EF          % Element end-forces                        {i,j}(k,l)
        Reaction    % Reaction Forces at fixed points           {i,j}
        
        % --------------------------   Coupling Parameters ------------------------------%
        
        w           % Relaxation factor                         {i,j}{dof}
        maxDev      % Max deviation between consec. Airfoils    (i,j)
        meanDev     % Mean deviation between consec. Airfoils   (i,j)
        
    end
    
    methods
        
        %-------------------------- Class Constructor and Saving ------------------------%
        
        function this = Results
            % Create Dialog Box to name and describe current simulation
            % (Automatically creates a unique test case name for fast runs)
%             [y, m, d, h, mn, sec] = datevec(now);      %get Current Date and Time
%             
%             prompt      = {'File Name:','Subfolder:','Description:'};
%             dlg_title   = 'Name and Comments';
%             num_lines   = 1;
%             def = {['testFSI_',num2str(y),'-',num2str(m),'-',num2str(d),'_',num2str(h),...
%                     'h',num2str(mn),'m',num2str(round(sec)),'s'],'Results-Tests',...
%                     'FSI test case'};
%             this.descriptor  = inputdlg(prompt,dlg_title,num_lines,def,'on');
%             
%             this.path        = strcat('./',...
%                                 this.descriptor{2,1},'/',this.descriptor{1,1});
           
        end
        
        function TitlePath(this)
            [y, m, d, h, mn, sec] = datevec(now);      %get Current Date and Time
            
            prompt      = {'File Name:','Subfolder:','Description:'};
            dlg_title   = 'Name and Comments';
            num_lines   = 1;
            def = {['testFSI_',num2str(y),'-',num2str(m),'-',num2str(d),'_',num2str(h),...
                    'h',num2str(mn),'m',num2str(round(sec)),'s'],'Results-Tests',...
                    'FSI test case'};
            this.descriptor  = inputdlg(prompt,dlg_title,num_lines,def,'on');
            
            this.path        = strcat('./',...
                                this.descriptor{2,1},'/',this.descriptor{1,1});
        end

        %----------------------------- Read Module Outputs ------------------------------%
        
        function readConstants(this,V,rho,mu,Re,Ncrit,Cp_int,v,Area,E,chord)
            % Static air properties
            this.V      = V;    this.rho    = rho;   this.mu  = mu;     this.Re  = Re;
            this.Cp_int = Cp_int;   
            this.Ncrit = Ncrit;
            % Static sailcloth properties
            this.v      = v;    this.Area_Cross   = Area;               this.E   = E;
            this.chord  = chord;
        end
        
        function readGeometry(this,i,j,Geom)
            this.Mast{i,1}      = Geom.Mast{1,1};
            this.Mast{i,2}      = Geom.Mast{1,2};
            this.Layers{i,j}    = Geom.Layers;
            this.Airfoil{i,j}   = Geom.Airfoil;
            this.Sep_Point{i,j} = Geom.Sep_Point;
            this.conform{i,j}   = Geom.conform;
            this.Normals{i,j}   = Geom.Normals;
            this.Tangents{i,j}  = Geom.Tangents;
            this.N{i,j}         = zeros(length(this.Airfoil{i,j}),1);
        end
        
        function readXFOIL(this,i,j,alph,Couple)
            fid = fopen('Polar.txt');
            pdata=textscan(fid,'%f%f%f%f%f%f%f','HeaderLines',12);
            fclose(fid);
            
            if isempty(pdata{1,1})
                this.Unconv{i,j} = 'The CFD solver did not converge';
                disp(['The CFD solver did not converge at i=',num2str(i),...
                    ' and j=',num2str(j),'.'])
                this.Alphas(i,j)        = alph;
            else
                this.Unconv{i,j}        = 0;
                this.Alphas(i,j)        = pdata{1,1}(2,1);
                this.CL(i,j)            = pdata{1,2}(2,1);
                this.CD(i,j)            = pdata{1,3}(2,1);
                this.CDp(i,j)           = pdata{1,4}(2,1);
                this.CM(i,j)            = pdata{1,5}(2,1);
                this.Top_Xtr(i,j)       = pdata{1,6}(2,1);
                this.Bot_Xtr(i,j)       = pdata{1,7}(2,1);
                
                fid = fopen('Pressure.txt');
                ppdata=textscan(fid,'%f%f','HeaderLines',1);
                fclose(fid);
                
                this.Cp{i,j}        = ppdata{1,2};
                
                [l, ~]              = size(this.Cp{i,j});
                
                fid = fopen('Dump.txt');
                ddata=textscan(fid,'%f%f%f%f%f%f%f%f','HeaderLines',1);
                fclose(fid);
                
                %Nd = size(ddata{1,1},1);
                
                this.s{i,j}                 = ddata{1,1}(1:l,1);      %
                this.AirfoilXF{i,j}(:,1)    = ddata{1,2}(1:l,1); % is this necessary?
                this.AirfoilXF{i,j}(:,2)    = ddata{1,3}(1:l,1); % is this necessary?
                this.Vrat{i,j}              = ddata{1,4}(1:l,1);%
                this.Theta{i,j}             = ddata{1,6}(1:l,1);
                this.Cf{i,j}                = ddata{1,7}(1:l,1);
                this.H{i,j}                 = ddata{1,8}(1:l,1);
                
                % Unlikely to need this:
                % this.Dstar{i}     = ddata{1,5};
                
                if nargin>4
                    Couple.genForces(this,i,j)
                end
            end
            
        end
        
        function readFEMBeam(this,FEM,Results,i,j,layer)
            % Read x- & y- displacements as well as rotations of each node respectively.

            if layer == 1
                [len, ~] = size(Results.Airfoil{i,j});
                this.U_fem{i,j}  = zeros(len,3);
                this.U_fem{i,j}(1:Results.Sep_Point{i,j}{1},:)  = ...
                                                            [FEM.F(1:3:FEM.NN*3-2,1), ...
                                                             FEM.F(2:3:FEM.NN*3-1,1), ...
                                                             FEM.F(3:3:FEM.NN*3  ,1) ];
            else
                this.U_fem{i,j}(Results.Sep_Point{i,j}{2}:end,:)  = ...
                                                            [FEM.F(1:3:FEM.NN*3-2,1), ...
                                                             FEM.F(2:3:FEM.NN*3-1,1), ...
                                                             FEM.F(3:3:FEM.NN*3  ,1) ];
            end
            
            % Read member end-forces (forces/moments acting on the ends of each element)
            this.EF{i,j}{layer} = FEM.EF;
            
            % Read reaction forces (trailing edge is second)
            if layer == 1
                this.Reaction{i,j}{1} = [FEM.REACT(1,[3,1]) FEM.REACT(1,[4,2])];
            else
                this.Reaction{i,j}{2} = [FEM.REACT(1,[1,3]) FEM.REACT(1,[2,4])];
            end
                       
        end
        
        function readFEMTruss(this,FEM,Results,i,j,layer)
            % Read x- & y- displacements as well as rotations of each node respectively.

            if layer == 1
                [len, ~] = size(Results.Airfoil{i,j});
                this.U_fem{i,j}  = zeros(len,2);
                this.U_fem{i,j}(1:Results.Sep_Point{i,j}{1},:)  = ...
                                                            [FEM.F(1:2:FEM.NN*2-1,1), ...
                                                             FEM.F(2:2:FEM.NN*2  ,1)];
            else
                this.U_fem{i,j}(Results.Sep_Point{i,j}{2}:end,:)  = ...
                                                            [FEM.F(1:2:FEM.NN*2-1,1), ...
                                                             FEM.F(2:2:FEM.NN*2  ,1)];
            end
            
            % Read member end-forces (forces/moments acting on the ends of each element)
            this.Stress{i,j}{layer} = FEM.STRESS;
            
            % Read reaction forces (trailing edge is second)
            if layer == 1
                this.Reaction{i,j}{1} = [FEM.REACT(1,[3,1]) FEM.REACT(1,[4,2])];
            else
                this.Reaction{i,j}{2} = [FEM.REACT(1,[1,3]) FEM.REACT(1,[2,4])];
            end
                       
        end
    
        %--------------------- Log Errors and Unconverged Solutions ---------------------%
        
        function logUnconv(this,i,j,String)
            this.Unconv{i,j} = String;
            
        end
        
        function clearOpPoint(this,i)
            
            % Geometry results
            this.Airfoil{i}        = [];
            this.Sep_Point{i}      = [];
            this.Unconv{i}         = [];
            this.s{i}              = [];
            this.Areas{i}          = [];
            this.Tangents{i}       = [];
            this.Normals{i}        = [];
            this.Length{1,i}       = [];
            
            % Fluid results
            this.CL(i,1:end)       = 0;
            this.CD(i,1:end)       = 0;
            this.CDp(i,1:end)      = 0;
            this.CM(i,1:end)       = 0;
            this.Top_Xtr(i,1:end)  = 0;
            this.Bot_Xtr(i,1:end)  = 0;
            this.s{i}              = [];
            this.AirfoilXF{i}      = [];
            this.AirfoilXF{i}      = []; 
            this.Vrat{i}           = [];
            this.Theta{i}          = [];
            this.Cf{i}             = [];
            this.H{i}              = [];
            
            % Structure results
            this.F{i}              = [];
            this.U{i}              = [];
            this.U_fem{i}          = [];
            this.EF{i}             = [];
            this.Reaction{i}       = [];
            
            % FSI results
            this.w{i}              = [];
            this.meanDev(i,1:end)  = 0;
            this.maxDev(i,1:end)   = 0;
            
        end
        %---------------------- Infer Forces Acting on each Node ------------------------%
        
        function genForces(this,i,j)
            % Determine Area for each corresponding node. (The first and
            % last are excluded because they are fixed in the model.)
            Areas_temp = ((this.Airfoil{i,j}(2:end,1) - ...
                           this.Airfoil{i,j}(1:end-1,1)).^2 + ...
                          (this.Airfoil{i,j}(2:end,2) - ...
                           this.Airfoil{i,j}(1:end-1,2)).^2).^0.5;
                            
            this.Areas{i,j}(:,1) = 1/2*(Areas_temp(1:end-1,1) + Areas_temp(2:end,1));
            
            % Convert Pressure and Friction Coefficients into Forces
            [len, ~]            = size(this.Cp{i,j}(2:end-1));
            this.rho            = 1.225;%kg/m^3
            Normal_F            = this.Cf{i,j}(2:end-1) .* this.Areas{i,j} ...
                                    * 0.5 * this.rho * this.V^2 ;
            Tangen_F            = (this.Cp{i,j}(2:end-1) + this.Cp_int*ones(len,1))...
                                    .* this.Areas{i,j} * 0.5 * this.rho * this.V^2 ;
                                
            % Write Tangent and Normal Forces in terms of Cartesian
            % Coordinate system
            this.F{i,j}(:,1)    = Normal_F .* this.Normals{i,j}(2:end-1,1) + ...
                                  Tangen_F .* this.Tangents{i,j}(2:end-1,1);
            this.F{i,j}(:,2)    = Normal_F .* this.Normals{i,j}(2:end-1,2) + ...
                                  Tangen_F .* this.Tangents{i,j}(2:end-1,2);
        end
        
        function FEMInput(this,i,j,layer)
            SP = this.Sep_Point{i,j}{layer};
            
            if layer == 1
                [NN, ~]     = size(this.Airfoil{i,j}(1:SP,:));
                NODE        = [1:1:NN; this.Airfoil{i,j}(1:SP,:)']';
                %LOAD = zeros(2,NN*2);
                LOAD(:,1:2:NN*2-5) = [3:2:NN*2-3; this.F{i,j}(2:SP-1,1)'];
                LOAD(:,2:2:NN*2-4) = [4:2:NN*2-2; this.F{i,j}(2:SP-1,2)'];
            else
                [NN, ~]     = size(this.Airfoil{i,j}(SP:end,:));
                NODE        = [1:1:NN; this.Airfoil{i,j}(SP:end,:)']';
                %LOAD = zeros(2,NN*2);
                LOAD(:,1:2:NN*2-5) = [3:2:NN*2-3; this.F{i,j}(SP+1:end-1,1)'];
                LOAD(:,2:2:NN*2-4) = [4:2:NN*2-2; this.F{i,j}(SP+1:end-1,2)'];
            end
            
            
            
            NE = NN-1;
            NL = NN*2-4;
            
            First_Line  = [NN NE 1 2 2 2];
            Second_Line = [4 NL 0];
            
            ELEM        =  [1:1:NE
                            1:1:NN-1
                            2:1:NN    
                            ones(1,NE) 
                            this.Area_Cross * ones(1,NE)
                            zeros(1,NE)]';
                        
            PRE_DIS     = [1 2 NN*2-1 NN*2; zeros(1,4)]';
            
            PROP        = [1 this.E 12E-6];
                
            
            
            % Write intput to file 'Frame.in'
            fid = fopen('Frame.in','w');
            fprintf(fid,'<<2-D Truss Analysis>> \r\n');
            fprintf(fid,'Sailwing FSI \r\n');
            fprintf(fid,'NN  NE  NM   NDIM  NEN  NDN \r\n');
            fprintf(fid,'%3d %3d %3d %3d %3d %3d \r\n',First_Line');
            
            fprintf(fid,'ND NL  NMPC \r\n');
            fprintf(fid,'%2d %3d %4d \r\n',Second_Line');
            
            fprintf(fid,'Node#       X         Y\r\n');
            fprintf(fid,'%2d % 8.5g % 8.5g \r\n',NODE');
            
            fprintf(fid,'Elem#  N1   N2  Mat#  Area TempRise\r\n');
            fprintf(fid,'%2d   %2d  %2d  %2d   %8.3g   % 2d\r\n',ELEM');
            
            fprintf(fid,'DOF#   Displacement\r\n');
            fprintf(fid,'%2d %8.5d \r\n',PRE_DIS');
            
            fprintf(fid,'DOF#      Load \r\n');
            fprintf(fid,'%2d % 5e\n',LOAD);
            
            fprintf(fid,'MAT#   E       Alpha \r\n');
            fprintf(fid,'%2d % 5G % 5G \r\n',PROP');
            
            fprintf(fid,'B1 i   B2  j   B3  (Multi-point constr. B1*Qi+B2*Qj=B3) \r\n');
            fclose(fid);
        end
        
    end
end