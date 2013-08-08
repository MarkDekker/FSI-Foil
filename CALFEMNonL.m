%-------------------------------------------------------------------------%
% Class CALFEMNonL based on Non-Linear 2D Beams in open source CALFEM     %
% distribution.  Stiffness matrix is a function of internal loading.      %
% Problem solved iteratively, updating the stiffness matrix constantly,   % 
% with the option of updating applied forces too (perhaps to remain       %
% orthogonal).  System of equations solved with internal MATLAB system    %
% solvers (i.e "\" operation).                                            %
%-------------------------------------------------------------------------%

classdef CALFEMNonL < handle
    properties
        Edof    % Node Connectivity Matrix (Defines Elements)
        Ex      % X-locations of Nodes per Element
        Ex1
        Ey      % Y-locations of Nodes per Element
        Ey1
        Ep      % Section Properties
        Ed      % Calculated element displacement
        Es      % Internal loads of every node (local coord)
        g       % Internal Loads in global coordinates 
        r       % Difference between applied and internal force
        Q       % Cumulative displacement
        f       % Forces Applied at Nodes
        N       % Internal normal forces
        Bc      % Boundary conditions (prescribed displacements)
        Kmat    % Log Stiffness Matrix
        R       % Reaction Forces at Nodes
        REACT   % Reaction Forces as applied to Airfoil
        U       % Displacements as applied to the Airfoil
        SeP     % Separation Point
    end
    
    methods (Static)
        function   [K,f] = assem(edof,K,Ke,f,fe)
            % K=assem(edof,K,Ke)
            % [K,f]=assem(edof,K,Ke,f,fe)
            %-------------------------------------------------------------
            % PURPOSE
            %  Assemble element matrices Ke ( and fe ) into the global
            %  stiffness matrix K ( and the global force vector f )
            %  according to the topology matrix edof.
            %
            % INPUT: edof:       dof topology matrix
            %          K :       the global stiffness matrix
            %          Ke:       element stiffness matrix
            %          f :       the global force vector
            %          fe:       element force vector
            % OUTPUT:  K :       the new global stiffness matrix
            %          f :       the new global force vector
            %-------------------------------------------------------------

            [nie,n]=size(edof);
            t=edof(:,2:n);
            for i = 1:nie
                K(t(i,:),t(i,:)) = K(t(i,:),t(i,:))+Ke;
                if nargin==5
                    f(t(i,:))=f(t(i,:))+fe;
                end
            end
        end
        
        function [Ke,fe] = beam2g(ex,ey,ep,N,eq)
            % Ke=beam2g(ex,ey,ep,N)
            % [Ke,fe]=beam2g(ex,ey,ep,N,eq)
            %-------------------------------------------------------------
            %    PURPOSE
            %       Compute the element stiffness matrix for a two dimensional
            %       beam element with respect to geometric nonlinearity.
            %
            %    INPUT:  ex = [x1 x2]   element node coordinates
            %            ey = [y1 y2]
            %            ep = [E A I]   element properties;
            %                  E: Young's modulus
            %                  A: Cross section area
            %                  I: Moment of inertia
            %            N:  axial force in the beam.
            %            eq: distributed transverse load
            %    OUTPUT: Ke : element stiffness matrix (6 x 6)
            %            fe : element load vector (6 x 1)
            %-------------------------------------------------------------

            if nargin==4; eq=0; end
            
            if length(eq)>1
                disp('eq should be a scalar!!!')
                return
            end
            
            b=[ ex(2)-ex(1); ey(2)-ey(1) ];
            L=sqrt(b'*b);
            n=b/L;
            
            E=ep(1); A=ep(2); I=ep(3); rho=-N*L^2/(pi^2*E*I);
            
            kL=pi*sqrt(abs(rho));%+eps;
            
            if rho>0
                f1=(kL/2)/tan(kL/2);
                f2=(1/12)*kL^2/(1-f1);
                f3=f1/4+3*f2/4;
                f4=-f1/2+3*f2/2;
                f5=f1*f2;
                
                h=6*(2/kL^2-(1+cos(kL))/(kL*sin(kL)));
                
            elseif rho<0
                f1=(kL/2)/tanh(kL/2);
                f2=-(1/12)*kL^2/(1-f1);
                f3=f1/4+3*f2/4;
                f4=-f1/2+3*f2/2;
                f5=f1*f2;
                
                h=-6*(2/kL^2-(1+cosh(kL))/(kL*sinh(kL)));
            else
                f2=1;f3=1;f4=1;f5=1;h=1;%f1=1;
            end
            
            Kle=[E*A/L  0            0        -E*A/L      0          0 ;
                0  12*E*I*f5/L^3   6*E*I*f2/L^2  0 -12*E*I*f5/L^3 6*E*I*f2/L^2;
                0  6*E*I*f2/L^2    4*E*I*f3/L    0  -6*E*I*f2/L^2   2*E*I*f4/L;
                -E*A/L  0            0         E*A/L      0          0 ;
                0  -12*E*I*f5/L^3 -6*E*I*f2/L^2  0  12*E*I*f5/L^3 -6*E*I*f2/L^2;
                0  6*E*I*f2/L^2    2*E*I*f4/L    0  -6*E*I*f2/L^2   4*E*I*f3/L];
            
            fle=eq*L*[0 1/2 L*h/12 0 1/2 -L*h/12]';
            
            
            G=[n(1) n(2)  0    0    0   0;
                -n(2) n(1)  0    0    0   0;
                0    0    1    0    0   0;
                0    0    0   n(1) n(2) 0;
                0    0    0  -n(2) n(1) 0;
                0    0    0    0    0   1];
            
            Ke=G'*Kle*G;  fe=G'*fle;
        end
        
        function      es = beam2gs(ex,ey,ep,ed,N,eq)
            % es=beam2gs(ex,ey,ep,ed,N,eq)
            %-------------------------------------------------------------
            %    PURPOSE
            %      Calculate section forces in a two dimensional nonlinear
            %      beam element.
            %
            %    INPUT:  ex = [x1 x2]
            %            ey = [y1 y2]          element node coordinates
            %            ep = [E A I]          element properties,
            %                  E: Young's modulus
            %                  A: cross section area
            %                  I: moment of inertia
            %            ed = [u1 ... u6]       element displacement vector
            %            N			            axial force
            %            eq = [qy]              distributed transverse load
            %    OUTPUT: es = [N1 V1 M1 ;
            %                  N2 V2 M2 ]       element forces, local directions
            %-------------------------------------------------------------

            b=[ex(2)-ex(1); ey(2)-ey(1)];
            L=sqrt(b'*b);  n=b/L;
            %
            if nargin==5; eq=0; end
            %
            E=ep(1); A=ep(2); I=ep(3); rho=-N*L^2/(pi^2*E*I);
           
            kL=pi*sqrt(abs(rho))+eps;
            
            if rho>0
                f1=(kL/2)/tan(kL/2);
                f2=(1/12)*kL^2/(1-f1);
                f3=f1/4+3*f2/4;
                f4=-f1/2+3*f2/2;
                f5=f1*f2;
                
                h=6*(2/kL^2-(1+cos(kL))/(kL*sin(kL)));
                
            elseif rho<0
                f1=(kL/2)/tanh(kL/2);
                f2=-(1/12)*kL^2/(1-f1);
                f3=f1/4+3*f2/4;
                f4=-f1/2+3*f2/2;
                f5=f1*f2;
                
                h=-6*(2/kL^2-(1+cosh(kL))/(kL*sinh(kL)));
            else
                f2=1;f3=1;f4=1;f5=1;h=1;%f1=1;
            end
            
            Kle=[E*A/L  0            0        -E*A/L      0          0 ;
                0  12*E*I*f5/L^3   6*E*I*f2/L^2  0 -12*E*I*f5/L^3 6*E*I*f2/L^2;
                0  6*E*I*f2/L^2    4*E*I*f3/L    0  -6*E*I*f2/L^2   2*E*I*f4/L;
                -E*A/L  0            0         E*A/L      0          0 ;
                0  -12*E*I*f5/L^3 -6*E*I*f2/L^2  0  12*E*I*f5/L^3 -6*E*I*f2/L^2;
                0  6*E*I*f2/L^2    2*E*I*f4/L    0  -6*E*I*f2/L^2   4*E*I*f3/L];
            
            fle=eq*L*[0 1/2 L*h/12 0 1/2 -L*h/12]';
            
            G=[n(1) n(2)  0    0    0   0;
                -n(2) n(1)  0    0    0   0;
                0    0    1    0    0   0;
                0    0    0   n(1) n(2) 0;
                0    0    0  -n(2) n(1) 0;
                0    0    0    0    0   1];
            
            u=ed';
            P=(Kle*G*u-fle);
            es=[-P(1) -P(2) -P(3)
                 P(4)  P(5)  P(6)];
        end
        
        function    [ed] = extract(edof,a)
            %-------------------------------------------------------------
            % PURPOSE
            %  Extract element displacements from the global displacement
            %  vector according to the topology matrix edof.
            %
            % INPUT:     a:  the global displacement vector
            %         edof:  topology matrix
            % OUTPUT:   ed:  element displacement matrix
            %-------------------------------------------------------------
            
            [nie,n]=size(edof);

            t=edof(:,2:n);
            ed = zeros(nie,n-1);
            for i = 1:nie
                ed(i,1:(n-1))=a(t(i,:))';
            end
        end
        
        function   [d,Q] = solveq(K,f,bc)
            % a=solveq(K,f)
            % [a,Q]=solveq(K,f,bc)
            %-------------------------------------------------------------
            % PURPOSE
            %  Solve static FE-equations considering boundary conditions.
            %
            % INPUT:   K : global stiffness matrix, dim(K)= nd x nd
            %          f : global load vector, dim(f)= nd x 1
            %         bc : boundary condition matrix
            %              dim(bc)= nbc x 2, nbc : number of b.c.'s
            % OUTPUT:  a : solution including boundary values
            %          Q : reaction force vector
            %              dim(a)=dim(Q)= nd x 1, nd : number of dof's
            %-------------------------------------------------------------
            
            if nargin==2;
                d=K\f ;
            elseif nargin==3;
                [~,nd]=size(K);
                fdof = 1:nd;

                d = zeros(size(fdof'));
                %Q = zeros(size(fdof'));

                pdof       = bc(:,1);
                dp         = bc(:,2);
                fdof(pdof) = [];

                s = K(fdof',fdof')\(f(fdof')-K(fdof',pdof)*dp);

                d(pdof)  = dp;
                d(fdof') = s;
            end
            Q=K*d-f;
            
        end
    end
    
    methods
        function this = CALFEMNonL
        end
        
        function run(this,Results,i,j)
            this.U      = zeros(length(Results.Airfoil{i,j}),2);
            this.SeP(1) = Results.Sep_Point{i,j}{1};
            this.SeP(2) = Results.Sep_Point{i,j}{2};
            for layer = 1:2
                this.InputData(Results,layer,i,j)
                this.NewtonRhapson
                this.Output(Results,layer,i,j)
            end
        end
        
        function InputData(this,Results,layer,i,j)
            % Read in top and bottom sail-mast separation points
            % Use base geometry where no internal stresses are present (iter 1)
            SP = this.SeP(layer);
            
            x = 1;
            y = 1;
            
            % -------  Define Geometry -------- %
            % Coordinates
            if layer == 1
                NN      = length(Results.Airfoil{x,y}(1:SP,:));
                NDOF = NN*3;        % Number of Degrees of Freedom
                NE   = NN-1;        % Number of Elements
                this.Ex = [Results.Airfoil{x,y}(1:SP-1,1),...
                           Results.Airfoil{x,y}(2:SP,1)];
                this.Ey = [Results.Airfoil{x,y}(1:SP-1,2),...
                           Results.Airfoil{x,y}(2:SP,2)];
                       
                % Allocate forces
                this.f             = zeros(NDOF,1);
                this.f(4:3:NDOF-2) = Results.F{i,j}(1:SP-1,1);
                this.f(5:3:NDOF-1) = Results.F{i,j}(1:SP-1,2);
                
                % Boundary Conditions
                this.Bc = [1 2 4 5 NDOF-2 NDOF-1
                           0 0 0 0 0      0     ]';
                
                this.Es = zeros(NE,6);
            else
                NN      = length(Results.Airfoil{x,y}(SP:end,:));
                NDOF    = NN*3;        % Number of Degrees of Freedom
                NE      = NN-1;        % Number of Elements
                this.Ex = [Results.Airfoil{x,y}(SP:end-1,1),...
                           Results.Airfoil{x,y}(SP+1:end,1)];
                this.Ey = [Results.Airfoil{x,y}(SP:end-1,2),...
                           Results.Airfoil{x,y}(SP+1:end,2)];
                       
                % Allocate forces
                this.f             = zeros(NDOF,1);
                this.f(1:3:NDOF-5) = Results.F{i,j}(SP-1:end,1);
                this.f(2:3:NDOF-4) = Results.F{i,j}(SP-1:end,2);
                
                % Boundary Conditions
                this.Bc = [1 2 NDOF-2 NDOF-1
                           0 0 0      0     ]';
                       
                this.Es = zeros(NE,6);
            end
            
            %Capture initial geometry
            this.Ex1 = this.Ex;
            this.Ey1 = this.Ey;

            % Elements, connectivity and node locations
            
            COL  = 1:3:NDOF-5;  % DOF in first column of Edof
            
            this.Edof = [1:1:NE
                         COL
                         COL+1
                         COL+2
                         COL+3
                         COL+4
                         COL+5 ]';
                     
            % Element properties
            A       = Results.Area_Cross;
            I       = 1/12*(Results.Area_Cross)^3;
            this.Ep = [Results.E A I];
            
            % Initialise internal forces to prevent initial divergence
            this.N = ones(NE,1);
        end
        
        function Output(this,Results,layer,i,j)
            % Prepare variables for output
            if     layer == 1
                this.U(2:this.SeP(1)-1,:) = [this.Ex(2:end,1),  this.Ey(2:end,1)] - ...
                                            [this.Ex1(2:end,1), this.Ey1(2:end,1)];
                Results.N{i,j}(1:this.SeP(1)-1,:) = this.N(:,1);
            elseif layer == 2
                this.U(this.SeP(2)+1:end-1,:) = [this.Ex(2:end,1), this.Ey(2:end,1)] - ...
                                                [this.Ex1(2:end,1),this.Ey1(2:end,1)];          
                Results.U_fem{i,j}            = this.U;
                Results.N{i,j}(this.SeP(2):end-1,:) = this.N(:,1);
            end
            
        end
        
        function G = InternalLoads(this,i)
            % Determine inner load vectors to check balance of forces
            
            % Obtain neighbouring node locations descriptions based on
            % current node

            actv1   = 0;
            actv2   = 0;
            
            if i>1 && i<length(this.Ex)+1
                p1 = [this.Ex(i-1,1) this.Ey(i-1,1)]; 
                p2 = [this.Ex(i,1)   this.Ey(i,1)] ;
                p3 = [this.Ex(i,2)   this.Ey(i,2)];
                
                L1      = sqrt((p2(1)-p1(1))^2+(p2(2)-p1(2))^2); 
                L2      = sqrt((p3(1)-p2(1))^2+(p3(2)-p2(2))^2); 
                
                % Obtain force component vectors for normal forces
                comp_x1 = (p1(1)-p2(1))/L1;
                comp_y1 = (p1(2)-p2(2))/L1;
                comp_x2 = (p3(1)-p2(1))/L2;
                comp_y2 = (p3(2)-p2(2))/L2;
                
                actv1   = 1;
                actv2   = 1;
                
                u_i1 = this.Ed(i-1,:);%[0, 0, 0, this.Ed(i-1,4:6)-this.Ed(i-1,1:3)];
                u_i  = this.Ed(i,:);%[0, 0, 0, this.Ed(i,4:6)-this.Ed(i,1:3)];
                
                es1 = this.beam2gs(this.Ex(i-1,:),this.Ey(i-1,:),this.Ep,u_i1,this.N(i-1));
                es2 = this.beam2gs(this.Ex(i,:),  this.Ey(i,:),  this.Ep,u_i,  this.N(i));
                
                this.Es(i,:) = [es2(1,:) es2(2,:)];
                
            elseif i==1
                p2 = [this.Ex(i,1)   this.Ey(i,1)]; 
                p3 = [this.Ex(i,2)   this.Ey(i,2)];
                 
                L2      = sqrt((p3(1)-p2(1))^2+(p3(2)-p2(2))^2); 
                
                % Obtain force component vectors for normal forces
                comp_x1 = 0;
                comp_y1 = 0;
                comp_x2 = (p3(1)-p2(1))/L2;
                comp_y2 = (p3(2)-p2(2))/L2;
                
                actv2   = 1;
                
                u_i  = this.Ed(i,:);%[0, 0, 0, this.Ed(i,4:6)-this.Ed(i,1:3)];
                
                es1 = [0 0 0
                       0 0 0];
                es2 = this.beam2gs(this.Ex(i,:),this.Ey(i,:),this.Ep,u_i,this.N(i));
                
                this.Es(i,:) = [es2(1,:) es2(2,:)];
                
            elseif i==length(this.Ex)+1
                p1 = [this.Ex(i-1,1) this.Ey(i-1,1)]; 
                p2 = [this.Ex(i-1,2) this.Ey(i-1,2)]; 
                
                L1      = sqrt((p2(1)-p1(1))^2+(p2(2)-p1(2))^2); 
                
                % Obtain force component vectors for normal forces
                comp_x1 = (p1(1)-p2(1))/L1;
                comp_y1 = (p1(2)-p2(2))/L1;
                comp_x2 =  0;
                comp_y2 =  0;
                
                actv1   = 1;
                
                u_i1 = this.Ed(i-1,:);%[0, 0, 0, this.Ed(i-1,4:6)-this.Ed(i-1,1:3)];
                
                es1 = this.beam2gs(this.Ex(i-1,:),this.Ey(i-1,:),this.Ep,u_i1,this.N(i-1));
                es2 = [0 0 0
                       0 0 0];
                
            end
            % Rotate force components to determine shear force
            % contributions
            compc_x1  = comp_x1 * cos(pi/2) - comp_y1 * sin(pi/2);
            compc_x2  = comp_x2 * cos(pi/2) - comp_y2 * sin(pi/2);
            compc_y1  = comp_x1 * sin(pi/2) + comp_y1 * cos(pi/2);
            compc_y2  = comp_x2 * sin(pi/2) + comp_y2 * cos(pi/2);
            
            gx  = comp_x1*es1(1,1) + comp_x2*es2(1,1) + (compc_x1*es1(1,2) + compc_x2*es2(1,2));
            gy  = comp_y1*es1(1,1) + comp_y2*es2(1,1) + (compc_y1*es1(1,2) + compc_y2*es2(1,2));
            gm  = - es1(2,3)*actv1 + es2(1,3)*actv2;
            
            G   = [gx gy gm];
            
        end
        
        %  ---------------------- Iterative Schemes -----------------------------------  %
        
        function Iterate(this)
            NE = length(this.Ex); NN = NE+1; NDOF = NN*3;
            
            this.Ed = zeros(NE,3*2);
            this.Q  = zeros(NE-1,3);
            
            for n=1:100
                if n<20; F = this.f/20*n;
                else F = this.f;
                end
                
                Ed0 = this.Ed;
                
                % Determine and Assemble Stiffness Matrix
                K=zeros(NDOF,NDOF);
                for i = 1:NE
                    Ke = this.beam2g(this.Ex(i,:),this.Ey(i,:),this.Ep,this.N(i));
                    K  = this.assem(this.Edof(i,:),K,Ke);
                end
                this.Kmat = K;
                
                [a,re]   = this.solveq(K,F,this.Bc); 
                
                this.Ed = this.extract(this.Edof,a);
                this.R  = this.extract(this.Edof,re);
                
                % Update internal loads
                Loads = zeros(NE);
                for i = 1:NE
                    es = this.beam2gs(this.Ex(i,:),this.Ey(i,:),this.Ep,....
                          this.Ed(i,:),this.N(i));
                      Loads(i) = es(1,1);
                end
                
                this.N = Loads;
                
                % Check for convergence
                if sum(sum(abs(this.Ed-Ed0)))<0.00001;
                    return;
                elseif n==100; disp('The solution has not converged')
                end
            end
            
            % Log Reaction Forces
            this.R = [re(1:3:NDOF-5), re(2:3:NDOF-4), re(3:3:NDOF-3), re(4:3:NDOF-2), ...
                re(5:3:NDOF-1), re(6:3:NDOF)];
        end
        
        function NewtonRhapson(this)
            NE = length(this.Ex);
            NN = NE+1; NDOF = NN*3;
            
            % Initialise variables
            this.Ed = zeros(NE,3*2); this.Q  = zeros(NE-1,3); F_appl  = zeros(NDOF,1);
            
            % Specify number of force increments up to full loading
            steps = 20;     

            for n=1:steps;
                % Increase force slightly every iteration but subtract the
                % previously applied force to prevent unbounded and
                % unphysical extension of the material.  (every new iteration 
                % proceeds from the previous state discarding the previous force 
                % and displacement information) 
                F = this.f * (n/steps)^3 - F_appl;
                F_appl = F + F_appl;        % The total applied force up to this point
                
                % Reset displacement
                this.Ed = 0*this.Ed;
                
                if n==1
                    % Determine and Assemble Stiffness Matrix
                    K=zeros(NDOF,NDOF);
                    for i = 1:NE
                        Ke = this.beam2g(this.Ex(i,:),this.Ey(i,:),this.Ep,this.N(i));
                        K  = this.assem(this.Edof(i,:),K,Ke);
                    end
                    
                    [a,~]   = this.solveq(K,F,this.Bc);
                    this.Ed = this.extract(this.Edof,a);
                end
                
                this.r = ones(NE*3+3,1);
                count = 0;
                

                while norm(this.r(7:end-3,1))>10^-7 || count<5 
                    count = count+1;
                    
                    % Determine internal loading 
                    this.g  = zeros(NE*3+3,1);
                    for i = 1:NE+1
                        G = this.InternalLoads(i);
                        this.g(i*3-2) = G(1);     
                        this.g(i*3-1) = G(2);
                        if ~isnan(G(3));  this.g(i*3) = G(3); end
                    end
                    
                    % Compute Global components of internal forces
                    this.N = this.Es(:,1);
                    
                    % Determine and Assemble Stiffness Matrix
                    K=zeros(NDOF,NDOF);
                    for i = 1:NE
                        Ke = this.beam2g(this.Ex(i,:),this.Ey(i,:),this.Ep,this.N(i));
                        K  = this.assem(this.Edof(i,:),K,Ke);
                    end
                    
                    % Calculate force residual (difference between internal
                    % and external forces)
                    this.r = F + this.g;
                    
                    % Adjust displacement for force residual
                    [a, ~] = this.solveq(K,this.r,this.Bc);
                    this.Ed = this.Ed + this.extract(this.Edof,a);
                    
                    % Interrupt simulation in case of divergence
                    if norm(this.r)>10^9
                        disp('FEM divergence')
                        return
                    end
                end
                
                % Update node locations
                this.Ex = this.Ex + 0.1*this.Ed(:,[1 4]);
                this.Ey = this.Ey + 0.1*this.Ed(:,[2 5]);
            end
            
            this.Ex = this.Ex + this.Ed(:,[1 4]);
            this.Ey = this.Ey + this.Ed(:,[2 5]);
            
            Geom1    = [[this.Ex1(:,1);this.Ex1(end,2)],[this.Ey1(:,1);this.Ey1(end,2)]];
            Lengths1 = ((Geom1(1:end-1,1)-Geom1(2:end,1)).^2+...
                        (Geom1(1:end-1,2)-Geom1(2:end,2)).^2).^0.5;
            
            Geom     = [[this.Ex(:,1);this.Ex(end,2)],[this.Ey(:,1);this.Ey(end,2)]];
            Lengths  = ((Geom(1:end-1,1)-Geom(2:end,1)).^2+...
                        (Geom(1:end-1,2)-Geom(2:end,2)).^2).^0.5;
                   
            %Determine internal loading for total extension from rest
            this.N = (Lengths-Lengths1) * this.Ep(1) * this.Ep(2) ./ Lengths1;
        end

        %  ---------------------- Catenary Specific -----------------------------------  %
        
        function CatenaryIn(this,Shape,F,Properties)
            % -------  Define Geometry -------- %
            % Coordinates
            NN      = length(Shape);
            NDOF    = NN*3;            % Number of Degrees of Freedom
            NE      = NN-1;           % Number of Elements
            this.Ex = [Shape(1:end-1,1),Shape(2:end,1)];
            this.Ey = [Shape(1:end-1,2),Shape(2:end,2)];
            this.Ex1 = this.Ex;
            this.Ey1 = this.Ey;
            
            % Allocate forces
            this.f             = zeros(NDOF,1);
            this.f(1:3:NDOF-2) = F(1:end,1);
            this.f(2:3:NDOF-1) = F(1:end,2);
            
            % Boundary Conditions
            this.Bc = [1 2 NDOF-2 NDOF-1
                       0 0 0      0     ]';
            
            this.Es = zeros(NE,6);

            % Elements, connectivity and node locations
            COL  = 1:3:NDOF-5;  % DOF in first column of Edof
            
            this.Edof = [1:1:NE
                         COL
                         COL+1
                         COL+2
                         COL+3
                         COL+4
                         COL+5 ]';
                     
            % Element properties
            A       = Properties(1);
            I       = 1/12*Properties(1)^3;
            E       = Properties(2);
            this.Ep = [E A I];
            
            % Initialise internal forces to prevent initial divergence
            this.N = ones(NE,1);
        end
        
        function CatenaryRun(this,Shape,F,Properties)
            this.CatenaryIn(Shape,F,Properties)
            this.NewtonRhapson
            %this.CatenaryOut(Results,layer,i,j)
        end
 
    end
    
end