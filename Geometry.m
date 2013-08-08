 classdef Geometry < handle
    properties
        Length          % Non-normalised chord length
        Mast            % Mast Profile {top or bottom}(x or y,:)
        Mast_fun        % Mathematical Description of Mast Shape
        Sep_Point       % Point at which cloth separates from Mast {1}[val] - {layer}[ind]
        Shift_SP        % Shift in calculated top and bottom cloth separation points
        LayerLengths    % Top and Bottom Layer Lengths (in terms of normalised chord)          
        Layers          % Top and Bottom Layer Profiles {1}(1,:) - {layer}(x or y,:)
        Layers_fun      % Top and Bottom Layer Profiles
        Airfoil         % XFOIL compatible geometry
        TProfile        % Complete Sail Profile
        VdBMast         % Van der Borne Mast Profile
        VdBLayers       % Van der Borne Layer Lengths
        NoseAngle       % Angle of Mast w.r.t. Alpha (in rad)
        TECoord         % Trailing edge location
        physical        % Is the case possible? (Default - true)
        conform         % Normalising factor and rotation to obtain chord = 1 & alpha = 0              
        Tangents        % Tangents of discretised surface (based on Airfoil)
        Normals         % Normals of discretised surface (based on Airfoil)
        
    end %properies
    
    methods (Static)
        %------------------------- Discretisation techniques ----------------------------%
        function Discret = DiscrCos(y,x_min,x_max,Nodes)
            %***************************************%
            %   Cosine-based Shape Discretisation   %
            %***************************************%

            % Generate equi-spaced nodes in x
            x_equi  = x_min:(x_max-x_min)/(Nodes):x_max;
            
            %Set up function to stretch elements for optimum local resolution
            syms x_n d
            
            x_n  = 0.53*(1-cos(0.847071748*pi*d));
            %x_n  = 0.5* (1-cos(pi*d));
            
            %Create Coordinates Variable (All profiles assumed symmetrical)
            Discret(1,:) =  subs(x_n,'d',(x_equi-x_min)/(x_max-x_min))...
                                *(x_max-x_min)+x_min;
            Discret(2,:) =  subs(y,'x', Discret);
        end
        
        function Discret = DiscrHalfCos(y1,y2,x_min,x_mid,x_max,Nodes)
            %***************************************%
            %   Cosine-based Shape Discretisation   %
            %***************************************%
            
            %   NB!  requires symmetrical profile   %
            % Designed for partial functions y1(x) in front and y2(x) at the
            % rear
            
            % Generate equi-spaced nodes in x
            x_equi  = [x_min:2*(x_mid-x_min)/(Nodes):x_mid, ...
                       x_mid:2*(x_max-x_mid)/(Nodes):x_max];

            
            %Set up function to stretch elements for optimum local resolution
            syms x_n d
            
            x_n  = 0.5*(1-cos(pi*d));

            %Create Coordinates Variable (All profiles assumed symmetrical)
            Discret(1,:) =  subs(x_n,'d',(x_equi-x_min)/(x_max-x_min))...
                                *(x_max-x_min)+x_min;
            Discret(2,:) = [subs(y1,'x', Discret(1:Nodes/2+1)), ...
                            subs(y2,'x', Discret(Nodes/2+2:end))];
        end
        
        function Discret = DiscrSqrtCos(y,x_min,x_max,Nodes)
            %***************************************%
            %   Cosine-based Shape Discretisation   %
            %***************************************%
            
            %   NB!  requires symmetrical profile   %

            % Generate equi-spaced nodes in x
            x_equi  = x_min:(x_max-x_min)/(Nodes-1):x_max;
            
            %Set up function to stretch elements for optimum local resolution
            syms x_n d
            
            x_n  = cos(pi*(d/2.7-0.5))*1/cos(pi*(1/2.7-0.5));
            
            %Create Coordinates Variable (All profiles assumed symmetrical)
            Discret(1,:) =  subs(x_n,'d',(x_equi-x_min)/(x_max-x_min))...
                                *(x_max-x_min)+x_min;
            Discret(2,:) =  subs(y,'x',Discret);
        end
        
        function Discret = DiscrUniform(y,x_min,x_max,nnodes)
            %***************************************%
            %  Uniform Distribution Discretisation  %
            %***************************************%
            
            %Create Coordinates Variable (All profiles assumed symmetrical)
            Discret(1,:) =  x_min:(x_max-x_min)/(nnodes-1):x_max;
            Discret(2,:) =  subs(y,'x', Discret);
        end
    end
    
    methods

        %----------------------------- Class Constructor --------------------------------%

        function this = Geometry(nose_angle,chord,top_layer,...
                                    bot_layer,Shift_SP)
            %this.NoseAngle = nose_angle;
            this.Length    = chord;
            this.LayerLengths = [top_layer bot_layer];
            
            % Define location of trailing edge
            this.NoseAngle = nose_angle;
            this.TECoord   = [ cos(nose_angle/180*pi) 
                              -sin(nose_angle/180*pi)];
                          
            this.physical = true;
            
            % Record manual shift in cloth separation points
            this.Shift_SP = Shift_SP;
            
        end %class constructor

        %--------------------------- Create Custom Geometry -----------------------------%
         
        function Generate(this,MastType,Nodes,NACA,Scale)    
            % this.LayerLengths = [TopLength; BotLength];
            
            if nargin<5
                Scale = 1;
            end
            
            switch MastType
                case 'Elliptic'
                    this.Elliptic(NACA,Scale)       % NACA and Scale are actually a and b
                    %Generate top and bottom layer profiles with pretension
                    %technique.
                    
                    this.Pretension(1); %Top
                    this.Pretension(2); %Bottom
                    
                    this.Coarsen(Nodes);
                case 'Conformed NACA'
                    this.NACAMastGen(NACA,Scale)
                    %Generate top and bottom layer profiles with pretension
                    %technique.
                    
                    this.Pretension(1); %Top
                    this.Pretension(2); %Bottom
                    
                    this.Coarsen(Nodes);
                case 'Create other cases'
                    disp('You still need to define this technique.')
                otherwise
                    disp(['No technique defined to generate this mast',...
                        ' geometry.'])
            end
            
            
            
            % Prepare geometry for computation
            this.Conform;
            this.Export;
            this.TanNorm;
            
        end %Custom
        
        function Shape = Catenary(this,p1,p2,rad,Nodes)
            % Arc subtended from 
            x1 = p1(1); y1 = p1(2);
            x2 = p2(1); y2 = p2(2);
            
            % Choose initial radius
            r = rad;
            
            if r < sqrt((x1-x2)^2+(y1-y2)^2)/2
                disp('chosen radius is too small for endpoints')
            end
            
            % Find circle centrepoint
            xc = 0;         % Lies on the y-axis for simplicity
            yc = -sqrt(r^2-(xc-x1)^2)+y1;
            
            % Define semicircle function as a starting geometry for
            % catenary
            
            syms y x
            y = -sqrt(r^2-(xc-x)^2)+yc;

            %Discretise shape
            Shape = this.DiscrUniform(y,x1,x2,Nodes);
        end
        
        function NACA0012Gen(this,NACA,Nodes,Uni)
            % Decypher NACA 4-Digit Foil
            m = floor(NACA/1000)/100;
            p = floor((NACA-m*100000)/100)/10;
            t = (NACA-m*100000-p*1000)/100;
            
            % Generate Thickness Function
            syms y yt yc1 yc2 yvdB yU yL x xU xL theta1 theta2
            yt  = t/(0.2*1.0089305) * (0.2969*(x*1.0089305)^.5 - 0.126*(x*1.0089305) - ...
                0.3516*(x*1.0089305)^2 + 0.2843*(x*1.0089305)^3 - 0.1015*(x*1.0089305)^4);
            
            % Discretise
            if Uni == 1
                Discret = this.DiscrUniform(yt,0,1,Nodes/2);
            elseif Uni == 2
                Discret = this.DiscrCos(yt,0,1,Nodes/2);
            end
            
            % Rearrange coordinates (counterclockwise around airfoil from 
            % top of tail to bottom tail)
            this.Airfoil = [fliplr(Discret(1,:))  Discret(1,2:end)
                            fliplr(Discret(2,:)) -Discret(2,2:end)]';
            
            if p>0 || m>0
                disp('Not possible to generate cambered profiles')
            end
            
            this.TanNorm;
            this.Export;
            
        end
        
        %--------------------------- Create Mast Geometries -----------------------------%
        
        function Elliptic(this,a,b) 
            %***************************************%
            % Define Shape of Elipse Mathematically %
            %***************************************%
            syms y x
            y = sqrt(b^2-((x-b)^2)/1);

            roots = double(solve(y));
            x_min = min(roots); x_max = max(roots);
            
            Discret = this.DiscrCos(y,x_min,x_max,1000);
            this.Mast{1}(1,:) =  Discret(1,:);
            this.Mast{1}(2,:) =  Discret(2,:);  %Top Surface
            this.Mast{2}(1,:) =  Discret(1,:);
            this.Mast{2}(2,:) = -Discret(2,:);  %Bottom Surface
            
            this.Mast_fun = y;
            
        end %Elliptic Mast
        
        function NACAMastGen(this,NACA,Scale)
            % Decypher NACA 4-Digit Foil
            m = floor(NACA/1000)/100;
            p = floor((NACA-m*100000)/100)/10;
            t = (NACA-m*100000-p*1000)/100;
            
            % Generate Thickness Function
            syms y yt yc1 yc2 yvdB yU yL x xU xL theta1 theta2
            yt  = t/0.2 * (0.2969*x^.5 - 0.126*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4);
            
            %van der Borne rear polynomial
            xB   = 0.15;
            x1   = 0.07;
            
            y0   = subs(yt,x,xB);
            dy0  = subs(diff(yt,x),x,xB);
            ddy0 = subs(diff(yt,x,2),x,xB);
            
            a    = 0.07;
            b    = 1/(3*x1^2)*y0 + 1/(3*x1)*dy0 + 2/3*ddy0;
            c    = -2/x1*y0 - 3*dy0 - 2*x1*ddy0;
            d    = 8/(3*x1^.5)*y0 + 8/3*x1^.5*dy0 + 4/3*x1^1.5*ddy0;
            
            yvdB = b*(a-(x-xB))^2 + c*(a-(x-xB)) + d*(a-(x-xB))^.5;
            
            % Discretise uncambered mast
            x_min = 0;  x_mid = xB; x_max = x1+xB;
            Discret = this.DiscrHalfCos(yt,yvdB,x_min,x_mid,x_max,2000);
            
            % Add camber to profile (Unlikely that mast includes position of max camber)
            if m>0
                yc = m/p^2 * (2*p*x - x^2);
                theta = atan(diff(yc,x));
            else
                yc = 0*x;
                theta = 0*x;
            end

            xU = x - y*sin(theta);  yU = yc + y * cos(theta);
            xL = x + y*sin(theta);  yL = yc - y * cos(theta);

            this.Mast{1}(1,:) =  subs(xU,{x,y},{Discret(1,:),Discret(2,:)})*Scale;
            this.Mast{1}(2,:) =  subs(yU,{x,y},{Discret(1,:),Discret(2,:)})*Scale; %Top Surface
            this.Mast{2}(1,:) =  subs(xL,{x,y},{Discret(1,:),Discret(2,:)})*Scale;
            this.Mast{2}(2,:) =  subs(yL,{x,y},{Discret(1,:),Discret(2,:)})*Scale; %Bot Surface
            
        end
       
         %------------------------- Discretisation techniques ----------------------------%
        
        function Coarsen(this,NNodes) % Reduce number of nodes (Call before Conform)
            % Determine coarsening factor
            [~, len]        = size([this.Mast{1}(1,:) this.Layers{1}(1,:)...
                                    this.Mast{2}(1,:) this.Layers{2}(1,:)]);
            reduct          = floor(len*3/4/NNodes);
            
            %Coarsen Mast
            this.Mast{1}    = [this.Mast{1}(1:2,1:reduct:this.Sep_Point{1}) ...
                               this.Mast{1}(1:2,this.Sep_Point{1}:reduct:end) ...
                               this.Mast{1}(1:2,end)];
            this.Mast{2}    = [this.Mast{2}(1:2,1:reduct:this.Sep_Point{2}) ...
                               this.Mast{2}(1:2,this.Sep_Point{2}:reduct:end) ...
                               this.Mast{2}(1:2,end)];
            %Coarsen Layers
            this.Layers{1}  = [this.Layers{1}(:,1:reduct:end) this.Layers{1}(:,end)];
            this.Layers{2}  = [this.Layers{2}(:,1:reduct:end) this.Layers{2}(:,end)];
            
            %Recalculate where separation point is in new discretisation
            this.Sep_Point{1} =  ceil(this.Sep_Point{1}/reduct);     
            this.Sep_Point{2} =  ceil(this.Sep_Point{2}/reduct)+1; 
            
        end %Coarsen
            
        %--------------- Generate and Pretension Top and Bottom Layers ------------------%
        
        function Pretension(this,layer)
            % Read in trailing edge location
            b1      = this.TECoord(1);  b2      = this.TECoord(2);

            Cost(1) = 1;                          % Initialise Cost Vector
            
            % Determine Alebraic expression for membrane
            
            for c = 2:length(this.Mast{layer}(1,:))-1   
                a1      = this.Mast{layer}(1,c);
                a2      = this.Mast{layer}(2,c);
                
                att(c)  = sum(((this.Mast{layer}(1,1:c-1)-this.Mast{layer}(1,2:c)).^2 + ...
                    (this.Mast{layer}(2,1:c-1)- ...
                    this.Mast{layer}(2,2:c)).^2).^0.5);
                
                len(c)  = this.LayerLengths(layer) - att(c);
                
                dMast_dx = (this.Mast{layer}(2,c+1)-this.Mast{layer}(2,...
                    c-1))/(this.Mast{layer}(1,c+1) - ...
                    this.Mast{layer}(1,c-1));
                
                
                % Point between (a1,a2) and (b1,b2)
                m1      = a1 - (a1-b1)/2;
                m2      = a2 - (a2-b2)/2;
                
                % Line orthogonal to tangent at point (a1,a2)
%                 l1      = -1/dMast_dx * x + (a2 + 1/dMast_dx * a1);
                
                % Line orthogonal to A-B through point M
%                 l2      = -(a1-b1)/(a2-b2) * x + (m2 + (a1-b1)/(a2-b2) * m1);
                
                
                % Calculate the centre of the circle
                c1      = (m2-a2 + (a1-b1)/(a2-b2) * m1 - 1/dMast_dx * a1)/...
                    (-1/dMast_dx + (a1-b1)/(a2-b2));
                c2      = -1/dMast_dx * c1 + (a2 + 1/dMast_dx * a1);
                
                C(c,:)  = [c1 c2];
                
                % Calculate circle radius and distance between A and B
                r(c)       = sqrt((b1-c1)^2+(b2-c2)^2);
                d         = sqrt((a1-b1)^2+(a2-b2)^2);
                
                % Calculate arc length of current circle definition
                Arc_tem = r(c) * (2*asin((d/2)/r(c)));
                
                % Cost function to be minimised                
                Cost(c) = abs(Arc_tem - len(c));  
                
                if -1/dMast_dx > (a2-b2)/(a1-b1)
%                     Cost(c) = 1.0;
%                     continue
                elseif len(c)<d
                    Cost(c) = 1.1;
                    continue
                end
                
            end
            
            %Uncomment lines below to track convergence
            figure(layer)
            clf
            plot(Cost)
            
            figure
            plot(len+att)
            
            % Find point where fabric separates from Mast
            if layer == 2
                [~, index, ~] = find(Cost<0.0008,1,'last');
            elseif layer == 1
                [~, index, ~] = find(Cost<0.00015,1,'first');
            end
                
            if Cost(index)>0.0028
                this.physical = false;
                disp(['Layer number ',num2str(layer),' is not long enough',...
                    ' to satisfy the pretension conditions.'])
            end

            % Save mathematical description of profile
            syms x 
            this.Layers_fun{layer} = sqrt(r(index)^2 - (x-C(index,1))^2) ...
                                        + C(index,2);

            % Discretise profile with uniform distribution
             Nodes   = 4000;
             Discret = this.DiscrSqrtCos(this.Layers_fun{layer},...
                         this.Mast{layer}(1,index),b1,Nodes);
                    
            this.Layers{layer}(1,:) = Discret(1,:);
            this.Layers{layer}(2,:) = Discret(2,:);
            
            this.Sep_Point{layer}   = index;
            
        end % Pretension
        
        %--------------------- Rotate, Normalise and Export Airfoil ---------------------%
        
        function Conform(this)
            points = 30;

            % Read in trailing edge location
            b1      = this.TECoord(1);  b2      = this.TECoord(2);

            d = zeros(2,points/2); % Distance Cost Vector
            

            for i=1:points/2;
                a1      = this.Mast{1}(1,i);
                a2u     = this.Mast{1}(2,i);
                a2b     = this.Mast{2}(1,i);
                
                d(1,i)  = sqrt((a1-b1)^2+(a2u-b2)^2);
                d(2,i)  = sqrt((a1-b1)^2+(a2b-b2)^2);
                
            end
            
            [dis_t, ind_t] = max(d(1,:));
            [dis_b, ind_b] = max(d(2,:));
            
            if dis_t<dis_b
                ind       = ind_b;
                dist      = dis_b;
                shift     = [-this.Mast{1}(1,ind) -this.Mast{2}(1,ind)];
            else
                ind       = ind_t;
                dist      = dis_t;
                shift     = [-this.Mast{1}(1,ind) -this.Mast{1}(2,ind)];
            end
            
            %Shift Mast
            this.Mast{1}(1,:)   = this.Mast{1}(1,:)   + shift(1);
            this.Mast{2}(1,:)   = this.Mast{2}(1,:)   + shift(1);
            this.Mast{1}(2,:)   = this.Mast{1}(2,:)   + shift(2);
            this.Mast{2}(2,:)   = this.Mast{2}(2,:)   + shift(2);
            
            % Shift Layers
            this.Layers{1}(1,:) = this.Layers{1}(1,:) + shift(1);
            this.Layers{1}(2,:) = this.Layers{1}(2,:) + shift(2);
            this.Layers{2}(1,:) = this.Layers{2}(1,:) + shift(1);
            this.Layers{2}(2,:) = this.Layers{2}(2,:) + shift(2);

            % Rotate geometry
            theta = - asin(this.Layers{1}(2,end)/dist);
            
            transfert(1,:)      = this.Mast{1}(1,:) * cos(theta) - ...
                                  this.Mast{1}(2,:) * sin(theta);
            transfert(2,:)      = this.Mast{1}(1,:) * sin(theta) + ...
                                  this.Mast{1}(2,:) * cos(theta);
            transferb(1,:)      = this.Mast{2}(1,:) * sin(theta) + ...
                                  this.Mast{2}(2,:) * cos(theta);
            transferb(2,:)      = this.Mast{2}(1,:) * cos(theta) - ...
                                  this.Mast{2}(2,:) * sin(theta);
                              
            this.Mast{1}(1,:)     =  transfert(1,:);    
            this.Mast{1}(2,:)     =  transfert(2,:);  
            this.Mast{2}(2,:)     =  transferb(1,:);  
            this.Mast{2}(1,:)     =  transferb(2,:);  
            
            
            transfer2(1,:)     = this.Layers{1}(1,:)  * cos(theta) - ...
                                  this.Layers{1}(2,:)  * sin(theta);
            transfer3(1,:)     = this.Layers{2}(1,:)  * cos(theta) - ...
                                  this.Layers{2}(2,:)  * sin(theta);
            transfer2(2,:)     = this.Layers{1}(1,:)  * sin(theta) + ...
                                  this.Layers{1}(2,:)  * cos(theta);
            transfer3(2,:)     = this.Layers{2}(1,:)  * sin(theta) + ...
                                  this.Layers{2}(2,:)  * cos(theta);
                              
            this.Layers{1}(1,:) = transfer2(1,:);
            this.Layers{2}(1,:) = transfer3(1,:);
            this.Layers{1}(2,:) = transfer2(2,:);
            this.Layers{2}(2,:) = transfer3(2,:);
            
            % Normalise chord
            norm           = 1/dist;
            this.Mast{1}        = this.Mast{1}   * norm;
            this.Mast{2}        = this.Mast{2}   * norm;
            this.Layers{1}      = this.Layers{1} * norm;
            this.Layers{2}      = this.Layers{2} * norm;
            
            disp(' ');
            disp(repmat('-',1,70));
            disp(['The original airfoil with a chord length of ',...
                num2str(dist),' has been normalised to 1.']);
            disp(['It has also been rotated by ', num2str(theta/pi*180),...
                  ' degrees to an angle of attack of zero.'])
            disp(['This means that the relative top layer length is now ',...
                num2str(this.LayerLengths(1))])
            disp(['and the bottom layer length is now ', ...
                num2str(this.LayerLengths(2)),'.'])
            disp(repmat('-',1,70));
            disp(' ');
            
            % Record conforming parameters
            this.conform.len = norm;
            this.conform.ang = theta;  % Defined counterclockwise is positive
            
            % Rearrange coordinates (counterclockwise around airfoil from 
            % top of tail to bottom tail)
            this.Airfoil = [flipud(this.Layers{1,1}')
                            flipud(this.Mast{1}(:,1:this.Sep_Point{1})')
                            this.Mast{2}(:,2:this.Sep_Point{2}-1)'
                            this.Layers{1,2}']; 
            
            % Update separation point
            [~, ll1] = size(this.Layers{1}(1,:));
            [~, ll2] = size(this.Layers{2}(1,:));
            [la, ~] = size(this.Airfoil);
            
            this.Sep_Point{1} = ll1 + this.Shift_SP(1);
            this.Sep_Point{2} = la - ll2 + this.Shift_SP(2);
            
        end %Conform
        
        function Export(this)
            % Write XFoil export '.dat' file
            delete Current_Airfoil.dat
            fid          = fopen('Current_Airfoil.dat','w+');
            fprintf(fid,'%s\r\n','Current Airfoil');
            fprintf(fid,' %f %f \r\n',[this.Airfoil(:,1) this.Airfoil(:,2)]');
            fclose(fid);
        end % Export
        
        %----------------- Find tangents and normals of airfoil surface -----------------%
        
        function TanNorm(this)
            [len, ~]         = size(this.Airfoil);
            
            % Determine tangent vectors
            this.Tangents      = (this.Airfoil(3:end,2)-this.Airfoil(1:end-2,2))./...
                                 (this.Airfoil(3:end,1)-this.Airfoil(1:end-2,1));
            this.Tangents      = [this.Tangents(1,1)
                                  this.Tangents
                                  this.Tangents(end,1)];
            this.Tangents      = [ones(len,1) this.Tangents];
            norm_t             = (sum(this.Tangents.^2,2)).^-0.5;
            this.Tangents      = this.Tangents.*[norm_t norm_t];
            
            % Determine normal vectors
            this.Normals       = -(this.Airfoil(3:end,1)-this.Airfoil(1:end-2,1))./ ...
                                  (this.Airfoil(3:end,2)-this.Airfoil(1:end-2,2));
            this.Normals       = [this.Normals(1,1)
                                  this.Normals
                                  this.Normals(end,1)];
            this.Normals       = [ones(len,1) this.Normals];
            norm_n             = (sum(this.Normals.^2,2)).^-0.5;
            this.Normals       = this.Normals.*[norm_n norm_n];
            
            % Correct normal vectors to ensure they face outward
            nose               = find(this.Airfoil(:,1) == 0);
            err_t              = find(this.Normals(1:nose,2) < 0);
            err_b              = find(this.Normals(nose:end,2) > 0)+nose-1;
            this.Normals(err_t,:) = this.Normals(err_t,:)*-1;
            this.Normals(err_b,:) = this.Normals(err_b,:)*-1;
            
            % Uncomment to check that tangents and normals are correct:
%             figure
%             hold on
%                 plot(this.Airfoil(:,1),this.Airfoil(:,2))
%                 quiver(this.Airfoil(:,1),this.Airfoil(:,2),...
%                        this.Tangents(:,1),this.Tangents(:,2),0.1)
%                 quiver(this.Airfoil(:,1),this.Airfoil(:,2),...
%                        this.Normals(:,1),this.Normals(:,2),0.1)
%             hold off
%             set(gca,'DataAspectRatio',[1 1 1])
            
        end 
        
        %----------------- Reevaluate where Separation Point should be ------------------%
        
        function Separation(this,Solution,i,j)
            % Determine the expected separation point
            % Check the difference in the gradient before and after the
            % separation point
            SP = Solution.Sep_Point{i,j}{1};
            Grads = (Solution.Airfoil{i,j}(SP:SP+1,2) - Solution.Airfoil{i,j}(SP-1:SP,2)) ./ ...
                    (Solution.Airfoil{i,j}(SP:SP+1,1) - Solution.Airfoil{i,j}(SP-1:SP,1));
            
            % Shift the separation point gradually based on the magnitude
            % of the difference in gradients.  (Also account for stagnation
            % pressure where kink would be)
            
            
            if     (Grads(1,1)-Grads(2,1)) > 0.15
                Solution.Sep_Point{i,j}{1} = SP+10;
                Solution.F{i,j}(SP,:) = Solution.F{i,j}(SP-1,:);
            elseif (Grads(1,1)-Grads(2,1)) > 0.09
                Solution.Sep_Point{i,j}{1} = SP+2;
                Solution.F{i,j}(SP,:) = Solution.F{i,j}(SP-1,:);
            elseif (Grads(1,1)-Grads(2,1)) > 0.03
                Solution.Sep_Point{i,j}{1} = SP+1;
                Solution.F{i,j}(SP,:) = Solution.F{i,j}(SP-1,:);
                
            elseif (Grads(2,1)-Grads(1,1)) > 0.15
                Solution.Sep_Point{i,j}{1} = SP-5;
            elseif (Grads(2,1)-Grads(1,1)) > 0.07
                Solution.Sep_Point{i,j}{1} = SP-1;
            elseif (Grads(2,1)-Grads(1,1)) > 0.05
                Solution.Sep_Point{i,j}{1} = SP;
            end
            
            this.Sep_Point{1} = Solution.Sep_Point{i,j}{1};
        end
    end
end