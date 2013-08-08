classdef Coupling < handle
    methods
        %-------------------------- Class Constructor and Saving ------------------------%
        
        function this = Coupling
            
        end  
    end
    
    methods(Static)
        
        %----------------------------- Read Module Outputs ------------------------------%
        
        function updateAirfoil(Results,Geom,i,j)
            Results.Airfoil{i,j+1} = Results.Airfoil{1,1} + Results.U{i,j};
            Geom.Airfoil = Results.Airfoil{i,j+1};
            Geom.TanNorm
            Results.Normals{i,j+1}      = Geom.Normals;
            Results.Tangents{i,j+1}     = Geom.Tangents;
            Results.Sep_Point{i,j+1}    = Geom.Sep_Point;
            Results.N{i,j+1}            = zeros(length(Results.Airfoil{i,j+1}),1);
            Geom.Export
        end
        
        function newAlpha(Results,Geom,i,j)
            Results.Airfoil{i+1,1} = Results.Airfoil{1,1};
            
            Geom.Airfoil = Results.Airfoil{i+1,1};
            Geom.TanNorm
            Results.Normals{i+1,1}      = Geom.Normals;
            Results.Tangents{i+1,1}     = Geom.Tangents;
            Results.Sep_Point{i+1,1}    = Geom.Sep_Point;
            Results.N{i+1,1}            = zeros(length(Results.Airfoil{i,j}),1);
            Geom.Export
        end
        
        %---------------------- Infer Forces Acting on each Node ------------------------%
        
        function genForces(Results,i,j)
            % Determine Area for each corresponding node. (The first and
            % last are excluded because they are fixed in the model.)
            Areas_temp = ((Results.Airfoil{i,j}(2:end,1) - ...
                           Results.Airfoil{i,j}(1:end-1,1)).^2 + ...
                          (Results.Airfoil{i,j}(2:end,2) - ...
                           Results.Airfoil{i,j}(1:end-1,2)).^2).^0.5;
                            
            Results.Areas{i,j}(:,1) = 1/2*(Areas_temp(1:end-1,1) + Areas_temp(2:end,1)) ...
                                        * 0.36;
            
            if length(Results.V)==1
                c = 1;
            else
                c = i;
            end
            
            % Convert Pressure and Friction Coefficients into Forces
            [len, ~]        = size(Results.Cp{i,j}(2:end-1));
            Results.rho     = 1.225;%kg/m^3
            Tangen_F        =   Results.Cf{i,j}(2:end-1) .* Results.Areas{i,j} ...
                                * 0.5 * Results.rho * Results.V(c)^2 * Results.chord;
            Normal_F        = -(Results.Cp{i,j}(2:end-1) - Results.Cp_int(i)*ones(len,1)) ...
                               .* Results.Areas{i,j} * 0.5 * Results.rho * Results.V(c)^2 ...
                                * Results.chord;
                                
            % Write Tangent and Normal Forces in terms of Cartesian
            % Coordinate system
            Results.F{i,j}(:,1)    = Normal_F .* Results.Normals{i,j}(2:end-1,1) + ...
                                        Tangen_F .* Results.Tangents{i,j}(2:end-1,1);
            Results.F{i,j}(:,2)    = Normal_F .* Results.Normals{i,j}(2:end-1,2) + ...
                                        Tangen_F .* Results.Tangents{i,j}(2:end-1,2);
                                    
        end  
        
        %------------------------- Determine Length of Layers ---------------------------%
        
        function MeasureLyrs(Results,i,j)
            [~, zero] = min(Results.Airfoil{i,j}(:,1));
            
            Results.Length{i}(j,1) = sum(((Results.Airfoil{i,j}(2:zero,1) - ...
                                           Results.Airfoil{i,j}(1:zero-1,1)).^2 + ...
                                          (Results.Airfoil{i,j}(2:zero,2) - ...
                                           Results.Airfoil{i,j}(1:zero-1,2)).^2).^0.5);
            Results.Length{i}(j,2) = sum(((Results.Airfoil{i,j}(zero+1:end,1) - ...
                                           Results.Airfoil{i,j}(zero:end-1,1)).^2 + ...
                                          (Results.Airfoil{i,j}(zero+1:end,2) - ...
                                           Results.Airfoil{i,j}(zero:end-1,2)).^2).^0.5);
            
        end
        
        %------------- Relax displacements calculated by structural solver --------------%
        
        function relaxU(Results,i,j)
            if j<3
                rel = 1;
                Results.w{i,j}{1} = rel; Results.w{i,j}{2} = rel; Results.w{i,j}{3} = rel;
                Results.U{i,j} = rel*Results.U_fem{i,j};
            else
                for x = 1:2
                    r_i1 = Results.U_fem{i,j-1}(:,x) - Results.U{i,j-2}(:,x);
                    r_i2 = Results.U_fem{i,j}(:,x) - Results.U{i,j-1}(:,x);
                    
                    Results.w{i,j}{x}  = - Results.w{i,j-1}{x} * ...
                        (r_i1'*(r_i2-r_i1))/((r_i2-r_i1)'*(r_i2-r_i1));
                    
                    Results.U{i,j}(:,x)  = Results.U{i,j-1}(:,x) + Results.w{i,j}{x}*r_i2;
                end
                
            end
        end
        
        function reEval(Geom,Results,i,j) % Also update airfoil in here
            Results.U_fem{i,j}    = Results.U_spare*0.9;
            if j == 1
                Results.Airfoil{i,j+1} = Results.Airfoil{1,1};
            elseif j<3 && j>1 
                rel = 0.2*j^2;
                Results.w{i,j}{1} = rel; Results.w{i,j}{2} = rel; Results.w{i,j}{3} = rel;
                Results.U{i,j}    = rel*Results.U_spare;
                
                Results.Airfoil{i,j+1} = Results.Airfoil{i,1} + Results.U{i,j}(:,1:2);
            else 
                for x = 1:2
                    rel = Results.w{i,j-1}{x} * random('Normal',0.9,0.15);
                    Results.w{i,j}{x} = rel;
                    Results.U{i,j}(:,x)  = Results.U_spare(:,x) * Results.w{i,j}{x};
                end
                
                Results.Airfoil{i,j+1} = Results.Airfoil{i,1} + Results.U{i,j}(:,1:2);
            end
            Geom.Airfoil = Results.Airfoil{i,j+1};
            Geom.TanNorm
            Results.Normals{i,j+1}      = Geom.Normals;
            Results.Tangents{i,j+1}     = Geom.Tangents;
            Results.Sep_Point{i,j+1}    = Geom.Sep_Point;
            Results.N{i,j+1}            = zeros(length(Results.Airfoil{i,j+1}),1);
            Geom.Export
        end
        
        %----------------------------- Measure Convergence ------------------------------%
        
        function convCheck(Results,i,j)
            % Determine maximum and mean deviations between consecutive geometries.
            Dist = ((Results.Airfoil{i,j}(:,1)-Results.Airfoil{i,j+1}(:,1)).^2 + ...
                    (Results.Airfoil{i,j}(:,2)-Results.Airfoil{i,j+1}(:,2)).^2).^0.5;
            Results.maxDev(i,j)     = max(Dist);
            Results.meanDev(i,j)    = mean(Dist);

        end
        
    end
end