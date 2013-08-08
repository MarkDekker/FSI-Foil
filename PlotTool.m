classdef PlotTool < handle
    properties
        fig_hndl        % Main Figure Handle
        x               % x-coordinates of nodes
        y               % y-coordinates of nodes
        Alpha           % Angles of Attack
        Cp              % Pressure Coefficients
        Cf              % Friction Coefficients
        Cl              % Lift Coefficiens (per angle of attack)
        Cd              % Drag Coefficients (per angle of attack)
        Defomation      % Last deformation of every node
        Batch           % Struct of results files to be compared directly 
        ncases = 1      % number of cases currently under consideration
        plotID          % individual plot handles (i,j) i - figure number, j - plot number
        
        % Ladson NACA 0012 validation data
        L0   = importdata('./ValidationData/Ladson0.mat'); 
        L80  = importdata('./ValidationData/Ladson80.mat');  
        L120 = importdata('./ValidationData/Ladson120.mat');
        L180 = importdata('./ValidationData/Ladson180.mat'); 
        
        % Gregory NACA 0012 validation data
        G0   = importdata('./ValidationData/Gregory0.mat');  
        G10  = importdata('./ValidationData/Gregory10.mat');
        G15  = importdata('./ValidationData/Gregory15.mat');
        
        % Is current plot a validation case? (default: false)
        validation = false 
        
        skip_alph       % Logical array that makes it easy to skip angles 
                        % of attack where the solution did not converge
    end %properies
    
    methods (Static)
        function hline(y,in1,in2)
            % function h=hline(y, linetype, label)
            %
            % Draws a horizontal line on the current axes at the location specified by 'y'.  Optional arguments are
            % 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
            % label appears in the same color as the line.
            %
            % The line is held on the current axes, and after plotting the line, the function returns the axes to
            % its prior hold state.
            %
            % The HandleVisibility property of the line object is set to "off", so not only does it not appear on
            % legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
            % return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be
            % overridden by setting the root's ShowHiddenHandles property to on.
            %
            % h = hline(42,'g','The Answer')
            %
            % returns a handle to a green horizontal line on the current axes at y=42, and creates a text object on
            % the current axes, close to the line, which reads "The Answer".
            %
            % hline also supports vector inputs to draw multiple lines at once.  For example,
            %
            % hline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
            %
            % draws three lines with the appropriate labels and colors.
            %
            % By Brandon Kuczenski for Kensington Labs.
            % brandon_kuczenski@kensingtonlabs.com
            % 8 November 2001
            
            if length(y)>1  % vector input
                for I=1:length(y)
                    switch nargin
                        case 1
                            linetype='k:';
                            label='';
                        case 2
                            if ~iscell(in1)
                                in1={in1};
                            end
                            if I>length(in1)
                                linetype=in1{end};
                            else
                                linetype=in1{I};
                            end
                            label='';
                        case 3
                            if ~iscell(in1)
                                in1={in1};
                            end
                            if ~iscell(in2)
                                in2={in2};
                            end
                            if I>length(in1)
                                linetype=in1{end};
                            else
                                linetype=in1{I};
                            end
                            if I>length(in2)
                                label=in2{end};
                            else
                                label=in2{I};
                            end
                    end
                    h(I)=hline(y(I),linetype,label);
                end
            else
                switch nargin
                    case 1
                        linetype='k:';
                        label='';
                    case 2
                        linetype=in1;
                        label='';
                    case 3
                        linetype=in1;
                        label=in2;
                end
                
                
                
                
                g=ishold(gca);
                hold on
                
                x=get(gca,'xlim');
                h=plot(x,[y y],linetype);
                if ~isempty(label)
                    yy=get(gca,'ylim');
                    yrange=yy(2)-yy(1);
                    yunit=(y-yy(1))/yrange;
                    if yunit<0.2
                        text(x(1)+0.02*(x(2)-x(1)),y+0.02*yrange,label,'color',get(h,'color'))
                    else
                        text(x(1)+0.02*(x(2)-x(1)),y-0.02*yrange,label,'color',get(h,'color'))
                    end
                end
                
                if g==0
                    hold off
                end
                set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends
            end % else
        end
        
        function vline(x,in1,in2)
            % function h=vline(x, linetype, label)
            %
            % Draws a vertical line on the current axes at the location specified by 'x'.  Optional arguments are
            % 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
            % label appears in the same color as the line.
            %
            % The line is held on the current axes, and after plotting the line, the function returns the axes to
            % its prior hold state.
            %
            % The HandleVisibility property of the line object is set to "off", so not only does it not appear on
            % legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
            % return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be
            % overridden by setting the root's ShowHiddenHandles property to on.
            %
            % h = vline(42,'g','The Answer')
            %
            % returns a handle to a green vertical line on the current axes at x=42, and creates a text object on
            % the current axes, close to the line, which reads "The Answer".
            %
            % vline also supports vector inputs to draw multiple lines at once.  For example,
            %
            % vline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
            %
            % draws three lines with the appropriate labels and colors.
            %
            % By Brandon Kuczenski for Kensington Labs.
            % brandon_kuczenski@kensingtonlabs.com
            % 8 November 2001
            
            if length(x)>1  % vector input
                for I=1:length(x)
                    switch nargin
                        case 1
                            linetype='k:';
                            label='';
                        case 2
                            if ~iscell(in1)
                                in1={in1};
                            end
                            if I>length(in1)
                                linetype=in1{end};
                            else
                                linetype=in1{I};
                            end
                            label='';
                        case 3
                            if ~iscell(in1)
                                in1={in1};
                            end
                            if ~iscell(in2)
                                in2={in2};
                            end
                            if I>length(in1)
                                linetype=in1{end};
                            else
                                linetype=in1{I};
                            end
                            if I>length(in2)
                                label=in2{end};
                            else
                                label=in2{I};
                            end
                    end
                    h(I)=vline(x(I),linetype,label);
                end
            else
                switch nargin
                    case 1
                        linetype='k:';
                        label='';
                    case 2
                        linetype=in1;
                        label='';
                    case 3
                        linetype=in1;
                        label=in2;
                end
                
                
                
                
                g=ishold(gca);
                hold on
                
                y=get(gca,'ylim');
                h=plot([x x],y,linetype);
                if length(label)
                    xx=get(gca,'xlim');
                    xrange=xx(2)-xx(1);
                    xunit=(x-xx(1))/xrange;
                    if xunit<0.8
                        text(x+0.01*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
                    else
                        text(x-.05*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
                    end
                end
                
                if g==0
                    hold off
                end
                set(h,'tag','vline','handlevisibility','off')
            end % else
            
        end
    end
          
    methods
         %--------------------------- class constructor ---------------------------------%

         function this = PlotTool
             %this.fig_hndl = figure(1);
         end %class constructor
         
         %-------------------------- Load data of relevant case -------------------------%
         
         function loadMultiple(this)
             % Instantiate Results.m to read objects
             Solution = Results;
             
             % Get names of all files contained in the batch folder
             fnames = dir([pwd,'/Results-Tests/Load/*.mat']);
             
             % From this, determine how many cases will be cmpared
             this.ncases = length(fnames);
             this.Batch  = cell(1,length(fnames));
             
             % Load all cases into the object 'Solution'
             for k=1:length(fnames)
                 fname = fnames(k).name;
                 this.Batch{k} = load([pwd,'/Results-Tests/Load/',fname],'Solution');
             end
         end % Load previous cases for direct comparison
             
         %---------------------- Define various plotting operations ---------------------%
         
         % Plot Global Properties
         
         function PlotCL(this,Solution,clr,figr,plotN)
             if nargin<6
                 figr = 1;
                 plotN = 1;
             end
             
             if nargin<5
                 clr = 0.7*[1, 1, 1];
             end
             
             Cases  = length(Solution.Iter);
             Alphas = zeros(1,Cases); CLs = zeros(1,Cases);
             % Get final converged values of Alpha and CL
             for i = 1:Cases
                j         = Solution.Iter(i)-1;
                Alphas(i) = Solution.Alphas(i,j);
                CLs(i)     = Solution.CL(i,j);
             end
             
             %Plot these results
             this.plotID(figr,plotN) = plot(Alphas,CLs,'.-','Color',clr);
             %this.plotID(figr,plotN+1) = plot(Alphas+1.2,CLs,'.-','Color',clr*0.5);
             xlabel('\alpha [deg]','FontSize',12)
             ylabel('C_l','rot',0,'Units','Normalized','Position',[-0.18, 0.5, 0],'FontSize',14)
             xlim([Alphas(1)-1 Alphas(end)+2])
         end 
         
         function PlotCD(this,Solution,clr,figr,plotN)
             if nargin<6
                 figr = 1;
                 plotN = 1;
             end
             
             if nargin<5
                 clr = 0.7*[1, 1, 1];
             end
             
             % Get final converged values of Alpha and CL
             for i = 1:length(Solution.Iter)
                j         = Solution.Iter(i)-1;
                Alphas(i) = Solution.Alphas(i,j);
                CDs(i)    = Solution.CD(i,j);
             end
             
             %Plot these results
             this.plotID(figr,plotN) = plot(Alphas,CDs,'.-','Color',clr);
              %this.plotID(figr,plotN+1) = plot(Alphas+1.2,CDs,'.-','Color',clr*0.5);
             xlabel('\alpha [deg]','FontSize',12)
             ylabel('C_d','rot',0,'Units','Normalized','Position',[-0.19, 0.5, 0],...
                 'FontSize',14)
             xlim([Alphas(1)-1 Alphas(end)+2])
         end
         
         function PlotLD(this,Solution,clr,figr,plotN)
             if nargin<6
                 figr = 1;
                 plotN = 1;
             end
             
             if nargin<5
                 clr = 0.7*[1, 1, 1];
             end
             
             % Get final converged values of Alpha and CL
             for i = 1:length(Solution.Iter)
                j         = Solution.Iter(i);
                Alphas(i) = Solution.Alphas(i,j);
                CDs(i)    = Solution.CD(i,j);
                CLs(i)    = Solution.CL(i,j); 
             end
             
             %Plot these results
             this.plotID(figr,plotN) = plot(Alphas,CLs./CDs,'Color',clr,'LineWidth',1.5);
             xlabel('\alpha [deg]')
             ylabel('L/D','rot',0)
             xlim([Alphas(1)-1 Alphas(end)+1])
             ylim([-70 140])
         end
         
         % Plot Local Properties
         
         function plotcp(this,Solution)
             init = find(round(Solution.Alphas)==0);
             plot(Solution.Airfoil{init,1}(:,1),-Solution.Cp{init,1}(:,1),...
                 'Color',0.7*[1, 1, 1],'LineWidth',1.5,'MarkerSize',12)
             title(['\alpha = ',num2str(Solution.Alphas(init)),'.'],'FontSize',16);
             xlabel('% chord','FontSize',16)
             ylabel('-C_p','FontSize',16)
             xlim([-0.01 1.01])
             
             % First get the figure's data-cursor mode, activate it, and set some of its properties
             cursorMode = datacursormode(gcf);
             set(cursorMode,'enable','on','UpdateFcn',{@setDataTipTxt,Solution,this},...
                 'NewDataCursorOnClick',false);
             
             function output_txt = setDataTipTxt(~,event_obj,Solution,this)
                 % Display the position of the data cursor
                 % obj          Currently not used (empty)
                 % event_obj    Handle to event object
                 % output_txt   Data cursor text string (string or cell array of strings).
                 
                 pos = get(event_obj,'Position');
                 output_txt = {['alpha: ',num2str(pos(1),4)],...
                     ['CL: ',num2str(pos(2),4)]};
                 
                 Case = find(Solution.Alphas==pos(1));
                 subplot(2,3,4);
                 cla reset;
                 hold on
                 plot(Solution.Airfoil{Case,1}(:,1),-Solution.Cp{Case,1}(:,1),...
                     'Color',0.7*[1, 1, 1],'LineWidth',1.5)
                 title(['\alpha = ',num2str(pos(1))...
                     ,'.'],'FontSize',16);
                 xlabel('% chord','FontSize',16)
                 ylabel('-C_p','FontSize',16)
                 xlim([-0.01 1.01])
                 
                 
                 % Check whether angle of attack is roughly in the region
                 % of the relevant validation data (and that this is a 
                 % validation case). If this is the case, plot it. 
                 appr = round(pos(1));
                 
                 if appr ==0 || appr == 10 || appr ==15
                     this.NACA0012validation_cp(appr)
                 end
                 
                 hold off
             end
         end % PlotCp
         
         function PlotCP(this,Solution,i,j,clr,figr,plotN)
             if nargin<6
                 figr = 1;
                 plotN = 1;
             end
             
             if nargin<5
                 clr = [0.5, 0.5, 0.5];
             end
             
             this.plotID(figr,plotN) = plot(Solution.AirfoilXF{i,j}(:,1),...
                                           -Solution.Cp{i,j}(:,1),'Color',clr);
             xlabel('% chord','FontSize',12)
             ylabel('-C_p','rot',0,'Units','Normalized','Position',[-0.14, 0.5, 0],...
                 'FontSize',14)
             xlim([-0.01 1.01])

         end % PlotCp
         
         % Geometry related Plotting
         
         function PlotProfile(this,Results,i,j,clr,figr,plotN,min_y)
             if nargin<7
                 figr = 1;
                 plotN = 1;
             end
             
             if nargin<5
                 clr = [0.5, 0.5, 0.5];
             end
             
             if nargin<8
                 min_y = -0.1;
             end
             
            offset = min_y + 0.1;
             
            % Calculate and plot profiles
            Profile = Results.Airfoil{i,j};
            this.plotID(figr,plotN) = plot(Profile(:,1),Profile(:,2)+offset,...
                'Color',clr,'LineWidth',0.5);
            
            % Format Profile Plot
            set(gca,'DataAspectRatio',[1 1 1])
            xlim([0 1])
            ylim([min_y 0.2])
            xlabel('% chord','fontsize',13)
            %ylabel('thickness (as % of Chord)','fontsize',11)
         end
         
         function PlotBaseProfile(this,Results,clr,figr,plotN)
             if nargin<6
                 figr = 1;
                 plotN = 1;
             end
             
             if nargin<5
                 clr = [0, 0, 0];
             end
             
            % Calculate and plot profiles
            Profile = Results.Airfoil{1,1};
            this.plotID(figr,plotN) = plot(Profile(:,1),Profile(:,2),...
                'Color',clr,'LineWidth',0.7,'Marker','.','MarkerSize',4);
            
            % Format Profile Plot
            set(gca,'DataAspectRatio',[1 1 1])
            xlim([0 1])
            ylim([-0.1 0.15])
            xlabel('% chord','fontsize',11,'Position',[0.5,-0.15,1])
            ylabel('thickness (as % of Chord)','fontsize',11)
         end
        
         function PlotMast(this,Results,figr,plotN,clr,min_y)
             if nargin<4
                 figr = 1;
                 plotN = 1;
             end
             
             if nargin<5
                 clr = [0.5, 0.5, 0.5];
             end
             
             if nargin<6
                 min_y = -0.1;
             end
             
             offset = min_y + 0.1;
             
             
             % Plot mast
             this.plotID(figr,plotN) = plot(Results.Mast{1,1}(1,:),Results.Mast{1,1}(2,:)+offset,...
                 'Color',clr,'LineWidth',0.5,'LineStyle','--');
             plot(Results.Mast{1,2}(1,:),Results.Mast{1,2}(2,:)+offset,...
                 'Color',clr,'LineWidth',0.5,'LineStyle','--');
             % Format Profile Plot
             set(gca,'DataAspectRatio',[1 1 1])
             xlim([0 1])
             ylim([-0.1 0.2])
             xlabel('% chord','fontsize',11)
             %ylabel('thickness (as % of Chord)','fontsize',11)
         end
         
         function PlotBL(this,Solution,i,j,clr,figr,plotN)
             if nargin<6
                 figr = 1;
                 plotN = 1;
             end
             
             if nargin<5
                 clr = [0.5, 0.5, 0.5];
             end
             
            % Calculate and plot profiles
            % Generate Boundary layer representation
            BL = [Solution.Theta{i,j}(1:length(Solution.Airfoil{i,j})) .* ...
                  Solution.H{i,j}(1:length(Solution.Airfoil{i,j})) .* ...
                  Solution.Normals{i,j}(:,1), ...
                  Solution.Theta{i,j}(1:length(Solution.Airfoil{i,j})) .* ...
                  Solution.H{i,j}(1:length(Solution.Airfoil{i,j})) .* ...
                  Solution.Normals{i,j}(:,2)]; 
            
            Profile = Solution.Airfoil{i,j} + BL;
            this.plotID(figr,plotN) = plot(Profile(:,1),Profile(:,2),...
                'Color',clr,'LineWidth',0.7,'LineStyle',':');
            
            % Format Profile Plot
            set(gca,'DataAspectRatio',[1 1 1])
            xlim([0 1])
            ylim([-0.1 0.2])
            xlabel('% chord','fontsize',14,'Position',[0.5,-0.15,1])
            ylabel('thickness (% chord)','fontsize',14)
         end
        
        
         %-------------------- Define functions to plot validation data -----------------%
         
         function NACA0012validation_cl(this)
             % Plot the NACA0012 validation Results + add Legend
             plot(this.L0(:,1),this.L0(:,2),'.','Color',[0 0 0],'MarkerSize',12);
         end 
         
         function NACA0012validation_cd(this)
             % Plot the NACA0012 validation Results + add Legend
             plot(this.L0(:,1),this.L0(:,3),'.','Color',[0 0 0],'MarkerSize',12);
         end 
         
         function NACA0012validation_LtoD(this)
             % Plot the NACA0012 validation Results + add Legend
             plot(this.L0(:,1),this.L0(:,2)./this.L0(:,3),'.','Color',[0 0 0],'MarkerSize',12);
         end 
         
         function NACA0012validation_cp(this,alpha)
             % Plot the NACA0012 validation Results + add Legend
             if alpha == 0
             plot(this.G0(:,1),-this.G0(:,2),'.','Color',[0 0 0],'MarkerSize',12);
             elseif alpha == 10
             plot(this.G10(:,1),-this.G10(:,2),'.','Color',[0 0 0],'MarkerSize',12);    
             elseif alpha == 15
             plot(this.G15(:,1),-this.G15(:,2),'.','Color',[0 0 0],'MarkerSize',12);
             end
         end
         
         function NACA0012_CL_error(this,Solution)
             %Plot the relative error of simulation results to validation
             %data
             CL_Lad0  = this.L0(:,2);
             CL_Lad120 = this.L120(:,2);
             CL_Lad180 = this.L180(:,2);
             
             CL_Error80  = -(CL_Lad0 - Solution.CL)*100;
             CL_Error120 = -(CL_Lad120 - Solution.CL)*100;
             CL_Error180 = -(CL_Lad180 - Solution.CL)*100;
             
             bar(Solution.Alphas,CL_Error80)
             ylabel('%Error (normalised w.r.t. C_L=1)')
             xlabel('\alpha [deg]')
             xlim([-5 20])
             ylim([-20 20])
         end
         
         function NACA0012_CD_error(this,Solution)
             %Plot the relative error of simulation results to validation
             %data
             CD_Lad0   = this.L0(:,3);
             CD_Lad120 = this.L120(:,3);
             CD_Lad180 = this.L180(:,3);
             
             CD_Error0  = -(CD_Lad0 - Solution.CD)./...
                            (CD_Lad0)*100;
             CD_Error120 = -(CD_Lad120 - Solution.CD)./...
                            (CD_Lad120)*100;
             CD_Error180 = -(CD_Lad180 - Solution.CD)./...
                            (CD_Lad180)*100;
             
             bar(Solution.Alphas, CD_Error0)
             ylabel('% Error (normalised)')
             xlabel('\alpha [deg]')
             xlim([-5 20])    
             ylim([-50 50])
         end
         
         function vdBProfiles(this,Alpha,min_y,max_y)
             
            if nargin<3
               min_y = -.1;
               max_y = 0.2; 
            end
             
            % replace with an image of your choice
            if round(Alpha) == 1
                    img = imread('./ValidationData/Alpha_1.png');
            elseif round(Alpha) == 3
                    img = imread('./ValidationData/Alpha_3.png'); 
            elseif round(Alpha) == 5
                    img = imread('./ValidationData/Alpha_5.png');
            elseif round(Alpha) >= 7
                    img = imread('./ValidationData/Alpha_7_up.png');
            end
            
            % set the range of the axes
            % The image will be stretched to this.
            min_x = 0;
            max_x = 1;
            
            % Flip the image upside down before showing it
            imagesc([min_x max_x], [min_y max_y], flipdim(img,1));
         end
         
         %------------ Define various visualisation layouts that can be called ----------%
         
         % NACA0012  
         
         function NACA0012Results(this,Solution)
             clf
             clr = 0.7*[1, 1, 1];
             % Plot the current case (hold statements overlay plots)
             subplot(2,3,1);
             hold on
             this.PlotCL(Solution,clr,1,1)
             this.NACA0012validation_cl
             legend('Simulation','Ladson Free Transition','Location','NorthWest')
             this.vline(0);
             this.hline(0);
             
             a(2) = subplot(2,3,2);
             hold on
             this.PlotCD(Solution,clr,2,1)
             this.NACA0012validation_cd
             legend('Simulation','Ladson Free Transition','Location','NorthWest')
             this.vline(0);
             
             a(3) = subplot(2,3,3);
             hold on
             this.PlotLD(Solution,clr,3,1)
             this.NACA0012validation_LtoD
             legend('Simulation','Ladson Free Transition','Location','SouthEast')
             this.vline(0);
             this.hline(0);
             
             a(4) = subplot(2,3,4);
             hold on
             this.plotcp(Solution)
             this.NACA0012validation_cp(0)
             hold off
             this.hline(0);
             
             a(5) = subplot(2,3,5);
             this.NACA0012_CD_error(Solution)
             legend('Ladson Free Transition','Location','SouthWest')
             
             a(6) = subplot(2,3,6);
             this.NACA0012_CL_error(Solution)
             legend('Ladson Free Transition','Location','SouthEast')
             
             hold off
             set(this.fig_hndl,'Units','normalized','Position',[0 0 1 1])
           
             %saveas(this.fig_hndl,path,'eps');
             %saveas(this.fig_hndl,path,'fig');
         end 
         
         function NACA0012_MeshConv(this)
             this.loadMultiple;
             
             figure(1)
             hold on
             for i = 1:this.ncases
                 clr = (1 - i*0.2)*[1 1 1];
                 this.PlotCL(this.Batch{1,i}.Solution,clr,1,i)
             end
             legend('35 Nodes','70 Nodes','140 Nodes',...
                 '280 Nodes','Location','NorthWest')
             this.hline(0);
             this.vline(0);
             hold off
             
             figure(2)
             hold on
             for i = 1:this.ncases
                 clr = (1 - i*0.2)*[1 1 1];
                 this.PlotCD(this.Batch{1,i}.Solution,clr,2,i)
             end
             legend('35 Nodes','70 Nodes','140 Nodes',...
                 '280 Nodes','Location','NorthWest')
             this.vline(0);
             hold off
             
             figure(3)
             hold on
             clr = 0.7*[1 1 1];
             this.PlotCL(this.Batch{1,3}.Solution,clr,3,1)
             this.NACA0012validation_cl
             legend('140 Nodes','Ladson Free Transition','Location','NorthWest')
             this.hline(0);
             this.vline(0);
             hold off
             
             figure(4)
             hold on
             clr = 0.7*[1 1 1];
             this.PlotCD(this.Batch{1,3}.Solution,clr,3,1)
             this.NACA0012validation_cd
             legend('140 Nodes','Ladson Free Transition','Location','NorthWest') 
             this.vline(0);
             hold off
         end
         
         function COMSOLConvergence(this)
             COMSOL_C  = importdata('./ValidationData/COMSOL_C.mat'); 
             COMSOL_BC = importdata('./ValidationData/COMSOL_B1.mat');
             COMSOL_M  = importdata('./ValidationData/COMSOL_M.mat');
             COMSOL_BF = importdata('./ValidationData/COMSOL_B2.mat');
             COMSOL_F  = importdata('./ValidationData/COMSOL_F.mat');
             
             COMSOL_Dom  = importdata('./ValidationData/COMSOL_Dom.mat');
             
             MrkrSz = 12;
             
             figure(1)
             hold on
             plot(COMSOL_C(:,1), COMSOL_C(:,3), '.','color',[1 1 1]*0.7,'MarkerSize',MrkrSz);
             plot(COMSOL_BC(:,1),COMSOL_BC(:,3),'o','color',[1 1 1]*0.5,'MarkerSize',MrkrSz/2);
             plot(COMSOL_C(:,1), COMSOL_M(:,3), '.','color',[1 1 1]*0.3,'MarkerSize',MrkrSz);
             plot(COMSOL_BC(:,1),COMSOL_BF(:,3),'o','color',[1 1 1]*0,  'MarkerSize',MrkrSz/2);
             plot(COMSOL_C(:,1), COMSOL_F(:,3), '.','color',[1 1 1]*0,  'MarkerSize',MrkrSz);
             this.hline(0)
             xlabel('\alpha [deg]','FontSize',12)
             ylabel('C_d','rot',0,'FontSize',12)
             legend({'coarse','intermediate coarse','medium','intermediate fine','fine'})
             hold off
              
             
             figure(2)
             hold on 
             plot(COMSOL_C(:,1), COMSOL_C(:,2), '.','color',[1 1 1]*0.7,'MarkerSize',MrkrSz);
             plot(COMSOL_BC(:,1),COMSOL_BC(:,2),'o','color',[1 1 1]*0.5,'MarkerSize',MrkrSz/2);
             plot(COMSOL_C(:,1), COMSOL_M(:,2), '.','color',[1 1 1]*0.3,'MarkerSize',MrkrSz);
             plot(COMSOL_BC(:,1),COMSOL_BF(:,2),'o','color',[1 1 1]*0,  'MarkerSize',MrkrSz/2);
             plot(COMSOL_C(:,1), COMSOL_F(:,2), '.','color',[1 1 1]*0,  'MarkerSize',MrkrSz);
             this.hline(0)
             xlabel('\alpha [deg]','FontSize',12)
             ylabel('C_l','rot',0,'FontSize',12)
             legend({'coarse','intermediate coarse','medium','intermediate fine','fine'})
             hold off
             
             for i=1:1:6
                 %Corresponding number of cells: C=20'000 B1=52'000 M=77'000 B2=204'000
                 %F= 302'000
                 %Corresponding number of DOF: C=52'900 B1=140'000 M=200'000 B2=520'000
                 %F= 770'000
                figure(2+i)
                
                subplot(1,2,1)
                cla
                plot(20000,COMSOL_C(i,2),  '.','color',[1 1 1]*0.7,'MarkerSize',MrkrSz);
                hold on
                plot(52000,COMSOL_BC(i,2), 'o','color',[1 1 1]*0.5,'MarkerSize',MrkrSz/2);
                plot(77000,COMSOL_M(i,2),  '.','color',[1 1 1]*0.3,'MarkerSize',MrkrSz);
                plot(204000,COMSOL_BF(i,2),'o','color',[1 1 1]*0,'MarkerSize',MrkrSz/2);
                plot(302000,COMSOL_F(i,2), '.','color',[1 1 1]*0,'MarkerSize',MrkrSz);
                xlabel('number of cells','FontSize',13)
                ylabel('C_l','rot',0,'FontSize',13)
                title(['solution convergence at \alpha =',num2str(COMSOL_C(i,1))],'FontSize',12)
                hold off
                
                subplot(1,2,2)
                cla
                plot(20000,COMSOL_C(i,3),  '.','color',[1 1 1]*0.7,'MarkerSize',MrkrSz);
                hold on
                plot(52000,COMSOL_BC(i,3), 'o','color',[1 1 1]*0.5,'MarkerSize',MrkrSz/2);
                plot(77000,COMSOL_M(i,3),  '.','color',[1 1 1]*0.3,'MarkerSize',MrkrSz);
                plot(204000,COMSOL_BF(i,3),'o','color',[1 1 1]*0,'MarkerSize',MrkrSz/2);
                plot(302000,COMSOL_F(i,3), '.','color',[1 1 1]*0,'MarkerSize',MrkrSz);
                xlabel('number of cells','FontSize',13)
                ylabel('C_d','rot',0,'FontSize',13)
                hold off
             end
             
             figure(9)
             plot(COMSOL_Dom(:,1),COMSOL_Dom(:,2),'.','color',[1 1 1]*0,'MarkerSize',MrkrSz);
             title('Re = 6''000''000 | \alpha = 10^o','FontSize',12)
             xlabel('domain size (normalised w.r.t. chord)','FontSize',12)
             ylabel('C_l','rot',0,'FontSize',12)
             
             figure(10)
             plot(COMSOL_Dom(:,1),COMSOL_Dom(:,3),'.','color',[1 1 1]*0,'MarkerSize',MrkrSz);
             title('Re = 6''000''000 | \alpha = 10^o','FontSize',12)
             xlabel('domain size (normalised w.r.t. chord)','FontSize',12)
             ylabel('C_d','rot',0,'FontSize',12)
             
         end
         
         function ClCdValidationComparison(this,Solution)
             COMSOL_M  = importdata('./ValidationData/COMSOL_B2.mat');

             subplot(1,2,1)
             hold on
             plot(COMSOL_M(:,1),COMSOL_M(:,2),'*','color',[1 1 1]*0.5);
             this.PlotCL(Solution,[1 1 1]*0.4,1,1)
             this.NACA0012validation_cl
             this.hline(0)
             this.vline(0)
             xlabel('\alpha [deg]','FontSize',16)
             ylabel('C_l','rot',0,'FontSize',16)
             hold off

             subplot(1,2,2)
             hold on 
             plot(COMSOL_M(:,1),COMSOL_M(:,3),'*','color',[1 1 1]*0.5);
             this.PlotCD(Solution,[1 1 1]*0.4,1,1)
             this.NACA0012validation_cd
             this.vline(0)             
             xlabel('\alpha [deg]','FontSize',16)
             ylabel('C_d','rot',0,'FontSize',16)
             legend({'COMSOL medium mesh','XFOIL 140 nodes','Ladson measurements'},'NorthEast')
             hold off
 
         end
         
         function NACACpCompare(this,Solution)
             clr = [0.7 * ones(1,3) 
                    0.5 * ones(1,3)
                    0.2 * ones(1,3)];
             
             for i = 2
                 figure(1)
                 if i==1; mesh = 'coarse'; elseif i == 2; mesh = 'medium'; elseif i == 3; 
                     mesh = 'fine'; end
                 
                 COMSOL_Cp = importdata(['./ValidationData/COMSOL_',mesh,'Cp_0.mat']);
                 hold on
                 plot(COMSOL_Cp(:,1),COMSOL_Cp(:,2),'Color',clr(i,:));
             end
             
             this.PlotCP(Solution,3,1,[0 0 0],1,5)
             this.NACA0012validation_cp(0)
             
             hold off
             title([mesh,' mesh at \alpha = 0']);
             xlabel('% chord','FontSize',14)
             ylabel('-C_p','rot',0,'FontSize',14)
             xlim([-0.01 1.01])
                 
             for i = 2
                 figure(2)
                 if i==1; mesh = 'coarse'; elseif i == 2; mesh = 'medium'; elseif i == 3; 
                     mesh = 'fine'; end
                 
                 COMSOL_Cp  = importdata(['./ValidationData/COMSOL_',mesh,'Cp_10.mat']);

                 hold on
                 plot(COMSOL_Cp(:,1),COMSOL_Cp(:,2),'Color',clr(i,:));
             end
             this.PlotCP(Solution,8,1,[0 0 0],1,5)
             this.NACA0012validation_cp(10)
             hold off
             title([mesh,' mesh at \alpha = 10']);
             xlabel('% chord','FontSize',14)
             ylabel('-C_p','rot',0,'FontSize',14)
             legend('COMSOL','XFOIL','Ladson experiments')
             xlim([-0.01 1.01])
         end
         
         % van den Borne

         function vdBProfileCompare(this,Solution,op_pt,fig)
             if nargin<3
                 op_pt = 1;
                 if nargin<4
                    fig = 1;
                 end
             end
             
             figure(fig)
             subplot(2,1,1) % Plot the series of airfoils to convergence
             hold on
             this.PlotMast(Solution,1,1,[0 0 0.7])
             this.PlotBaseProfile(Solution,[0 0 0.7],1,2)
             Cases = Solution.Iter(op_pt);
             for iter = 2:Cases
                 clr = [1, 1, 1]*(Cases-iter)/Cases;
                 this.PlotProfile(Solution,op_pt,iter,clr,1,iter+2)
             end
             this.hline(0);
             title(['\alpha = ',num2str(Solution.Alphas(op_pt,Cases)),'^{ o}'])
             hold off
             
             subplot(2,1,2)
             hold on
             this.vdBProfiles(Solution.Alphas(op_pt,1))
             this.PlotMast(Solution,2,1,[0 0 0.7])
             this.PlotProfile(Solution,op_pt,Cases,[1 1 1]*0,2,2)
             this.PlotProfile(Solution,op_pt,Cases,[1 1 1]*0.5,2,3)
             this.PlotProfile(Solution,op_pt,Cases,[0 0 0.7],2,4)
             %this.PlotBL(Solution,op_pt,Cases-1,[0.5 0.5 0.5],2,5)   
             
             h = this.plotID(2,1:5);
             legend(h,'mast profile','van den Borne measurement',...
                      'van den Borne simulation','FSIFoil simulation')%'FSIFoil boundary layer (\delta *)'
             this.hline(0);
             hold off
             
         end 
         
         function vdBForceCompare(this,Solution)
             vdB_Config1 = importdata('./ValidationData/vdB_Config1.mat');
             
             clr = [1 1 1]*0.7;
             
             subplot(1,2,1)
             hold on
             this.plotID(1,1) = plot(vdB_Config1(:,1),vdB_Config1(:,2),...
                                     '.-','Color',[0 0 0],'MarkerSize',6);
             this.PlotCL(Solution,clr,1,2)
             legend('van den Borne measurements','FSIFoil results (original)','FSIFoil results (offset)')
             hold off
             
             subplot(1,2,2)
             hold on
             this.plotID(2,1) = plot(vdB_Config1(:,1),vdB_Config1(:,3),...
                                     '.-','Color',[0 0 0],'MarkerSize',6);
             this.PlotCD(Solution,clr,2,2)
             %legend('van den Borne measurements','FSIFoil results')
             hold off
         end
         
         function ProfileOverlay(this,Solution)
             figure(1)
             for i = 1:5
                 Alfa = i*2-1; op_pt = find(round(Solution.Alphas(:,1))==Alfa);
                 Cases = Solution.Iter(op_pt);
                 min_y = -0.1 - (i-1)*0.3;
                 max_y =  0.2 - (i-1)*0.3;
                 hold on
                 this.vdBProfiles(Solution.Alphas(op_pt,1),min_y,max_y)
                 this.PlotMast(Solution,2,1,[0 0 1],min_y)
                 this.PlotProfile(Solution,op_pt,Cases,[0 0 1],2,4,min_y)
                 
                 this.hline(min_y+0.1);
                 hold off
             end
             set(gca,'YTick',[]) 
         end 
         
         % Sensitivity Analysis

         function ParameterVar(this,Names)
             this.loadMultiple;
             
             Fig1 = figure(1);
             hold on
             this.PlotMast(this.Batch{1,1}.Solution,2,1,[0 0 1])
             this.PlotMast(this.Batch{1,1}.Solution,2,1,[0 0 1],-0.3)

             for i = 1:this.ncases
                 Case1 = this.Batch{1,i}.Solution.Iter(2);
                 Case2 = this.Batch{1,i}.Solution.Iter(4);
                 clr = (1 - i/this.ncases)*[1 1 1];
                 
                 this.PlotProfile(this.Batch{1,i}.Solution,2,Case1,clr,1,i);
                 this.PlotProfile(this.Batch{1,i}.Solution,4,Case2,clr,1,i,-0.3);
                 if i == this.ncases; 
                     legend(this.plotID(1,:),Names,'Location','NorthEast');
                     this.hline(-0.2); this.hline(0);  set(gca,'YTick',[]);
                     set(gca,'Units','Normalized','Position',[0.01, 0.07, 0.98 0.98])
                 end 
                 
             end
             
             annotation('textbox','Units','Normalized','Position',[0.35 0.7 0.1 0.1],...
                        'FontName','Times New Roman','FontAngle','Italic','FontSize',16,...
                        'LineStyle','none','String','\alpha = 3^o')
                    
             annotation('textbox','Units','Normalized','Position',[0.35 0.33 0.1 0.1],...
                        'FontName','Times New Roman','FontAngle','Italic','FontSize',16,...
                        'LineStyle','none','String','\alpha = 7^o')
             hold off
             
             
             Fig2 = figure(2);
             for i = 1:this.ncases
                 Case1 = this.Batch{1,i}.Solution.Iter(2)-1;
                 Case2 = this.Batch{1,i}.Solution.Iter(5)-1;
                 clr = (1 - i/this.ncases)*[1 1 1];
                 subplot(1,2,1)
                 hold on
                 this.PlotCP(this.Batch{1,i}.Solution,2,Case1,clr,1,i);
                 
                 if i == this.ncases; 
                     this.hline(0);
                     set(gca,'Units','Normalized','Position',[0.09, 0.16, 0.39 0.78])
                 end
                 
                 subplot(1,2,2)
                 hold on
                 this.PlotCP(this.Batch{1,i}.Solution,5,Case2,clr,1,i);
                 if i == this.ncases; 
                     legend(Names,'Location','NorthEast'); 
                     this.hline(0);
                     set(gca,'Units','Normalized','Position',[0.6, 0.16, 0.39 0.78])
                 end 
                 
             end
             annotation('textbox','Units','Normalized','Position',[0.2 0.5 0.1 0.1],...
                        'FontName','Times New Roman','FontAngle','Italic','FontSize',16,...
                        'LineStyle','none','String','\alpha = 3^o')
                    
             annotation('textbox','Units','Normalized','Position',[0.67 0.44 0.1 0.1],...
                        'FontName','Times New Roman','FontAngle','Italic','FontSize',16,...
                        'LineStyle','none','String','\alpha = 9^o')
             hold off
            
             
             Fig3 = figure(3);
    
             for i = 1:this.ncases
                 clr = (1 - i/this.ncases)*[1 1 1];
                 
                 subplot(1,2,1)
                 hold on
                 this.PlotCL(this.Batch{1,i}.Solution,clr,2,i)
                 
                 if i == this.ncases; 
                     
                     set(gca,'Units','Normalized','Position',[0.09, 0.16, 0.39 0.78])
                 end

                 
                 subplot(1,2,2)
                 hold on
                 this.PlotCD(this.Batch{1,i}.Solution,clr,2,i)
                 
                 if i == this.ncases; 
                     legend(Names,'Location','NorthWest'); 
                     set(gca,'Units','Normalized','Position',[0.6, 0.16, 0.39 0.78])
                 end
                 
             end
             hold off
             % move figures to desired position and resize
             set(Fig1, 'Position', [10  400 600 330])
             set(Fig2, 'Position', [10  60  700 250])
             set(Fig3, 'Position', [750 60  700 250])
         end

    end %methods
    
end  %class