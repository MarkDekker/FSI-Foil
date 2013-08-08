 classdef XFOIL < handle
    properties
        Airfoil
        Actions    = {}
        Polars
        Dumps
        Pressures
        PolarFiles = {};
        PressureFiles = {};
        DumpFiles = {};
        Visible = true
        Process
        XFOILExecutable = 'xfoil.exe';
        KeepFiles = false;
        ID
    end
    
    properties (SetAccess = private)
        AirfoilFile = ''
    end
    
    
    methods (Static, Hidden)
        function ID = NewID()
            persistent LastID
            if isempty(LastID)
                LastID=0;
            end
            ID=LastID+1;
            LastID=ID;
        end
    end
    
    methods (Hidden)
        function CreateActionsFile(this)
            fid = fopen(this.ActionsFile,'wt+');
            if ~this.Visible
                fprintf(fid,'PLOP\nG\n\n');
            end
            fprintf(fid,'LOAD %s\n',this.AirfoilFile);
            fprintf(fid,'%s\n',this.Actions{:});
            fclose(fid);
        end       
    end
    
    methods (Static)
        function DownloadXFOIL
            h = waitbar(0,'Please wait, downloading XFOIL...');
            URL = 'http://web.mit.edu/drela/Public/web/xfoil/xfoil6.96.zip';
            [f, status] = urlwrite(URL, 'xfoil.zip');
            if status == 0
                warning('XFOIL:NotFound','XFOIL executable was not found in current MATLAB path, and I failed to download it. Please get it at http://web.mit.edu/drela/Public/web/xfoil/');
            end
            ClassPath = mfilename('fullpath');
            ClassPath = fileparts(fileparts(ClassPath)); %Get the directory above the class
            TargetPath = fullfile(ClassPath,'XFOIL');
            unzip(f,TargetPath);
            delete(f);
            addpath(genpath(TargetPath));
            delete(h);            
        end
 
        function AF = ActionsFile
            AF = 'Actions.txt';
        end
        
    end
    
    methods
        function this = XFOIL(XFOILExecutable)
            if nargin == 1 && ~isempty(XFOILExecutable)
                this.XFOILExecutable = XFOILExecutable;
            end
            this.ID=this.NewID;
            
            if this.ID == 1
                disp(repmat('-',1,70));
                disp(' XFOIL - MATLAB interface v1.0');
                disp(' Copyright (c) 2011 by Rafael Oliveira - Contact: rafael@rafael.aero')
                disp(repmat('-',1,70));
            end
            
            if isempty(which(this.XFOILExecutable))
                if ispc
                    ButtonName = questdlg('XFOIL executable not found, should I download it?', 'XFOIL','Yes','No','Yes');
                    if strcmp(ButtonName,'Yes')
                        XFOIL.DownloadXFOIL;
                    else
                        warning('XFOIL:NotFound','XFOIL executable was not found in current MATLAB path, Please get it at http://web.mit.edu/drela/Public/web/xfoil/');
                    end
                else
                    error('Unix version not yet implemented!')
                end
            end
        end
          
        function run(this)
            FoilFile = 'Current_Airfoil';
            this.AirfoilFile = [FoilFile '.dat'];

            if ~exist(this.AirfoilFile,'file')
                error('Airfoil file not found: %s', FoilFile)
            end
            
            this.CreateActionsFile;
            
            this.Actions = {};      % Make sure you clear actions for  
                                    % next calculation
                
            warning('off','MATLAB:DELETE:FileNotFound')
            for i=1:length(this.PolarFiles)
                delete(this.PolarFiles{1})
            end
            warning('on','MATLAB:DELETE:FileNotFound')
            
            if ispc
                if ~exist(fullfile(pwd,this.XFOILExecutable),'file')
                    xfEXE = which(this.XFOILExecutable);
                    [success,msg,~] = copyfile(xfEXE,fullfile(pwd,this.XFOILExecutable),'f');
                    if ~success
                        error('Error when copying XFOIL to current directory: %s',msg)
                    end
                end                
                arg = {'cmd', '/c',sprintf('"%s < %s"',this.XFOILExecutable,this.ActionsFile), '>','nul'};
            else
                error('Unix version not yet implemented!')
            end
            PB=java.lang.ProcessBuilder(arg);
            this.Process = PB.start;
            

        end
        
        function finished = wait(this,timeout)
            tStart = tic;
            finished=false;
            if nargin<2
                timeout=inf;
            end
            while toc(tStart)<timeout && finished==false
                try %#ok<TRYNC>
                    ev = this.Process.exitValue(); %#ok<NASGU>
                    this.Process=[];
                    finished=true;
                end
                pause(0.01)
            end
        end
        
        function kill (this)
            if ~isempty(this.Process)
                this.Process.destroy;
                this.Process =[];
            end
        end
        
        function addActions (this, NewActions)
            if ~iscell(NewActions)
                NewActions={NewActions};
            end
            
            this.Actions = cat(1,this.Actions, NewActions);
        end
                
        function addActionSet(this,last,Alpha,Re,Ma,N,...
                            Vacc,XTrTop,XTrBot,Iter) %Add all actions in a go

            % Assign default values if insuffcent input argument were given
            if nargin < 2
                last = true;
            end
            if nargin < 3
                Alpha = 0;
            end
            if nargin < 4
                Re    = 6000000;
            end
            if nargin < 5
                Ma    = 0;
            end
            if nargin < 6
                N     = 9;
            end
            if nargin < 7
                Vacc  = 0.01;
            end
            if nargin < 8
                XTrTop= 1;
            end
            if nargin < 9
                XTrBot=1;
            end
            if nargin < 10
                Iter  = 200;
            end
            
%             NewActions3 = {''; ''; ''; '';'' ;''; ''; ''; ...
%                 'PANE'; 'MDES'; 'FILT 1.00'; ''};
%             this.addActions (NewActions3);
            
            NewActions2 = {'OPER'; 'VPAR'; ...
                sprintf('N %1.2f',N); ...
                sprintf('VACC %1.4f',Vacc); ...
                'XTR'; sprintf('%1.4f',XTrTop); ...
                sprintf('%1.4f',XTrBot); ''; ...
                sprintf('VISC %1.4f',Re); ...
                sprintf('MACH %1.6f',Ma)};
            this.addActions (NewActions2);
            this.addActions(sprintf('ITER %i',Iter))
            
            PolarFile = 'Polar.txt';
            this.addActions({'PACC'; PolarFile; ''});
            this.PolarFiles{end+1} = PolarFile;

            this.addActions({sprintf('ALFA %2.4f',Alpha);'INIT'})
            
            this.addActions({sprintf('ALFA %2.4f',Alpha)})
            
            PressureFile = 'Pressure.txt';
            this.addActions(sprintf('CPWR %s', PressureFile));
            this.PressureFiles{end+1} = PressureFile;
            
            DumpFile = 'Dump.txt';
            this.addActions(sprintf('DUMP %s', DumpFile));
            this.DumpFiles{end+1} = DumpFile;
            
            this.addActions({'PACC';''})
            

            this.addActions({'';'';'';'';'';'';'';'QUIT';''});

        end  

    end    
 end