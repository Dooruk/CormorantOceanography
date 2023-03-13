clear;close all;
addpath('~/Documents/MATLAB/export_fig');
addpath('~/Documents/MATLAB/customcolormap');
addpath('~/Documents/MATLAB/rgb');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%*********** k is automatically converted to the -2 power (to linearize)*********

% This code combines assimilation variables from different sources and time periods.

% Requires some pre-processing of the variables and includes a lot of
% switches impacting the outputs and the plots.

% Automatically clears out the blewn members using Ens.memb variable that OUGHT TO be created
% during the processing (reading/interpolating) ensemble outputs stage (see below).

% A visual ensemble variable vs. ensemble parameter (e.g. u_surf vs. bathymetry) check is beneficial.

% Eliminates Hs observations inside the mouth. (east of longitude -124.05)

% Error levels are defined using the difference between observations
% and Prior

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
swch.Lc            =    1;                         % Omega_Vec correction (Localization on/off)
swch.prior         =    0;                         % use prior for the innovation [0=use ensemble average]
swch.twin          =    0;                         % Calculate Cov. for Hh [1] i.e. obs. points or h [0] i.e. the whole domain
swch.trmpoly       =    1;                         % Ignore points outside the polygon (Arbitrarily chosen)
swch.doprint       =    0;                         % print figures
swch.onscreen      =    1;                         % Display figure on screen or not
swch.filter        =    0;                         % Tidal filter (made things worse, avoid)
swch.filter2       =    0;                         % Adjoint Sensitivity filter (made things worse, avoid)
swch.filter3       =    1;                         % Std.Dev filter (made things better, embrace)
Trm.option         =    1;                         %[Trimming options: (1) only coord., (2) time and coord.,
%(3) multiple time and coord., (4) time and trim outside localization
%any option other than 1 requires unique obs. loop
swch.dostat        =   0;
% swch.multarea      =    0;                       % only relevant for dostat, divides the stat area into 4 regions
UVmax              =    7;                         % (m/s) absolute max. of the vel. variations (not active right now impacts folder names only)
swch.makeprior     =    0;
swch.savemat       =    0;                         % used for outputting .mat files for Kitchen_sink_plot_comb
Names.pert='pert2';                                % pert5 or pert2 (sigma_h (2m or 5m) for the ensemble members)

% Plot constants
Plt.ms      = 0.55;                  % Marker size
Plt.ms2     = 3;
Plt.covfac  = 12;
Plt.fs3     = 9;                   % Figure 3 Title font
Plt.fs      = 20;                   % axes titles
Plt.skl     = 0.6;                  % Only relevant for bird obs. Cuts down the low skill measurements (it was 0.75)
Plt.ulim    = 3;
Plt.res     = '-r150';            % Resolution of the plots
Plt.type    = '-png';            % Resolution of the plots
Plt.psize   = 20;                %pointsize for markers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Log:
%    - Variable names change: us,vs - surface velocities (u_surf)
%                             ub,vb - depth avg. velocities (u_bar)
%                             hs,l,p - sign. wave height, wave length, peak period
%                             dp    - peak wave direction

%    - To simplify combining different times, started writing file numbers in the ensemble variable
%    called Ens.memb.

%    - Check Ensemble members during combination


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rval is defined for each variable (u,v,hs and k) 4 times at differing levels. And then depending
% on the size of the Ensemble variables they are concatenated during the combination of the
% variables stage. So if you have 300 u and 100 k observations you end up with [300 R_u; 100 R_k]

% How to use this code (on top of the switch options above):

% 0.a) Eliminate observations outside the temporal and spatial scope of your project. Ideally you
% should have an idea (through Adjoint Sensitivities perhaps?) on when to use what type of
% observations but that's not imperative. Then BEFORE you start using this code, interpolate your
% Truth, Prior, Ensemble variables onto observation times and locations (Model * Interpolation_Matrix).
% Save them in .mat files and bring to a folder located in {$PWD}.

% 0.b) Make sure your interpolated .mat files are organized well, variable names are consistent, observation
% variables are named consistently (this will require some work but it will be worth it while using
% this code). At this time, I manually edited these files and put them into corresponding folders. Specify your
% original ensemble member sample size (N). All observation structures should be named 'Obs'.

N=300;
Time.basedate=datenum(2013,5,1,0,0,0);                 %model start time
Names.grid='../../romsprep/mcr_200m_2013_grd_v1.nc';
Names.polygon='Polygon_Pts6.mat';

% 1) Enter filenames manually. I originally wrote automated combinationS but manual typing will prevent
% file ordering issues from occuring.

% Truth and Prior files have 5 minute timesteps, Ensemble files have 15 minutes time steps so they
% are in different file structures. (Ofcourse, those times are application dependent.)

%n2 and n3 in the beginning of the bird files

birdavg=[];
% birdavg='n2/';

for stdvar=[6]
    for lll=[2] %localization loop
% 
        % Ens structure % ORDER IS IMPORTANT
        FileListEns = {...

                          '2013/SARJun3_ens_pert2.mat';...
                          '2013/SARJun3_13_ens_pert2.mat';...
%                         '2013/drifter_ens_pert2.mat';...
%                         ['AvgBirds2014/' birdavg 'Bird_longtide_pert2_15min_min_v_sur250m900s_time1.mat'];...
%                         ['AvgBirds2014/' birdavg 'Bird_longtide_pert2_15min_min_v_sur250m900s_time2.mat'];...
%                         ['AvgBirds2014/' birdavg 'Bird_longtide_pert2_15min_min_v_sur250m900s_time3.mat'];...
%                                                 '2019/2019_sound_ens_pert2.mat';...


            %real2:
            %iter1
%                         'real2/SARJun3_ens_real2iter1_pert2.mat';...
%                         'real2/SARJun3_13_ens_real2iter1_pert2.mat';...

%                                                 'real2/Bird_pert2_real2iter1_15min_250m900s_time1.mat';...
%                                                 'real2/Bird_pert2_real2iter1_15min_250m900s_time2.mat';...
%                                     'real2/Bird_pert2_real2iter1_15min_250m900s_time3.mat';...
%                                     'real2/2013_drifter_ens_real2iter1pert2.mat';...
            %                                             'real2/2019_sound_ens_pert2_30sec.mat';

            %            %iter2 (actual iter3)
%                         'real2/iter2/SARJun3_ens_real2iter2_pert2.mat';...
%                         'real2/iter2/SARJun3_13_ens_real2iter2_pert2.mat';...
%                         'real2/iter2/Bird_pert2_real2iter2_15min_250m900s_time1.mat';...
%                         'real2/iter2/Bird_pert2_real2iter2_15min_250m900s_time2.mat';...
%                         'real2/iter2/Bird_pert2_real2iter2_15min_250m900s_time3.mat';...
%                         'real2/iter2/2013_drifter_ens_real2iter2pert2.mat';...

            %            %iter3



%             'real2/iter3/2019_sound_ens_real2iter3_pert2_15sec_2obs.mat';...

};

        % True, Prior and Observation structures
        FileListTP = {...
                        '2013/SARJun3_trupri_pert2.mat';...
                        '2013/SARJun3_13_trupri_pert2.mat';...
%                         '2013/drifter_trupri_pert2.mat';...
%                         ['AvgBirds2014/' birdavg 'TrupriBird_comb_pert2_5min_v_sur250m900s_time1.mat'];...
%                         ['AvgBirds2014/' birdavg 'TrupriBird_comb_pert2_5min_v_sur250m900s_time2.mat'];...
%                         ['AvgBirds2014/' birdavg 'TrupriBird_comb_pert2_5min_v_sur250m900s_time3.mat'];...
%                                                 '2019/2019_sound_trupri_pert2.mat';...
            %real2:
            %iter1
%                         'real2/SARJun3_trupri_real2iter1_pert2.mat';...
%                         'real2/SARJun3_13_trupri_real2iter1_pert2.mat';...
%             
%                                     'real2/TrupriBird_pert2_real2iter1_5min_250m900s_time1.mat';...
%                         'real2/TrupriBird_pert2_real2iter1_5min_250m900s_time2.mat';...
%                         'real2/TrupriBird_pert2_real2iter1_5min_250m900s_time3.mat';...
%                         'real2/2013_drifter_trupri_real2iter1pert2.mat';...

            %                                 'real2/2019_sound_trupri_pert2_30sec.mat';...

            %iter2
%                         'real2/iter2/SARJun3_trupri_real2iter2_pert2.mat';...
%                         'real2/iter2/SARJun3_13_trupri_real2iter2_pert2.mat';...
%             
%                           'real2/iter2/TrupriBird_pert2_real2iter2_5min_250m900s_time1.mat';...
%                           'real2/iter2/TrupriBird_pert2_real2iter2_5min_250m900s_time2.mat';...
%                          'real2/iter2/TrupriBird_pert2_real2iter2_5min_250m900s_time3.mat';...
%                         'real2/iter2/2013_drifter_trupri_real2iter2pert2.mat';...
            %iter3
%             'real2/iter3/SARJun3_trupri_real2iter3_pert2.mat';...
%             'real2/iter3/SARJun3_13_trupri_real2iter3_pert2.mat';...
%             %             %
%             'real2/iter3/TrupriBird_pert2_real2iter3_5min_250m900s_time1.mat';...
%             'real2/iter3/TrupriBird_pert2_real2iter3_5min_250m900s_time2.mat';...
%             'real2/iter3/TrupriBird_pert2_real2iter3_5min_250m900s_time3.mat';...
%             'real2/iter3/2013_drifter_trupri_real2iter3pert2.mat';...

%             'real2/iter3/2019_sound_trupri_real2iter3_pert2_15sec_2obs.mat';...
           };
        %
        % 2) Decide what variables you would like to combine (refer to the variable list below). The code
        % will look for these fields

        %                             fields={'us_rho','vs_rho','hh_rho','k_rho'};
%         fields={'us_rho','vs_rho','hs_rho','k_rho'};
%         fields={'us_rho','vs_rho','hs_rho','k_rho','hh_rho'};
%                         fields={'us_rho','vs_rho'};
%                         fields={'us_rho','vs_rho','hs_rho','k_rho'};
%                                         fields={'hh_rho'};
                                        fields={'k_rho'};
%                                         fields={'hs_rho'};

        for ff=1:length(fields); Ens.(fields{ff})=[];end
        for ff=1:length(fields); True.(fields{ff})=[];end
        for ff=1:length(fields); Prior.(fields{ff})=[];end
        for ff=1:length(fields); Obs.(fields{ff})=[] ;end

        Rval1=[];Rval2=[];Rval3=[];Rval4=[];
        Rtitle={'Highest R','Higher but Realistic R','Lower but Realistic R','Lowest R'};
        close all
        if(swch.onscreen)
            %Set the subplot figure handles
            f3=figure(3);clf
            set(gcf,'units','centimeters','position',[5 5 40 40]); %paperposition is different
            [hf3,~]=tight_subplot(3,2,0.05,[0.05 0.1],[0.05 0.05]);
            fig3title = sgtitle('some initial title');

            f2=figure(2);clf
            set(gcf,'units','centimeters','position',[5 5 40 40]);
            [hf2,~]=tight_subplot(2,2,0.05,[0.05 0.05],[0.05 0.05]);
        end
        % ORIGINAL

        %v5
        Rvallist=[0.65,0.5,0.35,0.25;...             % Error options for us_rho
            0.4,0.3,0.2,0.15;...             % Error options for vs_rho
            0.75,0.5,0.1,0.05;...             % Error options for hs_rho
            80,60,40,20;...            % Error options for k (because k^2)
            %     0.025,0.01,0.005,0.0025;...            % Error options for k
%             3,2,1,0.5];           % Error options for hh_rho (v5)
    2.5,1.5,0.5,0.25];           % Error options for hh_rho %v6 (realistic is 2.5)

        % v4, these work well
        %                 Rvallist=[1,0.75,0.25,0.15;...             % Error options for us_rho
        %                     1.5,1,0.5,0.25;...             % Error options for vs_rho
        %                     0.5,0.5,0.1,0.05;...             % Error options for hs_rho
        %                     200,150,25,10;...            % Error options for k (because k^2)
        %                     %     0.025,0.01,0.005,0.0025;...            % Error options for k
        %                     1,0.5,0.25,0.1];           % Error options for hh_rho

        % Rval=[0.025,0.01,0.005,0.001;...             % Error options for u
        % Rval=[1;0.005;0.001;0.0005];...            % Error options for k

        % Localization lengths, lll loops over them
        Lc=[0.006;0.008;0.012;0.016;0.02;0.025;0.03;0.04;0.05;0.06;0.08];
        Transects={'Hor','Ver','Shoal','SChnl','Cakan'};                          % matters during the plotting stage (trttr)

        Names.folder       =    'multiall/';               % This changes automatically if twintest is 0
        Names.Lc           =    '_NoLc';                   % No Omega_Vec correction defaultname

        % Don't touch this part, all processed files were created by these subset limits
        % so only change it if you would like to create new set of assimilations!
        ltmin=110; ltmax=175;lnmin=120;lnmax=250;
        %     StrtTrm = [lnmin ltmin];
        % SizeTrm = [lnmax-lnmin ltmax-ltmin];
        coord=readgrid(Names.grid,ltmin,ltmax,lnmin,lnmax);
        clear ltmin ltmax lnmin lnmax
        %% READING AND COMBINING VARIABLES
        %  This part requires combining structure variables from multiple different .mat files. Thinking
        %  long term, combination of different time periods will also require blown-up members to be
        %  eliminated.

        %  To make it clear, if ensemble member 271 blew up you can't combine that with another time period.
        %  Hence I added these Ens.memb variables to show which ensemble member each Ens variables is from
        %  First should be finding missing Ens.memb files for each file being combined and then clearing those.

        nomemb =[];
        for iFile = 1:numel(FileListEns)
            FileData     = load(FileListEns{iFile});
            mlist=1:N;
            nomemb=[nomemb,mlist(~ismember(mlist,FileData.Ens.memb))];
        end
        nomemb= unique(nomemb(:))';

        % check whether any blew up members are included in the ensemble files and create cell arrays
        % for each file to clear in the next step
        for iFile = 1:numel(FileListEns)
            FileData     = load(FileListEns{iFile});
            indblw=[];
            for nn=1:length(nomemb)
                indblw=[indblw,find(nomemb(nn)==FileData.Ens.memb(:))];
            end
            Indblow{iFile}=indblw;
        end

        clear FileData field
        % Now lets combine the variables from multiple file sources

        % Ens, True and Prior are more straightforward but for observation only combine when the field
        % is there and I had to rename those files

        Obs.lon2=[]; Obs.lat2=[]; Obs.rms=[];Obs.timegmt2=[];Ens.hH=[];True.hH=[];Prior.hH=[];
        % This following line is created just for some tests
        Obs.rmswch=[];True.us_rho2=[];True.vs_rho2=[];Prior.us_rho2=[];Prior.vs_rho2=[];True.dp_rho2=[];
        Prior.dp_rho2=[];True.flag=[];Obs.noaa_f=[];
        for ff=1:length(fields)

            for iFile = 1:numel(FileListEns)

                FileDataEns    = load(FileListEns{iFile});
                FileDataTP     = load(FileListTP{iFile});

                FileDataEns.Ens=clearstr2(FileDataEns.Ens,Indblow{iFile});

                if isfield(FileDataEns.Ens,fields{ff}) % Only consider variables in the Ensemble

                    if(strcmp(fields{ff},'k_rho'))% | strcmp(fields{ff},'hs_rho')) %take inverse square of the k variables
                        FileDataEns.Ens.(fields{ff})=FileDataEns.Ens.(fields{ff}).^-2;
                        FileDataTP.True.(fields{ff})=FileDataTP.True.(fields{ff}).^-2;
                        FileDataTP.Prior.(fields{ff})=FileDataTP.Prior.(fields{ff}).^-2;
                        FileDataTP.Obs.(fields{ff})=FileDataTP.Obs.(fields{ff}).^-2;
                    end

                    Ens.(fields{ff}) = [Ens.(fields{ff});FileDataEns.Ens.(fields{ff})];
                    True.(fields{ff}) = [True.(fields{ff});FileDataTP.True.(fields{ff})];
                    Prior.(fields{ff}) = [Prior.(fields{ff});FileDataTP.Prior.(fields{ff})];
                    True.flag  = [ True.flag;ff*ones(length(FileDataTP.True.(fields{ff})),1)];
                    Obs.(fields{ff}) = [Obs.(fields{ff});FileDataTP.Obs.(fields{ff})];
                    Obs.lon2=[Obs.lon2;FileDataTP.Obs.lon];
                    Obs.lat2=[Obs.lat2;FileDataTP.Obs.lat];

                    % Combine Observation error if exists
                    % Define Obs.rmswch as either [1]: With Obs. RMSE or
                    % [0]: No Obs. RMSE. This switch is used in the later part of
                    % the code to redefine RMS values for Obs. error.
                    if isfield(FileDataTP.Obs,[fields{ff}(1:3),'rmse']) % Only consider variables in the Ensemble
                        Obs.rmswch=[Obs.rmswch;ones(length(FileDataTP.Obs.lat),1)];
                        Obs.rms=[Obs.rms;FileDataTP.Obs.([fields{ff}(1:3),'rmse'])];
                    else
                        Obs.rmswch=[Obs.rmswch;zeros(length(FileDataTP.Obs.lat),1)];
                        Obs.rms=[Obs.rms;zeros(length(FileDataTP.Obs.lat),1)];
                    end
                    % % %       Temporary % These are created just for some tests
%                                     Obs.noaa_f=[Obs.noaa_f;FileDataTP.Obs.noaa_f];
%                                     True.us_rho2=[True.us_rho2;FileDataTP.True.us_rho_ignore];
%                                     True.vs_rho2=[True.vs_rho2;FileDataTP.True.vs_rho_ignore];
%                                     True.dp_rho2=[True.dp_rho2;FileDataTP.True.dp_rho];
%                     
%                                     Prior.us_rho2=[Prior.us_rho2;FileDataTP.Prior.us_rho_ignore];
%                                     Prior.vs_rho2=[Prior.vs_rho2;FileDataTP.Prior.vs_rho_ignore];
%                                     Prior.dp_rho2=[Prior.dp_rho2;FileDataTP.Prior.dp_rho];
                    % % %       Temporary % These are created just for some tests

                    Obs.timegmt2=[Obs.timegmt2;FileDataTP.Obs.timegmt];
                    %Hh is repeated every time a new variable is considered
                    Ens.hH = [Ens.hH;FileDataEns.Ens.hH];
                    True.hH = [True.hH;FileDataTP.True.hH];
                    Prior.hH = [Prior.hH;FileDataTP.Prior.hH];
                end
            end
        end

        % This is same for all files
        Ens.h=FileDataEns.Ens.h;
        True.h=FileDataTP.True.h;
        Prior.h=FileDataTP.Prior.h;
        %
        %     figure(66);clf
        %     plot(True.k_rho,Obs.k_rho,'x');hold on;
        %     plot(Prior.k_rho,Obs.k_rho,'x');
        % %     for pp=1:length(Obs.k_rho)
        % %         text(True.k_rho(pp),Obs.k_rho(pp),num2str(pp))
        % %     end
        %     xlabel('Modeled k [using NOAA sigma, SWAN direction and model velocities]')
        %         limneg=0;limpos=0.2;grid on
        %         plot([limneg limpos],[limneg limpos])
        %     ylabel('Obs kmag');legend('True Model','Prior Model')

        % Now augment all the varibles
        Ens.var=[];   True.var=[]; Prior.var=[]; Obs.var =[]; %Obs.swchsum =0;
        Names.figname = [];
        if(swch.dostat); Names.statname = []; end

        for ff=1:length(fields)

            Names.figname=[Names.figname fields{ff}(1)];
            if(swch.dostat); Names.statname=[Names.statname fields{ff}(1:2)]; end

            Ens.var=[Ens.var;Ens.(fields{ff})];
            Prior.var=[Prior.var;Prior.(fields{ff})];
            True.var=[True.var;True.(fields{ff})];
            Obs.var=[Obs.var;Obs.(fields{ff})];

            % Concatenating the error values according to the field type
            if(fields{ff}(1:2)=='us'); rf=1;end
            if(fields{ff}(1:2)=='vs'); rf=2;end
            if(fields{ff}(1:2)=='hs'); rf=3;end
            if(fields{ff}(1)=='k');    rf=4;end
            if(fields{ff}(1:2)=='hh'); rf=5;end

            %             if(sum(Obs.rmswch)<1) %I.e., if no Obs.rms was given
            %                 sprintf('Prior and true error levels deneme %s',fields{ff})
            % define RVAL according to Error levels %v5

            err_rate=mean(abs(Prior.var-Obs.var))/mean(abs(Obs.var));
            %                 err_rate=mean(abs(True.var-Obs.var))/mean(abs(Obs.var));
            %
            %                 Rval1=[Rval1;2.5*err_rate*Obs.(fields{ff})];
            %                 Rval2=[Rval2;2*err_rate*Obs.(fields{ff})];
            %                 Rval3=[Rval3;1.5*err_rate*Obs.(fields{ff})];
            %                 Rval4=[Rval4;err_rate*Obs.(fields{ff})];
            %             else
            %             % define RVAL (seems to be the easiest way) %v4
            Rval1=[Rval1;repmat(Rvallist(rf,1),[length(Ens.(fields{ff})),1])];
            Rval2=[Rval2;repmat(Rvallist(rf,2),[length(Ens.(fields{ff})),1])];
            Rval3=[Rval3;repmat(Rvallist(rf,3),[length(Ens.(fields{ff})),1])];
            Rval4=[Rval4;repmat(Rvallist(rf,4),[length(Ens.(fields{ff})),1])];
            %             end

            Ens=rmfield(Ens,fields{ff});
            Prior=rmfield(Prior,fields{ff});
            True=rmfield(True,fields{ff});
            Obs=rmfield(Obs,fields{ff});
        end
        Ens.var_avg=mean(Ens.var,2);              %ensemble average
        Obs.np=numel(Obs.lat2);
        Rval= [Rval1,Rval2,Rval3,Rval4];
        clear Rval1 Rval2 Rval3 Rval4
        if(~swch.makeprior)
%             clear FileListEns FileDataEns FileDataTP FileListTP Indblow indblw
        end
        coord.ne=size(Ens.var,2);

        Unq.Unqtime        =    unique(Obs.timegmt2);
        Unq.Ut_N           =    numel(Unq.Unqtime);  % Finds how many timesteps in the code.

        % This part is nexessary for the unique time >>>>>>>
        % to simplify things, we will copy the original structure files before the unique loop and then
        % use them within the loop rather than reading everything again
        %     Ens2=Ens; Prior2=Prior; True2=True; Obs2=Obs; Rval2=Rval;
        %
        % %     for utt=[207:208] %Unq.Ut_N]  % unique time loop       % you can (and should) manually comment this out for bird observations, line ~1638
        %         Ens=Ens2; Prior=Prior2; True=True2; Obs=Obs2; Rval=Rval2;
        % <<<<<<< This part is nexessary for the unique time

        % %     % Find the rows corresponding to these logical statements
        % %     [indall,jndall]=find(abs(Ens.us_rho) > UVmax | abs(Ens.vs_rho) > UVmax);
        % %
        % %     indall=unique(indall);
        % %
        % %     %Now get rid of these outliers for Obs,True,Prior and Ens structures
        % %
        % %     True=clearstr(True,indall,'clear',Obs.np);
        % %     Prior=clearstr(Prior,indall,'clear',Obs.np);
        % %     Ens=clearstr(Ens,indall,'clear',Obs.np);
        % %
        % %     Obs=clearstr(Obs,indall,'clear',Obs.np);
        % %     Obs.np = numel(Obs.lat2);
        % %
        % %     clear indall jndall

        % Take out Obs observations outside the Area of Interest
        if(swch.trmpoly)
            % This is the optional part of trimming observations in irrelevant locations.
            % It uses inpolygon function which as the name suggests checks whether
            % given coordinates (i.e. Obs observations)are inside a predetermined
            % polygon. I chose those points manually and created a polygon with
            % 10 sides (10 points). And then plotted the inpolygon output for a visual
            % check

            % I also combined it with an option on time restrictions

            if(isempty(dir(Names.polygon)))
                disp('you are missing the polygon file')
                return
            else
                load(Names.polygon)
                IN = inpolygon(Obs.lon2, Obs.lat2,Poly.lon_rho, Poly.lat_rho );
            end

            disp(' ');
            disp(['You chose Trim option ' num2str(Trm.option) ', make sure this is what you want']);
            disp(' ');
            switch(Trm.option)
                case(1)
                    % Just the polygon (option 1)
                    indall=find(IN==1);
                    disp('No time limits, all observations are combined');
                    disp(' ');
                case(2)
                    disp(['Time limits according to the unique time steps, observations within' ...
                        ' a certain time window around the unique time step are considered']);
                    disp(' ');
                    % Polygon and time limit(option 2)

                    % Spring Multiple Timesteps
                    Sta.ts=Unq.Unqtime(utt)-0.000025/24;
                    Sta.te=Unq.Unqtime(utt)+0.000025/24;
                    Names.figname=[Names.figname num2str(utt)];
                    disp(['Unique time ' num2str(utt) ' of ' num2str(Unq.Ut_N) ' -> ' datestr(Unq.Unqtime(utt)) ])

                    indall=find((IN==1) & (Obs.timegmt2>Sta.ts & Obs.timegmt2<Sta.te));
                    if(isempty(indall))
                        continue
                    end

                case(3)
                    disp(['Combination of multiple time windows, this code is lazy to print them out so find them yourself']);
                    disp(' ');
                    % This step took a lot of tries to make it more precies but essentially it
                    % is combining two logical statements with an or statement
                    % (option 3)

                    % % % 6points birds dates
%                     indall=find((IN==1) & (Obs.timegmt2>datenum(2014,05,29,1,0,0) ...
%                         & Obs.timegmt2<datenum(2014,05,29,2,0,0))      ...
%                         |   (IN==1) & (Obs.timegmt2>datenum(2014,06,04,1,0,0) ...
%                         & Obs.timegmt2<datenum(2014,06,04,2,0,0))...
%                         );
 
                                        indall=find((IN==1) & (Obs.timegmt2>datenum(2014,05,27,1,0,0) ...
                                            & Obs.timegmt2<datenum(2014,05,27,19,0,0))      ...
                                            |   (IN==1) & (Obs.timegmt2>datenum(2014,05,27,18,0,0) ...
                                            & Obs.timegmt2<datenum(2014,05,28,2,0,0))...
                                            );

                case(4)
                    disp(['Depends on the bcbc location but this draws a circle' ...
                        ' around the mouth or upper estuary point ' newline 'and only considers' ...
                        ' those points for assimialtion']);
                    disp(' ');
                    %OPTION 4: Trim outside the localization area AND the unique time
                    %step, with the help of the viscircles function

%                     xo(1)=416501; yo(1)=5122176; desc{1}='mouth'; %%Cgdem1

                    xo(1)=416001; yo(1)=5123076; desc{1}='mouth'; %%Paper2
                    xo(2)=422000; yo(2)=5122583; desc{2}='upper estuary'; %%Cgdem2

                    [yo,xo]=utm2ll(xo,yo,10);
                    bcbc=1;
%                     [~,ppp_c1]=min(spheric_dist(yo(swch.ibc(bcbc)),Obs.lat,xo(swch.ibc(bcbc)),Obs.lon));                    [~,ppp_c1]=min(spheric_dist(yo(swch.ibc(bcbc)),Obs.lat,xo(swch.ibc(bcbc)),Obs.lon));

                    [~,ppp_c1]=min(spheric_dist(yo(bcbc),Obs.lat2,xo(bcbc),Obs.lon2));

                    %                     [~,ppp_c2]=min(spheric_dist(yo(2),Obs.lat,xo(2),Obs.lon));

                    f767=figure(767);clf;
                    pcolor(coord.lon_rho,coord.lat_rho,coord.h);hold on;
                    h2=viscircles([Obs.lon2(ppp_c1),Obs.lat2(ppp_c1)],Lc(lll));
%                     h2=viscircles([xo(1) yo(1)],Lc(lll));

                    plot(xo,yo,'yx')
                    Circle.X=h2.Children.XData; %lon
                    Circle.Y=h2.Children.YData; %lat

                    IN = inpolygon(Obs.lon2, Obs.lat2,Circle.X, Circle.Y );

                    %org
%                     Sta.ts=Unq.Unqtime(utt)-0.25/24;
%                     Sta.te=Unq.Unqtime(utt)+0.25/24;
%                     Names.figname=[Names.figname num2str(utt)];
                    %org

                    %June3 vematR
                    %     '03-Jun-2013 02:20:16'
                    %     '03-Jun-2013 13:47:27'
                    Sta.ts=datenum(2013,6,2,16,20,16)-0.25/24;
                    Sta.te=datenum(2013,6,3,5,47,27)+0.25/24;
                    indall=find((IN==1) & (Obs.timegmt2>Sta.ts & Obs.timegmt2<Sta.te));
                    close(f767)
                    clear h2


            end
        end
        % This part typically does not change
        kpclr='keep';
        True=clearstr(True,indall,kpclr,Obs.np);
        Prior=clearstr(Prior,indall,kpclr,Obs.np);
        Ens=clearstr(Ens,indall,kpclr,Obs.np);
        Rval=Rval(indall,:);
        Obs=clearstr(Obs,indall,kpclr,Obs.np);
        Obs.np = numel(Obs.lat2);
        clear indall IN

        % Eliminate Hsig inside the mouth
                if(sum(contains(fields,'hs_rho'))>0)
                    hsflag = find(contains(fields,'hs_rho') ==1);
                    %             fields
                    disp([' Eliminating Hsig inside the mouth '])
                    indall = find(Obs.lon2 > -124.05 & True.flag == hsflag);
        
                    kpclr='clear';
                    True=clearstr(True,indall,kpclr,Obs.np);
                    Prior=clearstr(Prior,indall,kpclr,Obs.np);
                    Ens=clearstr(Ens,indall,kpclr,Obs.np);
                    Rval(indall,:)=[];
                    Obs=clearstr(Obs,indall,kpclr,Obs.np);
                    Obs.np = numel(Obs.lat2);
                    clear indall IN hsflag
                end

        %%
        % Reshape truth and prior and Use observations if not a Twin Test

        Ens.h=reshape(Ens.h,[length(coord.lon_rho(:)) coord.ne]);

        % Use less ensemble testing part
        %     if(swch.Nens)
        %         fields = fieldnames(Ens);
        %         for ii=1:length(fields)
        %             Ens.(fields{ii}) = Ens.(fields{ii})(:,1:Trm.Nens);
        %         end
        %         clear fields
        %     end

        True.orgvar=True.var; %duplicate true model output for plotting purposes

        % Twin test or True test
        if(~swch.twin)
            %         for ff=1:length(fields)
            True.var= Obs.var;
            %             endr
            Names.folder='Realtest/';
        end

        %         % This part is to read in salt wedge >>>>>>>
        %         filehis='../Model_Outputs/2013_PRIOR_COUPLED/ocean_his.nc';
        %         ltmin=110; ltmax=175;lnmin=120;lnmax=250;
        %         StrtTrm = [lnmin ltmin];
        %         SizeTrm = [lnmax-lnmin ltmax-ltmin];
        %
        %         Pri2.time=ncread(filehis,'ocean_time');
        %         Pri2.time=Time.basedate+Pri2.time./86400;
        %         Time.timemin=findnearest(Unq.Unqtime(utt),Pri2.time);
        %         Pri2.s_b=ncread(filehis,'salt',[StrtTrm 1 Time.timemin],[SizeTrm 1 1]);
        %         % <<<<<<< This part is to read in salt wedge

        if(swch.filter) % tidal filter (using station output from ROMS) (didn't work well)

            dr='../Model_Outputs/2013_TRUTH_COUPLED/';file='ocean_sta.nc';
            vars={'u','v'};
            Sta = roms_readstation(dr,file,2013,vars,[2,3]);
            indfilt=[];
            for ff=1:length(Obs.timegmt2)

                iT=findnearest(Sta.timegmt,Obs.timegmt2(ff));
                %         ff
                %         disp(datestr(Obs.timegmt2(ff)) )
                %         disp(datestr(Sta.timegmt(iT)) )

                if(abs(Sta.MGT3.u(iT)) < 0.5)
                    indfilt=[indfilt ff];  %combine obs time steps when slack tide
                end
            end

            kpclr='clear';
            True=clearstr(True,indfilt,kpclr,Obs.np);
            Prior=clearstr(Prior,indfilt,kpclr,Obs.np);
            Ens=clearstr(Ens,indfilt,kpclr,Obs.np);
            Rval(indfilt,:)=[];
            Obs=clearstr(Obs,indfilt,kpclr,Obs.np);
            Obs.np = numel(Obs.lat2);
            clear indfilt
        end %filter1end
        %     Sta.MGT3.u_low=Sta.MGT3.u(abs(Sta.MGT3.u)<0.35);
        %
        %     [pks,locs]=findpeaks(abs(Sta.MGT3.u_low),'MinPeakHeight',0.3,'MinPeakDistance',2)
        %
        %     plot(abs(Sta.MGT3.u_low),'x-');hold on;
        %     plot(locs,pks,'o')
        %

        %                     disp(datestr(Obs.timegmt2(ff)) )
        %                     disp(datestr(Sta.timegmt(iT)) )
        %             figure(444)
        %     plot(Sta.timegmt,Sta.MGT3.u);hold on;
        %     plot(Obs.timegmt2,Obs.var*0,'x');
        %     plot(Obs.timegmt2,Sens./max(Sens(:)),'d');
        %     Sens = roms_ens_sens(Ens,Obs,utt,Prior,True);

        if(swch.filter2) % Adjoint filter, currently it is set for 0.05
            Sens = roms_ens_sens(Ens,Obs,1,Prior,True);

            indall=find(abs(Sens)>0.05);
            kpclr='keep';
            True=clearstr(True,indall,kpclr,Obs.np);
            Prior=clearstr(Prior,indall,kpclr,Obs.np);
            Ens=clearstr(Ens,indall,kpclr,Obs.np);
            Rval=Rval(indall,:);
            Sens=abs(Sens(indall));
            %     Sens = round2decimal(Sens,0.01);
            %     colors=turbo(length(unique(Sens)));
            Obs=clearstr(Obs,indall,kpclr,Obs.np);
            Obs.np = numel(Obs.lat2);
            clear indall
            %     colors=bluewhitered(Sens);
        end %filter2end

        if(swch.filter3) % std dev filter
            diffobs=abs(Prior.var-Obs.var);
            std_v=std(Ens.var')';
            indall=find(diffobs<stdvar*std_v);

            %             std2=std(True.orgvar-True.var);
            % %            std2= std(True.orgvar(True.flag==2)-True.var(True.flag==2));
            %             %
            %             fi=1;
            %
            %         limneg=min(True.orgvar(True.flag==fi));
            %         limpos=max(True.orgvar(True.flag==fi));
            %         figure(21);clf
            %         subplot(121)
            % %                 plot(True.orgvar(indall),True.var(indall),'x');hold on
            % %         plot(Obs.var(indall),Prior.var(indall),'x');hold on
            % %         plot(Obs.var(indall),True.orgvar(indall),'o')
            % %         plot(Obs.var(True.flag==fi),Prior.var(True.flag==fi),'x');hold on
            %         plot(Obs.var,True.orgvar,'o');hold on
            %         plot(Obs.var(indall),True.orgvar(indall),'x');hold on
            %
            %         plot([limneg limpos],[limneg limpos])
            %         plot([limneg+std2 limpos+std2],[limneg limpos])
            %         plot([limneg-std2 limpos-std2],[limneg limpos])
            %         xlabel('Observation');ylabel('model');legend('true','std')
            %         subplot(122)
            %         plot(Obs.var,Prior.var,'o');hold on
            %         plot(Obs.var(indall),Prior.var(indall),'x');
            %
            %         plot([limneg limpos],[limneg limpos])
            %         plot([limneg+std2 limpos+std2],[limneg limpos])
            %         plot([limneg-std2 limpos-std2],[limneg limpos])
            %         xlabel('Observation');ylabel('model');legend('prior','std')
            %
            %             %         plot(Prior.var(indall)+3*std(Ens.hH(indall,:)')')
            %             %         plot(Prior.var(indall)-3*std(Ens.hH(indall,:)')')
            % return
            kpclr='keep';
            True=clearstr(True,indall,kpclr,Obs.np);
            Prior=clearstr(Prior,indall,kpclr,Obs.np);
            Ens=clearstr(Ens,indall,kpclr,Obs.np);
            Rval=Rval(indall,:);
            Obs=clearstr(Obs,indall,kpclr,Obs.np);
            Obs.np = numel(Obs.lat2);
            clear indall diffobs std_v
        end
        % I tried to simplify this portion and combining variables if they are "active" using
        % the binary swch.var term.

        % Ens.var=Ens.var.^-2;    True.var=True.var.^-2;    Prior.var=Prior.var.^-2;    True.orgvar=True.orgvar.^-2;
        % Ens.var_avg=Ens.var_avg.^-2;
        % Onto the covariance part, C_uu represents the augmented variables
        Cov.C_uu=myCov(Ens.var,Ens.var);
        if(~swch.savemat)
            %             Cov.C_hh=myCov(Ens.h,Ens.h);  %Optional, so it's off for now
        end
        Cov.C_hu=myCov( Ens.h,Ens.var);

        % Localization
        if(swch.Lc)
            Cov.C_hu=Cov.C_hu.*omega_vec(coord.lat_rho, coord.lon_rho,...
                Obs.lat2,Obs.lon2,Lc(lll));
            Cov.C_uu=Cov.C_uu.*omega_vec(Obs.lat2,Obs.lon2,...
                Obs.lat2,Obs.lon2,Lc(lll));
            Names.Lc=['_Lc=' num2str(Lc(lll)*111.11*1000,4) 'm'];
        end

        if (swch.savemat) %does not save C.hh (takes too much disk space)
%                         save(['sinkout_using_' Names.figname],'True','Prior','Obs','Cov','Rval','Ens')
%             save(['sinkouts/circle/SARJun3_c1_' Names.figname],'True','Prior','Obs','Cov','Rval','Ens','Circle')
            save(['sinkouts/multivar_iter4_' Names.figname],'True','Prior','Obs','Cov','Rval','Ens')

            disp('Looks like swch.savemat is on, the code will stop running now.')
            return
        end
        %% Plots

%         Names.usace='usace_data/2019_filt_gridded_usace_bathy_200m.mat';  %Org
%         Names.usace='usace_data/gridded_usace_bathy_250m.mat'; % Dylan

        %                 Names.usace='usace_data/2019_gridded_usace_bathy_250m.mat';
%         usace2019=load(Names.usace);
%         usace2019.z=cell2mat(usace2019.z);
%         usace2019.z=usace2019.z{2};

        Names.usace='usace_data/2019_filt_gridded_usace_bathy_250m.mat';
        % Names.usace='usace_data/gridded_usace_bathy_250m.mat';
        usace2019=load(Names.usace);
        usace2019.z=cell2mat(usace2019.z);

        for trtr=[1,2]  %transect loop -> [1,2,3,4,5] or any combination
            % have to start from 1

            Trans.dir=Transects{trtr};
            swch.Hor=false; %default

            % Define Transect coordinates
            switch(Trans.dir)
                case('Hor')
                    % Horizontal Transect (Original is the top one)
                    %original old
                    %             Trans.xtrans(1)=-124.1; Trans.ytrans(1)=46.257;
                    %             Trans.xtrans(2)=-123.99; Trans.ytrans(2)=46.255;
%                     load('../Adjoint_density/SAR_trns_NChnl2.mat');
                    load('SAR_trns_NChnl_usace.mat');

                    Trans.xtrans=Poly.lon_rho;
                    Trans.ytrans=Poly.lat_rho;
                    clear Poly
                    swch.markJetA=true;                               %Dotted line for JettyA
                    swch.Hor=true;
                case('Ver')
                    % Vertical Transect
                    %CORRECT VALS
                                        Trans.xtrans(1)=-124.038; Trans.ytrans(1)=46.266;
                                        Trans.xtrans(2)=-124.039; Trans.ytrans(2)=46.235;
                    %Alternative
%                     Trans.xtrans(1)=-123.95; Trans.ytrans(1)=46.25;
%                     Trans.xtrans(2)=-123.97; Trans.ytrans(2)=46.215;
                    swch.markJetA=false;                               %Dotted line for JettyA
                    swch.Hor=false;
                case('Shoal')
                    % Transect at the shoal

                    Trans.xtrans(1)=-124.1; Trans.ytrans(1)=46.248;
                    Trans.xtrans(2)=-124.085; Trans.ytrans(2)=46.21;
                    swch.markJetA=false;
                    swch.Hor=false;

                case('SChnl')
                    % Transect at the shoal
                    %             Trans.xtrans(1)=-124.007; Trans.ytrans(1)=46.257;
                    %             Trans.xtrans(2)=-123.9; Trans.ytrans(2)=46.233;
                    %             swch.markJetA=false;
                    load('../Adjoint_density/SAR_trns_SChnl.mat');
                    Trans.xtrans=Poly.lon_rho;
                    Trans.ytrans=Poly.lat_rho;
                    clear Poly
                    swch.markJetA=true;                               %Dotted line for JettyA
                    swch.Hor=true;
                case('Cakan')
                    % show a transect (picked by hand)
                    %                     Trans.xtrans(1)=418704; Trans.ytrans(1)=5124211;
                    %                     Trans.xtrans(2)=424311; Trans.ytrans(2)=5122245;
                    Trans.xtrans(1)=-124.0550; Trans.ytrans(1)=46.2666;
                    Trans.xtrans(2)=-123.9819; Trans.ytrans(2)=46.2496;
            end

            % -124.038     46.266     "JETTY A" (kind of where "the mouth" is)

            Trans.JetLat=46.266;
            Trans.JetLon=-124.038;
            [Trans.eastJet,Trans.northJet] = ll2UTM(Trans.JetLon,Trans.JetLat,23,'10T');
            Trans.eastJet=Trans.eastJet./1000;                         %convert to km
            Trans.northJet=Trans.northJet./1000;                       %convert to km

            if(trtr==1)
                for ir=1:4
                    %                             ir=3;
                    %                 if(swch.twin)
                    clear R
                    % If this is not a twin test, use the existing rmse from the Observations.
                    if(swch.twin==0)
                        Rval(Obs.rmswch==1,ir)=Obs.rms(Obs.rmswch==1);
                    end
                    R=Rval(:,ir).^2.*eye(length(Ens.var_avg));
                    %  R=Rval*Rval'; %*eye(length(Ens.var_avg));
                    %                 else
                    %                     R=[Obs.u_rms;Obs.v_rms];
                    %                     R=diag(R.^2)./1;
                    %                 end

                    Ci=inv(R+Cov.C_uu);

                    inov_1=Ci*(True.var-Ens.var_avg);
                    if(swch.prior)
                        inov_1=Ci*(True.var-Prior.var);
                    end

                    inov=Cov.C_hu*inov_1;                inovvar=Cov.C_uu*inov_1;
                    %                     Cov.C_hh_a=Cov.C_hh-Cov.C_hu*Ci*Cov.C_hu'; %Posterior Uncertainty
                    %Cov.C_uu*inov_l is for U innovations

                    display (['sum of inovations=', num2str(sum(inov(:))) 'm for Rval =' num2str(Rval(ir))])
                    innovation_2d(:,:,ir)=reshape(inov,[size(coord.lat_rho)]);
                    %Plotting

                    % %                     innovation_2d(ir,:,:)(coord.maskR==0)=NaN;
                    % %                             subplot(4,1,ir)
%                     True.h( coord.maskR==0)=NaN;
%                     Prior.h( coord.maskR==0)=NaN;
                    %
                    if(swch.onscreen)
                        axes(hf3(ir));cla;
                        pcolor( coord.lon_rho', coord.lat_rho',( Prior.h+innovation_2d(:,:,ir))')
                        % shading flat
                        caxis([0 35]);% daspect([1 1 1]);
                        hcb=colorbar;
                        colormap(flipud(haxby(64)));shading interp
                        title(Rtitle{ir})

                        if(swch.Lc)
                            fig3title.String= [Names.figname,' Assimilation using L_c = ',num2str(Lc(lll)*111.11*1000,4),...
                                'm'];
                            %                     fig3title.String=['U velocity at MGT3 = ' num2str(Sta.MGT3.u(iT)) ...
                            %                         ', Sens U,V = ' num2str(Sens)];
                            drawnow
                        end
                    end % onscreen loop
                    % % %                             subplot(4,1,1)
                end    % ir loop
            end  % you only need to run this once
            %
            if(swch.onscreen)

                axes(hf3(5));cla;
%                 pcolor( coord.lon_rho', coord.lat_rho',True.h')
                pcolor(usace2019.lon,usace2019.lat,usace2019.z);hold on; caxis([0 35]);
                %                                     hcb=colorbar; colormap(flipud(haxby(64))); shading interp
                % shading flat
                caxis([0 35]); %daspect([1 1 1]);
                hcb=colorbar;
                colormap(flipud(haxby(64))); shading interp
                hold on;
                %         plot(Trans.xtrans,Trans.ytrans,'r-','LineWidth',1)
                plot(Obs.lon2,Obs.lat2,'ko','MarkerSize',Plt.ms*2,'MarkerFaceColor','k')
                %                 for pp=1:length(Obs.lon2)
                %                     text(Obs.lon2(pp),Obs.lat2(pp),num2str(pp));%datestr(Obs.timegmt(pp)))
                %                 end
                title(['True Bathymetry'],'FontSize',Plt.fs3)
                %         load(Names.polygon)
                %         plot(Poly.lon_rho,Poly.lat_rho)
                hold off

                axes(hf3(6));cla;
                pcolor( coord.lon_rho', coord.lat_rho',Prior.h')
                caxis([0 35]); %daspect([1 1 1]);
                hcb=colorbar;
                colormap(flipud(haxby(64))); shading interp
                hold on;       % plot(Sta.MGT3.lon,Sta.MGT3.lat,'x')
                %       gscatter(Obs.lon2,Obs.lat2,Sens,colors)
                plot(Trans.JetLon,Trans.JetLat,'ko','MarkerSize',Plt.ms+2,'MarkerFaceColor','k')
                plot(Trans.xtrans,Trans.ytrans,'r-','LineWidth',2)
                text(Trans.JetLon,Trans.JetLat+0.01,'Jetty A')

                title('Prior Bathymetry','FontSize',Plt.fs3)
                hold off
            end %onscreen loop
            %         % This part is to plot salt wedge >>>>>>>
            %
            %             axes(hf3(6));cla;
            %             pcolor(coord.lon_rho',coord.lat_rho',Pri2.s_b');shading interp;hold on;
            %             plot(Obs.lon2,Obs.lat2,'ko','MarkerSize',Plt.ms*3,'MarkerFaceColor','k')
            %
            %             caxis([0 35]);% daspect([1 1 1]);
            %             colormap(hf3(6),cm_thermal(128));colorbar;
            %             end
            %         % <<<<<<< This part is to plot salt wedge

            % Transects
            if(~swch.Hor)
                Trans.xtrans=linspace(Trans.xtrans(1),Trans.xtrans(2),50);
                Trans.ytrans=linspace(Trans.ytrans(1),Trans.ytrans(2),50);
            end
            % MODEL DATA
            [Trans.eastT,Trans.northT] = ll2UTM(Trans.xtrans,Trans.ytrans,23,'10T');
            Trans.eastT=Trans.eastT./1000;                                 %convert to km
            Trans.northT=Trans.northT./1000;                                %convert to km

            if(swch.onscreen)

                %                                     f11=figure(1)
                %                                     pcolor(usace2019.lon,usace2019.lat,usace2019.z);hold on; caxis([0 35]);
                %                                     hcb=colorbar; colormap(flipud(haxby(64))); shading interp
                %                                     load('../Adjoint_density/SAR_trns_NChnl2.mat');
                %                                     plot(Poly.lon_rho,Poly.lat_rho)
                %                                     clear Poly
                %                                     load('../Adjoint_density/SAR_trns_SChnl.mat');
                %                                     plot(Poly.lon_rho,Poly.lat_rho)
                %
                %                                     f22=figure(2)
                %                                     pcolor( coord.lon_rho', coord.lat_rho',True.h');caxis([0 35]);
                %                                     hcb=colorbar; colormap(flipud(haxby(64))); shading interp;hold on;
                %                                     load('../Adjoint_density/SAR_trns_NChnl2.mat');
                %                                     plot(Poly.lon_rho,Poly.lat_rho)
                %                                     clear Poly
                %                                     load('../Adjoint_density/SAR_trns_SChnl.mat');
                %                                     plot(Poly.lon_rho,Poly.lat_rho)
                %
                %                                     export_fig(f11,'2014usace',Plt.type,Plt.res,'-p0.01');
                %                                     export_fig(f22,'2014profile',Plt.type,Plt.res,'-p0.01');
                for ir=1:4
                    %                                             break
                    set(hf2(ir),'XLabel',[]);        %Plot then clear the x handl
                    Trans.hfi(ir,:)=griddata( coord.lon_rho, coord.lat_rho, Prior.h,Trans.xtrans,...
                        Trans.ytrans);
                    Trans.hai(ir,:)=griddata( coord.lon_rho, coord.lat_rho,( Prior.h+innovation_2d(:,:,ir)),...
                        Trans.xtrans,Trans.ytrans);
                    Trans.hti(ir,:)=griddata( coord.lon_rho, coord.lat_rho, True.h,Trans.xtrans,...
                        Trans.ytrans);
                    Trans.hui(ir,:)=griddata( usace2019.lon, usace2019.lat, usace2019.z,Trans.xtrans,...
                        Trans.ytrans);

                    axes(hf2(ir));cla;
                    % f2.WindowState = 'maximized';
                    hold on;
                    if(ir >2); xlabel('Distance along the Transect [km]'); end
                    switch(swch.markJetA)

                        case (false)
                            ltrans=sqrt((Trans.eastT-Trans.eastT(1)).^2+(Trans.northT-...
                                Trans.northT(1)).^2);
                            plot(ltrans,Trans.hfi(ir,:),'r','linewidth',1.5)
                            plot(ltrans,Trans.hai(ir,:),'b','linewidth',1.5)
                            %                             plot(ltrans,Trans.hti(ir,:),'k','linewidth',1.5)        %2013 model grid truth
                            plot(ltrans,Trans.hui(ir,:),'k','linewidth',1.5)        %USACE Truth

                            if(contains(Trans.dir,'Ver'))
                                if(ir >2);  xlabel('Distance from Jetty A [km]'); end
                            end
                        case(true)

                            %vertical line at JettyA
                            Trans.yJet = [0:35];
                            Trans.xJet=(Trans.eastJet-Trans.eastT(1))*ones(1,length(Trans.yJet));

                            plot(Trans.eastT-Trans.eastT(1),Trans.hfi(ir,:),'r','linewidth',1.5)
                            plot(Trans.eastT-Trans.eastT(1),Trans.hai(ir,:),'b','linewidth',1.5)
                            %                             plot(Trans.eastT-Trans.eastT(1),Trans.hti(ir,:),'k','linewidth',1.5)        %2013 model grid truth
                            plot(Trans.eastT-Trans.eastT(1),Trans.hui(ir,:),'k','linewidth',1.5)        %USACE Truth

                            plot(Trans.xJet,Trans.yJet,'.','linewidth',1.5)
                            text(Trans.xJet(1),Trans.yJet(end),'Jetty A')
                            %                         set(gcf,'paperposition',[0 0 8 3])
                    end

                    if(ir==1); legend('Prior','Posterior','Truth'); end
                    %                     if(ir==1); legend('Prior','Posterior','Truth'); end

                    box on
                    ylabel('Depth [m]');
                    set(gca,'fontsize',12)
                    %             title(['R=',num2str(Rval(ir)),'^2'],'FontSize',8)
                    title(Rtitle{ir})
                    hold off;

                end
            end %onscreen
            if(swch.doprint)

                %             %drifter outputs >>>>>>>>>>
                %             export_fig(f2,strcat(Names.folder,'drifter2/',...
                %                 Names.pert,'_',Trans.dir,'_Transect_',Names.figname,Names.Lc),... %num2str(zstr(utt,3))),...
                %                 Plt.type,Plt.res,'-p0.01');%,'-m2')
                %
                %             export_fig(f3,strcat(Names.folder,'drifter2/',...
                %                 Names.pert,'_Mouth_',Names.figname,Names.Lc),...         %num2str(zstr(utt,3))),...
                %                 Plt.type,Plt.res,'-p0.01');
                %             %<<<<<<<drifter outputs

                %             %org outputs >>>>>>>>>>
                export_fig(f2,strcat(Names.folder,'UV',num2str(UVmax),'/',birdavg,...
                    Names.pert,'_',Trans.dir,'_Transect_',Names.figname,Names.Lc,num2str(stdvar)),...
                    Plt.type,Plt.res,'-p0.01');%,'-m2')

                export_fig(f3,strcat(Names.folder,'UV',num2str(UVmax),'/',birdavg,...
                    Names.pert,'_Mouth_',Names.figname,Names.Lc,num2str(stdvar)),...
                    Plt.type,Plt.res,'-p0.01');
                %  <<<<<<<<<<<<<<    %org outputs

            end
            clear Trans           %necessary for the transect loop
            %                     pause
        end % End of the transect loop
        %     pause

        %% Statistics
        if(swch.dostat)
            Names.polyarea = 'RMSEtest.mat';

            Names.Nens =[ Names.statname '.txt'];

            if swch.twin==0
                Names.Nens =[ 'Real_' Names.Nens];
            end
            % RMSE for area
            load(Names.polyarea);
            %         f22=figure(22); clf;
            % %         subplot(211)
            %         pcolor(coord.lon_rho',coord.lat_rho',True.h'); hold on;
            %         caxis([0 35]); daspect([1 1 1]); hcb=colorbar;
            %         colormap(flipud(haxby(64)));shading flat
            %         plot(Area.Poly.lon_rho,Area.Poly.lat_rho,'b-','MarkerSize',3,'MarkerFaceColor','k')
            %         subplot(212)
            %         pcolor(coord.lon_rho',coord.lat_rho',Prior.h'); hold on;
            %         caxis([0 35]); daspect([1 1 1]); hcb=colorbar;
            %         colormap(flipud(haxby(64)));shading flat
            %         plot(Area.Poly.lon_rho,Area.Poly.lat_rho,'b-','MarkerSize',3,'MarkerFaceColor','k')

            %     if(swch.Nens==0);Trm.Nens=coord.ne;end
            for ir=1:4

                h_inov= Prior.h+innovation_2d(:,:,ir);
                IN = inpolygon(coord.lon_rho, coord.lat_rho,Area.Poly.lon_rho, Area.Poly.lat_rho );
                fmt='%4.0f   %5.2f   %5.2f   %5.2f %8.2f %5.2f %8.2f %4.0f %4.0f \n';

                %             for pp=1:numel(coord.lon_rho)
                %                 text(coord.lon_rho(pp),coord.lat_rho(pp),num2str(IN(pp)))
                %             end
                Area.indall=find(IN==1);

                %                 A=h_inov(Area.indall);
                %                 B=True.h(Area.indall);
                usace2019.h=griddata(usace2019.lon,usace2019.lat,usace2019.z,coord.lon_rho,coord.lat_rho);

                A=h_inov(Area.indall);
                % A=Prior.h(Area.indall);

                %         B=True.h(Area.indall);
                B=usace2019.h(Area.indall);
                [stat.r,stat.d,stat.b,stat.cc,stat.si,stat.ms]=rmse(A,B);

                if(isempty(dir(Names.Nens)))
                    fid = fopen(Names.Nens,'w');
                    fprintf(fid, 'Lc     RMSE     CC      Bias     D       si       MS   IR   std  \n');
                    fprintf(fid, fmt,...
                        Lc(lll)*111.11*1000,stat.r,stat.cc,stat.b,stat.d,stat.si,stat.ms,ir,stdvar);
                    fclose(fid);

                else
                    fid = fopen(Names.Nens, 'a+');
                    fprintf(fid, fmt,...
                        Lc(lll)*111.11*1000,stat.r,stat.cc,stat.b,stat.d,stat.si,stat.ms,ir,stdvar);
                    fclose(fid);

                end

            end
            %             R=corrcoef(A,B);
            %             [R,stat.d,stat.b,stat.cc,stat.si,stat.ms]=rmse(A,B);
            %
            %             disp([ ' rmse =  ' num2str(R)]);
            %             disp(' ')
        end

        %         Names.figname = [];%unique time loop
        %     end %unique time loop
    end % localization loop
end % std loop

%% Create Prior (Optional)
if(swch.makeprior)
    for i=2:2
        makeprior(Prior,innovation_2d,i,coord,1,'mcr_200m_real2_iter4',FileListTP);
    end
end

%%

% if(~swch.twin && ~swch.doprint && swch.onscreen)
if(~swch.twin && swch.onscreen)

    disp('Real Test -> Plotting Model vs. Observations')
    f44=figure(44); clf; fig44title = sgtitle('some initial title');

    set(gcf,'units','centimeters','position',[5 5 40 40]);
    [hf44,~]=tight_subplot(2,2,0.05,[0.05 0.1],[0.05 0.05]);
    for ff=1:length(fields)

        indall=find(True.flag==ff);
        axes(hf44(ff));cla;

        limneg=min(True.orgvar(indall));
        limpos=max(True.orgvar(indall));
%         plot(True.orgvar(indall),True.var(indall),'x');hold on

        plot(True.orgvar(indall),Prior.var(indall),'x');hold on
        plot(Obs.var(indall),Prior.var(indall)+inovvar(indall),'o')
        plot([limneg limpos],[limneg limpos])
        % xlabel('True k (model output from true bathy)')
        % ylabel('Prior and Posterior k error')
        % legend('Prior-True','Posterior-True');title('SAR inversion')

        % xlabel('True u,v and k (model output from true bathy)')
        % ylabel('Prior and Posterior u,v and k values')
        legend('Prior','Posterior');
        fig44title.String= 'Red Circles Indicate Model Improvements';
        %         title('SAR + Drifter inversion')

        % Model - Observation Comparison
        ylabel([fields{ff} ' (model output)'],'Interpreter','None')
        xlabel(['Observed ' fields{ff} ],'Interpreter','None');
        if(fields{ff}(1:3)=='k_r')
            ylabel('$$k^{-2} $$(model output)','Interpreter','Latex')
            xlabel('Observed $$k^{-2}$$','Interpreter','Latex');
        end
        [r,d,b,cc,~,~]=rmse(Obs.var(indall),Prior.var(indall));
        [rp,d,b,ccp,~,~]=rmse(Obs.var(indall),Prior.var(indall)+inovvar(indall));

        text(min(Obs.var(indall)),max(Prior.var(indall)),['Prior RMSE = ' num2str(r) newline 'Posterior RMSE = ' num2str(rp)]);
        set(gca,'fontsize',12);
        rf(ff)=r;ccf(ff)=cc;
        % legend('Prior-True','Posterior-True');title('SAR inversion')

    end

    %     f45=figure(45); clf; fig44title = sgtitle('some initial title');
    %
    %     set(gcf,'units','centimeters','position',[5 5 40 40]);
    %     [hf45,~]=tight_subplot(2,2,0.05,[0.05 0.1],[0.05 0.05]);
    %
    %     for ff=1:length(fields)
    %         indall=find(True.flag==ff);
    %         axes(hf45(ff));cla;
    %         plot(Obs.var(indall),'o');hold on
    %         plot(Prior.var(indall),'x');
    %         plot(Prior.var(indall)+3*std(Ens.var(indall,:)')')
    %         plot(Prior.var(indall)-3*std(Ens.var(indall,:)')')
    %
    %     end
    % for pp=1:length(Obs.lon2)
    %     text(True.orgvar(pp),True.var(pp),num2str(pp));%datestr(Obs.timegmt(pp)))
    % end

end

% export_fig(f44,['drifter_hs_' num2str(stdvar) '_improvement'],Plt.type,Plt.res,'-p0.01');
% save('Articleworthy/20132014_Comparisons/birds_iter3','True','Prior','Obs','fields','inovvar','FileListTP')

return
%% To plot linear vs nonlinear k/h comparison
% You need to comment the k-2 part and uncomment the vs_rho2 etc part while
% the the matrices are concatenated.
% For SAR, find corresponding wave numbers

%use only these two
%                         '2013/SARJun3_trupri_pert2.mat';...
%                         '2013/SARJun3_13_trupri_pert2.mat';...

% And uncomment this section and comment the k_Rho^2 part
%                                     Obs.noaa_f=[Obs.noaa_f;FileDataTP.Obs.noaa_f];

True.dp_rhon = 270 - True.dp_rho2;
kbarxp = cosd(True.dp_rhon);
kbaryp = sind(True.dp_rhon);
True.k_rho = disper(Obs.noaa_f*2*pi(),True.hH,True.us_rho2,True.vs_rho2,kbarxp,kbaryp);

Prior.dp_rhon = 270 - Prior.dp_rho2;
kbarxp_p = cosd(Prior.dp_rhon);
kbaryp_p = sind(Prior.dp_rhon);
% Prior.k_rho = disper(SAR.noaa_f*2*pi(),Prior.hH,Prior.us_rho2,Prior.vs_rho2,kbarxp,kbaryp);

% Ens.dp_rhon = 270 - Ens.dp_rho;
% kbarxp = cosd(Ens.dp_rhon);
% kbaryp = sind(Ens.dp_rhon);
%
% Ens.k_rho = disper(SAR.noaa_f*2*pi(),Ens.hH,Ens.us_rho,Ens.vs_rho,kbarxp,kbaryp);
% True.h( coord.maskR==0)=NaN;
% Prior.h( coord.maskR==0)=NaN;
Plt.msz=3.5;
Plt.lsz=10;
%quick visual check for the ensemble distribution, to make sure all is well
f5=figure(5);clf
set(gcf,'units','centimeters','position',[5 5 40 40]);
[hf5,~]=tight_subplot(3,1,0.05,[0.05 0.05],[0.05 0.05]);
        
% Good examples -> ppp_c1 =    21,370 (with june3 and june3_13

% for ppp_c1=1:length(Ens.hH)
for ppp_c1=[21,370]

    if(Obs.lat2(ppp_c1) < 46.25 & Obs.lat2(ppp_c1) > 46.235 & Obs.lon2(ppp_c1) > -124.08)
        ppp_c1
        dh=Ens.hH (ppp_c1,:);
        du=Ens.var (ppp_c1,:);
        p=polyfit(dh,du,1);
        fit=polyval(p,2:30);

        %     pkp=pade_disper(Obs.noaa_f(ppp_c1),2:30);
        pkt = disper(Obs.noaa_f(ppp_c1)*2*pi(),2:30,True.us_rho2(ppp_c1),True.vs_rho2(ppp_c1),...
            kbarxp(ppp_c1),kbaryp(ppp_c1));
        pkp = disper(Obs.noaa_f(ppp_c1)*2*pi(),2:30,Prior.us_rho2(ppp_c1),Prior.vs_rho2(ppp_c1),...
            kbarxp_p(ppp_c1),kbaryp_p(ppp_c1));
        % combined)
        axes(hf5(1));cla;
        pcolor( coord.lon_rho', coord.lat_rho',True.h')
        caxis([0 35]);cb=colorbar;xtickformat('degrees');ytickformat('degrees');cb.Label.String ='Depth [m]';    
        colormap(flipud(haxby(64))); shading interp; cb.Label.FontSize =10; hold on;
        plot(Obs.lon2(ppp_c1),Obs.lat2(ppp_c1),'ko','MarkerSize',2,'MarkerFaceColor','k');
        T1=text(min(coord.lon_rho(:))+0.005,max(coord.lat_rho(:))-0.01,'a)');
        T1.BackgroundColor='w';T1.FontWeight='bold';T1.FontSize=Plt.lsz;
        axes(hf5(2));cla;
        plot(dh,du,'.','MarkerSize',10);hold on;grid on;
        ss=scatter(True.hH(ppp_c1),True.orgvar(ppp_c1),'LineWidth',Plt.msz,'Marker','d','MarkerEdgeColor','k')
        scatter(Prior.hH(ppp_c1),Prior.var(ppp_c1),'LineWidth',Plt.msz,'Marker','d','MarkerEdgeColor','r');
        xlabel('Depth [m]'); ylabel('k [1/m]');
        T1=text(0.5,0.14,'b)');T1.FontWeight='bold';T1.FontSize=Plt.lsz;
        %     plot([dh True.orgvar(ppp_c1)],fit,'m--');
        plot(2:30,fit,'b--');
%         plot(2:30,pkt,'-m'); 
        pl=plot(2:30,pkp,'-r');pl.Color(4) = 0.7;
        legend('ens','true','pri','linfit','pridisp');
        xlim([0 30]);ylim([0.02 0.15]);
        %updated fit
        p=polyfit(dh,du.^-2,1);
        fit=polyval(p,2:30);
        axes(hf5(3));cla;
        plot(dh,du.^-2,'.','MarkerSize',10);hold on;grid on;
        scatter(True.hH(ppp_c1),True.orgvar(ppp_c1).^-2,'LineWidth',Plt.msz,'Marker','d','MarkerEdgeColor','k')
        scatter(Prior.hH(ppp_c1),Prior.var(ppp_c1).^-2,'LineWidth',Plt.msz,'Marker','d','MarkerEdgeColor','r');
        xlabel('Depth [m]'); ylabel('k^{-2} [m^2]');
        %     plot([dh True.orgvar(ppp_c1)],fit,'m--');
        plot(2:30,fit,'b--');
        T1=text(0.5,375,'c)');T1.FontWeight='bold';T1.FontSize=Plt.lsz;

%         plot(2:30,pkt,'-m'); pl=plot(2:30,pkp,'-k');pl.Color(4) = 0.3;
        legend('ens','true','pri','linfit');
        xlim([0 30]);ylim([50 400]);
%                 export_fig(f5,strcat('Articleworthy/kplots/',num2str(ppp_c1)),'-pdf','-p0.01');
a = annotation('rectangle',[0 0 1 1],'Color','w');
                exportgraphics(f5,'Articleworthy/kplots/asda.pdf','ContentType','vector')
        pause(2)
    end
end