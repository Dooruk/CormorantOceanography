clear;%close all;
addpath('~/Documents/MATLAB/export_fig');
addpath('~/Documents/MATLAB/customcolormap');
addpath('~/Documents/MATLAB/rgb');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code that is only meant for plotting the processed variables from the
% Kitchen_sink_assimilate script.

% Currently it is loading two different Priors. Original prior (iter1) for
% plotting purposes and Prior_iter4 for posterior calculation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Included ThreeAreas.mat and defined three regions to quantify
% improvements for these parts of the MCR (Outside,Mouth,Inside). Inside
% was chosen this way due to the lack of observations and vicinity to the
% smoothing domain.

% TODO: Include k_mod vs. k_obs with the following
% scatter(True.var.^(-1/2),True.orgvar.^(-1/2));hold on;
% plot([0.05,0.14],[0.05,0.14])

% Drifter u comp
% scatter(True.var(True.flag==1),True.orgvar(True.flag==1))
% v comp
% scatter(True.var(True.flag==2),True.orgvar(True.flag==2))
% also includes convergence plots

%*********** k is automatically converted to the -2 power (to linearize)*********\

% Original Prior is loaded before transects at the end for plotting purposes
% Plotted for irir error levels (1 to 4, high to low Obs_RMSE)
irir=[2:3];

Transects={'Hor','Ver','Shoal','SChnl','Cakan'};                          % matters during the plotting stage (trttr)
Lc=[0.006;0.008;0.012;0.016;0.02;0.025;0.03;0.04;0.05;0.06;0.08];
lll=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
swch.Lc            =    1;                         % Omega_Vec correction (Localization on/off)
swch.prior         =    0;                         % use prior for the innovation [0=use ensemble average]
swch.twin          =    0;                         % Calculate Cov. for Hh [1] i.e. obs. points or h [0] i.e. the whole domain
swch.trmpoly       =    1;                         % Ignore points outside the polygon (Arbitrarily chosen)
swch.doprint       =    0;                         % print figures
swch.onscreen      =    1;                         % Display figure on screen or not
swch.filter        =    0;                         % Tidal filter (under development)
swch.filter2       =    0;                         % Adjoint Sensitivity filter
swch.filter3       =    1;                         % Std.Dev filter
Trm.option         =    1;                         %[Trimming options: (1) only coord., (2) time and coord.,
%(3) multiple time and coord., (4) time and trim outside localization
%any option other than 1 requires unique loop
swch.dostat        =    1;
swch.naninterp     =    0;                         % naninterp during statistical comparisons (to compensate for USACE survey gaps)
swch.makeprior     =    0;
swch.panel5        =    1;                         % Plotting inside the mouth with 5 panels or not
Names.pert='pert2';                                % pert5 or pert2 (sigma_h (2m or 5m) for the ensemble members)

pos.post = [0.078171,0.46,0.4,0.2];
pos.true = [0.078171,0.728,0.4,0.2];
pos.ver = [0.5239,0.468,0.42139,0.45782];
if(swch.panel5)
    pos.hor = [0.5239,0.078,0.42139,0.3]; %new horizontal
else
    pos.hor = [0.078171,0.078,0.85,0.3]; %org horizontal
end
pos.shoal = [0.068171,0.078,0.4,0.3]; %new shoal

% Plot constants
Plt.ms      = 0.55;                  % Marker size
Plt.ms2     = 3;
Plt.covfac  = 12;
Plt.fs3     = 8;                   % Figure 3 Title font
Plt.fs      = 20;                   % axes titles
Plt.skl     = 0.6;                  % Only relevant for bird obs. Cuts down the low skill measurements (it was 0.75)
Plt.ulim    = 3;
Plt.res     = '-r600';            % Resolution of the plots
Plt.type    = '-png';            % Resolution of the plots
Plt.psize   = 20;                %pointsize for markers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Log:
%    - Variable names change: us,vs - surface velocities (u_surf)
%                             ub,vb - depth avg. velocities (u_bar)
%                             hs,l,p - sign. wave height, wave length, peak period
%                             dp    - peak wave direction

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
% this code). At time I manually edited these files and put them into corresponding folders. Specify your
% original ensemble member sample size (N). All observation structures should be named 'Obs'.

N=300;
Time.basedate=datenum(2013,5,1,0,0,0);                 %model start time
Names.grid='../../romsprep/mcr_200m_2013_grd_v1.nc';
Names.polygon='Polygon_Pts6.mat';

Rval1=[];Rval2=[];Rval3=[];Rval4=[];
Rtitle={'Highest R','Higher but Realistic R','Lower but Realistic R','Lowest R'};

%v5
% Rvallist=[0.65,0.5,0.35,0.25;...             % Error options for us_rho
%     0.4,0.3,0.2,0.15;...             % Error options for vs_rho
%     0.75,0.5,0.1,0.05;...             % Error options for hs_rho
%     80,60,40,20;...            % Error options for k (because k^2)
%     %     0.025,0.01,0.005,0.0025;...            % Error options for k
% %     3,2,1,0.5];           % Error options for hh_rho %v5 original

% Don't touch this part, all processed files were created by these limits
% so only change it if you would like to create new set of assimilations!
ltmin=110; ltmax=175;lnmin=120;lnmax=250;
%     StrtTrm = [lnmin ltmin];
% SizeTrm = [lnmax-lnmin ltmax-ltmin];
coord=readgrid(Names.grid,ltmin,ltmax,lnmin,lnmax);
clear ltmin ltmax lnmin lnmax

% Read USACE Survey
Names.usace='usace_data/2019_v4_gridded_usace_bathy_250m.mat';
% Names.usace='usace_data/gridded_usace_bathy_250m.mat';
usace2019=load(Names.usace);
usace2019.z=cell2mat(usace2019.z);

Names.usace='usace_data/2014_v3_gridded_usace_bathy_250m.mat';
usace2014=load(Names.usace);
% usace2019=load(Names.usace);

usace2014.z=cell2mat(usace2014.z);
% usace2019.z=usace2019.z{2}; %1 -> 2014, 2 -> 2019, 3 -> 2020
% usace2014.z=cell2mat(usace2014.z);

% % % % % Some temporary plotting >>>>>>>>>>>>>>>
% % % % % These are few observations during bird1 timeframe, I am using them to
% % % % % create a figure that shows the gradient descent process
% 
% Iter1=load('sinkouts/multivar_iter1_uvhk.mat');
% Iter2=load('sinkouts/multivar_iter2_uvhk.mat');
% Iter3=load('sinkouts/multivar_iter3_uvhk.mat');
% Iter4=load('sinkouts/multivar_iter4_uvhk.mat');
% 
% % %p_oe = observation error, p_me = model error
% % flags indicate the observation type(s): 1: u,2: v, 3:H_s, 4:k
% % Iteration 1
% 
% % Iter1.Rval(Iter1.Obs.rmswch==1,ir)=Iter1.Obs.rms(Iter1.Obs.rmswch==1);
% % R1=Iter1.Rval(:,ir).^2.*eye(length(Iter1.Obs.var));
% %
% %
% % flg=Iter1.True.flag;
% % Cost.p_oe(ir,ff,1)=(Iter1.True.var(flg==ff)-Iter1.Prior.var(flg==ff))'...
% %     *inv(R1(flg==ff,flg==ff))*(Iter1.True.var(flg==ff)-Iter1.Prior.var(flg==ff));
% 
% for ir=1:4
% for ff=1:4
% %
% Iter1.Rval(Iter1.Obs.rmswch==1,ir)=Iter1.Obs.rms(Iter1.Obs.rmswch==1);
% R1=Iter1.Rval(:,ir).^2.*eye(length(Iter1.Obs.var));
% 
% flg=Iter1.True.flag;
% Cost.p_oe(ir,ff,1)=(Iter1.True.var(flg==ff)-Iter1.Prior.var(flg==ff))'...
%     *inv(R1(flg==ff,flg==ff))*(Iter1.True.var(flg==ff)-Iter1.Prior.var(flg==ff));
% 
% Cost.full_p_oe(ir,1)=(Iter1.True.var-Iter1.Prior.var)'...
%     *inv(R1)*(Iter1.True.var-Iter1.Prior.var);
% 
% % % % Iteration 2
% Iter2.Rval(Iter2.Obs.rmswch==1,ir)=Iter2.Obs.rms(Iter2.Obs.rmswch==1);
% R2=Iter2.Rval(:,ir).^2.*eye(length(Iter2.Obs.var));
% 
% flg=Iter2.True.flag;
% Cost.p_oe(ir,ff,2)=(Iter2.True.var(flg==ff)-Iter2.Prior.var(flg==ff))'...
%     *inv(R2(flg==ff,flg==ff))*(Iter2.True.var(flg==ff)-Iter2.Prior.var(flg==ff));
% 
% Cost.full_p_oe(ir,2)=(Iter2.True.var-Iter2.Prior.var)'...
%     *inv(R2)*(Iter2.True.var-Iter2.Prior.var);
% % % % Iteration 3
% Iter3.Rval(Iter3.Obs.rmswch==1,ir)=Iter3.Obs.rms(Iter3.Obs.rmswch==1);
% R3=Iter3.Rval(:,ir).^2.*eye(length(Iter3.Obs.var));
% 
% flg=Iter3.True.flag;
% Cost.p_oe(ir,ff,3)=(Iter3.True.var(flg==ff)-Iter3.Prior.var(flg==ff))'...
%     *inv(R3(flg==ff,flg==ff))*(Iter3.True.var(flg==ff)-Iter3.Prior.var(flg==ff));
% 
% Cost.full_p_oe(ir,3)=(Iter3.True.var-Iter3.Prior.var)'...
%     *inv(R3)*(Iter3.True.var-Iter3.Prior.var);
% 
% % % % Iteration 4
% Iter4.Rval(Iter4.Obs.rmswch==1,ir)=Iter4.Obs.rms(Iter4.Obs.rmswch==1);
% R4=Iter4.Rval(:,ir).^2.*eye(length(Iter4.Obs.var));
% 
% flg=Iter4.True.flag;
% Cost.p_oe(ir,ff,4)=(Iter4.True.var(flg==ff)-Iter4.Prior.var(flg==ff))'...
%     *inv(R4(flg==ff,flg==ff))*(Iter4.True.var(flg==ff)-Iter4.Prior.var(flg==ff));
% 
% Cost.full_p_oe(ir,4)=(Iter4.True.var-Iter4.Prior.var)'...
%     *inv(R4)*(Iter4.True.var-Iter4.Prior.var);
% end
% end
% 
% f913=figure(913);clf;
% set(gcf,'units','centimeters','position',[5 5 8 8])
% plot(Cost.full_p_oe(2,:)./Cost.full_p_oe(2,1),'-o')
% xlabel('Iteration #')
% ylabel('$$J_{obs}$$ (Normalized by Iter. 1)', Interpreter='latex')
% xticks([1,2,3,4]);set(gca,'fontsize',12)
%         export_fig(f913,'Realtest/2023/obserror2',Plt.type,Plt.res,'-p0.025');%,'-m2')
% 
% return
% % % % Some temporary plotting <<<<<<<<<<<<<<<<<<<<

% Load the data that was read in Kitchen_sink_assimilate

% load('sinkouts/SARMay21_k.mat'); %

% Below files were used for Paper 2
% Names.sinkout='sinkouts/birds_uv.mat'; %Paper 2, 06/2022
% Names.sinkout='sinkouts/SARJun3_k.mat'; % Paper 2, 06/2022
Names.sinkout='sinkouts/drifters_uvh.mat'; % Paper 2, 06/2022
% %
% Names.sinkout='sinkouts/multivar_iter4_uvhk.mat'; % Paper 2, 06/2022
% Names.sinkout='sinkouts/multivar_iter3_uvhk.mat'; % Paper 2, 06/2022
% Names.sinkout='sinkouts/multivar_iter2_uvhk.mat'; % Paper 2, 06/2022
% Names.sinkout='sinkouts/multivar_iter1_uvhk.mat'; % Paper 2, SAR June3 only 06/2022

% circle (small outside area) comparisons for Hs vs. k
% Names.sinkout='sinkouts/circle/SARJun3_c1_k.mat'; % Paper 2, 06/2022
% Names.sinkout='sinkouts/circle/drifter_c3_h.mat'; % Paper 2, 06/2022

load(Names.sinkout);
% load('sinkouts/birds_h.mat'); %  %another iter 4(combining it with depth soundings, different from OS2022)

% Save Prior_iter4 for calculation purposes, however original prior
% (best guess linear profile)  will be used for plotting transects using:

Prior_iter4=Prior;

%% Plots
%         break
%     f22=figure(22);clf
%                 set(gcf,'units','centimeters','position',[5 5 40 40])
%     set(gcf,'units','inches','position',[0 0 7.8 5],'Color',[1 1 1]);
%     Trans.dir=Transects{trtr};
swch.Hor=false; %default

% Define Transect coordinates
%     load('../Adjoint_density/SAR_trns_NChnl2.mat');
load('SAR_trns_NChnl_usace.mat');

TransH.xtrans=Poly.lon_rho;
TransH.ytrans=Poly.lat_rho;
clear Poly
swch.markJetA=true;                               %Dotted line for JettyA
%                     case('Ver')
%     Vertical Transect (requires inpaint_nans)
% TransV.xtrans(1)=-124.038; TransV.ytrans(1)=46.2645;
% TransV.xtrans(2)=-124.039; TransV.ytrans(2)=46.235;
% swch.channel       =    0;   %used for N and S channel naming

% (works better for usace 2019)
%                 TransV.xtrans(1)=-124.03; TransV.ytrans(1)=46.2645;
%                 TransV.xtrans(2)=-124.031; TransV.ytrans(2)=46.235;
%         swch.channel       =    0;   %used for N and S channel naming

swch.markJetA=false;                               %Dotted line for JettyA
%orgshoal (looks good)
TransS.xtrans(1)=-124.101; TransS.ytrans(1)=46.26;
TransS.xtrans(2)=-124.085; TransS.ytrans(2)=46.23;

%alternative shoal
%             TransS.xtrans(1)=-124.12; TransS.ytrans(1)=46.248;
%             TransS.xtrans(2)=-124.105; TransS.ytrans(2)=46.21;

%Cross Channel track (org)
%     TransV.xtrans(1)=-123.95; TransV.ytrans(1)=46.25;
%     TransV.xtrans(2)=-123.97; TransV.ytrans(2)=46.215;

%Cross Channel(inside) track (New for 2014)
TransV.xtrans(1)=-123.965; TransV.ytrans(1)=46.255;
TransV.xtrans(2)=-123.985; TransV.ytrans(2)=46.218;
swch.channel       =    1;   %used for N and S channel naming

%Cross Channel(inside) track (New for 2014)

% -124.038     46.266     "JETTY A" (kind of where "the mouth" is)

Trans.JetLat=46.266;
Trans.JetLon=-124.038;
[Trans.eastJet,Trans.northJet] = ll2UTM(Trans.JetLon,Trans.JetLat,23,'10T');
Trans.eastJet=Trans.eastJet./1000;                         %convert to km
Trans.northJet=Trans.northJet./1000;                       %convert to km

f22=figure(22);clf
%                 set(gcf,'units','centimeters','position',[5 5 40 40])
set(gcf,'units','inches','position',[0 0 7.8 5],'Color',[1 1 1]);
for ir=irir
    %                             ir=3;
    clear R

    % If this is not a twin test, use the existing rmse from the Observations.
    if(swch.twin==0)
        Rval(Obs.rmswch==1,ir)=Obs.rms(Obs.rmswch==1);
    end

    R=Rval(:,ir).^2.*eye(length(Obs.var));
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

    %             Cov.C_hh_a=Cov.C_hh-Cov.C_hu*Ci*Cov.C_hu'; %Posterior Uncertainty
    %Cov.C_uu*inov_l is for U innovations

    display (['sum of inovations=', num2str(sum(inov(:))) 'm for Rval =' num2str(Rval(1,ir))])
    %     return
    innovation_2d(:,:,ir)=reshape(inov,[size(coord.lat_rho)]);
    %Plotting

    %                 innovation_2d(ir,:,:)(coord.maskR==0)=NaN;
    %             subplot(4,1,ir)
    trutemp=True.h;                    pritemp=Prior.h;
    True.h( coord.maskR==0)=NaN;
    Prior.h( coord.maskR==0)=NaN;
    Prior_iter4.h( coord.maskR==0)=NaN;

end    % ir loop

subplot('Position',pos.true)

pcolor(usace2014.lon,usace2014.lat,usace2014.z)
% shading flat
caxis([0 35]); xlim([min(coord.lon_rho(:)) max(coord.lon_rho(:))]);
ylim([min(coord.lat_rho(:)) max(coord.lat_rho(:))]);
hcb=colorbar;
colormap(flipud(haxby(64))); shading flat
hold on;

%     plot(TransH.xtrans,TransH.ytrans,'r-','LineWidth',1.3)
%     scatter(TransV.xtrans(1),TransV.ytrans(1))
%         plot(TransV.xtrans,TransV.ytrans,'k-','LineWidth',1.3)

%     plot(TransS.xtrans,TransS.ytrans,'b-','LineWidth',1.3)
%org
plot(Obs.lon2,Obs.lat2,'k.','MarkerSize',3,'MarkerFaceColor','k');
%org
% tind=True.flag<=2 & Obs.rmswch==1;
% % plot(Obs.lon2,Obs.lat2,'k.','MarkerSize',3,'MarkerFaceColor','k');
% plot(Obs.lon2(tind),Obs.lat2(tind),'r.','MarkerSize',3,'MarkerFaceColor','r');
% tind=True.flag<=2 & Obs.rmswch==0;
% plot(Obs.lon2(tind),Obs.lat2(tind),'g.','MarkerSize',3,'MarkerFaceColor','r');
% tind=True.flag=4 & Obs.rmswch==0;

pl_o = plot(nan, nan, 'k.', 'MarkerSize', 8, 'DisplayName', 'Observations');
legend(pl_o,'FontSize',6.25,'Location','southwest');
% ll1.Position = [0.2935 0.8860 0.0715 0.0292];

T1=text(min(coord.lon_rho(:))+0.01,max(coord.lat_rho(:))-0.02,'a)');xtickangle(0);
T1.BackgroundColor='w';T1.FontWeight='bold';T1.FontSize=8;xtickformat('degrees');ytickformat('degrees');
%                 for pp=1:length(Obs.lon2)
%                     text(Obs.lon2(pp),Obs.lat2(pp),num2str(pp));%datestr(Obs.timegmt(pp)))
%                 end
title(['USACE Survey'],'FontSize',Plt.fs3)
%         load(Names.polygon)
%         plot(Poly.lon_rho,Poly.lat_rho)
hold off

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
%             return
TransV.xtrans=linspace(TransV.xtrans(1),TransV.xtrans(2),50);
TransV.ytrans=linspace(TransV.ytrans(1),TransV.ytrans(2),50);

TransS.xtrans=linspace(TransS.xtrans(1),TransS.xtrans(2),50);
TransS.ytrans=linspace(TransS.ytrans(1),TransS.ytrans(2),50);

[TransH.eastT,TransH.northT] = ll2utm(TransH.ytrans,TransH.xtrans);
TransH.eastT=TransH.eastT./1000;                                 %convert to km
TransH.northT=TransH.northT./1000;                                %convert to km

[TransV.eastT,TransV.northT] = ll2utm(TransV.ytrans,TransV.xtrans);
TransV.eastT=TransV.eastT./1000;                                 %convert to km
TransV.northT=TransV.northT./1000;                                %convert to km

[TransS.eastT,TransS.northT] = ll2utm(TransS.ytrans,TransS.xtrans);
TransS.eastT=TransS.eastT./1000;                                 %convert to km
TransS.northT=TransS.northT./1000;                                %convert to km

load('sinkouts/smooth_prior.mat'); %original prior

for ir=irir
    F = inpaint_nans(usace2014.z);

    TransV.hfi(ir,:)=griddata( coord.lon_rho, coord.lat_rho, Prior.h,TransV.xtrans,...
        TransV.ytrans);
    TransV.hai(ir,:)=griddata( coord.lon_rho, coord.lat_rho,( pritemp+innovation_2d(:,:,ir)),...
        TransV.xtrans,TransV.ytrans);
    %         TransV.hti(ir,:)=griddata( coord.lon_rho, coord.lat_rho, trutemp,TransV.xtrans,...
    %             TransV.ytrans);
    TransV.hui(ir,:)=griddata( usace2019.lon, usace2019.lat, usace2019.z,TransV.xtrans,...
        TransV.ytrans);
    TransV.hu2i(ir,:)=griddata( usace2014.lon, usace2014.lat, F,TransV.xtrans,...
        TransV.ytrans);

    % Iter1 and Iter 2 for paper figure
    TransV.hP3i(ir,:)=griddata( coord.lon_rho, coord.lat_rho, Prior_iter4.h,TransV.xtrans,...
        TransV.ytrans);
    %     TransV.hP2i(ir,:)=griddata( coord.lon_rho, coord.lat_rho, Iter3.Prior.h,TransV.xtrans,...
    %         TransV.ytrans);
    %     TransV.hP1i(ir,:)=griddata( coord.lon_rho, coord.lat_rho, Iter2.Prior.h,TransV.xtrans,...
    %         TransV.ytrans);

    TransH.hfi(ir,:)=griddata( coord.lon_rho, coord.lat_rho, Prior.h,TransH.xtrans,...
        TransH.ytrans);
    TransH.hai(ir,:)=griddata( coord.lon_rho, coord.lat_rho,( pritemp+innovation_2d(:,:,ir)),...
        TransH.xtrans,TransH.ytrans);
    %         TransH.hti(ir,:)=griddata( coord.lon_rho, coord.lat_rho, trutemp,TransH.xtrans,...
    %             TransH.ytrans);

    %         TransH.hui(ir,:)=griddata( usace2019.lon, usace2019.lat, usace2019.z,TransH.xtrans,...
    %             TransH.ytrans);

    TransH.hu2i(ir,:)=griddata( usace2014.lon, usace2014.lat, F,TransH.xtrans,...
        TransH.ytrans);

    % Iter1 and Iter 2 for paper figure
    TransH.hP3i(ir,:)=griddata( coord.lon_rho, coord.lat_rho, Prior_iter4.h,TransH.xtrans,...
        TransH.ytrans);
    %     TransH.hP2i(ir,:)=griddata( coord.lon_rho, coord.lat_rho, Iter3.Prior.h,TransH.xtrans,...
    %         TransH.ytrans);
    %     TransH.hP1i(ir,:)=griddata( coord.lon_rho, coord.lat_rho, Iter2.Prior.h,TransH.xtrans,...
    %         TransH.ytrans);

    TransS.hfi(ir,:)=griddata( coord.lon_rho, coord.lat_rho, Prior.h,TransS.xtrans,...
        TransS.ytrans);
    TransS.hai(ir,:)=griddata( coord.lon_rho, coord.lat_rho,( pritemp+innovation_2d(:,:,ir)),...
        TransS.xtrans,TransS.ytrans);
    TransS.hui(ir,:)=griddata( usace2019.lon, usace2019.lat, usace2019.z,TransS.xtrans,...
        TransS.ytrans);
    TransS.hu2i(ir,:)=griddata( usace2014.lon, usace2014.lat, usace2014.z,TransS.xtrans,...
        TransS.ytrans);

    % Iter1 and Iter 2 for paper figure
    TransS.hP3i(ir,:)=griddata( coord.lon_rho, coord.lat_rho, Prior_iter4.h,TransS.xtrans,...
        TransS.ytrans);
    %     TransS.hP2i(ir,:)=griddata( coord.lon_rho, coord.lat_rho, Iter3.Prior.h,TransS.xtrans,...
    %         TransS.ytrans);
    %     TransS.hP1i(ir,:)=griddata( coord.lon_rho, coord.lat_rho, Iter2.Prior.h,TransS.xtrans,...
    %         TransS.ytrans);

    subplot('Position',pos.ver);cla
    ltrans=sqrt((TransV.eastT-TransV.eastT(1)).^2+(TransV.northT-...
        TransV.northT(1)).^2);
    plot(ltrans,TransV.hfi(ir,:),'r','linewidth',1.5);hold on;
    plot(ltrans,TransV.hai(ir,:),'b','linewidth',1.5) %usace

    %         plot(ltrans,TransV.hti(ir,:),'k-x','linewidth',1.5) %true model
    %         plot(ltrans,TransV.hui(ir,:),'k','linewidth',1.5) %usace

    plot(ltrans,TransV.hu2i(ir,:),'k','linewidth',1.5) %usace (2014 is truth for Paper2)
    %     plot(ltrans,TransV.hP3i(ir,:),'k--','linewidth',0.75)
    %     plot(ltrans,TransV.hP2i(ir,:),'k--','linewidth',0.75) %usace
    %     plot(ltrans,TransV.hP1i(ir,:),'k--','linewidth',0.75) %usace
    %         xlabel('Distance along the shoal transect [km]');xlim([0 2.5]);
    %                 xlabel('Distance along the shoal transect [km]');xlim([0 3.5]);

    xlabel('Distance along the cross-channel transect [km]');

    %                title('Cross-channel Transect','FontSize',Plt.fs2)


    if(swch.channel)
        T1=text(0.1,19,'c)');T1.BackgroundColor='w';T1.FontWeight='bold';T1.FontSize=8;
        ylim([5 20]); xlim([0 4.5]);
        T3=text(1,6,'North Channel','HorizontalAlignment','center');T3.FontSize=8;
        T4=text(3,6,'South Channel','HorizontalAlignment','center');T4.FontSize=8;
    else
        T1=text(0.1,37,'c)');T1.BackgroundColor='w';T1.FontWeight='bold';T1.FontSize=8;
        ylim([0 40]);xlim([0 3.5]);
        patch([1.75,1.75,3.5,3.5],[0,40,40,0],[0.5,0.5,0.5],'FaceAlpha',0.25,'LineStyle','none');
    end

    ll=legend('Prior','Posterior','USACE');ll.Location='north';ll.Box='off';
    box on;ylabel('Depth [m]');
    hold off; %set(gca,'fontsize',12);
    subplot('Position',pos.hor);cla
    %vertical line at JettyA
    Trans.yJet = [0:40];
    Trans.xJet=(Trans.eastJet-TransH.eastT(1))*ones(1,length(Trans.yJet));

    plot(TransH.eastT-TransH.eastT(1),TransH.hfi(ir,:),'r','linewidth',1.5);hold on;
    plot(TransH.eastT-TransH.eastT(1),TransH.hai(ir,:),'b','linewidth',1.5) %usace (iter3 for Paper2)

    %         plot(TransH.eastT-TransH.eastT(1),TransH.hti(ir,:),'k','linewidth',1.5)    %true model grid
    %         plot(TransH.eastT-TransH.eastT(1),TransH.hui(ir,:),'m','linewidth',1.5)
    plot(TransH.eastT-TransH.eastT(1),TransH.hu2i(ir,:),'k','linewidth',1.5) %usace (2014 is truth for Paper2)

    %     plot(TransH.eastT-TransH.eastT(1),TransH.hP3i(ir,:),'k--','linewidth',0.75)
    %     plot(TransH.eastT-TransH.eastT(1),TransH.hP2i(ir,:),'k--','linewidth',0.75) %usace (iter2 for Paper2)
    %     plot(TransH.eastT-TransH.eastT(1),TransH.hP1i(ir,:),'k--','linewidth',0.75) %usace (iter1 for Paper2)

    plot(Trans.xJet,Trans.yJet,'k.','linewidth',1.5)
    text(Trans.xJet(1)-0.25,Trans.yJet(1)+0.5,'Jetty A','FontWeight','bold','Rotation',90,'FontSize',6)
    %                 text(Trans.xJet(1)-0.25,Trans.yJet(end)-4,'Jetty A','FontWeight','bold','Rotation',90,'FontSize',6)

    %                         set(gcf,'paperposition',[0 0 8 3])
    xlabel('Distance along the channel transect [km]');
    %         legend('Prior','Posterior','USACE');
    box on;ylabel('Depth [m]');xlim([0 13]);
    %                 title('Along-channel Transect','FontSize',Plt.fs2)
    if(swch.panel5)
        T1=text(0.25,36,'e)');T1.BackgroundColor='w';T1.FontWeight='bold';T1.FontSize=8;
    else
        T1=text(0.25,36,'d)');T1.BackgroundColor='w';T1.FontWeight='bold';T1.FontSize=8;
    end
    hold off;%set(gca,'fontsize',12);

    if(swch.panel5)
        subplot('Position',pos.shoal);cla
        ltrans=sqrt((TransS.eastT-TransS.eastT(1)).^2+(TransS.northT-...
            TransS.northT(1)).^2);
        plot(ltrans,TransS.hfi(ir,:),'r','linewidth',1.5);hold on;
        plot(ltrans,TransS.hai(ir,:),'b','linewidth',1.5) %%usace (iter3 for Paper2)

        %         plot(ltrans,TransS.hui(ir,:),'k','linewidth',1.5) %usace
        plot(ltrans,TransS.hu2i(ir,:),'k','linewidth',1.5) %%usace (2014 is truth for Paper2)

        %         plot(ltrans,TransS.hP3i(ir,:),'k--','linewidth',0.75)
        %         plot(ltrans,TransS.hP2i(ir,:),'k--','linewidth',0.75) %%usace (iter2 for Paper2)
        %         plot(ltrans,TransS.hP1i(ir,:),'k--','linewidth',0.75) %%usace (iter1 for Paper2)

        xlabel('Outside the mouth transect [km]');ylabel('Depth [m]');
        box on; %ylabel('Depth [m]');
        %             set(gca, 'YAxisLocation', 'right')
        %                title('Cross-channel Transect','FontSize',Plt.fs2)
        %         T1=text(0.1,22.5,'c)');T1.BackgroundColor='w';T1.FontWeight='bold';T1.FontSize=8;
        T1=text(0.075,27.5,'d)');T1.BackgroundColor='w';T1.FontWeight='bold';T1.FontSize=8;
        ylim([10 30]);xlim([0 2.5]);
        hold off; %set(gca,'fontsize',12);
    end

    Prior.h( coord.maskR==0)=NaN;

    subplot('Position',pos.post);cla;
    posterior=Prior_iter4.h+innovation_2d(:,:,ir);
    %         posterior=Prior.h+innovation_2d(:,:,ir);
    posterior(posterior<2)=2;
    pcolor( coord.lon_rho', coord.lat_rho',posterior');hold on;
    pl(1)=plot(TransH.xtrans,TransH.ytrans,'k--','LineWidth',1.3);
    pl(2)=plot(TransV.xtrans,TransV.ytrans,'k-','LineWidth',1.3);
    pl(3)=plot(TransS.xtrans,TransS.ytrans,'k:','LineWidth',1.3);

    plot(Trans.JetLon,Trans.JetLat,'ko','MarkerSize',Plt.ms+1,'MarkerFaceColor','k')
    T2=text(Trans.JetLon,Trans.JetLat+0.01,'Jetty A');T2.FontSize=8;
    caxis([0 35]);hcb=colorbar;colormap(flipud(haxby(64)));shading interp
    title(['Posterior Bathymetry'],'FontSize',Plt.fs3)
    legend(pl(1:end),'Along Channel','Cross Channel','Out of the MCR','Location','southwest','FontSize',6.25)
    T1=text(min(coord.lon_rho(:))+0.01,max(coord.lat_rho(:))-0.02,'b)');xtickangle(0);
    T1.BackgroundColor='w';T1.FontWeight='bold';T1.FontSize=8;xtickformat('degrees');ytickformat('degrees');
    if(swch.Lc)
        fig3title.String= [' Assimilation using L_c = ',num2str(Lc(lll)*111.11*1000,4),...
            'm'];
        %                     fig3title.String=['U velocity at MGT3 = ' num2str(Sta.MGT3.u(iT)) ...
        %                         ', Sens U,V = ' num2str(Sens)];
        drawnow
    end

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
        fn=split(Names.sinkout,{'/','_'});
        export_fig(f22,['Realtest/2023/' fn{2} '_' fn{3} '_' num2str(ir)],Plt.type,Plt.res,'-p0.025');%,'-m2')
        %         export_fig(f22,['Realtest/drifters_only' num2str(ir)],Plt.type,Plt.res,'-p0.025');%,'-m2')
        %         export_fig(f22,['Realtest/2023/' fn{2} '_' num2str(ir)],Plt.type,Plt.res,'-p0.025');%,'-m2')

        %  <<<<<<<<<<<<<<    %org outputs

    end
end %ir loop, we want to plot each error level
clear Trans           %necessary for the transect loop
%                     pause
% end % End of the transect loop
%     pause

%% STATISTICAL COMPARISONS

% Uses the original Prior for comparisons

% I decided to divide the MCR into three regions: Outside, Mouth, and
% Inside. I make lot of statements in Paper #2 regarding the spatial so it
% makes sense looking at corrections for these different regions. I think
% the difference between the mouth and the inside of the domain occur due
% to the current profile difference one is like a channel one is a more
% complicated (meandering) domain.

% I also did a quick "Circle" comparison for a small region outside the
% domain.

if(swch.dostat)
    Names.polyarea = 'ThreeAreas_v2.mat';
    Locs={'Outside','Mouth','Inside'};

    % RMSE for area
    locs=load(Names.polyarea);

%     %
%         f35=figure(35); clf;
%         set(gcf,'units','centimeters','position',[5 5 12 3]); %paperposition is different
% 
% %         subplot(211)
% %         pcolor(coord.lon_rho',coord.lat_rho',usace2014.h'); hold on;
%         pcolor(usace2014.lon,usace2014.lat,usace2014.z); hold on;
% 
%         caxis([0 35]); colormap(flipud(haxby(64)));shading flat;
%         hcb=colorbar;hcb.Label.String='Depth [m]';hcb.Label.FontWeight='bold';
%         
%         pp(1)=plot(locs.Outside.lon,locs.Outside.lat,'r--',LineWidth=1.25);
%         pp(2)=plot(locs.Mouth.lon,locs.Mouth.lat,'k-',LineWidth=1.25);
%         pp(3)=plot(locs.Inside.lon,locs.Inside.lat,'b-.',LineWidth=1.25);
%         yticks([46.17:0.02:46.28])
%         xticks([-124.18:0.05:-123.86])
%         xtickformat('degrees');ytickformat('degrees');
%         legend(pp(1:end),'Outside','Mouth','Inside','Location','south','FontSize',5)
%         set(gca,'fontsize',6);
% 
%          export_fig(f35,'Realtest/2023/threearea',Plt.type,Plt.res,'-p0.025');
% 
%         subplot(212)
% %             pcolor(coord.lon_rho',coord.lat_rho',Prior.h'); hold on;
%         pcolor(coord.lon_rho',coord.lat_rho',posterior'); hold on;
%     
%         caxis([0 35]); colormap(flipud(haxby(64)));hcb=colorbar;shading flat
%         plot(locs.Outside.lon,locs.Outside.lat,'b-','MarkerSize',3,'MarkerFaceColor','k')
%         plot(locs.Mouth.lon,locs.Mouth.lat,'r-','MarkerSize',3,'MarkerFaceColor','k')
%         plot(locs.Inside.lon,locs.Inside.lat,'k-','MarkerSize',3,'MarkerFaceColor','k')

    fmt=['%5.2f   %5.2f   %5.2f %8.2f %5.2f %8.2f %4.0f \n' ...
        '%5.2f   %5.2f   %5.2f %8.2f %5.2f %8.2f %4.0f  \n'];
    fmt2=['%5.2f   %5.2f   %5.2f %8.2f %5.2f %8.2f %4.0f \n'];

    tmpname=split(Names.sinkout,{'/','_'}); %just a fancy way to name output

    for ll=1:3
        for ir=irir
            h_inov1= Prior_iter4.h+innovation_2d(:,:,ir);
            %             h_inov= Prior.h+innovation_2d(:,:,ir);
            if(swch.naninterp)
                usace2014.h=griddata(usace2014.lon,usace2014.lat,F,coord.lon_rho,coord.lat_rho);
            else
                %         usace2014.h=griddata(usace2014.lon,usace2014.lat,usace2014.z,coord.lon_rho,coord.lat_rho);
                h_inov=griddata(coord.lon_rho,coord.lat_rho,h_inov1,usace2014.lon,usace2014.lat);
                Prior.h2=griddata(coord.lon_rho,coord.lat_rho,Prior.h,usace2014.lon,usace2014.lat);

            end

            if(exist('Circle')) %only to plot mouth circle plots

                Names.Nens = ['stats/circle/' tmpname{3} '_' tmpname{2} '_062022.txt'];

                IN = inpolygon(coord.lon_rho, coord.lat_rho,Circle.X,Circle.Y);

                Area.indall=find(IN==1);

                % Calculate Prior vs. Truth and Posterior vs. Truth
                A=Prior.h(Area.indall);
                B=h_inov(Area.indall);
                C=usace2014.h(Area.indall);

                %                 f35=figure(35); clf;
                %                 subplot(211)
                %                 pcolor(coord.lon_rho,coord.lat_rho,usace2014.h);hold on;
                %                 plot(coord.lon_rho(Area.indall),coord.lat_rho(Area.indall),'k.');
                %                 caxis([0 35]); colormap(flipud(haxby(64)));shading flat
                [stat1.r,stat1.d,stat1.b,stat1.cc,stat1.si,stat1.ms]=rmse(A,C);
                [stat2.r,stat2.d,stat2.b,stat2.cc,stat2.si,stat2.ms]=rmse(B,C);

                if(isempty(dir(Names.Nens)))
                    fid = fopen(Names.Nens,'w');
                    fprintf(fid, 'RMSE     CC      Bias     D       si       MS   IR  \n');
                    fprintf(fid, fmt,...
                        stat1.r,stat1.cc,stat1.b,stat1.d,stat1.si,stat1.ms,ir, ...
                        stat2.r,stat2.cc,stat2.b,stat2.d,stat2.si,stat2.ms,ir);
                    fclose(fid);

                else
                    fid = fopen(Names.Nens, 'a+');
                    fprintf(fid, fmt2,...
                        stat2.r,stat2.cc,stat2.b,stat2.d,stat2.si,stat2.ms,ir);
                    fclose(fid);

                end

            else
                %             Names.Nens = ['stats/' Locs{ll} '_amultivariate_062022.txt'];
                if(swch.naninterp)
                    Names.Nens = ['stats/2023/' Locs{ll} '_' tmpname{2} '_interp.txt'];
                else
                    Names.Nens = ['stats/2023/' Locs{ll} '_' tmpname{2} '_nointerp.txt'];
                end
                % Locations are created with cells
                [usace2014.lon2,usace2014.lat2]= meshgrid(usace2014.lon,usace2014.lat);
                IN = inpolygon(usace2014.lon2,usace2014.lat2,locs.(Locs{ll}).lon, locs.(Locs{ll}).lat);

                %Temporary for finding obs count
                IN2 = inpolygon(Obs.lon2,Obs.lat2,locs.(Locs{ll}).lon, locs.(Locs{ll}).lat);
                disp([ ' Obs in =  ' Locs{ll} num2str(sum(IN2))]);
                %Temporary for finding obs count

                Area.indall=find(IN==1);


                % Calculate Prior vs. Truth and Posterior vs. Truth
                A=Prior.h2(Area.indall);
                B=h_inov(Area.indall);
                C=usace2014.z(Area.indall);
                %
                %                 f35=figure(35); clf;
                %                 subplot(211)
                %                 pcolor(coord.lon_rho,coord.lat_rho,usace2014.h);hold on;
                %                 plot(coord.lon_rho(Area.indall),coord.lat_rho(Area.indall),'k.');
                %                 caxis([0 35]); colormap(flipud(haxby(64)));shading flat
                %                 return %

                [stat1.r,stat1.d,stat1.b,stat1.cc,stat1.si,stat1.ms]=rmse(A,C);
                [stat2.r,stat2.d,stat2.b,stat2.cc,stat2.si,stat2.ms]=rmse(B,C);

                if(isempty(dir(Names.Nens)))
                    fid = fopen(Names.Nens,'w');
                    fprintf(fid, 'RMSE     CC      Bias     D       si       MS   IR  \n');
                    fprintf(fid, fmt,...
                        stat1.r,stat1.cc,stat1.b,stat1.d,stat1.si,stat1.ms,ir, ...
                        stat2.r,stat2.cc,stat2.b,stat2.d,stat2.si,stat2.ms,ir);
                    fclose(fid);

                else
                    fid = fopen(Names.Nens, 'a+');
                    fprintf(fid, fmt2,...
                        stat2.r,stat2.cc,stat2.b,stat2.d,stat2.si,stat2.ms,ir);
                    fclose(fid);
                end
            end %end circle loop
        end  %ir loop


    end %location loop
    %             R=corrcoef(A,B);
    %             [R,stat.d,stat.b,stat.cc,stat.si,stat.ms]=rmse(A,B);
    %
    %             disp([ ' rmse =  ' num2str(R)]);
    %             disp(' ')
end