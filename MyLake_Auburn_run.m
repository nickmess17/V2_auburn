%function[SS_DO]=MyLake_Auburn_run(par)

path(path,'C:\Users\Nick\Desktop\MyLake_v2_Vansjo-master\MyLake-v2.0') %path for MyLake model code
path(path,'C:\Users\Nick\Desktop\MyLake_v2_Vansjo-master\Sediment-v2.0') %path for MyLake model code
path(path,'C:\Users\Nick\Desktop\MyLake_v2_Vansjo-master\Sediment-v2.0\ph-module') %path for MyLake model code

tic
disp('Started at:')
disp(datetime('now'));

is_metrics = true; % print metrics in the end

lake='Lake Auburn';
year=2012;
m_start=[2012,3,23]; %year,month,day
m_stop=[2013,12,31];

save_initial_conditions = false; % save final concentrations as initial for the next run

[lake_params, sediment_params] = load_params();
name_of_scenario = 'IO/Scenarios/Auburn_historical.txt';

file_name = 'IO/Auburn.mat';
%physical parameters
lake_params{2} = 0.0185; % 0.0442401  % 2     open water diffusion parameter (-)
lake_params{5} = 0.5; % 0.5 % 5     wind shelter parameter (-)
lake_params{16} = 1.8911; % 1 % 16    scaling factor for inflow volume (-)
lake_params{17} = 3.6; % 0 % 17    adjusting delta for inflow temperature (-)
lake_params{39} = 0.3398; % 2.5 % 39     non-PAR light attenuation coeff. (m-1)
lake_params{40} = 1.5623; %1.05 %40    PAR light attenuation coefficient for water (m-1)
%biological parameters
% lake_params{10} = 0.00005; %     PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
% lake_params{20} = 1.5; %     scaling factor for inflow concentration of total P (-)
% lake_params{56} = 0.1; %     Settling velocity for Chl a (m day-1)
% lake_params{57} = 0.1; %    Loss rate (1/day) at 20 deg C
% lake_params{58} = 1; %     Specific growth rate (1/day) at 20 deg C
% lake_params{53} = 0.5;%     Half saturation growth P level (mg/m3)
%oxygen parameters
% lake_params{70}=par(1);%Q10 for reactions of respiration
% lake_params{71}=par(2);%wc factor
% sediment_params{72}=par(3); % effective depth
% sediment_params{6}=par(4);%Km_O2
% sediment_params{13}=par(4);%Kin_O2

% try
run_ID = 0;
clim_ID = 0;
run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

[MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params , lake_params, name_of_scenario, use_INCA, run_INCA, run_ID, clim_ID, save_initial_conditions); % runs the model and outputs obs and sim

disp('Saving results...')
save(file_name, 'MyLake_results','Sediment_results')%'Sediment_results'
disp('Finished at:')
disp(datetime('now'));

%load observations
[obs_temp]=xlsread('C:\Users\Nick\Desktop\MyLake_v2_Vansjo-master\MyLake_temp_profiles.xlsx');%year,month,day,depth,temperature
[Obs_TP]=xlsread('C:\Users\Nick\Desktop\MyLake_v2_Vansjo-master\MyLake_TP_profiles.xlsx');
[obs_DO]=xlsread('C:\Users\Nick\Desktop\MyLake_v2_Vansjo-master\MyLake_DO_profiles.xlsx');

tt_mod=MyLake_results.basin1.days-datenum(year,1,1);%time now scaled so that it begins from the 1 january of the "year"

%Temperature profile observations(tt_mod, z, T)
date_of_observation=[obs_temp(:,1),obs_temp(:,2),obs_temp(:,3)];
TempObs=[datenum(date_of_observation) - datenum(year,1,1), obs_temp(:,4), obs_temp(:,5)];

%=align temperature observations with model results
initial=[1;find(diff(TempObs(:,1))~=0)+1];
finish=[find(diff(TempObs(:,1))~=0); length(TempObs)];
for i=1:length(initial)
operation=find(tt_mod==TempObs(initial(i),1));
    if (isempty(operation)==0)
    TempMod(initial(i):finish(i))=interp1(MyLake_results.basin1.z,MyLake_results.basin1.T(:,operation),TempObs(initial(i):finish(i),2));
    else
    TempMod(initial(i):finish(i))=NaN;    
    end    
end

%calculate sum squares for temperature
ss=zeros(1,numel(TempMod));
for i=1:numel(TempMod)
    ss(1,i)=((TempMod(1,i)-TempObs(i,3)).^2);
end
SS_temp=nansum(ss);
disp(SS_temp)

%linear regression on temperature
mdl=fitlm(TempObs(:,3),TempMod);
disp(mdl)

%TP profile observations(tt_mod, z, T)
date_of_observationP=[Obs_TP(:,1),Obs_TP(:,2),Obs_TP(:,3)];
TPobs=[datenum(date_of_observationP) - datenum(year,1,1), Obs_TP(:,4), Obs_TP(:,5)];
TP=MyLake_results.basin1.concentrations.PP+MyLake_results.basin1.concentrations.P+MyLake_results.basin1.concentrations.DOP;

%=align TP observations with model results
initialP=[1;find(diff(TPobs(:,1))~=0)+1];
finishP=[find(diff(TPobs(:,1))~=0); length(TPobs)];
for i=1:length(initialP)
operationP=find(tt_mod==TPobs(initialP(i),1));
    if (isempty(operationP)==0)
    TPMod(initialP(i):finishP(i))=interp1(MyLake_results.basin1.z,TP(:,operationP),TPobs(initialP(i):finishP(i),2));
    else
    TPMod(initialP(i):finishP(i))=NaN;    
    end    
end

%DO profile observations(tt_mod, z, T)
date_of_observationO=[obs_DO(:,1),obs_DO(:,2),obs_DO(:,3)];
DOobs=[datenum(date_of_observationO) - datenum(year,1,1), obs_DO(:,4), obs_DO(:,5)];
MyLake_results.basin1.concentrations.O2=MyLake_results.basin1.concentrations.O2./1000;%unit conversion

%=align DO observations with model results
initialO=[1;find(diff(DOobs(:,1))~=0)+1];
finishO=[find(diff(DOobs(:,1))~=0); length(DOobs)];
for i=1:length(initialO)
operationO=find(tt_mod==DOobs(initialO(i),1));
    if (isempty(operationO)==0)
    DOMod(initialO(i):finishO(i))=interp1(MyLake_results.basin1.z,MyLake_results.basin1.concentrations.O2(:,operationO),DOobs(initialO(i):finishO(i),2));
    else
    DOMod(initialO(i):finishO(i))=NaN;    
    end    
end

%calculate sum squares for DO
ss=zeros(1,numel(DOMod));
for i=1:numel(DOMod)
    ss(1,i)=((DOMod(1,i)-DOobs(i,3)).^2);
end
SS_DO=nansum(ss);
disp(SS_DO)

%linear regression on temperature
mdl=fitlm(DOobs(:,3),DOMod);
disp(mdl)

F_OM=1e+6*0.012;    %mass fraction [mg kg-1] of P of dry organic matter (assuming 50% of C, and Redfield ratio)

%prepare for plotting
zlim = [0 max(MyLake_results.basin1.z)];
tlim = [min(tt_mod) max(tt_mod)];

thermo = MyLake_results.basin1.MixStat(12,:);%thermocline depth

% figure(1)%temp vs depth profiles
% clf
% fign=9;
% post_operation=find(isnan(TempMod(initial))==0);
% post_initial=initial(post_operation);
% post_finish=finish(post_operation);
% for i = 1:min(fign,length(post_initial))
%    subplot(3,3,i);
%    operation=find(tt_mod==TempObs(post_initial(i),1));
%    plot(MyLake_results.basin1.T(:,operation),MyLake_results.basin1.z,'-b',TempObs(post_initial(i):post_finish(i),3),TempObs(post_initial(i):post_finish(i),2),'.r');
%    axis([0 30 zlim])
%    axis ij
%    title(['Lake Auburn ' datestr(tt_mod(operation)+datenum(year,1,1)),],'fontsize',8); 
%    set(gca,'FontSize',8) 
% end
% subplot(337)
% ylabel('Depth (m)','fontsize',8)
% xlabel('Temperature (^{o}C)','fontsize',8)
% subplot(331)
% ylabel('Depth (m)','fontsize',8)
% legend('Modelled','Observed');
% subplot(334)
% ylabel('Depth (m)','fontsize',8)
% subplot(338)
% xlabel('Temperature (^{o}C)','fontsize',8)
% subplot(339)
% xlabel('Temperature (^{o}C)','fontsize',8)
% 
% figure(2)%temp vs depth profiles
% clf
% for i = fign+1:min(18,length(post_initial))
%    subplot(3,3,i-fign);
%    set(gca,'FontSize',8) 
%    operation=find(tt_mod==TempObs(post_initial(i),1));
%    plot(MyLake_results.basin1.T(:,operation),MyLake_results.basin1.z,'-b',TempObs(post_initial(i):post_finish(i),3),TempObs(post_initial(i):post_finish(i),2),'.r');
%    axis([0 30 zlim])
%    axis ij
%    title(['Lake Auburn ' datestr(tt_mod(operation)+datenum(year,1,1)),],'fontsize',8); 
% end
%    subplot(337)
%    ylabel('Depth (m)','fontsize',8)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%     subplot(331)
%    ylabel('Depth (m)','fontsize',8)
%    legend('Modelled','Observed');
%     subplot(338)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%      subplot(339)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%    subplot(334)
% ylabel('Depth (m)','fontsize',8)
%    
% figure(3)%temp vs depth profiles
% clf
% for i = fign+10:min(27,length(post_initial))
%    subplot(3,3,i-(fign+9));
%    set(gca,'FontSize',8) 
%    operation=find(tt_mod==TempObs(post_initial(i),1));
%    plot(MyLake_results.basin1.T(:,operation),MyLake_results.basin1.z,'-b',TempObs(post_initial(i):post_finish(i),3),TempObs(post_initial(i):post_finish(i),2),'.r');
%    axis([0 30 zlim])
%    axis ij
%    title(['Lake Auburn ' datestr(tt_mod(operation)+datenum(year,1,1)),],'fontsize',8); 
% end
%    subplot(337)
%    ylabel('Depth (m)','fontsize',8)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%     subplot(331)
%    ylabel('Depth (m)','fontsize',8)
%    legend('Modelled','Observed');
%     subplot(338)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%      subplot(339)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%    subplot(334)
% ylabel('Depth (m)','fontsize',8)
%    
%    figure(4)%temp vs depth profiles
% clf
% for i = fign+19:min(36,length(post_initial))
%    subplot(3,3,i-(fign+18));
%    set(gca,'FontSize',8) 
%    operation=find(tt_mod==TempObs(post_initial(i),1));
%    plot(MyLake_results.basin1.T(:,operation),MyLake_results.basin1.z,'-b',TempObs(post_initial(i):post_finish(i),3),TempObs(post_initial(i):post_finish(i),2),'.r');
%    axis([0 30 zlim])
%    axis ij
%    title(['Lake Auburn ' datestr(tt_mod(operation)+datenum(year,1,1)),],'fontsize',8); 
% end
%    subplot(337)
%    ylabel('Depth (m)','fontsize',8)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%     subplot(331)
%    ylabel('Depth (m)','fontsize',8)
%    legend('Modelled','Observed');
%     subplot(338)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%      subplot(339)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%    subplot(334)
% ylabel('Depth (m)','fontsize',8)
% 
%    figure(5)%temp vs depth profiles
% clf
% for i = fign+28:min(45,length(post_initial))
%    subplot(3,3,i-(fign+27));
%    set(gca,'FontSize',8) 
%    operation=find(tt_mod==TempObs(post_initial(i),1));
%    plot(MyLake_results.basin1.T(:,operation),MyLake_results.basin1.z,'-b',TempObs(post_initial(i):post_finish(i),3),TempObs(post_initial(i):post_finish(i),2),'.r');
%    axis([0 30 zlim])
%    axis ij
%    title(['Lake Auburn ' datestr(tt_mod(operation)+datenum(year,1,1)),],'fontsize',8); 
% end
%    subplot(337)
%    ylabel('Depth (m)','fontsize',8)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%     subplot(331)
%    ylabel('Depth (m)','fontsize',8)
%    legend('Modelled','Observed');
%     subplot(338)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%      subplot(339)
%    xlabel('Temperature (^{o}C)','fontsize',8)
%    subplot(334)
% ylabel('Depth (m)','fontsize',8)
% 
% figure(6)%observed vs modelled scatter plot
% clf
% plot([0 25],[0 25],':b', TempObs(:,3), TempMod, '.r');
% set(gca,'fontsize',9);
% ylabel('Modelled temperature (^oC)');
% xlabel('Measured temperature (^oC)');
% axis([0 25 0 25]);
% axis square;
% title([lake ' ' datestr(datenum(m_start),28) '--' datestr(datenum(m_stop),28)]); 
% grid on;

%  tlims=[datenum([2013,1,1]) datenum([2014,12,31])];
% 
% figure(7)%time series
% clf
% subplot(511)
% inx=find(round(TempObs(:,2))==0);
% H=plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r.','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),MyLake_results.basin1.T(1,:),'b-');
%  set(gca,'ylim',[0 30]);
%  ylabel('^oC','fontsize',9)
%  title('Temperature  0m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xticklabel',[]);
% set(gca,'xlim',tlims)
%  
% subplot(512)
% inx=find((round(TempObs(:,2))==5));
% plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r.','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),MyLake_results.basin1.T(6,:),'b-');
%  set(gca,'ylim',[0 30]);
%  ylabel('^oC','fontsize',9)
%  title('Temperature  5m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xticklabel',[]);
%  set(gca,'xlim',tlims)
%  
% subplot(513)
% inx=find((round(TempObs(:,2))==9));
% plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r.','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),MyLake_results.basin1.T(10,:),'b-');
%  set(gca,'ylim',[0 30]);
%  ylabel('^oC','fontsize',9)
%  title('Temperature  9m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xticklabel',[]);
%  set(gca,'xlim',tlims)
%  
% subplot(514)
% inx=find((round(TempObs(:,2))==15));
% plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r.','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),MyLake_results.basin1.T(16,:),'b-');
%  set(gca,'ylim',[0 25]);
%  ylabel('^oC','fontsize',9)
%  title('Temperature  15m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xticklabel',[]);
%  set(gca,'xlim',tlims)
%  
%  subplot(515)
% inx=find((round(TempObs(:,2))==33));
% plot(TempObs(inx,1)+datenum(year,1,1),TempObs(inx,3),'r.','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),MyLake_results.basin1.T(34,:),'b-');
%  set(gca,'ylim',[0 25]);
%  ylabel('^oC','fontsize',9)
%  xlabel('year','fontsize',9)
%  title('Temperature 33m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xlim',tlims)
 
%  figure(8)%area vs depth
%  plot(MyLake_results.basin1.params.Az./1e+6,MyLake_results.basin1.z)
%  axis ij
%  xlabel('Horizontal Area x 10^6(km^2)')
%  ylabel('Depth (m)')
%  title('Bathymetry for Lake Auburn');
%  set(gca,'Ylim',[0 35])
%  grid on
%  
%  figure(9)%time vs depth temperature color plots
% clf
% pcolor(tt_mod,MyLake_results.basin1.z,MyLake_results.basin1.T)
% shading interp
% axis ij
% hold on
% plot(tt_mod,thermo,'k','LineWidth',1);
% hold off
% datetick('x','yy');
% set(gca,'ylim',zlim);
% caxis([0 24]);
% colorbar;
% set(gca,'fontsize',9);
% ylabel('Depth (m)')
% set(gca,'TickDir','out')
% grid on;
% 
% figure(10)%time series for TP
% clf
% subplot(211)
% inx=find(round(TPobs(:,2))<2);
% G=plot(TPobs(inx,1)+datenum(year,1,1),TPobs(inx,3),'r.','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),TP(2,:),'-');
%  set(gca,'ylim',[0 30]);
%  ylabel('ppb','fontsize',9)
%  xlabel('year','fontsize',9)
%  title('TP  1m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
% set(gca,'xlim',tlims)
% 
% subplot(212)
% inx=find((round(TPobs(:,2))==32)|(round(TPobs(:,2))==33));
% plot(TPobs(inx,1)+datenum(year,1,1),TPobs(inx,3),'r.','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),TP(34,:),'b-');
%  set(gca,'ylim',[0 60]);
%  ylabel('ppb','fontsize',9)
%  xlabel('year','fontsize',9)
% title('TP  33m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xlim',tlims)
%  
% figure(11)%time series for DO
% clf
% subplot(311)
% inx=find(round(DOobs(:,2))==1);
% I=plot(DOobs(inx,1)+datenum(year,1,1),DOobs(inx,3),'r.','linewidth',2);
% hold on
% plot(tt_mod+datenum(year,1,1),MyLake_results.basin1.concentrations.O2(2,:),'b-');
%  set(gca,'ylim',[0 20]);
%  ylabel('mg/L','fontsize',9)
%  title('DO  1m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xticklabel',[]);
% set(gca,'xlim',tlims)
% 
% subplot(312)
% inx=find((round(DOobs(:,2))==14));
% plot(DOobs(inx,1)+datenum(year,1,1),DOobs(inx,3),'r.','linewidth',2);
% hold on
% y=zeros(1,numel(tt_mod));
% y(1,:)=2;
% plot(tt_mod+datenum(year,1,1),y,"--");
% hold on
% plot(tt_mod+datenum(year,1,1),MyLake_results.basin1.concentrations.O2(15,:),'b-');
%  set(gca,'ylim',[0 20]);
%  ylabel('mg/L','fontsize',9)
%  title('DO 14m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xlim',tlims)
%  
% subplot(313)
% inx=find((round(DOobs(:,2))==29));
% plot(DOobs(inx,1)+datenum(year,1,1),DOobs(inx,3),'r.','linewidth',2);
% hold on
% y=zeros(1,numel(tt_mod));
% y(1,:)=2;
% plot(tt_mod+datenum(year,1,1),y,"--");
% hold on
% plot(tt_mod+datenum(year,1,1),MyLake_results.basin1.concentrations.O2(30,:),'b-');
%  set(gca,'ylim',[0 20]);
%  ylabel('mg/L','fontsize',9)
%  xlabel('year','fontsize',9)
%  title('DO 29m','fontweight','bold')
%  datetick('x','yy')
%  grid on
%  set(gca,'fontsize',9)
%  set(gca,'xlim',tlims)
 
%end