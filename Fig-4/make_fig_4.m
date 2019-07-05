
clear

addpath('../../../../../Matlab_Utilities/Plotting-Scripts/'); 
horvat_colors; 


load('../Fig-3/pseudo-pond-stats-100new.mat')

I_crit = 4.5;
k_w = .11;

% Model for transmission to ice base - based on values from Light (2008)
alpha_w = .15;
alpha_i = .75;
thick = 1;
k_i = 1; % Incl. SSL
k_p = 0.7; % Melt ponds
I_p = 0.75; 
I_b = 0.4; 

SW = 500/pi; % Radiance - for these computations. 


I_p = SW * I_p*exp(-thick*k_p)*(1-alpha_w); 
I_i = SW * I_b*exp(-thick*k_i)*(1-alpha_i); 


[mpfrac,Irr_net,Dstar,depth_i,depth_p,depth_net, ...
    tot_irr,tot_ice,tot_mp, ...
    abs_net,I0_p,I0_i] = ...
    deal(nan(length(pond_save),1));


depthfun = @(D,Irr) D*k_w - (Irr/I_crit)*(1-exp(1).^(-D*k_w)); 

fun = @(D) depthfun(D,I_p); 
depth_theo_p = fzero(fun,I_p/(k_w*I_crit));

fun = @(D) depthfun(D,I_i); 
depth_theo_i = fzero(fun,I_i/(k_w*I_crit));


%%

for i = 1:length(pond_save)
    
    dz = diff(z_save{i});
    
    edge_mp = pond_save{i};
    mpfrac(i) = sum(edge_mp(:))/numel(edge_mp);
    
    if (mpfrac(i) > 0.05)&&(mpfrac(i) < .95)
        
        % Compute the critical depth in general. 
        
        Irr_net = mpfrac(i)*I_p + (1-mpfrac(i))*I_i;
        
        fun = @(D) depthfun(D,Irr_net); 

        depth_net(i) = fzero(fun,depth_theo_p);
        
        % Compute the irradiance in ponded or unponded regions
        tot_irr(i) = sum(I_net_save{i}(2:end).*dz'); % Intergral of net I
        tot_ice(i) = sum(I_ic_save{i}(2:end).*dz'); % Under ice
        tot_mp(i) = sum(I_mp_save{i}(2:end).*dz'); % Under ponds
        
        abs_net(i) = Irr_net; % Total absorbed sunlight
        I0_i(i) = abs_net(i).*tot_ice(i)./tot_irr(i); %Fraction of irradiance in ice
        I0_p(i) = abs_net(i).*tot_mp(i)./tot_irr(i); % Ditto for ponds
        
        % Compute the critical depth for ponded regions incl. scattering. 
               
        
        fun = @(D) depthfun(D,I0_i(i)); 
        depth_i(i) = fzero(fun,10*depth_theo_i); 
     
        fun = @(D) depthfun(D,I0_p(i)); 
        
        depth_p(i) = fzero(fun,depth_theo_p); 
        
%% Compute fractal dimension
        [n,r] = boxcount(pond_save{i});
        
        n = n(1:end);
        r = r(1:end);
        
        dx = x_save{i}(2) - x_save{i}(1);
        r = r * dx^2;
        
        bc_grad = -gradient(log10(n))./gradient(log10(r));
        
        Dstar(i) = sum((n./sum(n)).*bc_grad);
        
    end
    
    
    
end

%% 
close all

pos = [6 2]; 
figure('windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
drawnow
cplots = [228,26,28
    55,126,184]/256;

change_depth_i = depth_i - depth_theo_i;
change_depth_p = depth_p - depth_theo_p;

% Ax{1} = subplot(131);
% 
% scatter(Dstar,depth_net,10,'filled','markerfacecolor',[0 0 0]);
% 
% box on
% grid on
% set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
% ylabel('m');
% titles{1} = 'Critical Depth';
% xlim([min(Dstar) max(Dstar)]);

Ax{1} = subplot(121);

scatter(Dstar,depth_i,10,'filled','markerfacecolor',cplots(1,:));
hold on
scatter(Dstar,depth_p,10,'filled','markerfacecolor',cplots(2,:));
% scatter(Dstar,depth_net,10,'filled','markerfacecolor','k');


xx = linspace(min(Dstar),max(Dstar),100);
plot(xx,0*xx + depth_theo_p,'--','color',cplots(2,:));


xx = linspace(min(Dstar),max(Dstar),100);
plot(xx,0*xx + depth_theo_i,'--','color',cplots(1,:));

legend('Ice','Ponds','location','northwest')


box on
grid on
set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
xlabel('PDD'); ylabel('Critical Depth','Fontsize',9);
xlim([min(Dstar) max(Dstar)]);
ylim([0 depth_theo_p*1.1])
title('Critical Depth','interpreter','latex');
    set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)

% Add CICE Runs
Ax{2} = subplot(122);
cla
LAT = double(ncread('std-cice-ponds.cice.h.ymonavg.var.1991-2015.nc','TLAT'));
LON = double(ncread('std-cice-ponds.cice.h.ymonavg.var.1991-2015.nc','TLON'));
% LON = 2*(LON + 90); 
% LAT(abs(LAT) > 1000) = nan; 
% LON(abs(LON) > 1000) = nan; 


std_area = ncread('std-cice-ponds.cice.h.ymonavg.var.1991-2015.nc','aice'); 
std_thick = ncread('std-cice-ponds.cice.h.ymonavg.var.1991-2015.nc','hi'); 
 
corr_area = ncread('std-cice-ponds-corr.cice.h.ymonavg.var.1991-2015.nc','aice'); 
corr_thick = ncread('std-cice-ponds-corr.cice.h.ymonavg.var.1991-2015.nc','hi'); 

corr_area_2 = ncread('std-cice-ponds-corr-2.cice.h.ymonavg.var.1991-2015.nc','aice'); 
corr_thick_2 = ncread('std-cice-ponds-corr-2.cice.h.ymonavg.var.1991-2015.nc','hi'); 


diff = corr_thick - std_thick;

axesm('mapprojection','ortho','maplatlimit',[70 90],'maplonlimit',[0,360],...
    'Frame','on','Grid','on','MLineLocation',90,'MeridianLabel','on',...
    'MLabelParallel',70,'ParallelLabel','On','PLabelMeridian',180,...
    'labelformat','compass','labelrotation','on','glinestyle',':','glinewidth',0.5,...
    'gcolor',[0.500 0.500 0.500],'grid','on');

set(gca,'visible','off')
tightmap

load('/Users/chorvat/Dropbox (Brown)/Research Projects/Published/Phytoplankton Blooms/Revision-2/Scripts/Fig-3/data_samegrid_bymonth','lat','lon','continents')
continents = double(continents);
contourm(lat,lon,continents,1,'k')
pcolorm(LAT,LON,diff(:,:,9));
contourm(LAT,LON,corr_area(:,:,9),[1 .15],'r'); 

set(gca,'clim',.05*[-1 1]); 
colormap(cmap_symm)
posser = get(gca,'position'); 
title('\Delta Ice Thickness','interpreter','latex')
tightmap

posy = get(gca,'position'); 
posy(2) = posy(2)+posy(4); 
posy(4) = .085; 
annotation('textbox',posy, ...
        'String','Thickness Change','Interpreter','latex','LineStyle','none','FontName','Helvetica', ...
        'Fontsize',9.9,'Tag','legtag','horizontalalignment','center');
    
colorbar('position',[posser(1) + posser(3)   posser(2)    .01    posser(4)])
% contourm(lat,lon,continents,1,'k')
    % surfm(lat,lon,plot)


%%


% Ax{2} = subplot(122);
% 
% 
% I = imread('melt-pond-light-sit.png'); 
% 
% I_sit = I(408:1074,805:1461,:); 
% 
% imshow(I_sit)
% title('Change in Thickness (m)','interpreter','latex','fontsize',9); 
% colorbar
% set(gca,'clim',[-.2 .2]); 
% colormap(cmap_symm)
% 
%     set(gca,'layer','top','fontname','helvetica','Fontsize',9)
% 


letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};


for i = 1:length(Ax)
    
    grid(Ax{i},'on');
    box(Ax{i},'on');
%    xlabel(Ax{i},'Fractal Dimension','Fontsize',9); 
%    title(Ax{i},titles{i},'interpreter','latex')
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1)-.05 posy(2)+posy(4)+.05 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'Fontsize',9,'Tag','legtag');
    
end


%%
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
saveas(gcf,'Fig-4.pdf')
saveas(gcf,'Fig-4.fig')


