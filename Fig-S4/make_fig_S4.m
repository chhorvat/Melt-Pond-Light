clear all

nx = 100; % Number of lateral grid points
dx = 5; % Grid spacing in meters

load(['../Fig-S3/pseudo-pond-stats-' num2str(nx) '-old.mat'])
%
SW_orig = 500;
SW = 350;
k_i_orig = 0.7; % original extinction coefficient
k_i = 0.7; % better extinction coefficient
thick = 1;

rel_I = (SW/SW_orig)*exp(-thick*k_i)/exp(-thick*k_i_orig);

[mpfrac,pfrac,abs_net,abs_net_theo,abs_ice_theo,abs_mp_theo, ...
    tot_irr,tot_ice,tot_mp,max_vals,max_depths,pondastat,pondrstat,pondpstat,pondnstat, ...
    D_1,D_bc,D_bc_weighted,pond_below] = ...
    deal(nan(length(A_surf_save),1));


kappa_w = .11;

a_threshold = (10)^2;

%%

%
for i = 1:length(A_pond_save)
    
    % Compute the melt pond fraction
    mpfrac(i) = A_surf_save(i) / (max(max(x_save{i}))^2);
    pfrac(i) = P_surf_save(i)/(kappa_w * max(max(x_save{i}))^2);
    
    % Compute the absolute forcing to the surface
    abs_net(i) = rel_I*I_net_save{i}(1); % Through ice
    abs_net_theo(i) = SW*(mpfrac(i)*0.4131 + (1-mpfrac(i))*.1241); % Theoretical through ice
    % abs_net and abs_net_theo should be the same (+/- a little due to radiation scheme)
    
    % Theoretical forcing to unponded-ponded regions
    abs_ice_theo(i) = rel_I*(1-mpfrac(i))*I_ic_save{i}(1);
    abs_mp_theo(i) = rel_I*mpfrac(i)*I_mp_save{i}(1);
    
    
    dz = diff(z_save{i});
    
    % The total irradiance summed over the depth - useful to compare the
    % fraction of irradiance disposed of in each category
    tot_irr(i) = rel_I*sum(I_net_save{i}(2:end).*dz');
    tot_ice(i) = rel_I*sum(I_ic_save{i}(2:end).*dz');
    tot_mp(i) = rel_I*sum(I_mp_save{i}(2:end).*dz');
    
    % Location and value of maximum irradiance under ice
    [max_vals(i),max_depths(i)] = max(smooth(I_ic_save{i}));
    max_depths(i) = z_save{i}(max_depths(i));
    
    
    %% Do fractal dimension calculation
    if (mpfrac(i) > 0.05).*(mpfrac(i) < .95)
        
        
        pond_below(i) = sum(A_pond_save{i}(A_pond_save{i} <= (a_threshold)))/(max(max(x_save{i}))^2);
        
        
        
    end
    
    pondastat(i) = mean(A_pond_save{i});
    pondpstat(i) = mean(P_pond_save{i});
    % Box-counting dimension
    % Remove all small areas because they bias the result.
    % pond_agg = bwareaopen(pond_save{i},60);
    % [n,r] = boxcount(pond_agg);
    [n,r] = boxcount(pond_save{i});
    
    n = n(1:end);
    r = r(1:end);
    
    dx = x_save{i}(2) - x_save{i}(1);
    r = r * dx^2;
    
    bc_grad = -gradient(log10(n))./gradient(log10(r));
    
    D_bc(i) = mean(bc_grad);
    
    D_bc_weighted(i) = sum((n./sum(n)).*bc_grad);
    
    
end
%% Do aggregate stats.

fracdim = D_bc_weighted;%  + 1 - min(D_bc_weighted(:));
% fracdim = D_bc;

abs_ice = (1-mpfrac).*abs_net.*tot_ice./tot_irr;
abs_mp = mpfrac.*abs_net.*tot_mp./tot_irr;

rel_mp = abs_mp ./ abs_net;
rel_ice = abs_ice ./ abs_net;

rel_mp_theo = abs_mp_theo./abs_net;
rel_ice_theo = abs_ice_theo./abs_net;

rel_change_mp = (rel_mp-rel_mp_theo)./rel_mp_theo;
rel_change_ice = (rel_ice-rel_ice_theo)./rel_ice_theo;


%%
close all

pos = [7 2.25];
figure('windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
drawnow


cplots = [228,26,28
    55,126,184]/256;

Ax{1} = subplot(121);
scatter(pond_below,fracdim,10,'filled','markerfacecolor','k');

hold on
box on
grid on
set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
xlabel('Small Pond Fraction','interpreter','latex');
ylabel('PDD'); 



Ax{2} = subplot(122);


scatter(pond_below,rel_change_mp,10,'filled','markerfacecolor',cplots(1,:));
hold on
scatter(pond_below,rel_change_ice,10,'filled','markerfacecolor',cplots(2,:));

hold on
box on
grid on
set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
xlabel('Small Pond Fraction','interpreter','latex');
ylabel('Enhancement','interpreter','latex')
legend('Ponds','Ice');
%%

letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};

for i = 1:length(Ax)
    
    posser = get(Ax{i},'position');
    %    xlabel(Ax{i},xlabs{i});
    %    set(Ax{i},'position',posser + [0 .03 0 -.02]);
    set(Ax{i},'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
    grid(Ax{i},'on');
    box(Ax{i},'on');
    %    title(Ax{i},titles{i})
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1)-.05 posy(2)+posy(4)+.05 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'Fontsize',9,'Tag','legtag');
    
end
% Plotting tools

%%
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');

saveas(gcf,'Fig-S4.pdf')
saveas(gcf,'Fig-S4.fig')

%%

