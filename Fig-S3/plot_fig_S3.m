% make_fig_2
%%
SW_orig = 500;
SW = 350;


rel_I = (SW/SW_orig);

[mpfrac,pfrac,abs_net,abs_net_theo,abs_ice_theo,abs_mp_theo, ...
    tot_irr,tot_ice,tot_mp,max_vals,max_depths] = ...
    deal(zeros(length(A_surf_save),1));

[pondastat,pondpstat,D_0,D_1,D_mean,D_horv,D_bc,D_bc_weighted, ...
    D_bc_weighted_2,D_bc_weighted_3] = deal(zeros(1,length(A_pond_save)));


%%
for i = 1:length(A_pond_save)
    
    % Compute the melt pond fraction
    mpfrac(i) = A_surf_save(i) / (max(max(x_save{i}))^2);
    pfrac(i) = P_surf_save(i) / (max(max(x_save{i}))^2);
    
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
    if (mpfrac(i) > 0.05)&&(mpfrac(i) < .95)
        
        %  inds = logical((A_pond_save{i} > 10).*(A_pond_save{i} < 1000));
        
        % Fit A-P info
        
        if length(A_pond_save{i}) > 1
            
            [fstat,stat_fit] = polyfit(log10(A_pond_save{i}),log10(P_pond_save{i}),1);
            
            D_0(i) = fstat(1);
            
        else
            
            D_0(i) = nan;
            
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
        
        
        D_bc_weighted_2(i) = sum((n./sum(n)).*log10(n)./log10(r));
        
        D_bc_weighted_3(i) = sum((n./sum(n)).*log10(n./sum(n))./log10(r));
        
        % -(n(end) - n(1))/sum(n);
        
        % Hausdorff Dimension
        
        D_1(i) = hausDim(pond_save{i});
        
        
        % Surface dimension
        D_horv(i) = 2*log10(sum(P_pond_save{i}))/log10(sum(A_pond_save{i}));
        
        
        % Mean dimension
        D_mean(i) = mean(D_pond_save{i});
        
        
        
    else
        
        pondastat(i) = nan;
        pondpstat(i) = nan;
        D_0(i) = nan;
        D_1(i) = nan;
        D_mean(i) = nan;
        D_horv(i) = nan;
        D_bc(i) = nan;
        D_bc_weighted(i) = nan;
        D_bc_weighted_2(i) = nan;
        D_bc_weighted_3(i) = nan;
        
    end
    
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

X = [ones(size(A_surf_save)); fracdim; mpfrac'; pfrac']';
[reg_coeff,~,~,~,~] = regress(rel_change_mp,X);

model_mp = sum(bsxfun(@times,reg_coeff,X'),1)';

X = [ones(size(A_surf_save)); fracdim]';
[reg_coeff_1,~,~,~,~] = regress(rel_change_mp,X);

model_mp_1 = sum(bsxfun(@times,reg_coeff_1,X'),1)';

madeup_coeff_1 = [-.6 .3]';

model_mp_madeup = sum(bsxfun(@times,madeup_coeff_1,X'),1)';


alpha_w = .15;
alpha_i = .75;

prefac = ((1-alpha_w)/(1-alpha_i));


model_ice = -model_mp.*(rel_mp_theo./rel_ice_theo);
model_ice_1 = -model_mp_1.*(rel_mp_theo./rel_ice_theo);

Qscat = abs_mp-abs_mp_theo;

thick = 1;
k_i = 0.9;

Qscat_model = SW*exp(-k_i*thick)*mpfrac.*model_mp;
Qscat_model_1 = SW*exp(-k_i*thick)*mpfrac.*model_mp_1;
Qscat_madeup = SW*exp(-k_i*thick)*mpfrac.*model_mp_madeup;


X = [ones(size(A_surf_save)); fracdim; mpfrac'; pfrac']';
[reg_coeff_2,~,~,~,~] = regress(rel_change_ice,X,.01);

frac_cutoff = .95;


Rsq = comp_rsq(rel_change_mp,model_mp);

fprintf('The full model for E_p fits with R^2 = %d \n',Rsq);

Rsq = comp_rsq(rel_change_mp,model_mp_1);

fprintf('The D model for E_p fits with R^2 = %d \n',Rsq);


% I busted something here!!!

Rsq = comp_rsq(rel_change_ice,model_ice);

fprintf('The full model for E_i fits with R^2 = %d \n',Rsq);

Rsq = comp_rsq(rel_change_ice(mpfrac < frac_cutoff),model_ice_1(mpfrac < frac_cutoff));

fprintf('The D model for E_i below mpfrac = %d fits with R^2 = %d \n',frac_cutoff,Rsq);

Rsq = comp_rsq(Qscat,Qscat_model);

fprintf('The model for Qscat  fits with R^2 = %d \n',Rsq);

Rsq = comp_rsq(Qscat,Qscat_model_1);

fprintf('The D model for Qscat  fits with R^2 = %d \n',Rsq);

Rsq = comp_rsq(Qscat,Qscat_madeup);

fprintf('The easy model for Qscat  fits with R^2 = %d \n',Rsq);
%%
close all

pos = [7 2.25];
figure('windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
drawnow


cplots = [228,26,28
    55,126,184]/256;


% Ax{1} = subplot(231);
%
% scatter(D_1,fracdim,10,'filled','markerfacecolor',cplots(2,:));
% hold on
% box on
% grid on
% set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
% xlab{1}='Hausdorff Dim'; ylabel('D');
% hold off
% ylim([0.5 2]); xlim([1 2])
% titles{1} = '';

Ax{1} = subplot(131);

scatter(fracdim,mpfrac,10,'filled','markerfacecolor','k');
hold on
box on
grid on
set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
ylabel('\phi'); xlabel('PDD');
% title('Pond Fractal Dimension')
hold off
% ylim([1 2]); xlim([0 1])
title('Pond Fraction','interpreter','latex');

Ax{2} = subplot(132);
%
drawnow
scatter(fracdim,rel_change_ice,10,'filled','markerfacecolor',cplots(1,:));
hold on
scatter(fracdim,rel_change_mp,10,'filled','markerfacecolor',cplots(2,:));
% scatter(fracdim,model_ice,'filled');
plot(fracdim,model_mp_1,'--k','linewidth',1)
% scatter(fracdim,model_ic)

box on
grid on
set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
ylabel('Multiple','interpreter','latex'); xlabel('PDD','interpreter','latex');
title('Enhancement Factor','interpreter','latex');
hold off
% legend('Ice','Ponds','Linear Model')

%
Ax{3} = subplot(133);

scatter(fracdim,(abs_ice-abs_ice_theo),10,'filled','markerfacecolor',cplots(1,:));
hold on
scatter(fracdim,(abs_mp-abs_mp_theo),10,'filled','markerfacecolor',cplots(2,:));
box on
grid on
set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
ylabel('Q_{scat} (W/m^2)'); xlabel('PDD','interpreter','latex');
title('Change in Solar Flux','interpreter','latex');
hold off


% plot(fracdim,abs_mp_theo.*(1+model_mp'),'--k','linewidth',1)
% scatter(fracdim,model_ic)


% subplot(133)
%
% scatter(fracdim,max_depths,'filled','markerfacecolor',cplots(2,:));
% hold on
% box on
% grid on
% set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
% ylabel('Depth (m)'); xlabel('Fractal Dimension');
% title('Location of sub-ice max')
% hold off
% xlim([1 2]);
% ylim([0 10]);
% xlim([1 2])

% subplot(134)
%
% scatter(mpfrac,rel_ice./rel_ice_theo);
% hold on
% box on
% grid on
% set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
% ylabel('Increase in I'); xlabel('Fractal Dimension');
% title('Increase in heat under ice')
% hold off
% xlim([1 2])

letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};

%%
for i = 1:length(Ax)
    
    posser = get(Ax{i},'position');
    %    xlabel(Ax{i},xlabs{i});
    %    set(Ax{i},'position',posser + [0 .1 0 -.1]);
    set(Ax{i},'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
    grid(Ax{i},'on');
    box(Ax{i},'on');
    %    title(Ax{i},titles{i})
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1)-.06 posy(2)+posy(4)+.06 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'Fontsize',9,'Tag','legtag');
    
end

%legend(Ax{3},'Ice','Ponds','position',[0.5307    0.0017    0.1280    0.1512],'orientation','horizontal')

% Plotting tools

%
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');

saveas(gcf,'Fig-S3.pdf')
saveas(gcf,'Fig-S3.fig')
