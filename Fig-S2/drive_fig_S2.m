%drive_fig_2
clear
close all


kappa_vals = [.055 .11 .22]; % kappa value to use. Default is .11
dx_vals = [2.5 5 10];
nx = 100; % Number of lateral grid points

A_vals = .1 + .6*rand(25,1);
%
%

for ind_k = 1:length(kappa_vals)
    
    dx = 5;
    
    R_vals = round(logspace(log10(5),log10(nx*dx),25));

    L = nx*dx;
    
    num_it = 0;
    
    kappa_pert = kappa_vals(ind_k);
    
    outstr = ['L-' num2str(nx*dx) '-k-' num2str(round(100*kappa_pert)) '-dx-' num2str(round(dx))];
    
    for ind_R = 1:length(R_vals)
        
        for ind_A = 1:length(A_vals)
            
            
            try
                
                fprintf('Iter %d of %d \n',(ind_R-1)*length(A_vals)+ind_A,length(A_vals)*length(R_vals));
                clearvars -except Ax ind_* *_vals nx dx* *_save num_it kappa* outstr
                
                base_radius = R_vals(ind_R);
                rad_mult = R_vals(ind_R);
                target_a = A_vals(ind_A);
                
                process_fig_S2_nosave_sens;
                
            catch err_proc
                
                disp(err_proc.message);
                
            end
            
        end
        
    end
    
    %% Make fig panels
    SW_orig = 500;
    SW = 350;
    k_i_orig = 0.7; % original extinction coefficient
    k_i = 0.7; % better extinction coefficient
    thick = 1;
    
    rel_I = (SW/SW_orig)*exp(-thick*k_i)/exp(-thick*k_i_orig);
    
    [mpfrac,pfrac,abs_net,abs_net_theo,abs_ice_theo,abs_mp_theo, ...
        tot_irr,tot_ice,tot_mp,max_vals,max_depths] = ...
        deal(zeros(length(A_surf_save),1));
    
    %
    for i = 1:length(A_pond_save)
        
        % Compute the melt pond fraction
        mpfrac(i) = A_surf_save(i) / (max(max(x_save{i}))^2);
        pfrac(i) = P_surf_save(i) / (max(max(x_save{i})));
        
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
        
        %% Do fractal dimension calculation
        if (mpfrac(i) > 0.05).*(mpfrac(i) < .95)
            
            pondastat(i) = mean(A_pond_save{i});
            pondpstat(i) = mean(P_pond_save{i});
            
            [n,r] = boxcount(pond_save{i});
            
            n = n(1:end);
            r = r(1:end);
            
            dx = x_save{i}(2) - x_save{i}(1);
            r = r * dx^2;
            
            bc_grad = -gradient(log10(n))./gradient(log10(r));
            
            D_bc_weighted(i) = sum((n./sum(n)).*bc_grad);
            
            
            
        else
            
            pondastat(i) = nan;
            pondpstat(i) = nan;
            D_bc_weighted(i) = nan;
            
        end
        %
        
    end
    
    fracdim = D_bc_weighted;
    
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
    
    Rsq = comp_rsq(rel_change_mp,model_mp);
    fprintf('The full model for E_p fits with R^2 = %d \n',Rsq);
    Rsq = comp_rsq(rel_change_mp,model_mp_1);
    fprintf('The D model for E_p fits with R^2 = %d \n',Rsq);
    
    fprintf('The slope of the D model is %d \n',reg_coeff_1(2));
    
    
    %%
    %%
    
    
    cplots = [228,26,28
        55,126,184]/256;
    
    Ax{ind_k} = subplot(2,3,ind_k);
    scatter(fracdim,rel_change_mp,10,'filled','markerfacecolor',cplots(1,:));
    hold on
    scatter(fracdim,rel_change_ice,10,'filled','markerfacecolor',cplots(2,:));
    plot(fracdim,model_mp_1,'--k','linewidth',1)
    box on
    grid on
    set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
    xlabel('Pond D','interpreter','latex');
    ylabel('$E_p$','interpreter','latex')
    
    
end

%%

for ind_x = 1:length(dx_vals)
    
    dx = dx_vals(ind_x);
    kappa_pert = .11;
    
    
    L = nx*dx;
    
    num_it = 0;
    
    
    outstr = ['L-' num2str(nx*dx) '-k-' num2str(round(100*kappa_pert)) '-dx-' num2str(round(dx))];
    
    for ind_R = 1:length(R_vals)
        
        for ind_A = 1:length(A_vals)
            
            
            try
                
                fprintf('Iter %d of %d \n',(ind_R-1)*length(A_vals)+ind_A,length(A_vals)*length(R_vals));
                clearvars -except Ax ind_* *_vals nx dx* *_save num_it kappa* outstr
                
                base_radius = R_vals(ind_R);
                rad_mult = R_vals(ind_R);
                target_a = A_vals(ind_A);
                
                process_fig_S2_nosave_sens;
                
            catch err_proc
                
                disp(err_proc.message);
                
            end
            
        end
        
    end
    
    %% Make fig panels
    SW_orig = 500;
    SW = 350;
    k_i_orig = 0.7; % original extinction coefficient
    k_i = 0.7; % better extinction coefficient
    thick = 1;
    
    rel_I = (SW/SW_orig)*exp(-thick*k_i)/exp(-thick*k_i_orig);
    
    [mpfrac,pfrac,abs_net,abs_net_theo,abs_ice_theo,abs_mp_theo, ...
        tot_irr,tot_ice,tot_mp,max_vals,max_depths] = ...
        deal(zeros(length(A_surf_save),1));
    
    %
    for i = 1:length(A_pond_save)
        
        % Compute the melt pond fraction
        mpfrac(i) = A_surf_save(i) / (max(max(x_save{i}))^2);
        pfrac(i) = P_surf_save(i) / (max(max(x_save{i})));
        
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
        
        %% Do fractal dimension calculation
        if (mpfrac(i) > 0.05).*(mpfrac(i) < .95)
            
            pondastat(i) = mean(A_pond_save{i});
            pondpstat(i) = mean(P_pond_save{i});
            
            [n,r] = boxcount(pond_save{i});
            
            n = n(1:end);
            r = r(1:end);
            
            dx = x_save{i}(2) - x_save{i}(1);
            r = r * dx^2;
            
            bc_grad = -gradient(log10(n))./gradient(log10(r));
            
            D_bc_weighted(i) = sum((n./sum(n)).*bc_grad);
            
            
            
        else
            
            pondastat(i) = nan;
            pondpstat(i) = nan;
            D_bc_weighted(i) = nan;
            
        end
        %
        
    end
    
    fracdim = D_bc_weighted;
    
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
    
    Rsq = comp_rsq(rel_change_mp,model_mp);
    fprintf('The full model for E_p fits with R^2 = %d \n',Rsq);
    Rsq = comp_rsq(rel_change_mp,model_mp_1);
    fprintf('The D model for E_p fits with R^2 = %d \n',Rsq);
    
    fprintf('The slope of the D model is %d \n',reg_coeff_1(2));
    
    
    %%
    
    
    
    cplots = [228,26,28
        55,126,184]/256;
    
    Ax{ind_x+3} = subplot(2,3,ind_x+3);
    scatter(fracdim,rel_change_mp,10,'filled','markerfacecolor',cplots(1,:));
    hold on
    scatter(fracdim,rel_change_ice,10,'filled','markerfacecolor',cplots(2,:));
    plot(fracdim,model_mp_1,'--k','linewidth',1)
    box on
    grid on
    set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
    xlabel('PDD','interpreter','latex');
    ylabel('$E_p$','interpreter','latex')
    
    
end

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

%%
pos = [7 4.5];
% figure('windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
drawnow

saveas(gcf,'Fig-S2.pdf')
saveas(gcf,'Fig-S2.fig')

