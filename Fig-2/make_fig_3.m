clear
close all

load_file = 0;

%% Initialize variables
nx = 500; % Number of lateral grid points
dx = 1; % Grid spacing in meters
Lx = dx*nx; % Size of domain in meters
edge_mp = zeros(nx,nx); % melt pond surface

pond_X = dx*(1:nx); % The grid
[pond_Y,pond_X] = meshgrid(pond_X,pond_X); % Grid mesh

z = [0:.5:10];% depths where we evaluate the light field, in meters.
SW = 350; % Surface SW irradiance in W/m^2

baserad = [.25 2.5 25];

%%

load lighties; 
load pondies; 

%%

for pond_i = 1:3
   
    edge_mp = ponder{pond_i};
    % Get the statistics of the pond surface
    % CC = bwconncomp(edge_mp);
    % stats = table2array(regionprops('table',CC,'Area','Perimeter'));
    
    [n,r] = boxcount(edge_mp);
    r = r * dx^2;
    bc_grad = -gradient(log10(n))./gradient(log10(r));
    D_bc(pond_i) = mean(bc_grad);
    D_bc_weighted(pond_i) = sum((n./sum(n)).*bc_grad);
    
    fprintf('BC Dim is %d, Weighted is %d \n',D_bc(pond_i),D_bc_weighted(pond_i));

end

%% 
close all

for pond_i = 1:3
    %% Now plot things
    
    edge_mp = ponder{pond_i};
    % Plot the pond surface and the light field at depth
    
    depth_ind = 13;
    
    Ax{pond_i} = subplot(2,3,pond_i);
    imagesc([0 Lx],[0 Lx],interp2(edge_mp',2));
    hold on
    title(sprintf('PDD = %.1f',D_bc_weighted(pond_i)));
    set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
    grid on
    box on
    xlabel('m')
    ylabel('m')
    hold off
    pbaspect([1 1 1]); 
    % climmer = get(gca,'clim');
    % set(gca,'clim',[0 max(climmer)]);
    
    Ax{pond_i+3} = subplot(2,3,3+pond_i);
    
    
    % Irradiance as a function of depth
    clines = [228,26,28
        55,126,184]/256;
    
    I_ic_plot = I_ic_pred{pond_i} + smooth(I_ic{pond_i} - I_ic_pred{pond_i}',3)'; 
    I_mp_plot = I_mp_pred{pond_i} + smooth(I_mp{pond_i} - I_mp_pred{pond_i}',5)'; 
    
    
    plot((I_ic_plot),-z,'linewidth',1,'color',clines(1,:))
    hold on
    plot((I_mp_plot),-z,'linewidth',1,'color',clines(2,:));
    plot(I_net{pond_i},-z,'k','linewidth',1)
    
    
    plot(I_ic_pred{pond_i},-z,'--','color',clines(1,:));
    plot(I_mp_pred{pond_i},-z,'--','color',clines(2,:));
    % plot(I_net_pred,-z,'--k')
    
    hold off
%    legend('Ice','Pond','Net','location','southeast')
    grid on
    box on
    set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
    ylim([min(-z) 0]);
    xlim([0 max(I_mp{pond_i}(:))]);
    xlabel('W/m^2')
    if pond_i == 1
    ylabel('Depth (m)')
    end
    title('Average Irradiance')
    
    
end

cmap = [256 256 256
    247,251,255
222,235,247
198,219,239
158,202,225
107,174,214
66,146,198
33,113,181
8,81,156
8,48,107]/256; 

colormap(cmap)

pos = [6 3.5];

set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

%

letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};

for i = 1:length(Ax)
    
  %  figure; subplot(2,3,i); pos = get(gca,'position'); close; set(Ax{i},'position',pos);
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


%
for i = 4:6
        posser = get(Ax{i},'position'); 
        set(Ax{i},'position',posser + [0 .1 0 -.1]); 
end

legend('Ice','Pond','Net','position',[.325 .02 .35 .05],'orientation','horizontal')




%%
% saveas(gcf,'Fig-2.pdf')
% saveas(gcf,'Fig-2.fig')

