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

z = [0:1:10 12.5 15 17.5 20];% depths where we evaluate the light field, in meters.

% z = [0 15]; 
SW = 350; % Surface SW irradiance in W/m^2

depth_ind = 13; 

if load_file
    
    load('lightfile')
    
else
    disp('Making Synthetic MP surface')
    % Now create the random circle radii and locations.
    target_a = 0.3; % Target area, if no ponds overlap
    base_radius = 1; % Radius of unperturbued circular in meters
    
    % This is the number of ponds in each direction if they do not overlap.
    nponds = round(sqrt(target_a*Lx^2 / (pi*base_radius^2)));
    
    % The amount over which the radii and positions can vary
    rad_mult = 4*base_radius;
    lattice_position_mult = 4*Lx/nponds;
    
    %% Create the lattice of random circle ponds on the domain.
    x_lattice = linspace(1,Lx,nponds); % Evenly space ponds on the grid
    y_lattice = x_lattice;
    [y_centers, x_centers] = meshgrid(x_lattice,x_lattice);
    
    radii = base_radius + rad_mult*(rand(length(x_lattice)*length(y_lattice),1)-.5);
    % radii = base_radius + rad_mult*(randn(length(x_lattice)*length(y_lattice),1)-.5);
    radii(radii < 0) = 0;
    pos_circ = [x_centers(:) y_centers(:)] + ...
        lattice_position_mult*(rand(nponds^2,2)-.5);
    
    %% Place the randomly adjusted circles on the domain
    for i = 1:length(pos_circ)
        
        cp = (pond_X - pos_circ(i,1)).^2 ...
            + (pond_Y - pos_circ(i,2)).^2 <= radii(i).^2;
        
        edge_mp(cp) = 1;
        
    end
    
    disp(sum(edge_mp(:))/numel(edge_mp))
    
end
% Get the statistics of the pond surface
CC = bwconncomp(edge_mp);
stats = table2array(regionprops('table',CC,'Area','Perimeter'));

addpath('../Fig-2'); 
[n,r] = boxcount(edge_mp);
r = r * dx^2; 
bc_grad = -gradient(log10(n))./gradient(log10(r));
D_bc = mean(bc_grad);
D_bc_weighted = sum((n./sum(n)).*bc_grad);

fprintf('BC Dim is %d, Weighted is %d',D_bc,D_bc_weighted); 

%% Now take the resulting data and fit it.
if size(stats,1) > 5
    
    %Area and Perimeter of ponds.
    A_data = dx^2*stats(:,1);
    P_data = dx*stats(:,2);
    
    % Sort it for doing gradients
    [A_data,I] = sort(A_data);
    P_data = P_data(I);
    P_data(P_data == 0) = 1;
        
else
    
    error('No ponds')
    
end

%% Now fit the data
% Logarithmically spaced to do line fitting.
x_data=log10(A_data);
y_data=log10(P_data);


% Make logarithmically spaced data to handle fitting
x_plot = linspace(min(log10(A_data)),max(log10(A_data)),200);

% Fit the logarithmically spaced data to the model P = k A^(D/2)
% Without assuming there is a bend
[D_0,intercept] = polyfit(x_data,y_data,1);

% Fit the a-p data to a tanh model as suggested by Bowen
fit_fun=fittype('x*(a1*tanh(a2*(x-a3))+a4)+c');
[f,gof,output] = fit(x_data(2:numel(x_data)),y_data(2:numel(x_data)),fit_fun,'StartPoint',[0.25,1,6,1.25,1],'Lower',[0,-Inf,-Inf,1,-Inf],'Upper',[0.5,+Inf,+Inf,1.5,+Inf]);

% Perimeter values obtained from the tanh fit.
logP_Brady = feval(f,x_plot);
logP_Line = polyval(D_0(2),x_plot);

Dstar = 2*log10(sum(P_data))/log10(sum(A_data));


% Two fittings with D = 1 and D = 2 and the total D

% Perimeter values obtained in each case.
P_Obs = P_data;
P_Fit = 10.^(logP_Brady);
P_Line = 10.^(logP_Line);

% Now compute the fractal dimensions
D_Line = 0*A_data + D_0(1); % Linear fit
D_Brady = f.a1*tanh(f.a2*(x_plot-f.a3))+f.a4;
D_Gradient = gradient(logP_Brady,x_plot);
D_data = 2*y_data./x_data; 


%% Next Compute the light field


if ~exist('light','var')
    

    disp('Computing Light Field')
    light = SW*comp_light_field_CH(pond_X,pond_Y,z,edge_mp);
    
end

%% Get some statistics of the light field
[I_mp,I_ic,I_net] = deal(zeros(length(z),1));  

for i = 1:length(z)
    
    temp = light(:,:,i);
    I_mp(i) = mean(mean(temp(logical(edge_mp))));
    I_ic(i) = mean(mean(temp(logical(1-edge_mp))));
    I_net(i) = mean(mean(temp)); 
    
end

I_net_pred = I_net(1)*exp(-z*.1);
I_mp_pred = I_mp(1)*exp(-z*.1);
I_ic_pred = I_ic(1)*exp(-z*.1);


%% Now plot things

% Plot the pond surface and the light field at depth
subplot(221)
imagesc([0 Lx],[0 Lx],interp2(light(:,:,1)',2));
hold on
title('Irradiance at Ice Base')
set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
grid on
box on
xlabel('m')
ylabel('m')
hold off
climmer = get(gca,'clim'); 
set(gca,'clim',[0 max(climmer)]); 

subplot(222)
imagesc([0 Lx],[0 Lx],interp2(light(:,:,depth_ind)',2));
title(sprintf('Depth %d m',z(depth_ind)));
set(gca,'yticklabel',[],'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
grid on
box on
xlabel('m')
% ylabel('km')

climmer = get(gca,'clim'); 

set(gca,'clim',[0 max(climmer)]); 

% Now do irradiance as a function of depth
clines = [228,26,28
55,126,184]/256; 

subplot(223)

plot((I_ic),-z,'linewidth',1,'color',clines(1,:))
hold on
plot((I_mp),-z,'linewidth',1,'color',clines(2,:));
plot(I_net,-z,'k','linewidth',1)


plot(I_ic_pred,-z,'--','color',clines(1,:)); 
plot(I_mp_pred,-z,'--','color',clines(2,:)); 
% plot(I_net_pred,-z,'--k')

hold off
legend('Under Ice','Under Pond','Net','location','southeast')
grid on
box on
set(gca,'ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
ylim([min(-z) 0]);
xlim([0 max(I_mp(:))]);
xlabel('W/m^2')
ylabel('Depth (m)')
title('Average Irradiance')

% Now include pond stats

subplot(224)

scatter(A_data,P_Obs,10)
hold on
box on
grid on

logP_1 = polyval([.5 .55],x_plot);
logP_2 =  polyval([1 -.5],x_plot);
logP_D = polyval([Dstar/2 0],x_plot);

% plot(A_plot,P_Fit,'--','linewidth',1)
% loglog(A_plot,P_Line,'-','linewidth',1)
loglog(10.^(x_plot),10.^(logP_1),'--','linewidth',1);
loglog(10.^(x_plot),10.^(logP_2),'--','linewidth',1);
loglog(10.^(x_plot),10.^(logP_D),'--','linewidth',1);
xlim([5 max(A_data)]);
ylim([5 5e3])
set(gca,'xscale','log','xtick',[1 10 100 1e3 1e4 1e5],'yscale','log','ydir','normal','layer','top','fontname','helvetica','Fontsize',9)
ylabel('P (m)'); xlabel('A (m^2)');
legend('Ponds','D=1','D=2','PDD','location','northwest')
title('Pond Geometry')
hold off

% Prep for plotting

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

pos = [6 3];


%
subplot(221)
cpos = get(gca,'position'); 
cpos(1) = cpos(1)+cpos(3) + .01; 
cpos(3) = .025; 
colorbar('position',cpos); 

subplot(222)
cpos = get(gca,'position'); 
cpos(1) = cpos(1)+cpos(3) + .01; 
cpos(3) = .025; 
colorbar('position',cpos); 

letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(e)','(c)'};

delete(findall(gcf,'Tag','legtag'))

for i = 1:4
    subplot(2,2,i)
    set(gca,'fontname','helvetica','fontsize',9,'xminortick','on','yminortick','on')
    posy = get(gca,'position');
    annotation('textbox',[posy(1)-.035 posy(2)+posy(4)+.05 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',10,'Tag','legtag');
    
end

%%
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

% saveas(gcf,'Fig-1.pdf')
% saveas(gcf,'Fig-1.fig')

