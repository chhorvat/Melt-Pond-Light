%% Initialize Domain
Lx = dx*nx; % Size of domain in meters
edge_mp = zeros(nx,nx); % melt pond surface

pond_X = dx*(1:nx); % The grid
[pond_Y,pond_X] = meshgrid(pond_X,pond_X); % Grid mesh

z = [0:.5:5 10 15 20 25 30];
SW = 500; % Surface SW irradiance in W/m^2

%% Begin by generating a synthetic pond surface

% This is the number of ponds in each direction if they do not overlap.
nponds = round(sqrt(target_a*Lx^2 / (pi*base_radius^2)));

lattice_position_mult = Lx/nponds;


%% Create the lattice of random circle ponds on the domain.
x_lattice = linspace(1,Lx,nponds+2);
% Evenly space ponds on the grid. Don't want one at the start or end
x_lattice = x_lattice(2:end-1);
y_lattice = x_lattice;
[y_centers, x_centers] = meshgrid(x_lattice,x_lattice);

radii = base_radius + rad_mult*2*(rand(length(x_lattice)*length(y_lattice),1)-.5);
% radii(radii < 0) = 0;
pos_circ = [x_centers(:) y_centers(:)] + ...
    lattice_position_mult*2*(rand(nponds^2,2)-.5);

%% Place the randomly adjusted circles on the domain
for i = 1:size(pos_circ,1)
    
    cp = (pond_X - pos_circ(i,1)).^2 ...
        + (pond_Y - pos_circ(i,2)).^2 <= radii(i).^2;
    
    edge_mp(cp) = 1;
    
end

%% Compute the light field and geometry stats

CC = bwconncomp(edge_mp);
stats = table2array(regionprops('table',CC,'Area','Perimeter'));

if ~isempty(stats)
    
    %Area and Perimeter of ponds.
    A_data = dx^2*stats(:,1);
    P_data = dx*stats(:,2);
    
    % Sort it for doing gradients
    [A_data,I] = sort(A_data);
    P_data = P_data(I);
    P_data(P_data == 0) = 1;
    
    D_data = 2*log10(P_data)./log10(A_data);
    
    A_surface = sum(A_data);
    P_surface = sum(P_data);
    D_surface = 2*log10(P_surface)/log10(4*pi*A_surface);
      
    num_it = num_it + 1;
       
    pond_save{num_it} = edge_mp;
    
    A_pond_save{num_it} = A_data;
    P_pond_save{num_it} = P_data;
    D_pond_save{num_it} = D_data;
    
    A_surf_save(num_it) = A_surface;
    P_surf_save(num_it) = P_surface;
    D_surf_save(num_it) = D_surface;
    

    x_save{num_it} = pond_X;
    y_save{num_it} = pond_Y;
    z_save{num_it} = z;
        
else
    
    disp('No ponds')
    
end