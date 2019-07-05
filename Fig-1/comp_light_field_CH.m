function light = comp_light_field(X,Y,Z,meltpond)
% This code takes the melt pond meltpond mp over coordinates X, Y, and Z and
% mp_ij = 0 if it is bare ice, 1 if it is melt pond.
% Heavily borrowed from C. Katlein, Alfred-Wegener-Institut
% Bremerhaven, 2015

if min(Z) ~= 0
    Z = [0 Z];
end

leb = getLebedevSphere(110); %Generate Lebedev Sphere of low order
% leb = getLebedevSphere(5810); %Generate Lebedev Sphere of highest order

%Define light fluxes /transmittances

% Model for transmission to ice base - based on values from Light (2008)
alpha_w = .15;
alpha_i = .75;
thick = 1;
k_i = 0.7;


meltpond(meltpond==1)=exp(-thick*k_i)*(1-alpha_w); % 0.22; from before
meltpond(meltpond==0)=exp(-thick*k_i)*(1-alpha_i); % 0.04; from before

% Extinction coefficient of light
extinction_coefficient = .11; % 1/m - decay length scale is 10 m.

%% Define the grid for the new light field
nx = length(X);
ny = length(Y);

Lx = max(X(:));
Ly = max(Y(:));

dx = Lx/nx; 
dy = Ly/ny; 

mid_X = round(Lx/2);
mid_Y = round(Ly/2);

light=zeros([nx ny length(Z)]);

% Because we want to ensure that we can grab light from all the other
% points, we used to shift the matrix to put the location at the middle,
% i.e. assuming periodicity. It turns out to be much easier to make a
% larger meltpond map and just take chunks out of it.
meltpond_enhanced = repmat(meltpond,[3 3]);

%% Initialize LebedevSphere
% The Sphere is necessary to integrate light over the upper half sphere
%construct the sphere using 5810 points
%for increased speed choose a smaller value from the functions help
% leb = getLebedevSphere(5810);
%convert beam coordinated to azimuth and elevation angles (phi, theta)
[azimuth,elevation,~] = cart2sph(leb.x,leb.y,leb.z);
%select upward looking
upidx=find(elevation>0);

for d_ind = 1:length(Z)
    %%
    % At each depth level, we ask: suppose I have a sensor in the middle
    % of the domain. Which are the surface points that we will accumulate to
    % determine the light field there?
    
    % Depth of sensor
    z = Z(d_ind);
    %calculate ray length between this depth and all 50 points on the
    %Lebedev Sphere
    ray_length=z./sin(elevation);
    
    % Get the x and y coordinates of the surface points at the elevation
    % angle and the ray length. 
    [x_loc,y_loc,~]=sph2cart(azimuth,elevation,ray_length); % Units of meters
    
    % Now shift these points to the middle of the domain
    x_loc = x_loc+mid_X; % Centered location in meters
    y_loc = y_loc+mid_Y;
    
    % Next grab the melt pond cells that will send light here. 
    % Get indices of meltpond cells that give light to the midpoint one

    % Eliminate all that point outside of the domain
    x_loc(x_loc<=0 | x_loc > Lx) = NaN; 
    y_loc(y_loc<=0 | y_loc > Ly) = NaN; 
    
    % Now find the indices of these points
    x_ind = round(x_loc/dx);
    y_ind = round(y_loc/dy); 
    
    % Keep the ones that are above the current point and aren't nans
    idx=find(elevation>0 & ~isnan(x_ind) & ~isnan(y_ind));
    % A vector which will hold the surface values
    surf_val=zeros(size(x_ind));
    
    % Now for every x,y point
    for Ix = 1:nx
        
        for Iy = 1:ny
            %% 
            % Will shift the melt pond matrix to center it over (mid_X,mid_Y);
            % I used to do this explicitly and it was super expensive, so
            % now we just change the coordinates
            shift_X = nx/2 + Ix;
            shift_Y = ny/2 + Iy;
 
            % For every point that will contribute to the light field here
            for k=1:length(idx)
                % Grab the surface irradiance at that point shifted
                % appropriately
                
                % surf_val(idx(k))= ...
                %   meltpond_enhanced(x_ind(idx(k))+shift_X,y_ind(idx(k))+shift_Y);

                surf_val(idx(k))= ...
                    meltpond_enhanced(x_ind(idx(k))+shift_X,y_ind(idx(k))+shift_Y);

                
            end
            
            % Isotropic Scattering
            %  em=1; %real isotropic.
 
            % Use Katlein et al 
            gamma=0.6;
            em=((1/3)+(2/3).*cos((pi/2)-elevation)).*cos((pi/2)-elevation).*(1-gamma)...
                +gamma.*exp(((pi/2)-elevation).*0.05681);%see Katlein et al 2014, JGR
            em=em.*1.1770; %needs to be rescaled to be converted into proper escape function depending on gamma!
            
            % surf_em = surf_val.*em; % Irradiance
            surf_em=surf_val.*em./(pi); %divide by pi for irradiance to radiance conversion

            
            % Scalar Irradiance
            % sens_cell=surf_em; % will always give slightly too low estimates, as horizontal rays are not included 
            
            
            % Planar Irradiance
            % Reduce by the angle over which the light comes in
            sens_cell=surf_em.*cos((pi/2)-elevation); %
            

            
            % Extinguish along the depth z. 
            abs_sens_cell=sens_cell.*exp(-z*extinction_coefficient);
            
            % Sum up along all upward looking angles. 
            light(Ix,Iy,d_ind) =sum(abs_sens_cell(upidx).*leb.w(upidx)); 
            
        end
    end
    
end

end
