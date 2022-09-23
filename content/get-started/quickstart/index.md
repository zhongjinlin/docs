% Zhongjin Lin UBC,2022-02-22
% matlab code for the layout of speckle spectrometer



%%
clear all
close all
clc
tic
rng('shuffle');
% create a structure to hold elements (cell name)
global gs;
gs = gds_structure('rand_10');

load('r_for_random.mat');
n_r=length(result);

R=50; 
Rmin=5;
N_array=15;
gap=0.2;
si_w=0.5;

r_used=zeros(N_array,N_array);
L=zeros(N_array,N_array);

%% Create structures
for i=1:1:N_array
    for j=1:1:N_array
        a=1;
        while a==1
            nnn=floor(rand()*(n_r-1))+1;
            theta=result(nnn,4);
            r=result(nnn,2)*10^6;
            if r>Rmin
                r2=sqrt(2)/2*(R-r)/sin(theta-pi/4)-r;
                delta=sqrt(2)/2*(R-r)*(1-tan(theta-pi/4-pi/2));
                x_point=(sqrt(2)/2*(R-r)+(r+r2)*cos(theta-pi/4))*sqrt(2)/2;
                if r2>Rmin && delta>r2 %&& (theta<pi/2||x_point>r2+gap/2+si_w/2)
                    a=0;
                    r_used(i,j)=r;
                    L(i,j)=8*(r*theta+r2*(theta-pi/4));
                 %%
                    port   = [0 0 0];   % start at (0,0) point with 0 degrees
                    width  =0.5;       % arc waveguide width
                    rotate = theta;        % rotation angle of the arc
                    N      = 1200;       % number of points per arc line
                    layer  = 1;         % layer number to draw the arc on
                    x_position=(2*i-1)*(R+gap/2+width/2);
                    y_position=(2*j-1)*(R+gap/2+width/2);
                    
                    [ports,compA] = Arc_for_speckle_spectrometer(port,width,r,r2,R,rotate,N,x_position,y_position,layer);
                    
                    layer  = 3;
                    width  =4.5;
                    N      = 120;
                    [ports,compA] = Arc_for_speckle_spectrometer(port,width,r,r2,R,rotate,N,x_position,y_position,layer);
                    
                    %% Create structures
                    r_hole=sqrt(2)/2*(R-r)+(r+r2)*cos(theta-pi/4)-r2;
                    for k=1:1:4
                        center  = [x_position+r_hole*cos(pi/2*k+pi/4) y_position+r_hole*sin(pi/2*k+pi/4)]; % center of the disc
                        radius  = 0.075;    % disk radius
                        N       = 30;   % number of points per circumference
                        layer   = 4;     % layer number to draw the disc on                       
                        [compA] = Disc(center,radius,N,layer);
                    end                   
                end
            end
        end
    end
end

L_u=L/min(min(L));

%% GDS conversion

% create a library to hold the structure
glib = gds_library('example_library', 'uunit',1e-6, 'dbunit',1e-9,gs);

% finally write the library to a file
write_gds_library(glib,'GDS_rand10_R50_test_without_hole.gds');
toc
