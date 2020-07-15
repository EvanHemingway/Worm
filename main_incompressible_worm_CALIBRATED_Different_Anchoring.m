clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GEOMETRY
global ref_radius ell N
ref_radius = 0.0076/2; % 7.6 mm diameter based on Kurth and Kier's scaling of a 10 g worm (2014)
ell = 0.25; % 250 mm long based on Kurth and Kier's scaling of 10 g worm (2014)
N = 500; % Number of xi points
skips = 5; % Numvber of cross-sections to skip plotting in animation

%TIMES
t_start = 0; % start time in seconds
t_step = 0.01; % step time in seconds
t_pause = 0; % pause between cycles in seconds
t_pro = 3.15; %protrusion time for 10g worm, Quillen, kinematic scaling (1999)
t_stance = 1.75; %stance time for 10g worm, Quillen, kinematic scaling (1999)
t_end = 2*(t_pro+t_stance); % end time in seconds

%STIFFNESSES
E1 = 2E3;
A0 = pi*ref_radius^2;
global kI
kI = 2/3*E1*A0;
global kG
kG = E1*pi*ref_radius^4/12;

%INERTIAS
mass = 0.01; %Adult worms are 10 grams (Kurth and Kier, 2014)
rho_0_star = mass/pi/ref_radius^2/ell; %kg/m^3
g = 9.81; %m/s^2

%FORCING
global Dbar
peakfactor1 = 3/5;
Dbar = peakfactor1*0.0436; % 0.0436 N is the peak contractile force
peakfactor2 = 1.19/5;
p = peakfactor2*5700; % 5700 kPa is the pressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%END INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GET rho_0
rho_0 = rho_0_star*pi*ref_radius^2; % linear mass density in kg/m

%NONDIMENSIONALIZATIONS
global itilde deltatilde ctilde mtilde gamtilde etatilde
itilde = kI*ell^2/kG; deltatilde=Dbar*ell^2/kG; ctilde=p*pi*ref_radius^2*ell^2/kG; mtilde=ref_radius*rho_0*g*ell^2/kG;

%DETERMINE TIME PARAMETERS
t_tot = t_end-t_start; %total elapsed time
num_t_points = floor(t_tot/t_step); %get an integer from a continuum
t_step_discrete = t_tot/num_t_points;
t = linspace(t_start,t_end,num_t_points);
M = length(t); % total number of time indices

%BUILD XI
xi = linspace(0,ell,N);
xitilde = xi/ell;

%GRAB CRITICAL INDICES
I=floor(t_pause/t_step_discrete); % How many time indices from start until cycle initiation
J=floor((t_start+t_stance)/t_step_discrete); % How many time indices from start to end of eta
K=floor((t_pause+t_pro)/t_step_discrete); % How many time indices from start until start of new eta
L=floor((t_pause+t_stance+t_pro)/t_step_discrete); % How many time indices from start until end of a cycle

%%%%%%%%%%%%%%%%%%%BUILD BASE FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INITIALIZE GAMMA AND ETA BASE FUNCTIONS
base_gams=nan(L,1);
base_etas=nan(L,1);

k=1;
for i=I+1:J
    dist = (0.25*(xi(end-1)-xi(2))/t_stance)*(t(k)-t_start); %for a given elapsed time, find distance propagated
    [~,S] = min(abs(xi-(xi(end-1)-xi(2)-dist))); %find closest space value at that distance
    if S == 1 %ensure that we are not close to the endpoints of the space mesh
        base_gams(i) = xi(2);
    elseif S == N 
        base_gams(i) = xi(end-1);
    else
        base_gams(i) = xi(S);
    end
    k=k+1;
end
dist_end = dist;

t_final = t(k-1);

for i=J+1:L
    dist = (0.75*(xi(end-1)-xi(2))/t_pro)*(t(k)-t_final)+dist_end; %for a given elapsed time, find distance propagated
    [~,S] = min(abs(xi-(xi(end-1)-xi(2)-dist))); %find closest space value at that distance
    if S == 1 %ensure that we are not close to the endpoints of the space mesh
        base_gams(i) = xi(2);
    elseif S == N 
        base_gams(i) = xi(end-1);
    else
        base_gams(i) = xi(S);
    end
    k=k+1;
end

k=1;
for i=I+1:J
    dist = (0.5*(xi(end-1)-xi(2))/t_stance)*(t(k)-t_start); %for a given elapsed time, find distance propagated
    [~,S] = min(abs(xi-ell/2+dist)); %find closest space value at that distance
    if S == 1 %ensure that we are not close to the endpoints of the space mesh
        base_etas(i) = xi(2);
    elseif S == N 
        base_etas(i) = xi(end-1);
    else
        base_etas(i) = xi(S);
    end
    k=k+1;
end

k=1;
for i=K:L-2
    dist = (0.5*(xi(end-1)-xi(2))/t_stance)*(t(k)-t_start); %for a given elapsed time, find distance propagated
    [~,S] = min(abs(xi-(xi(end-1)-xi(2)-dist))); %find closest space value at that distance
    if S == 1 %ensure that we are not close to the endpoints of the space mesh
        base_etas(i) = xi(2);
    elseif S == N 
        base_etas(i) = xi(end-1);
    else
        base_etas(i) = xi(S);
    end
    k=k+1;
end

%%%%%%%%%%%%%%%%%%%END BUILD BASE FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%BUILD GAMMA AND ETA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INITIALIZE GAMMA AND ETA IN TIME
gams = nan(M,1);
etas = nan(M,1);

for i=1:L-4:M
    if i<M-L
    etas(i:i-1+L) = base_etas;
    end
end

j=1;
for k=i:M
    etas(k) = base_etas(j);
    j=j+1;
end

for i=1:L-4:M
    if i<M-L
    gams(i:i-1+L) = base_gams;
    end
end

j=1;
for k=i:M
    gams(k) = base_gams(j);
    j=j+1;
end

%%%%%%%%%%%%%%%%%%%END BUILD GAMMA AND ETA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NONDIMENSIONALIZE
gamtildes = gams/ell;
etatildes = etas/ell;

%PLOT THE PRESCRIBED PERISTALSIS
f3 = figure(3);
plot(t,gamtildes,'go')
hold on
plot(t,etatildes,'bo')
xlim([t_start t_end])
legend('gamtildes','etatildes')
xlabel('Time')
ylabel('Nondimensional Location')
movegui(f3,[40 0])
hold off


%INITIALIZE POSITIONS VECTORS AND DIRECTORS IN SPACE
r = zeros(N,3);
d1 = zeros(N,3);
d2 = zeros(N,3);

%INITIALIZE VOLUME MEASURE
vol = zeros(N,1);

f2 = figure(2);
xlabel('meters')
movegui(f2,[600,0])
drawnow
hold on
view(0,90)
camorbit(-90,0)

%COMPARISON WITH QUILLIN'S (1998 and 1999) SEGMENT 50 MEASUREMENT 
[~,ind_segment50] = min(abs(xi-3*ell/5));
mu1_segment50 = zeros(length(t),1);
mu3_segment50 = zeros(length(t),1);
pressure = zeros(length(t),1);
anchoring_point = zeros(length(t),1);

%GET INDICIES FOR DIRECTOR TAPER FOR WORM ANIMATION
[~,ind_tail] = min(abs(xi-49/1090*ell));
x_nondim = linspace(0,1,ind_tail);
sin_func = sin(pi/2*x_nondim)'; %sine function for tapering tail
[~,ind_neck_start] = min(abs(xi-734/1090*ell));
neck_factor = 1; %Jump for "neck"
[~,ind_neck_end] = min(abs(xi-783/1090*ell));
[~,ind_head] = min(abs(xi-1002/1090*ell));
x_nondim = linspace(0,1,N-ind_head);
cos_func = cos(pi/2*x_nondim)'; %cosine function for tapering nose

writerObj = VideoWriter('Worm_Two_Cycles_Different_Anchoring.avi'); % Name it.
writerObj.FrameRate = 1/t_step_discrete; % How many frames per second.
open(writerObj);


ind_t = 1; %Time index

[~,ind_t_pro_one_third] = min(abs(t-t_pro/3));
[~,ind_t_pro_two_third] = min(abs(t-2*t_pro/3));
[~,ind_t_pro] = min(abs(t-t_pro));
[~,ind_t_stance_half] = min(abs(t-t_pro-t_stance/2));
[~,ind_t_stance] = min(abs(t-t_pro-t_stance));

for C=1:floor(M/L)
    for i=(1+(C-1)*L):(L*C)
        sph1 = gobjects(1); %initialize singularity sphere
        sph2 = gobjects(2); %initialize second singularity sphere
        gamtilde = gamtildes(i);
        etatilde = etatildes(i);
        if isnan(gamtilde) && isnan(etatilde) %resting worm
            [mu1_temp,mu1p_temp,mu2_temp,mu2p_temp,xi_temp]=resting_incom();
            lam_temp = 0*mu1_temp;
            I = 1;
            title(strcat('Fully resting, t=',sprintf('%1.1f',t(i))))
        elseif isnan(etatilde) && ~isnan(gamtilde) %we have just one gamma: circ right, long left
            flag = 1;
            [~,J] = min(abs(xi-ell*gamtildes(i))); %get space index of the gamma
            [lam_temp,mu1_temp,mu1p_temp,mu2_temp,mu2p_temp,xi_temp]=just_gamma_incom();
            I = 1;
            title(strcat('PROTRUSION, t=',sprintf('%1.1f',t(i))))
        elseif ~isnan(gamtilde) && ~isnan(etatilde) 
            if gamtilde < etatilde %we have gamma and eta: a circular pulse in the middle
                flag = 2;
                [~,J] = min(abs(xi-ell*gamtildes(i))); %get space index of the gamma
                [~,K] = min(abs(xi-ell*etatildes(i))); %get space index of the eta
                [lam_temp,mu1_temp,mu1p_temp,mu2_temp,mu2p_temp,xi_temp]=gamma_left_eta_right_incom();
                title(strcat('STANCE, t=',sprintf('%1.1f',t(i))))
                I = N; %See Quillen (1999), she assumes anterior end to be fixed during stance phase.
            else %we have gamma and eta: a longitudinal pulse in the middle
                flag = 4;
                [~,J] = min(abs(xi-ell*gamtildes(i))); %get space index of the gamma
                [~,K] = min(abs(xi-ell*etatildes(i))); %get space index of the eta
                [lam_temp,mu1_temp,mu1p_temp,mu2_temp,mu2p_temp,xi_temp]=eta_left_gamma_right_incom();
                title(strcat('PROTRUSION, t=',sprintf('%1.1f',t(i))))
                if 1-gamtildes(i) < etatildes(i)
                    I = N;
                else
                    I = 1;
                end
            end
        else %we have just one eta: circ left, long right
            flag = 3;
            [~,K] = min(abs(xi-ell*etatildes(i))); %get space index of the eta
            [lam_temp,mu1_temp,mu1p_temp,mu2_temp,mu2p_temp,xi_temp]=just_eta_incom();
            I = N;
            title(strcat('PROTRUSION, t=',sprintf('%1.1f',t(i))))
        end
        
        mu1 = interp1(xi_temp, mu1_temp, xi);
        mu1p = interp1(xi_temp, mu1p_temp, xi);
        mu2 = interp1(xi_temp, mu2_temp, xi);
        mu2p = interp1(xi_temp, mu2p_temp, xi);
        mu3 = 1./mu1./mu2;
        
        %Track Quillin's segment 50
        mu1_segment50(ind_t) = mu1(ind_segment50);
        mu3_segment50(ind_t) = mu3(ind_segment50);
        
        
        lam = interp1(xi_temp, lam_temp, xi);
        lamreal = lam - kI*(2-mu1-mu2);
        pcoelom = -lamreal.*mu1.*mu2/pi/ref_radius^2;
        pressure(ind_t) = pcoelom(ind_segment50);
        
        anchoring_point(ind_t) = xi(I)/ell;
        
        vol = mu1.*mu2.*mu3;
        
        theta = 180/pi*asin(ref_radius*mu1p./mu3);
        
        %For an azimuth and elevation of 30 degrees, Matlab's basis vectors, M1 M2
        %and M3 are such that M1 is directed rightward, M2 into the page, and M3
        %up. In our theory, this view has E3 rightward, E2 out of the page, and E1
        %up. Therefore, M1 = E3, M2 = -E2, M3 = E1. For plotting purposes, the
        %vectors are constructed on the MATLAB basis.
        d1_temp = ref_radius*mu1'; %d1 in E1 in theory, M3 for MATLAB
        d2_temp = -ref_radius*mu2'; %d2 in E2 in theory, -M2 for MATLAB
        
        %Taper the directors for worm animation
        d1(1:ind_tail,3) = sin_func.*d1_temp(1:ind_tail);
        d1(ind_tail+1:ind_neck_start,3) = d1_temp(ind_tail+1:ind_neck_start);
        d1(ind_neck_start+1:ind_neck_end,3) = neck_factor*d1_temp(ind_neck_start+1:ind_neck_end);
        d1(ind_neck_end+1:ind_head,3) = d1_temp(ind_neck_end+1:ind_head);
        d1(ind_head+1:end,3) = cos_func.*d1_temp(ind_head+1:end);
        
        d2(1:ind_tail,2) = sin_func.*d2_temp(1:ind_tail);
        d2(ind_tail+1:ind_neck_start,2) = d2_temp(ind_tail+1:ind_neck_start);
        d2(ind_neck_start+1:ind_neck_end,2) = neck_factor*d2_temp(ind_neck_start+1:ind_neck_end);
        d2(ind_neck_end+1:ind_head,2) = d2_temp(ind_neck_end+1:ind_head);
        d2(ind_head+1:end,2) = cos_func.*d2_temp(ind_head+1:end);
        
        rhold=r(I,1);
        
        if i==1 %first iteration, integrate from the nose rearward
            r(1,1)=0;
            for j=2:length(xi)
                r(j,1)=r(j-1,1)+sqrt(mu3(j)^2-ref_radius^2*mu1p(j)^2)*(xi(j)-xi(j-1)); %build r along M1
            end
        else
            r(I,1)=rhold;
            if I==1
                for j=2:N
                    r(j,1)=r(j-1,1)+sqrt(mu3(j)^2-ref_radius^2*mu1p(j)^2)*(xi(j)-xi(j-1)); %build r along M1
                end
            elseif I==N
                for j=I-1:-1:1
                    r(j,1)=r(j+1,1)-sqrt(mu3(j)^2-ref_radius^2*mu1p(j)^2)*(xi(j+1)-xi(j)); %build r along M1
                end
            else
                for j=I+1:N
                    r(j,1)=r(j-1,1)+sqrt(mu3(j)^2-ref_radius^2*mu1p(j)^2)*(xi(j)-xi(j-1)); %build r along M1
                end
                for j=I-1:-1:1
                    r(j,1)=r(j+1,1)-sqrt(mu3(j)^2-ref_radius^2*mu1p(j)^2)*(xi(j+1)-xi(j)); %build r along M1
                end
            end
        end
        r(:,3)=ref_radius*(mu1'-1);
        
        cross_sec_points = 100;
        
        [V,rod] = circular_rod_plot(r,d1,d2,100);
        hold on
        
        ind_base = 1:cross_sec_points;
        
        for ii=1:floor(N/skips)
            inds_temp = (ii*skips-1)*cross_sec_points+ind_base;
            indflip1 = inds_temp(floor(cross_sec_points/2)+1:end);
            indflip2 = inds_temp(1:floor(cross_sec_points/2));
            inds((ii-1)*cross_sec_points+1:ii*cross_sec_points) = [indflip1, indflip2];
        end
        
        circles = plot3(V(inds,1),V(inds,2),V(inds,3),'k');
        
        %q1 = quiver3(0,0,0,1,0,0,'b'); %M1 = E3
        %q2 = quiver3(0,0,0,0,-1,0,'g'); %M2 = -E2
        %q3 = quiver3(0,0,0,0,0,1,'r'); %M3 = E1

        
        XXXX = [0 2*ell 2*ell 0];
        YYYY = [-4*ref_radius -4*ref_radius 4*ref_radius 4*ref_radius];
        ZZZZ = [-ref_radius -ref_radius -ref_radius -ref_radius];
        plane = patch(XXXX,YYYY,ZZZZ,'white');
        
        
%         if flag == 1
%             [XXX,YYY,ZZZ]=sphere;
%             sph1 = surf(XXX*0.01+r(J,1), YYY*0.01, ZZZ*0.01-ref_radius,'FaceColor', [0 1 0]);
%         elseif flag == 2
%             [XXX,YYY,ZZZ]=sphere;
%             sph1 = surf(XXX*0.01+r(J,1), YYY*0.01, ZZZ*0.01-ref_radius,'FaceColor', [0 1 0]);
%             [XXXX,YYYY,ZZZZ]=sphere;
%             sph2 = surf(XXXX*0.01+r(K,1), YYYY*0.01, ZZZZ*0.01-ref_radius,'FaceColor', [0 0 1]);
%         elseif flag == 3
%             [XXX,YYY,ZZZ]=sphere;
%             sph1 = surf(XXX*0.01+r(K,1), YYY*0.01, ZZZ*0.01-ref_radius,'FaceColor', [0 0 1]);
%         elseif flag == 4
%             [XXX,YYY,ZZZ]=sphere;
%             sph1 = surf(XXX*0.01+r(J,1), YYY*0.01, ZZZ*0.01-ref_radius,'FaceColor', [0 1 0]);
%             [XXXX,YYYY,ZZZZ]=sphere;
%             sph2 = surf(XXXX*0.01+r(K,1), YYYY*0.01, ZZZZ*0.01-ref_radius,'FaceColor', [0 0 1]);
%         end
        
        %[XX,YY,ZZ]=sphere;
        %sph = surf(XX*ref_radius/2+rhold, YY*ref_radius/2, ZZ*ref_radius/2-ref_radius,'FaceColor', [1 0 0]); %anchoring sphere
        axis([0 2*ell -ell/2 ell/2 -ell/2 ell/2])
        %    axis equal
        axis manual
        %    axis off
%         [XX,YY,ZZ]=sphere;
%         sph_seg50 = surf(XX*ref_radius/2+r(ind_segment50,1), YY*ref_radius/2+r(ind_segment50,2), ZZ*ref_radius/2+r(ind_segment50,3)+2*ref_radius,'FaceColor', [1 1 1]);
%         
%         if ind_t == 1 || ind_t == ind_t_pro_one_third || ind_t == ind_t_pro_two_third || ind_t == ind_t_pro || ind_t == ind_t_stance_half || ind_t == ind_t_stance
%             title('')
%             axis off
%             grid off
%             delete(plane)
%             pause()
%         end
        view(59,17)
        zoom(2)
axis off
grid off
pause(0.01)
frame = getframe(gcf);
        writeVideo(writerObj, frame);
        delete(rod)
        %delete(sph)
        delete(sph1)
        delete(sph2)
        delete(circles)
        delete(plane)
        %delete(q1)
        %delete(q2)
        %delete(q3)
        ind_t = ind_t+1;
    end
end


%%%%%%%%%%%%%%%%%COMPLETE THE FINAL CYCLE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if floor(M/L) == 0
    C=0;
end

for i=(L*C+1):length(t)
    sph1 = gobjects(1); %initialize singularity sphere
    sph2 = gobjects(2); %initialize second singularity sphere
    gamtilde = gamtildes(i);
    etatilde = etatildes(i);
    if isnan(gamtilde) && isnan(etatilde) %resting worm
        [mu1_temp,mu1p_temp,mu2_temp,mu2p_temp,xi_temp]=resting_incom();
        lam_temp = 0*mu1_temp;
        I = 1;
        title(strcat('Fully resting, t=',sprintf('%1.1f',t(i))))
    elseif isnan(etatilde) && ~isnan(gamtilde) %we have just one gamma: circ right, long left
        flag = 1;
        [~,J] = min(abs(xi-ell*gamtildes(i))); %get space index of the gamma
        [lam_temp,mu1_temp,mu1p_temp,mu2_temp,mu2p_temp,xi_temp]=just_gamma_incom();
        I = 1;
        title(strcat('PROTRUSION, t=',sprintf('%1.1f',t(i))))
    elseif ~isnan(gamtilde) && ~isnan(etatilde)
        if gamtilde < etatilde %we have gamma and eta: a circular pulse in the middle
            flag = 2;
            [~,J] = min(abs(xi-ell*gamtildes(i))); %get space index of the gamma
            [~,K] = min(abs(xi-ell*etatildes(i))); %get space index of the eta
            [lam_temp,mu1_temp,mu1p_temp,mu2_temp,mu2p_temp,xi_temp]=gamma_left_eta_right_incom();
            title(strcat('STANCE, t=',sprintf('%1.1f',t(i))))
            I = N; %See Quillen (1999), she assumes anterior end to be fixed during stance phase.
        else %we have gamma and eta: a longitudinal pulse in the middle
            flag = 4;
            [~,J] = min(abs(xi-ell*gamtildes(i))); %get space index of the gamma
            [~,K] = min(abs(xi-ell*etatildes(i))); %get space index of the eta
            [lam_temp,mu1_temp,mu1p_temp,mu2_temp,mu2p_temp,xi_temp]=eta_left_gamma_right_incom();
            title(strcat('PROTRUSION, t=',sprintf('%1.1f',t(i))))
            if 1-gamtildes(i) < etatildes(i)
                I=N;
            else
                I = 1;
            end
        end
    else %we have just one eta: circ left, long right
        flag = 3;
        [~,K] = min(abs(xi-ell*etatildes(i))); %get space index of the eta
        [lam_temp,mu1_temp,mu1p_temp,mu2_temp,mu2p_temp,xi_temp]=just_eta_incom();
        I = N;
        title(strcat('PROTRUSION, t=',sprintf('%1.1f',t(i))))
    end
    
    mu1 = interp1(xi_temp, mu1_temp, xi);
    mu1p = interp1(xi_temp, mu1p_temp, xi);
    mu2 = interp1(xi_temp, mu2_temp, xi);
    mu2p = interp1(xi_temp, mu2p_temp, xi);
    mu3 = 1./mu1./mu2;
    
    %Track Quillin's segment 50
    mu1_segment50(ind_t) = mu1(ind_segment50);
    mu3_segment50(ind_t) = mu3(ind_segment50);
    
    
    lam = interp1(xi_temp, lam_temp, xi);
    lamreal = lam - kI*(2-mu1-mu2);
    pcoelom = -lamreal.*mu1.*mu2/pi/ref_radius^2;
    pressure(ind_t) = pcoelom(ind_segment50);
    
    vol = mu1.*mu2.*mu3;
    
    theta = 180/pi*asin(ref_radius*mu1p./mu3);
    
    %For an azimuth and elevation of 30 degrees, Matlab's basis vectors, M1 M2
    %and M3 are such that M1 is directed rightward, M2 into the page, and M3
    %up. In our theory, this view has E3 rightward, E2 out of the page, and E1
    %up. Therefore, M1 = E3, M2 = -E2, M3 = E1. For plotting purposes, the
    %vectors are constructed on the MATLAB basis.
    d1_temp = ref_radius*mu1'; %d1 in E1 in theory, M3 for MATLAB
    d2_temp = -ref_radius*mu2'; %d2 in E2 in theory, -M2 for MATLAB
    
    %Taper the directors for worm animation
    d1(1:ind_tail,3) = sin_func.*d1_temp(1:ind_tail);
    d1(ind_tail+1:ind_neck_start,3) = d1_temp(ind_tail+1:ind_neck_start);
    d1(ind_neck_start+1:ind_neck_end,3) = neck_factor*d1_temp(ind_neck_start+1:ind_neck_end);
    d1(ind_neck_end+1:ind_head,3) = d1_temp(ind_neck_end+1:ind_head);
    d1(ind_head+1:end,3) = cos_func.*d1_temp(ind_head+1:end);
    
    d2(1:ind_tail,2) = sin_func.*d2_temp(1:ind_tail);
    d2(ind_tail+1:ind_neck_start,2) = d2_temp(ind_tail+1:ind_neck_start);
    d2(ind_neck_start+1:ind_neck_end,2) = neck_factor*d2_temp(ind_neck_start+1:ind_neck_end);
    d2(ind_neck_end+1:ind_head,2) = d2_temp(ind_neck_end+1:ind_head);
    d2(ind_head+1:end,2) = cos_func.*d2_temp(ind_head+1:end);
    
    rhold=r(I,1);
    
    if i==1 %first iteration, integrate from the nose rearward
        r(1,1)=0;
        for j=2:length(xi)
            r(j,1)=r(j-1,1)+sqrt(mu3(j)^2-ref_radius^2*mu1p(j)^2)*(xi(j)-xi(j-1)); %build r along M1
        end
    else
        r(I,1)=rhold;
        if I==1
            for j=2:N
                r(j,1)=r(j-1,1)+sqrt(mu3(j)^2-ref_radius^2*mu1p(j)^2)*(xi(j)-xi(j-1)); %build r along M1
            end
        elseif I==N
            for j=I-1:-1:1
                r(j,1)=r(j+1,1)-sqrt(mu3(j)^2-ref_radius^2*mu1p(j)^2)*(xi(j+1)-xi(j)); %build r along M1
            end
        else
            for j=I+1:N
                r(j,1)=r(j-1,1)+sqrt(mu3(j)^2-ref_radius^2*mu1p(j)^2)*(xi(j)-xi(j-1)); %build r along M1
            end
            for j=I-1:-1:1
                r(j,1)=r(j+1,1)-sqrt(mu3(j)^2-ref_radius^2*mu1p(j)^2)*(xi(j+1)-xi(j)); %build r along M1
            end
        end
    end
    r(:,3)=ref_radius*(mu1'-1);
    
    cross_sec_points = 100;
    
    [V,rod] = circular_rod_plot(r,d1,d2,100);
    hold on
    
    ind_base = 1:cross_sec_points;
    
    for ii=1:floor(N/skips)
        inds_temp = (ii*skips-1)*cross_sec_points+ind_base;
        indflip1 = inds_temp(floor(cross_sec_points/2)+1:end);
        indflip2 = inds_temp(1:floor(cross_sec_points/2));
        inds((ii-1)*cross_sec_points+1:ii*cross_sec_points) = [indflip1, indflip2];
    end
    
    circles = plot3(V(inds,1),V(inds,2),V(inds,3),'k');
    
    %q1 = quiver3(0,0,0,1,0,0,'b'); %M1 = E3
    %q2 = quiver3(0,0,0,0,-1,0,'g'); %M2 = -E2
    %q3 = quiver3(0,0,0,0,0,1,'r'); %M3 = E1
    
    
    XXXX = [0 2*ell 2*ell 0];
    YYYY = [-4*ref_radius -4*ref_radius 4*ref_radius 4*ref_radius];
    ZZZZ = [-ref_radius -ref_radius -ref_radius -ref_radius];
    plane = patch(XXXX,YYYY,ZZZZ,'white');
    
    
    %         if flag == 1
    %             [XXX,YYY,ZZZ]=sphere;
    %             sph1 = surf(XXX*0.01+r(J,1), YYY*0.01, ZZZ*0.01-ref_radius,'FaceColor', [0 1 0]);
    %         elseif flag == 2
    %             [XXX,YYY,ZZZ]=sphere;
    %             sph1 = surf(XXX*0.01+r(J,1), YYY*0.01, ZZZ*0.01-ref_radius,'FaceColor', [0 1 0]);
    %             [XXXX,YYYY,ZZZZ]=sphere;
    %             sph2 = surf(XXXX*0.01+r(K,1), YYYY*0.01, ZZZZ*0.01-ref_radius,'FaceColor', [0 0 1]);
    %         elseif flag == 3
    %             [XXX,YYY,ZZZ]=sphere;
    %             sph1 = surf(XXX*0.01+r(K,1), YYY*0.01, ZZZ*0.01-ref_radius,'FaceColor', [0 0 1]);
    %         elseif flag == 4
    %             [XXX,YYY,ZZZ]=sphere;
    %             sph1 = surf(XXX*0.01+r(J,1), YYY*0.01, ZZZ*0.01-ref_radius,'FaceColor', [0 1 0]);
    %             [XXXX,YYYY,ZZZZ]=sphere;
    %             sph2 = surf(XXXX*0.01+r(K,1), YYYY*0.01, ZZZZ*0.01-ref_radius,'FaceColor', [0 0 1]);
    %         end
    
    %[XX,YY,ZZ]=sphere;
    %sph = surf(XX*ref_radius/2+rhold, YY*ref_radius/2, ZZ*ref_radius/2-ref_radius,'FaceColor', [1 0 0]); %anchoring sphere
    axis([0 2*ell -ell/2 ell/2 -ell/2 ell/2])
    %    axis equal
    axis manual
            view(59,17)
            zoom(2)
    %    axis off
    axis off
    grid off
    pause(0.01)
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    if i< length(t)
        delete(rod)
        %delete(sph)
        delete(sph1)
        delete(sph2)
        delete(circles)
        delete(plane)
        %delete(q1)
        %delete(q2)
        %delete(q3)
    end
    ind_t = ind_t+1;
end

close(writerObj)

%%%%%%%%%%%%%%%%%%%%%%%%%PLOT DIRECTORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quiver3(r(:,1),r(:,2),r(:,3),d1(:,1),d1(:,2),d1(:,3),'r')
% quiver3(r(:,1),r(:,2),r(:,3),d2(:,1),d2(:,2),d2(:,3),'g')

%%%%%%%%%%%%%%%%%%%%%%PLOT INERTIAL BASIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quiver3(0,0,0,1,0,0,'b') %M1 = E3
% quiver3(0,0,0,0,-1,0,'g') %M2 = -E2
% quiver3(0,0,0,0,0,1,'r') %M3 = E1

%%%%%%%%%%%%%%%%%%%%%%%PLOT GROUND PLANE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = [-ell/2 2*ell 2*ell -ell/2];
% y = [-4*ref_radius -4*ref_radius 4*ref_radius 4*ref_radius];
% z = [-ref_radius -ref_radius -ref_radius -ref_radius];
% patch(x,y,z,'white')

%%%%%%%%%%%%%%%%%%%%%%%PLOT STRETCHES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(xi,mu1,'r*')
% hold on
% plot(xi,mu2,'go')
% plot(xi,mu3,'b-')
% legend('mu1','mu2','mu3')

%%%%%%%%%%%%%%%%%%%PLOT STRETCHES GRADIENTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(xi,mu1p,'r*')
% hold on
% plot(xi,mu2p,'go')
% legend('mu1p','mu2p')