
%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 SIMULATION SETUP %%%%%%%%%%%%%%%%%%%%%%%%
close all
rng(156, 'twister')
clear

%% Step 0: Diamond as a Medium
%Since we use a lot of E&M calculations, we need
DiaPerm = 5.66*8.854e-12;       %Permittivity in Diamond
EC = 1.602e-19/(4*pi*DiaPerm);  %Charge / 4*pi*epsilom

%In diamond there are two types of chargetraps

%P1 traps which can be neutral or negatively charged
P1occ = 13;            %Suspected Precent that are charged negative
P1dens = 5e21;         %Suspected Density of traps

%Vacancies which can be neutral or positively charged
VCocc = 13;            %Suspected Precent that are charged positive
VCdens = 5e21;         %Suspected Density of traps

%Finally, there are NV centers within Diamond which can be neutral or
%negatively charged, and are treated as dipoles.
NVocc = 82;            %Suspected Precent that are charged negative
NVdens = 0.2e21;         %Suspected Density of NV's


%% Step 1: Creating a Volume
%This step is pretty straight forward, for our simulation we need a volume
%of diamond to be used in the simulation. Thus we consider a box shape

BoxDim = 1e-6*[2;2;2];  %Box Size

%Here's a visualization of the box%
figure(1)  
BoxX = BoxDim(1)*[0,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0];   %Mess of points
BoxY = BoxDim(2)*[0,0,1,1,0,0,1,1,1,1,1,1,0,0,0,0];   %that create
BoxZ = BoxDim(3)*[0,0,0,0,0,1,1,0,1,1,0,1,1,0,1,1];   %the box shape
hold on
plot3(BoxX,BoxY,BoxZ,'b-','LineWidth',2);
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
axis equal


%% Step 2: Add Charge Traps Within Box
%We then estimate how many of each trap exists within the volume
P1num = round(P1dens*BoxDim(1)*BoxDim(2)*BoxDim(3));
VCnum = round(VCdens*BoxDim(1)*BoxDim(2)*BoxDim(3));

%Then we place the charge traps within the volume
P1pos(:,1:P1num) = [BoxDim(1)*rand(1,P1num);BoxDim(2)*rand(1,P1num);BoxDim(3)*rand(1,P1num)];
VCpos(:,1:VCnum) = [BoxDim(1)*rand(1,VCnum);BoxDim(2)*rand(1,VCnum);BoxDim(3)*rand(1,VCnum)];

%Here's what that looks like
scatter3(P1pos(1,:),P1pos(2,:),P1pos(3,:),2,'red','filled');
scatter3(VCpos(1,:),VCpos(2,:),VCpos(3,:),2,'blue','filled');


%% Step 3: Add NV's into Box
%Next we have to add our NV's into the box, same steps as before
NVnum = round(NVdens*BoxDim(1)*BoxDim(2)*BoxDim(3));
NVpos(:,1:NVnum) = [BoxDim(1)*rand(1,NVnum);BoxDim(2)*rand(1,NVnum,1);BoxDim(3)*rand(1,NVnum)];

%However, NV's are dipoles and thus have a random orientation
%Luckily because diamond is crystalline, there are only 4 possible
%orientations
n_poss = 1/sqrt(3)*[1,1,-1,-1;1,-1,1,-1;1,-1,-1,1]; %1/sqrt(3) is the normalization factor)
NVori = n_poss(:,randi([1,4],1,NVnum));

%Here's what the NV's look like with their orientations
scatter3(NVpos(1,:),NVpos(2,:),NVpos(3,:),5,'green','filled');
quiver3(NVpos(1,:),NVpos(2,:),NVpos(3,:),NVori(1,:),NVori(2,:),NVori(3,:),0.5,'green');

%we also have one more inherant constraint in NV's. NV's have inherent
%strains caused by defects within diamonds crystalline structure. As of
%right now the strain orientation is seemingly random. The magnitude is
%close to 10 GHZ in our spectrum, we convert this into SI using muE.
muE = 6.3*1e-6 ;    %GHz/(V/m) this specific value is from Tamarat et al

STmag = 10 / muE;  %10 GHZ to V/m
STvect = randi([-100,100],3,NVnum); %Creates a random vector per NV
STvect = STmag * STvect ./ vecnorm(STvect); %Normalizes the vector and applies the magnitude to it


%% Step 4: Reducing our scope
%Obviously this entire situation is a tad intense, So let's reduce it down
%Firstly lets focus on one NV, call this our main NV. Let's look at an NV
%close to the center that has 3 adjacent neighbors
%Firstly find the NV closest to the center
[~,A] = min(vecnorm(NVpos - BoxDim/2));

[~,B] = mink(vecnorm(NVpos - NVpos(:,A)),NVnum);
NVpos = NVpos(:,B);
NVori = NVori(:,B);

%Label our center NV in gold, the next 3 closest make larger in green
scatter3(NVpos(1,1),NVpos(2,1),NVpos(3,1),40,[1,0.7,0],'filled','MarkerEdgeColor','black');
scatter3(NVpos(1,2:4),NVpos(2,2:4),NVpos(3,2:4),40,'g','filled','MarkerEdgeColor','black');

%now overwrite our NV_pos and center everything about our main NV
P1pos = P1pos - NVpos(:,1);
[~,B] = mink(vecnorm(P1pos),P1num);
P1pos = P1pos(:,B);

VCpos = VCpos - NVpos(:,1);
[~,B] = mink(vecnorm(VCpos),P1num);
VCpos = VCpos(:,B);

NVpos = NVpos - NVpos(:,1);


%we will also remove any charges that are too far that they cause minimal
%interaction
MinDist = 0.5e-6;
P1pos = P1pos(:,(vecnorm(P1pos) < MinDist));
P1num = size(P1pos,2);
VCpos = VCpos(:,(vecnorm(VCpos) < MinDist));
VCnum = size(VCpos,2);
NVpos = NVpos(:,(vecnorm(NVpos) < MinDist));
NVori = NVori(:,(vecnorm(NVpos) < MinDist));
NVnum = size(NVpos,2);

%Now we can look at our new situation
figure(2)
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
axis(MinDist*[-1,1,-1,1,-1,1])

hold on
%Place Main NV
scatter3(NVpos(1,1),NVpos(2,1),NVpos(3,1),40,[1,0.7,0],'filled','MarkerEdgeColor','black');
quiver3(NVpos(1,1),NVpos(2,1),NVpos(3,1),NVori(1,1),NVori(2,1),NVori(3,1),1e-7,'color',[1,0.7,0]);

%Place closest 3 NV's
scatter3(NVpos(1,2:4),NVpos(2,2:4),NVpos(3,2:4),40,'g','filled','MarkerEdgeColor','black');
quiver3(NVpos(1,2:4),NVpos(2,2:4),NVpos(3,2:4),NVori(1,2:4),NVori(2,2:4),NVori(3,2:4),0.5,'green');

%Place all other NV's & charges
scatter3(NVpos(1,5:end),NVpos(2,5:end),NVpos(3,5:end),8,'g','filled');
scatter3(P1pos(1,:),P1pos(2,:),P1pos(3,:),8,'blue','filled');
scatter3(VCpos(1,:),VCpos(2,:),VCpos(3,:),8,'red','filled');



%% Step 4: Randomize Charge Per loop
%Now we will look at what happens every loop
LPnum = 100;
%Our First Loop is initialized randomly based on our occupancy values
P1stat = (100*rand(P1num,1) < P1occ);
VCstat = (100*rand(VCnum,1) < VCocc);
NVstat = (100*rand(NVnum,1) < NVocc);
%Just make the first 4 NV's active to start with, makes life easier
NVstat(1:4) = 1;

%This creates a vector where the columns represent the individual particles,
%and the rows represent what happens on the first loop. A value of 0 means the
%particle is neutral, and a value of 1 means the particle is charged,
%negatively or positively depending on the type.

%this is mainly useful for the mathematics, but it is interesting to see
%how it changes our situation for one loop
figure(3)
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
axis(MinDist*[-1,1,-1,1,-1,1])

hold on
%Place all NV's and Charges, colored based on if they're on or off
scatter3(NVpos(1,NVstat),NVpos(2,NVstat),NVpos(3,NVstat),8,'green','filled');
scatter3(NVpos(1,~NVstat),NVpos(2,~NVstat),NVpos(3,~NVstat),1,'black');
scatter3(P1pos(1,P1stat),P1pos(2,P1stat),P1pos(3,P1stat),8,'blue','filled');
scatter3(P1pos(1,~P1stat),P1pos(2,~P1stat),P1pos(3,~P1stat),1,'black','filled');
scatter3(VCpos(1,VCstat),VCpos(2,VCstat),VCpos(3,VCstat),8,'red','filled');
scatter3(VCpos(1,~VCstat),VCpos(2,~VCstat),VCpos(3,~VCstat),1,'black','filled');
%Place Main NV
scatter3(NVpos(1,1),NVpos(2,1),NVpos(3,1),40,[1,0.7,0],'filled','MarkerEdgeColor','black');
quiver3(NVpos(1,1),NVpos(2,1),NVpos(3,1),NVori(1,1),NVori(2,1),NVori(3,1),1e-7,'color',[1,0.7,0]);
%Place closest 3 NV's
scatter3(NVpos(1,2:4),NVpos(2,2:4),NVpos(3,2:4),40,'g','filled','MarkerEdgeColor','black');
quiver3(NVpos(1,2:4),NVpos(2,2:4),NVpos(3,2:4),NVori(1,2:4),NVori(2,2:4),NVori(3,2:4),0.5,'green');


%However We still have LPnum-1 situations remaining. For this we're
%going to consider random amounts of change on NVstat, P1stat and VCstat
%First let's initialize the remaining loops
NVstat = [NVstat,zeros(NVnum,LPnum-1)];
VCstat = [VCstat,zeros(VCnum,LPnum-1)];
P1stat = [P1stat,zeros(P1num,LPnum-1)];
Changes = zeros(2,LPnum);
%Changes is a matrix that describes every change
%the first row is dedicated to knowing what stat we're changing (NV,VC,P1)
%the second row is dedicated to knowing exactly what particle got changed

for ll = 2:LPnum
Changes(1,ll) = randi([1,3],1);
%Each loop can change a value from each of the 3 options.
%However, there is only 1 change per loop for simplicity
%So firsty, copy over the previous loops statistics.
NVstat(:,ll) = NVstat(:,ll-1);
P1stat(:,ll) = P1stat(:,ll-1);
VCstat(:,ll) = VCstat(:,ll-1);

switch Changes(1,ll)
    case 1  %An NV changes charge
        Changes(2,ll) = randi([5,5+floor(NVnum/50)],1);
        %Only look at the first eighth of our total NV's. These are the
        %closest to our actual NV and will be more apparent when a change
        %happens, all other NV's can be considered static noise.
        %Decide which NV's will be changed
        NVstat(Changes(2,ll),ll) = ~NVstat(Changes(2,ll),ll);
        %Flip those NV values

    case 2  %A P1 trap changes charge
        %Same step for the charges
        Changes(2,ll) = randi([1,5],1);
        %Again only look at the first 8th
        P1stat(Changes(2,ll),ll) = ~P1stat(Changes(2,ll),ll);

    case 3  %A Vacancy changes charge
        Changes(2,ll) = randi([1,5],1);
        VCstat(Changes(2,ll),ll) = ~VCstat(Changes(2,ll),ll);
end %Switch Case
end %ll Loop

%And Add a little thing to make sure all the stats are seen as logical
NVstat = (NVstat == 1);
P1stat = (P1stat == 1);
VCstat = (VCstat == 1);



%% Step 5: Electric Fields
%Now we look at how each loop effects the percieved electric fields at the
%point of each of our NV's

Efield = zeros(3,LPnum,NVnum);
for vv = 1:NVnum
%on a NV by NV Basis

%Firstly we add strain from each individual NV, strain doesnt
%change per loop so we can add it to all loops simultaneously
Efield(:,:,vv) = Efield(:,:,vv) + STvect(:,vv);

for ll = 1:LPnum

%Imporantly, An NV only reads out the electricfield if it's active,
%otherwise it is given a value of NaN
if NVstat(vv,ll) == 0
    Efield(:,ll,vv) = [NaN;NaN;NaN];
else

%Add contributions from our P1 Traps & Vacancies
Efield(:,ll,vv) = Efield(:,ll,vv) + EC*sum((P1pos(:,P1stat(:,ll)) - NVpos(:,vv)) ./ vecnorm(P1pos(:,P1stat(:,ll)) - NVpos(:,vv)).^3,2);
Efield(:,ll,vv) = Efield(:,ll,vv) - EC*sum((VCpos(:,VCstat(:,ll)) - NVpos(:,vv)) ./ vecnorm(VCpos(:,VCstat(:,ll)) - NVpos(:,vv)).^3,2);

%Add contributions from other NV's
NV_pospos = NVpos(:,[1:vv-1,vv+1:NVnum]) - NVpos(:,vv); %Relative NV Position Difference Per NV
Efield(:,ll,vv) = Efield(:,ll,vv) - EC*sum((NV_pospos(:,NVstat([1:vv-1,vv+1:NVnum],ll))) ./ vecnorm(NV_pospos(:,NVstat([1:vv-1,vv+1:NVnum],ll))).^3,2);

end %NV Active Test
end %ll loop
end %vv loop



%% Step 6: Longitudinal and Transversal Fields
%In our actual readouts we know the position of the spectral lines are the
%sum and difference of the projection and rejection fields respect to the
%NV orientation. We call the vector projection the Longitudinal field, and
%the rejection field the transversal field.
Elong = zeros(LPnum,NVnum);
Etran = zeros(LPnum,NVnum);
for vv = 1:NVnum
for ll = 1:LPnum
Elong(ll,vv) = dot(Efield(:,ll,vv) , NVori(:,vv));
Etran(ll,vv) = norm(Efield(:,ll,vv) - NVori(:,vv)*Elong(ll,vv));
end
end

%Here we'll display the values of the longitudinal field against the
%transversal field per NV
figure(4)
tiledlayout(2,2)
for vv = 1:4
nexttile(vv,[1,1])
scatter(Elong(:,vv)*muE,Etran(:,vv)*muE,'+')
axis square
xlabel('Longitudinal (GHz)')
ylabel('Transversal (GHz)')
title("NV " + vv)
end


%% Step 7: 1 NV Localization
%Now we have to locate where our change was. Firstly we'll consider
%locating our change with only our primary NV.

%A big issue with this localization, is that the transversal field
%condenses 2 dimensions into 1. In doing so the math gets a little
%complicated.

%Firstly let's look at what we have.
% Elong(1,1) == The longitudinal field of NV1 during the first loop
% Etran(1,1) == The transversal field of NV1 during the first loop
% Elong(2,1) == The longitudinal field of NV1 during the second loop
% Etran(2,1) == The transversal field of NV1 during the second loop

%There is some charge that changes that turns our first loop into our
%second loop.
%Let's call the field without this charge S for Strain, the field the
%charge gives C for charge, and the total field S + C = E

%Of course all these values are 3d-vectors. Lets refer to these 3
%dimensions as such
% n = Orientation of NV1
% s = Direction of S made perpendicular to n
% a = n × a     (the × is the cross product)

%In doing this we have a complete orthonormal basis. Thus our fields take
%the form
% E = n*n*E + s*s*E + a*a*E
% C = n*n*C + s*s*C + a*a*C
% S = n*n*S + s*s*S        (the way s is defined leads to S having no 'a' component)

%Our longitudinal fields are easy enough
% n(n*E) = n(n*C + n*S)
% n*E = n*C + n*S

%Our transversal fields not so much, as the dimensions get broken down
% magET = sqrt( (s*E)² + (a*E)² )
%Luckily instead of magnitude we can say that Eperp is 
%along the transverse plane in an arbitrary direction
% vec_E⊥ = E⊥*(s*cosθ + a*sinθ)
%In total our electric field looks like
% E = n*E∥ + E⊥*(s*cosθ + a*sinθ)

%If we solve for C, we find that
% C = E - S
% C = n*E∥ + E⊥(s*cosθ + a*sinθ) - n*S∥ - s*S⊥
% C = n(E∥ - S∥) + s(E⊥cosθ - S⊥) + a(E⊥sinθ)

%Now we have to convert this field into a position.
%To do this we use the function
% vec_E = Ec * vec_R / R^3
%We can remove one R by just looking at the unit vector
% vec_E = Ec * unit_R / R²
%If we solve this for direction and magnitude seperatly we see
% unit_R = ±unit_E = ±vec_E / mag_E
% mag_R = sqrt( Ec / mag_E )
%Where the ± is the charge of the particle and whether it was added or removed
%Thus in total
% vec_R = ±vec_E / mag_E * sqrt( Ec / mag_E )
% vec_R = ±vec_E * sqrt( Ec / mag_E^3/2 )
% vec_R = ±vec_E * sqrt(Ec) * mag_E^(-3/4)

%If we now plug in C for our E value, we find
% vec_R = ±( n(E∥ - S∥) + s(E⊥cosθ - S⊥) + a(E⊥sinθ) ) * sqrt(Ec)
%       * (E∥² + E⊥² + S∥² + S⊥² - 2E∥S∥ - 2E⊥S⊥cosθ )^(-3/4)

%Now for some annoyance, Although we know some unit vector s exists
%we don't know where it is exactly. Thus we have to rotate along the 
%s,a plane until we find the s value that matches
%Thus to initialize the setup we have an arbitrary s' value, which for
%simplicity I discribe as n cross [1;0;0], then normalized. we can then get
% a' = n' × s'
%thus our true s value is actually 
% s = s'cosφ + a'sinφ
%and our true a value is
% a = s'sinφ - a'cosφ
%and our E value is now
%To make this easier on ourselves, we'll rotate it all, after we find the
%values, speaking of, let's start writing some code

ThetaRes = 360; %How many test points for theta

%En = 4.6488e+03;    %E∥
%Ep = 1.7156e+04;    %E⊥
%Sn = 2.5892e+03;    %S∥
%Sp = 1.5530e+04;    %S⊥

En = Elong(2,1);    %E∥
Ep = Etran(2,1);    %E⊥
Sn = Elong(1,1);    %S∥
Sp = Etran(1,1);    %S⊥

Theta = linspace(0,2*pi,ThetaRes);
R = zeros(3,ThetaRes);
for tt = 1:ThetaRes
T = Theta(tt);      %this makes this look nicer, but is useless

R(:,tt) = [Ep*cos(T) - Sp ; Ep*sin(T) ; En - Sn] * sqrt(EC);        %±( n(E∥ - S∥) + s(E⊥cosθ - S⊥) + a(E⊥sinθ) ) * sqrt(Ec)
R(:,tt) = R(:,tt) * (En^2 + Ep^2 + Sn^2 + Sp^2 - 2*En*Sn - 2*Ep*Sp*cos(T))^(-3/4); %(E∥² + E⊥² + S∥² + S⊥² - 2E∥S∥ - 2E⊥S⊥cosθ )^(-3/4)
    
end

%Now to add the dimensions to take R from the n,s,a basis into x,y,z
n = NVori(:,1);
fs = cross(n,[1;0;0]);
fs = fs / norm(fs);     %s' == fraud s == fs
fa = cross(n,fs);

%Now we impliment some phi dependancy
PhiRes = 45;
Phi = linspace(0,2*pi,PhiRes);
R2 = zeros(3,ThetaRes,PhiRes);

for pp = 1:PhiRes
P = Phi(pp);

s = fs*cos(P) + fa*sin(P);
a = fs*sin(P) - fa*cos(P);

R2(:,:,pp) = R(1,:).*s + R(2,:).*a - R(3,:).*n;
%R2 is our final R value, with both θ & φ dependance

end

%Now for some graphs
figure(5)
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
axis square
hold on

%Firstly place our Main NV into the graph
scatter3(NVpos(1,1),NVpos(2,1),NVpos(3,1),40,[1,0.7,0],'filled','MarkerEdgeColor','black');
quiver3(NVpos(1,1),NVpos(2,1),NVpos(3,1),NVori(1,1),NVori(2,1),NVori(3,1),1e-8,'color',[1,0.7,0]);

%Now we have to place the particle that got changed
Type = Changes(1,2);
Index = Changes(2,2);
switch Type
    case 1
        scatter3(NVpos(1,Index),NVpos(2,Index),NVpos(3,Index),20,'g','filled','MarkerEdgeColor','black')
    case 2
        scatter3(P1pos(1,Index),P1pos(2,Index),P1pos(3,Index),20,'b','filled','MarkerEdgeColor','black')
    case 3
        scatter3(VCpos(1,Index),VCpos(2,Index),VCpos(3,Index),20,'r','filled','MarkerEdgeColor','black')
end

%Finally, we plot the R values
%But we need some colors to make it distinguishable
colors = interp1(linspace(0, 1, size(colormap(hsv), 1)), colormap(hsv), (1:PhiRes)/PhiRes);
for pp = 2:PhiRes
plot3(R2(1,:,pp),R2(2,:,pp),R2(3,:,pp), 'color', colors(pp,:))
end
plot3(R2(1,:,1),R2(2,:,1),R2(3,:,1), 'color', "k", 'LineWidth',2)
%test
%plot3(R(1,:),R(2,:),R(3,:))
axis square



hold off























