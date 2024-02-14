% This is a test for PCA method using the longitudinal fields of NV centers
% arranged in a grid to locate a point charge above

clear
clf
rng(14, 'twister')
c1 = datetime("now","Format","HH:mm:ss");
c = datetime("now","Format","HH:mm:ss");

% Constants %
P1occ = 13;                     %Suspected Precent that are charged negative
%Since we use a lot of E&M calculations, we need
DiaPerm = 5.66*8.854e-12;       %Permittivity in Diamond
EC = 1.602e-19/(4*pi*DiaPerm);  %Charge / 4*pi*epsilom

%Parameters
Nv = 5^2;                       %Total Number of NV s
Np = 10;                        %Total Number of P1 Charge Traps
Nl = 1000;                     %Total Number of Charge Instances
Nt = 150^2;                     %Total Number of Test Points
Ns = ceil(Nl*P1occ/100);        %Number of New Loops
BoxDim = 1e-7*[4;4;1];          %BoxDimensions
Indent = .20;                   %Precent indent


%% Charges, NV's, and Test Points
P1pos = [rand(1,Np)*BoxDim(1);rand(1,Np)*BoxDim(2);ones(1,Np)*BoxDim(3)]; %Our charges are xy random in the box and on the z roof
[~,sorter] = maxk(P1pos(2,:),10);
P1pos = P1pos(:,sorter);   %Sort our p1 Charge Indexes based on y position


%Place NVs
NVrow = ceil(sqrt(Nv));
IndentSize = Indent/2*[BoxDim(1),BoxDim(2)];
[y, x] = meshgrid(linspace(1, 0, NVrow), linspace(0, 1, NVrow));
NVpos = [IndentSize(1)+BoxDim(1)*(1-Indent)*x(:), IndentSize(2)+BoxDim(2)*(1-Indent)*y(:)];
NVpos = NVpos(1:Nv, :);
NVpos = [NVpos , zeros(Nv,1)]';
% The positions are in an xy grid all at the floor

%n_poss = 1/sqrt(3)*[1,1,-1,-1;1,-1,1,-1;1,-1,-1,1]; %1/sqrt(3) is the normalization factor)
%NV_ori = n_poss(:,randi([1,4],1,Nv));
NVori = repmat([0;0;1],1,Nv);

%place test points
Testrow = ceil(sqrt(Nt));
[y, x] = meshgrid(linspace(1, 0, Testrow), linspace(0, 1, Testrow));
Testpos = [IndentSize(1)+BoxDim(1)*(1-Indent)*x(:), IndentSize(2)+BoxDim(2)*(1-Indent)*y(:)];
Testpos = Testpos(1:Nt, :);
Testpos = [Testpos , BoxDim(3)*ones(Nt,1)]';
%TestPos = [TestPos , zeros(Nt,1)]';


%% Electric Fields
% Randomize charge instance per loop
P1stat = (100*rand(Np,Nl) < P1occ);

%Now we look at how each loop effects the percieved electric fields at the
%point of each of our NV's

Efield = zeros(3,Nl,Nv);
for vv = 1:Nv
for ll = 1:Nl
%Add contributions from our P1 Traps & Vacancies
Efield(:,ll,vv) = Efield(:,ll,vv) - EC*sum((P1pos(:,P1stat(:,ll)) - NVpos(:,vv)) ./ vecnorm(P1pos(:,P1stat(:,ll)) - NVpos(:,vv)).^3,2);

end %ll loop
end %vv loop

%% Longitudinal Fields and format
%In our actual readouts we know the position of the spectral lines are the
%sum and difference of the projection and rejection fields respect to the
%NV orientation. We call the vector projection the Longitudinal field, and
%the rejection field the transversal field.
DataMatrix = zeros(Nv,Nl);
for vv = 1:Nv
for ll = 1:Nl
DataMatrix(vv,ll) = dot(Efield(:,ll,vv) , NVori(:,vv));
end
end

%We also need to know what a suspected charge at every test point would
%look like
TestLong = zeros(Nv,Nt);
for tt = 1:Nt
for vv = 1:Nv
    TestField = -EC*(Testpos(:,tt) - NVpos(:,vv)) ./ vecnorm(Testpos(:,tt) - NVpos(:,vv)).^3;
    TestLong(vv,tt) = dot(TestField , NVori(:,vv));
end
end



Method = "Tom";
if Method == "Tom"
%% Tom Method
sigsqr = var(DataMatrix,0,2);
%Varience along the instances dimension

S = zeros(Nt,Nl);
for tt = 1:Nt
    S(tt,:) = sum(TestLong(:,tt)./sigsqr(:).*DataMatrix(:,:),1) / sum(TestLong(:,tt).^2 ./sigsqr(:));
    %First indexes are all vv
end

Test = var(S,0,2);
TestName = "Variance";



elseif Method == "Corr2"
%% Correlation 2
S = zeros(Nt,Nl);
Test = zeros(1,Nt);
for tt = 1:Nt

for ll = 1:Nl
[TempR,~] = corrcoef(DataMatrix(:,ll),TestLong(:,tt));
S(tt,ll) = TempR(2,1);
if isnan(S(tt,ll))
   S(tt,ll) = 0;
end
end

Test(tt) = var(S(tt,:));
TestName = "Variance of Correlations";
if mod(nnz(Test), 1000) == 0
disp(datetime("now","Format","HH:mm:ss") - c)
disp(floor((Nt - nnz(Test))/1000))
c = datetime("now","Format","HH:mm:ss");
end
end




end

c2 = datetime("now","Format","HH:mm:ss");
disp("Total Time Taken")
disp(c2-c1)

%% Graphing %%
figure(1)
hold on

colormap_custom = [0, 0, 1; 1, 0, 0]; % Blue to red
colors = interp1(linspace(min(Test), max(Test), size(colormap_custom, 1)), colormap_custom, Test);

scatter(Testpos(1,:),Testpos(2,:),50,colors,'filled')
scatter(P1pos(1,:),P1pos(2,:),50,'black','filled')
scatter(NVpos(1,:),NVpos(2,:),50,'green','filled')
%hold off

axis square
legend(TestName,'Charge Trap Positions','NV Sensor Positions')
title('Charge Trap Localization')
subtitle("Max(red) = " + max(Test) + "  Min(blue) = " + min(Test));
%scatter(P1pos(1,P1stat(:,53)),P1pos(2,P1stat(:,53)),50,'yellow','filled')

%% Analysis
figure(2)
tiledlayout(2,3);

nexttile
hold on
[~,t1] = max(Test);
[~,Index1] = maxk(S(t1,:),Nl);

plot(S(t1,Index1),'black')
scatter(linspace(0,Nl,Nl),2.5*P1stat(1,Index1),5,'filled','r')
scatter(linspace(0,Nl,Nl),2.4*P1stat(2,Index1),5,'filled','y')
scatter(linspace(0,Nl,Nl),2.3*P1stat(3,Index1),5,'filled','g')
scatter(linspace(0,Nl,Nl),2.2*P1stat(4,Index1),5,'filled','c')
scatter(linspace(0,Nl,Nl),2.1*P1stat(5,Index1),5,'filled','b')
ylim([0 5])
legend('S Sorted','Charge1 State','Charge2 State','Charge3 State','Charge4 State','Charge5 State')

nexttile
hold on
[~,t2] = min(vecnorm(P1pos(:,4)-Testpos(:,:),1));
[~,Index2] = maxk(S(t2,:),Nl);
plot(S(t2,Index2),'black')
scatter(linspace(0,Nl,Nl),2.5*P1stat(1,Index2),5,'filled','r')
scatter(linspace(0,Nl,Nl),2.4*P1stat(2,Index2),5,'filled','y')
scatter(linspace(0,Nl,Nl),2.3*P1stat(3,Index2),5,'filled','g')
scatter(linspace(0,Nl,Nl),2.2*P1stat(4,Index2),5,'filled','c')
scatter(linspace(0,Nl,Nl),2.1*P1stat(5,Index2),5,'filled','b')
ylim([0 5])
legend('S Sorted','Charge1 State','Charge2 State','Charge3 State','Charge4 State','Charge5 State')

nexttile
hold on
[~,t3] = min(vecnorm([2.63e-7;2.35e-7;0]-Testpos(:,:),1));
[~,Index3] = maxk(S(t3,:),Nl);
plot(S(t3,Index3),'black')
scatter(linspace(0,Nl,Nl),2.5*P1stat(2,Index3),5,'filled','r')
scatter(linspace(0,Nl,Nl),2.4*P1stat(4,Index3),5,'filled','y')
scatter(linspace(0,Nl,Nl),2.3*P1stat(6,Index3),5,'filled','g')
scatter(linspace(0,Nl,Nl),2.2*P1stat(8,Index3),5,'filled','c')
ylim([0 5])
legend('S Sorted','Charge2 State','Charge4 State','Charge6 State','Charge8 State')

nexttile
hold on
[~,t4] = min(vecnorm(P1pos(:,7)-Testpos(:,:),1));
[~,Index4] = maxk(S(t4,:),Nl);
plot(S(t4,Index4),'black')
scatter(linspace(0,Nl,Nl),2.5*P1stat(5,Index4),5,'filled','r')
scatter(linspace(0,Nl,Nl),2.4*P1stat(7,Index4),5,'filled','y')
scatter(linspace(0,Nl,Nl),2.3*P1stat(8,Index4),5,'filled','g')
scatter(linspace(0,Nl,Nl),2.2*P1stat(10,Index4),5,'filled','c')
ylim([0 5])
legend('S Sorted','Charge5 State','Charge7 State','Charge8 State','Charge10 State')

nexttile
hold on
[~,t5] = min(vecnorm([2.011e-7;6.792e-8;0]-Testpos(:,:),1));
[~,Index5] = maxk(S(t5,:),Nl);
plot(S(t5,Index5),'black')
scatter(linspace(0,Nl,Nl),2.5*P1stat(4,Index5),5,'filled','r')
scatter(linspace(0,Nl,Nl),2.4*P1stat(8,Index5),5,'filled','y')
scatter(linspace(0,Nl,Nl),2.3*P1stat(9,Index5),5,'filled','g')
scatter(linspace(0,Nl,Nl),2.2*P1stat(10,Index5),5,'filled','c')
ylim([0 5])
legend('S Sorted','Charge4 State','Charge8 State','Charge9 State','Charge10 State')

nexttile
hold on
[~,t6] = min(vecnorm([3.6e-7;1.6e-7;0]-Testpos(:,:),1));
[~,Index6] = maxk(S(t6,:),Nl);
plot(S(t6,Index6),'black')
scatter(linspace(0,Nl,Nl),2.5*P1stat(6,Index6),5,'filled','r')
scatter(linspace(0,Nl,Nl),2.4*P1stat(8,Index6),5,'filled','y')
scatter(linspace(0,Nl,Nl),2.3*P1stat(9,Index6),5,'filled','g')
ylim([0 5])
legend('S Sorted','Charge6 State','Charge8 State','Charge9 State')



%figure(2)
%hold on
%plot(S(1551,:))
%plot(P1stat(10,:))

%figure(3)
%hold on
%plot(DataMatrix(Index(1551,:),:)',"black")
%plot((TestLong(Index(1551,:),1551)*P1stat(10,:))',"red")



