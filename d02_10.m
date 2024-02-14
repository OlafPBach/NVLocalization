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
Nl = Np*20;                     %Total Number of Charge Instances
Nt = 98^2;                      %Total Number of Test Points
Ns = ceil(Nl*P1occ/100);        %Number of New Loops
BoxDim = 1e-7*[4;4;1];          %BoxDimensions
Indent = .20;                   %Precent indent


%% Charges, NV's, and Test Points
P1pos = [rand(1,Np)*BoxDim(1);rand(1,Np)*BoxDim(2);ones(1,Np)*BoxDim(3)]; 
%Our charges are xy random in the box and on the z roof

%Place NVs
NVrow = ceil(sqrt(Nv));
IndentSize = Indent/2*[BoxDim(1),BoxDim(2)];
[x, y] = meshgrid(linspace(0, 1, NVrow), linspace(0, 1, NVrow));
NVpos = [IndentSize(1)+BoxDim(1)*(1-Indent)*x(:), IndentSize(2)+BoxDim(2)*(1-Indent)*y(:)];
NVpos = NVpos(1:Nv, :);
NVpos = [NVpos , zeros(Nv,1)]';
% The positions are in an xy grid all at the floor

%n_poss = 1/sqrt(3)*[1,1,-1,-1;1,-1,1,-1;1,-1,-1,1]; %1/sqrt(3) is the normalization factor)
%NV_ori = n_poss(:,randi([1,4],1,Nv));
NVori = repmat([0;0;1],1,Nv);

%place test points
NVrow = ceil(sqrt(Nt));
[x, y] = meshgrid(linspace(0, 1, NVrow), linspace(0, 1, NVrow));
TESTpos = [IndentSize(1)+BoxDim(1)*(1-Indent)*x(:), IndentSize(2)+BoxDim(2)*(1-Indent)*y(:)];
TESTpos = TESTpos(1:Nt, :);
TESTpos = [TESTpos , BoxDim(3)*ones(Nt,1)]';
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
    TestField = -EC*(TESTpos(:,tt) - NVpos(:,vv)) ./ vecnorm(TESTpos(:,tt) - NVpos(:,vv)).^3;
    TestLong(vv,tt) = dot(TestField , NVori(:,vv));
end
end



Method = "Corr2+";
if Method == "Tom"
%% Tom Method
sigsqr = var(DataMatrix,0,2);
%Varience along the instances dimension

S = zeros(Nt,Nl);
for tt = 1:Nt
    S(tt,:) = sum(TestLong(:,tt)./sigsqr(:).*DataMatrix(:,:),1) / sum(TestLong(:,tt).^2 ./sigsqr(:));
    %First indexes are all vv
end

switch 1
    case 1
    Test = var(S,0,2);
    TestName = "Variance";
    case 2
    Test =  max(S,2)-min(S,2) ;
    TestName = "Max - Min";
    case 3
    Test = std(abs(S),2);
    TestName = "Standard Deviation";
end

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

elseif Method == "Corr2+"
%% Correlation 2
S = zeros(Nt,Nl);
Test = zeros(1,Nt);
Index = zeros(Nt,4);
for tt = 1:Nt

dist = vecnorm(TESTpos(:,tt) - NVpos(:,:),2);
[~,Index(tt,:)] = mink(dist,4);

for ll = 1:Nl
[TempR,~] = corrcoef(DataMatrix(Index(tt,:),ll),TestLong(Index(tt,:),tt));
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

elseif Method == "WeightCorr"
%% Weighted Correlations
SNR = zeros(Nt,Nv);
Weight = zeros(Nt,Nv);
for tt = 1:Nt
for vv = 1:Nv
STDenom = std(DataMatrix(vv,:));
SNR(tt,vv) = TestLong(vv,tt) / STDenom;
end
for vv = 1:Nv
Weight(tt,vv) = SNR(tt,vv) / sum(SNR(tt,:));
end
end

S = zeros(Nt,Nl);
Test = zeros(1,Nt);
for tt = 1:Nt

for ll = 1:Nl
%S(t,ll) = sum(Weight(t,:)' .* Elong(:,ll));
[TempR,TempP] = corrcoef(DataMatrix(:,ll),Weight(tt,:)'.*TestLong(:,tt));
%[TempR,TempP] = corrcoef(Elong(:,ll),TestFieldLong(:,t));
S(tt,ll) = TempR(2,1);
if isnan(S(tt,ll))
   S(tt,ll) = 0;
end
end

Test(tt) = var(S(tt,:));
TestName = "Varience of Wgt Correlations";
if mod(nnz(Test), 1000) == 0
disp(datetime("now","Format","HH:mm:ss") - c)
disp(floor((Nt - nnz(Test))/1000))
c = datetime("now","Format","HH:mm:ss");
end
end


elseif Method == "Olaf"
%% Olaf Method
Test = zeros(1,Nt);
S = zeros(Nt,Nl);
S2 = zeros(4,Nt);
Index = zeros(Nt,4);
Index2 = zeros(Nt,Ns);
for tt = 1:Nt
%Find the first 4 closest NV's to each test point
dist = vecnorm(TESTpos(:,tt) - NVpos(:,:),2);
[~,Index(tt,:)] = mink(dist,4);

for ll = 1:Nl
[TempR,~] = corrcoef(DataMatrix(Index(tt,:),ll),TestLong(Index(tt,:),tt));
if isnan(TempR(2,1))
    S(tt,ll) = 0;
else
    Ratio = abs( norm( TestLong(Index(tt,:),tt) ./ DataMatrix(Index(tt,:),ll) ) );
    if Ratio < 1
    S(tt,ll) = TempR(2,1) * Ratio;
    else
    S(tt,ll) = TempR(2,1) / Ratio;   
    end
%S(tt,ll) = TempR(2,1);
end
end

[~,Index2(tt,:)] = maxk(S(tt,:),Ns);
NewField = DataMatrix(Index(tt,:),Index2(tt,:));

S2(:,tt) = var(NewField,0,2);

Test(tt) = sum(S2(:,tt));
TestName = "Sum of Variance";
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
colormap_custom = [0, 0, 1; 1, 0, 0]; % Blue to red
colors = interp1(linspace(min(Test), max(Test), size(colormap_custom, 1)), colormap_custom, Test);
hold on
scatter(TESTpos(1,:),TESTpos(2,:),50,colors,'filled')
scatter(P1pos(1,:),P1pos(2,:),50,'black','filled')
scatter(NVpos(1,:),NVpos(2,:),50,'green','filled')
%hold off
axis square
legend(TestName,'Charge Trap Positions','NV Sensor Positions')
title('Charge Trap Localization')
subtitle("Max(red) = " + max(Test) + "  Min(blue) = " + min(Test));
%scatter(P1pos(1,P1stat(:,53)),P1pos(2,P1stat(:,53)),50,'yellow','filled')



%figure(2)
%hold on
%plot(S(1551,:))
%plot(P1stat(10,:))

%figure(3)
%hold on
%plot(DataMatrix(Index(1551,:),:)',"black")
%plot((TestLong(Index(1551,:),1551)*P1stat(10,:))',"red")



