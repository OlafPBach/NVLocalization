% This is a test for PCA method using the longitudinal fields of NV centers
% arranged in a grid to locate a point charge above

clear
clf
rng(14, 'twister')

%% Constants %%
%Since we use a lot of E&M calculations, we need
DiaPerm = 5.66*8.854e-12;       %Permittivity in Diamond
EC = 1.602e-19/(4*pi*DiaPerm);  %Charge / 4*pi*epsilom

%We only need p1 density as we'll say the NV's are always neutral and all
%charges are P1s. Also no density as all charges will be placed randomly
P1_occ = 13;            %Suspected Precent that are charged negative
BoxDim = 1e-7*[4;4;1];

%% Charges and NV's
P1_num = 250;
P1_pos = [rand(1,P1_num)*BoxDim(1);rand(1,P1_num)*BoxDim(2);ones(1,P1_num)*BoxDim(3)]; 
%P1_pos = [rand(1,P1_num)*BoxDim(1);rand(1,P1_num)*BoxDim(2);zeros(1,P1_num)]; 

%Our one charge is xy random in the box and on the z roof

NV_num = 5^2;
num_rows_cols = ceil(sqrt(NV_num));
[x, y] = meshgrid(linspace(0, 1, num_rows_cols), linspace(0, 1, num_rows_cols));
NV_pos = [BoxDim(1)*x(:), BoxDim(2)*y(:)];
NV_pos = NV_pos(1:NV_num, :);
NV_pos = [NV_pos , zeros(NV_num,1)]';
% The positions are in an xy grid all at the floor

n_poss = 1/sqrt(3)*[1,1,-1,-1;1,-1,1,-1;1,-1,-1,1]; %1/sqrt(3) is the normalization factor)
%NV_ori = n_poss(:,randi([1,4],1,NV_num));
NV_ori = repmat([0;0;1],1,NV_num);

% All NV's point along crystalline axis

%% Randomizing Charge Per loop
%As a test we only have 1 charge but still we need this
LoopNum = 3;
P1_stat = (100*rand(P1_num,LoopNum) < P1_occ);


%% Electric Fields
%Now we look at how each loop effects the percieved electric fields at the
%point of each of our NV's

Efield = zeros(3,LoopNum,NV_num);
for vv = 1:NV_num
for ll = 1:LoopNum
%Add contributions from our P1 Traps & Vacancies
Efield(:,ll,vv) = Efield(:,ll,vv) - EC*sum((P1_pos(:,P1_stat(:,ll)) - NV_pos(:,vv)) ./ vecnorm(P1_pos(:,P1_stat(:,ll)) - NV_pos(:,vv)).^3,2);

end %ll loop
end %vv loop

%% Longitudinal Fields and format
%In our actual readouts we know the position of the spectral lines are the
%sum and difference of the projection and rejection fields respect to the
%NV orientation. We call the vector projection the Longitudinal field, and
%the rejection field the transversal field.
Elong = zeros(NV_num,LoopNum);
for vv = 1:NV_num
for ll = 1:LoopNum
Elong(vv,ll) = dot(Efield(:,ll,vv) , NV_ori(:,vv));
end
end

TestNum = 150^2;
num_rows_cols = ceil(sqrt(TestNum));
[x, y] = meshgrid(linspace(0, 1, num_rows_cols), linspace(0, 1, num_rows_cols));
TestPos = [BoxDim(1)*x(:), BoxDim(2)*y(:)];
TestPos = TestPos(1:TestNum, :);
TestPos = [TestPos , BoxDim(3)*ones(TestNum,1)]';
%TestPos = [TestPos , zeros(TestNum,1)]';

%% Tom Method
figure(1)
SNR = zeros(TestNum,NV_num);
Weight = zeros(TestNum,NV_num);

for t = 1:TestNum
for vv = 1:NV_num
TestField = -EC*(TestPos(:,t) - NV_pos(:,vv)) ./ vecnorm(TestPos(:,t) - NV_pos(:,vv)).^3;
TestFieldLong = dot(TestField , NV_ori(:,vv));
STDenom = std(Elong(vv,:));
SNR(t,vv) = TestFieldLong / STDenom;
end
for vv = 1:NV_num
Weight(t,vv) = SNR(t,vv) / sum(SNR(t,:));
end
end

S = zeros(TestNum,LoopNum);
Variance = zeros(1,TestNum);
for t = 1:TestNum
for ll = 1:LoopNum
S(t,ll) = sum(Weight(t,:)' .* Elong(:,ll));
end

%Variance(t) = log10(std(abs(S(t,:))));
%Variance(t) = log10(max(abs(S(t,:))));
Variance(t) = var(S(t,:));
if mod(nnz(Variance), 1000) == 0
disp(floor((TestNum - nnz(Variance))/1000))
end
end

colormap_custom = [0, 0, 1; 1, 0, 0]; % Blue to red
colors = interp1(linspace(min(Variance), max(Variance), size(colormap_custom, 1)), colormap_custom, Variance);
hold on
scatter(TestPos(1,:),TestPos(2,:),50,colors,'filled')
scatter(P1_pos(1,:),P1_pos(2,:),50,'black','filled')
scatter(NV_pos(1,:),NV_pos(2,:),50,'green','filled')
hold off
axis square
legend('Variance','Charge Trap Positions','NV Sensor Positions')
title('Charge Trap Localization using Tom Method')
subtitle("Max(red) = " + max(Variance) + "  Min(blue) = " + min(Variance));

%% Correlation 1 %%
figure(2)
LongField = zeros(NV_num,1);
for t = 1:TestNum
for vv = 1:NV_num
%Add contributions from our P1 Traps & Vacancies
TestField = -EC*sum((TestPos(:,t) - NV_pos(:,vv)) ./ vecnorm(TestPos(:,t) - NV_pos(:,vv)).^3,2);
LongField(vv,1) = dot(TestField , NV_ori(:,vv));
end %vv loop

for ll = 1:LoopNum
[TempR,TempP] = corrcoef(Elong(:,ll),LongField);

HeatMapR(t,ll) = TempR(2,1);
HeatMapP(t,ll) = log10(abs(TempP(2,1))+1e-6);

end
end

colormap_custom = [0, 0, 1; 1, 0, 0]; % Blue to red
figure(1)
tiledlayout(2,3)

for ll = 1:LoopNum
nexttile
colors = interp1(linspace(-1, 1, size(colormap_custom, 1)), colormap_custom, HeatMapR(:,ll));
hold on
scatter(TestPos(1,:),TestPos(2,:),50,colors,'filled')
scatter(P1_pos(1,P1_stat(:,ll)),P1_pos(2,P1_stat(:,ll)),50,'black','filled')
scatter(NV_pos(1,:),NV_pos(2,:),20,'green','filled')
hold off
axis square
legend('Correlations','Charge Trap Positions','NV Sensor Positions')
title('Charge Trap Localization using Tom Method')
subtitle("Max(red) = " + max(HeatMapR(:,ll)) + "  Min(blue) = " + min(HeatMapR(:,ll)));
end


