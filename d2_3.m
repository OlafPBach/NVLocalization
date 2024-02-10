% This is a test for PCA method using the longitudinal fields of NV centers
% arranged in a grid to locate a point charge above

clear
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
P1_num = 25;
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
LoopNum = P1_num*50;
P1_stat = (100*rand(P1_num,LoopNum) < 1*P1_occ);


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


Method = "13% Correlations";
if Method == "Weighted Correlations"
%% Weighted Correlations
TestNum = 120^2;
num_rows_cols = ceil(sqrt(TestNum));
[x, y] = meshgrid(linspace(0, 1, num_rows_cols), linspace(0, 1, num_rows_cols));
TestPos = [BoxDim(1)*x(:), BoxDim(2)*y(:)];
TestPos = TestPos(1:TestNum, :);
TestPos = [TestPos , BoxDim(3)*ones(TestNum,1)]';
%TestPos = [TestPos , zeros(TestNum,1)]';

SNR = zeros(TestNum,NV_num);
Weight = zeros(TestNum,NV_num);

for t = 1:TestNum
for vv = 1:NV_num
TestField = -EC*(TestPos(:,t) - NV_pos(:,vv)) ./ vecnorm(TestPos(:,t) - NV_pos(:,vv)).^3;
TestFieldLong(vv,t) = dot(TestField , NV_ori(:,vv));
STDenom = std(Elong(vv,:));
SNR(t,vv) = TestFieldLong(vv,t) / STDenom;
end
for vv = 1:NV_num
Weight(t,vv) = SNR(t,vv) / sum(SNR(t,:));
end
end



S = zeros(TestNum,LoopNum);
Test = zeros(1,TestNum);
for t = 1:TestNum

for ll = 1:LoopNum
%S(t,ll) = sum(Weight(t,:)' .* Elong(:,ll));
[TempR,TempP] = corrcoef(Elong(:,ll),Weight(t,:)'.*TestFieldLong(:,t));
%[TempR,TempP] = corrcoef(Elong(:,ll),TestFieldLong(:,t));
S(t,ll) = TempR(2,1);
if isnan(S(t,ll))
   S(t,ll) = 0;
end
end

%Test(t) = log10(std(abs(S(t,:))));
%Test(t) = log10( max(S(t,:))-min(S(t,:)) );
%Test(t) = log10(var(S(t,:)));
Test(t) = var(S(t,:));

if mod(nnz(Test), 1000) == 0
disp(floor((TestNum - nnz(Test))/1000))
end
end

colormap_custom = [0, 0, 1; 1, 0, 0]; % Blue to red
colors = interp1(linspace(min(Test), max(Test), size(colormap_custom, 1)), colormap_custom, Test);
hold on
scatter(TestPos(1,:),TestPos(2,:),50,colors,'filled')
scatter(P1_pos(1,:),P1_pos(2,:),50,'black','filled')
scatter(NV_pos(1,:),NV_pos(2,:),50,'green','filled')
hold off
axis square
legend('Variance Log10','Charge Trap Positions','NV Sensor Positions')
title('Charge Trap Localization')
subtitle("Max(red) = " + max(Test) + "  Min(blue) = " + min(Test));


elseif Method == "13% Correlations"
TestNum = 150^2;
num_rows_cols = ceil(sqrt(TestNum));
[x, y] = meshgrid(linspace(0, 1, num_rows_cols), linspace(0, 1, num_rows_cols));
TestPos = [BoxDim(1)*x(:), BoxDim(2)*y(:)];
TestPos = TestPos(1:TestNum, :);
TestPos = [TestPos , BoxDim(3)*ones(TestNum,1)]';

for t = 1:TestNum
for vv = 1:NV_num
TestField = -EC*(TestPos(:,t) - NV_pos(:,vv)) ./ vecnorm(TestPos(:,t) - NV_pos(:,vv)).^3;
TestFieldLong(t,vv) = dot(TestField , NV_ori(:,vv));
SNR(t,vv) = TestFieldLong(t,vv) / std(Elong(vv,:));
end
for vv = 1:NV_num
Weight(t,vv) = SNR(t,vv) / sum(SNR(t,:));
end
end

%Now break it down to only look at the strongest 13% of instances per NV

for t = 1:TestNum
%WeightedFields = abs(sum(Weight(t,:)'.*Elong(:,:),1));
WeightedFields = vecnorm(Weight(t,:)'.*Elong(:,:));
for sl = 1:floor(LoopNum*(P1_occ)/110)
    [~,A(sl)] = max(WeightedFields);
    new(:,sl) = Elong(:,A(sl));
    WeightedFields(A(sl)) = 0;
end

Test(t) = -norm(var(new,0,2));
if mod(nnz(Test), 1000) == 0
    disp(floor((TestNum - nnz(Test))/1000))
end
end

colormap_custom = [0, 0, 1; 1, 0, 0]; % Blue to red
colors = interp1(linspace(min(Test), max(Test), size(colormap_custom, 1)), colormap_custom, Test);
hold on
scatter(TestPos(1,:),TestPos(2,:),50,colors,'filled')
scatter(P1_pos(1,:),P1_pos(2,:),50,'black','filled')
scatter(NV_pos(1,:),NV_pos(2,:),50,'green','filled')
hold off
axis square
legend('Variance Log10','Charge Trap Positions','NV Sensor Positions')
title('Charge Trap Localization using Tom Method')
subtitle("Max(red) = " + max(Test) + "  Min(blue) = " + min(Test));


end