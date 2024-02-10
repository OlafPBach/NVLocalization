% This is a test for PCA method using the longitudinal fields of NV centers
% arranged in a grid to locate a point charge above

clear
rng(14, 'twister')
c1 = clock;

%% Constants %%
%Since we use a lot of E&M calculations, we need
DiaPerm = 5.66*8.854e-12;       %Permittivity in Diamond
EC = 1.602e-19/(4*pi*DiaPerm);  %Charge / 4*pi*epsilom

%We only need p1 density as we'll say the NV's are always neutral and all
%charges are P1s. Also no density as all charges will be placed randomly
P1_occ = 13;            %Suspected Precent that are charged negative
BoxDim = 1e-7*[4;4;1];

%% Charges and NV's
P1_num = 5;
P1_pos = [rand(1,P1_num)*BoxDim(1);rand(1,P1_num)*BoxDim(2);ones(1,P1_num)*BoxDim(3)]; 
%P1_pos = [rand(1,P1_num)*BoxDim(1);rand(1,P1_num)*BoxDim(2);zeros(1,P1_num)]; 

%Our one charge is xy random in the box and on the z roof

NV_num = 25^2;
Indent = .20; %Precent indent
num_rows_cols = ceil(sqrt(NV_num));
IndentSize = Indent/2*[BoxDim(1),BoxDim(2)];
[x, y] = meshgrid(linspace(0, 1, num_rows_cols), linspace(0, 1, num_rows_cols));
NV_pos = [IndentSize(1)+BoxDim(1)*(1-Indent)*x(:), IndentSize(2)+BoxDim(2)*(1-Indent)*y(:)];
NV_pos = NV_pos(1:NV_num, :);
NV_pos = [NV_pos , zeros(NV_num,1)]';
% The positions are in an xy grid all at the floor

n_poss = 1/sqrt(3)*[1,1,-1,-1;1,-1,1,-1;1,-1,-1,1]; %1/sqrt(3) is the normalization factor)
%NV_ori = n_poss(:,randi([1,4],1,NV_num));
NV_ori = repmat([0;0;1],1,NV_num);

% All NV's point along crystalline axis

%% Randomizing Charge Per loop
LoopNum = P1_num*20;
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

TestNum = 100^2;
num_rows_cols = ceil(sqrt(TestNum));
[x, y] = meshgrid(linspace(0, 1, num_rows_cols), linspace(0, 1, num_rows_cols));
TestPos = [IndentSize(1)+BoxDim(1)*(1-Indent)*x(:), IndentSize(2)+BoxDim(2)*(1-Indent)*y(:)];
TestPos = TestPos(1:TestNum, :);
TestPos = [TestPos , BoxDim(3)*ones(TestNum,1)]';
%TestPos = [TestPos , zeros(TestNum,1)]';

Method = "Olaf3";
if Method == "Tom"
%% Tom Method
TestFieldLong = zeros(NV_num,TestNum);
for t = 1:TestNum
for vv = 1:NV_num
    TestField = -EC*(TestPos(:,t) - NV_pos(:,vv)) ./ vecnorm(TestPos(:,t) - NV_pos(:,vv)).^3;
    TestFieldLong(vv,t) = dot(TestField , NV_ori(:,vv));
end
end

sigsqr = var(Elong,0,2);
%Varience along the instances dimension

S = zeros(TestNum,LoopNum);
for t = 1:TestNum
    S(t,:) = sum(TestFieldLong(:,t)./sigsqr(:).*Elong(:,:),1) / sum(TestFieldLong(:,t).^2 ./sigsqr(:));
    %First indexes are all vv
end

Test = var(S,0,2);


%Test(t) = log10(std(abs(S(t,:))));
%Test(t) =  max(S(t,:))-min(S(t,:)) ;
%Test(t) = log10(var(S(t,:)));
%Test(t) = var(S(t,:));


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



elseif Method == "Olaf1"
%% Olaf Method

TestFieldLong = zeros(NV_num,TestNum);
for t = 1:TestNum
for vv = 1:NV_num
    TestField = -EC*(TestPos(:,t) - NV_pos(:,vv)) ./ vecnorm(TestPos(:,t) - NV_pos(:,vv)).^3;
    TestFieldLong(vv,t) = dot(TestField , NV_ori(:,vv));
end
end

NewNum = floor(LoopNum*P1_occ/101);
S = zeros(TestNum,LoopNum);
Test = zeros(TestNum,1);
for t = 1:TestNum
for ll = 1:LoopNum
    [TempR,~] = corrcoef(Elong(:,ll),TestFieldLong(:,t));
    Base(ll) = TempR(2,1);
end
for sl = 1:NewNum
    [~,A] = max(Base);
    Base(A) = 0;
    NewLong(:,sl) = Elong(:,A);
end
    Test(t) = norm(var(NewLong,0,2));
    if mod(nnz(Test), 1000) == 0
    c = clock;
    disp(floor((TestNum - nnz(Test))/1000)+ "  " + c(4) + ":" + c(5))
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
legend('Norm of Varience','Charge Trap Positions','NV Sensor Positions')
title('Charge Trap Localization')
subtitle("Max(red) = " + max(Test) + "  Min(blue) = " + min(Test));

elseif Method == "Olaf2"
%% Olaf 2
TestFieldLong = zeros(NV_num,TestNum);
for t = 1:TestNum
for vv = 1:NV_num
    TestField = -EC*(TestPos(:,t) - NV_pos(:,vv)) ./ vecnorm(TestPos(:,t) - NV_pos(:,vv)).^3;
    TestFieldLong(vv,t) = dot(TestField , NV_ori(:,vv));
end
end

NewNum = ceil(LoopNum*P1_occ/101);
for t = 1:TestNum
CheckField = vecnorm(Elong - TestFieldLong(:,t));
for sl = 1:NewNum
[~,A] = min(CheckField);
CheckField(A) = NaN;
NewField(:,sl) = Elong(:,A);
end
Test(t) = var(NewField,0,"all");
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


elseif Method == "Olaf3"
%% Olaf 3
TestFieldLong = zeros(NV_num,TestNum);
for t = 1:TestNum
for vv = 1:NV_num
    TestField = -EC*(TestPos(:,t) - NV_pos(:,vv)) ./ vecnorm(TestPos(:,t) - NV_pos(:,vv)).^3;
    TestFieldLong(vv,t) = dot(TestField , NV_ori(:,vv));
end
end

NewNum = floor(LoopNum*P1_occ/101);
S = zeros(TestNum,LoopNum);
Test = zeros(TestNum,1);
for t = 1:TestNum
    [TempR,~] = corrcoef(Elong(:,:),TestFieldLong(:,t)*ones(1,LoopNum));
    Base = TempR(2,1);
    Compare = zeros(1,LoopNum);

    for ll = 1:LoopNum
        [TempR,~] = corrcoef(Elong(:,[1:ll-1,ll+1:end]),TestFieldLong(:,t)*ones(1,LoopNum-1));
        Compare(ll) = TempR(2,1) - Base;
    end
    
    NewLong = zeros(NV_num,NewNum);
    for sl = 1:NewNum
        [~,A] = min(Compare);
    NewLong(:,sl) = Elong(:,A);
    Compare(A) = 1;
    end

    [TempR,~] = corrcoef(NewLong(:,:),TestFieldLong(:,t)*ones(1,NewNum));
    Test(t) = TempR(2,1);

    if mod(nnz(Test), 1000) == 0
    c = clock;
    disp(floor((TestNum - nnz(Test))/1000)+ "  " + c(4) + ":" + c(5))
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
c2 = clock;
TimeDiff = c2-c1;
disp("Total time taken " + TimeDiff(4)*60 + TimeDiff(5) + " minutes")


