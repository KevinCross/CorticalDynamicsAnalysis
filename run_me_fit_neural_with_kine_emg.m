%this script regresses the emg and kinematics signals on to the MC rotational dynamics
close ALL

% monkey can be either: Ayur, Patch, Puran or Xander
monkey='Patch';
proj_name=strcat('figures/MC_data/',monkey,'_Data.mat/Projection_Neural');
summ_name=strcat('figures/MC_data/',monkey,'_Data.mat/Summary_Neural');
kine_name=strcat('data_kinematics/',monkey,'_Joint_Data_Set');
EMG_name=strcat('data_EMG/',monkey,'_EMG_Data_set');

%% load data
load(proj_name)
load(summ_name)

%load kinematic and emg data and store in Data_Sens structure
tmp_kine=load(kine_name);
tmp_emg=load(EMG_name);
Data_Sens=struct;
for i=1:8
    Data_Sens(i).A=[tmp_kine.Data_struct(i).A(1:30,:),tmp_emg.Data_struct(i).A];
    Data_Sens(i).times=tmp_emg.Data_struct(i).times;
end

%% get data into matrix form and fit with linear equation
%matrices to hold the neural and sensory data
mat_neur=[];
for i=1:length(Proj_Neural)
    mat_neur=[mat_neur,Proj_Neural(i).proj.'];
end
[mat_sens]=preprocessing(Data_Sens);
mat_sens=mat_sens.';
%linear fit
W=mat_neur/mat_sens;

%% calculate the VAF by the fit
total_VAR=norm(mat_neur,'fro');
mat_neural_pred=W*mat_sens;  %predicted variance by fit
resd=norm(mat_neur-mat_neural_pred,'fro'); %residual error btwn predicted and observed
VAF=1-resd.^2/total_VAR.^2;   %VAF

VAF_plane=[];
for i=1:2:5
    total_VAR=norm(mat_neur(i:i+1,:),'fro');
    mat_neural_pred=W*mat_sens;  %predicted variance by fit
    resd=norm(mat_neur(i:i+1,:)-mat_neural_pred(i:i+1,:),'fro'); %residual error btwn predicted and observed
    VAF_plane=[VAF_plane 1-resd.^2/total_VAR.^2];   %VAF
end

%% plot jPCA planes for original and fitted data
Projection_Data_Pred=Proj_Neural; %create predicted Projection structure 
%replace observed projection values with predicted projection values
for i=1:length(Proj_Neural)
    Projection_Data_Pred(i).proj=mat_neural_pred(:,1+(i-1)*30:30*i).';
end

%plot planes
plotParams.planes2plot=[1,2,3];
phaseSpace(Proj_Neural,Summary_Neural,plotParams);
phaseSpace(Projection_Data_Pred,Summary_Neural,plotParams);
h=figure(4);
suptitle(strcat('Model R2 value=',num2str(VAF),' Plane R2 value=',num2str(VAF_plane(1))))
h=figure(5);
suptitle(strcat('Model R2 value=',num2str(VAF),' Plane R2 value=',num2str(VAF_plane(2))))
h=figure(6);
suptitle(strcat('Model R2 value=',num2str(VAF),' Plane R2 value=',num2str(VAF_plane(3))))

dir_name=strcat('figures/rotation_pred_from_neural/',monkey);
mkdir(dir_name)
for i=1:3
    h=figure(i);
    file_name=strcat(dir_name,'/jPCA_plot_original_',num2str(i));
    saveas(h,file_name)
    saveas(h,file_name,'epsc')
end
for i=4:6
    h=figure(i);
    file_name=strcat(dir_name,'/jPCA_plot_predicted_',num2str(i-3));
    saveas(h,file_name)
    saveas(h,file_name,'epsc')
end



function [bigA]=preprocessing(Data)
    %processing to follow same normalization and mean-subtraction procedure
    %found in jPCA_v2 applied to Data_Sens
    softenNorm=0.0000;
    numConds = length(Data);
    bigA = vertcat(Data.A);
    ranges = range(bigA);  % For each neuron, the firing rate range across all conditions and times.
    normFactors = (ranges+softenNorm);
    bigA = bsxfun(@times, bigA, 1./normFactors);  % normalize
    
    sumA = 0;
    for c = 1:numConds
        sumA = sumA + bsxfun(@times, Data(c).A, 1./normFactors);  % using the same normalization as above
    end
    meanA = sumA/numConds;
    bigA = bigA-repmat(meanA,numConds,1);
    
    
end