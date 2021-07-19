% this script applies jPCA analysis to the kinematic and muscle signals
% separately for each monkey
% calculates a null distribution based on down-sampling neural activity to
% match the number of kinematic/muscle signals

close ALL
%which monkey to analyze: can be Ayur, Manny, Patch, Puran or Xander
monkey_name='Ayur';
EMG_flag=1; %1= analyze emg, 0= don't analyze emg; set = 1 for all monkeys except Manny as no EMG was collected
numb_it=1000; %number of samples to generate the null distribution


%% import data
[KINE_data,M1_data,EMG_data]=import_data(monkey_name, EMG_flag);

%% select times of interest for kinematic data (emg and neural signals are already trimmed to the correct time epoch)
KINE_data_temp=struct;
for j=1:8
    KINE_data_temp(1,j).A =KINE_data(1,j).A(1:30,:);
end

%% bstrap data
[Proj_Kine,Summ_Kine]=bstrap_neural_signals(KINE_data_temp,M1_data,4,numb_it);
if EMG_flag
    [Proj_EMG,Summ_EMG]=bstrap_neural_signals(EMG_data,M1_data,6,numb_it);
end



function [Projection,Summary_data] = bstrap_neural_signals(data,M1_data,ndim,niters)
    %this function performs jPCA analysis on the data (kinematic or emg) and compares the
    %results to a null distribution generated from the neural data

    %assume data and M1_data are structs with fieldnames A and times
    % A should be timepoints x signals where signals is the number of
    % kinematic/muscle signals or neurons
    % A should also contain only the timepoints that are to be analzyed (no
    % preperturbation/reach preparation or steady-state activity)

    %output array of R2 from bstrap
    R2_bstrap_skew=zeros(1,niters);
    R2_bstrap_best=zeros(1,niters);
    %get number of signals in data
    [~,nsig]=size(data(1).A);
    [~,nneurons]=size(M1_data(1).A);


    %perform jPCA on kinematic/muscle signals and get the R2
    jPCA_params.softenNorm=0;
    jPCA_params.suppressBWrosettes=true;
    jPCA_params.suppressHistograms=true;
    jPCA_params.numPCs=ndim;
    jPCA_params.suppressText=true;
    
    for j=1:8
        data(1,j).times=1:10:300;
    end
    [Projection,Summary_data]=jPCA_v2(data,1:10:300,jPCA_params);
    R2_data_skew=Summary_data.R2_Mskew_kD;
    R2_data_best=Summary_data.R2_Mbest_kD;
    
    %% plot jPC planes
    plotParams.planes2plot=1;
    phaseSpace(Projection,Summary_data,plotParams);

    %% bootstrap neural data and perform jPCA analysis
    %jPCA parameters
    clear jPCA_params
    jPCA_params.softenNorm=5;
    jPCA_params.suppressBWrosettes=true;
    jPCA_params.suppressHistograms=true;
    jPCA_params.numPCs=ndim;
    jPCA_params.suppressText=true;

    for it=1:niters
        rand_ind=randperm(nneurons,nsig);
        for j=1:8
            temp(1,j).A    =M1_data(1,j).A(:,rand_ind);
            temp(1,j).times=1:10:300;
        end
        [~,Summary]=jPCA_v2(temp,1:10:300,jPCA_params);
        R2_bstrap_skew(1,it)=Summary.R2_Mskew_kD;
        R2_bstrap_best(1,it)=Summary.R2_Mbest_kD;
    end
    
    
    %% plot data as a bplot
    figure
    hold on 
    subplot(1,2,1)
    hold on 
    bplot(R2_bstrap_skew,1,'color','k');
    ylim([0,1])
    xlim([0,2])
    hline(R2_data_skew,'k')
    ylabel('Mskew R2 ')
    pval=sum(R2_bstrap_skew<R2_data_skew)/niters;
    temp_str=strcat('R2_data:',num2str(R2_data_skew),' pval:',num2str(pval));
    title(temp_str);
    
    subplot(1,2,2)
    hold on 
    bplot(R2_bstrap_best,1,'color','k');
    ylim([0,1])
    xlim([0,2])
    hline(R2_data_best,'k')
    ylabel('Mbest R2')
    pval=sum(R2_bstrap_best<R2_data_best)/niters;
    temp_str=strcat('R2_data:',num2str(R2_data_best),' pval:',num2str(pval));
    title(temp_str);
 
end



function [KINE_data,M1_data,EMG_data]=import_data(monkey_name, EMG_flag)
    nchar=length(monkey_name);
    %% EMG data

if EMG_flag
    cd('data_EMG')
    files=dir;
    for i=3:length(files)
        if strcmp(files(i).name(1:nchar),monkey_name)
            break;
        end
    end
    files(i).name
    EMG_data=load(files(i).name);
    EMG_data=EMG_data.Data_struct;
    cd('../')
else 
    EMG_data=0;
end

%% Neural Data
cd('data_neural/MC_data')
files=dir;
for i=3:length(files)
    if strcmp(files(i).name(1:nchar),monkey_name)
        break;
    end
end
files(i).name
M1_data=load(files(i).name);
M1_data=M1_data.Data_struct;
cd('../../')

%% kinematic data
cd('data_kinematics')
files=dir;
for i=3:length(files)
    if strcmp(files(i).name(1:nchar),monkey_name)
        break;
    end
end
files(i).name
KINE_data=load(files(i).name);
KINE_data=KINE_data.Data_struct;
cd('../')
end







