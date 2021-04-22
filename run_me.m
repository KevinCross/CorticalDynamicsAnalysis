close ALL
task='Posture';
offset=110; %time offset for when to start calculating jPCA
% task='Reach';
% offset=87;
net_names={'REC','NOREC'};
lay_names={'Input_layer','Output_layer','Muscle','Sensory_feedback','Kinematics'};

%% analysis of neural activities
%loop over the network configurations and layer types (only Input_layer and
%Output_layer)
for net_type=3:2
    for lay_type=1:2
        %get directory that the data are stored in
        curr_dir=strcat('data/',task,'/',net_names{net_type},'/',lay_names{lay_type});
        
        %run jPCA analysis on data
        [Projection_Data,Summary_Data]=calculate_jPCA(curr_dir,offset);
        
        %create directory for saving figures
        dir_name=strcat('figures/',task,'/',net_names{net_type},'/',lay_names{lay_type});
        mkdir(dir_name)
        for j=1:3
            h=figure(j);
            file_name=strcat(dir_name,'/jPCA_plot_',num2str(j));
            saveas(h,file_name)
            saveas(h,file_name,'epsc')
        end
        h=figure(4);
        saveas(h,strcat(dir_name,'/Metrics vs TME Null Distribution'))
        saveas(h,strcat(dir_name,'/Metrics vs TME Null Distribution'),'epsc')
        save(strcat(dir_name,'/Projection_Struct'),'Projection_Data')
        save(strcat(dir_name,'/Summary_Struct'),'Summary_Data')
        close ALL
    end
end

%% analysis of muscle and kinematic signals 
for net_type=1:2
    %get directors that store the data
    Neural_dir=strcat('data/',task,'/',net_names{net_type},'/',lay_names{2});
    
    Sens_dir=strcat('data/',task,'/',net_names{net_type},'/',lay_names{4});%note the kinematic signals are pulled from the sensory feedback
    %run jPCA analysis on the data
    bstrap_kine_and_muscle_signals(Neural_dir,Sens_dir,offset,task);
    %create directory to save figures to
    assert(false)
    %muscle figures
    dir_name=strcat('figures/',task,'/',net_names{net_type},'/',lay_names{3});
    mkdir(dir_name)
    h=figure(1);
    file_name=strcat(dir_name,'/jPCA_plot_1');
    saveas(h,file_name)
    saveas(h,file_name,'epsc')
    h=figure(2);
    file_name=strcat(dir_name,'/R2_vs_null_distribution');
    saveas(h,file_name)
    saveas(h,file_name,'epsc')
    
    
%     dir_name=strcat('figures/',task,'/',net_names{net_type},'/',lay_names{5});
%     mkdir(dir_name)
%     h=figure(3);
%     file_name=strcat(dir_name,'/jPCA_plot_1');
%     saveas(h,file_name)
%     saveas(h,file_name,'epsc')
%     h=figure(4);
%     file_name=strcat(dir_name,'/R2_vs_null_distribution');
%     saveas(h,file_name)
%     saveas(h,file_name,'epsc')
%     
    %kinematic signals
    dir_name=strcat('figures/',task,'/',net_names{net_type},'/',lay_names{4});
    mkdir(dir_name)
    save(strcat(dir_name,'/Projection_Struct'),'Proj_Sens_All')
    save(strcat(dir_name,'/Summary_Struct'),'Summary_Sens_All')

    close ALL
end



function []=bstrap_kine_and_muscle_signals(Neural_dir,Sens_dir,offset,task)
    niters=1000;
    [Data_Neural]=get_data_from_directory(Neural_dir,offset);
    [Data_Sens]=get_data_from_directory(Sens_dir,offset);
        
    %% extract the kinematic signals from the sensory feedback
    if strcmp(task,'Posture')
        Data_Kin=struct;
        for c=1:length(Data_Sens)
            Data_Kin(c).A=Data_Sens(c).A(:,3:6);
            Data_Kin(c).times=Data_Sens(c).times;
        end
        
    else
        %for reaching
        Data_Kin=struct;
        for c=1:length(Data_Sens)
            Data_Kin(c).A(:,1)=Data_Sens(c).A(:,4)+Data_Sens(c).A(:,2);
            Data_Kin(c).A(:,2)=Data_Sens(c).A(:,5)+Data_Sens(c).A(:,3);
            Data_Kin(c).A(:,3)=Data_Sens(c).A(:,6);
            Data_Kin(c).A(:,4)=Data_Sens(c).A(:,7);
            Data_Kin(c).times=Data_Sens(c).times;
        end
    end
    
    %% extract the muscle activities from the sensory feedback
    Data_Mus=struct;
    for c=1:length(Data_Sens)
        Data_Mus(c).A=Data_Sens(c).A(:,end-6:end);
        Data_Mus(c).times=Data_Sens(c).times;
    end
        
    %perform jPCA and bootstrap analysis
    [Proj_Musc] = bstrap_neural_signals(Data_Mus,Data_Neural,4,niters,[1]); 
    [Proj_Kine] = bstrap_neural_signals(Data_Kin,Data_Neural,2,niters,[1]);
    
end

function state_vec=convert_cart_to_jt(state_args)

    a1 = 0.145;
    a2 = 0.284;

    x = squeeze(state_args(1));
    y = squeeze(state_args(2));
    xvel = squeeze(state_args(3));
    yvel = squeeze(state_args(4));


    t2 = acos((x.^2 + y.^2 - a1.^2 - a2.^2)/(2*a1*a2));
    t1 = atan2(y, x) - atan2(a2*sin(t2), (a1 + a2*cos(t2)));

    J = [-a1*sin(t1)-a2*sin(t1+t2), -a2*sin(t1+t2); a1*cos(t1)+a2*cos(t1+t2), a2*cos(t1+t2)];
    jvel = inv(J) * [xvel; yvel];

    state_vec = [t1, t2, jvel.'];

end
        
        
        
        



function [Projection,Summary_data] = bstrap_neural_signals(data,M1_data,ndim,niters,port_seq)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    %assume data and M1_data are a struct with fieldname A
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
    jPCA_params.softenNorm=0.001;
    jPCA_params.suppressBWrosettes=true;
    jPCA_params.suppressHistograms=true;
    jPCA_params.numPCs=ndim;
    
%     for j=1:8
%         data(1,j).times=1:10:300;
%     end
    [Projection,Summary_data]=jPCA_v2(data,0:10:300,jPCA_params);
    R2_data_skew=Summary_data.R2_Mskew_kD;
    R2_data_best=Summary_data.R2_Mbest_kD;
    
    %% plot jPC planes
    plotParams.planes2plot=[port_seq];
    phaseSpace(Projection,Summary_data,plotParams);

    %% bootstrap neural data and perform jPCA analysis
    jPCA_params.suppressText=true;
    for it=1:niters
        rand_ind=randperm(nneurons,nsig);
        for j=1:8
            temp(1,j).A    =M1_data(1,j).A(:,rand_ind);
            temp(1,j).times=0:10:300;
        end
        [~,Summary]=jPCA_v2(temp,0:10:300,jPCA_params);
        R2_bstrap_skew(1,it)=Summary.R2_Mskew_kD;
        R2_bstrap_best(1,it)=Summary.R2_Mbest_kD;
    end
    
    
    %% plot data as a bplot
    figure
    hold on 
    subplot(1,2,1)
    hold on 
    bplot(R2_bstrap_skew,1,'color','k');
    hline(R2_data_skew,'k')
    ylim([0,1])
    xlim([0,2])
    ylabel('Mskew R2 ')
    pval=sum(R2_bstrap_skew<R2_data_skew)/niters;
    temp_str=strcat('R2_data:',num2str(R2_data_skew),' pval:',num2str(pval));
    title(temp_str);
    
    subplot(1,2,2)
    hold on 
    bplot(R2_bstrap_best,1,'color','k');
    hline(R2_data_best,'k')
    ylim([0,1])
    xlim([0,2])
    ylabel('Mbest R2')
    pval=sum(R2_bstrap_best<R2_data_best)/niters;
    temp_str=strcat('R2_data:',num2str(R2_data_best),' pval:',num2str(pval));
    title(temp_str);
 
end







function [Projection_Data,Summary_Data]=calculate_jPCA(task,offset)
    
    %% import data
    curr_dir=pwd;
    cd(task)
    files=dir;
    times=0:10:300;
    numb_it=1000; %changed from originally 1000->10
    clear Data_struct
    Data_struct=struct;
    count=1;
    for i=3:size(files,1)
        %load data
        tmp=load(files(i).name);
        fname=fieldnames(tmp);
        files(i).name
        %check to see if times are present in Data_pert struct
        for tp=1:size(tmp.(fname{1}),2)
            Data_struct(count).A=squeeze(tmp.(fname{1})(offset:offset+30,tp,:));
            Data_struct(count).times=times;
            count=count+1;
        end
    end
    %% jPCA parameters
    clear jPCA_params
    jPCA_params.softenNorm=0.0000;
    jPCA_params.suppressBWrosettes=true;
    jPCA_params.suppressHistograms=true;
    jPCA_params.numPCs=6;

    %run jPCA analysis
    [Projection_Data,Summary_Data,bigA]=jPCA_v2(Data_struct,times,jPCA_params);    

    %plot jPC planes
    plotParams.planes2plot=[1,2,3];
    [color]=phaseSpace(Projection_Data,Summary_Data,plotParams);
    
    
    %% calculate null distribution for goodness of fit
    %tensor Maximum entropy
    %convert to a data tensor 
    data_tensor=[];
    for j=1:length(Data_struct)
        data_tensor(j,:,:)=bigA(1+(j-1)*31:j*31,:);
    end

    %find maximum entropy of tensor
    [targetSigmaT, targetSigmaN, targetSigmaC, M] = extractFeatures2(data_tensor);
    params = [];
    params.margCov{1} = targetSigmaT;
    params.margCov{2} = targetSigmaN;
    params.margCov{3} = targetSigmaC;
    params.meanTensor = M.TNC; 
    maxEntropy = fitMaxEntropy(params);

    
    %sample from tensor entropy
    
    jPCA_params.suppressText=true;
    null_R_mskew=zeros(numb_it,1);
    null_R_mbest=zeros(numb_it,1);
    for it=1:numb_it
       [surrTensor1] = sampleTME(maxEntropy);  
       null_Data=struct;
       for k=1:length(Data_struct)
           null_Data(k).A=squeeze(surrTensor1(k,:,:));
           null_Data(k).times=times;
       end
       [Projection,null_Summary]=jPCA_v2(null_Data,times,jPCA_params);
       null_R_mskew(it)=null_Summary.R2_Mskew_kD;
       null_R_mbest(it) =null_Summary.R2_Mbest_kD;
    end
    
    %% plot data
    %plot null distribution as a box plot
    figure
    hold on 
    subplot(1,2,1)
    hold on 
    bplot(null_R_mbest,1,'color','k');
    hline(Summary_Data.R2_Mbest_kD,'k')
    ylim([0,1])
    xlim([0,2])
    title('R_MBest_all')
    text(0,0,strcat('pval: ',num2str(sum(null_R_mbest>Summary_Data.R2_Mbest_kD)/numb_it)))
    
    subplot(1,2,2)
    hold on 
    bplot(null_R_mskew,1,'color','k');
    hline(Summary_Data.R2_Mskew_kD,'k')
    ylim([0,1])
    xlim([0,2])
    title('R_MSkew_all')
    text(0,0,strcat('pval: ',num2str(sum(null_R_mskew>Summary_Data.R2_Mskew_kD)/numb_it)))

    
    cd(curr_dir)
end
