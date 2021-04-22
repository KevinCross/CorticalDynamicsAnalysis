close ALL

%get file names
area='MC_data';
fold=strcat('../data/',area);
data_sets=dir(fold);
numb_it=1000;

%% set up jPCA parameters
clear jPCA_params
jPCA_params.softenNorm=5;
jPCA_params.suppressBWrosettes=true;
jPCA_params.suppressHistograms=true;
jPCA_params.numPCs=6;

%% arrays to hold the relevant metrics

%% loop over files and perform jPCA analysis
for i=4:4%size(files,1)
    %% load data
    data_struct=load_data(data_sets(i).name,fold);
    times=data_struct(1).times;

    
    %% run jPCA analysis and plot
    [Proj_Neural,Summary_Neural,bigA]=jPCA_v2(data_struct,times,jPCA_params);    
    
    %plot jPC planes
    plotParams.planes2plot=[1,2,3];
    phaseSpace(Proj_Neural,Summary_Neural,plotParams);
    
    
    %% calculate null distribution for goodness of fit
    [null_R_mskew,null_R_mbest]=calculate_null_distribution(data_struct,bigA,numb_it,times,jPCA_params);
        
    %% plot data
    %plot null distribution as a box plot
    figure
    hold on 
    subplot(1,2,1)
    hold on 
    bplot(null_R_mbest,1,'color','k');
    hline(Summary_Neural.R2_Mbest_kD,'k')
    title('R_MBest_all')
    text(0,0,strcat('pval: ',num2str(sum(null_R_mbest>Summary_Neural.R2_Mbest_kD)/numb_it)))
    
    subplot(1,2,2)
    hold on 
    bplot(null_R_mskew,1,'color','k');
    hline(Summary_Neural.R2_Mskew_kD,'k')
    title('R_MSkew_all')
    text(0,0,strcat('pval: ',num2str(sum(null_R_mskew>Summary_Neural.R2_Mskew_kD)/numb_it)))
    
    for f=1:2
        subplot(1,2,f)
        ylim([0,1])
        xlim([0,2])
    end
        
    %% save figures
    fold_name=strcat('figures_test/',area,'/',files(i).name);
    mkdir(fold_name)
    for j=1:3
        h=figure(j);
        file_name=strcat(fold_name,'/',data_sets(i).name(1:6),'_jPCA_plot_',num2str(j));
        saveas(h,file_name)
        saveas(h,file_name,'epsc')
    end
    h=figure(4);
    saveas(h,strcat(fold_name,'/','Metrics vs TME Null Distribution'))
    saveas(h,strcat(fold_name,'/','Metrics vs TME Null Distribution'),'epsc')
    save(strcat(fold_name,'/','Projection_Neural'),'Proj_Neural')
    save(strcat(fold_name,'/','Summary_Neural'),'Summary_Neural')
    close ALL
    
end




function data_struct=load_data(fil_name,fold)
    %load data
    temp=load(strcat(fold,'/',fil_name));
    fil_name
    %check to see if times are present in Data_pert struct
    if isfield(temp.Data_pert,'times')
        times=Data_pert(1).times(20:49);
    else
        times=1:10:300;
        for j=1:8
            temp.Data_pert(1,j).times=times;
        end
    end
    data_struct=temp.Data_pert;
end

function [null_R_mskew,null_R_mbest]=calculate_null_distribution(data_struct,bigA,numb_it,times,jPCA_params)

%% calculate null distribution for goodness of fit
    %tensor Maximum entropy
    %convert to a data tensor 
    data_tensor=[];
    for j=1:length(data_struct)
        data_tensor(j,:,:)=bigA(1+(j-1)*30:j*30,:);
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
    null_R_mskew=zeros(numb_it,1);
    null_R_mbest=zeros(numb_it,1);
    for it=1:numb_it
       [surrTensor1] = sampleTME(maxEntropy);  
       null_Data=struct;
       for k=1:length(data_struct)
           null_Data(k).A=squeeze(surrTensor1(k,:,:));
           null_Data(k).times=times;
       end
       [~,null_Summary]=jPCA_v2(null_Data,times,jPCA_params);
       null_R_mskew(it)=null_Summary.R2_Mskew_kD;
       null_R_mbest(it) =null_Summary.R2_Mbest_kD;
    end
end