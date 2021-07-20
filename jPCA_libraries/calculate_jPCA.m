function [Projection_Data,Summary_Data]=calculate_jPCA(Data_struct,times,jPCA_params)
    numb_it=1000; %changed from originally 1000->10
    

    %run jPCA analysis
    [Projection_Data,Summary_Data,bigA]=jPCA_v2(Data_struct,times,jPCA_params);    

    %plot jPC planes
    plotParams.planes2plot=[1,2,3];
    [color]=phaseSpace(Projection_Data,Summary_Data,plotParams);
    
    
        %% calculate null distribution for goodness of fit
    %tensor Maximum entropy
    %convert to a data tensor 
    nTimes=size(Data_struct(1).A,1);
    data_tensor=[];
    for j=1:length(Data_struct)
        data_tensor(j,:,:)=bigA(1+(j-1)*nTimes:j*nTimes,:);
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

end