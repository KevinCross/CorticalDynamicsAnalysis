% this script performs jPCA analysis for each monkey/area separately and
% saves the resulting figures

close ALL
% area can be one of: MC_data, S1_data, A5_data, A2_data
area='MC_data';
curr_dir=strcat('data_neural/',area);
files=dir(curr_dir);

%% set up jPCA parameters
clear jPCA_params
jPCA_params.softenNorm=5;
jPCA_params.suppressBWrosettes=true;
jPCA_params.suppressHistograms=true;
jPCA_params.numPCs=6;

%% loop over files and perform jPCA analysis
for i=3:size(files,1)
    load(strcat(curr_dir,'/',files(i).name));
    [Proj_Neural,Summary_Neural]=calculate_jPCA(Data_struct,Data_struct(1).times,jPCA_params);
    
    %save files
    mkdir(strcat('figures/',area))
    cd(strcat('figures/',area))
    mkdir(files(i).name)
    cd(files(i).name)
    for j=1:3
        h=figure(j);
        file_name=strcat(files(i).name(1:6),'_jPCA_plot_',num2str(j));
        saveas(h,file_name)
        saveas(h,file_name,'epsc')
    end
    h=figure(4);
    saveas(h,'Metrics vs TME Null Distribution')
    saveas(h,'Metrics vs TME Null Distribution','epsc')
    save('Projection_Neural','Proj_Neural')
    save('Summary_Neural','Summary_Neural')
    close ALL
    cd('../../../')
end
    