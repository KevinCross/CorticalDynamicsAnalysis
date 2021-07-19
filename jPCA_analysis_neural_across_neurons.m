%this scripts pools neurons across monkeys and areas and performs jPCA
%analysis on the pooled data
close ALL

%% set up jPCA parameters
clear jPCA_params
jPCA_params.softenNorm=5;
jPCA_params.suppressBWrosettes=true;
jPCA_params.suppressHistograms=true;
jPCA_params.numPCs=6;

%% loop over files and add them to Data_struct 
%get file names
% file_names={'data_neural/S1_data/Patches_Data_S1','data_neural/A5_data/Patch_A5','data_neural/A2_data/Ayur_A2','data_neural/A5_data/Ayur_A5'};
% area='SC_all_collapsed';
file_names={'data_neural/MC_data/Patches_Data','data_neural/MC_data/Puran_Data','data_neural/MC_data/Xander_Data','data_neural/MC_data/Manny_Data_set1','data_neural/MC_data/Ayur_MC_Data'};
area='MC_all_collapsed';

Data_struct=struct('A',[]);
ncount=0;
for n=1:length(file_names)
    temp=load(file_names{n});
    nNeuronCurr=size(temp.Data_struct(1).A,2);
    for cond=1:length(temp.Data_struct)
        Data_struct(cond).A(:,ncount+1:ncount+nNeuronCurr)=temp.Data_struct(cond).A;
        Data_struct(cond).times=temp.Data_struct(cond).times;
    end
    ncount=ncount+nNeuronCurr;
end

%perform jPCA analysis
[Projection_Data,Summary_Data]=calculate_jPCA(Data_struct,Data_struct(1).times,jPCA_params);

%% save figures
% assert(false)
mkdir(strcat('figures/',area,'/'))
cd(strcat('figures/',area,'/'))
for j=1:3
    h=figure(j);
    file_name=strcat('jPCA_plot_',num2str(j));
    saveas(h,file_name)
    saveas(h,file_name,'epsc')
end
h=figure(4);
saveas(h,'Metrics vs TME Null Distribution')
saveas(h,'Metrics vs TME Null Distribution','epsc')
