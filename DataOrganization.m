%% delete sub
clear;clc
info_path='F:\Projects\Movie_graph_theory_2023_9_19\result\delete sub.xlsx';
[~,~,raw]=xlsread(info_path);
delete_name=string(raw);
sup_path='F:\Projects\Movie_graph_theory_2023_9_19\result\FIX-Denoised_200';
run_dir=dir(sup_path);
run_dir(1:2)=[];
for irun=1:4
    data_path=fullfile(sup_path,run_dir(irun).name,'GraphTheory_sparsity');
    sub_dir=dir(data_path);
    sub_dir(1:2)=[];
    for isub=1:length(sub_dir)
        if sum(strcmp(delete_name,sub_dir(isub).name(1:end-4)))
            file=fullfile(data_path,sub_dir(isub).name);
            recycle('on') % 设置'off'，则不进入回收站，直接删除。
            delete(file);
        end
    end
end

for irun=5:8
    clip_dir=dir(fullfile(sup_path,run_dir(irun).name,'clips_originvalue'));
    clip_dir(1:2)=[];
    for iclip=1:length(clip_dir)
        sub_dir=dir(fullfile(clip_dir(iclip).folder,clip_dir(iclip).name,'GraphTheory_sparsity'));
        sub_dir(1:2)=[];
        for isub=1:length(sub_dir)
            if sum(strcmp(delete_name,sub_dir(isub).name(1:end-4)))
                file=fullfile(sub_dir(isub).folder,sub_dir(isub).name);
                recycle('on') % 设置'off'，则不进入回收站，直接删除。
                delete(file);
            end
        end
    end
end


%% organize the head movement data to one file
clear;clc
% get the subjects' name which is analysysed
sub_dir=dir('F:\Projects\20230919_MovieGraphtheory\result\FIX-Denoised_200\tfMRI_MOVIE1_7T_AP\clips_originvalue\clip1\GraphTheory_sparsity');
sub_dir(1:2)=[];

run_names=["tfMRI_MOVIE1_7T_AP","tfMRI_MOVIE2_7T_PA","tfMRI_MOVIE3_7T_PA","tfMRI_MOVIE4_7T_AP",];
movedata_names=["Movement_AbsoluteRMS.txt","Movement_RelativeRMS.txt","Movement_Regressors.txt","Movement_Regressors_dt.txt"];
clips_label1=[[21,274];[285,516];[527,724];[735,808];[819,911]]+10;  % The start and end time_point of the clips of movie1
clips_label2=[[21,257];[268,536];[547,805];[816,898]]+10;    % The start and end time_point of the clips of movie2
clips_label3=[[21,211];[222,416];[427,640];[651,802];[813,905]]+10;    % The start and end time_point of the clips of movie3
clips_label4=[[21,263];[274,513];[524,788];[799,891]]+10;  
clips_num=[5,4,5,4];   % the clips number of each movie

for irun=1:4
    for isub=1:length(sub_dir)
        data_path=fullfile('F:\DataBase\HCP_MOVIE',sub_dir(isub).name(1:end-4),'MNINonLinear','Results',run_names(irun));
        % AbsoluteRMS_mean
        movedata_path=fullfile(data_path,movedata_names(1));
        data_AbsoluteRMS=load(movedata_path);
        % RelativeRMS_mean
        movedata_path=fullfile(data_path,movedata_names(2));
        data_RelativeRMS=load(movedata_path);
        for iclip=1:clips_num(irun)
            s=strcat('label1=','clips_label',num2str(irun),'(',num2str(iclip),',1);');
            eval(s);
             s=strcat('label2=','clips_label',num2str(irun),'(',num2str(iclip),',2);');
            eval(s);
        HeadMovement.AbsoluteRMS_mean(isub,irun,iclip)=mean(data_AbsoluteRMS(label1:label2));
        HeadMovement.RelativeRMS_mean(isub,irun,iclip)=mean(data_RelativeRMS(label1:label2));
        end
    end
end
RelativeRMS_mean_validation=(HeadMovement.RelativeRMS_mean(:,1,5)+HeadMovement.RelativeRMS_mean(:,1,5)...
    +HeadMovement.RelativeRMS_mean(:,1,5)+HeadMovement.RelativeRMS_mean(:,1,5))/4;   % validation clip

result_path='F:\Projects\20230919_MovieGraphtheory\result';
cd(result_path)
save HeadMovement.mat HeadMovement RelativeRMS_mean_validation
