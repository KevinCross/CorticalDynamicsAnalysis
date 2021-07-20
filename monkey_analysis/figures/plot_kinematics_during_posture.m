close ALL
% task='Posture';
% offset=100;
% time=(offset*-1:(300-offset-1))*10-6;

% task='Reach';
% offset=87;
% time=(offset*-1:(250-offset-1))*10-60;

task='Tracking';
offset=76;
time=(offset*-1:(250-offset-1))*10-60;

net_names={'REC','NOREC'};
lay_names={'Input_layer','Output_layer','Muscle','Sensory_feedback','Kinematics'};
net_type=1;
colours={[0.937500000000000,0,0]	[0.750000000000000,0,0]	[0.687500000000000,0,0]	[0.250000000000000,0,0]	[0,0.250000000000000,0]	[0,0.750000000000000,0]	[0,0.812500000000000,0]	[0,0.937500000000000,0]	[0,1,0]	[0,0.875000000000000,0]	[0,0.687500000000000,0]	[0.0625000000000000,0,0]	[0.625000000000000,0,0]	[0.812500000000000,0,0]	[1,0,0]	[0.875000000000000,0,0]	[0,0.125000000000000,0]	[0.500000000000000,0,0]	[0.375000000000000,0,0]	[0.187500000000000,0,0]	[0,0,0]	[0,0.187500000000000,0]	[0,0.375000000000000,0]	[0,0.500000000000000,0]	[0,0.625000000000000,0]	[0,0.562500000000000,0]	[0,0.437500000000000,0]	[0,0.312500000000000,0]	[0,0,0]	[0.125000000000000,0,0]	[0.312500000000000,0,0]	[0.437500000000000,0,0]	[0.562500000000000,0,0]	[0,0.0625000000000000,0]};
net_colour={'k','r'};
fig_dir=strcat('figures/kinematic_analysis/',task,'/kinematic_figures');
mkdir(fig_dir)



for net=1:2
    kin_dir=strcat('data_network_update_15_07_2021/',task,'/',net_names{net},'/',lay_names{5});
    kindata=get_data(kin_dir);
    sens_dir=strcat('data_network_update_15_07_2021/',task,'/',net_names{net},'/',lay_names{4});
    sensdata=get_data(sens_dir);

    theta=linspace(0,2*pi,16);
    theta=[theta, theta];
    %% loop over individual reach conditions
    for cond=1:size(kindata,2)
        tmp_kin=squeeze(kindata(:,cond,:));
        tmp_sens=squeeze(sensdata(:,cond,:));
%         tmp_kin=joint_to_cart(tmp_kin);
    %% plot kinematic figures
        figure(1)
        hold on 
        subplot(1,2,net)
        hold on
        handx=tmp_kin(:,1)-tmp_kin(offset-10,1);
        handy=tmp_kin(:,2)-tmp_kin(offset-10,2);
        plot(handx,handy,'color',colours{cond})
        
        if strcmp('Posture',task)
            scatter(handx(offset+30,1),handy(offset+30,1),'filled','k')
        end
        title(net_names{net})
%         assert(false)
        axis square
        
        
        %% plot shoulder angle and velocity
        if strcmp(task,'Posture')
            sho_ang=tmp_sens(:,3);
            sho_vel=tmp_sens(:,5);
        else
            sho_ang=tmp_sens(:,4)+tmp_sens(:,2);
            sho_vel=tmp_sens(:,6);
        end
        figure(2)
        hold on 
        subplot(2,2,net)
        hold on 
        plot(time,sho_ang,'color',colours{cond})
        xlim([-200,800])
        vline(0,'k')
        axis square
        ylabel('Shoulder Angle (rads)')
        title(net_names{net})

        subplot(2,2,net+2)
        hold on 
        plot(time,sho_vel,'color',colours{cond})
        xlim([-200,800])
        vline(0,'k')
        axis square
        ylabel('Shoulder Angle Velocity (rads/s)')
        title(net_names{net})
        
        
       %% plot elbow angle and velocity
        if strcmp(task,'Posture')
            sho_ang=tmp_sens(:,4);
            sho_vel=tmp_sens(:,6);
        else
            sho_ang=tmp_sens(:,5)+tmp_sens(:,3);
            sho_vel=tmp_sens(:,7);
        end
        figure(3)
        hold on 
        subplot(2,2,net)
        hold on 
        plot(time,sho_ang,'color',colours{cond})
        xlim([-200,800])
        vline(0,'k')
        axis square
        ylabel('Elbow Angle (rads)')
        title(net_names{net})

        subplot(2,2,net+2)
        hold on 
        plot(time,sho_vel,'color',colours{cond})
        xlim([-200,800])
        vline(0,'k')
        axis square
        ylabel('Elbow Angle Velocity (rads/s)')
        title(net_names{net})
        
        %% plot shoulder muscle
        sho_flx=tmp_sens(:,7);
        sho_ext=tmp_sens(:,8);
        figure(4)
        hold on 
        subplot(2,2,net)
        hold on 
        plot(time,sho_flx,'color',colours{cond})
        xlim([-200,800])
        vline(0,'k')
        axis square
        ylabel('Shoulder Flexor (a.u)')
        title(net_names{net})

        subplot(2,2,net+2)
        hold on 
        plot(time,sho_ext,'color',colours{cond})
        xlim([-200,800])
        vline(0,'k')
        axis square
        ylabel('Shoulder Extensor (a.u)')
        title(net_names{net})
        
        
        %% plot elbow muscle
        sho_flx=tmp_sens(:,9);
        sho_ext=tmp_sens(:,10);
        figure(5)
        hold on 
        subplot(2,2,net)
        hold on 
        plot(time,sho_flx,'color',colours{cond})
        xlim([-200,800])
        vline(0,'k')
        axis square
        ylabel('Elbow Flexor (a.u)')
        title(net_names{net})

        subplot(2,2,net+2)
        hold on 
        plot(time,sho_ext,'color',colours{cond})
        xlim([-200,800])
        vline(0,'k')
        axis square
        ylabel('Elbow Extensor (a.u)')
        title(net_names{net})
        
        
        %% plot biarticular muscles
        sho_flx=tmp_sens(:,11);
        sho_ext=tmp_sens(:,12);
        figure(6)
        hold on 
        subplot(2,2,net)
        hold on 
        plot(time,sho_flx,'color',colours{cond})
        xlim([-200,800])
        vline(0,'k')
        axis square
        ylabel('Elbow Flexor (a.u)')
        title(net_names{net})

        subplot(2,2,net+2)
        hold on 
        plot(time,sho_ext,'color',colours{cond})
        xlim([-200,800])
        vline(0,'k')
        axis square
        ylabel('Elbow Extensor (a.u)')
        title(net_names{net})
    end

end
% assert(false)
h=figure(1);
for i=1:2
    subplot(1,2,i)
    xlim([-0.06,0.06])
    ylim([-0.06,0.06])
    axis square
end

assert(false)
saveas(h,strcat(fig_dir,'/handpaths'))
saveas(h,strcat(fig_dir,'/handpaths'),'epsc')

h=figure(2);
saveas(h,strcat(fig_dir,'/handspeed'))
saveas(h,strcat(fig_dir,'/handspeed'),'epsc')
        
h=figure(3);
saveas(h,strcat(fig_dir,'/shoulder_angle_velocity'))
saveas(h,strcat(fig_dir,'/shoulder_angle_velocity'),'epsc')

h=figure(4);
saveas(h,strcat(fig_dir,'/elbow_angle_velocity'))
saveas(h,strcat(fig_dir,'/elbow_angle_velocity'),'epsc')

h=figure(5);
saveas(h,strcat(fig_dir,'/shoulder_muscles'))
saveas(h,strcat(fig_dir,'/shoulder_muscles'),'epsc')

h=figure(6);
saveas(h,strcat(fig_dir,'/elbow_muscles'))
saveas(h,strcat(fig_dir,'/elbow_muscles'),'epsc')

h=figure(7);
saveas(h,strcat(fig_dir,'/movement_time'))
saveas(h,strcat(fig_dir,'/movement_time'),'epsc')

function state_vec=joint_to_cart(state_args)

    l1 = 0.145;
    l2 = 0.284;

    sho_ang=state_args(:,1);
    elb_ang=state_args(:,2);
    sho_vel=state_args(:,3);
    elb_vel=state_args(:,4);
    x=(l1*cos(sho_ang))+l2*cos(sho_ang+elb_ang);
    y=l1*sin(sho_ang)+l2*sin(sho_ang+elb_ang);
    
    xvel=-l1.*sho_vel.*sin(sho_ang)-l2.*(sho_vel+elb_vel).*sin(sho_ang+elb_ang);
    yvel=l1.*sho_vel.*cos(sho_ang)+l2.*(sho_vel+elb_vel).*cos(sho_ang+elb_ang);
    
    
    state_vec = [x, y,xvel,yvel];

end

function data=get_data(dir_name)

    files=dir(dir_name);
    data=[];
    for i=3:length(files)
        tmp=load(strcat(dir_name,'/',files(i).name));
        fname=fieldnames(tmp);
        if i==3
            data=tmp.(fname{1});
        else
            data=cat(2,data,tmp.(fname{1}));
        end
    end

end


function [counter,location]=calc_time_to_target(kindata,target_loc,offset,task)
%assume kindata
    %calculate distance between hand and target
    if strcmp(task,'Reach')
        dist=((kindata(:,1)-target_loc(1)).^2+(kindata(:,2)-target_loc(2)).^2).^0.5;
        counter=offset;
        while dist(counter)>0.005
            counter=counter+1;
        end
        location=kindata(counter,1:2);
    else
        %calculate distance between hand and target
        dist=((kindata(:,1)-target_loc(1)).^2+(kindata(:,2)-target_loc(2)).^2).^0.5;
        [dmax_dist,tmax_dist]=max(dist);
        max_vals=[dmax_dist,tmax_dist];
        counter=tmax_dist;
        while dist(counter)>0.005
            counter=counter+1;
        end
        return_time=counter;
    end
end
    

