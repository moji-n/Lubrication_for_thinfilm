% Copyright: Mojtaba Norouzisadeh

%  e-mail: norouzi.mojtaba.sade@gmail.com
%
% date: 08/2023
% 
%
%
%This file is part of Lubrication_for_thinfilm.
%
%    Lubrication_for_thinfilm is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    Lubrication_for_thinfilm is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Diffmicro.  If not, see <https://www.gnu.org/licenses/>.
%

clear 
clc
clf
    
case_dir=pwd;


%% Open and read inputs


parameters=read_input(case_dir,'parameters');

time_step=read_input(case_dir,'timesteps');
tt=size(time_step,1)+1;

dis_parameters=read_input(case_dir,'dis_parameters');

    
%% Initialization of parameters
cluster=parameters(1); n_mobility=parameters(2); b=parameters(3);
plot_graph=parameters(4); laplacian_loop_mx=parameters(5); roughness_multiplyer=parameters(6);
alpha=parameters(7); del_t=parameters(8); max_del_t=parameters(9); end_time=parameters(10);
gravity=parameters(11); surf_tension=parameters(12); density=parameters(13);
xc=parameters(14); h_c=parameters(15); viscosity=parameters(16);
ex_f=parameters(17); x_dim=parameters(18); num_threads=parameters(19); start_from_zero=parameters(20);
del_x_init=parameters(21); init_profile=parameters(22); theta=parameters(23);potential_var=parameters(24:28);
if_refinement=parameters(29);period_bc=parameters(30);mobility_slop=parameters(31);slip_length=parameters(32);del_x_min=parameters(33);

g_coef=xc^2/(surf_tension/(density*gravity));
time_c=(3.0*viscosity*xc^4)/(surf_tension*h_c^3);
eps_ed=b^(10.0/3.0);

ham_cte1=dis_parameters(1); ham_cte2=dis_parameters(2); born_length=dis_parameters(3);
kb=dis_parameters(4); T=dis_parameters(5); epsilon=dis_parameters(6)*dis_parameters(7);
elec_charg=dis_parameters(8); avagad_num=dis_parameters(9);
nb=dis_parameters(10); zeta1=dis_parameters(11); zeta2=dis_parameters(12);
Lambda_2=dis_parameters(13); w0_h=dis_parameters(14);

psi_1=elec_charg*zeta1/(kb*T);
psi_2=elec_charg*zeta2/(kb*T);
psi_ave=(psi_1+psi_2)*0.5;
c=psi_ave^2;
d=((psi_1-psi_2))^2;
debye=sqrt(2*elec_charg^2*nb*avagad_num/(epsilon*kb*T));


dis_coef=xc^2/(surf_tension);

if cluster==1
    LASTN = maxNumCompThreads(num_threads); %on cluster
    fprintf('change max number of threads from %d into %d \n',LASTN,maxNumCompThreads)
end

%% Start from zero or not
qq=2;

if start_from_zero==1
    fprintf('solution from time zero, deleting previous results \n')
    delete(strcat(case_dir,'/output/*'))
        
        if mobility_slop==0
            prcurser=b;
        elseif mobility_slop==1
            prcurser=b;
        end
    
    
    
        rc=1;
        r_end=1.1;
        x=linspace(0,x_dim,2*floor(x_dim/del_x_init));

        % generate first grid
        [grid_neigh, x]=grid_gen_refin(del_x_init,x_dim,rc,r_end,x,del_x_min,ex_f,if_refinement);
        num_p=length(x);

        del_x=zeros(1,num_p);
        del_mid_x=zeros(1,num_p);
        for i=1:num_p-1
            del_x(i)=x(i+1)-x(i);
        end
        i=num_p;
        del_x(i)= (x(i)-x(i-1));
        for i=2:num_p-1
            del_mid_x(i)=(x(i)+x(i+1))/2- (x(i)+x(i-1))/2;
        end
        i=1;
        del_mid_x(i)=(x(i)+x(i+1));
        i=num_p;
        del_mid_x(i)= (x(i)-x(i-1));


        current_time=0;
        
        h_t=zeros(1,num_p);
        h_n=zeros(1,num_p);
        h_temp_1=zeros(1,num_p);
        mob_inter=zeros(num_p-1,6);
        mob_inter_temp=zeros(num_p-1,6);
        err=zeros(1,num_p)';
        del_t_old=del_t;
        xr_1=1;
        xr_2=1;
        xr_old_1=1;
        xr_old_2=xr_1;


        if init_profile==1          %circle for gracity and curvature
            for i=1:num_p
                %x(i)=del_x*(i-1);
                if (abs(x(i))<1.0)
                    h_t(i)=(1-x(i)^2)^0.5+prcurser;
                    %h_i(i)=(1-x(i)^2)^2+b;
                else
                    h_t(i)=prcurser;
                end
                
            end
        elseif init_profile==2      % for curvature
            for i=1:num_p
                %x(i)=del_x*(i-1);
                if (abs(x(i))<1.0)
                    %h_t(i)=(1-x(i)^2)^0.5+b;
                    h_t(i)=(1-x(i)^2)^2+prcurser;
                else
                    h_t(i)=prcurser;
                end
                
            end
        elseif init_profile==3      % for VDW
            h_dim_ana=1;%ham_cte1*(h_c/xc)/(6*pi()*density*(viscosity/density)^2*h_c)/2;
            vdw_coef=h_dim_ana;
            bc_vdw_an=3;
            h_c=1;
            xc=1;
            indice_r=1;
            r_indice_end=1;
            r_indice_init=1;
            clear x h_t h_n h_temp_1 mob_inter mob_inter_temp err grid_neigh;
            x=zeros(1,del_x_init);
            % generate first grid
            [grid_neigh, x]=grid_gen_refin(del_x_init,5*2*pi()*sqrt(2),rc,r_end,x,b,ex_f,if_refinement);
            num_p=length(x);
    
            del_x=zeros(1,num_p);
            del_mid_x=zeros(1,num_p);
            for i=1:num_p-1
                del_x(i)=x(i+1)-x(i);
            end
            i=num_p;
            del_x(i)= (x(i)-x(i-1));
            for i=2:num_p-1
                del_mid_x(i)=(x(i)+x(i+1))/2- (x(i)+x(i-1))/2;
            end
            i=1;
            del_mid_x(i)=(x(i)+x(i+1));
            i=num_p;
            del_mid_x(i)= (x(i)-x(i-1));
            
            h_t=zeros(1,num_p);
            h_n=zeros(1,num_p);
            h_temp_1=zeros(1,num_p);
            mob_inter=zeros(num_p-1,6);
            mob_inter_temp=zeros(num_p-1,6);
            err=zeros(1,num_p)';
    
    
    
                for i=1:num_p
                    
                    h_t(i)=1+0.1*sin((1/sqrt(2))*x(i));
    
                end
        elseif init_profile==4 % arc for small initilization
            theta_i=80;
            xr=1;
            rad_curv=xr/sin(theta_i/180*pi());
            h_out=rad_curv*(cos(theta_i/180*pi()));

            for i=1:num_p
                %x(i)=del_x*(i-1);
                if (abs(x(i))<1.0)
                    h_t(i)=(rad_curv^2-x(i)^2)^0.5-h_out+prcurser;
                    %h_i(i)=(1-x(i)^2)^2+b;
                else
                    h_t(i)=prcurser;
                end
                
            end

        end
    

        h_temp=h_t;
        h_temp_1=h_t;
        % Roughness generator
        h_orr=zeros(1,num_p);%roughness_gen(num_p,alpha,roughness_multiplyer);
        h_rr=h_orr;
        for i=1:length(x)
            if x(i)>1
                r_indice_init=i-5;
                r_indice_end=i+5;
                break
            end
        end

        x_r_indice_end_old=x(r_indice_end);
        x_r_indice_init_old=x(r_indice_init);

    else
    
        %% reading last time step and one before it
        fprintf('reading previous result to continue from last run \n')
        output_folder = dir('output'); 
        [output_pre_f,~]=size(output_folder);
        temp_name={output_folder.name};

        temp_name=natsort(temp_name,'\d+\.?\d*([eE][-+]?\d+)?');
    
        file_address_sum=temp_name{output_pre_f-2};
        
        output_1=read_input(strcat(case_dir,'/output/',file_address_sum),'output');
        
        if mobility_slop==0
            prcurser=b;
        elseif mobility_slop==1
            prcurser=b;
        end

        del_t=output_1{1};
        xr_1=output_1{2};
        xr_2=output_1{8};
            
        xr_old_1=xr_1;
        xr_old_2=xr_2;
        rc=xr_1;

        h_t=output_1{3};
        h_temp=h_t;
        x=output_1{4};
        current_time=output_1{5};
        num_p=length(h_t);

        h_orr=output_1{10};
        h_rr=h_orr;
        h_n=zeros(num_p);

        del_x=zeros(1,num_p);
        del_mid_x=zeros(1,num_p);
        for i=1:num_p-1
            del_x(i)=x(i+1)-x(i);
        end
        i=num_p;
        del_x(i)= (x(i)-x(i-1));
        for i=2:num_p-1
            del_mid_x(i)=(x(i)+x(i+1))/2- (x(i)+x(i-1))/2;
        end
        i=1;
        del_mid_x(i)=(x(i)+x(i+1));
        i=num_p;
        del_mid_x(i)= (x(i)-x(i-1));

        grid_neigh=[output_1{9}{1}, output_1{9}{2} ,output_1{9}{3} ,output_1{9}{4}];
        for i=1:length(x)
            if x(i)>xr_1
                r_indice_init=i-5;
                r_indice_end=i+5;
                break
            end
        end

        r_indice_end=1;
        r_indice_init=1;
        indice_r=1;
        x_r_indice_end_old=x(r_indice_end);
        x_r_indice_init_old=x(r_indice_init);
        

        file_address_sum=temp_name{output_pre_f-3};
        output_1=read_input(strcat(case_dir,'/output/',file_address_sum),'output');

        del_t_old=output_1{1};
        h_temp_1=output_1{3};
        x_old_timestep=output_1{4};

        h_temp_1=interp1(x_old_timestep,h_temp_1,x);

        
end


if plot_graph==1 %&& abs((current_time)/time_step(qq) -1.0) <1e-2
    fig1=set(gcf,'Position',[100 100 1240 480]);
    subplot(2,2,1);
    plot(x,h_t,'.-','DisplayName', strcat('h at t= 0.0 '));
    legend;
    xlim([0 5])
    ylim([0 1.2]);
%pause()
end


current_time=current_time+del_t;

timer_counter=0;
nn=0;
nn_old=0;
    
stop_negativity=0;



tot_vol=0.0;
    
for i=1:num_p-1
    tot_vol=tot_vol+abs(x(i)-x(i+1))*(h_t(i)+h_t(i+1))/2;
end
% finding the correct output entry
i=1;
while current_time>time_step(i)
    i=i+1;
    qq=qq+1;
end
del_temp=0;
h_t=h_t';
l=1;
while current_time<end_time && qq < tt && del_t>1e-30

    
    fprintf("time step size = %e, current time= %e , spreading radius= %e \n",del_t,current_time, xr_1)

    [min_temp_h, ind_min]=min(h_t);
    h_min(l)=min_temp_h;
    time_index(l)=current_time;
    del_x_min=max(min_temp_h,1e-4);
    
    if  (  abs(x_r_indice_init_old-x(r_indice_init))>5*ex_f*del_x_min)&& if_refinement==1

        fprintf("refining grid\n")
        %pause();
        x_old=x;
        h_orr=h_rr;
        clear x;
        clear h_rr;
        
        x_r_indice_end_old=x_old(r_indice_end);
        x_r_indice_init_old=x_old(r_indice_init2);

        
        r_start=x_old(r_indice_init);
        r_end=x_old(r_indice_init2);

        xr_old_1=xr_1;
        xr_old_2=xr_2;

        num_p_old=num_p;
        
        x=linspace(0,x_dim,floor(x_dim/del_x_init));
        [grid_neigh, x]=grid_gen_refin(del_x_init,x_dim,r_start,r_end,x,del_x_min,ex_f,if_refinement);
        
        num_p=length(x);

        clear del_x del_mid_x
        del_x=zeros(1,num_p);
        del_mid_x=zeros(1,num_p);
        for i=1:num_p-1
            del_x(i)=x(i+1)-x(i);
        end
        i=num_p;
        del_x(i)= (x(i)-x(i-1));


        for i=2:num_p-1
            del_mid_x(i)=(x(i)+x(i+1))/2- (x(i)+x(i-1))/2;
        end
        i=1;
        del_mid_x(i)=(x(i)+x(i+1));
        i=num_p;
        del_mid_x(i)= (x(i)-x(i-1));


      


        if interpolate_refine(h_t,num_p,num_p_old,x,x_old) ~= -1
            h_t=interpolate_refine(h_t,num_p,num_p_old,x,x_old);
        else
            break
        end
        if interpolate_refine(h_temp_1,num_p,num_p_old,x,x_old) ~= -1 
            h_temp_1=interpolate_refine(h_temp_1,num_p,num_p_old,x,x_old);
        else
            break
        end
        if interpolate_refine(h_orr,num_p,num_p_old,x,x_old) ~= -1
            h_rr=interpolate_refine(h_orr,num_p,num_p_old,x,x_old);
        else
            break
        end

        
        
        
        clear  h_n   err mob_inter mob_inter_temp
        clear Coef_matrix const_matrix;
        
        h_n=zeros(1,num_p);
        mob_inter=zeros(num_p-1,6);
        mob_inter_temp=zeros(num_p-1,6);
        err=zeros(1,num_p)';

        h_temp=h_t;
        i=1;
        while x(i)<x_r_indice_init_old && i< num_p
           i=i+1;
        end
        r_indice_init=i;

        while i< num_p && x(i)<x_r_indice_end_old 
           i=i+1;
        end
        r_indice_end=i;



    end
    
        
    
    timer_counter=timer_counter+1;
    
    x_t=(1.0+120.0*current_time)^(0.2);
    if init_profile==2
        for i=1:num_p
            
            if (abs(x(i))<x_t)
                
                h_n(i)=(1.0+120.0*current_time)^(-0.2)*(1-(x(i)/x_t)^2)^2+b;
            else
                h_n(i)=b;
            end
    
        end
    elseif init_profile==3

    else
        h_n(:)=b;
    end

    
    
    
    error_laplacian_phase=1.0;
    laplacian_loop_i=0;
    k=0;
    
    stop_negativity_inside=0;
    while error_laplacian_phase>1e-8 && laplacian_loop_i<laplacian_loop_mx && stop_negativity_inside~= 1 
        %error_laplacian_phase
        laplacian_loop_i=laplacian_loop_i+1;
        Coef_matrix=zeros(num_p,num_p);
        const_matrix=zeros(num_p,1);
        k=k+1;
        h_t_2=h_t;
        

        h_grad=gradient(h_t,x(1:num_p));
        h_grad2=gradient(h_grad,x(1:num_p));
        h_grad3=gradient(h_grad2,x(1:num_p));
        h_grad4=gradient(h_grad3,x(1:num_p));
        h_grad5=gradient(h_grad4,x(1:num_p));


        h_grad_t=gradient(h_temp,x(1:num_p));

        h_grad2_t=gradient(h_grad_t,x(1:num_p));

        clear phi_der phi_der_temp


        phi_der=zeros(num_p-1,3);
        phi_der_temp=zeros(num_p-1,3);



        for i_der=1:length(phi_der)
            phi_der(i_der,:)= Disjoin_der(0.5*h_t(i_der)+0.5*h_t(i_der+1),h_c,xc,ham_cte1,w0_h,nb,zeta1,zeta2,kb,T,epsilon,elec_charg,avagad_num,Lambda_2,potential_var,0.5*h_grad(i_der)+0.5*h_grad(i_der+1) ,0.5*h_grad2(i_der)+0.5*h_grad2(i_der+1) );
            phi_der_temp(i_der,:)= Disjoin_der(0.5*h_temp(i_der)+0.5*h_temp(i_der+1),h_c,xc,ham_cte1,w0_h,nb,zeta1,zeta2,kb,T,epsilon,elec_charg,avagad_num,Lambda_2,potential_var,0.5*h_grad_t(i_der)+0.5*h_grad_t(i_der+1) ,0.5*h_grad2_t(i_der)+0.5*h_grad2_t(i_der+1) );
        end



        sum_phi_der=zeros(num_p,1);
        for i_der=1:length(phi_der)
            sum_phi_der(i_der)=phi_der(i_der,1)+phi_der(i_der,2)+phi_der(i_der,3);
        end



        % interpolate mobility on the face
        for i=1:num_p-1
%             if h_t(i)==0
%                 i
%             end
            mob_inter(i,:)=interpolate_Mobility(h_t(i),h_t(i+1),n_mobility,prcurser,potential_var,phi_der(i,:),mobility_slop,(dis_coef*ham_cte1)/(8*pi()*h_c^2*xc^2),h_grad(i),h_grad(i+1),slip_length,born_length,h_c,xc);
            mob_inter_temp(i,:)=interpolate_Mobility(h_temp(i),h_temp(i+1),n_mobility,prcurser,potential_var,phi_der_temp(i,:),mobility_slop,(dis_coef*ham_cte1)/(8*pi()*h_c^2*xc^2),h_grad_t(i),h_grad_t(i+1),slip_length,born_length,h_c,xc);
        end

        
        [Coef_matrix, const_matrix]=populate_matrixes(num_p,grid_neigh,h_temp,x,del_x,del_mid_x,theta,mob_inter,mob_inter_temp,del_t,g_coef,dis_coef,Coef_matrix, const_matrix,period_bc);
            

        h_t=Coef_matrix\const_matrix;
        cond_num=cond(Coef_matrix);

        if cond_num > 5e4
            stop_negativity_inside=1;
        end
        

        error_laplacian_phase=0.0;
        for i=1:num_p
            if h_t_2(i)~=0
                error_laplacian_phase=error_laplacian_phase+ ((abs(h_t_2(i)-h_t(i))))/h_t_2(i);
            end
        end

        for i=1:num_p
            if h_t(i)<0 
                
                stop_negativity_inside=1;
                break
            end
        end
    %break    
    end
    
    error_vol(l)=0.0;
    
%laplacian_loop_i,error_laplacian_phase
    
    for i=1:num_p-1
        error_vol(l)=error_vol(l)+abs(x(i)-x(i+1))*(h_t(i)+h_t(i+1))/2;
        if h_temp(i)>1e-9 && h_t(i)>1e-9
            err(i)=(2*del_t/del_t_old)*(del_t_old*h_t(i)+del_t*h_temp_1(i)-(del_t_old+del_t)*h_temp(i))/((del_t_old+del_t)*h_temp(i));
        end
    end
    
    error_vol(l)=((error_vol(l)-tot_vol))/tot_vol*100;
    


    for i=1:num_p
        if h_t(i)<0
            if mobility_slop==0
                stop_negativity_inside=1;
            end
            fprintf("negative film thickness at index= %d, x= %e \n ",i,x(i))

        end

        if isnan(h_t(i))
            stop_negativity_inside=1;
            fprintf("NAN film thickness at index= %d, x= %e \n ",i,x(i))

            
        end


    end
    
    %max(err)
    if max(err)<1e-3  && stop_negativity_inside==0
        
        if plot_graph==1 %&& abs((current_time)/time_step(qq) -1.0) <1e-2
        subplot(2,3,1);
        plot(x(1:1:num_p),h_t(1:1:num_p),'.-','DisplayName', strcat('h at t= ',num2str(current_time)));
        hold on;
        plot(x(1:num_p),h_n,'r','DisplayName', strcat('h_{analytical} at t= ',num2str(current_time)));
        
        end



        h_grad=gradient(h_t,x(1:num_p));
        h_grad2=gradient(h_grad,x(1:num_p));
        h_grad3=gradient(h_grad2,x(1:num_p));
        h_grad4=gradient(h_grad3,x(1:num_p));

        h_grad_mult=h_grad.*h_grad2.*h_grad3.*h_grad4;





        if mobility_slop==0


            [~, r_indice_init]=max((h_grad2));
            indice_r=r_indice_init;
            r_indice_init2=r_indice_init+5;


                    xr_1=x(indice_r);

        elseif mobility_slop==1

            for i=1:num_p-2
                if h_grad2(i)>0
                    break
                end
            end

            r_indice_init=i;
            indice_r=r_indice_init;
            xr_1=x(r_indice_init);
        end


     


        if plot_graph==1

            plot(x(indice_r),h_t(indice_r),'or','MarkerSize',6,'MarkerFaceColor','none')

            plot(x(r_indice_init),h_t(r_indice_init),'og','MarkerSize',6,'MarkerFaceColor','g')


        hold off;

        subplot(2,3,2);


        plot(x(1:num_p),h_t,'.-','DisplayName', strcat('h at t= ',num2str(current_time)));
        hold on;
        plot(x(1:num_p),h_n,'DisplayName', strcat('h_{analytical} at t= ',num2str(current_time)));

        xlim([x(r_indice_init-5) x(r_indice_init+5)]);
        ylim([0 abs(h_t(r_indice_init)+b)+eps]);
        grid on
        plot(x(indice_r),h_t(indice_r),'or','MarkerSize',6,'MarkerFaceColor','r')

        plot(x(r_indice_init),h_t(r_indice_init),'og','MarkerSize',6,'MarkerFaceColor','g')

        hold off
        pause(0.0001) 
        end


        if abs((current_time-time_step(qq))/time_step(qq)) <1e-8
            

            [fid_s,msg] = fopen(sprintf('output/%d.txt',current_time),'wt');
            assert(fid_s>=3,msg)
            fprintf(fid_s,'%.18e %.18e %.18e\n',del_t,xr_1,xr_2);
            for i=1:num_p
                fprintf(fid_s,'%.18e %.18e %.18e %.18e %d %d %d %d\n',h_n(i),h_t(i),x(i),h_rr(i),grid_neigh(i,:));
            end
            
            fclose(fid_s);
            qq=qq+1;
            %break
        end
        
        
        
        
        
        h_temp_1=h_temp;
        h_temp=h_t;
        
        if timer_counter> 30 && nn_old < 10 && del_t<max_del_t && del_temp==0 && max(err) <0.008
            timer_counter=timer_counter-3;
            del_t_old=del_t;
            del_t=1.1*del_t; 
        elseif del_temp~=0
            del_t_old=del_t;
            del_t=del_temp; 
            del_temp=0;
        elseif nn_old>9
            nn_old=nn_old-1;

        end
        
        
         if qq < tt && (current_time+del_t)>time_step(qq)
            del_temp=del_t;
            del_t=time_step(qq)-current_time;
            current_time=time_step(qq);
        else
            current_time=current_time+del_t;
         end
        

        
        

        clear h_grad h_grad2;

        if plot_graph==1 
            subplot(2,3,3);
            hold on;
            yyaxis left
            plot (current_time,del_t,'b.');
            yyaxis right
            plot (current_time,max(err),'r.');
            


            subplot(2,3,4);
            hold on;
            %plot (current_time,max(h_t),'bo');
            plot (current_time,min(h_t),'b.');
        
        
            subplot(2,3,5);
            hold on;
            %plot (current_time,max(h_t),'b.');
            plot (current_time,cond_num,'b.');


            subplot(2,3,6);            

            plot(x,sum_phi_der,'.-')
        
            



        
        
        end
        l=l+1;
        nn=0;
    else
        del_t=0.9*del_t;
        nn=nn+1;
        nn_old=nn;
        timer_counter=0;
        h_t=h_temp;
        fprintf("number of failed iterations = %d \n",nn)

        
    end
    
        
    
    
    
    
end


[fid_s,msg] = fopen(sprintf('output/%d.txt',current_time),'wt');
assert(fid_s>=3,msg)
fprintf(fid_s,'%.18e %.18e %.18e\n',del_t,xr_1,xr_2);
for i=1:num_p
    fprintf(fid_s,'%.18e %.18e %.18e %.18e %d %d %d %d\n',h_n(i),h_t(i),x(i),h_rr(i),grid_neigh(i,:));
end

fclose(fid_s);

[fid_sg,msg] = fopen(sprintf('output_general.txt'),'wt');
assert(fid_sg>=3,msg)

for i=1:length(time_index)
    fprintf(fid_sg,'%.18e %.18e %.18e\n',time_index(i),h_min(i),error_vol(i));
end

fclose(fid_sg);






