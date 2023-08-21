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

case_path='/home/ISTO/mnorouzisadeh/cluster/Article_models/loc_s/b3/';
case_folder='hs_b4';
case_dir=strcat(case_path,case_folder);

%case_path='/home/ISTO/mnorouzisadeh/cluster/Article_models/ratio/b10/';
case_folder='hs_b10';
%case_dir=strcat(case_path,case_folder);

%% read input parameters

parameters=read_input(case_dir,'parameters');
dis_parameters=read_input(case_dir,'dis_parameters');



cluster=parameters(1); n_mobility=parameters(2); b=parameters(3);
plot_graph=parameters(4); laplacian_loop_mx=parameters(5); roughness_multiplyer=parameters(6);
alpha=parameters(7); del_t=parameters(8); max_del_t=parameters(9); end_time=parameters(10);
gravity=parameters(11); surf_tension=parameters(12); density=parameters(13);
xc=parameters(14); h_c=parameters(15); viscosity=parameters(16);
ex_f=parameters(17); x_dim=parameters(18); num_threads=parameters(19); start_from_zero=parameters(20);
del_x_init=parameters(21); init_profile=parameters(22); theta=parameters(23);potential_var=parameters(24:28);
if_refinement=parameters(29);period_bc=parameters(30);mobility_slop=parameters(31);slip_length=parameters(32);del_x_min=parameters(33);

dimen_case=1;

g_coef=xc^2/(surf_tension/(density*gravity));
time_c=(3.0*viscosity*xc^4)/(surf_tension*h_c^3);


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
% Check if it is correct? probably not!
if abs(zeta1+zeta2)<1e-6
    edl_coef=1*debye*xc^2/surf_tension;
else
    edl_coef=1*nb*avagad_num*kb*T*debye*xc^2/surf_tension;
end
vdw_coef=1*ham_cte1*xc^2/(2*pi()*surf_tension*h_c^4);
hyd_coef=1*w0_h*xc^2/(Lambda_2^2*surf_tension);


%% plot parameters
font_size=22;
line_width=2;

%% Reading the output results

fprintf('reading previous result to continue from last run \n')
output_folder = dir(strcat(case_dir,'/output/')); 
        [output_pre_f,~]=size(output_folder);
        temp_name={output_folder.name};

        temp_name=natsort(temp_name,'\d+\.?\d*([eE][-+]?\d+)?');
    
        num_timesteps=length(temp_name)-2;

        for i=1:num_timesteps
            file_address_sum=temp_name{i};
            
            output_1=read_input(strcat(case_dir,'/output/',file_address_sum),'output');
            
            dt(i)=time_c*output_1{1};
            x_r_lub_1(i)=xc*output_1{2};
            x_r_lub_2(i)=xc*output_1{8};
            
            h_lub_model{i}=h_c*output_1{3};
            
            
            %h_analytical{i}=h_c*output_1{7};


            
            x{i}=xc*output_1{4};

            j=1;
            while h_lub_model{i}(j)/h_c>b && j<length(h_lub_model{i}(j))-1
                j=j+1;
            end


            if h_lub_model{i}(j)/h_c <= b
                x_r_lub_1(i)=x{i}(j);
            end
            %plot(x{i}/xc,h_lub_model{i}/h_c,'.-')

            time(i)=time_c*output_1{5};
            num_p(i)=length(h_lub_model{i});

            x_t=(1.0+120.0*time(i)/time_c)^(0.2);
            h_analytical{i}=zeros(num_p(i),1);
            h_analytical{i}(:)=b*h_c;
            j=1;
%             while  abs(x{i}(j)/xc)<x_t
%                 h_analytical{i}(j)=((1.0+120.0*time(i)/time_c)^(-0.2)*(1-((x{i}(j)/xc)/x_t)^2)^2+b)*h_c;
%                 j=j+1;
% 
%             end

        end
        
        


        x_r_lub=x_r_lub_1;%(x_r_lub_1+x_r_lub_2)/2;
        %%
        indic=length(num_p);
        [min_temp_h, ind_min]=min(h_lub_model{indic}(:));

        [max_temp_h, ind_max]=max(h_lub_model{indic}(:));
    
        set(gcf,'Position',[100 100 1063 656])

        %plot(x{indic}(:)/x{indic}(ind_min),h_lub_model{indic}(:),'.-')
        plot(x{indic}(:),h_lub_model{indic}(:),'.-')
        %plot(x{indic}(:),gradient(h_lub_model{indic}(:),x{indic}(:)),'.-')

        xlab=xlabel('$x/x_{radius}$',interpreter='Latex');

        ylab=ylabel('$h ~ (m)$',interpreter='Latex');
        %ylim([0 2e-5])
        %xlim([0 3.5])

        grid(gca,'minor');




        %% read openfoam results 
if isequal(potential_var,[1 1 0 0 0]) || (isequal(potential_var,[1 1 1 1 1]))  % potential is only curvature and gradient
    openfoam_direct='/home/ISTO/mnorouzisadeh/cluster/VOF_lub_valid/scalesR/best_grid/';
    openfoam_case='ny100_film';

    file_address_openfoam_film=strcat(openfoam_direct,openfoam_case,'/film.mat');

    
    if exist(file_address_openfoam_film, 'file')
         load(file_address_openfoam_film);
    
         fprintf('film.mat is loaded in the memory  \n')
         % find the matching time
         time_match=[];
         time_match_index_vof=[];
         time_match_index_lub=[];
         for i=1:length(time_)
            ii=1;
            while abs(time_(i)-time(ii))/time_(i) >1e-3 && ii<length(time)
                ii=ii+1;
            end
            if abs(time_(i)-time(ii))/time_(i) <1e-3
                time_match=[time_match time_(i)];
                time_match_index_vof=[time_match_index_vof i];
                time_match_index_lub=[time_match_index_lub ii];
            end
         end

         % volume error
         for i=1:length(time_)
            vol_vof(i)=trapz(x_g,h(:,i));
            
            %vol_vof(i)=vol_vof(i)+(x_g(jjk)-x_g(jjk-1))*(h(,i)+h(,i))/2;   % This is one 2 case
         end
         
    else
        fprintf('info about openfoam film had not been created \n')
    end

    
    
end
%% plot del_X
for j=1:length(time)
    delx_=zeros(num_p(j)-1,1);
    
    
    x_select=x{j}/xc;
    for i=1:num_p(j)-1
        delx_(i)=x_select(i+1)-x_select(i);
    end
    delx_(num_p(j))=delx_(num_p(j)-1);
    
    plot(x_select,delx_,'.-')
    ylim([0.0001 0.02])
    xlim([0 4])
    grid on
    %pause
end
% figure()
% hold on
% grid on
% 
% plot(x_select,gradient(delx_,x_select),'.')
% pause()
% plot(x_select,gradient(gradient(delx_,x_select),x_select),'.')
% pause()
% plot(x_select,gradient(gradient(gradient(delx_,x_select),x_select),x_select),'.')
% pause()
% plot(x_select,gradient(gradient(gradient(gradient(delx_,x_select),x_select),x_select)),'.')


%% finding criteria
% grad1=gradient(h_lub_model{i}/h_c,x{i}/xc);
% grad2=gradient(grad1,x{i}/xc);
% grad3=gradient(grad2,x{i}/xc);
% grad4=gradient(grad3,x{i}/xc);
% 
% figure()
% plot(x{i}/xc,grad1.*grad2.*grad3.*grad4,'.-');
% grid on

%% error calculation
% 
% figure()
% 
% set(gcf,'Position',[100 100 800 800])
% if isequal(potential_var,[1 1 0 0 0])
%     error=zeros(1,length(time_match));
%     for i=1:length(time_match)
%         %error calculation based on VOF
%         ii_vof=time_match_index_vof(i);
%         ii_lub=time_match_index_lub(i);
%         vq = interp1(x{ii_lub}(:),h_lub_model{ii_lub}(:),x_g);
%         jjj=1;
%         while x{ii_lub}(jjj)<max(x_g)
%             error(i)=error(i)+(h_lub_model{i}(jjj)-vq(jjj))^2;
%             jjj=jjj+1;
%         end
%         error(i)=error(i)/(jjj-1);
%     end
% else
%     error=zeros(1,num_timesteps);
%     for i=1:num_timesteps
%         %error calculation based on analytical solution
% %         x_sel_error=linspace(0,x_dim*xc,100);
% %         y_lub_sel=spline(x{i},h_lub_model{i},x_sel_error);
% %         y_ana_sel=spline(x{i},h_analytical{i},x_sel_error);
% % 
% %         for jjj=1:100
% %             error(i)=error(i)+(y_lub_sel(jjj)-y_ana_sel(jjj))^2;
% %         end
% 
% %         for jjj=1:length(x{i})
% %             error(i)=error(i)+(h_lub_model{i}(jjj)-h_analytical{i}(jjj))^2;
% %         end
% %         error(i)=error(i)/length(x{i});
% 
%         
%         %error(i)=error(i)/100;
% 
%         error(i)=error(i)+(max(h_lub_model{i})-max(h_analytical{i}))^2;
%         
%     end
% end
% 
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% %yyaxis left
% if isequal(potential_var,[1 1 0 0 0])
%     error_plot=plot(time_match,error,'DisplayName', 'Means squere error');
% else
%     error_plot=plot(time,error,'DisplayName', 'Means squere error');
% end
% 
% 
% 
% error_plot.MarkerEdgeColor='[0.8500, 0.3250, 0.0980]	';
%     %vol_plot.Color=fluid_rough_color;
%     error_plot.LineStyle='-';
%     error_plot.LineWidth=3;
%     %vol_plot.MarkerSize=marker_size_1;
%     %vol_plot.Marker='s';
% 
% 
% xlab2=xlabel('time $(s) $','Interpreter','latex');
% %xlab1.FontSize=20;
% ylab2=ylabel('$\frac{1}{n}\sum_{i=1}^{n}\left(h_{{i}_{model}}-h_{{i}_{analytical}}\right)^2$','Interpreter','latex');
% %ylab2=ylabel('MSE','Interpreter','latex');
% 
% hold on
% 
% 
% % yyaxis right
% % dt_plot=plot(time,(dt),'DisplayName', 'Time-stepping');
% % 
% % dt_plot.MarkerEdgeColor='[0.8500, 0.3250, 0.0980]	';
% %     %dt_plot.Color=fluid_rough_color;
% %     dt_plot.LineStyle='-';
% %     dt_plot.LineWidth=3;
% %     %dt_plot.MarkerSize=marker_size_1;
% %     %dt_plot.Marker='s';
% % 
% % 
% % %  min_y_dt=round(log10(min(dt))); 
% % %  max_y_dt=round(log10(max(dt)));
% % %  dt_y_labels=min_y_dt:2:max_y_dt;
% % %  for i=1:length(dt_y_labels)
% % %      dt_yticks(i)=dt_y_labels(i);
% % %      yticklabels_s{i}=strcat('$ 10 ^{',num2str(dt_y_labels(i)),'}$');
% % %  end
% % 
% % %dt_yticks=[10^(min_y_dt) 10 100 ];
% % %yticks(dt_yticks);
% % %yticklabels(yticklabels_s)
% % ylab3=ylabel('$\Delta_t  \hspace{1mm} (s)$','Interpreter','latex');
% % ylim=[min(dt) max(dt)];
% grid on;
% 
% grid(gca,'minor');
% 
% %ylab.FontSize=40;
% set(gca,'FontSize',font_size,'LineWidth',line_width);
% 
% set(gca, 'XScale','linear')
% legend('Location', 'Best')
% 
% %savefig("../mass_balanc_time_step")
% %exportgraphics(gca,'results_pic/error.eps','Resolution',300) 
% %saveas(gcf,'results_pic/error_time_step','epsc')
% saveas(gcf,'results_pic/error','fig')
% 
% 
% 

        

        %% plot the results


%%% energy multiplier
e_dis=1;

y_res_energy=20;
my_marker_size=8;

vel_v=cell(1,num_timesteps);
vel_u=cell(1,num_timesteps);


%% velocity on the interface
% for i=1:num_timesteps
%     x_normal{i}=x{i}/x_r_lub(i);
% 
% end
% 
% 
% 
% for i=1:num_timesteps-1
%     %vel_v{i}(1)=(h_lub_model{i+1}(1)-h_lub_model{i}(1))/(time(i+1)-time(i));
%     %vel_u{i}(1)=0;
% 
%     h_2=interp1(x_normal{i+1}(:),h_lub_model{i+1}(:),x_normal{i}(:));
% 
%     vel_v{i}(:)=(h_2-h_lub_model{i}(:))/(time(i+1)-time(i));
%     vel_u{i}(:)= x_normal{i}(:)*(x_r_lub(i+1)-x_r_lub(i))/(time(i+1)-time(i));
% 
%     %for j=2:num_p(i)-1
%        
% 
% %         jj=1;
% %         if x_normal{i}(j)<1.5
% %             while x_normal{i+1}(jj) < x_normal{i}(j)
% %                jj=jj+1;
% %             end
% %             
% % 
% %             A=[x_normal{i+1}(jj) x_normal{i+1}(jj-1)];
% %             B=[h_lub_model{i+1}(jj) h_lub_model{i+1}(jj-1)];
% %             h_2=interp1(A,B,x_normal{i}(j));
% %             vel_v{i}(j)=(h_2-h_lub_model{i}(j))/(time(i+1)-time(i));
% %             vel_u{i}(j)= x_normal{i}(j)*(x_r_lub(i+1)-x_r_lub(i))/(time(i+1)-time(i));
% %         else
% %             vel_v{i}(j)=0;
% %             vel_u{i}(j)=0;
% %         end
%     %end
% end
% max_vel_v=0;
% max_vel_u=0;
% 
% for i=1:num_timesteps-1
%     if max(vel_v{i})>max_vel_v
%         max_vel_v=max(vel_v{i});
%     end
% 
%     if max(vel_u{i})>max_vel_u
%         max_vel_u=max(vel_u{i});
%     end
% 
% end
% for i=1:num_timesteps-1
%     vel_v{i}=(vel_v{i}./(1));
%     vel_u{i}=(vel_u{i}./(1));
% end
% figure()
%  for i=1:num_timesteps-1
%      plot(x{i}(1:num_p(i)-1),h_lub_model{i}(1:num_p(i)-1))
%      hold on
%      quiver(x{i}(1:num_p(i))',h_lub_model{i}(1:num_p(i))',vel_u{i},vel_v{i},0)
%      xlim([0 x_dim*xc]);
%      ylim([0 h_c]);
%      pause(0.1)
%      hold off
%  end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% PLOT the energy 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number_curve=length(num_timesteps);
clr_mp=hot(num_timesteps);

% my_map = {'#32144f', '#82054a', '#ae334f', '#c06b57', '#d09c5d', '#dfcc64', '#eefa6a'};
% clr_mp = validatecolor(my_map, 'multiple');


colormap(clr_mp);

marker_shape='.';


plot(0,0)
set(gcf,'Position',[100 100 1500 500])
parentAxis=gca;

xlab=xlabel('$ X \hspace(1mm)(m)$','Interpreter','latex');
xlab.FontSize=font_size;
ylab=ylabel('$h \hspace(1mm) (m)$','Interpreter','latex');
%ylab.FontSize=40;
set(gca,'FontSize',font_size,'LineWidth',line_width);
set(gca, 'FontName', 'Times New Roman')




for i=1:num_timesteps
    
    i
    

    vol(i)=0;
    vol_an(i)=0;
    for j=1:num_p(i)
        if x{i}(j)==x_r_lub(i)
            break
        end
    end


%     for jjk=2:num_p(i)
%         if dimen_case==2
%             x2=x{i}(jjk);
%             x1=x{i}(jjk-1);
%             y2=h_lub_model{i}(jjk);
%             y1=h_lub_model{i}(jjk-1);
% 
%             slop_l=(y2-y1)/(x2-x1);
%             coef_line=y2-slop_l*x2;
%             %if slop_l < 1e-6
%                 vol(i)=vol(i)+abs(pi*(x2^2+x1^2)*(y2-y1)/2);  % This is axisymmetrica case
%             %else
%             %    vol(i)=vol(i)+pi/slop_l^2*(x2-x1)*((x2^2+x1^2+x1*x2)/3-coef_line*(x2+x1)+coef_line^2);  % This is axisymmetrica case
%             %end
%         elseif dimen_case==1
%             
%             vol(i)=vol(i)+(x{i}(jjk)-x{i}(jjk-1))*(h_lub_model{i}(jjk)+h_lub_model{i}(jjk-1))/2;   % This is one 2 case
%         end
% 
%     end
    vol(i)=trapz(x{i}(:),h_lub_model{i}(:));
    
    for jjk=2:num_p(i)
        if dimen_case==2
            x2=x{i}(jjk);
            x1=x{i}(jjk-1);
            y2=h_analytical{i}(jjk);
            y1=h_analytical{i}(jjk-1);

            slop_l=(y2-y1)/(x2-x1);
            coef_line=y2-slop_l*x2;
            %if slop_l < 1e-6
                vol_an(i)=vol(i)+abs(pi*(x2^2+x1^2)*(y2-y1)/2);  % This is axisymmetrica case
            %else
            %    vol(i)=vol(i)+pi/slop_l^2*(x2-x1)*((x2^2+x1^2+x1*x2)/3-coef_line*(x2+x1)+coef_line^2);  % This is axisymmetrica case
            %end
        elseif dimen_case==1
            vol_an(i)=vol_an(i)+(x{i}(jjk)-x{i}(jjk-1))*(h_analytical{i}(jjk)+h_analytical{i}(jjk-1))/2;   % This is one 2 case
        end

    end


        h_x{i}=gradient(h_lub_model{i},x{i});
        h_max(i)=max(h_lub_model{i});
        [mag_thta_infl, indice_inflc]=max(atand(-h_x{i}));
        theta(i)=mag_thta_infl;
        inflection_point(i)=indice_inflc;
        %theta_loc(i,:)=atand(-h_x);

        h_xx{i}=gradient(h_x{i},x{i});

        
        
    
    
    lub_model=plot(x{i},h_lub_model{i},marker_shape);
    hold on
    points_top_sym=plot(-x{i},h_lub_model{i},marker_shape);
    lub_model.MarkerEdgeColor='k';
    %lub_model.MarkerSize=my_marker_size;

    points_top_sym.MarkerEdgeColor='k';
    %points_top_sym.MarkerSize=my_marker_size;


    analytical_=plot(x{i},h_analytical{i},'-');
    analytical_.Color=clr_mp(i,:);
    analytical_sym=plot(-x{i},h_analytical{i},'-');
    analytical_sym.Color=clr_mp(i,:);


    xlim([-x_dim*xc x_dim*xc])
    ylim([0 1*h_c]);

    set(gca,'DataAspectRatio',[9 1 1])

    %plot(x,h_lub_model_an,'.b','DisplayName', name);
    %legend;

    


    M=h_lub_model{i}(j);
    magnification_x=40;
    magnification_y=5;
   
    
    %b_h_init_p=plot(x_init,h_init);
    %b_h_init_p.LineWidth=line_size;

    %x_min_1=x{i}(j)*(1-2.0/(magnification_x*10));
    %x_max_1=x{i}(j)*(1+2.0/(magnification_x*10));

    x_min_1=x_r_lub(i)*(1-2.0/(magnification_x*10));
    x_max_1=x_r_lub(i)*(1+2.0/(magnification_x*10));


    y_min_1=0;
    y_max_1=M*(1.5+0.1/magnification_y);
    
    r_zoom=rectangle('Position',[x_min_1 y_min_1 (x_max_1-x_min_1) (y_max_1-y_min_1)]);

     quiver(x_min_1, y_min_1,(2*0.65-x_min_1), (0.35*1.5-y_min_1),0,'r');

     xlab=xlabel(' X (m)');
    xlab.FontSize=font_size;
    ylab=ylabel('h (m)');
    %ylab.FontSize=40;
    set(gca,'FontSize',font_size,'LineWidth',line_width);
    set(gca, 'FontName', 'Times New Roman')
    leg_str=strcat('OpenFOAM result at t=',num2str(time(i)),' s');
    %l_legend=legend([points_top_sym  ener_plot_2d],{leg_str, 'Lubrication model'});
    l_legend.Location='northwest';

    hold off

    %%%%%%%
    f=axes('position',[.62 .5 .22 .25]);
    box on % put box around new pair of axes


    hold on

    points_top_z=plot(x{i},h_lub_model{i},'r.-');
    %analytical_z=plot(x{i},h_analytical{i},'b-');
        
    
                   

    points_top_z.MarkerEdgeColor='k';
    points_top_z.MarkerSize=my_marker_size;




    xlim([x_min_1 x_max_1])
    ylim([0 y_max_1+eps]);

    %hold on;
    grid on;
    grid(gca,'minor');




    %%%%%%%%


    %pause()
    %set(gca,'FontSize',30,'LineWidth',1.5);
    %if mod(ind_t,1)==0
    %    ll=ll+1;
        Frame(i) = getframe(gcf) ;

    hold off
    if (i ~= num_timesteps)
        delete(points_top_z)
        delete(f)
    end


end

% %% plot vdw validation
% 
% marker_shape='.-';
% 
% 
% plot(0,0)
% set(gcf,'Position',[100 100 1500 500])
% parentAxis=gca;
% 
% xlab=xlabel('$ X \hspace(1mm)(m)$','Interpreter','latex');
% xlab.FontSize=font_size;
% ylab=ylabel('$h \hspace(1mm) (m)$','Interpreter','latex');
% %ylab.FontSize=40;
% set(gca,'FontSize',font_size,'LineWidth',line_width);
% set(gca, 'FontName', 'Times New Roman')
% 
% 
% 
% 
% for i=1:num_timesteps
%     
%     lub_model=plot(x{i},h_lub_model{i},marker_shape);
%     h_min(i)=min(h_lub_model{i});
%     hold on
%     lub_model.MarkerEdgeColor='k';
%     %lub_model.MarkerSize=my_marker_size;clear cle
% 
%     %points_top_sym.MarkerSize=my_marker_size;
%     %xlim([-x_dim*xc x_dim*xc])
%     ylim([0 2]);
% 
%     
%      xlab=xlabel(' X ');
%     xlab.FontSize=font_size;
%     ylab=ylabel('h ');
%     %ylab.FontSize=40;
%     set(gca,'FontSize',font_size,'LineWidth',line_width);
%     set(gca, 'FontName', 'Times New Roman')
%     leg_str=strcat('t=',num2str(time(i)/time_c));
%     l_legend=legend(lub_model,{leg_str});
%     l_legend.Location='northwest';
% 
%     hold off
% 
%     
%     grid on;
%     grid(gca,'minor');
% 
% 
% 
%         Frame(i) = getframe(gcf) ;
% 
%     hold off
%     if (i ~= num_timesteps)
%         delete(points_top_z)
%         delete(f)
%     end
% 
% 
% end
% figure()
% plot(-log(4-time/time_c),-log(h_min))
% grid on
%% Whole droplet image multiple timesteps
figure()
%set(gcf,'position',[100,100,1600,800])
%set(gcf,'Position',[100 100 874 656])
set(gcf, 'Units', 'Inches', 'Position', [2, 2, 13.5, 5], 'PaperUnits', 'Inches', 'PaperSize', [13.5, 5])
%set(gcf, 'Units', 'Inches', 'Position', [2, 2, 6.75, 5], 'PaperUnits', 'Inches', 'PaperSize', [6.75, 5])

num_curve_plot=5;



clr_mp=summer(7);

my_map = {'#32144f', '#8e0349', '#b64b52', '#ca895b', '#dcc262', '#eefa6a'};
%clr_mp = validatecolor(my_map, 'multiple');

%colormap(flipud(clr_mp))
colormap(clr_mp);



time_log=log10(time);
% y__=1;
% i=9;
% while i<num_timesteps
%     y__= [y__,i ];
%     i=i+9;
% end
%y__=[90 101 103  104  ];
%y__=[1 50 82  99];
%y__=[18 45 50 55 61 64];
y__=[1 4  9 12 14 16];
%y__=[1 4 5 9 11];
for i=1:length(y__)
    

    if isequal(potential_var,[1 1 0 0 0])   % potential is only curvature and gradient
        vof_i=time_match_index_vof(y__(i));
        lub_i=time_match_index_lub(y__(i));
        analytical_=plot(x_g,h(:,vof_i),'-');
         hold on
         lub_model=scatter(x{lub_i}(1:5:end),h_lub_model{lub_i}(1:5:end),'o',"filled");
         analytical_sym=plot(-x_g,h(:,vof_i),'-');
        
        lub_model_sym=scatter(-x{lub_i}(1:5:end),h_lub_model{lub_i}(1:5:end),'o',"filled");
        hl=legend('VOF simulation','Lubrication model');

        analytical_.Color=clr_mp(i,:);
    
    analytical_sym.Color=clr_mp(i,:);

    analytical_.LineWidth=2.5;
    analytical_sym.LineWidth=2.5;
    elseif isequal(potential_var,[1 1 1 1 1]) || isequal(potential_var,[1 1 2 1 1])  % potential is only curvature and gradient
        
        lub_i=y__(i);
       
         hold on
         lub_model=scatter(x{lub_i}(1:10:end),h_lub_model{lub_i}(1:10:end),'o',"filled");

        
        
        lub_model_sym=scatter(-x{lub_i}(1:10:end),h_lub_model{lub_i}(1:10:end),'o',"filled");

        hl=legend('VOF simulation','Lubrication model');

         lub_model_l=plot(x{lub_i},h_lub_model{lub_i},'-.');
         lub_model_sym_l=plot(-x{lub_i},h_lub_model{lub_i},'-.');
         lub_model_l.Color=clr_mp(i,:);
         lub_model_sym_l.Color=clr_mp(i,:);

         lub_model_l.LineWidth=1.5;
         lub_model_sym_l.LineWidth=1.5;
    else
        analytical_=plot(x{y__(i)},h_analytical{y__(i)},'-');
         hold on
         lub_model=scatter(x{y__(i)}(1:10:end),h_lub_model{y__(i)}(1:10:end),'o',"filled");
         analytical_sym=plot(-x{y__(i)},h_analytical{y__(i)},'-');
        
        lub_model_sym=scatter(-x{y__(i)}(1:10:end),h_lub_model{y__(i)}(1:10:end),'o',"filled");
        hl=legend('Analytical solution','Lubrication model');

        analytical_.Color=clr_mp(i,:);
    
    analytical_sym.Color=clr_mp(i,:);

    analytical_.LineWidth=3;
    analytical_sym.LineWidth=3;
    end

   
    

         
    

    
    alpha_marker=1;
    

    lub_model.MarkerEdgeColor='w';
    lub_model.MarkerFaceAlpha=alpha_marker;
    lub_model.MarkerEdgeAlpha=alpha_marker;
    lub_model.MarkerFaceColor='k';
    lub_model.LineWidth=1.5;
    lub_model.SizeData=7;
    distfromzero = abs(1-x{y__(i)}(1:10:end)/(x_r_lub(y__(i)))).^.2;
    lub_model.AlphaData = distfromzero;
    lub_model.MarkerFaceAlpha = 'flat';
    lub_model.MarkerEdgeAlpha = 'flat';


    lub_model_sym.MarkerEdgeColor='w';
    lub_model_sym.MarkerFaceColor='k';
    lub_model_sym.MarkerFaceAlpha=alpha_marker;
    lub_model_sym.MarkerEdgeAlpha=alpha_marker;
    lub_model_sym.LineWidth=1.5;
    lub_model_sym.SizeData=7;
    distfromzero = abs(1-x{y__(i)}(1:10:end)/(x_r_lub(y__(i)))).^.2;
    lub_model_sym.AlphaData = distfromzero;
    lub_model_sym.MarkerFaceAlpha = 'flat';
    lub_model_sym.MarkerEdgeAlpha = 'flat';


    
    %pause
end
xlim([-2.5*xc 2.5*xc])
ylim([0 1.1*h_c]);

set(findall(gcf,'-property','FontSize'),'FontSize',14)


aa=legend('Location','northeast' );
aa.FontSize=18;
set(hl, 'Interpreter','latex')

xlab=xlabel(' X (m)');
xlab.FontSize=18;
ylab=ylabel('h (m)');
ylab.FontSize=18;


grid on;
%grid(gca,'minor');
hcb=colorbar();
caxis([min(time) max(time)])

colorTitleHandle = get(hcb,'Title');
titleString = 'Time(s)';
set(colorTitleHandle ,'String',titleString);
%axis equal

set(gca,'DataAspectRatio',[10 1 1])
set(gca,'LineWidth',1.5);
set(gca, 'FontName', 'Times New Roman')
%set(gca,'xminorgrid','on','yminorgrid','on')

ax = gca;



exportgraphics(gca,'results_pic/whole_spreading.eps','Resolution',300) 
saveas(gcf,strcat('results_pic/whole_spreading ',case_folder,'.fig'),'fig')
%%
for i=1:num_timesteps
    h_t=h_lub_model{i};
    x_t=x{i};

    h_grad=gradient(h_t,x_t);
    h_grad2=gradient(h_grad,x_t);
    h_grad3=gradient(h_grad2,x_t);
    h_grad4=gradient(h_grad3,x_t);


    [~, r_indice_init1]=max(abs(h_grad3));
    [~, r_indice_init2]=max(abs(h_grad));
    r_indice_init=max(r_indice_init1,r_indice_init2);
    x_r_lub(i)=x_t(r_indice_init);
end



%% spreading only

figure()
%set(gcf,'Position',[100 100 800 800])
set(gcf, 'Units', 'Inches', 'Position', [2, 2, 6.75, 5], 'PaperUnits', 'Inches', 'PaperSize', [6.75, 5])
%plot(time,x_r_vof,'sk','DisplayName','VOF');
grid on;
grid(gca,'minor');
xlab2=xlabel('time $(s) $','Interpreter','latex');
%xlab1.FontSize=20;
ylab2=ylabel('x_r (m)','Interpreter','tex');
%ylab.FontSize=40;
set(gca,'FontSize',font_size,'LineWidth',1.5);
hold on
xr_lub=plot(time,x_r_lub,'.-','DisplayName',strcat('Lubrication',case_folder));


if isequal(potential_var,[1 1 0 0 0]) || isequal(potential_var,[1 1 1 1 1])
    %xr_ana=plot(time_,x_r_vof,'DisplayName','VOF');
else
    time_=0:0.1:200;
    %xr_ana=plot(time_,xc*(1.0+120.0.*time_/time_c).^(0.2),'DisplayName','Analytical');
end


xr_lub.Marker='.';
xr_lub.MarkerSize=my_marker_size;
xr_lub.LineWidth=1.5;
xr_lub.MarkerEdgeColor='#6cb0d6';
xr_lub.MarkerFaceColor='#226e9c';
xr_lub.LineStyle='-';

%xr_ana.LineWidth=2;
%xr_ana.Color='#0D4A70';
%xr_ana.LineStyle='-';

xlim([0 max(time)])
legend('Location','northwest')
set(findall(gcf,'-property','FontSize'),'FontSize',14)


%set(gca,'FontSize',font_size,'LineWidth',line_width);
set(gca, 'FontName', 'Times New Roman')
saveas(gcf,strcat('results_pic/spreading_radius_only_',case_folder,'.fig'),'fig')





%% Spreading
for i=1:num_timesteps
    h_lun_0(i)=h_lub_model{i}(1);
end

figure()
set(gcf,'Position',[100 100 800 800])
%plot(time,x_r_vof,'sk','DisplayName','VOF');
grid on;
grid(gca,'minor');
xlab2=xlabel('time $(s) $','Interpreter','latex');
%xlab1.FontSize=20;
ylab2=ylabel('$x_{droplet}$ (m)','Interpreter','latex');
%ylab.FontSize=40;
set(gca,'FontSize',font_size,'LineWidth',1.5);
hold on
yyaxis left
xr_lub=plot(time,x_r_lub,'DisplayName','Lubrication');


if isequal(potential_var,[1 1 0 0 0]) % || (isequal(potential_var,[1 1 1 1 1]))
    xr_ana=plot(time_,x_r_vof,'DisplayName','VOF');
else
    time_=0:0.1:200;
    xr_ana=plot(time_,xc*(1.0+120.0.*time_/time_c).^(0.2),'DisplayName','Analytical');
end



set(gca,'YColor','#32144f')

yyaxis right
hr_lub=plot(time,h_lun_0,'DisplayName','Lubrication');
if isequal(potential_var,[1 1 0 0 0])
    hr_ana=plot(time_,h_vof_0,'DisplayName','VOF');
    legend('Lubrication model','VOF model')
else
    time_=0:0.1:200;
    hr_ana=plot(time_,h_c*(1.0+120.0.*time_/time_c).^(-0.2),'DisplayName','Analytical');
    legend('Lubrication model','Analytical solution')
end

ylab3=ylabel('$h_{0}$ (m)','Interpreter','latex');

xr_lub.Marker='s';
xr_lub.MarkerSize=my_marker_size;
xr_lub.LineWidth=1.5;
xr_lub.MarkerEdgeColor='#6cb0d6';
xr_lub.MarkerFaceColor='#32144f';
xr_lub.LineStyle='none';

hr_lub.Marker='s';
hr_lub.MarkerSize=my_marker_size;
hr_lub.LineWidth=1.5;
hr_lub.MarkerEdgeColor='#f9d8e6';
hr_lub.MarkerFaceColor='#ae334f';
hr_lub.LineStyle='none';


xr_ana.LineWidth=2;
xr_ana.Color='#32144f';
xr_ana.LineStyle='-';

hr_ana.LineWidth=2;
hr_ana.Color='#ae334f';
hr_ana.LineStyle='-';

%set(gca, 'XScale', 'log')
set(gca,'XMinorTick','on')
%min_t=min(time);
%pw_min = ceil(log10(min_t));
%res_log = 10^(pw_min-1);
%max_t=max(time);
%for i_tick=1:floor(log10(max_t/min_t))
%    xlog_tic(i_tick)=floor(min_t/res_log)*res_log*10^(i_tick);
%end

xlim([0 20])
%xlim([0 80])
legend('Location','northwest')

%set(gca, 'XTick', xlog_tic);
set(gca,'YColor','#ae334f')

set(gca,'FontSize',font_size,'LineWidth',line_width);
set(gca, 'FontName', 'Times New Roman')

%savefig("../spreading_radius")
exportgraphics(gca,'results_pic/spreading_radius_height.eps','Resolution',300) 

saveas(gcf,'results_pic/spreading_radius','epsc')
saveas(gcf,'results_pic/spreading_radius','fig')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%% MASS BALANCE plot  +++++ Deta t plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

figure()
set(gcf,'Position',[100 100 800 800])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
yyaxis left
vol_plot=semilogx(time,(vol)./vol(1),'DisplayName', 'Relative volume');
hold on
if isequal(potential_var,[1 1 0 0 0])
    vol_plot_an=semilogx(time_,(vol_vof)./vol_vof(1),'DisplayName', 'Relative volume analytical');
else
    vol_plot_an=semilogx(time,(vol_an)./vol_an(1),'DisplayName', 'Relative volume analytical');
end

vol_plot.MarkerEdgeColor='[0.8500, 0.3250, 0.0980]	';
    %vol_plot.Color=fluid_rough_color;
    vol_plot.LineStyle='-';
    vol_plot.LineWidth=3;
    %vol_plot.MarkerSize=marker_size_1;
    %vol_plot.Marker='s';


xlab2=xlabel('time $(s) $','Interpreter','latex');
%xlab1.FontSize=20;
ylab2=ylabel('$\frac{V}{V_{initial}}$','Interpreter','latex');

yyaxis right
dt_plot=semilogx(time,log10(dt),'DisplayName', 'Time-stepping');

dt_plot.MarkerEdgeColor='[0.8500, 0.3250, 0.0980]	';
    %dt_plot.Color=fluid_rough_color;
    dt_plot.LineStyle='-';
    dt_plot.LineWidth=3;
    %dt_plot.MarkerSize=marker_size_1;
    %dt_plot.Marker='s';


 min_y_dt=round(log10(min(dt))); 
 max_y_dt=round(log10(max(dt)));
 dt_y_labels=min_y_dt:2:max_y_dt;
 for i=1:length(dt_y_labels)
     dt_yticks(i)=dt_y_labels(i);
     yticklabels_s{i}=strcat('$ 10 ^{',num2str(dt_y_labels(i)),'}$');
 end

%dt_yticks=[10^(min_y_dt) 10 100 ];
yticks(dt_yticks);
yticklabels(yticklabels_s)
ylab3=ylabel('$\Delta_t  \hspace{1mm} (s)$','Interpreter','latex');

grid on;

grid(gca,'minor');

%ylab.FontSize=40;
set(gca,'FontSize',font_size,'LineWidth',1.2);


legend('Location', 'Best')

%savefig("../mass_balanc_time_step")

saveas(gcf,'results_pic/mass_balanc_time_step','epsc')
saveas(gcf,'results_pic/mass_balanc_time_step','fig')



%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%%%%        Contact angle vs velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
 %yyaxis left

%theta_p=semilogx(time,theta,'.-','DisplayName', 'Contact angle');
theta_p=plot(time,theta,'.-','DisplayName', strcat('Contact angle ',case_folder));

theta_p.MarkerEdgeColor='[0.8500, 0.3250, 0.0980]	';
    %theta_p.Color=fluid_rough_color;
    theta_p.LineStyle='-';
    theta_p.LineWidth=3;
    theta_p.MarkerSize=20;
    %theta_p.Marker='s';


xlab2=xlabel('time $(s) $','Interpreter','latex');
%xlab1.FontSize=20;
ylab2=ylabel('$\theta$','Interpreter','latex');
xlim([0 max(time)])
%hold on
%yyaxis right

saveas(gcf,strcat('results_pic/CA_',case_folder,'.fig'),'fig')


figure()
clear x_r_lub2 time2 theta_fit2 vel
j=0;
for i=1:length(x_r_lub)
    if time(i)>1e-12
        j=j+1;
        x_r_lub2(j)= x_r_lub(i);
        time2(j)=time(i);
        theta_fit2(j)=theta(i);
    end

end

vel=smooth(gradient(smooth(x_r_lub2),time2));
%vel=smooth(diff(x_r_lub2)./diff(time2));

%vel_plot=semilogy(time2,(vel),'.-','DisplayName', 'Velocity');
vel_plot=plot(time2,(vel),'o-','DisplayName', strcat('Velocity ',case_folder));
%vel_plot=plot(time(1:end-1),(vel),'.-','DisplayName', 'Velocity');
vel_plot.MarkerEdgeColor='[0.8500, 0.3250, 0.0980]	';
    %vel_plot.Color=fluid_rough_color;
    vel_plot.LineStyle='-';
    vel_plot.LineWidth=3;
    vel_plot.MarkerSize=20;
    %vel_plot.Marker='s';


min_y_dv=round(log10(min(vel(vel>0)))); 
 max_y_dv=round(log10(max(vel)));
 dt_y_v_labels=min_y_dv:2:max_y_dv;
 for i=1:length(dt_y_v_labels)
     dt_v_yticks(i)=dt_y_v_labels(i);
     yticklabels_v{i}=strcat('$ 10 ^{',num2str(dt_y_v_labels(i)),'}$');
 end

%dt_yticks=[10^(min_y_dt) 10 100 ];
%yticks(dt_v_yticks);
yticklabels(yticklabels_v)
ylab3=ylabel('$ u_{r}  \hspace{1mm} (m/s)$','Interpreter','latex');



grid on;

grid(gca,'minor');

%ylab.FontSize=40;
set(gca,'FontSize',font_size,'LineWidth',line_width);



legend('Location', 'Best')


%savefig("../CA_velocity")
saveas(gcf,'results_pic/CA_velocity','epsc')
saveas(gcf,strcat('results_pic/CA_velocity',case_folder,'.fig'),'fig')

vel=(viscosity/surf_tension)*vel;
clear vel_fit theta_fit time_fit
j=0; 
% for i=1:length(vel)
%     if time2(i)>1e-6
%         if (vel(i)>0)
%             j=j+1;
%             vel_fit(j)= vel(i);
%             theta_fit(j)=theta_fit2(i);
%             time_fit(j)=time2(i);
%             
%         end
%     end
% 
% end

vel_fit= vel';
theta_fit=theta_fit2;
time_fit=time2;


theta_s=min(theta);

myfittype = fittype('a.*(theta_fit).^b',...
    'dependent',{'vel_fit'},'independent',{'theta_fit'},...
    'coefficients',{'a','b'});




% myfittype = fittype('a*exp(b*(theta_fit-theta_s))',...
%     'dependent',{'vel_fit'},'independent',{'theta_fit'},...
%     'coefficients',{'a','b','theta_s'});

%coeffnames(myfittype);
myfit = fit(theta_fit',vel_fit',myfittype,'StartPoint',[1e-8,3],'Lower',[1e-14,2.7],'Upper',[1e-1,3.2]);

% myfittype_2 = fittype('(a.*h_max+c)*(theta).^b',...
%     'dependent',{'vel'},'independent',{'theta','h_max'},...
%     'coefficients',{'a','b','c'});
% myfit2 = fit([theta', h_max'],vel_fit',myfittype_2,'StartPoint',[1e-5,3,1e-9],'Lower',[1e-6,2.8,1e-10],'Upper',[1000,3.3,1e-5]);
% 
%  myfittype_3 = fittype('(a.*(theta_fit).^b)/(1+c.*(theta_fit).^d)',...
%      'dependent',{'vel_fit'},'independent',{'theta_fit'},...
%      'coefficients',{'a','b','c','d'});
% 
% myfit3 = fit(theta_fit',vel_fit',myfittype_3,'StartPoint',[1e-4,3,1e-2,2],'Lower',[1e-11,2.5,1e-7,0.5],'Upper',[1e-1,3.5,2e-1,2.5]);

% figure()
% 
% plot(myfit,theta_fit',vel_fit')
% set(gca, 'YScale', 'log')
% %set(gca, 'XScale', 'log')
% 
% xlab2=xlabel('$\theta$','Interpreter','latex');
% %xlab1.FontSize=20;
% ylab2=ylabel('$u_{r}$ $(\frac{m}{s})$','Interpreter','latex');
% 
% grid on;
% 
% grid(gca,'minor');
% 
% set(gca,'FontSize',15,'LineWidth',1.2);
% 
% legend('Location', 'Best')
% hl=legend('simulation',['$u_{r}=',num2str(myfit.a),'\theta^{',num2str(myfit.b),'}$']);
% 
% set(hl, 'Interpreter','latex')



%x = linspace(0,2*pi,50);
%y = sin(x) + randi(50,1,50);
%c = linspace(1,10,length(theta_fit'));
figure()
set(gcf,'Position',[100 100 800 800])
s=scatter(theta_fit,vel_fit,[],log10(time_fit),'filled','o','MarkerEdgeColor','k','DisplayName',case_folder);
s.SizeData = 100;
colorbar
colormap summer(256)
hold on
x_fit=min(theta')*0.01:0.1:max(theta');
y_fit=myfit.a*x_fit.^myfit.b;
%y_fit2=(myfit2.a.*h_max).*(theta.^myfit2.b);
%y_fit3=(myfit3.a*x_fit.^myfit3.b)./(1-myfit3.c*x_fit.^myfit3.d);
%plot(theta_fit,vel_fit)
plot(x_fit,y_fit,'DisplayName',strcat(case_folder, 'Ca=',num2str(myfit.a),'\times \theta^{',num2str(myfit.b),'}'))
%plot(theta,y_fit2,'DisplayName',strcat('u_{r}=',num2str(myfit2.a),'\theta^{',num2str(myfit2.b),'}'))
%plot(x_fit,y_fit3,'DisplayName',strcat('u_{r}=',num2str(myfit3.a),'\theta^{',num2str(myfit3.b),'}/1-',num2str(myfit3.c),'theta^{',num2str(myfit3.d),'}'))


hl=legend('Location', 'Best');
set(hl, 'Interpreter','tex')

xlab2=xlabel('\theta','Interpreter','tex');
%ylab2=ylabel('$u_{r}$ $(\frac{m}{s})$','Interpreter','tex');
ylab2=ylabel('Ca','Interpreter','tex');
grid on;
%grid(gca,'minor');

set(gca,'FontSize',font_size,'LineWidth',line_width);
set(gca, 'YScale', 'log')
%savefig("../CA_velocity_tanner")

saveas(gcf,strcat('results_pic/CA_velocity_tanner_',case_folder),'epsc')
saveas(gcf,strcat('results_pic/CA_velocity_tanner_',case_folder,'.fig'),'fig')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5



% figure()
% %Plot contact angle versus distance to contact line
% plot_ra_b=10;       %plot range before inflection point
% plot_ra_a=10;
% dist_ax=-plot_ra_b:plot_ra_a;
% for i=1:length(time)
%     for j=1:length(dist_ax)
%         thta_3d(i,j)=theta_loc(i,inflection_point(i)+dist_ax(j));
%     end
% end
% for i=1:length(time)
%     dist_3d(i,:)=dist_ax;
%     for j=1:length(dist_ax)
%         vel_3d(i,j)=vel(i);
%     end
% end
% 
% s=surf(vel_3d,dist_3d,thta_3d);
% colorbar
% %s=surf(Xq,Yq,Vq,'FaceAlpha',0.5);
% s.EdgeColor = 'none';
% s.FaceColor='interp';
% s.FaceLighting='gouraud';
% xlabel('Contact line velocity');
% %caxis([1 700])
% ylabel('distance to inlfection point');
% zlabel('Contact angle');
% 
% set(gca, 'XScale', 'log')
% hold on
% y_gen=zeros(1,length(vel_fit));
% scatter3(vel_fit',y_gen',theta_fit',[],log(time_fit'*time_c),'filled','o','MarkerEdgeColor','k')
% plot3(y_fit,zeros(1,length(y_fit)),x_fit,'DisplayName',strcat('u_{r}=',num2str(myfit.a),'\theta^{',num2str(myfit.b),'}'))

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% writerObj = VideoWriter(strcat(case_folder,'evolution.avi'),'Uncompressed AVI');
%  %writerObj.Quality=100;
%  %writerObj.VideoCompressionMethod='H.264';
%   writerObj.FrameRate = 10;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:num_timesteps
%     % convert the image to a frame
%     frame = Frame(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);

% 
