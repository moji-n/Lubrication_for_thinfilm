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
case_path='/home/ISTO/mnorouzisadeh/cluster/Article_models/ratio/b10/';
case_folder='hs_b10_s';
case_dir=strcat(case_path,case_folder);

%% plot parameters
font_size=22;
line_width=2;
my_marker_size=7;

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
if_refinement=parameters(29);%;period_bc=parameters(30);mobility_slop=parameters(31);slip_length=parameters(32);del_x_min=parameters(33);

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
set(gcf, 'Units', 'Inches', 'Position', [2, 2, 13.5, 3], 'PaperUnits', 'Inches', 'PaperSize', [13.5, 3])
parentAxis=gca;

xlab=xlabel('$ X \hspace(1mm)(m)$','Interpreter','latex');
xlab.FontSize=font_size;
ylab=ylabel('$h \hspace(1mm) (m)$','Interpreter','latex');
%ylab.FontSize=40;
set(gca,'FontSize',font_size,'LineWidth',line_width);
set(gca, 'FontName', 'Times New Roman')




for i=1:num_timesteps
%for j=1:length(time_match_index_lub)
    
    

    %i=time_match_index_lub(j);
    %ii=time_match_index_vof(j);
    
    ii=i;
    

for j=1:num_p(i)
        if x{i}(j)==x_r_lub(i)
            break
        end
    end

        h_x{i}=gradient(h_lub_model{i},x{i});
        h_max(i)=max(h_lub_model{i});
        [mag_thta_infl, indice_inflc]=max(atand(-h_x{i}));
        theta(i)=mag_thta_infl;
        inflection_point(i)=indice_inflc;
        %theta_loc(i,:)=atand(-h_x);

        h_xx{i}=gradient(h_x{i},x{i});

        area(x{i},h_lub_model{i},FaceColor="#0072BD",FaceAlpha=0.7)
        hold on
    area(-x{i},h_lub_model{i},FaceColor="#0072BD",FaceAlpha=0.7)
        
    
    analytical_=plot(x{i},h_analytical{i},'-');
    %analytical_=plot(x_g(1:10:end),h(1:10:end,ii),'r-');
    analytical_.Color="#A2142F";
    analytical_sym=plot(-x{i},h_analytical{i},'-');
    %analytical_sym=plot(-x_g(1:10:end),h(1:10:end,ii),'r-');
    analytical_sym.Color="#A2142F";

    analytical_.LineWidth=3;
    analytical_sym.LineWidth=3;

    
    
    lub_model=plot(x{i}(1:10:end),h_lub_model{i}(1:10:end),marker_shape);
    hold on
    points_top_sym=plot(-x{i}(1:10:end),h_lub_model{i}(1:10:end),marker_shape);
    lub_model.MarkerEdgeColor='k';
    lub_model.MarkerSize=my_marker_size;

    points_top_sym.MarkerEdgeColor='k';
    points_top_sym.MarkerSize=my_marker_size;


    

    xlim([-x_dim*xc x_dim*xc])
    ylim([0 1*h_c]);

    %set(gca,'DataAspectRatio',[9 1 1])

    %plot(x,h_lub_model_an,'.b','DisplayName', name);
    str_lg=strcat('time = ',num2str(time(i)), ' (s)');
     hl=legend(str_lg,Location='northwest',AutoUpdate='off');
    %legend;

    


    M=h_lub_model{i}(j);
    magnification_x=2;
    magnification_y=5;
   
    
    %b_h_init_p=plot(x_init,h_init);
    %b_h_init_p.LineWidth=line_size;

    %x_min_1=x{i}(j)*(1-2.0/(magnification_x*10));
    %x_max_1=x{i}(j)*(1+2.0/(magnification_x*10));

    x_min_1=x_r_lub(i)*(1-2.0/(magnification_x*10));
    x_max_1=x_r_lub(i)*(1+2.0/(magnification_x*10));


    y_min_1=0;
    y_max_1=M*(1.5+0.1/magnification_y);
    
    %r_zoom=rectangle('Position',[x_min_1 y_min_1 (x_max_1-x_min_1) (y_max_1-y_min_1)]);
    
     aro_=quiver(x_min_1, y_min_1,(x_dim*xc*0.67-x_min_1), (0.42*1*h_c-y_min_1),0,'--r');
     aro_.LineWidth=3;
     %aro_.LineSpec="--";


    


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
    f=axes('position',[.72 .42 .17 .38]);
    box on % put box around new pair of axes


    hold on

    points_top_z=plot(x{i},h_lub_model{i},'r.-');
    analytical_z=plot(x{i},h_analytical{i},'b-');

    %analytical_z=plot(x_g(:),h(:,ii),'-r');    


     analytical_z.Color="#A2142F";
     analytical_z.LineWidth=3;
    area(x{i},h_lub_model{i},FaceColor="#0072BD",FaceAlpha=0.7)               

    points_top_z.MarkerEdgeColor='k';
    points_top_z.MarkerSize=12;




    xlim([x_min_1 x_max_1])
    ylim([0 4*b*h_c+1e-9]);

    %hold on;
    grid on;
    %grid(gca,'minor');




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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%num_timesteps=length(time_match);
%%
for i=1:12
    Frame(num_timesteps+i)=Frame(num_timesteps+i-1);
end
%%

writerObj = VideoWriter(strcat('results_vid/evolution_',case_folder,'.avi'),'Uncompressed AVI');
 %writerObj.Quality=100;
 %writerObj.VideoCompressionMethod='H.264';
  writerObj.FrameRate = 12;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:num_timesteps+12
    % convert the image to a frame
    frame = Frame(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);




