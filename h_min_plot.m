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


clc
clear 
folder_path='/home/ISTO/mnorouzisadeh/cluster/Article_models/ratio/b10/';
case_name='curv_g2_fluc';
data_read = readmatrix(strcat(folder_path,case_name,'/output_general.txt'));

%data_read = readmatrix(strcat('output_general.txt'));

%h_min=smooth(data_read(:,2),10,'moving');
h_min=data_read(:,2);

hmin_grad=gradient(h_min,data_read(:,1));
hmin_grad2=gradient(gradient(h_min,data_read(:,1)),data_read(:,1));


ii=1;
ii2=1;
for i=2:length(hmin_grad)-1
    if hmin_grad(i)*hmin_grad(i+1)<0 
        if hmin_grad2(i)>0 %%&& abs(hmin_grad(i))<1
            time_select1(ii)=data_read(i,1);
            h_min_select1(ii)=data_read(i,2);
            ii=ii+1;
        elseif hmin_grad2(i)<0
            time_select2(ii2)=data_read(i,1);
            h_min_select2(ii2)=data_read(i,2);
            ii2=ii2+1;

        end
    end
end


for i=1:length(time_select1)-1
    period_osc(i)=time_select1(i+1)-time_select1(i);
end

for i=1:min(length(h_min_select2),length(h_min_select1))
    ampl_osc(i)=h_min_select2(i)-h_min_select1(i);
end
%%
figure()
plot(data_read(:,1),h_min,'-b','DisplayName',strcat('hmin ',case_name))
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold on
%plot(time_select1,h_min_select1,'.g')

%plot(time_select2,h_min_select2,'.r')

saveas(gcf,strcat('h_min_',case_name),'fig')


figure()
plot(time_select1(1:end-1),period_osc)

figure()
plot(time_select1(1:length(ampl_osc)),ampl_osc)
%%
%g = fittype('a*x^(b)+c');
%fitobject = fit(time_select1',h_min_select1',g);
%figure()
%plot(fitobject,time_select1,h_min_select1)
