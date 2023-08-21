
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



function [output] = Disjoin_der(h1,hc,xc,ham_cte1,w0_h,nb,zeta1,zeta2,kb,T,epsilon,elec_charg,avagad_num,Lambda_2,s,h_grad,h_grad2)



    
    
    psi_1=elec_charg*zeta1/(kb*T);
    psi_2=elec_charg*zeta2/(kb*T);
    
    psi_ave=(psi_1+psi_2)*0.5;
    c=psi_ave^2;
    d=((psi_1-psi_2))^2;
    debye=sqrt(2*elec_charg^2*nb*avagad_num/(epsilon*kb*T));
    

    edl_coef=nb*avagad_num*kb*T*debye;

    vdw_coef=ham_cte1/(2*pi());
    hyd_coef=w0_h/(Lambda_2^2);
    born_length=2e-11;
    d_stern=1e-12;
    born_coef=0;

    if h1*hc<2e-10
        h1=0;
        %fprintf("negative value for h? turned to zero for disjoining  ")
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% solvation Props


    if h1*hc<1e-5 && s(4)==1 && h1~=0

       EDL_CC_der=-edl_coef*(-(c*coth(debye*(h1*hc+d_stern)/2)*csch(debye*(h1*hc+d_stern)/2)^2)/sqrt(c*csch(debye*(h1*hc+d_stern)/2)^2+1) ...
                                                        +d*exp(-debye*(h1*hc+d_stern))/(c*csch(debye*(h1*hc+d_stern)/2)^2+1) ...
                                                       -(c*d*exp(-debye*(h1*hc+d_stern))*coth(debye*(h1*hc+d_stern)/2)*csch(debye*(h1*hc+d_stern)/2)^2)/(c*csch(debye*(h1*hc+d_stern)/2)^2+1)^2);
     
                    
    else
        EDL_CC_der=0;
    end
    
    if h1~=0 && h1*hc<1e-5 && s(5)==1 && h1~=0
        Hydration_p_der=hyd_coef*exp(-(h1)*hc/Lambda_2);
    else 
        Hydration_p_der=0;
    end
    

    if h1~=0 && h1*hc<1e-5 && s(3)==1 
        vdw_p_der=-(vdw_coef)*(1/((h1*hc)^4));

     elseif s(3)==2 && h1~=0
        vdw_p_der=-(vdw_coef/(6*born_length^2))*(6*born_length^2/((h1*hc+born_length)^4)-72*born_coef*born_length^8/((h1*hc+born_length)^10))...
            +vdw_coef*(0.75*(hc/xc)^2*h_grad(1)^2/((h1*hc+born_length)^4)-(hc/(xc^2))*h_grad2(1)/((h1*hc+born_length)^3));


        
        
    else
        vdw_p_der=0;
    end


    output=[vdw_p_der, EDL_CC_der, Hydration_p_der];
        
end
