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

function output=roughness_gen(num_p,alpha,roughness_multiplyer)

    nx_rough=2.^nextpow2(num_p);
    eta=randn(2*nx_rough+1);
    for i=nx_rough/2:nx_rough+nx_rough/2
        %for j=ny/2:ny+ny/2
            sum=0.0;
            for k=-nx_rough/2+1:nx_rough/2-1
                %for l=-ny/2+1:ny/2-1
                    if (k==0 )
                        sum=sum+alpha^2*pi()*eta(i+k);
                    else
                        
                        sum=sum+alpha*besselj(1,2*pi()*alpha*sqrt(k^2))/sqrt(k^2)*eta(i+k);
                    end
                    
                    
                %end
            end
            z(i)=sum;
        %end
    end
    
    h_rough=z(nx_rough/2:nx_rough+nx_rough/2);
    if min(h_rough)<0
       h_rough=h_rough- min(h_rough);
    end
    h_rough=h_rough*roughness_multiplyer;
    output=h_rough(1:num_p);

end
