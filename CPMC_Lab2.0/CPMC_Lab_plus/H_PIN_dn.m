function H=H_PIN_dn(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz)
% function H=H_K(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz) 只有一个粒子的时候的哈密顿量动能项
% Generate the one-body kinetic term of the Hubbard Hamiltonian with the given parameters
% Input:
%   Lx: The number of lattice sites in the x direction.
%   Ly: The number of lattice sites in the y direction.
%   Lz: The number of lattice sites in the z direction.
%   kx: The x component of the twist angle in TABC (twist-averaging boundary conditions)
%   ky: The y component of the twist angle in TABC
%   kz: The z component of the twist angle in TABC
%   tx: The hopping amplitude between nearest-neighbor sites in the x direction
%   ty: The hopping amplitude between nearest neighbor sites in the y direction
%   tz: The hopping amplitude between nearest neighbor sites in the y direction
% Output
%   H: The one-body kinetic Hamiltonian in the form of a square matrix of size (Lx*Ly*Lz) 
% 
% Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
% ?014 v1.0
% Package homepage: http://cpmc-lab.wm.edu
% Distributed under the <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">Computer Physics Communications Non-Profit Use License</a>
% Any publications resulting from either applying or building on the present package 
%   should cite the following journal article (in addition to the relevant literature on the method):
% "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" Comput. Phys. Commun. (2014)

    r=0;
    N_sites=Lx*Ly*Lz;
    kx=sqrt(-1)*pi*kx;
    ky=sqrt(-1)*pi*ky;
    kz=sqrt(-1)*pi*kz;
    H=zeros(N_sites,N_sites);

    for mz=1:Lz
        for iy=1:Ly
            for jx=1:Lx
                r=r+1;      % r=(iy-1)*Lx+jx;
                if Lx~=1
                    if jx==1
                        H(r,r)=H(r,r)+(-1)^Ly*0;
                    elseif jx==Lx
                        H(r,r)=H(r,r)+(-1)^Ly*0;
                    end
                end
            end
        end
    end

end