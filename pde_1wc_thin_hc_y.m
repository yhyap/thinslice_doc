function yt=pde_1wc_thin_hc_y(t,y)

%% Graphical representation of grids
%       u(nu) uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
%       u ^
%       u 2                SOLID CORDIERITE PHASE (Tku)
%       u 1 
%       s(ns) sssssssssssssssssssssssssssssssssssssssssssssssssssssssss
%       s ^
%       s 3                   SOLID WASHCOAT PHASE
%       s 2                        (cas Tks)
%       s 1 --------------------------------------------------- surface
%                                  GAS PHASE
%       ->  r 1 2 3 4 5 6 >         (ca Tk)
%       zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

%% Global area
global pi nr nz ns nu ds du dz dss dus dzs r r0 s u u0 z P Mco Mair Mo2 Mco2 Mc3h8 Vco Vair Vco2 Vc3h8 re por tau Ya Yb Yc Tk Yas Ybs Ycs Tks Tku Yae Ybe Yce Tke fcv fcCO fcO fcOO RCO HisRa HisRb HisRc HisRd HisRe HisRf HisTks HisTke G H Av r_gtc rho_wc rho_cord R vmean ncall History_t

% ignition: 543K, extinction: 513 - 463K
%% Kinetic parameters (sourced from Salomons et al)
% Temperature exponent and reaction order correction are assumed zero
% First assume CO2 does not adsorb back onto active sites
% a) Parameters used for adsorption of CO onto vacant sites (Chatterjee et al & Salomons)
    S0_CO=0.85;      % Sticking coefficient: 0.84
    %AaB=1.8e13;       % Pre-exponential factor in mol cm s (Chatterjee et al)\: 1e13 then 1.8e9
    AaB=1.0e8;%0.1e-6;        % Previous: 0.02e13(1), 0.0072e13(2)
    % Mind desorption term has no G term, at 500K, adsorption >> desorption
% b) Parameters used for adsorption of O2 onto vacant sites (Chatterjee et al & Salomons)
    %S_O2=0.047;%0.045;      % Sticking coefficient - Engel & Ertl said 0.1<x<1.0, default: 0.07   
    S_O2=0.15;
% c) Parameters used for bimolecular surface reaction of O* and CO* (Chatterjee et al & Salomons)
    %AcF=8.156e9;%5.560e9;       % Taken from pg.217 in mol mol-1 s (Salomons's): 2.56e9
    AcF=8.2e8;  
% d) Parameters used for compressed 2O* sites by CO
    AdF=1e14;       % Pre-exponential factor in mol mol-1 s (Salomons pg.229&251):1e12**********
    EdF=-50;        % Activation energy in kJ mol-1 (independent of coverage): -50, -115
% e) Parameters used for bimolecular surface reaction of OO and CO (Salomons only)
    AeF=1e15;%1e15;          % Pre-exponential factor in mol mol-1 s (Salomons pg.251): 1e15
    EeF=-115;       % Activation energy in kJ mol-1 (independent of coverage): -115
% f) Parameters used for surface reaction of OO onto two O sites (Salomons only)
    AfF=1e3;%1e20;       % Pre-exponential factor in mol mol-1 s (Salomons pg.251): 1e15
    EfF=-105;       % Activation energy in kJ mol-1 (independent of coverage)
% Parameters used for propane oxidation
    Ahc=3e12;       % Pre-exponential term for propane oxidation, mol mol-1 s, :1e13
    Ec3h8=-1.08e4;  % including R, activation energy for propane oxidation 
    Ahci=0.208;     % Pre-exponential factor for inhibition term
    Ec3h8i=380;     % From Voltz, scaled accordingly compared with H&K, pls check with Siemund et al
    
%% Other parameters
YO2=0.0933;     % % of O2
dHR=-282.55e3;  % J/mol(CO)
DH=1.20e-3;     % hydraulic diameter (m)= 5.45e-4 x 2, 1.208688mm from SEM
ff=1;

%% 1D to 2D matrices
for i=1:nz
    % Gas phase
        ij=i;
        Ya(i)=y(ij);
        Yb(i)=y(ij+nz);
        Yc(i)=y(ij+2*nz);
        Tk(i)=y(ij+3*nz);
    % Solid washcoat phase
    for jj=1:ns
        ij=4*nz+(i-1)*ns+jj;
        Yas(i,jj)=y(ij);
        Ybs(i,jj)=y(ij+ns*nz);
        Ycs(i,jj)=y(ij+2*ns*nz);
        Tks(i,jj)=y(ij+3*ns*nz);
        fcv(i,jj)=y(ij+4*ns*nz);
        fcCO(i,jj)=y(ij+5*ns*nz);
        fcO(i,jj)=y(ij+6*ns*nz);
        fcOO(i,jj)=y(ij+7*ns*nz);
    end
    % Solid cordierite phase
    for jjj=1:nu
	  ij=4*nz+8*ns*nz+(i-1)*nu+jjj;
	  Tku(i,jjj)=y(ij);
    end
end

%% Step through the grid points in r and z
for i=1:nz
    %% Gas phase
        %% Varying temperature to see light off and hysteresis    
            rt1=200;                % linear ramp time in second
            rt2=500;                % linear ramp time in second
            rt3=300;                % 300 s of constant temperature
            rt4=100;                % linear ramp time in second
            rt5=200;                % linear ramp time in second
            rt6=500;                % linear ramp time in second
            Tkef1=505;              % temperature point 1
            Tkef2=545;              % temperature point 2, max
            Tkef4=485;              % temperature point 3
            Tkef5=445;              % temperature point 4
            mp1=(Tkef1-Tke)/rt1;    % multiplier 1
            mp2=(Tkef2-Tkef1)/rt2;  % multiplier 2
            % mp3 is flat therefore no equation
            mp4=(Tkef4-Tkef2)/rt4;  % multiplier 4
            mp5=(Tkef5-Tkef4)/rt5;  % multiplier 5
            mp6=(Tke-Tkef5)/rt6;    % multiplier 6
            if t<rt1                
                Tkee(i)=Tke+mp1*t;
            elseif t<(rt1+rt2)      
                Tkee(i)=Tkef1+mp2*(t-rt1);
            elseif t<(rt1+rt2+rt3)       
                Tkee(i)=Tkef2;
            elseif t<(rt1+rt2+rt3+rt4)                
                Tkee(i)=Tkef2+mp4*(t-(rt1+rt2+rt3));
            elseif t<(rt1+rt2+rt3+rt4+rt5)
                Tkee(i)=Tkef4+mp5*(t-(rt1+rt2+rt3+rt4));
            elseif t<(rt1+rt2+rt3+rt4+rt5+rt6)
                Tkee(i)=Tkef5+mp6*(t-(rt1+rt2+rt3+rt4+rt5));
            else
                Tkee(i)=Tke;
            end
            
        %% Parameters & variables
            % Calculation of bulk diffusion coefficient and dispersion coefficient (m2/s)
            v(i)=vmean*Tk(i)/298;
            Dca(i)=(1.013e-2*Tk(i)^1.75*((1/Mco)+(1/Mair))^(1/2))/(P*(Vco^(1/3)+Vair^(1/3))^2);
            Dcb(i)=(1.013e-2*Tk(i)^1.75*((1/Mco2)+(1/Mair))^(1/2))/(P*(Vco2^(1/3)+Vair^(1/3))^2);
            Dcc(i)=(1.013e-2*Tk(i)^1.75*((1/Mc3h8)+(1/Mair))^(1/2))/(P*(Vc3h8^(1/3)+Vair^(1/3))^2);
            Dia(i)=Dca(i)+(v(i)*r0)^2/(48*Dca(i));      % Taylor Aris
            Dib(i)=Dcb(i)+(v(i)*r0)^2/(48*Dcb(i));      % Taylor Aris
            Dic(i)=Dcc(i)+(v(i)*r0)^2/(48*Dcc(i));      % Taylor Aris
        
            % Calculation for thermal diffusion coefficient in the fluid (assume all is air)
            k_a(i)=1.679e-2+5.073e-5*Tk(i);             % W/mK (thermal conductivity of air) from H&K
            rho_a(i)=P*Mair/(1000*R*Tk(i));           % kg/m3 (fluid density) - of air at 20C
            Cp_a(i)=(28.09+1.965e-3*Tk(i)+4.799e-6*(Tk(i))^2-1.965e-9*(Tk(i))^3)/.02896; % J/kg K (air specific heat)
            Dt(i)=k_a(i)/(rho_a(i)*Cp_a(i));            % m2/s (thermal diffusivity) (typical:1.5e-5 m2/s)
            Dit(i)=Dt(i)+(v(i)*r0)^2/(48*Dt(i));        % Taylor Aris
    
        %% Reaction terms 
        for jj=1:ns
            %% Mol fraction in the solid phase
            YCO(i,jj)=Yas(i,jj); 
            %% Heat of reaction for propane (from H&K pg. 321)
            dHRhc(i,jj)=(-2.059e6+72.3*Tks(i,jj)-9.69e-2*(Tks(i,jj))^2+4.34e-5*(Tks(i,jj))^3+7.56e-9*(Tks(i,jj))^4);
            %% Self inhibition term for propane reaction
            RHCi(i,jj)=(1+(Ahci*exp(Ec3h8i/Tks(i,jj))*Ycs(i,jj)))^2;
                        
            %% Elementary reaction terms on active sites (follow Salomons et al)
            % UNIT: mol (CO) mol-1 (Pt) s-1
            % Reaction a: CO + * <--> CO*
            EaB(i,jj)=-126.4+33*fcCO(i,jj);       % Activation energy in kJ mol-1 (desorption depends on CO coverage)
            S_CO=S0_CO*((1-fcCO(i,jj))/(1-0.85*fcCO(i,jj)));    % from Engel and Ertl (1979)
            Ra(i,jj)=S_CO/(G*R*Tks(i,jj))*P*YCO(i,jj)*fcv(i,jj)*(R*1000*Tks(i,jj)/(2*pi*Mco))^0.5-AaB*fcCO(i,jj)*exp(EaB(i,jj)*1000/(R*Tks(i,jj)));
            % Reaction b: O2 + 2* -> 2O*
            Rb(i,jj)=2*S_O2/(G*R*Tks(i,jj))*P*YO2*fcv(i,jj)^2*(R*1000*Tks(i,jj)/(2*pi*Mo2))^0.5;
            % Reaction c: O* + CO* -> CO2 + 2*
            %EcF(i,jj)=-108.0+33*fcCO(i,jj);       % Activation energy in kJ mol-1 (depends on CO & NO coverage)
            EcF(i,jj)=-104.6+46.024*fcCO(i,jj);     % from Engel and Ertl pg.47
            Rc(i,jj)=AcF*fcCO(i,jj)*fcO(i,jj)*exp(EcF(i,jj)*1000/(R*Tks(i,jj)));
            % Reaction d: CO + 2O* -> CO* + OO* to CO + O* -> CO* + OO*
            %EdF(i,jj)=-108.0+33*fcCO(i,jj);        % My own modification
            Rd(i,jj)=AdF*P*YCO(i,jj)/(R*Tks(i,jj))*fcO(i,jj)^2*exp(EdF*1000/(R*Tks(i,jj)));
            % Reaction e: CO* + OO* -> CO2 + O* + * (This step is implausible)
            Re(i,jj)=AeF*fcCO(i,jj)*fcOO(i,jj)*exp(EeF*1000/(R*Tks(i,jj)));
            % Reaction f: OO* + * -> 2O* to OO* -> O*
            Rf(i,jj)=AfF*fcOO(i,jj)*fcv(i,jj)*exp(EfF*1000/(R*Tks(i,jj)));
        
            %% Species reaction terms 
            % Rate of disappearance of CO
            RCO(i,jj)=(Ra(i,jj)+Rd(i,jj));
            % Rate of disappearance of O2
            RO2(i,jj)=Rb(i,jj);
            % Rate of appearance of CO2
            RCO2(i,jj)=(Rc(i,jj)+Re(i,jj));
            % Rate of disappearance of propane
            RHC(i,jj)=Ahc*YO2*Ycs(i,jj)*exp(Ec3h8/Tks(i,jj))/RHCi(i,jj);
            
            %% To obtain Damkohler number
            % Overall reaction rate of CO oxidation across the monolith
            %for i=1:nz
            %    RCOjj(i)=(RCO(i,1)+RCO(i,2)+RCO(i,3)+RCO(i,4))*H;
            %end
            %RCOjj=sum(RCO);             % Sum of RCO(i) at each jj
            %RCOtot=sum(RCOjj)*H;        % Total reaction rate
            % I couldn't calculate RCOtot based on each (i) and hence calc the Da number for entire monolith
            
            % Convective mass transport rate
            %cadiff(i)=abs(cas(i,1));                % absolute value
            %Dam(i)=RCOtot*DH/(cadiff(i)*Dca(i));    % Damkohler number
            %Dam(i)=RCO(i,1)*DH/(cas(i,1)*Dca(i));    % Damkohler number

        end

            %% Back to gas phase
            % Dimensionless number (H&K pg.315)
            mu(i)=7.701e-6+4.166e-8*Tk(i)-7.531e-12*(Tk(i))^2;  % Pa.s, from H&K pg.323 & 660
            Rey(i)=rho_a(i)*vmean*DH/mu(i);                     % Local Reynolds Number    
            Pr(i)=Cp_a(i)*mu(i)/k_a(i);%0.7;                    % Prandtl Number for air, pg. 315
            L(i)=i*dz;                                          % Distance from entrance in dirn of gas flow
            Gz(i)=Rey(i)*Pr(i)*DH/L(i);                         % Graetz Number
            NuT(i)=3.657+8.827*(1000/Gz(i))^-0.545*exp(-48.2/Gz(i));    % Nu for constant wall temperature
            NuH(i)=4.364+13.18*(1000/Gz(i))^-0.524*exp(-60.2/Gz(i));    % Nu for constant wall flux
            %Nu(i)=0.5*(NuH(i)-(Dam(i)*NuH(i)/NuT(i))+((NuH(i)-Dam(i)*NuH(i)/NuT(i))^2+4*Dam(i)*NuH(i))^0.5);
            Nu(i)=((NuT(i)+NuH(i))/2)+1;
            % Transfer coefficients and assume mass-heat transfer analogy
            Sh(i)=Nu(i);                        % Sherwood number  
            kma(i)=Sh(i)*Dca(i)/DH;             % MT coefficient for CO
            kmb(i)=Sh(i)*Dcb(i)/DH;             % MT coefficient for CO2
            kmc(i)=Sh(i)*Dcc(i)/DH;             % MT coefficient for propane
            hm(i)=Nu(i)*k_a(i)/DH;              % HT coefficient
                         
        %% Yazz, Ybzz, Yczz, Tkzz, diffusion term in axial direction
        if (i==1)
            Yazz(i)=(Ya(i+1)-2.0*Ya(i)+Yae)/dzs;
            Ybzz(i)=(Yb(i+1)-2.0*Yb(i)+Ybe)/dzs;
            Yczz(i)=(Yc(i+1)-2.0*Yc(i)+Yce)/dzs;
            Tkzz(i)=(Tk(i+1)-2.0*Tk(i)+Tkee(i))/dzs;
        elseif (i==nz)
            Yazz(i)=2.0*(Ya(i-1)-Ya(i))/dzs;
            Ybzz(i)=2.0*(Yb(i-1)-Yb(i))/dzs;
            Yczz(i)=2.0*(Yc(i-1)-Yc(i))/dzs;
            Tkzz(i)=2.0*(Tk(i-1)-Tk(i))/dzs;
        else
            Yazz(i)=(Ya(i+1)-2.0*Ya(i)+Ya(i-1))/dzs;
            Ybzz(i)=(Yb(i+1)-2.0*Yb(i)+Yb(i-1))/dzs;
            Yczz(i)=(Yc(i+1)-2.0*Yc(i)+Yc(i-1))/dzs;
            Tkzz(i)=(Tk(i+1)-2.0*Tk(i)+Tk(i-1))/dzs;
        end
        
        %% Yaz, Ybz, Ycz, Tkz, advection term in axial direction
        if (i==1)
            Yaz(i)=(Ya(i)-Yae)/dz;
            Ybz(i)=(Yb(i)-Ybe)/dz;
            Ycz(i)=(Yc(i)-Yce)/dz;
            Tkz(i)=(Tk(i)-Tkee(i))/dz;
        else
            Yaz(i)=(Ya(i)-Ya(i-1))/dz;
            Ybz(i)=(Yb(i)-Yb(i-1))/dz;
            Ycz(i)=(Yc(i)-Yc(i-1))/dz;
            Tkz(i)=(Tk(i)-Tk(i-1))/dz;
        end
        
        %% PDEs
        Yat(i)=Dia(i)*(Yazz(i))-v(i)*Yaz(i)-4/DH*kma(i)*(Ya(i)-Yas(i,1));
        Ybt(i)=Dib(i)*(Ybzz(i))-v(i)*Ybz(i)-4/DH*kmb(i)*(Yb(i)-Ybs(i,1));
        Yct(i)=Dic(i)*(Yczz(i))-v(i)*Ycz(i)-4/DH*kmc(i)*(Yc(i)-Ycs(i,1));
        Tkt(i)=Dit(i)*(Tkzz(i))-v(i)*Tkz(i)+4/DH*hm(i)/(rho_a(i)*Cp_a(i))*(Tks(i,1)-Tk(i));  
    
    
    %% Solid washcoat phase
    for jj=1:ns
        %% Parameters & variables
            % Calculation for effective diffusivity of washcoat
            Dca(i,jj)=(1.013e-2*Tks(i,jj)^1.75*((1/Mco)+(1/Mair))^(1/2))/(P*(Vco^(1/3)+Vair^(1/3))^2);
            Dcb(i,jj)=(1.013e-2*Tks(i,jj)^1.75*((1/Mco2)+(1/Mair))^(1/2))/(P*(Vco2^(1/3)+Vair^(1/3))^2);
            Dcc(i,jj)=(1.013e-2*Tks(i,jj)^1.75*((1/Mc3h8)+(1/Mair))^(1/2))/(P*(Vc3h8^(1/3)+Vair^(1/3))^2);
            Dka(i,jj)=97*re*(Tks(i,jj)/Mco)^0.5;            % Knudsen diffusion coefficient for CO
            Dkb(i,jj)=97*re*(Tks(i,jj)/Mco2)^0.5;           % Knudsen diffusion coefficient for CO2
            Dkc(i,jj)=97*re*(Tks(i,jj)/Mc3h8)^0.5;          % Knudsen diffusion coefficient for propane
            Dva(i,jj)=1/((1/Dca(i,jj))+(1/Dka(i,jj)));      % Diffusion coefficient in the pore for CO
            Dvb(i,jj)=1/((1/Dcb(i,jj))+(1/Dkb(i,jj)));      % Diffusion coefficient in the pore for CO2
            Dvc(i,jj)=1/((1/Dcc(i,jj))+(1/Dkc(i,jj)));      % Diffusion coefficient in the pore for propane
            Dea(i,jj)=ff*por*Dva(i,jj)/tau;                    % m2/s (effective diffusivity for CO)
            Deb(i,jj)=ff*por*Dvb(i,jj)/tau;                    % m2/s (CO2)           
            Dec(i,jj)=ff*por*Dvc(i,jj)/tau;                    % m2/s (propane)
            
            % Calculation for thermal effective diffusivity in the washcoat
            k_s(i,jj)=0.9558-2.09e-4*Tks(i,jj);           % W/mK (solid thermal conductivity)
            rho_s(i,jj)=rho_wc/1000;                      % kg/m3 (washcoat density)
            Cp_s(i,jj)=948+0.2268*Tks(i,jj);              % J/kg K (solid specific heat) from H&K wall
            Ds(i,jj)=k_s(i,jj)/(rho_s(i,jj)*Cp_s(i,jj));  % m2/s (thermal diffusivity) (typical:1.6e-6 m2/s)
        
        %% (1/s)*Yass, (1/s)*Ybss, (1/s)*Ycss, (1/s)*Tkss, first derivative groups in radial direction
        if (jj==1)  % Is this correct, similar to above 
            % ***Does l'Hopital's rule applicable here since s=0 when jj=1
            % Equations below account for flux from fluid to surface to solid
            % Front: Borrows the idea from Yasss when jj=ns, Back: H&K pg.339 & 340, need to add to balance unit
            % dw, eff. wall thickness is c.s.a./wetted perimeter, H&K pg.308
            % DH or ds or dr for denominator of MT?
            Yass(i,jj)=1.0*(1/por)*((Dea(i,jj)*(Yas(i,jj+1)-Yas(i,jj))/dss+kma(i)*(Ya(i)-Yas(i,jj))/ds)); 
            Ybss(i,jj)=1.0*(1/por)*((Deb(i,jj)*(Ybs(i,jj+1)-Ybs(i,jj))/dss+kmb(i)*(Yb(i)-Ybs(i,jj))/ds)); 
            Ycss(i,jj)=1.0*(1/por)*((Dec(i,jj)*(Ycs(i,jj+1)-Ycs(i,jj))/dss+kmc(i)*(Yc(i)-Ycs(i,jj))/ds)); 
            Tkss(i,jj)=1.0*(Ds(i,jj)*(Tks(i,jj+1)-Tks(i,jj))/dss)-(hm(i)*(Tks(i,jj)-Tk(i))/(ds*rho_s(i,jj)*Cp_s(i,jj)));
        elseif (jj==ns) % Boundary
            Yass(i,jj)=0;       % As dc/ds=0
            Ybss(i,jj)=0;  
            Ycss(i,jj)=0;
            %Tkss(i,jj)=0;
            Tkss(i,jj)=(1.0/s(jj))*(Tku(i,1)-Tks(i,jj-1))/(2.0*ds);
            %Tkss(i,jj)=(1.0/s(jj))*((Tku(i,1)-Tks(i,jj))/((du+ds)/2));
        else
            Yass(i,jj)=(1.0/(s(jj)))*(Yas(i,jj+1)-Yas(i,jj-1))/(2.0*ds);
            Ybss(i,jj)=(1.0/(s(jj)))*(Ybs(i,jj+1)-Ybs(i,jj-1))/(2.0*ds);
            Ycss(i,jj)=(1.0/(s(jj)))*(Ycs(i,jj+1)-Ycs(i,jj-1))/(2.0*ds);
            Tkss(i,jj)=(1.0/(s(jj)))*(Tks(i,jj+1)-Tks(i,jj-1))/(2.0*ds);
        end
        
        %% Yasss, Ybsss, Ycsss, Tksss, second derivative diffusion term in radial direction
        if (jj==1)
            Yasss(i,jj)=0;%((cas(i,jj+2)-2.0*cas(i,jj+1)+cas(i,jj))*drs/dss+Dc/De*(ca(i,nr-1)-2.0*ca(i,nr)+cas(i,jj))*dss/drs)/2;
            Ybsss(i,jj)=0;
            Ycsss(i,jj)=0;
            Tksss(i,jj)=0;%((Tks(i,jj+2)-2.0*Tks(i,jj+1)+Tks(i,jj))*drs/dss+Dt/Du*(Tk(i,nr-1)-2.0*Tk(i,nr)+Tks(i,jj))*dss/drs)/2;
        elseif (jj==ns) % Boundary
            Yasss(i,jj)=2.0*(Yas(i,jj-1)-Yas(i,jj))/dss;
            Ybsss(i,jj)=2.0*(Ybs(i,jj-1)-Ybs(i,jj))/dss;
            Ycsss(i,jj)=2.0*(Ycs(i,jj-1)-Ycs(i,jj))/dss;
            %Tksss(i,jj)=2.0*(Tks(i,jj-1)-Tks(i,jj))/dss;
            Tksss(i,jj)=(Tku(i,1)-2.0*Tks(i,jj)+Tks(i,jj-1))/dss;
        else
            Yasss(i,jj)=(Yas(i,jj+1)-2.0*Yas(i,jj)+Yas(i,jj-1))/dss;
            Ybsss(i,jj)=(Ybs(i,jj+1)-2.0*Ybs(i,jj)+Ybs(i,jj-1))/dss;
            Ycsss(i,jj)=(Ycs(i,jj+1)-2.0*Ycs(i,jj)+Ycs(i,jj-1))/dss;
            Tksss(i,jj)=(Tks(i,jj+1)-2.0*Tks(i,jj)+Tks(i,jj-1))/dss;
        end
        
        %% Yaszz, Ybszz, Ycszz, Tkszz, diffusion term in axial direction
        if (i==1)
            Yaszz(i,jj)=2.0*(Yas(i+1,jj)-Yas(i,jj))/dzs;
            Ybszz(i,jj)=2.0*(Ybs(i+1,jj)-Ybs(i,jj))/dzs;
            Ycszz(i,jj)=2.0*(Ycs(i+1,jj)-Ycs(i,jj))/dzs;
            Tkszz(i,jj)=2.0*(Tks(i+1,jj)-Tks(i,jj))/dzs;
        elseif (i==nz)
            Yaszz(i,jj)=2.0*(Yas(i-1,jj)-Yas(i,jj))/dzs;
            Ybszz(i,jj)=2.0*(Ybs(i-1,jj)-Ybs(i,jj))/dzs;
            Ycszz(i,jj)=2.0*(Ycs(i-1,jj)-Ycs(i,jj))/dzs;
            Tkszz(i,jj)=2.0*(Tks(i-1,jj)-Tks(i,jj))/dzs;
        else
            Yaszz(i,jj)=(Yas(i+1,jj)-2.0*Yas(i,jj)+Yas(i-1,jj))/dzs;
            Ybszz(i,jj)=(Ybs(i+1,jj)-2.0*Ybs(i,jj)+Ybs(i-1,jj))/dzs;
            Ycszz(i,jj)=(Ycs(i+1,jj)-2.0*Ycs(i,jj)+Ycs(i-1,jj))/dzs;
            Tkszz(i,jj)=(Tks(i+1,jj)-2.0*Tks(i,jj)+Tks(i-1,jj))/dzs;
        end
        
                  
        %% PDEs
        if (jj==1)
            Yast(i,jj)=Yass(i,jj)+(1/por)*(Dea(i,jj)*Yaszz(i,jj)-RCO(i,jj)*H*Av*R*Tks(i,jj)/P);%/r_gtc);
            Ybst(i,jj)=Ybss(i,jj)+(1/por)*(Deb(i,jj)*Ybszz(i,jj)+RCO2(i,jj)*H*Av*R*Tks(i,jj)/P);%/r_gtc);
            Ycst(i,jj)=Ycss(i,jj)+(1/por)*(Dec(i,jj)*Ycszz(i,jj)-RHC(i,jj)*H*Av*R*Tks(i,jj)/P);%/r_gtc);
            Tkst(i,jj)=Tkss(i,jj)+Ds(i,jj)*Tkszz(i,jj)+(-RCO(i,jj)*dHR-RHC(i,jj)*dHRhc(i,jj))*H*Av/(rho_s(i,jj)*Cp_s(i,jj));
            fcvt(i,jj)=-Ra(i,jj)-2*Rb(i,jj)+2*Rc(i,jj)+Re(i,jj)-Rf(i,jj);
            fcCOt(i,jj)=Ra(i,jj)-Rc(i,jj)+Rd(i,jj)-Re(i,jj);
            fcOt(i,jj)=2*Rb(i,jj)-Rc(i,jj)-2*Rd(i,jj)+Re(i,jj)+2*Rf(i,jj);
            fcOOt(i,jj)=Rd(i,jj)-Re(i,jj)-Rf(i,jj);
        else
            Yast(i,jj)=(1/por)*(Dea(i,jj)*(Yass(i,jj)+Yasss(i,jj)+Yaszz(i,jj))-RCO(i,jj)*H*Av*R*Tks(i,jj)/P);%/r_gtc);
            Ybst(i,jj)=(1/por)*(Deb(i,jj)*(Ybss(i,jj)+Ybsss(i,jj)+Ybszz(i,jj))+RCO2(i,jj)*H*Av*R*Tks(i,jj)/P);%/r_gtc);
            Ycst(i,jj)=(1/por)*(Dec(i,jj)*(Ycss(i,jj)+Ycsss(i,jj)+Ycszz(i,jj))-RHC(i,jj)*H*Av*R*Tks(i,jj)/P);%/r_gtc);
            Tkst(i,jj)=Ds(i,jj)*(Tkss(i,jj)+Tksss(i,jj)+Tkszz(i,jj))+(-RCO(i,jj)*dHR-RHC(i,jj)*dHRhc(i,jj))*H*Av/(rho_s(i,jj)*Cp_s(i,jj));
            fcvt(i,jj)=-Ra(i,jj)-2*Rb(i,jj)+2*Rc(i,jj)+Re(i,jj)-Rf(i,jj);
            fcCOt(i,jj)=Ra(i,jj)-Rc(i,jj)+Rd(i,jj)-Re(i,jj);
            fcOt(i,jj)=2*Rb(i,jj)-Rc(i,jj)-2*Rd(i,jj)+Re(i,jj)+2*Rf(i,jj);
            fcOOt(i,jj)=Rd(i,jj)-Re(i,jj)-Rf(i,jj);
        end

    end

    %% Solid cordierite phase
    for jjj=1:nu
	    %% Parameters & variables  
            % Calculation for thermal effective diffusivity in the washcoat
            k_u(i,jjj)=5.423e-9*Tku(i,jjj)^3-1.276e-5*Tku(i,jjj)^2+8.107e-3*Tku(i,jjj)+7.975e-1; % W/mK Solid thermal conductivity (H&K pg. 657 from Corning)
            rho_u(i,jjj)=rho_cord/1000;                     % kg/m3 (density of substrate)
            Cp_u(i,jjj)=1071+0.156*Tku(i,jjj)-3.437e7*Tku(i,jjj)^-2; % J/kg K (solid specific heat) from H&K pg.655
            Du(i,jjj)=k_u(i,jjj)/(rho_u(i,jjj)*Cp_u(i,jjj));% m2/s (thermal diffusivity) (typical:1.6e-6 m2/s)
    
    %% (1/u)*Tkuu, first derivative groups in radial direction
        if (jjj==1)  % Is this correct, similar to above 
            % Equations below account for flux from fluid to surface to solid
            % Front: Borrows the idea from Yasss when jj=ns, Back: H&K pg.339 & 340, need to add to balance unit
            % dw, eff. wall thickness is c.s.a./wetted perimeter, H&K pg.308
            % DH or ds or dr for denominator of MT?
            Tkuu(i,jjj)=0;%(2.0*Du(i,jjj)*(Tku(i,jjj+1)-Tku(i,jjj))/dus)-(2.0*(Ds(i,ns)*(Tku(i,jjj)-Tks(i,ns))/((dus+dss)/2)));
        elseif (jjj==nu) % Boundary
            Tkuu(i,jjj)=0;
        else
            Tkuu(i,jjj)=(1.0/(u(jjj)))*(Tku(i,jjj+1)-Tku(i,jjj-1))/(2.0*du);
        end

    %% Tkuuu, second derivative thermal diffusion term in radial direction
        if (jjj==1)
            %Tkuuu(i,jjj)=0;
            Tkuuu(i,jjj)=1*(Du(i,jjj)*(Tku(i,jjj+1)-Tku(i,jjj))/dus)-((Ds(i,ns)*(Tku(i,jjj)-Tks(i,ns))/dss));
        elseif (jjj==nu) % Boundary
            Tkuuu(i,jjj)=2.0*(Tku(i,jjj-1)-Tku(i,jjj))/dus;
        else
            Tkuuu(i,jjj)=(Tku(i,jjj+1)-2.0*Tku(i,jjj)+Tku(i,jjj-1))/dus;
        end
 
    %% Tkuzz, diffusion term in axial direction
        if (i==1)
            Tkuzz(i,jjj)=2.0*(Tku(i+1,jjj)-Tku(i,jjj))/dzs;
        elseif (i==nz)
            Tkuzz(i,jjj)=2.0*(Tku(i-1,jjj)-Tku(i,jjj))/dzs;
        else
            Tkuzz(i,jjj)=(Tku(i+1,jjj)-2.0*Tku(i,jjj)+Tku(i-1,jjj))/dzs;
        end
 
        %% PDEs
        if (jjj==1)
            Tkut(i,jjj)=Tkuu(i,jjj)+Tkuuu(i,jjj)+Du(i,jjj)*Tkuzz(i,jjj); 
        else
            Tkut(i,jjj)=Du(i,jjj)*(Tkuu(i,jjj)+Tkuuu(i,jjj)+Tkuzz(i,jjj)); 
        end

     end

end

            %% Record of rate and temperature for hysteresis study
            for i=1:nz
                for jj=1:ns
            if (t==0)
                History_t=t;
                HisTks=Tks(2,2);
                HisTke=Tkee(1);
                HisRa=Ra(2,2);
                HisRb=Rb(2,2);
                HisRc=Rc(2,2);
                HisRd=Rd(2,2);
                HisRe=Re(2,2);
                HisRf=Rf(2,2);
            else
                History_t=[History_t;t];
                HisTks=[HisTks;Tks(2,2)];
                HisTke=[HisTke;Tkee(1)];
                HisRa=[HisRa;Ra(2,2)];
                HisRb=[HisRb;Rb(2,2)];
                HisRc=[HisRc;Rc(2,2)];
                HisRd=[HisRd;Rd(2,2)];
                HisRe=[HisRe;Re(2,2)];
                HisRf=[HisRf;Rf(2,2)];
            end
                end
            end

%% 2D to 1D matrices
for i=1:nz
    % Gas phase
        ij=i;
        yt(ij)=Yat(i);
        yt(ij+nz)=Ybt(i);
        yt(ij+2*nz)=Yct(i);
        yt(ij+3*nz)=Tkt(i);
    % Solid washcoat phase
    for jj=1:ns
        ij=4*nz+(i-1)*ns+jj;
        yt(ij)=Yast(i,jj);
        yt(ij+ns*nz)=Ybst(i,jj);
        yt(ij+2*ns*nz)=Ycst(i,jj);
        yt(ij+3*ns*nz)=Tkst(i,jj);
        yt(ij+4*ns*nz)=fcvt(i,jj);
        yt(ij+5*ns*nz)=fcCOt(i,jj);
        yt(ij+6*ns*nz)=fcOt(i,jj);
        yt(ij+7*ns*nz)=fcOOt(i,jj);
    end
    % Solid cordierite phase
    for jjj=1:nu
        ij=4*nz+8*ns*nz+(i-1)*nu+jjj;
        yt(ij)=Tkut(i,jjj);
    end
end

%% Transpose and steps count
yt=yt';
ncall=ncall+1;
fprintf('\n t = %5d\n',t);
