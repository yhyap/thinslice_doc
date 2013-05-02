function yt=pde_1wc_thin_hc_y2D(t,y)

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
global pi nr nz ns nu dr ds du dz drs dss dus dzs r r0 s u u0 z P Mco Mair Mo2 Mco2 Mc3h8 Vco Vair Vco2 Vc3h8 re por tau Ya Yb Yc Tk Yas Ybs Ycs Tks Tku Yae Ybe Yce Tke fcv fcCO fcO fcOO RCO HisRa HisRb HisRc HisRd HisRe HisRf HisTks HisTke G H Av r_gtc rho_wc rho_cord R vmean ncall History_t

% ignition: 543K, extinction: 513 - 463K
%% Kinetic parameters (sourced from Salomons et al)
% Temperature exponent and reaction order correction are assumed zero
% First assume CO2 does not adsorb back onto active sites
% a) Parameters used for adsorption of CO onto vacant sites (Chatterjee et al & Salomons)
    S0_CO=0.85;      % Sticking coefficient: 0.84
    %AaB=1.8e13;       % Pre-exponential factor in mol cm s (Chatterjee et al)\: 1e13 then 1.8e9
    AaB=1e8;        % Previous: 0.02e13(1), 0.0072e13(2)
    % Mind desorption term has no G term, at 500K, adsorption >> desorption
% b) Parameters used for adsorption of O2 onto vacant sites (Chatterjee et al & Salomons)
    S_O2=0.15;      % Sticking coefficient - Engel & Ertl said 0.1<x<1.0, default: 0.07
    
% c) Parameters used for bimolecular surface reaction of O* and CO* (Chatterjee et al & Salomons)
    AcF=8.2e8;       % Taken from pg.217 in mol mol-1 s (Salomons's): 1.006e12, previous was 3.2e11** 
                % 2.5e12
% d) Parameters used for compressed 2O* sites by CO
    AdF=0;        % Pre-exponential factor in mol mol-1 s (Salomons pg.229&251):1e12**********
    EdF=-50;        % Activation energy in kJ mol-1 (independent of coverage): -50, -115
% e) Parameters used for bimolecular surface reaction of OO and CO (Salomons only)
    AeF=0;          % Pre-exponential factor in mol mol-1 s (Salomons pg.251): 1e15
    EeF=-115;       % Activation energy in kJ mol-1 (independent of coverage): -115
% f) Parameters used for surface reaction of OO onto two O sites (Salomons only)
    AfF=1e15;       % Pre-exponential factor in mol mol-1 s (Salomons pg.251): 1e15
    EfF=-105;       % Activation energy in kJ mol-1 (independent of coverage)
% Parameters used for propane oxidation
    Ahc=3e12;       % Pre-exponential term for propane oxidation, mol mol-1 s, :1e13
    Ec3h8=-1.08e4;  % including R, activation energy for propane oxidation 
    Ahci=0.208;     % Pre-exponential factor for inhibition term
    Ec3h8i=380;     % From Voltz, scaled accordingly compared with H&K, pls check with Siemund et al
    
%% Other parameters
YO2=0.0933;     % % of O2
dHR=-282.55e3;  % J/mol(CO)
ff=1;
gg=1;

%% 1D to 2D matrices
for i=1:nz
    % Gas phase
    for j=1:nr
        ij=(i-1)*nr+j;
        Ya(i,j)=y(ij);
        Yb(i,j)=y(ij+nz*nr);
        Yc(i,j)=y(ij+2*nz*nr);
        Tk(i,j)=y(ij+3*nz*nr);
    end
    % Solid washcoat phase
    for jj=1:ns
        ij=4*nz*nr+(i-1)*ns+jj;
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
	  ij=4*nz*nr+8*ns*nz+(i-1)*nu+jjj;
	  Tku(i,jjj)=y(ij);
    end
end

%% Step through the grid points in r and z
for i=1:nz
    %% Gas phase
    for j=1:nr
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
                Tkee(1,j)=Tke+mp1*t;
            elseif t<(rt1+rt2)      
                Tkee(1,j)=Tkef1+mp2*(t-rt1);
            elseif t<(rt1+rt2+rt3)       
                Tkee(1,j)=Tkef2;
            elseif t<(rt1+rt2+rt3+rt4)                
                Tkee(1,j)=Tkef2+mp4*(t-(rt1+rt2+rt3));
            elseif t<(rt1+rt2+rt3+rt4+rt5)
                Tkee(1,j)=Tkef4+mp5*(t-(rt1+rt2+rt3+rt4));
            elseif t<(rt1+rt2+rt3+rt4+rt5+rt6)
                Tkee(1,j)=Tkef5+mp6*(t-(rt1+rt2+rt3+rt4+rt5));
            else
                Tkee(1,j)=Tke;
            end
            
        %% Parameters & variables
            % Calculation of bulk diffusion coefficient and dispersion coefficient (m2/s)
            v(i,j)=2*vmean*(1-((r(j))^2/r0^2))*Tk(i,j)/298;
            Dca(i,j)=gg*(1.013e-2*Tk(i,j)^1.75*((1/Mco)+(1/Mair))^(1/2))/(P*(Vco^(1/3)+Vair^(1/3))^2);
            Dcb(i,j)=gg*(1.013e-2*Tk(i,j)^1.75*((1/Mco2)+(1/Mair))^(1/2))/(P*(Vco2^(1/3)+Vair^(1/3))^2);
            Dcc(i,j)=gg*(1.013e-2*Tk(i,j)^1.75*((1/Mc3h8)+(1/Mair))^(1/2))/(P*(Vc3h8^(1/3)+Vair^(1/3))^2);
        
            % Calculation for thermal diffusion coefficient in the fluid (assume all is air)
            k_a(i,j)=1.679e-2+5.073e-5*Tk(i,j);                 % W/mK (thermal conductivity) from H&K
            rho_a(i,j)=P*Mair/(1000*R*Tk(i,j));                 % kg/m3 (fluid density) - of air at 20C
            Cp_a(i,j)=(28.09+1.965e-3*Tk(i,j)+4.799e-6*(Tk(i,j))^2-1.965e-9*(Tk(i,j))^3)/.02896; % J/kg K (air specific heat)
            Dt(i,j)=gg*k_a(i,j)/(rho_a(i,j)*Cp_a(i,j));            % m2/s (thermal diffusivity) (typical:1.5e-5 m2/s)
    
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
        end
                         
        %% (1/r)*Yar, (1/r)*Ybr, (1/r)*Ycr, (1/r)*Tkr, first derivative groups in radial direction
        if (j==1)   % l'Hopital's rule applies because r=0 and dc/dr=0
            Yar(i,j)=2.0*(Ya(i,j+1)-Ya(i,j))/drs;
            Ybr(i,j)=2.0*(Yb(i,j+1)-Yb(i,j))/drs;
            Ycr(i,j)=2.0*(Yc(i,j+1)-Yc(i,j))/drs;
            Tkr(i,j)=2.0*(Tk(i,j+1)-Tk(i,j))/drs;
        elseif (j==nr)  % Last point in gas phase before surface
            % There is no boundary setting and should be treated as continuum
            % Therefore the below equations are correct
            % Should have Dea here instead in the solid phase, above
            % assumed bulk diff coeff to the first layer
            % but above only apply when dc/dr=0, so...
            %Yar(i,j)=(1.0/r(j))*(Yas(i,1)-Ya(i,j))/((dr+ds)/2);
            %Ybr(i,j)=(1.0/r(j))*(Ybs(i,1)-Yb(i,j))/((dr+ds)/2);
            %Ycr(i,j)=(1.0/r(j))*(Ycs(i,1)-Yc(i,j))/((dr+ds)/2);
            %Tkr(i,j)=(1.0/r(j))*(Tks(i,1)-Tk(i,j))/((dr+ds)/2);
            Yar(i,j)=(1.0/r(j))*(Yas(i,1)-Ya(i,j-1))/(2.0*dr);
            Ybr(i,j)=(1.0/r(j))*(Ybs(i,1)-Yb(i,j-1))/(2.0*dr);
            Ycr(i,j)=(1.0/r(j))*(Ycs(i,1)-Yc(i,j-1))/(2.0*dr);
            Tkr(i,j)=(1.0/r(j))*(Tks(i,1)-Tk(i,j-1))/(2.0*dr);
            % (dr+ds)/2 can be approximated by dr?
            %Yar(i,j)=(1.0/r(j))*(Yas(i,1)-Ya(i,j-1))/(2.0*dr);%-(Ya(i,j)-Ya(i,j-1))/dr);
            %Ybr(i,j)=(1.0/r(j))*(Ybs(i,1)-Yb(i,j-1))/(2.0*dr);%-(Yb(i,j)-Yb(i,j-1))/dr);
            %Ycr(i,j)=(1.0/r(j))*(Ycs(i,1)-Yc(i,j-1))/(2.0*dr);%-(Yc(i,j)-Yc(i,j-1))/dr);
            %Tkr(i,j)=(1.0/r(j))*(Tks(i,1)-Tk(i,j-1))/(2.0*dr);%-(Tk(i,j)-Tk(i,j-1))/dr);
            %Yar(i,j)=((Yas(i,1)-Ya(i,j))/((drs+dss)/2))-((Ya(i,j)-Ya(i,j-1))/drs);
            %Ybr(i,j)=((Ybs(i,1)-Yb(i,j))/((drs+dss)/2))-((Yb(i,j)-Yb(i,j-1))/drs);
            %Ycr(i,j)=((Ycs(i,1)-Yc(i,j))/((drs+dss)/2))-((Yc(i,j)-Yc(i,j-1))/drs);
            %Tkr(i,j)=((Tks(i,1)-Tk(i,j))/((drs+dss)/2))-((Tk(i,j)-Tk(i,j-1))/drs); % This needs change??
        else
            Yar(i,j)=(1.0/r(j))*(Ya(i,j+1)-Ya(i,j-1))/(2.0*dr);
            Ybr(i,j)=(1.0/r(j))*(Yb(i,j+1)-Yb(i,j-1))/(2.0*dr);
            Ycr(i,j)=(1.0/r(j))*(Yc(i,j+1)-Yc(i,j-1))/(2.0*dr);
            Tkr(i,j)=(1.0/r(j))*(Tk(i,j+1)-Tk(i,j-1))/(2.0*dr);
        end
        
        %% Yarr, Ybrr, Ycrr, Tkrr, second derivative groups in radial direction
        if (j==1)
            Yarr(i,j)=2.0*(Ya(i,j+1)-Ya(i,j))/drs;
            Ybrr(i,j)=2.0*(Yb(i,j+1)-Yb(i,j))/drs;
            Ycrr(i,j)=2.0*(Yc(i,j+1)-Yc(i,j))/drs;
            Tkrr(i,j)=2.0*(Tk(i,j+1)-Tk(i,j))/drs;
        elseif (j==nr)  % treat as continuum, no MTC or HTC
            %Yarr(i,j)=0;%2.0*(ca(i,j-1)-ca(i,j))/drs; %(cas(i,1)-2.0*ca(i,j)+ca(i,j-1))/drs; %
            %Ybrr(i,j)=0;
            %Ycrr(i,j)=0;
            %Tkrr(i,j)=0;%(Tks(i,1)-2.0*Tk(i,j)+Tk(i,j-1))/drs;
            % approximate the distance below
            Yarr(i,j)=(Yas(i,1)-2.0*Ya(i,j)+Ya(i,j-1))/drs;
            Ybrr(i,j)=(Ybs(i,1)-2.0*Yb(i,j)+Yb(i,j-1))/drs;
            Ycrr(i,j)=(Ycs(i,1)-2.0*Yc(i,j)+Yc(i,j-1))/drs;
            Tkrr(i,j)=(Tks(i,1)-2.0*Tk(i,j)+Tk(i,j-1))/drs;
            %Yarr(i,j)=((Yas(i,1)-Ya(i,j))/((drs+dss)/2))-((Ya(i,j)-Ya(i,j-1))/drs);
            %Ybrr(i,j)=((Ybs(i,1)-Yb(i,j))/((drs+dss)/2))-((Yb(i,j)-Yb(i,j-1))/drs);
            %Ycrr(i,j)=((Ycs(i,1)-Yc(i,j))/((drs+dss)/2))-((Yc(i,j)-Yc(i,j-1))/drs);
            %Tkrr(i,j)=((Tks(i,1)-Tk(i,j))/((drs+dss)/2))-((Tk(i,j)-Tk(i,j-1))/drs);
        else
            Yarr(i,j)=(Ya(i,j+1)-2.0*Ya(i,j)+Ya(i,j-1))/drs;
            Ybrr(i,j)=(Yb(i,j+1)-2.0*Yb(i,j)+Yb(i,j-1))/drs;
            Ycrr(i,j)=(Yc(i,j+1)-2.0*Yc(i,j)+Yc(i,j-1))/drs;
            Tkrr(i,j)=(Tk(i,j+1)-2.0*Tk(i,j)+Tk(i,j-1))/drs;
        end
            
        %% Yazz, Ybzz, Tkzz, diffusion term in axial direction
        if (i==1)
            Yazz(i,j)=(Ya(i+1,j)-2.0*Ya(i,j)+Yae)/dzs;
            Ybzz(i,j)=(Yb(i+1,j)-2.0*Yb(i,j)+Ybe)/dzs;
            Yczz(i,j)=(Yc(i+1,j)-2.0*Yc(i,j)+Yce)/dzs;
            Tkzz(i,j)=(Tk(i+1,j)-2.0*Tk(i,j)+Tkee(1,j))/dzs;
        elseif (i==nz)
            Yazz(i,j)=2.0*(Ya(i-1,j)-Ya(i,j))/dzs;
            Ybzz(i,j)=2.0*(Yb(i-1,j)-Yb(i,j))/dzs;
            Yczz(i,j)=2.0*(Yc(i-1,j)-Yc(i,j))/dzs;
            Tkzz(i,j)=2.0*(Tk(i-1,j)-Tk(i,j))/dzs;
        else
            Yazz(i,j)=(Ya(i+1,j)-2.0*Ya(i,j)+Ya(i-1,j))/dzs;
            Ybzz(i,j)=(Yb(i+1,j)-2.0*Yb(i,j)+Yb(i-1,j))/dzs;
            Yczz(i,j)=(Yc(i+1,j)-2.0*Yc(i,j)+Yc(i-1,j))/dzs;
            Tkzz(i,j)=(Tk(i+1,j)-2.0*Tk(i,j)+Tk(i-1,j))/dzs;
        end
        
        %% Yaz, Ybz, Tkz, advection term in axial direction
        if (i==1)
            Yaz(i,j)=(Ya(i,j)-Yae)/dz;
            Ybz(i,j)=(Yb(i,j)-Ybe)/dz;
            Ycz(i,j)=(Yc(i,j)-Yce)/dz;
            Tkz(i,j)=(Tk(i,j)-Tkee(1,j))/dz;
        else
            Yaz(i,j)=(Ya(i,j)-Ya(i-1,j))/dz;
            Ybz(i,j)=(Yb(i,j)-Yb(i-1,j))/dz;
            Ycz(i,j)=(Yc(i,j)-Yc(i-1,j))/dz;
            Tkz(i,j)=(Tk(i,j)-Tk(i-1,j))/dz;
        end
        
        %% PDEs 
            Yat(i,j)=Dca(i,j)*(Yar(i,j)+Yarr(i,j)+Yazz(i,j))-v(i,j)*Yaz(i,j);
            Ybt(i,j)=Dcb(i,j)*(Ybr(i,j)+Ybrr(i,j)+Ybzz(i,j))-v(i,j)*Ybz(i,j);
            Yct(i,j)=Dcc(i,j)*(Ycr(i,j)+Ycrr(i,j)+Yczz(i,j))-v(i,j)*Ycz(i,j);
            Tkt(i,j)=Dt(i,j)*(Tkr(i,j)+Tkrr(i,j)+Tkzz(i,j))-v(i,j)*Tkz(i,j);
    end
    
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
            % Beware of Dca(i,jj) follows and dc/dr not = 0
            %Yass(i,jj)=(1/por)*((Dea(i,jj)*(Yas(i,jj+1)-Yas(i,jj))/dss)-(Dca(i,nr)*(Yas(i,jj)-Ya(i,nr))/((drs+dss)/2))); %% drs+dss/2??
            %Ybss(i,jj)=(1/por)*((Deb(i,jj)*(Ybs(i,jj+1)-Ybs(i,jj))/dss)-(Dcb(i,nr)*(Ybs(i,jj)-Yb(i,nr))/((drs+dss)/2)));
            %Ycss(i,jj)=(1/por)*((Dec(i,jj)*(Ycs(i,jj+1)-Ycs(i,jj))/dss)-(Dcc(i,nr)*(Ycs(i,jj)-Yc(i,nr))/((drs+dss)/2)));
            %Tkss(i,jj)=(Ds(i,jj)*(Tks(i,jj+1)-Tks(i,jj))/dss)-(Dt(i,nr)*(Tks(i,jj)-Tk(i,nr))/((drs+dss)/2));
            % Following equation don't work bcos s(jj=1)=0
            %Yass(i,jj)=(1/por)*(1.0/s(jj))*((Dea(i,jj)*(Yas(i,jj+1)-Yas(i,jj))/ds)-(Dca(i,nr)*(Yas(i,jj)-Ya(i,nr))/((dr+ds)/2)));
            %Ybss(i,jj)=(1/por)*(1.0/s(jj))*((Deb(i,jj)*(Ybs(i,jj+1)-Ybs(i,jj))/ds)-(Dcb(i,nr)*(Ybs(i,jj)-Yb(i,nr))/((dr+ds)/2)));
            %Ycss(i,jj)=(1/por)*(1.0/s(jj))*((Dec(i,jj)*(Ycs(i,jj+1)-Ycs(i,jj))/ds)-(Dcc(i,nr)*(Ycs(i,jj)-Yc(i,nr))/((dr+ds)/2)));
            %Tkss(i,jj)=(1.0/s(jj))*((Ds(i,jj)*(Tks(i,jj+1)-Tks(i,jj))/ds)-(Dt(i,nr)*(Tks(i,jj)-Tk(i,nr))/((dr+ds)/2)));
            %Yass(i,jj)=(1/por)*(1.0/s(jj))*((Dea(i,jj)*(Yas(i,jj+1)-Yas(i,jj))/ds)-(Dca(i,nr)*(Yas(i,jj)-Ya(i,nr))/((dr+ds)/2)));
            %Ybss(i,jj)=(1/por)*(1.0/s(jj))*((Deb(i,jj)*(Ybs(i,jj+1)-Ybs(i,jj))/ds)-(Dcb(i,nr)*(Ybs(i,jj)-Yb(i,nr))/((dr+ds)/2)));
            %Ycss(i,jj)=(1/por)*(1.0/s(jj))*((Dec(i,jj)*(Ycs(i,jj+1)-Ycs(i,jj))/ds)-(Dcc(i,nr)*(Ycs(i,jj)-Yc(i,nr))/((dr+ds)/2)));
            %Tkss(i,jj)=(1.0/s(jj))*((Ds(i,jj)*(Tks(i,jj+1)-Tks(i,jj))/ds)-(Dt(i,nr)*(Tks(i,jj)-Tk(i,nr))/((dr+ds)/2)));
            Yass(i,jj)=0;%(1/por)*Dca(i,nr)*(1.0/s(jj))*(Yas(i,jj)-Ya(i,nr))/((dr+ds)/2);
            Ybss(i,jj)=0;%(1/por)*Dcb(i,nr)*(1.0/s(jj))*(Ybs(i,jj)-Yb(i,nr))/((dr+ds)/2);
            Ycss(i,jj)=0;%(1/por)*Dcc(i,nr)*(1.0/s(jj))*(Ycs(i,jj)-Yc(i,nr))/((dr+ds)/2);
            Tkss(i,jj)=0;%Dt(i,nr)*(1.0/s(jj))*(Tks(i,jj)-Tk(i,nr))/((dr+ds)/2);
            % Can't go with 1/0 = indeterminant
        elseif (jj==ns) % Boundary
            Yass(i,jj)=0;       % As dc/ds=0
            Ybss(i,jj)=0;  
            Ycss(i,jj)=0;
            %Tkss(i,jj)=((Tku(i,1)-Tks(i,jj))/((dss+dus)/2))-((Tks(i,jj)-Tks(i,jj-1))/dss);
            Tkss(i,jj)=(1.0/s(jj))*(Tku(i,1)-Tks(i,jj-1))/(2.0*ds);
        else
            Yass(i,jj)=(1.0/s(jj))*(Yas(i,jj+1)-Yas(i,jj-1))/(2.0*ds);
            Ybss(i,jj)=(1.0/s(jj))*(Ybs(i,jj+1)-Ybs(i,jj-1))/(2.0*ds);
            Ycss(i,jj)=(1.0/s(jj))*(Ycs(i,jj+1)-Ycs(i,jj-1))/(2.0*ds);
            Tkss(i,jj)=(1.0/s(jj))*(Tks(i,jj+1)-Tks(i,jj-1))/(2.0*ds);
        end
        
        %% Yasss, Ybsss, Ycsss, Tksss, second derivative diffusion term in radial direction
        if (jj==1)
            Yasss(i,jj)=2.0*(1/por)*((Dea(i,jj)*(Yas(i,jj+1)-Yas(i,jj))/dss)-(Dca(i,nr)*(Yas(i,jj)-Ya(i,nr))/drs)); %% drs+dss/2??
            Ybsss(i,jj)=2.0*(1/por)*((Deb(i,jj)*(Ybs(i,jj+1)-Ybs(i,jj))/dss)-(Dcb(i,nr)*(Ybs(i,jj)-Yb(i,nr))/drs));
            Ycsss(i,jj)=2.0*(1/por)*((Dec(i,jj)*(Ycs(i,jj+1)-Ycs(i,jj))/dss)-(Dcc(i,nr)*(Ycs(i,jj)-Yc(i,nr))/drs));
            Tksss(i,jj)=2.0*(Ds(i,jj)*(Tks(i,jj+1)-Tks(i,jj))/dss)-(Dt(i,nr)*(Tks(i,jj)-Tk(i,nr))/drs);
        elseif (jj==ns) % Boundary
            Yasss(i,jj)=2.0*(Yas(i,jj-1)-Yas(i,jj))/dss;
            Ybsss(i,jj)=2.0*(Ybs(i,jj-1)-Ybs(i,jj))/dss;
            Ycsss(i,jj)=2.0*(Ycs(i,jj-1)-Ycs(i,jj))/dss;
            %Tksss(i,jj)=((Tku(i,1)-Tks(i,jj))/dus)-((Tks(i,jj)-Tks(i,jj-1))/dss);
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
            Yast(i,jj)=Yass(i,jj)+Yasss(i,jj)+(1/por)*(Dea(i,jj)*Yaszz(i,jj)-RCO(i,jj)*H*Av*R*Tks(i,jj)/P);%/r_gtc
            Ybst(i,jj)=Ybss(i,jj)+Ybsss(i,jj)+(1/por)*(Deb(i,jj)*Ybszz(i,jj)+RCO2(i,jj)*H*Av*R*Tks(i,jj)/P);%/r_gtc
            Ycst(i,jj)=Ycss(i,jj)+Ycsss(i,jj)+(1/por)*(Dec(i,jj)*Ycszz(i,jj)-RHC(i,jj)*H*Av*R*Tks(i,jj)/P);%/r_gtc
            Tkst(i,jj)=Tkss(i,jj)+Tksss(i,jj)+Ds(i,jj)*(Tkszz(i,jj))+(-RCO(i,jj)*dHR-RHC(i,jj)*dHRhc(i,jj))*H*Av/(rho_s(i,jj)*Cp_s(i,jj));
            fcvt(i,jj)=-Ra(i,jj)-2*Rb(i,jj)+2*Rc(i,jj)+Re(i,jj)-Rf(i,jj);
            fcCOt(i,jj)=Ra(i,jj)-Rc(i,jj)+Rd(i,jj)-Re(i,jj);
            fcOt(i,jj)=2*Rb(i,jj)-Rc(i,jj)-2*Rd(i,jj)+Re(i,jj)+2*Rf(i,jj);
            fcOOt(i,jj)=Rd(i,jj)-Re(i,jj)-Rf(i,jj);
        else
            Yast(i,jj)=(1/por)*(Dea(i,jj)*(Yass(i,jj)+Yasss(i,jj)+Yaszz(i,jj))-RCO(i,jj)*H*Av*R*Tks(i,jj)/P);%/r_gtc
            Ybst(i,jj)=(1/por)*(Deb(i,jj)*(Ybss(i,jj)+Ybsss(i,jj)+Ybszz(i,jj))+RCO2(i,jj)*H*Av*R*Tks(i,jj)/P);%/r_gtc
            Ycst(i,jj)=(1/por)*(Dec(i,jj)*(Ycss(i,jj)+Ycsss(i,jj)+Ycszz(i,jj))-RHC(i,jj)*H*Av*R*Tks(i,jj)/P);%/r_gtc
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
        if (jjj==1)  
            %Tkuu(i,jjj)=(Du(i,jjj)*(Tku(i,jjj+1)-Tku(i,jjj))/dus)-((Ds(i,ns)*(Tku(i,jjj)-Tks(i,ns))/((dus+dss)/2)));
            Tkuu(i,jjj)=0;%Ds(i,ns)*(1.0/u(jjj))*(Tku(i,jjj)-Tks(i,ns))/((du+ds)/2);
            % Can't go with 1/0 = indeterminant
        elseif (jjj==nu) % Boundary
            Tkuu(i,jjj)=0;
        else
            Tkuu(i,jjj)=(1.0/(u(jjj)))*(Tku(i,jjj+1)-Tku(i,jjj-1))/(2.0*du);
        end

    %% Tkuuu, second derivative thermal diffusion term in radial direction
        if (jjj==1)
            %Tkuuu(i,jjj)=0;
            Tkuuu(i,jjj)=2*(Du(i,jjj)*(Tku(i,jjj+1)-Tku(i,jjj))/dus)-((Ds(i,ns)*(Tku(i,jjj)-Tks(i,ns))/dss));
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
                %for jj=1:ns
            if (t==0)
                History_t=t;
                HisTks=Tks(2,2);
                HisTke=Tkee(1,1);
                HisRa=Ra(2,2);
                HisRb=Rb(2,2);
                HisRc=Rc(2,2);
                HisRd=Rd(2,2);
                HisRe=Re(2,2);
                HisRf=Rf(2,2);
            else
                History_t=[History_t;t];
                HisTks=[HisTks;Tks(2,2)];
                HisTke=[HisTke;Tkee(1,1)];
                HisRa=[HisRa;Ra(2,2)];
                HisRb=[HisRb;Rb(2,2)];
                HisRc=[HisRc;Rc(2,2)];
                HisRd=[HisRd;Rd(2,2)];
                HisRe=[HisRe;Re(2,2)];
                HisRf=[HisRf;Rf(2,2)];
            end
                %end
            end

%% 2D to 1D matrices
for i=1:nz
    % Gas phase
    for j=1:nr    
        ij=(i-1)*nr+j;
        yt(ij)=Yat(i,j);
        yt(ij+nz*nr)=Ybt(i,j);
        yt(ij+2*nz*nr)=Yct(i,j);
        yt(ij+3*nz*nr)=Tkt(i,j);
    end
    % Solid washcoat phase
    for jj=1:ns
        ij=4*nz*nr+(i-1)*ns+jj;
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
        ij=4*nz*nr+8*ns*nz+(i-1)*nu+jjj;
        yt(ij)=Tkut(i,jjj);
    end
end

%% Transpose and steps count
yt=yt';
ncall=ncall+1;
fprintf('\n t = %5d\n',t);
