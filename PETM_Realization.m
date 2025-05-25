clear;clc;
%Calculate the period estimation verification of the actual XPNAV-1 observation data
%Set the number of observation groups according to the actual situation
Z=121;
for j= 1:Z
    %Import the flag bit information required for preprocessing data
    load('flag.mat')
    %Import the actual period of the JB Observatory to facilitate the subsequent calculation of error values
    load('JB.mat')
    %Initialize the matrix variables
    TOA1=[];
    TOA2=[];
    TOA=[];
    %Automatically import the data corrected by Einstein one by one
    str1='tdb_time11.txt';
    num=flag(j,1);
    str11='11';
    str12=num2str(num);
    str1=strrep(str1, str11, str12);
    a1=importdata(str1);
    %Perform MJD time conversion
    DATA=MJD(a1);
    %Automatically import the shapiro delay correction information one by one
    str2='r_shapiroDE432s11.mat';
    str21='11';
    str22=num2str(num);
    str2=strrep(str2, str21, str22);
    load(str2);
    %Automatically import the roemer delay correction information one by one
    str3='r_roemerDE432s11.mat';
    str31='11';
    str32=num2str(num);
    str3=strrep(str3, str31, str32);
    load(str3);
    %Perform delay information correction
    str41='r_roemerDE432s';
    str42=num2str(num);
    str4=strcat(str41, str42);
    str51='r_shapiroDE432s';
    str52=num2str(num);
    str5=strcat(str51, str52);
    DATA=DATA+eval(str4)+eval(str5);
    str61='';
    str62=num2str(flag(j,2));
    str63='.txt';
    str6=strcat(+str61, str62,str63);
    a2=importdata(str6);
    %Calculate the size of the preprocessed data
    [m,n] = size(DATA);
    %Extract the arrival time of photons
    TOA1(:,1)=DATA(:,1)-DATA(1,1)*ones(m,1);
    %Extract the energy carried by photons
    TOA1(:,2)=a2(flag(j,3):flag(j,4),2);
    %Based on the photon energy probability distribution criterion, the energy probability information is generated
    TOA1(:,3)=Energy_probability(TOA1);
    TOA=TOA1;
    [m,n] = size(TOA);
    %Set the number of bins
    bin_num = BIN;
    chi_range = 1;
    %Lower bound of the search cycle
    eP_lowerbound = P1;
    %Upper bound of the search cycle
    eP_upperbound = P2;
    %Search step size
    eP_step = 0.0000000001;
    %Set the proportion coefficient of template fusion, and combine K1 and K2 through the optimization algorithm to achieve the optimal and most stable estimation effect
    K1=10;
    K2=10;
    for P = eP_lowerbound:eP_step:eP_upperbound
        Tb = P/bin_num;
        %Perform the photon epoch folding operation
        profile = Epoch_folding(bin_num,P,Tb,m,TOA);
        %Align the folded contour with the template in phase
        move_profile=move_phase(profile);
        %Standardized folding profile
        normalized_profile=(move_profile-min(move_profile))/(max(move_profile)-min(move_profile));
        %Import the Hilbert conversion template
        load('muban.mat')
        hilbert(16,16)=0;
        for i=1:16
            for b=1:16
                for k=1:256
                    if muban(i,b)==k
                        hilbert(i,b)=normalized_profile(k,1);
                    end
                end
            end
        end
        %Import the standard profile template of the pulsar
        load('matchingstd1.mat')
        load('matchingstd2.mat')
        %Extract the folded profile template of the pulsar
        n1=9;n2=12;n3=3;n4=8;
        n5=1;n6=16;n7=1;n8=16;
        matching1=hilbert(n1:n2,n3:n4);
        matching2=hilbert(n5:n6,n7:n8);
        %Calculate the PCC coefficient
        chi2(chi_range,1) = K1*PCC(matchingstd1,matching1)+K2*PCC(matchingstd2,matching2);
        chi_range =  chi_range+1;
    end
    %Search for the minimized PCC and its index
    [maxchi1,p_pos] = max(chi2(:,1));
    P = (p_pos-1)*eP_step+eP_lowerbound;
    error(j,1)=P;
    error(j,2)=abs(P-JB(j,1))*1e9
end