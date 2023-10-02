clc;clear;close all;%designed by hzLiu, Apr.20th,2023

MAXOBS=60;%PRN<=60
% X(m),Y(m),Z(m),pseudorange(m)
obsdata=zeros([MAXOBS,5]);
satid=nan([MAXOBS,1]);
X=zeros([1,4]);
nepoch=0;
c=299792458;

m0sum=[];
KXsum=[];
Xsum=[];
m0vec=[];
Km0vec=[];
neuPos=[];
KneuPos=[];
clock=[];
Kclock=[];
timeSeq=[];
kfQx0=[];
kfX0=zeros(1,4);
kfx=zeros(4,1);
kfm0=[];

filename='../txtDir/CUSV_20212220_BDS_sgnac_no_eph_error.txt';
outfile='../txtDir/outpos.txt';
outfile1='../txtDir/Koutpos.txt';
fid=fopen(filename,'r');
fop=fopen(outfile,'w+');
fop1=fopen(outfile1,'w+');

disp(filename);
if(fid==-1 || fop==-1)
    disp('can not find or open the file: wrong!');
else
    %% read file
    while ~feof(fid)
        nepoch=nepoch+1;
        tline=fgetl(fid);
        firstline=tline;%copy the first line
        if(strcmp(tline(1),'#'))
            %epoch, gpsweek, gpssecond, number of obs
            obshead=sscanf(tline(2:end),"%f");
            for i=1:obshead(4)
                tline=fgetl(fid);
                % satid=tline(1:3); %omit prn
                % X(m),Y(m),Z(m),pseudorange(m),variance(m)
                obsdata(i,:)=sscanf(tline(4:end),"%f",[1,5]);
            end
        end
        timeSeq=[timeSeq;obshead(3)];

        %% calculate
        %Static Kalman Filter
        if nepoch>=2
            if nepoch==2
                X_est=kfX0.';
                Qxk=kfQx0;
                x_est=kfx0;
                m0k=kfm0;
                allsatsum=kfsatnum;
                Km0vec=[Km0vec;m0k];
            end
            m=0;
            kvecbx=zeros(obshead(4),1);
            kvecby=zeros(obshead(4),1);
            kvecbz=zeros(obshead(4),1);
            kvecbc=-1.*ones(obshead(4),1);
            z=zeros(obshead(4),1);
            ks0=zeros(obshead(4),1);
            %default:sigma0=1,then Q=D
            q=zeros(obshead(4));
            kmaxCorr=1;
            for i=1:obshead(4)
                    Xs=obsdata(i,1);
                    Ys=obsdata(i,2);
                    Zs=obsdata(i,3);
                    ks0(i,1)=sqrt((Xs-X_est(1))*(Xs-X_est(1))+(Ys-X_est(2))*(Ys-X_est(2))+(Zs-X_est(3))*(Zs-X_est(3)));
                    kvecbx(i,1)=-(Xs-X_est(1))/ks0(i,1);
                    kvecby(i,1)=-(Ys-X_est(2))/ks0(i,1);
                    kvecbz(i,1)=-(Zs-X_est(3))/ks0(i,1);
                    z(i,1)=obsdata(i,4)-ks0(i,1)+X_est(4);
                    q(i,i)=obsdata(i,5);
            end
            w=inv(q);
            allsatsum=allsatsum+obshead(4);
            %iteration: one time is enough
            while abs(kmaxCorr)>1e-9&&m<1 
                h=[kvecbx kvecby kvecbz kvecbc];
                deltaz=z-h*x_est;
                %deltaz=z-h*X_est;%wrong
                K=Qxk*h.'/(inv(w)+h*Qxk*h.');
                xk1=x_est+K*deltaz;
                %Xk1=X_est+K*deltaz;%wrong
                Xk1=X_est+xk1;
                Qk1=Qxk-K*h*Qxk;          
                vk1=h*xk1-z;
                m0k1=sqrt(1/(allsatsum-4)*(m0k*m0k*(allsatsum-obshead(4)-4)+deltaz.'*K.'*Qxk*K*deltaz+vk1.'*w*vk1));
                Dk1=m0k1*m0k1*Qk1;
                %renewal
                x_est=xk1;
                %X_est=Xk1;%no need for the estimated value renewal
                Qxk=Qk1;
                kmaxCorr=max(abs(K*deltaz));
                m=m+1;
            end
            disp("KalmanEpoch "+string(nepoch)+" Iteration Times: "+string(m));
            KXsum=[KXsum;Xk1.'];
            Km0vec=[Km0vec;m0k1];
            KXmean=mean(KXsum(:,1));
            KYmean=mean(KXsum(:,2));
            KZmean=mean(KXsum(:,3));
            %StaPos=[Xmean Ymean Zmean];
            KStaPos=[-1132915.01648681 6092528.50388968 1504633.16777129];
            KstdX=std(KXsum(:,1));
            KstdY=std(KXsum(:,2));
            KstdZ=std(KXsum(:,3));

            outdata1=[Xk1(1),Xk1(2),Xk1(3),-1.*Xk1(4)];
            fprintf(fop1,string(firstline)+'\r\n');
            fprintf(fop1,'%14s %14s %14s %14s\r\n','X(m)','Y(m)','Z(m)','T(m)');
            fprintf(fop1,'%14.4f %14.4f %14.4f %14.4f\r\n',outdata1);
            fprintf(fop1,"Posterior Unit Weight Error(m):%13.8f\r\n",m0k1);
            fprintf(fop1,"Posterior Estimated Variance(m^2):\r\n");
            fprintf(fop1,'%14.4f %14.4f %14.4f %14.4f\r\n',Dk1(1,:));
            fprintf(fop1,'%14.4f %14.4f %14.4f %14.4f\r\n',Dk1(2,:));
            fprintf(fop1,'%14.4f %14.4f %14.4f %14.4f\r\n',Dk1(3,:));
            fprintf(fop1,'%14.4f %14.4f %14.4f %14.4f\r\n',Dk1(4,:));
        end

        %Single Epoch Cauculation by Least Squares
        vecbx=zeros(obshead(4),1);
        vecby=zeros(obshead(4),1);
        vecbz=zeros(obshead(4),1);
        vecbc=-1.*ones(obshead(4),1);
        l=zeros(obshead(4),1);
        s0=zeros(obshead(4),1);
        %default:sigma0=1,then Q=D
        Q=zeros(obshead(4));
        x=zeros(4,1);
        n=0;
        maxCorr=1;
        while abs(maxCorr)>1e-9&&n<10
            for i=1:obshead(4)
                Xs=obsdata(i,1);
                Ys=obsdata(i,2);
                Zs=obsdata(i,3);
                s0(i,1)=sqrt((Xs-X(1))*(Xs-X(1))+(Ys-X(2))*(Ys-X(2))+(Zs-X(3))*(Zs-X(3)));
                vecbx(i,1)=-(Xs-X(1))/s0(i,1);
                vecby(i,1)=-(Ys-X(2))/s0(i,1);
                vecbz(i,1)=-(Zs-X(3))/s0(i,1);
                l(i,1)=obsdata(i,4)-s0(i,1)+X(4);
                Q(i,i)=obsdata(i,5);
            end
            B=[vecbx vecby vecbz vecbc];
            %P=inv(Q);
            x=((B.')/Q*B)\((B.')/Q*l);
            X(1)=X(1)+x(1);
            X(2)=X(2)+x(2);
            X(3)=X(3)+x(3);
            X(4)=X(4)+x(4);
            maxCorr=max(x(1:3));
            n=n+1;
        end
        V=B*x-l;
        m0=sqrt((V.'/Q*V)/(obshead(4)-4));
        Qxx=inv(B.'/Q*B);
        Dxx=m0*m0.*Qxx;
        Xsum=[Xsum;X];
        m0vec=[m0vec;m0];

        Xmean=mean(Xsum(:,1));
        Ymean=mean(Xsum(:,2));
        Zmean=mean(Xsum(:,3));
        %StaPos=[Xmean Ymean Zmean];
        StaPos=[-1132915.01648681 6092528.50388968 1504633.16777129];
        stdX=std(Xsum(:,1));
        stdY=std(Xsum(:,2));
        stdZ=std(Xsum(:,3));

        if nepoch==1 %give the initialized INFO for the static kalman filtering
            kfX0=X;
            kfQx0=Qxx;
            kfx0=x;
            kfm0=m0;
            kfDxx=(m0*m0).*kfQx0;
            kfsatnum=obshead(4);
            allsatsum=kfsatnum;
            outdata1=[X(1),X(2),X(3),-1.*X(4)];
            KXsum=[KXsum;X];
            
            fprintf(fop1,string(firstline)+'\r\n');
            fprintf(fop1,'%14s %14s %14s %14s\r\n','X(m)','Y(m)','Z(m)','T(m)');
            fprintf(fop1,'%14.4f %14.4f %14.4f %14.4f\r\n',outdata1);
            fprintf(fop1,"Posterior Unit Weight Error(m):%13.8f\r\n",kfm0);
            fprintf(fop1,"Posterior Estimated Variance(m^2):\r\n");
            fprintf(fop1,'%14.4f %14.4f %14.4f %14.4f\r\n',kfDxx(1,:));
            fprintf(fop1,'%14.4f %14.4f %14.4f %14.4f\r\n',kfDxx(2,:));
            fprintf(fop1,'%14.4f %14.4f %14.4f %14.4f\r\n',kfDxx(3,:));
            fprintf(fop1,'%14.4f %14.4f %14.4f %14.4f\r\n',kfDxx(4,:));
        end
        disp("Epoch "+string(nepoch)+" Iteration Times: "+string(n));
        outdata=[X(1),X(2),X(3),X(4)];
        
        fprintf(fop,string(firstline)+'\r\n');
        fprintf(fop,'%14s %14s %14s %14s\r\n','X(m)','Y(m)','Z(m)','T(m)');
        fprintf(fop,'%14.4f %14.4f %14.4f %14.4f\r\n',outdata);
        fprintf(fop,"Posterior Unit Weight Error(m):%13.8f\r\n",m0);
        fprintf(fop,"Posterior Estimated Variance(m^2):\r\n");
        fprintf(fop,'%14.4f %14.4f %14.4f %14.4f\r\n',Dxx(1,:));
        fprintf(fop,'%14.4f %14.4f %14.4f %14.4f\r\n',Dxx(2,:));
        fprintf(fop,'%14.4f %14.4f %14.4f %14.4f\r\n',Dxx(3,:));
        fprintf(fop,'%14.4f %14.4f %14.4f %14.4f\r\n',Dxx(4,:));
        
        %% Clear the data and prepare the next epoch
        obsdata = obsdata * 0; 

    end

    %% draw imgs
    Xsum(:,4)=-1.*Xsum(:,4);
    for i=1:length(Xsum)
        neuPos=[neuPos;xyz2neu(StaPos,Xsum(i,1:3))];
        clock=[clock;Xsum(i,4)];
    end
    KXsum(:,4)=-1.*KXsum(:,4);
    for i=1:length(KXsum)
        KneuPos=[KneuPos;xyz2neu(KStaPos,KXsum(i,1:3))];
        Kclock=[Kclock;KXsum(i,4)];
    end

    Nmean=mean(neuPos(:,1));
    Emean=mean(neuPos(:,2));
    Umean=mean(neuPos(:,3));
    rmsN=rms(neuPos(:,1));
    rmsE=rms(neuPos(:,2));
    rmsU=rms(neuPos(:,3));
    neuPos=[timeSeq neuPos];
    xyzPos=[timeSeq Xsum];
    clock=[timeSeq clock];
    m0sum=[timeSeq m0vec Km0vec];

    neuName="neuPos";
    neuDraw(neuPos,neuName);
    clockName="Receiver clock difference";
    RclockdiffDraw(clock,clockName);
    m0name="Posterior Unit Weight Error";
    m0Draw(m0sum,m0name);

    KNmean=mean(KneuPos(:,1));
    KEmean=mean(KneuPos(:,2));
    KUmean=mean(KneuPos(:,3));
    KrmsN=rms(KneuPos(:,1));
    KrmsE=rms(KneuPos(:,2));
    KrmsU=rms(KneuPos(:,3));
    KneuPos=[timeSeq KneuPos];
    KxyzPos=[timeSeq KXsum];
    Kclock=[timeSeq Kclock];
    KneuName="KneuPos";
    neuDraw(KneuPos,KneuName);
    KclockName="KReceiver clock difference";
    RclockdiffDraw(Kclock,KclockName);

    xyzName="xyzPos";
    xyzDraw(xyzPos,xyzName);
    KxyzName="KxyzPos";
    xyzDraw(KxyzPos,KxyzName);

    pos3dname="blhPos3d";
    blhpos3d(Xsum(:,1),Xsum(:,2),Xsum(:,3),pos3dname,StaPos);
    Kpos3dname="KblhPos3d";
    blhpos3d(Xsum(:,1),Xsum(:,2),Xsum(:,3),Kpos3dname,KStaPos);

    cats6=['stdX';'stdY';'stdZ';'rmsN';'rmsE';'rmsU';];
    vals6=[stdX,stdY,stdZ,rmsN,rmsE,rmsU; ...
          KstdX,KstdY,KstdZ,KrmsN,KrmsE,KrmsU];
    histogramname6="Hist1_Kalman_SingleLeastSquares";
    histogram6(cats6,vals6,histogramname6);

    cats3=['Xmean';'Ymean';'Zmean';'Nmean';'Emean';'Umean';];
    vals3=[Xmean,Ymean,Zmean,Nmean,Emean,Umean; ...
          KXmean,KYmean,KZmean,KNmean,KEmean,KUmean];
    histogramname3="Hist2_Kalman_SingleLeastSquares";
    histogram3(cats3,vals3,histogramname3);

    %% console display
    disp('Xmean='+string(Xmean)+'m');
    disp('Ymean='+string(Ymean)+'m');
    disp('Zmean='+string(Zmean)+'m');
    disp('stdX='+string(stdX)+'m');
    disp('stdY='+string(stdY)+'m');
    disp('stdZ='+string(stdZ)+'m');
    disp('Nmean='+string(Nmean)+'m');
    disp('Emean='+string(Emean)+'m');
    disp('Umean='+string(Umean)+'m');
    disp('rmsN='+string(rmsN)+'m');
    disp('rmsE='+string(rmsE)+'m');
    disp('rmsU='+string(rmsU)+'m');

    disp('KXmean='+string(KXmean)+'m');
    disp('KYmean='+string(KYmean)+'m');
    disp('KZmean='+string(KZmean)+'m');
    disp('KstdX='+string(KstdX)+'m');
    disp('KstdY='+string(KstdY)+'m');
    disp('KstdZ='+string(KstdZ)+'m');
    disp('KNmean='+string(KNmean)+'m');
    disp('KEmean='+string(KEmean)+'m');
    disp('KUmean='+string(KUmean)+'m');
    disp('KrmsN='+string(KrmsN)+'m');
    disp('KrmsE='+string(KrmsE)+'m');
    disp('KrmsU='+string(KrmsU)+'m');

    disp('process over');
end
fclose(fid); fclose(fop);fclose(fop1);