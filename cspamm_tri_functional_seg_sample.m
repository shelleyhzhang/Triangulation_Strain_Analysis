function cspamm_tri_functional_segment
%%%%%%%%%%
%%Functional Description:
%%% 1. calculates the mechanical strains of the myocardial wall based on
%%% trianglular finite elements
%%% 2. divide the wall into 4 segment and delineate iron labeled cells
%%% author: shelleyzhang@gmail.com
%%% date: Aug 2010
%%%%%%%%%%
%% Load data
folder = ['C:\Varian\HUPC\CardiacDense\nrat_baseline\082809_nrat3_leftrightear2hole_\']; pss='pss-0.3_01.fid';
nti = read_procpar([folder 'dir1tag0_' pss ],'ne ');
COx = (reform(load_fid([folder 'dir1tag0_' pss])-load_fid([folder 'dir1tag2_' pss]),nti));
COy = (reform(load_fid([folder 'dir2tag0_' pss])-load_fid([folder 'dir2tag2_' pss]),nti));
rekx =(reform(load_fid([folder  'dir1tag2r_' pss]),nti));
reky = (reform(load_fid([folder 'dir1tag2r_' pss]),nti)); 
spacing = read_procpar([folder 'dir1tag0_' pss],'dtag ');
lro = read_procpar([folder 'dir1tag0_' pss],'lro ');
XI = size(COx,1); YI = size(COx,2);

kshift=-5;
COx=circshift(COx,[kshift 0 0]); COy=circshift(COy,[kshift 0 0]);
rekx=circshift(rekx,[kshift 0 0]);reky=circshift(reky,[kshift 0 0]);
cut=50;
COx(1:cut,:,:)=0;COx(cut,:,:)=COx(cut+1,:,:)*2/3; COx(cut-1,:,:)=COx(cut+1,:,:)/3;
COy(1:cut,:,:)=0;COy(cut,:,:)=COy(cut+1,:,:)*2/3; COy(cut-1,:,:)=COy(cut+1,:,:)/3;
COy = imrotate(COy,-90); 
rotflagy=0; rotflagx=90;
reky = imrotate(reky,rotflagy);
rekx = imrotate(rekx,rotflagx);
% fft k-space
for n=1:nti
    Dx(:,:,n) = fftshift(fft2((COx(:,:,n)),XI,YI));
    Dy(:,:,n) = fftshift(fft2((COy(:,:,n)),XI,YI));
    refx(:,:,n) = fftshift(fft2(rekx(:,:,n),XI,YI));
    refy(:,:,n) = fftshift(fft2(reky(:,:,n),XI,YI));

    cc = normxcorr2(abs(Dx(:,:,n))/max(max(abs(Dx(:,:,n)))), abs(Dy(:,:,n))/max(max(abs(Dy(:,:,n)))));
    [max_cc, imax] = max(abs(cc(:))); [ypeak, xpeak] = ind2sub(size(cc),imax(1));
    corr_offset = [ (ypeak-XI) (xpeak-YI) ];
    fprintf(['corr Dx and Dy: ' num2str(corr_offset) '\n'])
    Dy(:,:,n) = circshift(Dy(:,:,n),[-corr_offset(1) -corr_offset(2) 0]);
end

yshift = 0;
Dx = circshift(Dx,[0 yshift 0]);
Dy = circshift(Dy,[0 yshift 0]);
refx = circshift(refx,[0 yshift 0]);
refy = circshift(refy,[0 yshift 0 ]);
Dx    = abs(Dx).*exp(1i*(angle(Dx)-angle(imrotate(refx,-rotflagx)))); 
Dy    = abs(Dy).*exp(1i*(angle(Dy)-angle(imrotate(refx,-rotflagx))));
end

refy=imrotate(refy,-rotflagy);
refx=imrotate(refx,-rotflagx);

[cropped rect]= imcrop(abs(refy));
Dxc = imcrop(real(Dx),rect)+1i*imcrop(imag(Dx),rect);
Dyc = imcrop(real(Dy),rect)+1i*imcrop(imag(Dy),rect);
figure,subplot(2,1,1),montager(angle(Dxc),[]),title('angle of Dx')
subplot(2,1,2),montager(angle(Dyc),[-pi pi]),title('angle of Dy')

%% Contouring 
scrsz = get(0,'ScreenSize');
h=figure('Units','pixels','Position',[0, 0,scrsz(3),scrsz(4)]);
a0=reform(load_fid([folder 'dir1tag0_' pss]),nti,'f',XI,YI);
for n=[1:nti]
    m=max(abs(Dx(:)));
    figure(h),
    subplot(1,2,1),imshow(abs(refy(:,:,n)),[0 max(abs(refy(:)))/1.3]),title(['cardiac phase ' num2str(n)])
    subplot(1,2,2), imshow(abs(Dx(:,:,n)),[0 m/1.5],'initialmagnification','fit');
    disp('####   Hint 1    #####')
    disp('Select epi and endo ROI, press q to confirm ROI or any other key to repeat drawing, followed by the Enter key')
    user_entry = '';  
    while (~strcmp(user_entry,'q')) && (~strcmp(user_entry,'r'))
        user_entry = input('Press q or Q after selecting ROI: ','s'); 
        if (user_entry=='r')
            subplot(1,2,1);
        end
        [epi(:,:,n),xi,yi] = roipoly;
        hold on; plot(xi,yi,'r');
        endo(:,:,n)= roipoly;
        
    end            
    disp('####   end of Hint 1    ####');
    hcontour(:,:,n)=epi(:,:,n)-endo(:,:,n);
end
cine(hcontour,4,2)
%%  Locating iron spot
cine(abs(refx),3,5)
for n=[1:nti]
    imshow(abs(refx(:,:,n)),[],'initialmagnification','fit'),title(num2str(n),'fontsize',12)
    ironspot(:,:,n)=roipoly;        
end
cine(cat(4,ironspot.*hcontour,abs(refx)/max(abs(refx(:))),zeros(size(refx))),3)
%% Region Growing phase unwrapping
clear sy sx
figure,subplot(2,1,1),montager(angle(Dxc)),subplot(2,1,2),montager(angle(Dyc))
if ~exist('rect')
       figure
       [cropped rect] = imcrop(abs(refy(:,:,1)),[]);
       hbw=imcrop(hcontour,rect);
else
    hbw=circshift(imcrop(hcontour,rect),[0 0 0]);
end
cropped = imcrop(abs(refy),rect);
Dxc = imcrop(real(Dx),rect)+1i*imcrop(imag(Dx),rect);
Dyc = imcrop(real(Dy),rect)+1i*imcrop(imag(Dy),rect);
hbw=imcrop(hcontour,rect);
for n=[1:nti]
       chestwall=zeros(size(cropped(:,:,1)));
       if ~exist('xpoint')
           [syy(:,:,n),xpoint,ypoint]= flood_fill_dense_mask(Dyc(:,:,n),hbw(:,:,n)); 
       else
           [syy(:,:,n)]= flood_fill_dense_mask(Dyc(:,:,n),hbw(:,:,n),xpoint,ypoint);
       end
       xp(:,n)=[xpoint ypoint]';
       figure
       subplot(1,3,1),imshow(angle(Dyc(:,:,n)),[-pi pi]),title('orig','FontSize',12)
       subplot(1,3,2),imshow(syy(:,:,n),[]),colorbar,title('unwrapped','FontSize',12)
       subplot(1,3,3),imshow(syy(:,:,n)-angle(Dyc(:,:,n)),[]),colorbar,title('Diff','FontSize',12)
       if ~exist('xpointx')
            [sxx(:,:,n),xpointx,ypointx] = flood_fill_dense_mask(Dxc(:,:,n),hbw(:,:,n));
       else
            [sxx(:,:,n)] = flood_fill_dense_mask(Dxc(:,:,n),hbw(:,:,n),xpointx,ypointx);
       end
       figure
       subplot(1,3,1),imshow(angle(Dxc(:,:,n)),[-pi pi]),title('orig','FontSize',12)
       subplot(1,3,2),imshow(sxx(:,:,n),[]),colorbar,title('unwrapped phase')
       subplot(1,3,3),imshow(sxx(:,:,n)-angle(Dxc(:,:,n)),[]),colorbar,title('difference between the two')
end
figure,subplot(2,1,1),montager(sxx),colorbar,subplot(2,1,2),montager(syy),colorbar
%% Quiver displacement maps
clear COx COy rekx reky
hcontour=epi-endo;hbw=imcrop(hcontour,rect);
cropped = imcrop(abs(refy),rect);
scrsz = get(0,'ScreenSize');
ke=lro/XI/(spacing)  %cycle per pixel

clear ux uy
for j=1:nti
    n=mod(j-1,nti)+1;
    ux = sxx(:,:,n)/(ke*(2*pi)); uy = syy(:,:,n)/(ke*(2*pi));
    figure('Position',[scrsz(3)/10 scrsz(4)/10 scrsz(3)/2 scrsz(3)/2])
    imshow(cropped(:,:,n),[],'InitialMagnification','fit'),title(num2str(n),'FontSize',12)
    hold on
    ur=size(ux);
    uxo=(ux-mean(nonzeros(ux.*hbw(:,:,n)))).*hbw(:,:,n);
    uyo=(uy-mean(nonzeros(uy.*hbw(:,:,n)))).*hbw(:,:,n);
    dispvec=complex(uxo,uyo); %in vnmrj2.3C
    dispvec=complex(-uyo,uxo);
    init_pos=(complex(meshgrid(1:ur(2),1:ur(1)),meshgrid(1:ur(1),1:ur(2))')-dispvec);
    %  ------>real(init_pos)
    %  |
    %  |
    % \|/ imag(init_pos)
    quiver(real(init_pos),imag(init_pos),real(dispvec),imag(dispvec),0,'y');
    hold off
end

%%  Strain calculation and divide functional sectors
%phase = displacement * Ke, Ke= delta_x/spacing=40mm/128pixel /2mm =0.15625 cycle/pixel
hcontour=epi-endo; hbw=imcrop(hcontour,rect);
if ~exist('ke')
    ke=lro/XI/(spacing);
end
if ~exist('ironspotc')
     ironspotc = zeros(size(hbw));
end
ur=size(sxx(:,:,1));
cropped = imcrop(abs(refy),rect);

optflag=1; %0 - linear strain;1 - finite strain
ironflag =0; %0 - without iron labeled cells; 1 - with iron labeled cells
nsec=6; %number of phases
Err_reg =zeros(nsec,nti);
Ecc_reg =zeros(nsec,nti);
disp_reg=zeros(nsec,nti);
E1_reg  =zeros(nsec,nti);
E2_reg  =zeros(nsec,nti);
section = zeros(ur(1),ur(2),nsec);
E1  = zeros(1500,nti);
E2  = E1;
Ecc = E1;
Err = E1;
E1_theta=E1;

for n=1:nti
    LVcen1 = regionprops(double(hbw(:,:,1)),'Centroid'); LVcen1 = LVcen1.Centroid;
    LVcen  = regionprops(double(hbw(:,:,n)),'Centroid'); LVcen = round(LVcen.Centroid);    
    ux = sxx(:,:,n)/(ke*(2*pi));  uy = syy(:,:,n)/(ke*(2*pi));
    %  <----syy----|     
    %              |
    %             \|/  sxx
    uxo=(ux-mean(nonzeros(ux.*hbw(:,:,n)))).*hbw(:,:,n);
    uyo=(uy-mean(nonzeros(uy.*hbw(:,:,n)))).*hbw(:,:,n);
    dispvec = complex(uxo,-uyo);
  
    [px,py] = find(hbw(:,:,n));
    dt = DelaunayTri(px,py);
    ic = incenters(dt);
    c=1;
    pz = zeros(numel(px),1);
    while c<=size(ic,1) 
        mid1 = [mean(px(dt.Triangulation(c,1:2))) mean(py(dt.Triangulation(c,1:2)))];
        mid2 = [mean(px(dt.Triangulation(c,2:3))) mean(py(dt.Triangulation(c,2:3)))];
        mid3 = [mean(px(dt.Triangulation(c,1:2:3))) mean(py(dt.Triangulation(c,1:2:3)))];
        tf   = or(ismember(floor([mid1;mid2;mid3]),[px py],'rows'),...
                  ismember( ceil([mid1;mid2;mid3]),[px py],'rows'));
        while sum(tf)<2
            %fprintf('out of myocardium zone %dth triangle\n',c)
            ic(c,:)=[1 1];
            Ecc(c,n)=-100;
            c = c+1;
            mid1 = [mean(px(dt.Triangulation(c,1:2))) mean(py(dt.Triangulation(c,1:2)))];
            mid2 = [mean(px(dt.Triangulation(c,2:3))) mean(py(dt.Triangulation(c,2:3)))];
            mid3 = [mean(px(dt.Triangulation(c,1:2:3))) mean(py(dt.Triangulation(c,1:2:3)))];
            tf   = or(ismember(floor([mid1;mid2;mid3]),[px py],'rows'),...
                      ismember( ceil([mid1;mid2;mid3]),[px py],'rows'));
        end
       Uix = real(dispvec(px(dt.Triangulation(c,1)), py(dt.Triangulation(c,1))));
       Uiy = imag(dispvec(px(dt.Triangulation(c,1)), py(dt.Triangulation(c,1))));
       Ujx = real(dispvec(px(dt.Triangulation(c,2)), py(dt.Triangulation(c,2))));
       Ujy = imag(dispvec(px(dt.Triangulation(c,2)), py(dt.Triangulation(c,2))));
       Ukx = real(dispvec(px(dt.Triangulation(c,3)), py(dt.Triangulation(c,3))));
       Uky = imag(dispvec(px(dt.Triangulation(c,3)), py(dt.Triangulation(c,3))));
       U = [Uix Uiy Ujx Ujy Ukx Uky]';
       
       Xi = px(dt.Triangulation(c,1))-Uix;
       Xj = px(dt.Triangulation(c,2))-Ujx;
       Xk = px(dt.Triangulation(c,3))-Ukx;
       Yi = py(dt.Triangulation(c,1))-Uiy;
       Yj = py(dt.Triangulation(c,2))-Ujy;
       Yk = py(dt.Triangulation(c,3))-Uky;
       tri_area2 = Xi*(Yj-Yk)+Xj*(Yk-Yi)+Xk*(Yi-Yj); 
       
       B = [Yj-Yk   0   Yk-Yi   0   Yi-Yj   0
            Xk-Xj   0   Xi-Xk   0   Xj-Xi   0
              0   Yj-Yk   0   Yk-Yi   0   Yi-Yj
              0   Xk-Xj   0   Xi-Xk   0   Xj-Xi]/tri_area2;
    
       E = B*U;
       if optflag
            opt='finite';
            Em = [E(1)+0.5*(E(1)^2+E(3)^2) 0.5*(E(2)+E(3)+E(1)*E(2)+E(3)*E(4))
                  0.5*(E(2)+E(3)+E(1)*E(2)+E(3)*E(4)) E(4)+0.5*(E(2)^2+E(4)^2)];
            starow=12;
        else
            opt='linear';
            Em = [     E(1)        0.5*(E(2)+E(3))
                   0.5*(E(2)+E(3))       E(4)      ];       
            starow=21;
        end
       [v,d]=eigs(Em);
       d=diag(d);%descending in eigenvalue magnitudes
       Eig_max=d(d==max(d)); Eig_min=d(d==min(d));%
       Eig_max_vector = complex(v(1,d==max(d)), v(2,d==max(d)));
       Eig_min_vector = complex(v(1,d==min(d)), v(2,d==min(d)));
       radvec = complex(ic(c,1), ic(c,2))-complex(LVcen(2), LVcen(1));
       radvec = radvec/norm(radvec);
      
       theta = angle(radvec/1i);
       %http://mathworld.wolfram.com/RotationMatrix.html
       %rotm = [cos(theta) -sin(theta);sin(theta) cos(theta)];
       %here x is downward and y is rightward
       rotm = [cos(theta) sin(theta);-sin(theta) cos(theta)];
       %Mase-continuum mechanics
       Emr = rotm*Em*rotm';
       if length(Eig_max_vector)<2
           theta1 = angle(radvec/Eig_max_vector);theta2 = angle(-radvec/Eig_max_vector);
           theta = theta2;
           if abs(theta1)<abs(theta2)
               theta = theta1;
           end
       end
       if Eig_min>0.5
           fprintf('Eig_min %.2f',Eig_min)
       end
       if Eig_max<3&&Eig_max>-0.5&&Eig_min<0.5&&Eig_min>-3
           E1_theta(c,n) = theta*180/pi;
           E1(c,n) = Eig_max; E2(c,n) = Eig_min;  
           Err(c,n)= max(diag(Emr)); Ecc(c,n)= min(diag(Emr));  
       end
       c=c+1;
    end
    cp=cropped(:,:,n)'/max(cropped(:)); 
    cp=cat(3,cp,cp,cp);
    figure
    imshow(cp,[],'initialmagnification','fit'),hold on
    trisurf(dt.Triangulation,px(:),py(:),pz(:),Ecc(1:length(dt.Triangulation),n));
    colormap(jet)
    set(gca,'Clim',[-0.5 0.5])
    view(2)
    if ~exist('yp')
        figure,imshow(cropped(:,:,1),[0 max(cropped(:))/1.1],'initialmagnification','fit');
        impoint(gca,round(LVcen(1)),round(LVcen(2)));
        uiwait(msgbox('Point1 border-septum Point2 border-lateral'));
        [yp,xp] = ginput(2);
    end
    my=hbw(:,:,n);
    ypt = [LVcen(1);LVcen(1)]+6*(yp-[LVcen(1);LVcen(1)]);
    xpt = [LVcen(2);LVcen(2)]+6*(xp-[LVcen(2);LVcen(2)]) ;
    r1=complex(xpt(1)-LVcen(2),ypt(1)-LVcen(1));
    r2=complex(xpt(2)-LVcen(2),ypt(2)-LVcen(1));
    sepang=abs(angle(r2/r1));
    mypt= LVcen(1)+20*(mean(ypt)-LVcen(1));
    mxpt= LVcen(2)+20*(mean(xpt)-LVcen(1));

    radius=round(1.0*max(sqrt((yp(1)-LVcen(1))^2+(xp(1)-LVcen(2))^2),sqrt((yp(2)-LVcen(1))^2+(xp(2)-LVcen(2))^2)));
    h = fspecial('disk',radius); h(h>0)=1;
    pie=zeros(ur);
    pie(LVcen(2)-radius:LVcen(2)+radius,LVcen(1)-radius:LVcen(1)+radius)=h;
    section=zeros(ur(1),ur(2),5);
    %1=infarct
    section(:,:,1)=roipoly(pie,[LVcen(1) ypt(2) mypt ypt(1)],[LVcen(2) xpt(2) mxpt xpt(1)]).*my;
    %2=border-septum
    pt2=r1*exp(1i*pi/6);
    section(:,:,2)=roipoly(pie,[LVcen(1) ypt(1) LVcen(1)+imag(pt2)],[LVcen(2) xpt(1) LVcen(2)+real(pt2)]).*my;
    %3=border-lateral
    pt3=r2*exp(-1i*pi/6);
    section(:,:,3)=roipoly(pie,[LVcen(1) ypt(2) LVcen(1)+imag(pt3)],[LVcen(2) xpt(2) LVcen(2)+real(pt3)]).*my;
    %4=remote
    section(:,:,4)=my-sum(section(:,:,1:3),3).*my;
    %5=iron
    if (ironflag)
         figure,imshow(cropped(:,:,n),[0 max(cropped(:))/1.4],'initialmagnification','fit');
         title('Points to enclose IRON area');
         section(:,:,5)=roipoly.*my;
         ironspotc(:,:,n)=section(:,:,5);
    else
        section(:,:,5)=ironspotc(:,:,n);
    end
    %6=whole slice
    section(:,:,6)=my;
    croppedn=cropped(:,:,n)/max(cropped(:));
    figure,
    if ironflag 
        subplot(1,2,2),imshow(cat(3,section(:,:,5)+croppedn,croppedn,croppedn))
        subplot(1,2,1)
    end
    imshow(cat(3,section(:,:,4)/3+section(:,:,1)/3+croppedn,...
                 section(:,:,2)/3+section(:,:,1)/3+croppedn,...
                 section(:,:,3)/3+section(:,:,1)/3+croppedn),[],...
                 'initialmagnification','fit')
    title('infarct(w) bor-sep(g) bor-lat(b) remote(r)','fontsize',18);
    for sec=1:nsec
            [rx, ry]  = find(section(:,:,sec));
            rtf       = ismember(round(ic),[rx ry],'rows');
            f=2;
            if sec==1
               [Err_reg(sec,n),threErr]= filnstd(rtf.*Err(1:length(rtf),n),f);
               [Ecc_reg(sec,n),threEcc]= filnstd(rtf.*Ecc(1:length(rtf),n),f);
               [E1_reg(sec,n),threE1]  = filnstd(rtf.*E1(1:length(rtf),n),f);
               [E2_reg(sec,n),threE2]  = filnstd(rtf.*E2(1:length(rtf),n),f);
               [disp_reg(sec,n),thredp]= filnstd(section(:,:,sec).*abs(dispvec)*lro*10/XI,f);
               fprintf('\n threErr=[%.2f %.2f], threEcc=[%.2f %.2f], threE1=[%.2f %.2f], threE2=[%.2f %.2f], thredp=[%.2f %.2f]\n',...
                        threErr, threEcc,threE1,threE2,thredp);
            else
                Err_reg(sec,n) = filnstd(rtf.*Err(1:length(rtf),n),f,threErr);
                Ecc_reg(sec,n) = filnstd(rtf.*Ecc(1:length(rtf),n),f,threEcc);
                E1_reg (sec,n) = filnstd(rtf.*E1(1:length(rtf),n),f,threE1);
                E2_reg (sec,n) = filnstd(rtf.*E2(1:length(rtf),n),f,threE2);
                disp_reg(sec,n)= filnstd(section(:,:,sec).*abs(dispvec)*lro*10/XI,f,thredp);
           end
        end
    end




t=[1:nti]';
notation=['Infart '; 'Bor-Sep'; 'Bor-Lat'; 'Remote ';'iron   '; 'Average'];

figure,
subplot(1,4,1),plot(t,Ecc_reg(1,:),'-^',t,Ecc_reg(2,:),'->',t,Ecc_reg(3,:),'-<',t,Ecc_reg(4,:),'-v',t,Ecc_reg(5,:),'-*',t,Ecc_reg(6,:),'-'),title('Ecc','fontsize',18)
title([strrep(pss,'nt4_','') 'Ecc'],'fontsize',18); xlim([1 nti]);ylim([-0.5 0.05]);h = legend(notation,1);
subplot(1,4,2),plot(t,Err_reg(1,:),'-^',t,Err_reg(2,:),'->',t,Err_reg(3,:),'-<',t,Err_reg(4,:),'-v',t,Err_reg(5,:),'-*',t,Err_reg(6,:),'-'),title('Err','fontsize',18)
title([strrep(pss,'nt4_','') 'Err triangle'],'fontsize',18);  xlim([1 nti]);ylim([-0.05 0.5]);h = legend(notation,2);
subplot(1,4,3),plot(t,E2_reg(1,:),'-^',t,E2_reg(2,:),'->',t,E2_reg(3,:),'-<',t,E2_reg(4,:),'-v',t,E2_reg(5,:),'-*',t,E2_reg(6,:),'-'),title('E2','fontsize',18)
title([strrep(pss,'nt4_','') 'E2'],'fontsize',18); xlim([1 nti]);ylim([-0.5 0.05]);h = legend(notation,1);
subplot(1,4,4),plot(t,E1_reg(1,:),'-^',t,E1_reg(2,:),'->',t,E1_reg(3,:),'-<',t,E1_reg(4,:),'-v',t,E1_reg(5,:),'-*',t,E1_reg(6,:),'-'),title('E1','fontsize',18)
title([strrep(pss,'nt4_','') 'E1 triangle'],'fontsize',18);  xlim([1 nti]);ylim([-0.05 0.5]);h = legend(notation,2);


xlsname = 'C:\Varian\hupc\CardiacDense\result\CM_inf_strain2.xlsx';
sheetname = folder(strfind(folder,'nrat'):strfind(folder,'nrat')+5);
xlswrite(xlsname,{pss,['tri_' opt],'','Err','','','','','','',...
                            'Ecc','','','','','','',...
                            'E1','','','','','','',...
                            'E2','','','','','','',...
                            'displacement in mm'},sheetname,['A' num2str(starow)]);
xlswrite(xlsname,{'phase','Infart ','Bor-Sep','Bor-Lat','Remote ','iron','Average','',...
                          'Infart ','Bor-Sep','Bor-Lat','Remote ','iron','Average','',...
                          'Infart ','Bor-Sep','Bor-Lat','Remote ','iron','Average','',...
                          'Infart ','Bor-Sep','Bor-Lat','Remote ','iron','Average','',...
                          'Infart ','Bor-Sep','Bor-Lat','Remote ','iron','Average','',...
                  },sheetname,['A' num2str(starow+1)]);

save([folder sheetname '_' pss '_tri2.mat'],'epi','endo','yshift','rect',...
    'sxx','syy','rotflagx','rotflagy','ke','yp','xp','ironspotc','opt');
for n=1:nti
    xlswrite(xlsname,{num2str(n),...
        Err_reg(1,n),Err_reg(2,n),Err_reg(3,n),Err_reg(4,n),Err_reg(5,n),Err_reg(6,n),'',...
        Ecc_reg(1,n),Ecc_reg(2,n),Ecc_reg(3,n),Ecc_reg(4,n),Ecc_reg(5,n),Ecc_reg(6,n),'',...
        E1_reg(1,n) ,E1_reg(2,n) ,E1_reg(3,n), E1_reg(4,n), E1_reg(5,n),E1_reg(6,n),'',...
        E2_reg(1,n) ,E2_reg(2,n) ,E2_reg(3,n), E2_reg(4,n), E2_reg(5,n),E2_reg(6,n),'',...
        disp_reg(1,n),disp_reg(2,n),disp_reg(3,n),disp_reg(4,n),disp_reg(5,n),disp_reg(6,n),'',...
                     },sheetname,['A' num2str(starow+1+n)]);
end


