
% Example of code to run the flickering analysis
% This is for now mainly for confocal and epifluorescence images. For
% confocal images works very well. Maybe for epifluorescence it needs to be
% double checked. For movies taken with confocal the pixel dimension is
% extracted from the metadata.
% Could work also for the .movies - needs to be checked tough - (SERIES)

%clear all
%close all



addpath('C:\Users\vi211\Desktop\Nuclear shape fluctuation\Treatments\ATR_19_2_21\RAW_19_2_21_data')   
addpath('C:\Users\vi211\Desktop\Nuclear shape fluctuation')% put here the right path
% command to have the figures docked into the same window.
set(0,'DefaultFigureWindowStyle','docked')
warning('off','MATLAB:namelengthmaxexceeded')
name = 'ATR2.avi';
% few things to define for the analysis (put this into the for if the parameters change
% for the different files)


%zoom = 1;



%NSERIE = 1;              % if you have just one serie
Nrad=360;                   % number of points used to map the contour
%dpix = 0.1649*10^-6;% 0.164883988355282*10^-6% micron/px %micron in meters G1 Series004dx LaminGFP2 July2017       % set the value if the reader doesn't extract it from the metadata (in m)
%dpix = 0.114030363737983*10^-06;%meters for G2 Series002Lamin GFP3 July 2017 and S (I don't know its dpix)
%dpix = 0.109959392235374*10^-6;%G2Tlapses
%dpix = 0.0873615056460543*10^-6;%G2 LaminGFP2.lif - Series018
dpix = 0.140351565*10^-6; % All cells April2018
Temp=273.15+37;             % temperature of the experiment, in K
%WIDTH=40;   %G1                % width (in pixels) of the detection ring
%WIDTH = 30;%G2
%WIDTH = 60; %S
    WIDTH = 50;
plotting=true;

fluo=false;                 % false confocal, true epifluorescence

    %try
    
    movie = flickscript_nocircle_viSep2018(name, Nrad,Temp);
        %movie=flickering_fluo_mod(name,Nrad,Temp);%works Sep2018
        %movie=flickering_Davide_3([DIR,ves],Nrad,dpix,Temp,NSERIE,zoom);
        %(use this for brightfield)
   % catch
        %disp('There is no file')
       % break
    
    if movie.dpix == 0
        movie.dpix = dpix;
    end
    
    
        % initialization of the radius and center and contour finding
        
        movie.set_center('manual')
        movie.set_radius(WIDTH,'manual')
        movie.analyse_contour(true)
       
        % select frames without spikes and compute the spectrum
   %    
        close all
        ttFrames = [movie.Frames2Analyse(1):size(movie.contour_fine,1)];
        BADframes = checkCONTOUR(movie);
        GoodFrames = setdiff(ttFrames,BADframes);
        movie.get_spectrum(true,GoodFrames);
        %movie.get_spectrum(true);
      
        % fit
        
        movie.fit_vesicles_fluctuations(true);
        
        %save([name(1:end-4), '.mat']);
        
        %save([ves(1:end-4),'_flickering_WIDTH',num2str(WIDTH),'_VES_',num2str(count),'.mat'])
        
        % uncomment this to print fig in pdf (you can change the names of the files: now it gives the files a name which is ves_1, ves_2 for each file)
        %{
    set(gcf,'paperunits','centimeters');
    set(gcf,'paperposition', [1 1 11 7]); % left bottom width height
    set(gca,'ticklength',[0.02 0.01],'Linewidth',0.7)
    namefig_eps = ['ves_' num2str(FILES) '.eps'];
    namefig_pdf = ['ves_' num2str(FILES) '.pdf'];
    print(namefig_eps, '-depsc');
    system(sprintf(['ps2pdf -dEPSCrop ' namefig_eps ' ' namefig_pdf]));
        %}
        
        
        
        %%
        %%Calculate curvature and make videos
        %% now get cartesian coordinates from contour at each frame
        ang = linspace(-pi/2, 3*pi/2, size(movie.contour_fine,2));
        for i = size(movie.contour_fine,1):-1:1
            
            [xc(i,:), yc(i,:)] = pol2cart(ang, movie.contour_fine(i,:));
            xc(i,:) = xc(i,:) + movie.cen(i,1);
            yc(i,:) = yc(i,:) + movie.cen(i,2);
            

            curv(i,:) = calculate_signedcurvature(xc(i,:),yc(i,:));
            temp = ...
                calculate_signedcurvature(smooth(repmat(xc(i,:),1,3),30,'rloess'),...
                smooth(repmat(yc(i,:),1,3),30,'rloess'));
            sm_curv(i,:) = temp(size(curv,2) + 1 : 2*size(curv,2) );
        end %for
        
        %% plot for checking
        imagesc(sm_curv);
        h=colorbar
        t=get(h,'Limits');
        set(h,'Ticks',linspace(t(1),t(2),5))
        
        figure;
        ha = axes;
        ha.XLim = [min(xc(:)), max(xc(:))];
        ha.YLim = [min(yc(:)), max(yc(:))];
        ha.CLim = [-0.03 0.06];
        
        axis image;
        hold on;
        drawnow;
        
        v = VideoReader('MitosisLate_NE.avi');
      
         for i = size(movie.contour_fine,1):-1:1
            
            if exist('hp','var')
                delete(hp)
            end
            if exist('hi','var')
                delete(hi)
            end
            
            %frame = movie.video.read(i);
            frame = read(v,i);
            hi = imshow(frame);
            hold on        
            colormap('jet')

%             hp = scatter(xc(i,:), yc(i,:), 20, sm_curv(i,:), 'filled');
            hp = cline(xc(i,:), yc(i,:),[], sm_curv(i,:) );
            hp.LineWidth = 3;
            %pause(0.1);
            
         end
       
       %% to create videos
        imagesc(sm_curv);
        h=colorbar
        t=get(h,'Limits');
        set(h,'Ticks',linspace(t(1),t(2),5))
        
        figure;
        ha = axes;
        ha.XLim = [min(xc(:)), max(xc(:))];
        ha.YLim = [min(yc(:)), max(yc(:))];
        ha.CLim = [-0.03 0.06];
      
        axis image;
        hold on;
        drawnow;
        
        v = VideoReader('rot_transl240220_P15_Cell1.avi');
       
       %fl = dir('C:\Users\vi211\Desktop\Nuclear shape fluctuation\New Folder*.tif');
        
       %for fc = 1: numel(fl)
         
        for i = size(movie.contour_fine,1):-1:1
            
            if exist('hp','var')
                delete(hp)
            end
            if exist('hi','var')
                delete(hi)
            end
            
            %frame = movie.video.read(i);
            frame = read(v,i);
            hi = imshow(frame);
            hold on        
            colormap('jet')
           
            hp = cline(xc(i,:), yc(i,:),[], sm_curv(i,:) );
            hp.LineWidth = 3;%3;
            pause(0.1);
            eval(['print -dtiff Slice_' num2str(i) '.tiff'])         
            %saveas(figure(i),'prova','tiff')
        end
        
        %plottitle = fl(fc).name;
       % plottitle = plottitle(1:end-4);
       % plottitle = strcat(plottitle,'_curvature');
        
       %end 
       
        