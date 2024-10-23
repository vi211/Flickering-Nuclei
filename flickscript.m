%% FLICKERING CLASS: a Matlab algorithm for nuclei contour detection and analysis


classdef flickscript < handle
    % Class properties (variables) that can be set by the user in the
    % Matlab prompt, instead that using the class methods
    properties
        video               % contains the VideoReader object (with video infos)
        frame               % current frame
        current_frame       % ID of the current frame
        Frames2Analyse      % frames used for the analysis
        center              % center of the nucleus
        mean_radius         % mean value of the radius
        radius_limits       % limits of the ring region that contains the nuclear outline
        sigma               % parameter used to calculate radius_limits from mean_radius
        dt                  % Video Framerate
        dpix                % pixel to m conversion, in m
    end
    
    
    % Class properties (variables) that can be set only by the class methods
    properties (SetAccess = private)
        cen
        spectrum        % the spectra obtained from the contour
        results         % results of fitting
        Nrad            % Number of angular rays that maps the contour
        Temp            % temperature, in K
        ring            % contour region
        contour_fine    % logging contours (global fit) as 2D matrix
        xc              % X cartesian coordinate of the contour (for plotting)
        yc              % Y cartesian coordinate of the contour (for plotting)
        displacement
    end
    
    
    % Class methods (accessible by the users)
    methods
        
        function obj=flickscript(filename,Nrad,Temp,varargin)
            % Class constructor. It loads and initializes a
            % video file, creating the initial data structure.
            %
            % USAGE: obj=flickering(filename, Nrad, dpix, Temperature);
            %
            %   filename: the filename of the video to analyze;
            %
            %   Nrad:     number of points used to map the contour
            %             (usually, 360).
            %
            %   dpix:     pixel to meter conversion
            %             (spatial resolution of the camera).
            %
            %   Temperature: the temperature of the sample during the
            %                measurement. Required for correct normalization of the
            %                optput parameters of the flickering analysis
            
            % checking if the matlab pool of workers is active, to disable
            % some warnings.
%             if matlabpool('size')>0
%                 pctRunOnAll warning('off','MATLAB:nearlysingularmatrix')
%             else
%                 warning('off','MATLAB:nearlysingularmatrix')
%             end
            
            NSERIE=1;
            zoom=1;
            Frames2Analyse = [];
            
            if nargin == 6
                NSERIE=varargin{1};
                zoom=varargin{2};
            elseif nargin == 5 
                NSERIE=varargin{1};
            elseif nargin == 7
                NSERIE=varargin{1};
                zoom=varargin{2};
                Frames2Analyse = varargin{3};
            end
            

            obj.Nrad=Nrad;
            obj.Temp=Temp;
            obj.dpix=0;

            %reading the video file.
            
            if strfind(filename,'.movie')
                if which('moviereader');
                    obj.video=moviereader(filename);
                else
                    disp('Cannot open .movie videos');
                end
            elseif strfind(filename,'.lif')
                if which('bioreader')                  % for .lif and czi files same reader
                    obj.video=bioreader(filename);
                    obj.video.setSerie(NSERIE,zoom,0)   % this number here is the channel
                    obj.dpix=obj.video.dpix_zoom;
                else
                    disp('Cannot open .lif videos');
                end  
                
            elseif strfind(filename,'.czi')
                if which('bioreader');
                    obj.video=bioreader(filename);      % this is to read avi files with metadata 
                    obj.video.setSerie(NSERIE,zoom,0)   % this number here is the channel
                    obj.dpix=obj.video.dpix;
                else
                disp('Cannot open .avi with meta videos');
                end 
            else
                obj.video=VideoReader(filename);
            end
            
            if isempty(Frames2Analyse)==1
               Frames2Analyse = [1 obj.video.NumberOfFrames];
            end
            
            obj.Frames2Analyse = Frames2Analyse;
            
            %load first frame
            obj.frame = get_frame(obj.video,Frames2Analyse(1));
            obj.current_frame = Frames2Analyse(1);
            obj.dt=1./obj.video.FrameRate;
           
        end
        
        function obj=set_center(obj,varargin)
            % Function that sets the center of the nuclei.
            %
            % USAGE:
            %   obj.set_center;
            %   obj.set_center(method);
            %   obj.set_center(method,plotting);
            %
            % plotting: boolean value that activates plotting (default=true)
            %
            % method is a string that sets the center detection algorithm:
            %
            %   'manual':     center selected by a mouse click;
            %
            %   'basic_symm': (default) the center is desumed roughly as the center of
            %                 gravity of the grayscale images. The nuclei
            %                 has to be roughly at the center of the image;
            %
            %   'auto':       detected as the center of the largest object
            %                 in the automatically-thresholded image. Works
            %                 only with low-noise frames.
            
            method='basic_symm';
            plotting=true;
            
            if nargin==2;
                method=varargin{1};
                plotting=true;
            end
            
            if nargin==3
                method=varargin{1};
                plotting=varargin{2};
                if strcmp(method,'manual')
                    plotting=true;
                end
            end
            
            if plotting
                show_img(1,obj.frame)
            end
            if strcmp(method,'manual')
                title('Click on the center of the nucleus')
                obj.center=ginput(1);
                
            elseif strcmp(method,'basic_symm')
                xline=mean(obj.frame,1);
                yline=mean(obj.frame,2);
                x=1:length(xline);
                y=(1:length(yline))';
                obj.center=[sum(x.*xline)./sum(xline),sum(y.*yline)./sum(yline)];
                
            elseif strcmp(method,'auto')

                med = median(obj.frame(:));
                mad = median(abs(med -obj.frame(:)));
                m = max(0, med-6*mad);%min(obj.frame(:));
                M = med+6*mad;%M = max(obj.frame(:));
                normIM = (obj.frame-m)/(M-m+1);
                IMthres = graythresh(normIM);

                bw = im2bw(normIM,IMthres);
                
                bw=bwmorph(bw,'open',3);
                bw=bwmorph(bw,'clean');
                bw=bwmorph(bw,'close',3);
                %C=regionprops(bw,'FilledArea','Centroid');
                
                bw = imfill(bw,'holes');
                D = bwdist(~bw);
                [~, maxind] = max(D(:));
                [a,b] = ind2sub(size(D), maxind);
                
                obj.center = [b,a];
            end
            
            if plotting
                hold on
                plot(obj.center(1),obj.center(2),'r.')
            end
            
        end
        
        function obj=load_frame(obj,ind)
            % USAGE obj.load_frame( iFrame );
            % load the ith frame of the video, for debugging/control purposes.
            obj.current_frame=ind;
            obj.frame=get_frame(obj.video,ind);
        end
       
        function obj=set_radius(obj,delta,varargin)
            % Function that sets the radius of the nuclei.
            % USAGE:
            %   obj.set_radius(width);
            %   obj.set_radius(width,method);
            %   obj.set_radius(width,method,plotting);
            %
            % width: width in pixels of a ring, roughly centered on the
            %        contour of the nucleus, that contains the contour itself (usually set
            %        to 70)
            %
            % plotting: boolean value that activates plotting (default=true)
            %
            % method is a string that sets the radius detection algorithm:
            %
            %   'manual':     radius selected by a mouse click;
            %
            %   'basic_symm': (default) the radius is desumed roughly as the distance
            %                 of the maximum of the grayscale image center
            %                 with respect to the center of gravity of the
            %                 grayscale images. The nuclei (contour in white)
            %                 has to be roughly at the center of the image;
            %
            %   'auto':       detected as the radius of the largest object
            %                 in the automatically-thresholded image. Works
            %                 only with low-noise frames.
            
            method='basic_symm'; 
            plotting=true;
            
            if nargin==3;
                method=varargin{1};
                plotting=true;
            end
            
            if nargin==4
                method=varargin{1};
                plotting=varargin{2};
                if strcmp(method,'manual')%manual
                    plotting=true;
                end
            end
            
            if plotting
                show_img(1,obj.frame)
                hold on
                plot(obj.center(1),obj.center(2),'r.')
            end
            
            if strcmp(method,'manual')
                title('Select a point on the nuclear outline')
                p=ginput(1);
                
                obj.mean_radius=sqrt(sum((p-obj.center).^2));
                
            elseif strcmp(method,'basic_symm')
                a=abs(mean(obj.frame(round(obj.center(1))+(-10:1:10),:),1));
                [~,II]=max(a);
                r1=abs(II-round(obj.center(1)));
                
                b=abs(mean(obj.frame(:,round(obj.center(2))+(-10:1:10)),2));
                [~,II]=max(b);
                r2=abs(II-round(obj.center(2)));
                obj.mean_radius=mean([r1,r2]);
                
            elseif strcmp(method,'auto')           
                
                %It works well for confocal images with just one object in
                %the field of view. 
                med = median(obj.frame(:));
                mad = median(abs(med -obj.frame(:)));
                m = max(0, med-6*mad);
                M = med+6*mad;
                normIM = (obj.frame-m)/(M-m+1);
                IMthres = graythresh(normIM);

                bw = im2bw(normIM,IMthres);
                
                bw=bwmorph(bw,'open',3);
                bw=bwmorph(bw,'clean');
                bw=bwmorph(bw,'close',3);
                
                
                bw = imfill(bw,'holes');
                D = bwdist(~bw);
                
                obj.mean_radius = max(D(:));
                
            end
            
            % calculating radius and limits
            
            obj.sigma=delta/2/obj.mean_radius;
            obj.radius_limits=obj.sigma*obj.mean_radius*[1,-1]+obj.mean_radius;
            
            if plotting
                title('')
                ang=(0:0.1:360)';
                r=obj.mean_radius*[cosd(ang),sind(ang)];
                r1=obj.radius_limits(1)*[cosd(ang),sind(ang)];
                r2=obj.radius_limits(2)*[cosd(ang),sind(ang)];
                
                plot(obj.center(1)+r(:,1),obj.center(2)+r(:,2),'b--')
                plot(obj.center(1)+r1(:,1),obj.center(2)+r1(:,2),'g-.')
                plot(obj.center(1)+r2(:,1),obj.center(2)+r2(:,2),'g-.')
            end
        end
        
        function plot_contour(obj,fig)
            % function to plot the contour stored into obj.xc and obj.yc
            % (for debugging purposes)
            %
            % USAGE obj.plot_contour(figure_number)
            
            show_img(fig,obj.frame);
            hold on
            title(['Frame # ',num2str(obj.current_frame)])
            plot(obj.center(1),obj.center(2),'r.')
            plot(obj.xc+obj.center(1),obj.yc+obj.center(2),'r-');
        end
        
        
        
        function obj=analyse_contour(obj,varargin)
            % Function that detects the contour of the nucleus in every
            % frame of the movie. The video is analysed frame-by-frame
            % using a fast but raw fitting procedure.
            %
            % USAGE:
            %   obj.analyse_contour_raw;
            %   obj.analyse_contour_raw(plotting);
            %   obj.analyse_contour_raw(plotting,method);
            %
            % plotting: boolean variable to activate plotting (default: true)
            %
            % method:  "profile" (default) the contour is detected by
            %           correlating each radial profile with a previously
            %           determined template.

            fluo=0;
            plotting=0;
            
            if nargin==2
                plotting=varargin{1};
            elseif nargin==3
                plotting=varargin{1};
                fluo=varargin{2};
            elseif nargin==4
                plotting=varargin{1};
                fluo=varargin{2};
                NF=varargin{3};
            end
            
            FirstFrame = obj.Frames2Analyse(1);
            LastFrame = obj.Frames2Analyse(2);
            
            NF_Analyse = numel(FirstFrame:LastFrame);
            
            contour=zeros(NF_Analyse,obj.Nrad);
            cen_t=zeros(NF_Analyse,2);

            if fluo
                ga=fspecial('gaussian',double(round(obj.video.height/12)),double(round(obj.video.height/150)));
                bckg=imfilter(obj.frame,ga,'replicate');
                obj.frame=obj.frame-bckg+max(obj.frame(:));
            end
            
            %basic initialization for checks
            
            obj.ring=get_ring_interpolation(obj.frame,obj.Nrad,obj.radius_limits,obj.center);
            cc=find_contour(obj.ring,obj.radius_limits);
            obj.update_info(obj.ring,cc);
            
            
            %plotting the results every 60 frames
            if plotting
                figure(2)
                clf
                IM_handle=imagesc(obj.frame);
                colormap(gray)
                axis image
                hold on;
                CEN_handle=plot(obj.center(1),obj.center(2),'r.');
                CON_handle=plot(obj.xc+obj.center(1),obj.yc+obj.center(2),'r-');
                TITLE_handle=title('');
                xlim(round(obj.center(1)+obj.mean_radius*[-2,2]))
                ylim(round(obj.center(2)+obj.mean_radius*[-2,2]))
                
                BTN=uicontrol('Style','pushbutton','callback','break');
                
                drawnow
            end
            
            
            disp('Wait for contour calculation');
            tic
            t0=toc;
            
            
            for nn = FirstFrame:LastFrame;
                
                
                obj.frame=get_frame(obj.video,nn);
                
                if fluo
                    bckg=imfilter(obj.frame,ga,'replicate');
                    obj.frame=obj.frame-bckg+max(obj.frame(:));
                end
                
                obj.current_frame=nn;
                
                try
                    for fooind=1:3
                        obj.ring=get_ring_interpolation(obj.frame,obj.Nrad,obj.radius_limits,obj.center);
                        contour(nn,:)=find_contour(obj.ring,obj.radius_limits);
                        obj.update_info(obj.ring,contour(nn,:));
                        
                        cen_t(nn,:)=obj.center;
                    end
                    
                    if any(isnan(contour(nn,:))) || (mean(contour(nn,:)) > max(size(obj.frame))/2)
                        contour=contour(1:nn-3,:);         
                        cen_t=cen_t(1:nn-3,:);
                        break
                    end
                    
                catch err
                    contour=contour(1:nn-1,:);
                    cen_t=cen_t(1:nn-1,:);
                    disp(' ')
                    disp('Analysis interrupted by an undesired error.')
                    disp('The results calculated so far have been saved.')
                    break
                end
                
                str=['FPS: ',num2str(1/(toc-t0),'%.1f'),'  -  ',num2str(100*nn./NF_Analyse,'%.1f'),'%% completed'];
                t0=toc;
                
                if plotting
                    set(IM_handle,'CData', obj.frame);
                    set(CEN_handle,'XData',obj.center(1),'YData',obj.center(2));
                    set(CON_handle,'XData',obj.xc+obj.center(1),'YData',obj.yc+obj.center(2));
                    set(TITLE_handle,'string',strrep(str,'%%','%'));
                    drawnow;
                    xlim(obj.center(1)+[-1.5,1.5].*obj.mean_radius)
                    ylim(obj.center(2)+[-1.5,1.5].*obj.mean_radius)
                else
                    
                    if mod(nn,500)
                        fprintf(str);
                        Lstr=length(str)-1;
                    end
                end
            end
            
            obj.frame=get_frame(obj.video,size(contour,1));
            obj.current_frame=size(contour,1);
            
            disp(' ');
            disp('Done!')
            
            
            
            
            obj.contour_fine=contour;
            obj.cen=cen_t;
            if plotting
                figure(11)
                imagesc(obj.contour_fine')
                axis image
                xlabel('frame #')
                ylabel('azimuth angle #')
            end
        end
        

        function obj=get_spectrum(obj,varargin)
            %function that analyse the contour_fine property and extract
            %the relevant physical parameters of the membrane, and it stores
            %them in the "results" property.
            %
            % USAGE:
            %   obj.get_spectrum;
            %   obj.get_spectrum(plotting);
            %   obj.get_spectrum(plotting,frame_to_analyze);
            %
            % plotting: boolean value to trigger plotting (default=true)
            %
            % frame_to_analyze: array containing the frames that have to be
            %                   considered into this analysis. By default, the whole video is
            %                   considered. This parameter is useful if some contour shows
            %                   defects.
            
            plotting=false;
            cc=obj.contour_fine;
            frame_range = obj.Frames2Analyse(1):size(cc,1);
            
            if nargin==2
                plotting=varargin{1};
            elseif nargin==3
                plotting=varargin{1};
                frame_range=varargin{2};            % range of frames to use to build the spectrum
            end
            
            if length(frame_range)==1
                error('only 1 frame');
            end
            
           
            
            data.T=obj.Temp;

            dpix_t=obj.dpix;
            
            if dpix_t>0.1
                dpix_t=dpix_t*1e-6;
            end
            
            dt_t=obj.dt;
            clog=cc(frame_range,:);
            rr=clog*dpix_t;
            data.R = mean(rr(:));
            
            for i = 1: 360
                r0(i) = mean(rr(:,i));
            end
            
            data.displacement = rr - r0;%in meter
            dd = data.displacement;
            
            %clog1 = clog-repmat(mean(clog,2),1,size(clog,2));
            clear cc

            if numel(clog)==0
                error('Contour matrix is empty. Run the "analyse_contour" routine!');
            end
            
            %perimeter calculation
            ang=linspace(0,2*pi,obj.Nrad)';
            [x,y]=pol2cart(repmat(ang,1,size(rr,1)),rr');
            dx = x-circshift(x,[-1,0]);
            dy = y-circshift(y,[-1,0]);
            data.L=mean(sum(sqrt(dx.^2+dy.^2),1)); % real perimeter
            %data.L=2*pi*data.R;% perimeter as 2*pi*R (mean radius)
            
           
            % STATIC power spectrum, averaged over frames.
%       
            
            FFTdq=fft(dd,[],2);
            normdq=size(dd,2).^2;
            
            FFTdqm=mean(FFTdq.*conj(FFTdq)./normdq ,1);
            data.static.ps = FFTdqm(1:(obj.Nrad/2+1))';
            data.static.q = 2*pi*obj.Nrad/(data.L)*linspace(0,1,obj.Nrad+1)';
            data.static.q = data.static.q(1:(obj.Nrad/2+1));
            
            data.static.range=data.static.q>data.static.q(6)&data.static.q<data.static.q(36);
            data.static.err=std(FFTdq.*conj(FFTdq)./normdq,1,1)';
            data.static.err=data.static.err(1:(obj.Nrad/2+1));
            
            % DYNAMIC power spectrum, averaged over modes.
            
            FFTw=fft(rr,[],1);
            normw=size(rr,1).*(2*pi/dt_t);
            FFTwm=mean(FFTw.*conj(FFTw)./normw ,2);
            data.dynamic.ps=FFTwm(1:floor(length(FFTwm)/2));
            data.dynamic.omega=2*pi/dt_t*(1:length(data.dynamic.ps))'/length(data.dynamic.ps);
            
            %plot(data.dynamic.omega)
            data.dynamic.range=data.dynamic.omega>data.dynamic.omega(2)&data.dynamic.omega<data.dynamic.omega(3)*30;
            data.dynamic.err=std( FFTw.*conj(FFTw)./normw ,0,2);
            data.dynamic.err=data.dynamic.err(1:floor(length(FFTwm)/2));
            
            % Autocorrelation of the temporal evolution of the static spectrum's modes
            plot_autocor = true; %true for debugging
            tocorr=FFTdq.*conj(FFTdq)./normdq;
            data.corr.g2=autocor( tocorr',dt_t,plot_autocor);
            data.corr.g2=data.corr.g2(2:end,:); 
            
            %fitting with a simple exponential, and obtaining the relaxation time of
            %the modes.

 Fexp=fittype('a+b*exp(-x./tau)');
            options=fitoptions( 'Method','nonLinearLeastSquares',...
            'startpoint',[0,1,0.1],'lower',[-0.8,0.01,1e-3],'upper',[0.8,1.2,100]);
            
    
        
            qrange=(2:20);
            j=0;
            
            rt=data.corr.g2(:,1)<6;%keep only points tfit < 1.5
            for i = 1:length(rt)
                if rt(i) == 0
                   rt = rt(1:i-1);
                   break
                end
            end    
            data.corr.tau=zeros(length(qrange),1);
            data.corr.a=zeros(length(qrange),1);
            data.corr.b=zeros(length(qrange),1);
            data.corr.fit_g2=cell(length(qrange),1);
             figure(307)
            
            for ii=qrange
                j=j+1;
                tfit=data.corr.g2(rt,1);
                g2fit=data.corr.g2(rt,ii+1);
%                 
%                 % make it go between 1 and 0
                g2fit = g2fit - mean(g2fit(end-2:end));
                g2fit = g2fit ./ g2fit(1);
                
                HV=(max(g2fit)+min(g2fit))/2;
%                 
                [~,ind]=min(abs(g2fit-HV));
%                 %options.StartPoint(2) = g2fit(1)-options.StartPoint(1);
%                 %options.Upper(2) = 2*options.StartPoint(2);
                options.StartPoint(3)=tfit(ind);
                options.Upper(3)=tfit(ind)*2;
%                 
                [data.corr.fit_g2{j}, data.corr.gof_g2{j}]=fit(tfit,g2fit,Fexp,options);
                data.corr.tau(j)=data.corr.fit_g2{j}.tau;
                data.corr.a(j)=data.corr.fit_g2{j}.a;
                data.corr.b(j)=data.corr.fit_g2{j}.b;
                data.corr.fit_g2_score(j) = sum((data.corr.fit_g2{j}(tfit)-g2fit).^2);
                %check 
                cla
                semilogx(tfit, g2fit, 'o','MarkerSize',10)
                hold on
                semilogx(tfit,data.corr.fit_g2{j}(tfit),'g-')
                plot(tfit, g2fit, 'o')
                hold on
                plot(tfit,data.corr.fit_g2{j}(tfit),'g-')
                set(gca, 'XScale', 'linear')
                
                shg
                pause
            end
            
            data.corr.err=linspace(0.01,0.1,length(qrange))'.*data.corr.tau;
            data.corr.q=data.static.q(qrange);
            data.corr.range=(3:9);% to change the range 
            
            obj.spectrum=data;
            if plotting
                warning('OFF','MATLAB:Axes:NegativeDataInLogAxis');
                
                figure(3)
                clf
                errorbar(data.static.q(2:end),data.static.ps(2:end),data.static.err(2:end),'ro')
                set(gca,'yscale','log','xscale','log')
                xlabel('q (m^{-1})')
                ylabel('<|u(q_x)|^2> (m^2)')
                
                {
                figure(4)
                clf
                errorbar(data.dynamic.omega(3:end),data.dynamic.ps(3:end),data.dynamic.err(3:end),'ro')
                
                xlabel('\omega (rad/s)')
                ylabel('PSD (m^2/s)')
                
                figure(5)
                clf
                plot(data.corr.q,data.corr.tau,'ro')
                set(gca,'yscale','log','xscale','log')
                xlabel('q (m^{-1})')
                ylabel('\tau (s)')
                }
%                 
                drawnow
            end
        end
        
        function obj=fit_nuclei_fluctuations(obj,varargin)
            plotting=false;
            if nargin==2
                plotting=varargin{1};
            end
            data=obj.spectrum;
            
            % fitting the static spectrum
            start=[-6,-19];
            lower=[-13,-40];
            upper=[-2,-16];
           
            modelPSDq=fittype(@(s,k,x) log10(PSDq(s,k,0,data.T,data.L,x)));
            
            options=fitoptions('Method','nonlinearleastsquares',...
               'StartPoint',start,'lower',lower,'upper',upper);

            [fitPSDq,Gs]=fit(data.static.q(data.static.range),...
                log10(data.static.ps(data.static.range)),...
                modelPSDq,options);
            confidence_intervals=confint(fitPSDq,0.68);
            
            data.k=10.^fitPSDq.k;
            data.dk=mean(abs(10.^confidence_intervals(:,2)-data.k));
            data.s=10.^fitPSDq.s;
            data.ds=mean(abs(10.^confidence_intervals(:,1)-data.s));
            data.gamma=0;
            
            %             % fitting the dynamic spectrum
                        start=[fitPSDq.s,fitPSDq.k,-2];
                        lower=[fitPSDq.s-0.5,fitPSDq.k-0.5,-3];
                        upper=[fitPSDq.s+0.5,fitPSDq.k+0.5,-1];
                        modelPSDw=fittype(@(s,k,eta,x) log10(PSDw(s,k,eta,data.T,data.R,x)));
                        options=fitoptions('Method','nonlinearleastsquares',...
                                           'StartPoint',start,'lower',lower,'upper',upper,...
                                           'Weights',data.dynamic.ps(data.dynamic.range)./(data.dynamic.err(data.dynamic.range)));
                        [fitPSDw,Gw]=fit(data.dynamic.omega(data.dynamic.range),log10(data.dynamic.ps(data.dynamic.range)),...
                                  modelPSDw,options);
                        kd=10.^fitPSDw.k;
                        sd=10.^fitPSDw.s;
           
            
            

            modeltauq=fittype(@(eta,x) tauq(eta,log10(data.s),log10(data.k),0,data.R,x));
            options=fitoptions('Method','nonlinearleastsquares');
            options.Weights = 1./data.corr.fit_g2_score(data.corr.range);
            [fittau,Gt]=fit(data.corr.q(data.corr.range),data.corr.tau(data.corr.range),...
                modeltauq,options);
            data.eta=10.^fittau.eta;
            data.fit_goodness=Gs.sse+Gt.sse;
            data.fit_rsquare_static=Gs.rsquare;
            data.fit_rsquare_corr=Gt.rsquare;
            ms=PSDq(log10(data.s),log10(data.k),data.gamma,data.T,data.L,data.static.q);
            data.dev_from_model=sum(abs(data.static.ps(2:end)-ms(2:end))./ms(2:end))./length(ms);
            obj.results=data;
            
            % PLOTTING THE RESULTS % ######################################
            if plotting
                warning('OFF','MATLAB:Axes:NegativeDataInLogAxis');
                r=find(data.static.range);
                qs=logspace(log10(data.static.q(r(1))),log10(data.static.q(r(end))),1000)';
                modelstatic=PSDq(log10(data.s),log10(data.k),data.gamma,data.T,data.L,qs);
                
                om=logspace(log10(data.dynamic.omega(1)),log10(data.dynamic.omega(end)),1000)';
                modeldynamic=PSDw(log10(data.s),log10(data.k),log10(data.eta),data.T,data.R,om);
                
                qt=logspace(log10(data.corr.q(2)),log10(data.corr.q(8)),1000)';
                modelcorr=tauq(log10(data.eta),log10(data.s),log10(data.k),data.gamma,data.R,qt);
                
                figure(6)
                clf
                rr=find(data.static.range);
                err = data.static.err(rr)./sqrt(obj.dt*numel(find(obj.contour_fine(:,1),1,'first'):obj.current_frame));
                errorbar(data.static.q(rr),data.static.ps(rr),err,'ro')
                hold on
                plot(data.static.q,data.static.ps,'b.-')
                hold on
                plot(qs,modelstatic,'k-','linewidth',2)
                set(gca,'yscale','log','xscale','log')
                xlabel('q (m^{-1})')
                ylabel('<|u(q_x)|^2> (m^2)')
                
                
                              
                                figure(7)
                                clf
                                errorbar(data.dynamic.omega(3:end),data.dynamic.ps(3:end),data.dynamic.err(3:end),'r.')
                                hold on
                                plot(om,modeldynamic,'k-','linewidth',2)
                                set(gca,'yscale','log','xscale','log')
                                xlabel('\omega (rad/s)')
                                ylabel('PSD (m^2/s)')
                
                                figure(8)
                                clf
                                plot(data.corr.q,data.corr.tau,'ro')
                                hold on
                                plot(qt,modelcorr,'k-','linewidth',2)
                                set(gca,'xscale','log','yscale','log')
                
                                xlabel('q (m^{-1})')
                                ylabel('\tau (s)')
                
                                drawnow
                
            end
        end
        
      
        function obj = setT(obj,NewT)
            
            obj.Temp = NewT;
    
        end
        
    end

        
        
 
    
    methods (Hidden)
        % Hidden function that uptades the structure after each
        
        function obj=update_info(obj,ring,contour)
            obj.ring=ring;
            ang=linspace(-pi/2,3/2*pi,length(contour));
            [obj.xc,obj.yc]=pol2cart(ang,contour);
            
            obj.mean_radius=mean(contour);
            obj.radius_limits=obj.sigma*obj.mean_radius*[1,-1]+obj.mean_radius;
            obj.center=obj.center+[mean(obj.xc),mean(obj.yc)];
            
        end
    end
    
    
    
end

%% BASIC CONTOUR FITTING (PARABOLIC LOCAL FIT)
function frame=get_frame(movie,index)
%loading a frame from the movie object
frame=double(read(movie,index));
if length(size(frame))>2
    frame=squeeze(double(frame(:,:,2)));
end
end

function show_img(fig,img)
% a basic function to plot the image with the correct map and axis
if fig>0
    figure(fig)
    clf
    
end
imagesc(img)
colormap gray
axis image off
end

function II=clean_cont(II,lim,lmax)
% NB=smooth(abs(gradient(II)),3)'>2;
% II(NB)=round(mean(II(~NB)));
%
%
% PB_low=II>=lim+2;
% PB_high=II<=(lmax-(lim+2));
%
% II(~PB_low)=lim+2;
% II(~PB_high)=lmax-(lim+2);

for k=1:2
    NB=(abs(gradient(II))>2);
    NB=NB+(abs(gradient(gradient(II)))>1);
    NB=NB+(abs(gradient(gradient(gradient(II))))>0.5);
    NB=NB>0;
    
    II(NB)=ceil(nanmean(II(~NB)));
end
end

function [vx,vy,a]=lsq_parab_vertex(x,y)
%least square parabolic fitting
s0 = length(x); s1 = sum(x); s2 = sum(x.^2); s3 = sum(x.^3); s4 = sum(x.^4);
A = [s4,s3,s2;s3,s2,s1;s2,s1,s0];
d = [sum(x.^2.*y);sum(x.*y);sum(y)];
a = A\d;

vx=-a(2)./2./a(1);
vy = a(1)*vx.^2+a(2)*vx+a(3);
end


function ring=get_ring_interpolation(frame,Nrad,radius_limits,center)
%generating a coordinate array of the ring that contains the
%membrane
ang=linspace(-pi/2,3/2*pi,Nrad);
Np=round(abs(diff(radius_limits)));
rad=linspace(radius_limits(1),radius_limits(2),Np);
[theta,rho]=meshgrid(ang,rad);
[Xp,Yp]=pol2cart(theta,rho);

%generating an interpolant function of the whole frame
[X,Y]=meshgrid(1:size(frame,2),1:size(frame,1));
F=griddedInterpolant(X',Y',frame','cubic');

%interpolating (30 times faster than multiple line interpolation);
ring=F(Xp+center(1),Yp+center(2));

end



function contour=find_contour(rr,radius_limits)

rr(isnan(rr))=nanmean(nanmean(rr));

%{
% p=zeros(1,size(rr,2));
% x=(0:size(rr,1)-1)';
% for k=1:size(rr,2);
%     y=bsxfun(@minus,rr(:,k),nanmin(rr(:,k)));
%     C0=round(nansum(x.*y)./nansum(y));
% %     p(k)=lsq_parab_vertex(C0+(-5:5)',y(C0+(-5:5)));
%     p(k)=nansum(x(C0+(-3:3)).*y(C0+(-3:3)))./nansum(y(C0+(-3:3)));
% end
% isnan_parab=isnan(mean(p,1));
% p(:,isnan_parab)=repmat(nanmean(p,2),1,sum(isnan_parab));
%
% contour=radius_limits(1)-p-0.5;
% contour=clean_cont(contour,5,size(rr,1));
%}

 rr=bsxfun(@minus,rr,min(rr,[],1));          % this subtracts the background
 rr=bsxfun(@rdivide,rr,max(rr,[],1));        % it normalizes

p=zeros(size(rr,2),1);
v=zeros(size(rr,2),1);

for nn=1:size(rr,2);
    
    y=rr(:,nn);
    x=(1:length(y))';
    
    [~,p(nn)] = max(y);                             % just the maximum works fine for confocal images
    %centroid_ind = round(sum(x.*y)./sum(y));
    %[~,p(nn)] = max(y(1:centroid_ind));             % these two lines are for epifluorescence images
    
    
    v(nn)=3;                                        % this number sets the dimension of the window on which then the centroid is calculated
    
%     CCColor = hsv(5);    
%     if nn==4
%         figure
%         plot(x,y,'.-g')
       
        for iter = 1:5
%             hold on
%             plot(repmat(p(nn),1,2),[min(y),max(y)],'Color',CCColor(iter,:))
%             plot(repmat(max(p(nn)-v(nn),min(x)),1,2),[min(y),max(y)],'Color',CCColor(iter,:))
%             plot(repmat(min(p(nn)+v(nn),max(x)),1,2),[min(y),max(y)],'Color',CCColor(iter,:))
            indxx = x<=min(p(nn)+v(nn),max(x)) & x>=max(p(nn)-v(nn),min(x)); 
%             plot(x(indxx),ynI(indxx),'*b')
            p(nn)=sum(x(indxx).*y(indxx))./sum(y(indxx));
            %v(nn)=sqrt(nansum((x-p(nn)).^2.*y)./nansum(y));
        end
        
%         pause
%         close
%     end
    
end
contour=radius_limits(1)-p'+0.5;

end


%% DATA ANALYSIS ##########################################

function chisq = globalfit(par,data)
k=par(1);
kd=par(2);
s=par(3);
gamma=par(4);
eta=par(5);

modelstatic=PSDq(s,k,gamma,data.T,data.R,data.static.q(data.static.range));
modeldynamic=PSDw(s,kd,eta,data.T,data.R,data.dynamic.omega(data.dynamic.range));
modelcorr=tauq(eta,s,kd,gamma,data.R,data.corr.q(data.corr.range));

chiq=(log10(data.static.ps(data.static.range))-log10(modelstatic))./length(modelstatic);
chiw=(log10(data.dynamic.ps(data.dynamic.range))-log10(modeldynamic))./length(modeldynamic);
chit=(log10(data.corr.tau(data.corr.range))-log10(modelcorr))./length(modelcorr);

chiq=chiq./(data.static.err(data.static.range)./data.static.ps(data.static.range));
chiw=chiw./(data.dynamic.err(data.dynamic.range)./data.dynamic.ps(data.dynamic.range));
chit=chit./(data.corr.err(data.corr.range)./data.corr.tau(data.corr.range));

chisq=[chiq;chiw;chit];

end

function model=tauq(eta,s,k,g,R,qm)
eta=10.^eta;
s=10.^s;
k=10.^k;
etaM=1e-9;
etaext=1e-3;
model=0.8*2*(etaM/R^2+qm*(eta+etaext))./(2*g+s*qm.^2+k*qm.^4);
end

function spectrum=PSDq(s,k,g,T,L,x)
s=10.^s;
k=10.^k;
kb=1.38e-23; % SI
spectrum=kb*T/(L)*1/2/s*( 1./x - 1./sqrt(s/k+x.^2));%eq 1 Yoon2008
end

function spectrum=PSDw(s,k,eta,T,R,omega)
s=10.^s;
k=10.^k;
eta=10.^eta;

kb=1.38e-23; % SI

spectrum=zeros(size(omega));
for l=2:130
    vl=-l:l;
    Zlm=Zl(vl);
    sZl=sum(Zlm(~isinf(Zlm)));
    foo=kb.*T./R./eta./sZl.*(2*l+1)/2/pi./(wl(s,k,eta,R,l).^2+omega.^2);
    spectrum=spectrum+foo;
end
end

function wl=wl(s,k,eta,R,l)
wl=(k*(l+2).*(l-1).*l.*(l+1)+ s*R^2*(l+2).*(l-1))./(eta*R^3*Zl(l));
end

function Zl=Zl(l)
Zl=(2*l+1).*(2*l.^2+2.*l-1)./(l.*(l+1));
end


%% AUTOCOR ALGORITHM (multi-tau, vectorized algorithm)
function g2=autocor(I,dwell_time,show)
%	Function autocorrelate the data using multiple tau method.

global buf G num cts cur nobuf nolev

I=I./mean(I(:));

nobuf=8;  %must be even!
nolev=8;

timedelay=delays(dwell_time);
%initialize all the arrays
Ng2=size(I,1);
buf=zeros(nolev,nobuf,Ng2); %matrix of buffers
cts=zeros(nolev,1);
cur=nobuf*ones(nolev,1);
G=zeros((nolev+1)*nobuf/2,Ng2);
num=zeros(nolev,1);


if show
    figure(1)
    clf
    h=loglog(timedelay,G(:,1),'b.-');
    xlim([timedelay(2),timedelay(end)])
end
tic;
L=size(I,2);
for n=1:L
    insertimg(I(:,n))
    %timer,el1
    %et=el1-el
    if show
        fprintf('.:: processed frame: %d(%d)  - ',n,L)
        fprintf(' elapsed time: %.4f \r',toc);
    end
    if show && (mod(n,30)==0)
        set(h,'Ydata',G(:,1),'Xdata',timedelay);
        pause(0.01)
    end
    
end
if show
    disp('Done!');
    disp(' ');
end
norm=repmat(mean(I,2).^2,1,length(timedelay))';
g2=[timedelay,G./norm];
end

function dly=delays(time)
%    return array of delays.
%    KEYWORD:  time: scale delays by time ( should be time between frames)

global nolev nobuf

dly=zeros((nolev+1)*nobuf/2,1);
for i=1:nolev
    if i==1
        imin=2;
    else
        imin=nobuf/2+1;
    end
    ptr=(i-1)*nobuf/2+(imin:nobuf);
    dly(ptr)= ((imin-1):(nobuf-1)).*2^(i-1);
end
dly=dly*time;
end

function insertimg(In)
%   read and insert image, n, into accumulator buffers,
%   then increments appropriate correlators for new image.
global buf cts cur nobuf nolev


cur(1)=1+mod(cur(1),nobuf);  %increment buffer

buf(1,cur(1),:)=In;

process(1,cur(1));

processing=1;
lev=2;
while processing
    % either add previous two images then process or set flag next one
    if(cts(lev))
        prev=1+mod((cur(lev-1)-1-1+nobuf),nobuf);
        cur(lev)=1+mod(cur(lev),nobuf);
        buf(lev,cur(lev),:) = (buf(lev-1,prev,:)+buf(lev-1,cur(lev-1),:))./2;
        cts(lev)=0;
        process(lev,cur(lev));
        lev=lev+1;
        %Since this level finished, test if there is a next level for processing
        processing = (lev<=nolev);
    else
        cts(lev)=1; % set flag to process next time
        processing=0; % can stop until more images are accumulated
    end
end
end

function process(lev,bufno)
%	The buffer bufno at level lev has just been added so update
%	the correlation and intensity averages it affects
global buf G num nobuf


num(lev)=num(lev)+1; % one more in average
if (lev==1)
    imin=1;
else
    imin=1+nobuf/2;% no need for lower half of level > 1
end

for i=imin:min([num(lev),nobuf])%loop over delays
    ptr=(lev-1)*nobuf/2+i;
    delayno=1+mod((bufno-(i-1)-1+nobuf),nobuf); %cyclic buffers
    IP=squeeze(buf(lev,delayno,:));
    IF=squeeze(buf(lev,bufno,:));
    foo=(IF.*IP)';
    G(ptr,:) = G(ptr,:)+(foo-G(ptr,:))/(num(lev)-i+1);
end

end


