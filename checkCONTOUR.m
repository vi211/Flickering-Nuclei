

% check the contour: draw the contour to see if there are spikes etc...

function BADframes = checkCONTOUR(movie) 

BADframes = [];

for i = movie.Frames2Analyse(1):size(movie.contour_fine,1)
    
%       ang=linspace(0,2*pi,length(movie.contour_fine(i,:)));     
%       [xc,yc]=pol2cart(ang,movie.contour_fine(i,:));
    
    DIFF = diff(movie.contour_fine(i,:));
    CUT_OFF = DIFF>4;
    
    if i == movie.Frames2Analyse(1); R_1 = mean(movie.contour_fine(i,:));end 
    
    Ri = mean(movie.contour_fine(i,:)); 
   
    if sum(CUT_OFF)<1 && abs(R_1-Ri)<10 && nnz(DIFF)>200
%         figure(1)
%         plot(xc,yc,'.-b')
%         axis equal
%         xlim([-size(movie.frame,2)/2 size(movie.frame,2)/2])
%         ylim([-size(movie.frame,1)/2 size(movie.frame,1)/2])
%         pause(0.1)
        
    else
%         figure(1)
%         plot(xc,yc,'.-r')
%         axis equal
%         xlim([-size(movie.frame,2)/2 size(movie.frame,2)/2])
%         ylim([-size(movie.frame,1)/2 size(movie.frame,1)/2])
%         pause
        BADframes = [BADframes i];
        
    end
    
%     if round(i/60)==i/60
%         frame=read(movie.video,i);
%         figure(1)
%         imshow(frame,[])
%         hold on
%         plot(xc+movie.cen(i,1),yc+movie.cen(i,2),'r.-');
%         
%         pause(0.5)
%     end
    
end

end