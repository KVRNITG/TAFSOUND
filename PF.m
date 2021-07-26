 close all;
 clear all;
tic;


% ues data
% uex = x coordinate of the ue
% uey = y coordinate of the ue
% sbx = x coordinate of the s-BS
% sby = y coordinate of the s-BS
% sinrt;   % predefind SINR threshold value
% tp= transmission power of UE
%alpha=attenuation factor
%delta= thermal noise of the UE
% smax= mxmimum number ues supported by bs 
%pmax= Maximum transmission power of the s-BS
%pa= Circuit power of active s-BS
%ps= Circuit power of sleep s-BS
 %TR= Transmission range of s-bs
%ues= Read the coordinates of the UEs 
%SEE= small cell energy efficiency
uex=ues(:,1);  %intailzation points
 uey=ues(:,2); 
sue=size(ues,1);
%sbs=Read the coordinates of the small base station
ssb=size(sbs,1);
 sbsx=sbs(:,1);
 sbsy=sbs(:,2);
 
 

 Depysbs=1;
 srtpm=0;
 uel=0;
 while(Depysbs==1)
     [ssb,bsleast,sbsx,sbsy,Alluesconn,uelr,Depysbs]=connb(uex,uey,sbsx,sbsy,ssb,sue,srtpm,uel); % establish connected network 
   
     srtpk=bsleast;
     if(Alluesconn==1)
   uel=[sbsx,sbsy];
     else
        uel=uelr;
     end
     [ssb,xss,yss]=rbase(bsleast,sbsx,sbsy,ssb); % sleep the small base station (s-BS)
sbsx=0;sbsy=0;
sbsx=xss(:);
sbsy=yss(:);

 end
  toc;
 function [ssb,xss,yss]=rbase(bsleast,sbsx,sbsy,ssb)

xss=0;yss=0;
q=0;
 for l=1:ssb     
     if(l~=bsleast)
         q=q+1;
         xss(q)=sbsx(l);
         yss(q)=sbsy(l);
         
      
     end
 end
 ssb=ssb-1;
end
 function [ssb,bsleast,sbsx,sbsy,Alluesconn,uelr,yes]=connb(uex,uey,sbsx,sbsy,ssb,sue,srtpm,uel)

for i=1:ssb
    for j=1:sue        
        yDiff=sbsy(i)-uey(j);
        ySqr=yDiff*yDiff;
        xDiff=sbsx(i)-uex(j);
        xSqr=xDiff*xDiff;
        output(i,j)=sqrt(ySqr + xSqr);
        snr(i,j)=((tp*((output(i,j)^-alpha))/delta); %SINR caculation    
        
    end

end
 Assign=zeros(4,sue);
bsdegree=0;
for i=1:ssb r=0;
    for j=1:sue
       m=0;
        if(snr(i,j)>sinrt && Assign(1,j)==0) % ue not connectted any bs
            Assign(1,j)=1; %  UE j conneted to s-BS i
            Assign(2,j)=i; % with which bs
            Assign(3,j)=snr(i,j); % how much sinr between the ue to s-BS 
r=r+1;% count of each Bs connted to how many UEs
m=m+1;
        else if(snr(i,j)>sinrt&&Assign(1,j)==1)
                
                if(Assign(3,j)<snr(i,j))
                  kk=Assign(2,j); % find the BS which is minimum snr
                   bsdegree(kk)=bsdegree(kk)-1;
                    Assign(2,j)=i;
                    Assign(3,j)=snr(i,j);
                    Assign(4,j)=i;
                   if(m==0)
                    r=r+1 ;
                   end
                end
                
            end
        end 
    end
    bsdegree(i)=r;
end
Alluesconn=1;
for i=1:sue 
for j=1:ssb
    if(Assign(1,j)==0) 
        Alluesconn=0; % if any not connected to the s-BS
        srtp=srtpm; % back the ues and s-BS
        uelr=uel;
        break;
    end
end
end

if(Alluesconn==1)
    uelr=[sbsx,sbsy];
for i=1:ssb csc=1; snrst=0;
    
    while(bsdegree(i)>smax && bsdegree(i)>=csc)
        
        for j=1:sue
            if( Assign(2,j)==i)
                snrst(1,csc)=Assign(3,j); % know the ues which s-BS having greater smax ues
                snrst(2,csc)=j;
                csc=csc+1;
            end
        end
    end
 

if(bsdegree(i)>smax)
    
[temp,order]=sort(snrst(1,:)); % sort the associated ues of the s-bs and respective snr
snrst = snrst(:,order);
ueremove=bsdegree(i)-smax;
for k=1:ueremove
    % if the size is more than smax
         [a,b]=sort(snr(:,snrst(2,k)),'descend'); % first removed ue to all s-BSs snr prefernce is sorted  
         ueprefernce=[a,b];
         br=0;
   if(ueprefernce(1)==snrst(1,k)& br==0)
           uep=2;
       if(bsdegree(ueprefernce(2,uep))<smax)
          Assign(2,snrst(2,k))=ueprefernce(2,uep); 
          
        Assign(3,snrst(2,k))=ueprefernce(uep,1);
        br=1;
        bsdegree(ueprefernce(2,uep))=bsdegree(ueprefernce(2,uep))+1; % update degree of the ue assocaited s-BS
        bsdegree(i)=bsdegree(i)-1;  % update degree of the removed  ue assocaited s-BS
       else
           uep=uep+1;
       end
      
   end
end
end
   
end
end
for i=1:ssb    
    bsloads(i)=bsdegree(i)/smax;     % calculate the each s-BS load
    if(bsloads(i)>0)
        bsactive(i)=1;
    else
         bsactive(i)=0;
    end
end
[bsleastper,bsload]=sort(bsloads);
bsleast=bsload(1);
bsleastper=bsleastper(1);
bsleshred=0;
while(bsleshred==0)

count=0;
for i=1:ssb  
       yDiff=sbsy(bsleast)-sbsy(i);
        ySqr=yDiff*yDiff;
        xDiff=sbsx(bsleast)-sbsx(i);
        xSqr=xDiff*xDiff;
        bstbsd(i)=sqrt(ySqr + xSqr);
        if (bstbsd(i)<TR& bstbsd(i)~=0)
            count=count+1;
            nblist(1,count)=i; % 
            nblist(2,count)=bsloads(i);
            
        end
end
try
    Depysbs=1;
[temp,order]=sort(nblist(2,:));
snlist=nblist(:,order);
[temp,order]=sort(nblist(2,:),'descend');
sllist=nblist(:,order);
j=0;
fixs=(bsleastper/size(snlist,2));
flag=0;
for i=1:size(nblist,2)   
    
    snloads(1,i)= snlist(1,i); 
    snloads(2,i)=(sllist(2,i)/sum(sllist(2,:)))*bsleastper;
    
    if(snloads(2,i)>fixs & flag<bsdegree(bsleast))
        j=j+1;
        snloadsm(1,j)=snloads(1,i); % particular BS
        snloadsm(2,j)=snloads(2,i); % percentage of BS
        snloadsm(3,j)=ceil(snloadsm(2,j)/(1/smax)); % how many ues are possible to connect
        snloadsm(4,j)=0;
        flag=flag+snloadsm(3,j);
    end
    
    
end
while(flag~=bsdegree(bsleast))
      [a,b]=sort(snloadsm(3,:))
  snloadsm(3,b(1))=snloadsm(3,b(1))+1;
  flag=flag+1;
    end
bsleastcues=bsdegree(bsleast);
j=0;
for i=1:sue   
    if(Assign(2,i)==bsleast)
    j=j+1;
    bsleastcue(1,j)=i;
    bsleastcue(2,j)=0;
    end    
end
for i=1:size(snloadsm,2)  
    for j=1:bsleastcues
        if((snr(snloadsm(1,i),bsleastcue(1,j)))>sinrt)
    bsleastcuesnr(i,j)=snr(snloadsm(1,i),bsleastcue(1,j));
        else
         bsleastcuesnr(i,j)=0;
        end
    end    
end

for i=1:size(snloadsm,2)
    k=0;
  
      [a,b]=sort(bsleastcuesnr(i,:),'descend')
    while(snloadsm(4,i)~=1) % share the load of least loaded s-BS to neighbour s-BS based on the percentage
        
    if(bsdegree(snloadsm(1,i)<=smax))
      for j=1:size(bsleastcue,2)
         if((snloadsm(4,i)==0)&(bsleastcue(2,b(j))==0))
       bsleastcue(2,b(j))=1;
       degree=Assign(2,bsleastcue(1,b(j)));
       Assign(2,bsleastcue(1,b(j)))=snloadsm(1,i);
       Assign(3,bsleastcue(1,b(j)))=a(j);
       k=k+1; 
       bsdegree(snloadsm(1,i))=bsdegree(snloadsm(1,i))+1;
       bsdegree(degree)=bsdegree(degree)-1;
       if(snloadsm(3,i)==k)
           snloadsm(4,i)=1;
       end   
          
      end
      end
    end
    end
end
 
bsdegree(bsleast)=0;
bsleshred=1;
bsactive(bsleast)=0;

%%EE
%pa=power consumption of active s-BS;%
% ps=power consumption of active s-BS;

SNT=0;
NC=0;
ussb=ssb-1;
for c=1:ssb
    R_s=0;
     F_c=0;
      % upper part eqation of EE
    for v=1:sue
        if(Assign(2,v)==c)
  
    R(c,v)=(log2(1+(Assign(3,v)))*bsactive(c));
   R_s=R_s+R(c,v);
    F_c=F_c+1;
        end
        if(v==sue)
            R_ss(c)=R_s;
        end
    end
    T_s(c)=(R_ss(c))*(1);
    SNT=T_s(c)+SNT;
   
NC=NC+(bsactive(c)*(((pmax*F_c)/sm)+pa))+((1-bsactive(c))*ps);    

    if(c==ssb)
            SNC=NC;
        end
end
SEE=(SNT/SNC);
NoofAsbs=size(sbsx,2);

catch exception
    Depysbs=0;
break;
end


end
 end