Program Morphodynamic_gravelbed

!******************************** Morphodynamic_gravelbed *******************************
!************************ Copyright (C) 2019  Sadegh Jafarinik **************************



!****************************************************************************************
!****************************************************************************************
!*********** This program is free software: you can redistribute it and/or modify********
!*********** it under the terms of the GNU General Public License as published by********
!*********** the Free Software Foundation, either version 3 of the License  *************
!*********** , or (at your option) any later version.                       *************
!*********** This program is distributed in the hope that it will be useful, ************
!*********** but WITHOUT ANY WARRANTY; without even the implied warranty of *************
!*********** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the **************
!*********** GNU General Public License for more details.                  **************
!*********** You should have received a copy of the GNU General Public License **********
!*********** along with this program.  If not, see <https://www.gnu.org/licenses/>. *****
!****************************************************************************************
!****************************************************************************************

!----------------define parameters------------------------

implicit none
integer :: N=93,S=7, Mprint=30, Mtoprint=1500, i, j,jj,v,k,jjj,m,o

!----------------general parameters------------------------
real*8 :: dt,r, g,dx,alph,alphr,etaao

!----------------water and sediment parameters------------------------
real*8 :: qbf_meouge, qbf, lp,dg
real*8 :: cumulative,qbft=0,qbft_meouge=0
real*8, dimension(93,15):: qw,qww,b
real*8, dimension (93,15):: dsgbin,slopebin,wstarbin=0,taustarout=0
real*8, dimension(15):: qfreq,qw_meouge
real*8, dimension(13):: qbfi,qbfi_meouge
real*8, dimension(93):: x=0,etaa=0,h=0,slope=0,qbt=0,detaa=0,dgf=0 ,u, qvol=0,cf,catcharea
!----------------Bedload relation (Wilcock and Crowe 2003)------------------------

real*8,dimension(93,13)::bb,phi,gphi,wstar
real*8,dimension(93,13):: qb=0,fa,fl=.166,pb
real*8, dimension(93):: fss=0.0,taustarssrg,taustarsg,ustar,wstarg,taustart

!----------------Grain characteristics and D90,ks,la------------------------

real*8 :: fad=0,fau=0,d90i,psibar=0,psibarf=0,d50i
real*8, dimension(93):: la,d90, d50,ks
real*8, dimension(14):: db,psib,ff
real*8, dimension(13):: d,psi,f

!----------------sediment volume transported ------------------------
real*8,dimension(93,13)::qbi2019=0,qbi2028=0, qbi2031=0, qbi2082=0,qvoli=0

!----------------OUTPUTs ------------------------
real*8, dimension(93,102) :: etaaout=0, etabout=0, qbout=0, hout=0, dgout=0,slopeout=0,qvolout=0,d90out=0,d50out=0

!----------------STRATIGRAPHY PARAMETERS------------------------
integer , dimension (93):: p
real*8:: ls
real*8, dimension (93,13,1000) :: fs
real*8, dimension (93,1000) :: dg_sub
real*8, dimension (93) :: l_buffer
real*8, dimension (93,13):: f_buffer

!----------------Define constants------------------------
alph=0.2               !Hoey and Ferguson 1994 parameter for the fraction of sediment at interface between active layer and substrate
alphr=8.1              !constant for friction coefficient relation (Manning-Strickler formulation). 
dt=1.0/500.0           !time [per year], it should be small enough that the code doesn't crash. Depends on Dx.
dt=dt*365.25*24*60*60  !time [per second]
dx= 330.0              ! distance between computational nodes [m]
qbf=1790.0             ! sediment input volume per year at the upstream node (St sauveur dam) [m3/year]
qbf_meouge=0.0     ! Sediment input volume per year from any branch (Meouge river) [m3/year]
qbf= qbf/(365.25*24*60*60)  ! sediment input volume per second at the upstream node (St sauveur dam) [m3/sec]
qbf_meouge=qbf_meouge/(365.25*24*60*60)      ! Sediment input volume per second from Meouge river [m3/sec]
lp=0.4                 ! sediment porosity
r= 1.65                ! submerged specific gravity of sediment
g=9.81                 ! gravity acceleration [m/s^2]
ls=1.0                 ! thickness of the layers in the substrate for stratigraphy calculations [m]


!----------------grain size characteristics of the bedload------------------------
!---------------------------------------------------------------------------------

!-----Grain size distribution of the sediment at upstream/ meouge and the replenishment pile downstream of the dam----

ff(1)=0.023   !fraction  
db(1)=0.08    !Size [mm]         
ff(2)=0.045   !fraction
db(2)=0.2     !Size [mm] 
ff(3)=0.105	  !fraction 
db(3)=0.5     !Size [mm] 
ff(4)=0.13    !fraction		
db(4)=1.0     !Size [mm] 
ff(5)=0.148	  !fraction	
db(5)=2.0     !Size [mm] 
ff(6)=0.183	  !fraction	
db(6)=5.0     !Size [mm] 
ff(7)=0.24    !fraction		
db(7)=10.0    !Size [mm] 
ff(8)=0.353   !fraction		
db(8)=20.0    !Size [mm] 
ff(9)=0.49    !fraction	
db(9)=31.5    !Size [mm] 
ff(10)=0.655  !fraction		
db(10)=50.0   !Size [mm] 
ff(11)=0.823  !fraction		
db(11)=80.0   !Size [mm] 
ff(12)=0.878  !fraction		
db(12)=100.0  !Size [mm] 
ff(13)=0.955  !fraction		
db(13)=120.0  !Size [mm] 
ff(14)=1.0	  !fraction	
db(14)=150.0  !Size [mm] 

  
!----------------Geometric mean diameter of the input material and initial bed condition------------------------
psibar=0  
do i=1,14
  	psib(i)=log(db(i))/log(2.0)
end do
do i=1,13
  psi(i)=.5*(psib(i)+psib(i+1))
  d(i)=(db(i)*db(i+1))**.5 
  f(i)=ff(i+1)-ff(i)
end do

do i=1,13
  psibar=psi(i)*f(i)+psibar
end do
dg=2**psibar
dg=dg/1000    !geometric mean diameter of input sediment [m]
d=d/1000      !grain size class [m]   

!----------------D90 for input material and initial bed condition-----------------------
 do v=1,13
    fau=f(v)+fau
    if (fau>=.90)then
      fad=fau-f(v)
      d90i=2.72**((log(db(v+1))-log(db(v)))*(.9-fad)/(fau-fad)+log(db(v)))/1000
      exit 
    end if
 end do

la=2*d90i    !initial active layer thickness [m]

!----------------D50 for input material and initial bed condition-----------------------
do v=1,13
    fau=f(v)+fau
    if (fau>=.50)then
      fad=fau-f(v)
      d50i=2.72**((log(db(v+1))-log(db(v)))*(.5-fad)/(fau-fad)+log(db(v)))/1000
      exit 
    end if
 end do
 

!=======================================================================================
!--------------------------------Opening the I/O files----------------------------------
!=======================================================================================

!-------------------------------------INPUTS--------------------------------------------
                          open(41,file='input-eta')         !initial bed elevation, input
                          open(67,file='fa-input')                !fraction of sediment in the surface, Input
                          open(53,file='upstream-rating_curve_pre-dam')       !upstream rating curve for pre-dam condition, input
                          open(101,file='upstream-rating_curve_post-dam')     !upstream rating curve for post-dam condition, input
                          open(123, file='meouge_rating_curve')                 !Meouge rating curve                            
                          open(54,file='catchment_area')                      !contributing catchment area, input

!-------------------------------------OUTPUTS--------------------------------------------                          
                          open(1,file='etaa1')             !bed elevation, output
                          open(13,file='etaa2')            !bed elevation, output
                          open(14,file='etaa3')            !bed elevation, output
                          open(30,file='etaa4')            !bed elevation, output
                          open(31,file='etaa5')            !bed elevation, output
                          open(37,file='etaa6')            !bed elevation, output
                          open(38,file='etaa7')            !bed elevation, output
                          open(39,file='etaa8')            !bed elevation, output
                          open(4,file='dgf1')              !geometrice mean diameter of surface sediment, output
                          open(17,file='dgf2')             !geometrice mean diameter of surface sediment, output
                          open(18,file='dgf3')             !geometrice mean diameter of surface sediment, output
                          open(43,file='dgf4')             !geometrice mean diameter of surface sediment, output
                          open(44,file='dgf5')             !geometrice mean diameter of surface sediment, output
                          open(45,file='dgf6')             !geometrice mean diameter of surface sediment, output
                          open(46,file='dgf7')             !geometrice mean diameter of surface sediment, output
                          open(47,file='dgf8')             !geometrice mean diameter of surface sediment, output  
                          open(77,file='d90-1')                   !d90 of the surface sediment, output
                          open(78,file='d90-2')                   !d90 of the surface sediment, output
                          open(79,file='d90-3')                   !d90 of the surface sediment, output
                          open(80,file='d90-4')                   !d90 of the surface sediment, output
                          open(81,file='d90-5')                   !d90 of the surface sediment, output
                          open(82,file='d90-6')                   !d90 of the surface sediment, output
                          open(83,file='d90-7')                   !d90 of the surface sediment, output
                          open(84,file='d90-8')                   !d90 of the surface sediment, output
                          open(85,file='d50-1')                   !d50 of the surface sediment, output
                          open(86,file='d50-2')                   !d50 of the surface sediment, output
                          open(87,file='d50-3')                   !d50 of the surface sediment, output
                          open(88,file='d50-4')                   !d50 of the surface sediment, output
                          open(89,file='d50-5')                   !d50 of the surface sediment, output
                          open(90,file='d50-6')                   !d50 of the surface sediment, output
                          open(91,file='d50-7')                   !d50 of the surface sediment, output
                          open(92,file='d50-8')                   !d50 of the surface sediment, output                                           
                          open(5,file='water surface1')     !water depth, output
                          open(19,file='water surface2')    !water depth, output
                          open(20,file='water surface3')    !water depth, output
                          open(32,file='water surface4')    !water depth, output
                          open(33,file='water surface5')    !water depth, output
                          open(34,file='water surface6')    !water depth, output
                          open(35,file='water surface7')    !water depth, output
                          open(36,file='water surface8')    !water depth, output
                          open(7,file='slope1')             !bed slope, output
                          open(21,file='slope2')            !bed slope, output
                          open(22,file='slope3')            !bed slope, output
                          open(102,file='slope4')           !bed slope, output
                          open(103,file='slope5')           !bed slope, output
                          open(104,file='slope6')           !bed slope, output
                          open(105,file='slope7')           !bed slope, output
                          open(106,file='slope8')           !bed slope, output                         
                          open(68,file='fa-output')               !fraction of sediment in the surface, Output
                          open(69,file='Sediment_volume_1')       !total volume of sediment passing through each node, Output
                          open(70,file='Sediment_volume_2')       !total volume of sediment passing through each node, Output
                          open(71,file='Sediment_volume_3')       !total volume of sediment passing through each node, Output
                          open(72,file='Sediment_volume_4')       !total volume of sediment passing through each node, Output
                          open(73,file='Sediment_volume_5')       !total volume of sediment passing through each node, Output
                          open(74,file='Sediment_volume_6')       !total volume of sediment passing through each node, Output
                          open(75,file='Sediment_volume_7')       !total volume of sediment passing through each node, Output
                          open(76,file='Sediment_volume_8')       !total volume of sediment passing through each node, Output
                          open(107, file='qbiv2019-1')            !grain size specific volume passing through nodes in 2019, output
                          open(108, file='qbiv2019-2')            !grain size specific volume passing through nodes in 2019, output  
                          open(109, file='qbiv2019-3')            !grain size specific volume passing through nodes in 2019, output
                          open(110, file='qbiv2019-4')            !grain size specific volume passing through nodes in 2019, output
                          open(111, file='qbiv2028-1')            !grain size specific volume passing through nodes in 2028, output
                          open(112, file='qbiv2028-2')            !grain size specific volume passing through nodes in 2028, output
                          open(113, file='qbiv2028-3')            !grain size specific volume passing through nodes in 2028, output
                          open(114, file='qbiv2028-4')            !grain size specific volume passing through nodes in 2028, output
                          open(115, file='qbiv2031-1')            !grain size specific volume passing through nodes in 2031, output
                          open(116, file='qbiv2031-2')            !grain size specific volume passing through nodes in 2031, output
                          open(117, file='qbiv2031-3')            !grain size specific volume passing through nodes in 2031, output
                          open(118, file='qbiv2031-4')            !grain size specific volume passing through nodes in 2031, output
                          open(119, file='qbiv2082-1')            !grain size specific volume passing through nodes in 2082, output
                          open(120, file='qbiv2082-2')            !grain size specific volume passing through nodes in 2082, output
                          open(121, file='qbiv2082-3')            !grain size specific volume passing through nodes in 2082, output
                          open(122, file='qbiv2082-4')            !grain size specific volume passing through nodes in 2082, output

!----------------Reading flow discharge,catchment contributing areas, channel width and initial bed elevation from files-----------------------

 do m=1,s   !Depending the type of run, either pre-dam hydrograh or post-dam hydrograph should be chosen
     
!   read (53,*) qw(1,m), qfreq(m)         !pre-dam rating curve
    read(101,*) qw(1,m), qfreq(m)         !post-dam rating curve
    read(123,*) qw_meouge(m)   
 end do
 
 do i=1,n
   read (54,*) catcharea(i)   !  Catchment areas contributing to each computational node
   do m=1,s 
    if (i>56) then            !  Computing the discharge for each computational node based on contributing catchments. Meouge branch is added at node 57
      qw(i,m)= qw(1,m)*0.8*catcharea(i)+qw_meouge(m)+qw(1,m)
    else
      qw(i,m)= qw(1,m)*0.8*catcharea(i)+qw(1,m)
    end if

    b(i,m)= 0.8*(9.7238*qw(i,m)**0.5929)   ! This relation is obtained from hydraulic modeling of the study area using different discharges
   
    qww(i,m)=qw(i,m)/b(i,m)                ! Discharge per unit channel width for each node
   end do
 end do

 do i=1,n
   read(41,*) etaa(i)              ! Initial bed elevation, for pre-dam run it can be chosen arbitrarily because equilirbium elevation does not depend on initial bed elevation. For post-dam runs, it will be the equilirbium bed elevation obtained from the pre-dam run
   x(i)=i*dx
 end do

!----------------Sediment supply corresponding to each grain [m3/s]-----------------------


 Do i=1,13   
   if (i==1) then
     qbfi(i)=qbf*ff(i+1)                   
     qbfi_meouge(i)=qbf_meouge*ff(i+1)
   else
     qbfi(i)=qbf*(ff(i+1)-ff(i))
     qbfi_meouge(i)=qbf_meouge*(ff(i+1)-ff(i))
   end if
 end do   
   
!$$$$$$  qbfi(1)=qbf*ff(2)
!$$$$$$  qbfi(2)=qbf*(ff(3)-ff(2))
!$$$$$$  qbfi(3)=qbf*(ff(4)-ff(3))
!$$$$$$  qbfi(4)=qbf*(ff(5)-ff(4))
!$$$$$$  qbfi(5)=qbf*(ff(6)-ff(5))
!$$$$$$  qbfi(6)=qbf*(ff(7)-ff(6))
!$$$$$$  qbfi(7)=qbf*(ff(8)-ff(7))
!$$$$$$  qbfi(8)=qbf*(ff(9)-ff(8))
!$$$$$$  qbfi(9)=qbf*(ff(10)-ff(9))
!$$$$$$  qbfi(10)=qbf*(ff(11)-ff(10))
!$$$$$$  qbfi(11)=qbf*(ff(12)-ff(11))
!$$$$$$  qbfi(12)=qbf*(ff(13)-ff(12))
!$$$$$$  qbfi(13)=qbf*(ff(14)-ff(13))
!$$$$$$ !!!!!! Meouge!!!!!!!!!!!!!
!$$$$$$  qbfi_meouge(1)=qbf_meouge*ff(2)
!$$$$$$  qbfi_meouge(2)=qbf_meouge*(ff(3)-ff(2))
!$$$$$$  qbfi_meouge(3)=qbf_meouge*(ff(4)-ff(3))
!$$$$$$  qbfi_meouge(4)=qbf_meouge*(ff(5)-ff(4))
!$$$$$$  qbfi_meouge(5)=qbf_meouge*(ff(6)-ff(5))
!$$$$$$  qbfi_meouge(6)=qbf_meouge*(ff(7)-ff(6))
!$$$$$$  qbfi_meouge(7)=qbf_meouge*(ff(8)-ff(7))
!$$$$$$  qbfi_meouge(8)=qbf_meouge*(ff(9)-ff(8))
!$$$$$$  qbfi_meouge(9)=qbf_meouge*(ff(10)-ff(9))
!$$$$$$  qbfi_meouge(10)=qbf_meouge*(ff(11)-ff(10))
!$$$$$$  qbfi_meouge(11)=qbf_meouge*(ff(12)-ff(11))
!$$$$$$  qbfi_meouge(12)=qbf_meouge*(ff(13)-ff(12))
!$$$$$$  qbfi_meouge(13)=qbf_meouge*(ff(14)-ff(13))


!----------------Total sediment supply from upstream and Meouge [m3/s]-----------------------
do v=1,13
qbft=qbft+qbfi(v)
qbft_meouge=qbft_meouge+qbfi_meouge(v)
end do

!============================================================================================
!=========================Compute Initial Conditions==========================================
!============================================================================================


!----------------Quasi normal flow calculation-----------------------
ks=2*d90i           ! Initial roughness height of the bed
do i=1,n
  if (i<n) then
    slope(i)=(etaa(i)-etaa(i+1))/dx
  else
    slope(i)=(etaa(i-1)-etaa(i))/dx
  end if  
end do
do i=1,n  
   h(i)=(ks(i)**(1.0/3.0)*qww(i,1)**2.0/(alphr**2.0*g*slope(i)))**(3.0/10.0)  ! Manning-Strickler Formulation
   u(i)=qww(i,1)/h(i)
   cf(i)=(alphr**(-2.0))*(h(i)/ks(i))**(-1.0/3.0)
end do
   

!----------------Initial sediment fraction in active layer and subtrate-----------------------

do i=1,n
  do v=1,13
    fs(i,v,:)=f(v)              ! substrate grainsize fractions. It is assumed to be the same as the sediment supply
    read(67,*) fa(i,v)          ! active layer grain size fractions. For pre-dam run, it can be chosen arbitrarily. For post-dam runs, it will be the output of the pre-dam run.
  end do
  
  cumulative = 0.0
  
  do v=1,13                               !Normalizing the fractions to avoid unexpected divergence in the results
    cumulative = cumulative + fa(i,v)
  end do
  do v= 1,13
    fa(i,v) = fa(i,v)/cumulative
  end do
end do



!-----------------------Bedload transport (Wilcock and Crowe [2003])-----------------------

 qbt=0 
 do i=1,n
   fss(i)=0.0
   do v=1,13
    if (d(v)<0.0025) then
      fss(i)=f(v)+fss(i)
    end if
   end do
   taustarssrg(i)=(0.021+.015*exp(-20*fss(i)))
   ustar(i)=u(i)/cf(i)**(-.5)
   taustarsg(i)=ustar(i)**2/(r*g*dg)
  
   do v=1,13
    bb(i,v)=.67/(1.0+exp(1.5-d(v)/dg))
    phi(i,v)=(taustarsg(i)/taustarssrg(i))*(d(v)/dg)**(-bb(i,v))
    if (phi(i,v)<1.35) then
      gphi(i,v)=0.002*phi(i,v)**7.5
    else
      gphi(i,v)=14.0*(1.0-.894/phi(i,v)**.5)**4.5
    end if
    wstar(i,v)=gphi(i,v)
    if (taustarsg(i)< taustarssrg(i))  then
      qb(i,v)=0
    else
      qb(i,v)=(ustar(i)**3)*fa(i,v)*wstar(i,v)/(r*g)
    end if
    qbt(i)=qbt(i)+qb(i,v)
   end do
   do v=1,13
    if (qbt(i)/=0) then
      pb(i,v)= qb(i,v)/qbt(i)
    else
      pb(i,v)=0
    end if
   end do
 end do

!-----------------------Initialize the stratigraphy layer-----------------------

 do i=1,n
   p(i)= int((etaa(i)-la(i))/ls)+1
   l_buffer(i)=etaa(i)-la(i)-(p(i)-1)*ls
   do v=1,13
     f_buffer(i,v)=f(v)
     do k=1,p(i) 
      fs(i,v,k)=f(v)
     end do
   end do
 end do  

  

!-----------------------Print Initial Condition-----------------------

etaaout(:,1)=x       !distance from upstream is the first column of the output
etaaout(:,2)=etaa    !initial bed elevation is the second column of the output

slopeout(:,1)=x      !distance from upstream is the first column of the output
slopeout(:,2)=slope  !initial slope is the second column of the output

dgout(:,1)=x         !distance from upstream is the first column of the output
dgout(:,2)=dg        !initial geometric mean diameter of the surface material is the second column of the output

hout(:,1)=x          !distance from upstream is the first column of the output
hout(:,2)=h          !initial water depth is the second column of the output

d90out(:,1)=x        !distance from upstream is the first column of the output
d90out(:,2)=d90i     !initial d90 of the surface material is the second column of the output

d50out(:,1)=x        !distance from upstream is the first column of the output
d50out(:,2)=d50i     !initial d50 of the surface material is the second column of the output

qvolout(:,1)=x     !distance from upstream is the first column of the output

qbi2019(:,1)=x     !distance from upstream is the first column of the output
qbi2028(:,1)=x     !distance from upstream is the first column of the output
qbi2031(:,1)=x     !distance from upstream is the first column of the output
qbi2082(:,1)=x     !distance from upstream is the first column of the output
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!=========================Start the main calculations=========================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
                
do jjj=1,mprint                         ! time step for printing out results
  
  if (jjj==8) then                      ! one time reinjection in 2016, each jjj represents 3years
    etaa(1)=etaa(1)+2.45                ! for sediment reinjection (only post-dam model), we raise the bed elevation in the first two nodes of the model in such a way that the volume added due to the increas in the bed is equal to the reinjection volume
    etaa(2)=etaa(2)+1.25
  end if
  if (jjj>8) then                       ! periodic reinjection
    if (mod(jjj,2)==0) then             ! jjj is 3 years. The periodic reinjection can be adjusted in 3 years increments. It is now every 6 years.
      etaa(1)=etaa(1)+2.45
      etaa(2)=etaa(2)+1.25
    end if
  end if
  qvol=0
  qvoli=0
  do m=1,s                              !time steps for each flow discharge bin (discharge chosen from flow duration curve)
    
    do jj=1,mtoprint                   !time steps for the calculation between each printouts
   
!-----------------------Compute the new bed profile-----------------------
     do i=1,n
         if (i==1) then
           etaa(i)=etaa(i)-(qfreq(m)*(qbt(i)-qbf/b(1,m))/((1.0-lp)*dx))*dt
           detaa(i)=-qfreq(m)*(qbt(i)-qbf/b(1,m))/((1.0-lp)*dx)                            !detaa/dt is the time dependend bed elevation gradient. This term will be used in grain size specific conservation of mass equations 
         else
           etaa(i)=etaa(i)-(qfreq(m)*(qbt(i)-qbt(i-1))/((1.0-lp)*dx))*dt                   
           detaa(i)=-qfreq(m)*(qbt(i)-qbt(i-1))/((1.0-lp)*dx)                              !detaa/dt is the time dependend bed elevation gradient. This term will be used in grain size specific conservation of mass equations 
         end if   
     end do


!-----------------------Compute fl (the fractin of sediment at the interface of active layer and subtrate-----------------------   

   do i=1,n
     do v=1,13
       if (detaa(i)>0) then
         fl(i,v)=alph*fa(i,v)+(1-alph)*pb(i,v)
       else  
         fl(i,v)=f_buffer(i,v)
       end if
     end do
   end do

!-----------------------Compute fa i.e. sediment fractions in the active layer (grain size specific exner equations)-----------------------   
   do i=1,n
    do v=1,13
        if (i==1) then   
           fa(i,v)= (qfreq(m)*(-qb(i,v)+qbfi(v)/b(1,m))/(dx*(1-lp))-fl(i,v)*detaa(i))*(dt/la(i))+fa(i,v)
        else 
           fa(i,v)= (qfreq(m)*(-qb(i,v)+qb(i-1,v))/(dx*(1-lp))-fl(i,v)*detaa(i))*(dt/la(i))+fa(i,v)
        end if
    end do
!-----------------------Normalize the sediment fraction to avoid unexpected divergence-----------------------   
     cumulative = 0
     do v=1,13
	   cumulative = cumulative + fa(i,v)
     end do
     do v= 1,13
	   fa(i,v) = fa(i,v)/cumulative
     end do    
   end do
!----------------Geometric mean diameter of the surface sediment [m]------------------------                                
   do i=1,n
     psibarf=0               
     do v=1,13
       psibarf=psi(v)*fa(i,v)+psibarf
     end do
     dgf(i)=2**psibarf
     dgf(i)=dgf(i)/1000     
   end do
!--------------------------D90 of the surface sediment [m]------------------------  
   do i=1,n
     fau=0
     fad=0
     do v=1,14
       fau=fa(i,v)+fau
       if (fau>=.90)then
         fad=fau-fa(i,v)
         d90(i)=2.72**((log(db(v+1))-log(db(v)))*(0.9-fad)/(fau-fad)+log(db(v)))/1000.0
         exit 
       end if
     end do
   end do

!--------------------------D50 of the surface sediment [m]------------------------  
   do i=1,n
     fau=0
     fad=0
     do v=1,14
       fau=fa(i,v)+fau
       if (fau>=.50)then
         fad=fau-fa(i,v)
         d50(i)=2.72**((log(db(v+1))-log(db(v)))*(0.5-fad)/(fau-fad)+log(db(v)))/1000.0
         exit 
       end if
     end do
   end do
!--------------------------Active layer thickness and roughness height [m]------------------------  
   do i=1,n
     la(i)=2*d90(i)
     Ks(i)=2*d90(i)
   end do

!----------------------------------Quasi normal flow calculation------------------------------------

   do i=1,n
    if (i<n) then
     slope(i)=(etaa(i)-etaa(i+1))/dx               !bed slope
    else
     slope(i)=(etaa(i-1)-etaa(i))/dx
    end if  
   end do
  
   do i=1,n
    h(i)=(ks(i)**(1.0/3.0)*qww(i,m)**2.0/(alphr**2.0*g*slope(i)))**(3.0/10.0)   !Manning-Strickler formularion
    u(i)=qww(i,m)/h(i)
    cf(i)=(alphr**(-2.0))*(h(i)/ks(i))**(-1.0/3.0)
   end do
   
!------------------Store and access grain size distribution in the substrate (stratigraphy)-----------------------

  do i=1,n   
    if (detaa(i)< 0.0) then      !degradation
      if ((etaa(i)-la(i)-(p(i)-1)*ls)<0.0 ) then
        p(i)=p(i)-1
        l_buffer(i)= etaa(i)-la(i)-ls*(p(i)-1)
        do v=1,13
          f_buffer(i,v)=fs(i,v,p(i))
        end do
      else
        l_buffer(i)=etaa(i)-la(i)-ls*(p(i)-1)
        do v=1,13
          f_buffer(i,v)=fs(i,v,p(i))
        end do
      end if
    else                          !aggradation  
      if((etaa(i)-la(i)-(p(i)-1)*ls)>ls ) then          
        p(i)=p(i)+1
        l_buffer(i)=etaa(i)-la(i)-ls*(p(i)-1)
        do v=1,13
          fs(i,v,p(i))=fl(i,v)
          f_buffer(i,v)=fs(i,v,p(i))
        end do
      else
        do v=1,13
          f_buffer(i,v)=f_buffer(i,v)*l_buffer(i)/(l_buffer(i)+detaa(i))+fl(i,v)*detaa(i)/(l_buffer(i)+detaa(i))
       
        end do
          l_buffer(i)= etaa(i)-la(i)-ls*(p(i)-1)
      
       end if
    end if
  end do

!--------------------------geometric mean diameter of the substrate layers------------------------ 
  do i=1,n
      psibarf=0
     do k=1,p(i)
       psibarf=0
       do v=1,13
         psibarf=psi(v)*fs(i,v,k)+psibarf
       end do
       dg_sub(i,k)=2**psibarf
       dg_sub(i,k)=dg_sub(i,k)/1000  
     end do   
  end do      
!--------------------------Bedload transport (Wilcock and Crowe [2003])------------------------ 

   qbt=0 
   do i=1,n
    fss(i)=0.0
    do v=1,13
     if (d(v)<0.0025) then
      fss(i)=f(v)+fss(i)
     end if
    end do
    taustarssrg(i)=(0.021+.015*exp(-20*fss(i)))
    ustar(i)=u(i)/cf(i)**(-.5)
    taustarsg(i)=ustar(i)**2/(r*g*dg)
    
    do v=1,13
      bb(i,v)=.67/(1.0+exp(1.5-d(v)/dg))
      phi(i,v)=(taustarsg(i)/taustarssrg(i))*(d(v)/dg)**(-bb(i,v))
      if (phi(i,v)<1.35) then
        gphi(i,v)=0.002*phi(i,v)**7.5
      else
        gphi(i,v)=14.0*(1.0-.894/phi(i,v)**.5)**4.5
      end if
      wstar(i,v)=gphi(i,v)
      if (taustarsg(i)< taustarssrg(i))  then
        qb(i,v)=0
      else
        qb(i,v)=(ustar(i)**3)*fa(i,v)*wstar(i,v)/(r*g)             
      end if
!--------------------------Meouge branch sediment load------------------------ 
      if (i==57) then
        qb(i,v)=qb(i,v)+qbfi_meouge(v)/b(i,m)
      end if
!------------------------------------------------------------------------------ 
      qbt(i)=qbt(i)+qb(i,v)         !total sediment transport rate for each node
    end do     
      wstarg(i)=qbt(i)/taustarsg(i)
    do v=1,13
      if (qbt(i)/=0) then
        pb(i,v)= qb(i,v)/qbt(i)
      else
        pb(i,v)=0
      end if
    end do   
    qvol(i)=qvol(i)+qbt(i)*dt*b(i,m)*qfreq(m)   !total sediment volume transporting through each node
     do v=1,13
      qvoli(i,v)=qvoli(i,v)+qb(i,v)*dt*b(i,m)*qfreq(m)   !grain size specific sediment volume transporting through each node
     end do
   end do
 end do                                          ! end of the calculatin loop between each print outs (jj)
end do                                           ! end of the calculation loop for each discharge bin from flow duration curve

!---------------------------------Printouts-------------------------------- 

  etaaout(:,jjj+2)=etaa
  slopeout(:,jjj+2)=slope
  dgout(:,jjj+2)=dgf
  hout(:,jjj+2)=h
  qvolout(:,jjj+1)= qvol
  d90out(:,jjj+2)=d90
  d50out(:,jjj+2)= d50
  
!---------------------------------Print volume of each grain size passing though each node in years 2019,2028,2031 and 2082--------------------------------   

  if (jjj==9) then                 !2019
    do v=1,13
     qbi2019(:,v+1)=qvoli(:,v)/3.0
    end do 
  else if (jjj==12) then           !2028
    do v=1,13
     qbi2028(:,v+1)=qvoli(:,v)/3.0
    end do
  else if (jjj==13) then           !2031
    do v=1,13
     qbi2031(:,v+1)=qvoli(:,v)/3.0
    end do
  else if (jjj==30) then           !2082
    do v=1,13
     qbi2082(:,v+1)=qvoli(:,v)/3.0
    end do
  end if
!---------------------------------print grain size fraction of the active layer--------------------------------   
  
  do i=1,n                            !it is used to have the initial grain size distribution of the surface for the next run (use the pre-dam fa output as in input of the post-dam runs). Fa-output data should be copied/pasted in Fa-input
    do v=1,13
      write(68,*) fa(i,v)
    end do
  end do
 
end do                     ! end of the calculatin loop for printing out results (jjj)
 

 do i=1,n

!===============================================================================================================================================
!---------------------------each column of the following outputs represent one time step and each column is 3 years apart-----------------------
!===============================================================================================================================================

!---------------------------------bed elevation--------------------------------        
    write(1,*) etaaout(i,1),etaaout(i,2),etaaout(i,3),etaaout(i,4)
    write(13,*)etaaout(i,5),etaaout(i,6),etaaout(i,7),etaaout(i,8)
    write(14,*)etaaout(i,9),etaaout(i,10),etaaout(i,11),etaaout(i,12)
    write(30,*)etaaout(i,13),etaaout(i,14),etaaout(i,15),etaaout(i,16)
    write(31,*)etaaout(i,17),etaaout(i,18),etaaout(i,19),etaaout(i,20)
    write(37,*)etaaout(i,21),etaaout(i,22),etaaout(i,23),etaaout(i,24)
    write(38,*)etaaout(i,25),etaaout(i,26),etaaout(i,27),etaaout(i,28)
    write(39,*)etaaout(i,29),etaaout(i,30),etaaout(i,31),etaaout(i,32)

!---------------------------------volume of sediment transported through each node over the time step (3years)-------------------------------- 
   write(69,*) qvolout(i,1),qvolout(i,2),qvolout(i,3),qvolout(i,4)
   write(70,*) qvolout(i,5),qvolout(i,6),qvolout(i,7),qvolout(i,8)
   write(71,*) qvolout(i,9),qvolout(i,10),qvolout(i,11),qvolout(i,12)
   write(72,*) qvolout(i,13),qvolout(i,14),qvolout(i,15),qvolout(i,16)
   write(73,*) qvolout(i,17),qvolout(i,18),qvolout(i,19),qvolout(i,20)
   write(74,*) qvolout(i,21),qvolout(i,22),qvolout(i,23),qvolout(i,24)
   write(75,*) qvolout(i,25),qvolout(i,26),qvolout(i,27),qvolout(i,28)
   write(76,*) qvolout(i,29),qvolout(i,30)
   
!---------------------------------geometric mean diameter of the surface sediment-------------------------------- 
   write(4,*)dgout(i,1),dgout(i,2),dgout(i,3),dgout(i,4)
   write(17,*)dgout(i,5),dgout(i,6),dgout(i,7),dgout(i,8)
   write(18,*)dgout(i,9),dgout(i,10),dgout(i,11),dgout(i,12)
   write(43,*)dgout(i,13),dgout(i,14),dgout(i,15),dgout(i,16)
   write(44,*)dgout(i,17),dgout(i,18),dgout(i,19),dgout(i,20)
   write(45,*)dgout(i,21),dgout(i,22),dgout(i,23),dgout(i,24)
   write(46,*)dgout(i,25),dgout(i,26),dgout(i,27),dgout(i,28)
   write(47,*)dgout(i,29),dgout(i,30),dgout(i,31),dgout(i,32)
    
!---------------------------------d90 of the surface sediment-------------------------------- 
   write(77,*) d90out(i,1),d90out(i,2),d90out(i,3),d90out(i,4) 
   write(78,*) d90out(i,5),d90out(i,6),d90out(i,7),d90out(i,8)
   write(79,*) d90out(i,9),d90out(i,10),d90out(i,11),d90out(i,12)
   write(80,*) d90out(i,13),d90out(i,14),d90out(i,15),d90out(i,16)
   write(81,*) d90out(i,17),d90out(i,18),d90out(i,19),d90out(i,20)
   write(82,*) d90out(i,21),d90out(i,22),d90out(i,23),d90out(i,24)
   write(83,*) d90out(i,25),d90out(i,26),d90out(i,27),d90out(i,28)
   write(84,*) d90out(i,29),d90out(i,30),d90out(i,31)
!---------------------------------d50 of the surface sediment-------------------------------- 
   write(85,*) d50out(i,1),d50out(i,2),d50out(i,3),d50out(i,4) 
   write(86,*) d50out(i,5),d50out(i,6),d50out(i,7),d50out(i,8)
   write(87,*) d50out(i,9),d50out(i,10),d50out(i,11),d50out(i,12)
   write(88,*) d50out(i,13),d50out(i,14),d50out(i,15),d50out(i,16)
   write(89,*) d50out(i,17),d50out(i,18),d50out(i,19),d50out(i,20)
   write(90,*) d50out(i,21),d50out(i,22),d50out(i,23),d50out(i,24)
   write(91,*) d50out(i,25),d50out(i,26),d50out(i,27),d50out(i,28)
   write(92,*) d50out(i,29),d50out(i,30),d50out(i,31)
   
!---------------------------------water depth--------------------------------
   write(5,*) hout(i,1),hout(i,2),hout(i,3),hout(i,4)
   write(19,*)hout(i,5),hout(i,6),hout(i,7),hout(i,8)
   write(20,*)hout(i,9),hout(i,10),hout(i,11),hout(i,12)
   write(32,*)hout(i,13),hout(i,14),hout(i,15),hout(i,16)
   write(33,*)hout(i,17),hout(i,18),hout(i,19),hout(i,20)
   write(34,*)hout(i,21),hout(i,22),hout(i,23),hout(i,24)
   write(35,*)hout(i,25),hout(i,26),hout(i,27),hout(i,28)
   write(36,*)hout(i,29),hout(i,30),hout(i,31),hout(i,32)
    
!---------------------------------Bed slope--------------------------------
    write(7,*) slopeout(i,1),slopeout(i,2),slopeout(i,3),slopeout(i,4)
    write(21,*)slopeout(i,5),slopeout(i,6),slopeout(i,7),slopeout(i,8)
    write(22,*)slopeout(i,9),slopeout(i,10),slopeout(i,11),slopeout(i,12)
    write(102,*) slopeout(i,13),slopeout(i,14),slopeout(i,15),slopeout(i,16)
    write(103,*) slopeout(i,17),slopeout(i,18),slopeout(i,19),slopeout(i,20)
    write(104,*) slopeout(i,21),slopeout(i,22),slopeout(i,23),slopeout(i,24)
    write(105,*) slopeout(i,25),slopeout(i,26),slopeout(i,27),slopeout(i,28)
    write(106,*) slopeout(i,29),slopeout(i,30),slopeout(i,31),slopeout(i,32)
    
 

   
!---------------------------------volume transported thorugh each node in the year 2019--------------------------------
    
    write(107,*) qbi2019(i,1),qbi2019(i,2),qbi2019(i,3),qbi2019(i,4)
    write(108,*) qbi2019(i,5),qbi2019(i,6),qbi2019(i,7),qbi2019(i,8)
    write(109,*) qbi2019(i,9),qbi2019(i,10),qbi2019(i,11),qbi2019(i,12)
    write(110,*) qbi2019(i,13)
!---------------------------------volume transported thorugh each node in the year 2028--------------------------------
    write(111,*) qbi2028(i,9),qbi2028(i,10),qbi2028(i,11),qbi2028(i,12)
    write(112,*) qbi2028(i,5),qbi2028(i,6),qbi2028(i,7),qbi2028(i,8)
    write(113,*) qbi2028(i,9),qbi2028(i,10),qbi2028(i,11),qbi2028(i,12)
    write(114,*) qbi2028(i,13)
!---------------------------------volume transported thorugh each node in the year 2031--------------------------------
    write(115,*) qbi2031(i,9),qbi2031(i,10),qbi2031(i,11),qbi2031(i,12)
    write(116,*) qbi2031(i,5),qbi2031(i,6),qbi2031(i,7),qbi2031(i,8)
    write(117,*) qbi2031(i,9),qbi2031(i,10),qbi2031(i,11),qbi2031(i,12)
    write(118,*) qbi2031(i,13)
!---------------------------------volume transported thorugh each node in the year 2082--------------------------------
    write(119,*) qbi2082(i,9),qbi2082(i,10),qbi2082(i,11),qbi2082(i,12)
    write(120,*) qbi2082(i,5),qbi2082(i,6),qbi2082(i,7),qbi2082(i,8)
    write(121,*) qbi2082(i,9),qbi2082(i,10),qbi2082(i,11),qbi2082(i,12)
    write(122,*) qbi2082(i,13)
 end do

 
end program Morphodynamic_gravelbed
