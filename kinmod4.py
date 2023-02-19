''' Our Module '''

from sympy import solve
from sympy.physics.mechanics import cross, dot
from matplotlib import pyplot as plt
from math import sqrt,cos,sin,pi,radians

def draw_line_with_arrow(x1,x2,y1,y2,lname):
    xlen=sqrt((x1-x2)**2.0+(y1-y2)**2.0)
    c1=(x2-x1)/xlen
    s1=(y2-y1)/xlen
    alp=5*pi/6
    c2=c1*cos(alp)-s1*sin(alp)
    s2=s1*cos(alp)+c1*sin(alp)
    c3=c1*cos(alp)+s1*sin(alp)
    s3=s1*cos(alp)-c1*sin(alp)
    arrowlen=xlen/10
    x3=x2+arrowlen*c2
    y3=y2+arrowlen*s2
    x4=x2+arrowlen*c3
    y4=y2+arrowlen*s3

    plt.plot([x1,x2],[y1,y2],label=lname)
    plt.legend(loc='best')
    plt.plot([x2,x3],[y2,y3],'k')
    plt.plot([x2,x4],[y2,y4],'k')
    plt.plot([x3,x4],[y3,y4],'k')



      
    
class ourpoint:
    def __init__(self, name, *args):
        ''' Defines a point by its name and coordinates in a single or successive reference frames
            Example: B=ourpoint('B',N,[0,0])
                     C=ourpoint('C',N,[10,0],N1,[12,20]) # C.pos=10*N.x+ 12*N1.x +20*N1.y
        '''
        self.name=name
        self.pos=0*args[0].x #Initialization
        
        for i in range(0,len(args),2):
            self.pos=self.pos+args[i+1][0]*args[i].x+args[i+1][1]*args[i].y
            
        self.vt=[0*args[0].x + 0*args[0].y, \
                 0*args[0].x + 0*args[0].y, \
                 0*args[0].x + 0*args[0].y]
        self.at=[0*args[0].x + 0*args[0].y, \
                 0*args[0].x + 0*args[0].y, \
                 0*args[0].x + 0*args[0].y, \
                 0*args[0].x + 0*args[0].y, \
                 0*args[0].x + 0*args[0].y]
        self.vel=0*args[0].x
        self.accl=0*args[0].x
        self.namev=[]
        self.namea=[]
        print(self.name, 'has a position vector',self.pos)
        print()
        self.force = 0 * args[0].x
        self.moment = 0 * args[0].x
        
            
    def setvel_direct(self,ref,v):
        '''
            Example: B.setvel_direct(N,[10,20])
            Sets a velocity 10*N.x + 20*N.y to ourpoint object B
        '''
        self.vt[0]=v[0]*ref.x + v[1]*ref.y
        self.vt[1]=0*ref.x+0*ref.y
        self.vt[2]=0*ref.x+0*ref.y
        self.vel=self.vt[0]
        self.namev=['v'+ self.name,'','']
        print(self.name, 'has a velocity',self.vel)
        print()
        
        
    def calvel_rigid(self,omega,P):
        '''
            Example: B.calvel_rigid(10*N.z,A)
            B and A are ourpoint objects. Both points are rigidly attached to a body
            Velocity of A known.
            Angular velocity of the body 10*N.z
            Computes velocity of B
        '''
        self.vt[0]=P.vel
        self.vt[1]=cross(omega,self.pos-P.pos)
        self.vel=self.vt[0]+self.vt[1]
        self.namev=['v'+P.name,'v'+self.name+'/'+P.name,'']
        print(self.name,' rigidly connected to ', P.name, ' in body ',\
              P.name,'-',self.name,sep='')
        print(self.name, 'has a velocity ',self.vel)
        print()
         
    def calvel_sliding(self,omega,P,ref,vrel):
        '''
            Example: B.calvel_sliding(10*N.z,A,N1,[10,20])
            B and A are ourpoint objects. 
            Angular velocity of the body 10*N.z
            A rigidly attached to body.
            Velocity of A known.
            B has a relative velocity 10*N1.x + 20*N1.y w.r.t the body
            Computes velocity of B
        '''        
        self.vt[0]=P.vel
        self.vt[1]=cross(omega,self.pos-P.pos)        
        self.vt[2]=vrel[0]*ref.x + vrel[1]*ref.y
        self.vel=self.vt[0]+self.vt[1]+self.vt[2]
        self.namev=['v'+P.name,'v'+self.name+'c/'+P.name,'v_rel']
        print(self.name,' sliding wrt coincident point',self.name+'c')
        print(self.name+'c',' rigidly connected to ', P.name, ' in body ',\
              P.name,'-',self.name+'c',sep='')
        print(self.name, 'has a velocity ',self.vel)
        print()

    def setaccl_direct(self,ref,a):
        '''
            Example: B.setaccl_direct(N,[10,20])
            Sets an acceleration 10*N.x + 20*N.y to ourpoint object B
        '''               
        self.at[0]=a[0]*ref.x + a[1]*ref.y
        self.at[1]=0*ref.x+0*ref.y
        self.at[2]=0*ref.x+0*ref.y
        self.at[3]=0*ref.x+0*ref.y
        self.at[4]=0*ref.x+0*ref.y       
        self.accl=self.at[0]
        self.namea=['a'+self.name,'','','','']
        print(self.name, 'has an acceleration',self.accl)
        print()
        
    def calaccl_rigid(self,omega,alpha,P):
        
        '''
            Example: B.calaccl_rigid(10*N.z,2*N.z,A)
            B and A are ourpoint objects. Both points are rigidly attached to the body.
            Acceleration of A known
            Angular velocity of the body 10*N.z,Angular acceleration of the body 2*N.z 
            Computes acceleration of B
        '''        
        self.at[0]=P.accl
        self.at[1]=cross(alpha,self.pos-P.pos)
        self.at[2]=cross(omega,cross(omega,self.pos-P.pos))
        self.accl=self.at[0]+self.at[1]+self.at[2]
        self.namea=['a'+P.name,\
                     '(a'+self.name+'/'+P.name+')'+'t',\
                     '(a'+self.name+'/'+P.name+')'+'n','','']
        print(self.name,' rigidly connected to ', P.name, ' in body ',\
              P.name,'-',self.name,sep='')
        print(self.name, 'has a acceleration ',self.accl)
        print()
         
         
    def calaccl_sliding(self,omega,alpha,P,refv,vrel,refa,arel):
        '''
            Example: B.calaccl_sliding(10*N.z,2*N.z,A,N1,[10,20],N1,[1,2])
            B and A are ourpoint objects.
            Acceleration of A known
            Angular velocity of the body 10*N.z, angular acceleration of the body 2*N.z
            A rigidly attached to body.
            B has a relative velocity 10*N1.x + 20*N1.y and relative accl 1*N1.x + 2*N1.y w.r.t the body
            Computes acceleration of B
        '''        
        self.at[0]=P.accl
        self.at[1]=cross(alpha,self.pos-P.pos)
        self.at[2]=cross(omega,cross(omega,self.pos-P.pos))
        self.at[4]=arel[0]*refa.x+arel[1]*refa.y
        self.at[3]=2*cross(omega,vrel[0]*refv.x + vrel[1]*refv.y)     
        self.accl=self.at[0]+self.at[1]+self.at[2]+self.at[3]+self.at[4]      
        self.namea=['a'+P.name,\
                     '(a'+self.name+'c/'+P.name+')'+'t',\
                     '(a'+self.name+'c/'+P.name+')'+'n','a_Cor','a_rel']
        print(self.name,' sliding wrt coincident point',self.name+'c')
        print(self.name+'c',' rigidly connected to ', P.name, ' in body ',\
              P.name,'-',self.name+'c',sep='')
        print(self.name, 'has a acceleration ',self.accl)
        print()
         

    def setfor_xy(self, ref, f):
        '''
            Example: B.setfor_xy(N,[10,20])
            Sets a force 10*N.x + 20*N.y to ourpoint object B
        '''
        force = f[0] * ref.x + f[1] * ref.y

        self.force = self.force + force
        print(self.name, '-> applied current force', force)
        print(self.name, '-> total force', self.force)
        print()


    def setfor_mag_theta(self, ref, f):
        '''
            Example: B.setfor_mag_theta(N,[10,20])
            Sets a force 10*N.x + 20*N.y to ourpoint object B
        '''
        force = f[0] * cos(radians(f[1])) * ref.x + f[0] * sin(radians(f[1]))  * ref.y

        self.force = self.force + force
        print(self.name, '-> applied current force', force)
        print(self.name, '-> total force', self.force)
        print()


    def setmom_direct(self, ref, m):
        '''
            Example: B.setmom_direct(N,10)
            Sets a moment 10*N.z to ourpoint object B
        '''
        moment = m * ref.z

        self.moment = self.moment + moment
        print(self.name, '-> applied current moment', self.moment)
        print(self.name, '-> total moment', self.moment)
        print()


class oureqn:  
    def __init__(self):
        
        '''
            Defines an equation object which has methods for generating
            velocity and acceleration equations
            Example: e1=oureqn()
        '''
        
        print('Equation object created')
        
        
    def get_eqnv(self,p1,p2,ref):
        '''
           Example: e1.get_eqnv(B1,B2,N) -> Get 2 scalar velocity equations along N.x and N.y 
           by equating velocity of B computed using 2 different rigid bodies
        '''
        self.vel1=p1.vel
        self.vel2=p2.vel
        self.vt1=p1.vt.copy()
        self.vt2=p2.vt.copy()
        self.namev1=p1.namev.copy()
        self.namev2=p2.namev.copy()   

        v1=self.vel1.express(ref)
        v2=self.vel2.express(ref)
        self.eqnv1=dot((v1-v2),ref.x)
        self.eqnv2=dot((v1-v2),ref.y)
        print()
        print('Vel Eqn 1: ',self.eqnv1.n(3),'= 0')
        print('Vel Eqn 2: ',self.eqnv2.n(3),'= 0')
        print() 


    def subsv(self,**keyargs):
        '''
            Example: e1.subsv(a=0.70,b=-1.05) -> Substitute numerical values for symbols in vel equations
        '''
        for k in keyargs:
            self.eqnv1=self.eqnv1.subs(k,keyargs[k])
            self.eqnv2=self.eqnv2.subs(k,keyargs[k])
            print()
            print('Substitution of',k,'by',keyargs[k])
            print('Vel eqn 1 after subst:',self.eqnv1.n(3),'= 0')
            print('Vel eqn 2 after subst:',self.eqnv2.n(3),'= 0')
            print()
            

    def solvev(self,unk1,unk2):
        '''
           Solve velocity equations for 2 unknowns (symbolic variables)
           Example: e1.solvev('omBC','vrel')
        '''
        vsol=solve([self.eqnv1,self.eqnv2],[unk1,unk2])
        self.vsol=vsol
        print()
        print('Vel soln:',vsol)
        print()
        key_vsol=list(vsol.keys())
        val_vsol=list(vsol.values())
        for i in range(len(self.vt2)):
            #printrange((item)
            for j in range(len(key_vsol)):
                self.vt1[i]=self.vt1[i].subs(key_vsol[j],val_vsol[j])
                self.vt2[i]=self.vt2[i].subs(key_vsol[j],val_vsol[j])
                
                
    def get_eqna(self,p1,p2,ref):
        '''
           Example: e1.get_eqna(B1,B2,N) -> Get 2 scalar acceleration equations along N.x and N.y 
           by equating acceleration of B computed using 2 different rigid bodies
        '''
        self.accl1=p1.accl
        self.accl2=p2.accl
        self.at1=p1.at.copy()
        self.at2=p2.at.copy()
        self.namea1=p1.namea.copy()
        self.namea2=p2.namea.copy()   
               
        a1=self.accl1.express(ref)
        a2=self.accl2.express(ref)
        self.eqna1=dot((a1-a2),ref.x)
        self.eqna2=dot((a1-a2),ref.y)
        print()
        print('Accl Eqn 1: ',self.eqna1.n(3),'= 0')
        print('Accl Eqn 2: ',self.eqna2.n(3),'= 0')
        print()
        
        
    def subsa(self,**keyargs):
        '''
            Substitute numerical values for symbols in accl equations
            Examples :e1.subsa(om_AB=5,vrel=-1.05) 
                    e1.subsa(omBC=e2.vsol[omBC],vrel=e2.vsol[vrel])
         '''   
        for k in keyargs:
            self.eqna1=self.eqna1.subs(k,keyargs[k])
            self.eqna2=self.eqna2.subs(k,keyargs[k])
            print()
            print('Substitution of',k,'by',keyargs[k])
            print('Accl eqn 1 after subst:',self.eqna1.n(3),'= 0')
            print('Accl eqn 2 after subst:',self.eqna2.n(3),'= 0')
            print()
            
            
    def solvea(self,unk1,unk2):
        '''
           Solve acceleration equations for 2 unknowns (symbolic variables)
           Example: e1.solvea('alpha_BC','arel')
        '''
        asol=solve([self.eqna1,self.eqna2],[unk1,unk2])
        self.asol=asol
        print()
        print('Accl soln:',asol)
        print()
        key_asol=list(asol.keys())
        val_asol=list(asol.values())
        for i in range(len(self.at2)):
            #printrange((item)
            for j in range(len(key_asol)):
                self.at1[i]=self.at1[i].subs(key_asol[j],val_asol[j])
                self.at2[i]=self.at2[i].subs(key_asol[j],val_asol[j])
        key_asol=list(self.vsol.keys())
        val_asol=list(self.vsol.values())
        for i in range(len(self.at2)):
            #printrange((item)
            for j in range(len(key_asol)):
                self.at1[i]=self.at1[i].subs(key_asol[j],val_asol[j])
                self.at2[i]=self.at2[i].subs(key_asol[j],val_asol[j])
                
                
    def draw(self,ref,option):
        '''
            Example: e1.draw(N,'v') : -> Draw velocity diagram with components along N.x and N.y
                     e1.draw(N1,'a') : -> Draw acceleration diagram with components along N1.x and N1.y
        '''

        def callback(event):
            xp,yp=event.xdata,event.ydata
            try:
                for i in range(len(ldata1)):
                    xm,ym,xl,dx,dy=ldata1[i][0],ldata1[i][1],ldata1[i][2],\
                                    ldata1[i][3],ldata1[i][4]
                    if sqrt((xm-xp)**2+(ym-yp)**2)<0.05*xl:
                        print("%10s, Mag = %12.5e, X-comp = %12.5e, Y-comp = %12.5e"\
                              %(nlst1[i],xl,dx,dy))
                
                for i in range(len(ldata2)):
                    xm,ym,xl,dx,dy=ldata2[i][0],ldata2[i][1],ldata2[i][2],\
                                    ldata2[i][3],ldata2[i][4]
                    if sqrt((xm-xp)**2+(ym-yp)**2)<0.05*xl:
                        print("%10s, Mag = %12.5e, X-comp = %12.5e, Y-comp = %12.5e"\
                              %(nlst2[i],xl,dx,dy)) 
                               
            except:
                print('Probably mouse is outside plot window')
            
        #global ldata1, ldata1, nlst1, nlst2
        plt.close('all')
        fig=plt.figure(1)
        xst,yst=0,0
        ldata1,ldata2=[],[]
        if option=='v':
            lst1=self.vt1
            nlst1=self.namev1
            lst2=self.vt2
            nlst2=self.namev2            
        if option=='a':
            lst1=self.at1
            nlst1=self.namea1
            lst2=self.at2
            nlst2=self.namea2                 
        for i in range(len(lst1)) :
            item=lst1[i]
            delv=item.express(ref)
            delx=dot(delv,ref.x)
            dely=dot(delv,ref.y)
            #print(delx,dely)
            xfi=xst+delx
            yfi=yst+dely
            ldata1.append([xst+delx/2,yst+dely/2,sqrt(delx**2+dely**2),delx,dely])
            if delx !=0 or dely !=0:
                draw_line_with_arrow(xst,xfi,yst,yfi,nlst1[i])
            xst,yst=xfi,yfi
        xst,yst=0,0
        for i in range(len(lst2)) :
            item=lst2[i]
            delv=item.express(ref)
            delx=dot(delv,ref.x)
            dely=dot(delv,ref.y)
            #print(delx,dely)
            xfi=xst+delx
            yfi=yst+dely
            ldata2.append([xst+delx/2,yst+dely/2,sqrt(delx**2+dely**2),delx,dely])
            if delx !=0 or dely !=0:
                draw_line_with_arrow(xst,xfi,yst,yfi,nlst2[i])
            xst,yst=xfi,yfi
        
        fig.canvas.callbacks.connect('button_press_event', callback)
        plt.axis('equal')
        plt.show()
        #plt.hold()


class ourbody():

    def __init__(self,name,points,P_mc):
        '''
            Example : B1 = ourbody('Body1',[P1,P2,P3],P4)
            P1,P2,P3 are the points where forces are prescribed
            P_mc -> point about which moments are to be computed -> P4 here
        '''
        
        self.name = name
        self.points = points
        self.cm = 'Not defined'
        self.mc = P_mc
        self.dl=[]
        self.mass=0
        self.inertia=0
        self.accl=0
        self.ang_accl=0
    
    
    def set_for_dynamics(self,P_cm,mass,P_inertia,inertia,ref,accl,ang_accl):

        '''
            Example: B.set_for_dynamics(G,m,P1,(1/12.)*m*l**2,N,[0,alpha*(l/2)],alpha)
            
            P_cm -> center of mass -> G here
            mass -> mass of the body -> m
            P_inertia -> Point about which inertia value is supplied
            inertia -> Inertia value about point P_inertia, P1 here
            (program converts it to inertia about cm)
            ref -> reference in which linear acceleration is expressed -> N here
            accl -> list with two values, x and y components of acceleration of cm -> [0,alpha*(l/2)] here
            ang_accl -> value of angular acceleration -> ang_accl
        '''

        self.cm=P_cm
        self.mass=mass
        vec=(P_inertia.pos-self.cm.pos)
        xlen=vec.magnitude()
        self.inertia=inertia - mass*xlen*xlen
        self.accl=accl[0]*ref.x+accl[1]*ref.y
        self.ang_accl=ang_accl*ref.z
    
    
    
    def input_dist_load(self,dl_list):
        
        '''
        Example: B.input_dist_load([A1,A2,5,10])
        
        Right trapezoidal distributed load. A1 and A2 are 2 ourpoint objects.
        heights at A1 and A2 are 5 and 10 respectively. Load is perpendicular
        to line A1-A2. Positive direction is taken vectorially.
        
        Positive dir = Unit vector from A1 to A2 x k
        
        '''
        
        self.dl.append(dl_list)
        
        
    
    def generate_eqn(self,ref):
        
        '''
        Example: B.generate_eqn(N)
        
        Generates 2 force equations taking components along N.x and N.y
        and 1 moment equation for a body
        
        '''
        
        
        self.force = self.points[0].force
        for i in range (1,len(self.points)):
            self.force = self.force + self.points[i].force
        self.force = self.force - self.mass*self.accl
        
        self.mom = cross((self.points[0].pos-self.mc.pos),self.points[0].force)+self.points[0].moment
        for i in range (1,len(self.points)):
            self.mom = self.mom + cross((self.points[i].pos-self.mc.pos),self.points[i].force)
        if self.mass != 0:
            self.mom=self.mom - self.inertia*self.ang_accl  \
                      -cross((self.cm.pos-self.mc.pos),self.mass*self.accl)


        for i in range(len(self.dl)):
            h1=self.dl[i][2]
            h2=self.dl[i][3]
            fac=(1/3)*(h1+2*h2)/(h1+h2)
            vec=(self.dl[i][1].pos-self.dl[i][0].pos)
            dlcm=self.dl[i][0].pos + fac*vec
            dllen=vec.magnitude()
            vecn=cross(ref.z,vec)
            uvecn=vecn/vecn.magnitude()
            dlarea=(h1+h2)*dllen/2
            totforce=dlarea*uvecn
            print('Distributed Load',i+1,' - Total Force',totforce)
            print('Distributed Load',i+1,' - Centroid',dlcm)
            self.force=self.force+dlarea*uvecn
            self.mom=self.mom + cross((dlcm-self.mc.pos),totforce)


        self.eqn1 = dot(self.force,ref.x)
        self.eqn2 = dot(self.force,ref.y)
        self.eqn3 = dot(self.mom,ref.z)
        
        print('Equation 1: ',self.eqn1,'= 0')
        print('Equation 2: ',self.eqn2,'= 0')
        print('Equation 3: ',self.eqn3,'= 0')
        
        
    def eqn_substitute(self,**keyargs):
        
        '''      
        Substitute numerical values for symbols in equilibrium equations
        for a body
            Example : B.eqn_substitute(a=5,b=-1.05)           
        '''
        for k in keyargs:
            self.eqn1=self.eqn1.subs(k,keyargs[k])
            self.eqn2=self.eqn2.subs(k,keyargs[k])
            self.eqn3=self.eqn3.subs(k,keyargs[k])
            print()
            print('Substitution of',k,'by',keyargs[k])
            print('Eqn 1 after subst:',self.eqn1.n(3),'=0')
            print('Eqn 2 after subst:',self.eqn2.n(3),'=0')
            print('Eqn 3 after subst:',self.eqn3.n(3),'=0')
            print()
            
            
    def solve_eqn(self,unk):
        
        '''
        Bd.solve_eqn([P1x,P2y,P1y])
        Solve for the unknowns P1x, P2y, P1y
        
        '''    
        
        
        sol=solve([self.eqn1,self.eqn2,self.eqn3],unk)
        self.sol=sol
        for i in range(len(unk)):
            print(unk[i],' = ',self.sol[unk[i]])
            
        key_sol=list(self.sol.keys())
        val_sol=list(self.sol.values())
        
   
        for j in range(len(self.points)):
            for k in range(len(key_sol)):
                self.points[j].force = self.points[j].force.subs(key_sol[k],val_sol[k])
                self.points[j].moment = self.points[j].moment.subs(key_sol[k],val_sol[k])
    
            
    def draw_sfbm(self,ref):
        
        
        
        nlen=len(self.points)
        xc1 = self.points[0].pos
        xc2 = self.points[nlen-1].pos
        
        ndiv=100
        diff = (xc2-xc1)/ndiv
        
        sf = (ndiv+1) * [0*ref.x]
        bm = (ndiv+1) * [0*ref.x]
        xc = (ndiv+1) * [0*ref.x]
        
        for i in range(ndiv+1):
            xc[i] = xc1 + i*diff
            #print(xc[i])
            for j in range(len(self.points)):
                if xc[i].magnitude() > self.points[j].pos.magnitude():
                    #print(self.points[j].name)
                    bm[i] = bm[i] - cross((self.points[j].pos-xc[i]),self.points[j].force)
                    bm[i] = bm[i] - self.points[j].moment
                    sf[i] = sf[i] + self.points[j].force
            
            
            for j in range(len(self.dl)):
                if xc[i].magnitude() > self.dl[j][0].pos.magnitude():                       
                    h1=self.dl[j][2]
                    h2_act=self.dl[j][3]
                    if xc[i].magnitude() < self.dl[j][1].pos.magnitude():
                        fac1 = ((xc[i] - self.dl[j][0].pos).magnitude())/ \
                                 (self.dl[j][1].pos - self.dl[j][0].pos).magnitude()
                        vec=(xc[i]-self.dl[j][0].pos) 
                    else:
                        fac1=1
                        vec=(self.dl[j][1].pos-self.dl[j][0].pos) 
                    h2 = h1 + fac1*(h2_act-h1)
                    fac=(1/3)*(h1+2*h2)/(h1+h2)
                    dlcm=self.dl[j][0].pos + fac*vec
                    dllen=vec.magnitude()
                    vecn=cross(ref.z,vec)
                    uvecn=vecn/vecn.magnitude()
                    dlarea=(h1+h2)*dllen/2
                    totforce=dlarea*uvecn
                    #print('totforce',totforce)
                    sf[i]=sf[i]+totforce
                    bm[i]=bm[i] - cross((dlcm-xc[i]),totforce)
                    
            #print(sf[i].magnitude())
            #print(xc[i])
                
            sf[i]=float(dot(sf[i],ref.y))
            bm[i]=float(dot(bm[i],ref.z))
            xc[i]=float(xc[i].magnitude())
            #print(bm[i].magnitude())
            #print(xc[i].magnitude())
        fig,a=plt.subplots(2)
        #fig.suptitle('SF & BM diagrams')
        a[0].plot(xc,sf,'.',color='m')
        a[0].grid()
        a[1].plot(xc,bm,'.',color='c')
        a[1].grid()
        a[0].set_title('SF')
        a[1].set_title('BM')
        plt.show()
    
    
class group_of_bodies():    

    def __init__(self,name,bodies):
        
        
        '''
        Example : group1 = group_of_bodies('group1',[B1,B2])
        
        Generates a group group1 combining 2 bodies B1 & B2
        Equations are generated by appending equations for each
        constituent body
        
        '''
        self.name=name
        self.eqns=[]
        self.bodies=bodies.copy()
        
        for i in range(len(bodies)):
            self.eqns.append(bodies[i].eqn1)
            self.eqns.append(bodies[i].eqn2)
            self.eqns.append(bodies[i].eqn3)                
        #print(self.eqns)
        

    def gr_solve(self,unk):
        
        '''
        Example : group1.gr_solve([Bx,By,Cx,Cy,Dx,Dy])
        
        Solve equations of the group group1 for the unknowns
        Bx,By,Cx,Cy,Dx,Dy
        
        '''      
        sol=solve(self.eqns,unk)
        self.sol=sol
        for i in range(len(unk)):
            print(unk[i],' = ',self.sol[unk[i]])       

        key_sol=list(self.sol.keys())
        val_sol=list(self.sol.values())
        
   
        for i in range(len(self.bodies)):      
            for j in range(len(self.bodies[i].points)):
                for k in range(len(key_sol)):
                    self.bodies[i].points[j].force=\
                            self.bodies[i].points[j].force.subs(key_sol[k],val_sol[k])
                    self.bodies[i].points[j].moment=\
                            self.bodies[i].points[j].moment.subs(key_sol[k],val_sol[k])