'''
Created on Oct 28, 2009

@author: jakub
'''
'''
Created on Oct 28, 2009

@author: Andreas
'''

from scipy.spatial.distance import cdist
from numpy import delete,vstack, typename, array, append, reshape, sort, zeros, size, sum, transpose, swapaxes, mgrid, zeros_like, where, vdot, shape, min, average, repeat, insert, pi, arccos, dtype, sign
from numpy.linalg import solve
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, Dict, \
     Class  
from enthought.mayavi.mlab import surf, colorbar, show, mesh
from time import time

class LevelSet_discrete(HasTraits):
    
    out_grid = Array                #grid 
    work_grid = Array               #different type of order ([x0,y0],[x1,y1],[x2,y2]...[xn,yn])
    input_pts =Array
    scalars = Array                 #scalar distance to center line
    x_elements = Int(100.0)         #number of elements in x-direction used by set_grid, get_segment
    y_elements = Int(100.0)         #number of elements in y-direction used by set_grid, get_segment
    x_limit = array([-10.0,10.0])   #limits in x-direction used by set_grid, get_segment
    y_limit = array([-10.0,10.0])   #limits in y-direction used by set_grid, get_segment
    
    #values from preprocessor if used
    
    check_imput = False 
    line = Array
    length = Array
    normaverage = Array
    normnormal = Array
    
    def set_grid(self):
        """ defines a grid within a certain area with x_e elements in x and y_e elements in y direction. 
            values x_elements, y_elements, x_limits, y_limits should be defined in advance before starting set_grid"""        
        x_e=self.x_elements
        y_e=self.y_elements
        x_l=self.x_limit
        y_l=self.y_limit      
        a=(x_l[1]-x_l[0])/x_e
        b=(y_l[1]-y_l[0])/y_e
        self.out_grid=mgrid[x_l[0]:x_l[1]:a, y_l[0]:y_l[1]:b]
        x,y=self.out_grid
        self.work_grid=transpose(append(x.flatten(),y.flatten()).reshape(2,size(x)))       
        return self.out_grid
    
    def check_input_pts(self):
        """
        check_input_pts checks all input_pts for linear dependency of normal-average [i]&[i+1]vectors.
        if normal of line is identical to average normal of points, end point of line will be deleted
        if normal of line is unequal to average normal of points, new imput_pt will be added within line
        If line is closed, first and last point of input_pts have to be identical.
        """
        self.check_input=True
        input_pts=self.input_pts
        
        # calculation for line closed
        if input_pts[0].all()==input_pts[-1].all():
            line_closed=True
            input_pts=append(input_pts, input_pts[1]).reshape(-1,2)
            line=array(input_pts[1:]-input_pts[:-1])            
            normal=transpose(array((append(-line[:,1], line[:,0])).reshape(2,size(line)/2)))
            normnormal=normal*(sum(normal**2,axis=1)**-.5).reshape(len(normal),1)          
            average=(normnormal[:-1]+normnormal[1:])/2
            normaverage=average*(sum(average**2,axis=1)**-.5).reshape(len(average),1)
            normaverage=append(normaverage[-1],normaverage[:]).reshape(-1,2)
            # delete not needed input_pts which where only doubled to get normal average of connected point
            input_pts=delete(input_pts,-1,-2).reshape(-1,2)
            line=array(input_pts[1:]-input_pts[:-1]) 

        #calculation for line not closed
        else:
            line_closed=False      
            line=array(input_pts[1:]-input_pts[:-1])
            normal=transpose(array((append(-line[:,1], line[:,0])).reshape(2,size(line)/2)))
            normnormal=normal*(sum(normal**2,axis=1)**-.5).reshape(len(normal),1)          
            average=(normnormal[0:-1]+normnormal[1:])/2
            normaverage=average*(sum(average**2,axis=1)**-.5).reshape(len(average),1)
            # add normal of line for first and last point as average normal
            normaverage=append(append(normnormal[0],normaverage),normnormal[-1]).reshape(-1,2)
                
        # at least 3 points needed for calculation
        
        if len(input_pts)<3:
            print " Too few points need more to run discrete level_set (at least 3 points needed)"
        
        #    check for all three normals parallel
        else:
            check1 = where((normaverage[1:-1,0]==normaverage[2:,0])*(normaverage[2:,1]==normaverage[1:-1,1])
                            &(normaverage[2:,0]==normnormal[1:,0])*(normaverage[1:-1,1]==normnormal[1:,1]))[0]             
            check1=check1+1     # one has to be added because search is starting from point two on
            
            #delete input_pts, normed normal and average normal
            #recalculate line
            if len (check1)!=0:
                input_pts=delete(input_pts,check1,axis=0)
                line=array(input_pts[1:]-input_pts[:-1])            
                normaverage=delete(normaverage,check1,axis=0)
                normnormal=delete(normnormal,check1,axis=0)
                
        #    check for parallel average normal but not equivalent normal of line-> add point on line
            check2 = where((normaverage[:-1,0]==normaverage[1:,0])*(normaverage[1:,1]==normaverage[:-1,1])
                            &(normaverage[1:,0]!=normnormal[:,0])*(normaverage[:-1,1]!=normnormal[:,1]))[0]
            
            if len(check2)!=0:
                add_all_pts=[]
                add_all_normal=[]
            # add average of adjacent points of line to get new input points
            # first calculate all, and afterwards insert all, so we don`t stumble over indices
                for points in check2:
                    add_pt=(input_pts[points]+input_pts[points+1])/2
                    add_normal=(normnormal[points])
                    add_all_pts=append(add_all_pts,add_pt)
                    add_all_normal=append(add_all_normal,add_normal)
                
                input_pts=insert(input_pts,check2+1,add_all_pts.reshape(-1,2),axis=0)
                line=array(input_pts[1:]-input_pts[:-1])
                normaverage=insert(normaverage,check2+1,add_all_normal.reshape(-1,2),axis=0)
                normnormal=insert(normnormal,check2+1,add_all_normal.reshape(-1,2),axis=0)
            
            #save data
            
            self.input_pts=input_pts
            self.line=line
            self.length=sum(line**2,axis=1)**.5
            self.normnormal=normnormal
            self.normaverage=normaverage

            
    def get_value(self):        
        """
        calculates the signed distance for a discrete input_pointline
        three successive lines have to be  linear independent
        Because if so the average normal vectors of one line will be linear dependent with each other
        therefore the solution of their meeting point is not defined
        Also at least two lines are needed, for the same reason as above 
        It is recommended to use arrange_input_pts to check for correctness and delete unnecessary points
        If line of input_pts is not closed, level set of grid_pts outside of domain(=outside of edge normals), will be set to 0. 
        """

    # get variables
        work_grid=self.work_grid
        input_pts=self.input_pts
        
    #     set time    
        t_start = time()
        
    #     calculate variables: line, length, normal, average normal 
        if self.check_input==False:
            
            # calculation for line closed
            if input_pts[0].all()==input_pts[-1].all():
                line_closed=True
                input_pts=append(input_pts, input_pts[1]).reshape(-1,2)
                line=array(input_pts[1:]-input_pts[:-1])            
                normline=line*(sum(line**2,axis=1)**-.5).reshape(len(line),1)    
                normal=transpose(array((append(-line[:,1], line[:,0])).reshape(2,size(line)/2)))
                normnormal=normal*(sum(normal**2,axis=1)**-.5).reshape(len(normal),1)          
                average=(normnormal[:-1]+normnormal[1:])/2
                normaverage=average*(sum(average**2,axis=1)**-.5).reshape(len(average),1)
                normaverage=append(normaverage[-1],normaverage[:]).reshape(-1,2)
                # delete not needed input_pts which where only doubled to get normal average of connected point
                input_pts=delete(input_pts,-1,-2).reshape(-1,2)
                line=array(input_pts[1:]-input_pts[:-1]) 
    
            #calculation for line not closed
            else:
                line_closed=False      
                line=array(input_pts[1:]-input_pts[:-1])
                normline=line*(sum(line**2,axis=1)**-.5).reshape(len(line),1)    
                normal=transpose(array((append(-line[:,1], line[:,0])).reshape(2,size(line)/2)))
                normnormal=normal*(sum(normal**2,axis=1)**-.5).reshape(len(normal),1)          
                average=(normnormal[0:-1]+normnormal[1:])/2
                normaverage=average*(sum(average**2,axis=1)**-.5).reshape(len(average),1)
                # add normal of line for first and last point as average normal
                normaverage=append(append(normnormal[0],normaverage),normnormal[-1]).reshape(-1,2)
               
    # if preprocessing has been done just get it
        else:
            line=self.line
            normline=line*(sum(line**2,axis=1)**-.5).reshape(len(line),1)    
            length=self.length
            normnormal=self.normnormal
            normaverage=self.normaverage
        
    #    calculation of center_points, and angles where the average normals of each segment cross
    #    imput_pts & normaverage -> points (meeting point of average normals)
    #    imput_pts & points -> cosin(alpha) between meeting point and line points
        points=[] #center_points
        cosall=[]
        
        for i in range(0,len(normaverage)-1):
            Matrix=transpose(append(normaverage[i],-1*normaverage[i+1]).reshape(2,-1))
            value=solve(Matrix,line[i])
            point=input_pts[i]+normaverage[i]*value[0]
            cosalpha=vdot((input_pts[i]-point),(input_pts[i+1]-point))/(sum((input_pts[i]-point)**2,axis=0)*sum((input_pts[i+1]-point)**2,axis=0))**0.5
            points=append(points,point)
            cosall=append(cosall,cosalpha)             
        points=points.reshape(shape(line))
        
        #    vector from meeting point toward start point and end point of each line = vhelp
        #    will be needed later on for check if grid point lies within segment
 
        vhelp=[]
        for i in range(0,len(points)):
            vhelp=append(vhelp,[(input_pts[i]-points[i]),(input_pts[i+1]-points[i])])
        vhelp=vhelp.reshape(-1,2)
        normvtopoints=vhelp*(sum(vhelp**2,axis=1)**-.5).reshape(-1,1)

        
        """ start of algorithm, all following calculations are done over the hole domain of work_grid - points"""
 
        #    calculation of signed distance will only be evaluated, if point from work_grid lies within one segment
        result=zeros(len(work_grid)) #result array,where signeddistance will be saved in      
        for i in range(0,len(work_grid)): #for all work_grid points
            
            vtest=-points+work_grid[i]  #test vector
            normvtest=vtest*(sum(vtest**2,axis=1)**-.5).reshape(len(line),1)
            
            for j in range(0,len(normvtest)): # for all segments
                help1=vdot(normvtest[j],normvtopoints[2*j])     # cosin (alpha) between work_grid and point1 of line 
                help2=vdot(normvtest[j],normvtopoints[2*j+1])   # cosin (beta) between work_grid and point2 of line
                
                if min((help1,help2))>=cosall[j]:               # test if in between both vectors
                    
        #    check if projection on normal or minimum distance to points leads to minimum signed distance
                    check=vdot(normline[j],(work_grid[i]-input_pts[j])) 
                    if ((check<=length[j])* (check>=0.0))==True:
                        signeddistance=vdot(normnormal[j],(work_grid[i]-input_pts[j]))
                        
                    
                    else:
                        sig=sign(vdot(normnormal[j],(work_grid[i]-input_pts[j])))
                        length1=(sum((work_grid[i]-input_pts[j])**2,axis=0))**.5
                        length2=(sum((work_grid[i]-input_pts[j+1])**2,axis=0))**.5
                        signeddistance=sig*min(array([length1,length2]))
        
        #    some segments may still be overlapping and therefore absolute minimum value has to be saved 
        
                    if result[i]!=0:
                            result[i]=sign(signeddistance)*min(array([abs(signeddistance),abs(result[i])]))
                    else:
                        result[i]=signeddistance
                        
        print "Elapsed time:   %8.6f sec" % (time () - t_start) 
        self.scalars=array(result).reshape(shape(self.out_grid[0]))
    
    @show
    def show(self): 
        """Test contour_surf on regularly spaced co-ordinates like MayaVi."""
        x, y = self.out_grid
        s = surf(x, y, self.scalars, warp_scale = 0.5)
        colorbar(s)
        return s
 
 # not needed now....
     
    def get_segment(self,x,y):
        """evaluates segment, for equally spaced grid. In the following way:
            2    5    8
            1    4    7
            0    3    6
        points shouldn't lay on line otherwise inaccurate, because of round of errors.
        Hence otherwise point on outer line of 2 could become related to segment 3"""
        x_e=self.x_elements
        y_e=self.y_elements
        x_l=self.x_limit
        y_l=self.y_limit
        spacing_x=(x_l[1]-x_l[0])/x_e  
        spacing_y=(x_l[1]-x_l[0])/x_e 
        segment_y=y_e-1
        return int(int((x-x_l[0])/spacing_x)*segment_y+int((y-y_l[0])/spacing_y))

 # not needed now....
     
    def get_nodal_coordinates(self,number_of_segment):
        """ evaluate all 4 nodes, which are connected to segment n """
        
        x_e=self.x_elements
        y_e=self.y_elements
               
        node1=number_of_segment+int(number_of_segment/(x_e-1))
        node2=node1+1
        node3=int(node1+y_e)
        node4=int(node3+1)
        x=self.out_grid[0].flatten()
        y=self.out_grid[1].flatten()
        coordinates=array([x[node1],y[node1]],[x[node2],y[node2]],[x[node3],y[node3]],[x[node4],y[node4]])
        return coordinates
 
#test
if __name__ == '__main__':
    
    ls = LevelSet_discrete()
    ls.x_limit = array([-1.0,3.0])   
    ls.y_limit = array([-2.0,5.0])
    ls.x_elements = int(100.0)         #number of elements in x-direction used by set_grid, get_segment
    ls.y_elements = int(100.0)         #number of elements in y-direction used by set_grid, get_segment
    # steps
    ls.input_pts = array([[0.,0.],[1.0,0.0],[1.,1.0],[2.0,1.0],[2.0,2.0],[3.0,2.0],[3.0,3.0],[4.0,3.0],[4.0,4.0]])
    ls.set_grid()
    ls.check_input_pts()
    ls.get_value()
    ls.show()