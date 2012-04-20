#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Sep 8, 2011 by: matthias

from enthought.mayavi.core.api import PipelineBase
from enthought.mayavi.core.ui.api import MayaviScene, SceneEditor, \
    MlabSceneModel
from enthought.mayavi.modules.axes import Axes

from enthought.traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int, Bool, File, Array

from enthought.traits.ui.api import \
    View, Item, Group, ButtonEditor, RangeEditor, VGroup, HGroup, HSplit, Tabbed, \
    ViewSubElement, VGrid
from enthought.mayavi import mlab
from enthought.mayavi.core.api import Engine

import enthought.tvtk as tvtk

import tempfile
import os
import numpy as np
import string

# own Modules
from crease_pattern import CreasePattern

class CreasePatternView( HasTraits ):

    # data source
    data = Instance( CreasePattern )

    # plotting stuff
    scene = Instance( MlabSceneModel, () )
    plot = Instance( PipelineBase )
    scalefactor = 1.0

    # range of fold steps
    fold_step_min = Int( 0 )
    fold_step_max = Property
    def _get_fold_step_max( self ):
        return self.data.iteration_nodes.shape[0] - 1
    
    fold_step = Int( 0 )
    
    # constrain datas
    
    cnstr = None
    def get_cnstr( self ):
        
        #=======================================================================
        # This Method get the constrain-information from the actual crease-pattern
        # and divides it for easier calculation of the symbol-positions
        # 
        # The constrains are devided in three constrain-types:
        # - fixed constrains ( constrain is full fixed in his direction)
        # - connected constrains ( constrains are in an mathmatical abhaengigkeit
        #                        , e.g. constant or linear movementbehavior)
        # - load constrains ( constrains, which activates the numerical calculation)
        #
        # Different list for visualation:
        # fixed constrains:
        # - cn_f : indexarray of nodes of fixed constrains
        # - cd_f : direction (axes) of the fixed constrain [x, y, z]
        #
        # connected constrains:
        # - cn_c: indexarray of nodes of fixed constrains
        # - cc_c: connection of nodes as indexarrays, e.g. [[0, 1],
        #                                                   [2, 4]]
        #            each index represents a node
        # - cd_c: direction (axes) of the connected constrains [x, y, z]
        #
        # load constrains:
        # - cn_l : indexarray of nodes of load constrains
        # - cd_l : direction (axes) of the load constrain [x, y, z]
        #
        #=======================================================================
        
        # get constrain information of the creasepattern
            
        lhs = self.data.cnstr_lhs
        rhs = self.data.cnstr_rhs
        
        # get load constrains
                
        load_nodes = np.array( [] , dtype = int )
        load_dir = np.array( [] )
        count = 0
        while( count < len( rhs ) ):
            # all cell's in rhs which are different 0 represents a load and 
            # gives the direction 
            # the constrain on the same indexposition in lhs is the load constrain
            if ( rhs[count] != 0 ):
                node = lhs[count][0][0]
                dir_vec = np.array( [0, 0, 0] )
                dir_vec[lhs[count][0][1]] = 1
                
                if( rhs[count] < 0 ):
                    dir_vec *= -1
                    
                load_nodes = np.append( load_nodes, node )
                load_dir = np.append( load_dir, [dir_vec] )
                lhs.remove( lhs[count] ) # remove the constrain from lhs-list
                
            count += 1
        
        load_nodes = load_nodes.reshape( len( load_nodes ), 1 )
        load_dir = load_dir.reshape( len( load_nodes ), 3 )
        
        # divide all other constrains to fixed or connected constrains
        
        cnstr_fixed = lhs
        cnstr_connect = []
            
        count = 0
        while( count < len( cnstr_fixed ) ):
            # put all connected constrains out of cnstr_fixed into cnstr_connect    
            if len( cnstr_fixed[count] ) > 1:
                cnstr_connect.append( cnstr_fixed.pop( count ) )
                continue
            
            count += 1
            
        # create Cnstr Arrays
        
        fixed_nodes = np.array( [] , dtype = int )
        fixed_dir = np.array( [] )
        connect_nodes = np.array( [], dtype = int )
        connect_dir = np.array( [] )
        connect_connections = np.array( [], dtype = np.int )
        
        # build cn_f and cd_f
        
        for i in cnstr_fixed:
            fixed_nodes = np.append( [fixed_nodes], [i[0][0]] )
            dir_vec = np.array( [0, 0, 0] )
            dir_vec[i[0][1]] = 1
            fixed_dir = np.append( [fixed_dir], [dir_vec] )
        
        fixed_nodes = fixed_nodes.reshape( len( cnstr_fixed ), 1 )
        fixed_dir = fixed_dir.reshape( len( cnstr_fixed ), 3 ) 
        
        # get connections on reale node indexes  
        for i in cnstr_connect:
            con = np.array( [i[0][0], i[1][0]] )
            connect_nodes = np.append( connect_nodes, con )
       
        connect_nodes = connect_nodes.reshape( len( cnstr_connect ), 2 )

        # build an mask showing the nodes which are connected
        mask = np.ones( len( self.data.nodes ), dtype = np.int )
        mask[connect_nodes[:, 0]] = 0
        mask[connect_nodes[:, 1]] = 0
            
        # build cn_c and translate cc_c from global indexing to an array 
        # with index the effective cn_c array
            
        cc_c = np.array( connect_nodes )
        connect_nodes = np.array( [], dtype = np.int )
        countc = 0
        for i in range( len( mask ) ):
            if( mask[i] == 0 ):
                cc_c[cc_c == i] = countc
                connect_nodes = np.append( connect_nodes, i )
                countc += 1
        
        # build cd_c for all nodes in creasepattern
                
        connect_dir = np.zeros( ( len( mask ), 3 ) )
        for i in cnstr_connect:
            connect_dir[i[0][0]][i[0][1]] = 1
            connect_dir[i[1][0]][i[1][1]] = 1
        
        count = 0
        
        # delete all direction-array which are [0, 0, 0] and won't represent 
        # a connected constrain
        for i in range( len( mask ) ):
            if ( mask[i] == 1 ):
                connect_dir = np.delete( connect_dir, count, 0 )
            else:
                count += 1
        # return ( cn_f, cd_f, cn_c, cc_c, cd_c, cn_l, cn_d)
        return ( fixed_nodes, fixed_dir, connect_nodes, cc_c, connect_dir,
                load_nodes, load_dir )

    show_cnstr = Bool( True )

    # When the scene is activated
    @on_trait_change( 'scene.activated' )
    def create_plot( self ):
        
        # get the actual constrain information
        
        self.cnstr = self.get_cnstr()
        cn_f, cd_f, cn_c, cc_c, cd_c, cn_l, cd_l = self.cnstr
        
        nodes = self.data.nodes
        
        # setup an functional camera position
        fx, fy = np.average( nodes[:, ( 0, 1 )], axis = 0 )
        fz_arr = self.data.iteration_nodes[:, :, 2]
        fz_min = np.min( fz_arr )
        fz_max = np.max( fz_arr )
        fz = ( fz_min + fz_max ) / 2.0
        print 'fz', fz

        # Arrays of Point Data in axis
        x = nodes[:, 0]
        y = nodes[:, 1]
        z = nodes[:, 2]

        # scalar for Points
        scalar = np.ones( z.shape[0] )

        # Visualize the data 

        fig = self.scene.mlab.gcf()
        self.scene.mlab.figure( fig, fgcolor = ( 0, 0, 0 ),
                               bgcolor = ( 1, 1, 1 ) )

        # 3D Points
        #pts = self.scene.mlab.points3d(x, y, z, scalar,
        #                               scale_factor = 0.04, resolution = 10)

        # connections from creaseline Array
        #pts.mlab_source.dataset.lines = self.data.crease_lines

        # draw lines between nodes
        #tube = self.scene.mlab.pipeline.tube(pts, tube_radius = 0.015)

        triangles = self.scene.mlab.triangular_mesh( x, y, z, self.data.facets )
        triangles.mlab_source.dataset.lines = self.data.crease_lines

        tube = self.scene.mlab.pipeline.tube( triangles, tube_radius = 0.015 )
        #pts = self.scene.mlab.pipeline.points3d(x, y, z, scalar,
        #                               scale_factor = 0.04, resolution = 10)

        #triangles = self.scene.mlab.triangular_mesh(x, y, z, self.data.facets)
        #triangles = self.scene.mlab.pipeline.delaunay2d(pts)

        self.scene.mlab.pipeline.surface( tube, color = ( 1.0, 1.0, 0.9 ) )
        self.scene.mlab.pipeline.surface( triangles, color = ( 0.6, 0.6, 0.6 ) )

#        # connections from creaseline Array

        if self.show_cnstr:
            # old constrain visualiton
#            cnstr = self.data.get_cnstr_pos()
#            self.quiver3d = self.scene.mlab.quiver3d(*cnstr.T, color = (0.0, 0.0, 1.0))
            
            # autoscale on the min lenght of a creaseline
            minLenght = np.amin( self.data.c_lengths )
            self.scalefactor = 0.5 * minLenght
            # spacefactor is giving space between constrains an real node position
            spacefactor = 0.2 * self.scalefactor
            
            # fixed cnstr
            cp_f = np.array( [] )
            cp_f = nodes[cn_f]
            x, y, z = cp_f.T
            x, y, z = x[0], y[0], z[0]
            U, V, W = cd_f.T * self.scalefactor
            sU, sV, sW = cd_f.T * spacefactor
            
            x = x - U - sU
            y = y - V - sV
            z = z - W - sW
            
            self.cf_arrow = self.scene.mlab.quiver3d( x, y, z, U, V, W, mode = '2darrow', color = ( 0.0, 0.0, 1.0 ), scale_mode = 'vector', scale_factor = 1.0, line_width = self.scalefactor*0.5 )
            self.cf_cross = self.scene.mlab.quiver3d( x, y, z, U, V, W, mode = '2dcross', color = ( 0.0, 0.0, 1.0 ), scale_mode = 'vector', scale_factor = 1.0, line_width = self.scalefactor*0.5 )  
            
            self.scene.mlab.pipeline.surface( self.cf_cross )
            self.scene.mlab.pipeline.surface( self.cf_arrow )
            
            # alternative
#            pyramide = self.scene.mlab.quiver3d(x, y, z, U, V, W, mode = 'cone', resolution = 4,  scale_mode = 'vector', scale_factor = 1.0, color = (0.0, 0.0, 1.0) )
#            self.scene.mlab.pipeline.surface(pyramide)
 
            
            # load constrain
            
            cp_l = nodes[cn_l]
            
            x, y, z = cp_l.T
            x, y, z = x[0], y[0], z[0]
            U, V, W = cd_l.T * self.scalefactor
            sU, sV, sW = cd_l.T * spacefactor
            
            x = x - U - sU
            y = y - V - sV
            z = z - W - sW
            
            self.cl_arrow = self.scene.mlab.quiver3d( x, y, z, U, V, W, mode = 'arrow', color = ( 1.0, 0.0, 0.0 ), scale_mode = 'vector', scale_factor = 1.0 )
            self.scene.mlab.pipeline.surface( self.cl_arrow )
            
            # connected contrains
              
            cp_c = nodes[cn_c]
            
            x, y, z = cp_c.T
            
            U, V, W = cd_c.T * self.scalefactor
            sU, sV, sW = cd_c.T * spacefactor
            
            x = x - U - sU
            y = y - V - sV
            z = z - W - sW
            
            # position of connection line, a little under the constrain arrow
            #xL = x - sU
            #yL = y - sV
            #zL = z - sW
            
            
            self.cc_arrow = self.scene.mlab.quiver3d( x, y, z, U, V, W, mode = '2darrow', line_width = self.scalefactor*0.5, color = ( 0.0, 1.0, 0.0 ), scale_mode = 'vector', scale_factor = 1.0 )
            
            self.cc_arrow.mlab_source.dataset.lines = cc_c
            
            
            self.scene.mlab.pipeline.surface( self.cc_arrow, color = (0.0, 0.7, 0.0), line_width = self.scalefactor*0.5 )
            
            #self.cc_pts = self.scene.mlab.points3d( xL, yL, zL, mode = '2dvertex', color = ( 0.0, 1.0, 0.0 ), line_width = 0.1 )
            #self.cc_pts.mlab_source.dataset.lines = cc_c
            #self.scene.mlab.pipeline.surface( self.cc_pts, line_width = 0.01, color = ( 0.0, 0.5, 0.0 ) )
            
            # Draw cnstr (line/plane)
            # Outlines
            all_nodes = self.data.iteration_nodes
            Xmin = all_nodes[0][:,0].min()
            Ymin = all_nodes[0][:,1].min()
            Zmin = all_nodes[0][:,2].min()
            Xmax = all_nodes[0][:,0].max()
            Ymax = all_nodes[0][:,1].max()
            Zmax = all_nodes[0][:,2].max()
            for i in range(len(all_nodes)):
                tXmin =  all_nodes[i][:,0].min()
                tYmin =  all_nodes[i][:,1].min()
                tZmin =  all_nodes[i][:,2].min()
                tXmax =  all_nodes[i][:,0].max()
                tYmax =  all_nodes[i][:,1].max()
                tZmax =  all_nodes[i][:,2].max()
                if(Xmin > tXmin):
                    Xmin = tXmin
                if(Ymin > tYmin):
                    Ymin = tYmin
                if(Zmin > tZmin):
                    Zmin = tZmin
                if(Xmax < tXmax):
                    Xmax = tXmax
                if(Ymax < tYmax):
                    Ymax = tYmax
                if(Zmax < tZmax):
                    Zmax = tZmax
            
            
            print Xmin, Ymin, Zmin, Xmax, Ymax, Zmax
            
            hX = Xmax - Xmin
            hY = Ymax - Ymin
            hZ = Zmax - Zmin
            
            Xmin -= hX/2
            Ymin -= hY/2
            Zmin -= hZ/2
            Xmax += hX/2
            Ymax += hY/2
            Zmax += hZ/2   
            
            outlines = np.array([Xmin, Ymin, Zmin, Xmax, Ymax, Zmax], dtype = float)
            
            # calc Points
            pts_l = self.calc_points_line(self.data.cnstr[0],outlines)
            
            print pts_l
            self.line = self.scene.mlab.points3d(pts_l[:,0], pts_l[:,1], pts_l[:,2], mode = '2dvertex')
            self.line.mlab_source.dataset.lines = np.array([[0,1]])
            self.scene.mlab.pipeline.surface(self.line, line_width = 1, color = (1.0, 0.0, 0.0))
             

        self.scene.mlab.view( -60.0, 70.0, focalpoint = [fx, fy, fz] )


        self.plot = triangles


    
    def calc_points_line(self, cnstr_line, outlines):
        pts = np.array([], dtype =float)
        
        
        # Plane 1
        r = (outlines[0] - cnstr_line.base_v[0])/cnstr_line.a_dir_v[0]
        y = cnstr_line.base_v[1] + r*cnstr_line.a_dir_v[1]
        if( y < outlines[4] and y > outlines[1]):
            z = cnstr_line.base_v[2] + r*cnstr_line.a_dir_v[2]
            if( z < outlines[5] and z > outlines[2]):
                pts = np.append(pts, np.array([outlines[0],y,z]))
        
        # Plane 2
        r = (outlines[1] - cnstr_line.base_v[1])/cnstr_line.a_dir_v[1]
        x = cnstr_line.base_v[0] + r*cnstr_line.a_dir_v[0]
        if( x < outlines[3] and x > outlines[0]):
            z = cnstr_line.base_v[2] + r*cnstr_line.a_dir_v[2]
            if( z < outlines[5] and z > outlines[2]):
                pts = np.append(pts, np.array([x,outlines[1],z]))
                
        # Plane 3
        r = (outlines[2] - cnstr_line.base_v[2])/cnstr_line.a_dir_v[2]
        x = cnstr_line.base_v[0] + r*cnstr_line.a_dir_v[0]
        if( x < outlines[3] and x > outlines[0]):
            y = cnstr_line.base_v[1] + r*cnstr_line.a_dir_v[1]
            if( y < outlines[4] and y > outlines[1]):
                pts = np.append(pts, np.array([x,y,outlines[3]]))
                
        # Plane 4
        r = (outlines[3] - cnstr_line.base_v[0])/cnstr_line.a_dir_v[0]
        y = cnstr_line.base_v[1] + r*cnstr_line.a_dir_v[1]
        if( y < outlines[4] and y > outlines[1]):
            z = cnstr_line.base_v[2] + r*cnstr_line.a_dir_v[2]
            if( z < outlines[5] and z > outlines[2]):
                pts = np.append(pts, np.array([outlines[3],y,z]))
        
        # Plane 5
        r = (outlines[4] - cnstr_line.base_v[1])/cnstr_line.a_dir_v[1]
        x = cnstr_line.base_v[0] + r*cnstr_line.a_dir_v[0]
        if( x < outlines[3] and x > outlines[0]):
            z = cnstr_line.base_v[2] + r*cnstr_line.a_dir_v[2]
            if( z < outlines[5] and z > outlines[2]):
                pts = np.append(pts, np.array([x,outlines[4],z]))
                
        # Plane 6
        r = (outlines[5] - cnstr_line.base_v[2])/cnstr_line.a_dir_v[2]
        x = cnstr_line.base_v[0] + r*cnstr_line.a_dir_v[0]
        if( x < outlines[3] and x > outlines[0]):
            y = cnstr_line.base_v[1] + r*cnstr_line.a_dir_v[1]
            if( y < outlines[4] and y > outlines[1]):
                pts = np.append(pts, np.array([x,y,outlines[5]]))
        
        pts = pts.reshape(len(pts)/3, 3)
        return pts
    
    # when parameters are changed,
    # plot is updated
    @on_trait_change( 'fold_step, show_cnstr' )
    def update_plot( self ):

        # Array of current foldstep
        nodes = self.data.iteration_nodes[self.fold_step]
        x, y, z = nodes.T

        # Visualize the data 

        # set new position of 3D Points
        self.plot.mlab_source.set( x = x, y = y, z = z )

        if self.show_cnstr:
            # update constrain symbols
            
            cn_f, cd_f, cn_c, cc_c, cd_c, cn_l, cd_l = self.cnstr
            
            
            spacefactor = 0.2 * self.scalefactor
            
            # fixed cnstr
            cp_f = np.array( [] )
            cp_f = nodes[cn_f]
            x, y, z = cp_f.T
            x, y, z = x[0], y[0], z[0]
            U, V, W = cd_f.T * self.scalefactor
            sU, sV, sW = cd_f.T * spacefactor
            
            x = x - U - sU
            y = y - V - sV
            z = z - W - sW
                       
            self.cf_arrow.mlab_source.set( x = x, y = y, z = z )
            self.cf_cross.mlab_source.set( x = x, y = y, z = z )
            
            #load constrains
            
            cp_l = nodes[cn_l]
            
            x, y, z = cp_l.T
            x, y, z = x[0], y[0], z[0]
            U, V, W = cd_l.T * self.scalefactor
            sU, sV, sW = cd_l.T * spacefactor
            
            x = x - U - sU
            y = y - V - sV
            z = z - W - sW
            
            self.cl_arrow.mlab_source.set( x = x, y = y, z = z )
            
            # connected constrains
            
            cp_c = nodes[cn_c]
            
            x, y, z = cp_c.T
            
            U, V, W = cd_c.T * self.scalefactor
            sU, sV, sW = cd_c.T * spacefactor
                        
            x = x - U - sU
            y = y - V - sV
            z = z - W - sW
            
            #xL = x - sU
            #yL = y - sV
            #zL = z - sW
            
            # ToDo: cc_pts creates an vtk error on every FIRST update.
            #       this error won't treat the program behavior 
            #self.cc_pts.mlab_source.set( x = xL, y = yL, z = zL )
            self.cc_arrow.mlab_source.set( x = x, y = y, z = z )
            

    save_animation = Button
    animation_file = File
    def _animation_file_default( self ):
        return os.path.join( 'fig', 'bearigami.gif' )

    def _save_animation_fired( self ):

        #===========================================================================
        # Prepare plotting 
        #===========================================================================
        tdir = tempfile.mkdtemp()
        n_steps = len( self.data.iteration_nodes )

        steps_forward = range( n_steps )
        steps_backward = range( n_steps, 2 * n_steps )
        fnames_forward = [os.path.join( tdir, 'x%02d.jpg' % i )
                          for i in steps_forward ]
        fnames_backward = [os.path.join( tdir, 'x%02d.jpg' % i )
                           for i in steps_backward ]

        nodes_history = self.data.iteration_nodes
        for nodes, fname in zip( nodes_history, fnames_forward ):
            # Array of current foldstep
            x, y, z = nodes.T
            self.plot.mlab_source.set( x = x, y = y, z = z )

            if self.show_cnstr:
                # set new position of constraints
                cnstr = self.data.get_cnstr_pos( nodes )
                x, y, z = cnstr.T[:3]
                self.quiver3d.mlab_source.set( x = x, y = y, z = z )

            self.scene.mlab.savefig( fname, size = ( 300, 200 ) )

        for nodes, fname in zip( nodes_history[-1::-1], fnames_backward ):
            # Array of current foldstep
            x, y, z = nodes.T
            self.plot.mlab_source.set( x = x, y = y, z = z )

            if self.show_cnstr:
                # set new position of constraints
                cnstr = self.data.get_cnstr_pos( nodes )
                x, y, z = cnstr.T[:3]
                self.quiver3d.mlab_source.set( x = x, y = y, z = z )

            self.scene.mlab.savefig( fname, size = ( 300, 200 ) )

        fnames = fnames_forward + fnames_backward
        images = string.join( fnames, ' ' )
        destination = self.animation_file

        import platform
        if platform.system() == 'Linux':
            os.system( 'convert ' + images + ' ' + destination )
        else:
            raise NotImplementedError, 'film production available only on linux'
        print 'animation saved in', destination

    # The layout of the dialog created
    view = View( 
                HSplit( Group( 
                             Group( Item( 'show_cnstr' ), ),
                             Group( Item( 'save_animation', show_label = False ),
                                    Item( 'animation_file', show_label = False ),
                                    ),
                             id = 'creasepatternview.animation',
                             dock = 'tab'
                             ),
                              
                      VGroup( 
                             Item( 'scene', editor = SceneEditor( scene_class = MayaviScene ),
                                  show_label = False ),
                             Item( 'fold_step', editor = RangeEditor( low_name = 'fold_step_min',
                                                         high_name = 'fold_step_max',
                                                         format = '(%s)',
                                                         auto_set = False,
                                                         enter_set = False,
                                                         ),
                                    show_label = False
                                ),
                             id = 'creasepatternview.mayavi',
                             dock = 'tab'
                             ),
                      
                       
                        ),
                      
                dock = 'tab',
                resizable = True,
                id = 'creaspatternview',
                width = 1.0,
                height = 1.0
                )

if __name__ == '__main__':

    # initialise CreasePattern
    cp = CreasePattern()

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ],
                [ 1, 1, 0 ]]

    cp.crease_lines = [[ 0, 1 ], [1, 2]] # , [2, 0]]

    X = np.zeros( ( cp.n_dofs, ), dtype = float )

    cp.set_next_node( X )

    # initialise View
    my_model = CreasePatternView( data = cp )

    my_model.configure_traits()
