'''
Created on Jun 1, 2010

@author: alexander
'''
    data_array_ironed = Property( Array( float ),
                                  depends_on = 'data_array, +ironing_param, +axis_selection' )
    @cached_property
    def _get_data_array_ironed(self):
        '''remove the jumps in the displacement curves 
        due to resetting the displacement gauges. 
        '''
        print '*** curve ironing activated ***'
        
        # each column from the data array corresponds to a measured parameter 
        # e.g. displacement at a given point as function of time u = f(t))
        #
        data_array_ironed = copy( self.data_array )
        
        for idx in range( self.data_array.shape[1] ):
            
            # use ironing method only for columns of the displacement gauges.
            #
            if self.names_and_units[0][ idx ] != 'Kraft' and \
                self.names_and_units[0][ idx ] != 'Bezugskanal' and \
                self.names_and_units[0][ idx ] != 'Weg':
            
                # 1d-array corresponding to column in data_array
                data_arr = copy( data_array_ironed[:,idx] )

                # get the difference between each point and its successor
                jump_arr =  data_arr[1:] - data_arr[0:-1]
                
                #--------------
                
                jump_arr_idx = argsort( jump_arr )
                
                count_jumps = 0

                data_arr_range = max( data_arr ) - min( data_arr )
                print 'data_arr_range START', data_arr_range
                print 'jump_arr_idx', jump_arr_idx
                print 'jump_arr', jump_arr[ jump_arr_idx ][0:10]
                for jidx in jump_arr_idx:
                    shift = jump_arr[ jidx ]
                    print 'shift', shift
                    data_arr_trial = copy( data_arr )
                    data_arr_trial[ jidx+1: ] -= shift
                    # get the new range of the measured data 
                    data_arr_trial_range = max( data_arr_trial ) - min( data_arr_trial )
                    print 'data_arr_trial_range', data_arr_trial_range
                    print 'ratio', ( data_arr_trial_range - data_arr_range ) / data_arr_range
                    if abs( ( data_arr_trial_range - data_arr_range ) / data_arr_range ) > self.jump_rtol:
                        print 'in curve ironing: true jump received' 
                        data_arr = data_arr_trial 
                        data_arr_range = copy( data_arr_trial_range ) 
                        print 'data_arr_range', data_arr_range
                        count_jumps += 1
                    else:
                        print 'BREAK'
                        break
                    
                print 'number of jumps removed in data_array for column', self.names_and_units[0][ idx ], ': ', count_jumps
                
                #---------------
                 
#                # get the range of the measured data 
#                data_arr_range = max( data_arr ) - min( data_arr )
#
#                # determine the relevant criteria for a jump
#                # based on the data range and the specified tolerances:
#                jump_crit = self.jump_rtol * data_arr_range  
#    
#                # get the indexes in 'data_column' after which a 
#                # jump exceeds the defined tolerance criteria
#                jump_indexes = where( fabs(jump_arr) > jump_crit )
#                
#                print 'number of jumps removed in data_arr_ironed: ', jump_indexes[0].shape[0]
#                # glue the curve at each jump together
#                for n in range(jump_indexes[0].shape[0]):
#                    # get the offsets at each jump of the curve
#                    jidx = jump_indexes[0][n]
#                    shift = data_arr[jidx+1] - data_arr[jidx]
#                    # shift all succeeding values by the calculated offset
#                    data_arr[jidx+1:] -= shift

                data_array_ironed[:,idx] = data_arr[:]

        return data_array_ironed 