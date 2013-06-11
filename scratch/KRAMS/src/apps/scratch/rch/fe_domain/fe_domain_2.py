'''
Created on May 1, 2009

'''
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, \
     This, self, TraitError

from numpy import array, arange

class FEDomain(HasTraits):
    '''Test the state dependencies within the hierarchical domain representation.
    '''
    subdomains = List( domain_changed = True )
    @on_trait_change('subdomains[]')
    def _validate_subdomains(self):
        for domain in self.subdomains:
            if domain.parent != None:
                raise ValueError, 'only parentless subdomains can be inserted into domain'
    
    serialized_subdomains = List
    
    def _append_in_series( self, new_subdomain ):
        '''Link the new subdomain at the end of the series.
        '''
        if self.serialized_subdomains:
            last_subdomain = self.serialized_subdomains[-1]
            last_subdomain.next_domain = new_subdomain
            new_subdomain.previous_domain  = last_subdomain
        self.serialized_subdomains.append( new_subdomain )
        
    n_dofs = Property
    def _get_n_dofs(self):
        '''Return the total number of dofs in the domain.
        Use the last subdomain's: dof_offset + n_dofs 
        '''
        last_subdomain = self.serialized_subdomains[-1]
        return last_subdomain.dof_offset + last_subdomain.n_dofs

class FESubDomain(HasTraits):

    name = Str('<unnamed>')

    # local dof enumeration
    n_dofs = Int( 4, domain_changed = True )
    
    # dof offset within the global enumeration
    dof_offset = Property( Int, depends_on = 'previous_domain.dof_offset' )
    #cached_property
    def _get_dof_offset(self):
        if self.previous_domain:
            return self.previous_domain.dof_offset + self.previous_domain.n_dofs
        else:
            return 0

    def __repr__(self):
        if self.previous_domain:
            return self.name + ' <- ' + self.previous_domain.name
        else:
            return self.name
    
    # delete - this example does not need it
    elements = List(Int)
    
    # dependency link for sequential enumeration
    previous_domain = This( domain_changed = True )
    @on_trait_change('previous_domain')
    def _validate_previous_domain(self):
        if self.previous_domain == self:
            raise TraitError, 'cyclic reference for ' + self.name
    
    # dependency link for sequential enumeration
    next_domain = This( domain_changed = True )
    @on_trait_change('next_domain')
    def _validate_previous_domain(self):
        if self.next_domain == self:
            raise TraitError, 'cyclic reference for ' + self.name
    
class FERefinementLevel( HasTraits ):

    # container domain (why is this needed?)
    _domain = Instance(FEDomain)
    domain = Property
    def _set_domain( self, value ):
        'reset the domain of this domain'
        if self.parent != None:
            raise TraitError, 'child FESubDomain cannot be added to FEDomain'
        if self._domain:
            # unregister in the old domain
            raise NotImplementedError, 'FESubDomain cannot be relinked to another FEDomain'      

        self._domain = value
        # register in the domain as a subdomain
        self._domain.subdomains.append( self )
        self._domain._append_in_series( self )        
    def _get_domain(self):
        if self.parent != None:
            return self.parent.domain
        else:
            return self._domain
    
    # children domains: list of the instances of the same class
    children = List(This) 
    
    # parent domain
    _parent = This( domain_changed = True )
    parent = Property( This )
    def _set_parent( self, value ):
        'reset the parent of this domain'
        if self._parent:
            # check to see that the changed parent 
            # is within the same domain
            if value.domain != self._parent.domain:
                raise NotImplementedError, 'Parent change across domains not implemented'
            # unregister in the old parent
            self._parent.children.remove( self )
        else:
            # append the current subdomain at the end of the subdomain
            # series within the domain
            value.domain._append_in_series( self )
        # set the new parent
        self._parent = value
        # register the subdomain in the new parent
        self._parent.children.append( self )
    def _get_parent(self):
        return self._parent

class FERefinementGrid( FESubDomain, FERefinementLevel ):

    ldof_map = Property
    def _get_ldof_map(self):
        return arange( self.n_dofs )
    
    def get_dof_map(self):
        return self.ldof_map + self.dof_offset


class FEEnrDomain(HasTraits):
    
    pass

if __name__ == '__main__':

    domain1 = FEDomain()
        
    d1 = FERefinementGrid( name = 'd1', domain = domain1 )
    d2 = FERefinementGrid( name = 'd2', parent = d1 )
    d3 = FERefinementGrid( name = 'd3', parent = d2 )
    d4 = FERefinementGrid( name = 'd4', parent = d3 ) 
    d5 = FERefinementGrid( name = 'd5', domain = domain1 ) 

    print 'serialized domains'
    print domain1.serialized_subdomains
    # initial enumeration of the last domain
    print d4, d4.get_dof_map()
    # change the number of DOFs
    d2.n_dofs = 8
    # check that the domain gets larger
    print d2, d2.get_dof_map()
    # check the enumeration of some next domain
    print d4, d4.get_dof_map()
    # change the d2 once again
    d2.n_dofs = 12
    # did d4 change as well?
    print d4, d4.get_dof_map()
    # add new subdomain
    c1 = FERefinementGrid( name = 'c1', parent = d1 )
    # look at the serialized domains - is the c1 
    # domain correctly appended?
    print domain1.serialized_subdomains
    # check the enumeration of the new domain
    print c1, c1.get_dof_map()
    # what is now the enumeration of d3?
    print d3, d3.get_dof_map()
    # enlarge d2 once again 
    d2.n_dofs = 16
    # d3 should have been shifted
    print d3, d3.get_dof_map()
    # and c1 as well
    print c1, c1.get_dof_map()
    # change the parent of d4 from d3 to d2
    d4.parent = d2
    # look at the children of d2
    print 'd2 children', d2.children
    # look at the serialization
    print domain1.serialized_subdomains
    # show the enumeration of d4
    print d4, d4.get_dof_map()
    # construct a second domain
    domain2 = FEDomain()
    # try the exception - subdomain cannot have 
    # explicit link to a domain
    #d2.domain = domain2 
    # move d1 to domain2 !!! raises an exception
    #d1.domain = domain2
    # show the subdomains of domain1
    print 'subdomains of domain1', domain1.subdomains
    # show the serialized subdomains of domain1
    print 'serialized subdomains of domain1', domain1.serialized_subdomains
    print 'serialized subdomains of domain2', domain2.serialized_subdomains
    
    print 'n_dofs', domain1.n_dofs