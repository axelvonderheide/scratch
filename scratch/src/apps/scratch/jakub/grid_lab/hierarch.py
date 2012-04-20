'''
Created on May 1, 2009

@author: jakub
'''
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, \
     This, Self

from numpy import array, arange

class FeDomainList(HasTraits):
    
    _domains = List
    domains = Property(List)
    
    def _set_domains(self, list):
        for i, domain in enumerate(list):
            domain.d_list = self
            if i == 0:
                domain.previous_domain = FERefinementDomain()# this should be replaced by none
            else:
                domain.previous_domain = list[i-1]
        self._domains = list
    
    def _get_domains(self):
        return self._domains
    
    def add_domain(self, domain):
        domain.d_list = self #have to be here
        domain.previous_domain = self.domains[-1]
        self.domains.append(domain)


class FERefinementDomain(HasTraits):
    
    elements = List(Int)
    d_list = Instance(FeDomainList)
    previous_domain = This
    children = List(This) #list of the instances of the same class
    
    def add_child(self, c_domain):
        self.d_list.add_domain(c_domain)
        self.children.append(c_domain)
    
    ldof_map = Array(Int, value=[0,1,2,3])
    dof_offset = Property(Int, depends_on = 'previous_domain.dof_offset, previous_domain.n_dofs')
    @cached_property
    def _get_dof_offset(self):
        if self.previous_domain:
            return self.previous_domain.dof_offset + self.previous_domain.n_dofs
        else:
            return 0
    
    n_dofs = Property(Int, depends_on = 'ldof_map')
    @cached_property
    def _get_n_dofs(self):
        return self.ldof_map.shape[0]
    
    def get_dof_map(self):
        return self.ldof_map + self.dof_offset
    
class FEEnrDomain(HasTraits):
    
    pass

if __name__ == '__main__':
    
    d_list = FeDomainList()
    d1 = FERefinementDomain()
    d2 = FERefinementDomain()
    d3 = FERefinementDomain()
    d4 = FERefinementDomain() 
    
    d_list.domains = [d1,d2,d3] 
    d_list.add_domain(d4)
    
    print "ddm ", d4.get_dof_map()
    print "do ", d4.dof_offset
    
    d2.ldof_map = arange(8)
    
    print "ddm ", d4.get_dof_map()
    print "do ", d4.dof_offset
    
    c1 = FERefinementDomain()
    d1.add_child(c1)
    print "cdm ", c1.get_dof_map()
    print "co ", c1.dof_offset
    
    d2.ldof_map = arange(12)
    
    print "cdm ", c1.get_dof_map()
    print "co ", c1.dof_offset
#    
#    c2 = FERefinementDomain()
#    c1.add_child(c2)
#    print "cdm ", c2.get_dof_map()
#    print "co ", c2.dof_offset
#    
#    c3 = FERefinementDomain()
#    c1.add_child(c3)
#    print "cdm ", c3.get_dof_map()
#    print "co ", c3.dof_offset