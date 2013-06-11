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
# Created on Apr 20, 2010 by: rch

from has_traits_orm import HasTraitsORM

from enthought.traits.api import Str, HasTraits, Int, DelegatesTo, Property, \
    Enum, Range, on_trait_change

from sqlalchemy import Column, Integer, String, ForeignKey
from sqlalchemy.ext.declarative import _as_declarative, declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.interfaces import MapperProperty, InstrumentationManager
from sqlalchemy.schema import MetaData
import sqlalchemy as sa



class Hands(HasTraitsORM):
    
    user_id     = Int(None, sqldb=True, sqlpk=True, col_args=(ForeignKey('users.id'),))
    
    l_fingers     = Range(low = 0, high = 5, value = 5, sqldb=True)
    r_fingers     = Range(low = 0, high = 5, value = 5, sqldb=True)
    
    user        = sa.orm.relation('User', uselist=False)
    name        = DelegatesTo('user')
    
    def __init__(self, **traits):
        super(Hands, self).__init__(**traits)        
        self.orm_init()
    
    @sa.orm.reconstructor
    def orm_init(self):
        if self.l_fingers is None:
            self.l_fingers = self.trait('l_fingers').default
        
        if self.r_fingers is None:
            self.r_fingers = self.trait('r_fingers').default
        
    def __repr__(self):
        return "< Hands ('%s has %s fingers.') >" %\
             (self.name, self.l_fingers + self.r_fingers)
            

class User(HasTraitsORM):
 
    # If left out will default to lower case class name. 
    __tablename__ = 'users'   

    # Default value of autoinc fields needs to be set to None so SA knows to 
    # fetch it from a sequence via SQL.
    id = Int(None, sqldb=True, sqlpk=True)
    name = Str(sqldb=True)
    lastname = Str(sqldb=True)
    password = Str(sqldb=True)

    hands  = sa.orm.relation(Hands, uselist=False)
    
    fingers = Property(depends_on=['hands.[l_fingers, hands.r_fingers]'])
    
    fullname = Property(Str)
    
    def __init__(self, **traits):
        super(User, self).__init__(**traits)
        self.orm_init()
    
    @sa.orm.reconstructor
    def orm_init(self):
        if self.hands is None:
            self.hands = Hands(user=self)

    def _get_fullname(self):
        return '%s %s' % (self.name, self.lastname)
    
    
    def _get_fingers(self):
        return self.hands.l_fingers + self.hands.r_fingers
        
        
    def _name_changed(self, old, new):
        print "%s's name changed to %s." % (old, new) 
    
    @on_trait_change('hands.[l_fingers, hands.r_fingers]')
    def fingers_stauts(self, old, new):
        diff = new - old
        if diff < 0:
            print 'Poor %s just lost %s fingers.' % (self.name, abs(diff))
        else:
            print '%s just grew %s fingers.' % (self.name, abs(diff))

engine = sa.create_engine('postgres://postgres:postgres@localhost:5432/lynxgraphs', echo=False)
Session = sessionmaker(bind=engine)

HasTraitsORM.metadata.bind = engine
#HasTraitsORM.metadata.drop_all()
HasTraitsORM.metadata.create_all()
