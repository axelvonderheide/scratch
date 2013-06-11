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

##########################################################################
# has_traits_orm.py
##########################################################################

# Enthought library imports
from enthought.preferences.api import Preferences
from enthought.traits.api import \
    HasTraits, MetaHasTraits, Instance, Type, Int, Str, Bool, Float, Any, Class, \
    List, ListBool, ListFloat, ListInt, ListStr, String, Enum, DictStrAny

# Package imports
import sqlalchemy
from sqlalchemy import *
from sqlalchemy.orm import *
from sqlalchemy.exceptions import InvalidRequestError
from sqlalchemy.orm.interfaces import MapperProperty
from sqlalchemy.ext.declarative import *
from sqlalchemy.ext.declarative import _as_declarative

TRAIT_MAPPING = {
         Int : 'sqlalchemy.Integer',
         Str : 'sqlalchemy.Text',
         Enum : 'sqlalchemy.Text',
         String : 'sqlalchemy.Text',
         Float : 'sqlalchemy.Float',
         Bool : 'sqlalchemy.Boolean',
         }

class DeclarativeMetaTraits(MetaHasTraits):

    def __init__(cls, classname, bases, dict_):
        if '_decl_class_registry' in cls.__dict__:
            return MetaHasTraits.__init__(cls, classname, bases, dict_)

       # create sql columns from flagged traits
        if '__class_traits__' in cls.__dict__:
            traits = cls.__dict__['__class_traits__']
            for key, trait in traits.items():
                if getattr( trait, 'sqldb' ):
                    args = getattr( trait, 'col_args' ) or ()
                    kwargs = getattr( trait, 'col_kwargs' ) or {}
                    if 'name' not in kwargs:
                        kwargs['name'] = key
                    if 'type_' not in kwargs:
                        kwargs['type_'] = eval(TRAIT_MAPPING[type(trait.trait_type)])

                    c = Column(*args, **kwargs)
                    dict_[key] = c

        _as_declarative(cls, classname, dict_)
        return MetaHasTraits.__init__(cls, classname, bases, dict_)


    def __setattr__(cls, key, value):

        if '__mapper__' in cls.__dict__:
            if isinstance(value, Column):
                _undefer_column_name(key, value)
                cls.__table__.append_column(value)
                cls.__mapper__.add_property(key, value)
            elif isinstance(value, MapperProperty):
                cls.__mapper__.add_property(key, _deferred_relation(cls, value))
            else:
                MetaHasTraits.__setattr__(cls, key, value)
        else:
            MetaHasTraits.__setattr__(cls, key, value)



HasTraitsORM = declarative_base(name='HasTraitsORM', cls=HasTraits,
metaclass=DeclarativeMetaTraits, constructor=None)



###################################################


That works fine for simple classes and you can add normal sql Columns
in you class too without them being Traits.  The problem is though
that realations dont work. An quick example:

###############################################################################
# Some module
###############################################################################

from has_traits_orm import HasTraitsORM
from sqlalchemy import create_engine, Column, ForeignKey
from sqlalchemy.orm.collections import InstrumentedList
from sqlalchemy.orm import sessionmaker, relation, backref



class User(HasTraitsORM):
    __tablename__ = 'users'

    id = Column(Integer, primary_key=True)
    name = Str(sqldb=True)

    def __repr__(self):
        return "<User('%s')>" % (self.name)

class Address(HasTraitsORM):
    __tablename__ = 'addresses'

    id = Column(Integer, primary_key=True)
    email_address = Str(sqldb=True)
    user_id = Int(sqldb=True, col_args=(ForeignKey('users.id'),) )
    user = relation(User, backref=backref('addresses'),
collection_class=InstrumentedList)

    def __repr__(self):
        return "<Address('%s')>" \
            % (self.email_address)

engine = create_engine('sqlite:///:memory:', echo=True)

# bind the metadata
HasTraitsORM.metadata.bind = service.engine
HasTraitsORM.metadata.create_all()

Session = sessionmaker(bind=engine)
session = Session()

A1 = Address(email_address='john at gmail.com')
A3 = Address(email_address='jn at gmail.com')

U1 = User(name='John')
U1.addresses = [A1, A3]   #   <-----     bombs out right here

session.add(U1)
session.commit()




The bug lies with the collections classes.  After some serious
debugging I realised that when sqlalchemy sets up the mapper it must
somehow override object assigment, so normally (without Traits), when
assigning the addresses to U1.addresses sqlalchemy does some magic to
instrument the list passed in for its own puposes.  But as soon as we
inherit from HasTraits instead of sqlalchemy doing its magic execution
jumps to has_traits.py to the following method:

has_traits.py  line 3390
    #---------------------------------------------------------------------------
    #  Returns the trait definition for a specified name when there is no
    #  explicit definition in the class:
    #---------------------------------------------------------------------------

    def __prefix_trait__ ( self, name, is_set ):
