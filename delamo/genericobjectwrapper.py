# Copyright 2016-2018 Iowa State University Research Foundation, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# This module defines a function called generate_wrapper
# that, given  dispatch function
# creates a new class that uses the dispatch function to wrap
# all attribute accesses (including method lookups) on the class. 
# It is otherwise as invisible
# as possible. 
#
import sys
import os
import types
import functools



# set of names of Python magic methods
# magicnames omits __new__, __init__, __getattribute__,  
# otherwise this list is based on http://www.rafekettler.com/magicmethods.html    
magicnames=set(["__del__", "__cmp__", "__eq__","__ne__","__lt__","__gt__","__le__", "__ge__", "__pos__", "__neg__", "__abs__", "__invert__", "__round__", "__floor__", "__ceil__", "__trunc__", "__add__", "__sub__", "__mul__", "__floordiv__", "__div__", "__truediv__", "__mod__", "__divmod__", "__pow__", "__lshift__", "__rshift__", "__and__", "__or__", "__xor__", "__radd__", "__rsub__", "__rmul__", "__rfloordiv__", "__rdiv__", "__rtruediv__", "__rmod__", "__rdivmod__", "__rpow__", "__rlshift__", "__rrshift__", "__rand__", "__ror__", "__rxor__", "__iadd__", "__isub__", "__imul__", "__ifloordiv__", "__idiv__", "__itruediv__", "__imod__", "__ipow__", "__ilshift__", "__irshift__", "__iand__", "__ior__", "__ixor__", "__int__", "__long__", "__float__", "__complex__", "__oct__", "__hex__", "__index__", "__trunc__", "__coerce__", "__str__", "__repr__", "__unicode__", "__format__", "__hash__", "__nonzero__", "__dir__", "__sizeof__","__delattr__","__setattr__","__len__","__getitem__", "__setitem__","__delitem__","__iter__","__reversed__", "__contains__", "__missing__","__call__", "__getattr__","__enter__","__exit__","__get__","__set__","__delete__","__copy__","__deepcopy__","__getinitargs__","__getnewargs__","__getstate__","__setstate__","__reduce__","__reduce_ex__"])

if sys.version_info >= (2,7):
    magicnames.add("__subclasscheck__")  # cannot assign __subclasscheck__ prior to python 2.6
    magicnames.add("__instancecheck__") # cannot assign __instancecheck__ prior to python 2.6
    pass

# Not all magic functions are wrappable... for example we shouldn't  wrap __eq__ because it should return a boolean
# likewise __coerce__
magicnames_wrappable=set(["__abs__", "__invert__", "__round__", "__floor__", "__ceil__", "__trunc__", "__add__", "__sub__", "__mul__", "__floordiv__", "__div__", "__truediv__", "__mod__", "__divmod__", "__pow__", "__lshift__", "__rshift__", "__and__", "__or__", "__xor__", "__radd__", "__rsub__", "__rmul__", "__rfloordiv__", "__rdiv__", "__rtruediv__", "__rmod__", "__rdivmod__", "__rpow__", "__rlshift__", "__rrshift__", "__rand__", "__ror__", "__rxor__", "__iadd__", "__isub__", "__imul__", "__ifloordiv__", "__idiv__", "__itruediv__", "__imod__", "__ipow__", "__ilshift__", "__irshift__", "__iand__", "__ior__", "__ixor__", "__int__", "__long__", "__float__", "__complex__", "__oct__", "__hex__", "__index__", "__trunc__", "__format__", "__delattr__","__setattr__","__getitem__", "__setitem__","__delitem__","__reversed__", "__missing__","__call__", "__getattr__","__enter__","__exit__","__get__","__set__","__delete__","__copy__","__deepcopy__"])


def generate_object_wrapper(dispatch_function=None,wrappertype=None,wrapperinstance=None):
    """Generate a wrapper class for "object" uses the dispatch_function
    to perform all method calls (including __getattribute__). Specifically dispatch_function
    is called as  

       dispatch_function(wrapped_object,methodname,*args,**kwargs)
    
    where: 
      methodname is the name of the method being called
      args, kwargs   parameters of the call

    NOTE: if methodname=="__init__" then you should NOT return an object
"""
    
    class wrappedobject(object):
        _wrappedobj_classname=None
        _wrappertype=wrappertype
        _wrapperinstance=wrapperinstance
        def __new__(cls,*args, **kwargs):
            import traceback
            #retval=dispatch_function(cls,"__new__",*args,**kwargs)
            #sys.stderr.write("wrappedobject %s\n%s\n\n" % (str(retval),"\n".join(traceback.format_stack())));
            #return retval;
            return dispatch_function(cls,"__new__",*args,**kwargs)
        
        def __init__(self,*args,**kwargs):
            return dispatch_function(self,"__init__",*args,**kwargs)

        def __getattribute__(self,attrname):
            return dispatch_function(self,"__getattribute__",attrname)

        def __str__(self):
            return "wrap 0x%lx" % (id(self))

        def __iter__(self):
            raise ValueError("Can not iterate over a wrapped object (control flow may not depend on wrapped object")

        def __subclasscheck__(self):
            raise ValueError("Can not check subclass status of a wrapped object")
            
        def __instancecheck__(self):
            raise ValueError("Can not check instance status of a wrapped object")
            
        pass

    # override magic methods if present in original. Magic methods need
    # to be explicitly added because they cannot be overridden
    # with __getattribute__() 

    
    #setattr(wrappedclass,"__str__",lambda self, *args, **kwargs: self._wrap(object.__getattribute__(class_to_wrap,"__str__"),args,kwargs))
    
    for magicname in magicnames_wrappable:
        attrfunc=lambda magicname: lambda self, *args, **kwargs: dispatch_function(self,magicname,*args,**kwargs)            
        setattr(wrappedobject,magicname,attrfunc(magicname))
        pass


    return wrappedobject






if __name__=="__main__":
    # Wrapper diagnostics

    def dispatch(wrappedobj,methodname,*args,**kwargs):
        if methodname == "__init__":
            newobj=None
            pass
        elif methodname=="__new__":
            newobj=object.__new__(wrappedobj)  # wrappedobj is actually the class object 
            pass
        else:
            newobj=objwrap()
            pass

        print("dispatch %s.%s (%s,%s) -> %s" % (str(wrappedobj),methodname,",".join([str(arg) for arg in args]),",".join(["%s=%s" % (argname,str(kwargs[argname])) for argname in kwargs.keys()]),str(newobj)))
        return newobj

    objwrap = generate_object_wrapper(dispatch)

    myobj=objwrap(7)

    plus5 = myobj + 5
    
    plus5.myfun(32)
    
