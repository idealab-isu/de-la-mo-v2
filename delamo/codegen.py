import traceback
import inspect
import sys
import numpy as np

if sys.hexversion >= 0x03000000: # python 3.x
    # python 3 doesn't have long or unicode
    long=int
    unicode=str
    pass 

genericobjectwrapper=None

try:
    from . import genericobjectwrapper 
    pass
except ImportError:
    # From inside Abaqus we're unlikely to have genericobjectwrapper available
    # but we don't need it 
    pass

# Use for wrappers that need to be globals,
# but we don't want to monkeypatch globals everywhere,
# so we can create static references to namebindingwrappers
# in the module definition, then
# assign it to a wrapped object
# in codevariables.namedbindings
class namedbinding_wrapper(object):
    codestore=None
    wrapped=None
    name=None
    localvars=None # dictionary by id(codevariables object) of variables evaluated in that context
    
    def __init__(self,name):
        # reuse codestore here because it's pretty much the same as what we need. 
        self.name=name
        self.codestore=codestore(wrappertype=self.__class__,wrapperinstance=self)
        self.wrapped=self.codestore.preexisting_variable(name)
        self.localvars={}
        pass
    
    def evaluate(self,codevariables_id,wrapper_replacement,value):
        if codevariables_id not in self.localvars:
            self.localvars[codevariables_id]={}
            pass
        self.localvars[codevariables_id][self.name]=wrapper_replacement
        exec(self.codestore.get_code(),self.localvars[codevariables_id],self.localvars[codevariables_id])
        # now that we've executed the code so far once,
        # empty it so we don't re-execute
        self.codestore.code_lines=[]
        
        retval=eval(self.codestore.serialize_object(value),self.localvars[codevariables_id],self.localvars[codevariables_id])
        return retval
    pass



class codevariables(object):
    variabledict=None # dictionary by id of (wrapped object, assignment codestore,assignment_has_occurred(i.e. True or False))
    preexisting_vars=None # dictionary by id of preexisting variable name
    instancevariabledict=None # dictionary by id of (containing object, member name)
    classdict=None # dictionary by id of wrapped class object
    
    globals=None  # if not None, identify variable names from this global variable context

    namedmapping=None # dictionary by id of (human-interpretable variable name, wrapped object). Names are  _cg_<name>, _cg_<name>1, etc.... used as substitute for variabledict
    preexisting_names=None # Set of preexisting names

    namedbindings= None # dictionary by namedbinding_wrapper-id of (namedbinding_wrapper,wrapped object)

    def __init__(self,globals=None):
        self.variabledict={}
        self.preexisting_vars={}
        self.instancevariabledict={}
        self.classdict={}
        self.globals=globals
        self.namedmapping={}
        self.preexisting_names=set([])
        
        self.namedbindings={}
        

        pass


    def newnamedmapping(self,basename,wrappedobj):
        
        varid=id(wrappedobj)


        newname=basename
        cnt=1
        while newname in self.preexisting_names:
            newname="%s%d" % (basename,cnt)
            cnt+=1
            pass
            
        # ... and post a new namedmapping that's more meaningful
        self.namedmapping[varid]=(newname,wrappedobj)
        self.preexisting_names.add(newname)
        return newname

    def update_namedmapping_from_dict(self,vardict):
        # Find keys (ids) in variabledict that are also ids of 
        # named variables in vardict. 
        # Then use the name to give a new name for the 
        # variabledict entry


        vardict_names_by_id={ id(vardict[varname]): varname for varname in vardict if id(vardict[varname]) in self.variabledict }
        
        #preexisting_names=frozenset([ self.namedmapping[varid][0] for varid in self.namedmapping ])
        
        # Iterate over all global variables that have values in our variabledict
        for varid in vardict_names_by_id:
            (wrappedobj,codestoreobj,assignment_written)=self.variabledict[varid]
            if not assignment_written:
                continue # Cannot rename variables that haven't actually been written yet
            # remove from variabledict 
            del self.variabledict[varid]

            #if vardict_names_by_id[varid]=="newobj":
            #    import pdb
            #    pdb.set_trace()

            newname=self.newnamedmapping("_cg_%s" % (vardict_names_by_id[varid]),wrappedobj)


            # Generate named mapping in the same codestore 
            # where the wrapper was generated
            codestoreobj.append_code("%s=_codegen_%d" % (newname,varid),tbstriplayers=0,show_traceback=True)
            

            pass
        pass


    def update_namedmapping(self,ignorestackframes=1):
        # Find keys (ids) in variabledict that are also ids of globals
        # or local variables in the current stack backtrace
        # Then use the global name to give a new name for the 
        # variabledict entry

        # ... disabled if self.globals was not set on initialization
        if self.globals is None:
            return

        self.update_namedmapping_from_dict(self.globals)

        stacktrace=inspect.stack()
        for framecnt in range(len(stacktrace)-1,1+ignorestackframes-1,-1):
            self.update_namedmapping_from_dict(stacktrace[framecnt][0].f_locals)
            pass
        pass


    def add_named_binding(self,namedwrapper,realwrapper):
        self.namedbindings[id(namedwrapper)]=(namedwrapper,realwrapper)
        pass
    pass


class codestore(object):
    code_lines=None # list of lines
    codevariables=None # class codevariables (above)
    
    #wrapperclass_cache=None # dictionary by id of unwrapped class of wrapped class
    wrapped_object_class=None
    
    suppress_cnt=None
    

    
    def __init__(self,codevariableinstance=None,globals=None,wrappertype=None,wrapperinstance=None):
        self.code_lines=[]
        if codevariableinstance is None:
            codevariableinstance=codevariables(globals=globals)
            pass
        self.codevariables=codevariableinstance
        #self.wrapperclass_cache={}
        self.suppress_cnt=0
        self.wrapped_object_class=genericobjectwrapper.generate_object_wrapper(dispatch_function=self.method_wrapper,wrappertype=wrappertype,wrapperinstance=wrapperinstance)
        

        pass

    def suppress_code(self):  # suppress code generation 
        self.suppress_cnt+=1
        pass
    def permit_code(self):
        self.suppress_cnt-=1
        pass
    
    def get_code(self):
        return u"\n".join(self.code_lines)

    def save_code(self,filename):
        fh=file(filename,"w")
        fh.write(self.get_code())
        fh.close()
        pass
    

    def append_code(self,codestr,tbstriplayers=2,show_traceback=True):
        """tbstriplayers is number of most recent entries
        to strip from the traceback that is included with the code.
        Includes the call to append_code() itself"""

        if show_traceback:
            tbobj=traceback.extract_stack()[:-tbstriplayers]

            tracebackstr="".join(traceback.format_list(tbobj))

            commentedtracebackstr="#" + tracebackstr.replace('\n','\n#')+'\n'
            pass
        else:
            commentedtracebackstr=""
            pass
        self.code_lines.append(commentedtracebackstr+codestr+"\n\n")
        pass
    def assign_variable(self,name,value):
        if self.suppress_cnt==0:
            self.codevariables.update_namedmapping() # Capture newly defined globals
            pass

        self.append_code("%s=%s" % (name,self.serialize_object(value)))

        
        if id(value) in self.codevariables.variabledict:
            # simple rename of an existing variable
            self.codevariables.preexisting_vars[id(value)]=name
            del self.codevariables.variabledict[id(value)]
            obj=value
            pass
        elif id(value) in self.codevariables.namedmapping:
            # simple rename of an existing variable
            self.codevariables.preexisting_vars[id(value)]=name
            del self.codevariables.namedmapping[id(value)]
            obj=value
            pass
            
        else:
            
            obj=object()
            self.codevariables.preexisting_vars[id(obj)]=name
            pass
        
        return obj
        
    def add_named_binding(self,namedwrapper,realwrapper):
        self.codevariables.add_named_binding(namedwrapper,realwrapper)
        pass
    
    def preexisting_variable(self,name):
        self.suppress_code()
        obj=self.wrapped_object_class()
        self.permit_code()
        
        #self.codevariables.variabledict[id(obj)]=obj
        self.codevariables.preexisting_vars[id(obj)]=name
        
        return obj
    
    def wrap_function(self,function,function_name=None):
        if function_name is None:
            function_name=function.__name__
            pass
        return self.preexisting_variable(function_name)

    #def wrap_function(self,function,function_name=None):
    #
    #    if function_name is None:
    #        function_name=function.__name__
    #        pass
    #
    #    def function_wrapper(*args,**kwargs):
    #
    #        self.suppress_code() # prevent constructer from generating code, as we put the code in below
    #        retobj=self.wrapped_object_class()  
    #        self.permit_code()
    #        
    #        # store return variable as variable
    #        self.codevariables.variabledict[id(retobj)]=retobj
    #        
    #        # print("Writing code: Function name=%s" % (function.__name__))
    #
    #        self.append_code('%s=%s(*%s,**%s)' % (self.serialize_object(retobj),function_name,self.serialize_list(args),self.serialize_dict(kwargs)))
    #    
    #    
    #        #retval=origmethod(*args,**kwargs)
    #        
    #        #print("method_wrapper returning %s" % (str(retobj)))
    #        
    #        return retobj
    #    return function_wrapper

    def wrap_class(self,unwrappedclass):
        return self.preexisting_variable(unwrappedclass.__name__)
        

    #def wrap_class(self,unwrappedclass):
        # ***!!! TODO: unwrappedclass parameter no longer used or needed 
        #print("wrapping class %s" % (str(unwrappedclass)))
        # if id(unwrappedclass) in self.wrapperclass_cache:
        #     return self.wrapperclass_cache[id(unwrappedclass)]
        
        #wrappedclass=genericobjectwrapper.generate_object_wrapper(unwrappedclass,dispatch_function=self.method_wrapper)

        #self.wrapperclass_cache[id(unwrappedclass)]=wrappedclass
        #self.codevariables.classdict[id(wrappedclass)]=wrappedclass
        # print("wrapped class=%s" % (str(wrappedclass)))
    #    return self.wrapped_object_class
    
    def rewrapobj(self,prewrappedobj):
        """Rewrap an object wrapped by another codegen instance
        for use/execution by this codegen instance"""

        # Create a new wrapped object of this class
        self.suppress_code()  # suppress code generation 
        wrappedinstance=self.wrapped_object_class()
        self.permit_code()
        
        # store return variable as variable
        self.codevariables.variabledict[id(wrappedinstance)]=(wrappedinstance,self,False)

        # Then write assignment code so the two variables are equivalent
        self.append_code('%s=%s' % (self.serialize_object(wrappedinstance,unassigned_ok=True),self.serialize_object(prewrappedobj)))

        # Now mark it as written
        self.codevariables.variabledict[id(wrappedinstance)]=(wrappedinstance,self,True)

        return wrappedinstance
    
    #def accumulate_instanceclass_dict(self,dct,obj,cls):
    #    #if hasattr(cls,"__bases__"):
    #    for basecls in cls.__bases__:
    #        self.accumulate_instanceclass_dict(dct,obj,basecls)
    #        pass
    #    if hasattr(cls,"_%s__codegen_instanceclass_dict" % (cls.__name__)):
    #        dct.update(getattr(cls,"_%s__codegen_instanceclass_dict" % (cls.__name__)))
    #        pass
    #    pass

    #def accumulate_methodreturnclass_dict(self,dct,obj,cls):
    #    #if hasattr(cls,"__bases__"):
    #    for basecls in cls.__bases__:
    #        self.accumulate_methodreturnclass_dict(dct,obj,basecls)
    #        pass
    #    
    #    
    #    if hasattr(cls,"_%s__codegen_methodreturnclass_dict" % (cls.__name__)):
    #        dct.update(getattr(cls,"_%s__codegen_methodreturnclass_dict" % (cls.__name__)))
    #        pass
    #    pass

    
    def method_wrapper(self,wrappedobj,methodname,*args,**kwargs):
        #print("method_wrapper: methodname=%s" % (methodname))
        if self.suppress_cnt == 0:
            self.codevariables.update_namedmapping() # give any global variables better names
            pass

        if methodname=="__new__":
 
            newobj=object.__new__(wrappedobj)  # wrappedobj is the class object in this case
            
            if self.suppress_cnt == 0:
                self.codevariables.variabledict[id(newobj)]=(newobj,self,False)

                # !!!*** This next line is probably broken.... should be "%s=%s(*%s,**%s)" because serialize_object calls add the prefix
                # !!!*** so this probably hasn't been tested!!!***
                self.append_code('_codegen_%s=_codegen_%s(*%s,**%s)' % (self.serialize_object(newobj,unassigned_ok=True),self.serialize_object(wrappedobj),self.serialize_list(args),self.serialize_dict(kwargs)),tbstriplayers=2)

                # Mark variable as written
                self.codevariables.variabledict[id(newobj)]=(newobj,self,True)

                pass
            
            #import pdb
            #pdb.set_trace()
            return newobj
        elif methodname=="__init__":
            # generation handled by __new__ (above) 
            return None # __init__ always returns none
        elif methodname=="__getattribute__":
            self.suppress_code()  # suppress code generation 
            newobj=self.wrapped_object_class()
            self.permit_code()
            
            # store return variable as variable
            #self.codevariables.variabledict[id(newobj)]=(newobj,self)
            wrapped_serialization=self.serialize_object(wrappedobj)
            if wrapped_serialization.startswith("_cg_"):
                wrapped_serialization=wrapped_serialization[4:]
                pass

            # Give the new variable an intuitive name
            self.codevariables.newnamedmapping("_cg_%s_dot_%s" % (wrapped_serialization,args[0]),newobj)
            
            # print("Writing code: Methodname=%s" % (methodname))
            self.append_code('%s=%s.%s' % (self.serialize_object(newobj,unassigned_ok=True),self.serialize_object(wrappedobj),args[0]),tbstriplayers=3)
            return newobj

        elif methodname=="__call__":
            self.suppress_code()  # suppress code generation 
            newobj=self.wrapped_object_class()
            self.permit_code()
            
            # store return variable as variable
            self.codevariables.variabledict[id(newobj)]=(newobj,self,False)
            
            # print("Writing code: Methodname=%s" % (methodname))
            #self.append_code('%s=%s(*%s,**%s)' % (self.serialize_object(newobj),self.serialize_object(wrappedobj),self.serialize_list(args),self.serialize_dict(kwargs)),tbstriplayers=3)
            self.append_code('%s=%s%s' % (self.serialize_object(newobj,unassigned_ok=True),self.serialize_object(wrappedobj),self.serialize_params(args,kwargs)),tbstriplayers=3)

            # Mark variable as written
            self.codevariables.variabledict[id(newobj)]=(newobj,self,True)
            return newobj
        elif methodname=="__getitem__" and len(args)==1 and len(kwargs)==0:
            self.suppress_code()  # suppress code generation 
            newobj=self.wrapped_object_class()
            self.permit_code()
            
            # store return variable as variable
            self.codevariables.variabledict[id(newobj)]=(newobj,self,False)
            
            self.append_code('%s=%s[%s]' % (self.serialize_object(newobj,unassigned_ok=True),self.serialize_object(wrappedobj),self.serialize_object(args[0])),tbstriplayers=3)

            # Mark variable as written
            self.codevariables.variabledict[id(newobj)]=(newobj,self,True)
            
            return newobj

        assert(id(wrappedobj) in self.codevariables.variabledict or id(wrappedobj) in self.codevariables.preexisting_vars or id(wrappedobj) in self.codevariables.instancevariabledict or id(wrappedobj) in self.codevariables.classdict or id(wrappedobj) in self.codevariables.namedmapping) # can only operate on known objects

        # Generate return object
        
        self.suppress_code()  # suppress code generation 
        newobj=self.wrapped_object_class()
        self.permit_code()

        # print(methodname)
        # print(methodreturnclass_dict)
        
        # store return variable as variable
        self.codevariables.variabledict[id(newobj)]=(newobj,self,False)

        # print("Writing code: Methodname=%s" % (methodname))
        #self.append_code('%s=%s.%s(*%s,**%s)' % (self.serialize_object(newobj),self.serialize_object(wrappedobj),methodname,self.serialize_list(args),self.serialize_dict(kwargs)),tbstriplayers=3)
        self.append_code('%s=%s.%s%s' % (self.serialize_object(newobj,unassigned_ok=True),self.serialize_object(wrappedobj),methodname,self.serialize_params(args,kwargs)),tbstriplayers=3)

        # Mark variable as written
        self.codevariables.variabledict[id(newobj)]=(newobj,self,True)
        
        #retval=origmethod(*args,**kwargs)

        #print("method_wrapper returning %s" % (str(retobj)))
        
        return newobj

    def serialize_params(self,args,kwargs):
        retval='('
        for argcnt in range(len(args)):
            retval+=self.serialize_object(args[argcnt])
            if argcnt < len(args)-1 or len(kwargs) > 0:
                retval+=","
                pass
            pass
            
        keys=list(kwargs.keys())
        for keycnt in range(len(keys)):
            retval+="%s=%s" % (keys[keycnt],self.serialize_object(kwargs[keys[keycnt]]))
            if keycnt < len(keys)-1:
                retval+=","
                pass
            pass
        retval+=")"
        return retval

    def serialize_list(self,lst):
        retval='['
        for entry in lst:
            retval+=self.serialize_object(entry)+','
            pass
        retval+=']'
        return retval
    
    def serialize_tuple(self,tup):
        retval='('
        for entry in tup:
            retval+=self.serialize_object(entry)+','
            pass
        retval+=')'
        return retval

    def serialize_dict(self,dct):
        retval='{'
        for key in dct.keys():
            retval+=self.serialize_object(key)+':'+self.serialize_object(dct[key])+','
            pass
        retval+='}'
        return retval

    def serialize_set(self,setobj):
        retval='set(['
        for key in setobj:
            retval+=self.serialize_object(key)+','
            pass
        retval+='])'
        return retval


    def serialize_object(self,obj,unassigned_ok=False):
        # print("serialize_object: %s, __class__=%s" % (id(obj),obj.__class__.__name__))

        # is this a namedbinding_wrapper? if so unwrap and evaluate it
        # and serialize the evaluated wrapper object
        wt=None
        try: 
            wt=object.__getattribute__(obj,"_wrappertype")            
            pass
        except AttributeError:
            pass
        if wt is namedbinding_wrapper and len(list(self.codevariables.namedbindings.keys())) > 0:
            nbw_instance = object.__getattribute__(obj,"_wrapperinstance")
            
            obj=nbw_instance.evaluate(id(self.codevariables),self.codevariables.namedbindings[id(nbw_instance)][1],obj)
            pass

        
        if id(obj) in self.codevariables.namedmapping:
            return self.codevariables.namedmapping[id(obj)][0]
        if id(obj) in self.codevariables.variabledict:
            assert(unassigned_ok or self.codevariables.variabledict[id(obj)][2]) # variable must have been written unless unassigned_ok has been set
            return "_codegen_%d" % (id(obj))
        if id(obj) in self.codevariables.instancevariabledict:
            return "%s.%s" % (self.serialize_object(self.codevariables.instancevariabledict[id(obj)][0]),self.codevariables.instancevariabledict[id(obj)][1])
        # is this one of our preexisting variables?
        if id(obj) in self.codevariables.preexisting_vars:
            return self.codevariables.preexisting_vars[id(obj)] # return  variable name        

        if isinstance(obj,list):
            return self.serialize_list(obj)
        if isinstance(obj,set):
            return self.serialize_set(obj)
        if isinstance(obj,tuple):
            return self.serialize_tuple(obj)
        if isinstance(obj,dict):
            return self.serialize_dict(obj)
        if isinstance(obj,int) or isinstance(obj,long) or isinstance(obj,float):
            return repr(obj)
        if isinstance(obj,str) or isinstance(obj,unicode):
            return repr(obj)
        if isinstance(obj,np.ndarray):
            return "np.array(%s,dtype=np.%s)" % (np.array2string(obj,separator=',',suppress_small=False,threshold=np.inf,floatmode='unique'),repr(obj.dtype))

        
        ## is this one of our wrapped classes...
        #if hasattr(obj,"_wrappedobj"): # and obj._wrappedobj is None and id(obj._wrap_class_to_wrap) in self.wrapperclass_cache:
        #    #return self.wrapperclass_cache[id(obj._wrap_class_to_wrap)]._wrap_class_to_wrap.__name__
        #    return obj._wrap_class_to_wrap.__name__

        if obj is None:
            return "None"
        
        
        raise ValueError("Cannot serialize parameter %s" % (str(obj)))
    
    pass


# Decorator for classes that will have code generation
def codegen_class(classobj):
    # ***!!! This should no longer be necessary, so it is a no-op

    #classname=classobj.__name__
    #
    #mrcdict=getattr(classobj,"_%s__codegen_methodreturnclass_dict" % (classname))
    ## Fixup 'None' values in dictionary as referring to this class
    #
    #for key in mrcdict:
    #    if mrcdict[key] is None:
    #        mrcdict[key]=classobj
    #        pass
    #    pass
    return classobj







