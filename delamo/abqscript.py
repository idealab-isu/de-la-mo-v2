import sys
import os
import os.path
import ast

try:
    from io import StringIO
    pass
except ImportError:
    # Python 2.x
    from cStringIO import StringIO
    pass



class dummy(object):
    # dummy class used to identify this module and therefore this file
    pass


def build_abq_script(initinstrs,assemblyinstrs,bcinstrs,meshinstrs,runinstrs):
    # find installed directory
    thisfile=sys.modules[dummy.__module__].__file__
    thisdir=os.path.split(thisfile)[0]
    #pathdir=os.path.join(thisdir,"..")
    
    outbuf=StringIO(None)
    
    # load abqtemplate
    with open(os.path.join(thisdir,"abqtemplate.py"),"rb") as fh:
        line=fh.readline().decode('utf-8')

        while line != "":

            outbuf.write(line)
            
            #if line.startswith("###PATHLINE"):
            #    #outbuf.write(u"sys.path.append(r\"\"\"%s\"\"\")\n" % (pathdir))
            #    pass
            #elif line.startswith("###ABQFUNCS"):
            #    with file(os.path.join(thisdir,"abqfuncs.py")) as libs_fh:
            #        outbuf.write(libs_fh.read().decode('utf-8'))
            #        outbuf.write(u"\n")
            #        pass
            #    pass
            #elif line.startswith("###RUN_ABQPARAMS"):
                #for abqparams_filename in abqparams_filenames:
                #    with open(abqparams_filename,"rb") as abqparams_fh:
                #        abqparams_text=abqparams_fh.read().decode('utf-8')
                #        outbuf.write(abqparams_text)
                #        outbuf.write(u"\n")
                #        pass
                #    pass                
            #    pass
            if line.startswith("###RUN_INITINSTRS"):
                outbuf.write(initinstrs.get_code())
                outbuf.write(u"\n")
                pass
            elif line.startswith("###RUN_ASSEMBLYINSTRS"):
                outbuf.write(assemblyinstrs.get_code())
                outbuf.write(u"\n")
                pass
            elif line.startswith("###RUN_BCINSTRS"):
                outbuf.write(bcinstrs.get_code())
                outbuf.write(u"\n")
                pass
            elif line.startswith("###RUN_MESHINSTRS"):
                outbuf.write(meshinstrs.get_code())
                outbuf.write(u"\n")
                pass
            elif line.startswith("###RUN_RUNINSTRS"):
                outbuf.write(runinstrs.get_code())
                outbuf.write(u"\n")
                pass
            
            line=fh.readline().decode('utf-8')
            pass

        pass
    
    return outbuf.getvalue()
    
def write_abq_script(initinstrs,assemblyinstrs,bcinstrs,meshinstrs,runinstrs,outfilename):

    scriptstr=build_abq_script(initinstrs,assemblyinstrs,bcinstrs,meshinstrs,runinstrs)

    with open(outfilename,"w") as outfh:
        outfh.write(scriptstr)
        pass
    pass




def _capture_assignments_as_variables(syntax_tree_list,instrs,prevglobals,newglobals):
    # prevglobals are previously set globals
    # newglobals are updates from this set of assignments
    for syntax_element in syntax_tree_list:
        
        # A Function definition
        if isinstance(syntax_element, ast.FunctionDef):
            if syntax_element.name in prevglobals:
                raise ValueError("Attempt to override existing variable defining function %s" % (syntax_element.name))
                
            newglobals[syntax_element.name]=instrs.wrap_function(None,syntax_element.name)
            pass
                
        # An assignment to one or more targets
        if isinstance(syntax_element,ast.Assign) or isinstance(syntax_element,ast.AugAssign):
            if isinstance(syntax_element,ast.Assign):
                targets=syntax_element.targets
                pass
            else: 
                # AugAssign
                targets=[ syntax_element.target ]
                pass
                
            for target in syntax_element.targets:
                    
                if isinstance(target,ast.Name):
                    # simple assignments only (not assignments to attributes)
                    if target.id in prevglobals:
                        raise ValueError("Attempt to override existing variable when assigning %s" % (target.id))
                        
                    newglobals[target.id]=instrs.preexisting_variable(target.id)
                    pass
                pass
            pass
        if isinstance(syntax_element,ast.If) or isinstance(syntax_element,ast.For) or isinstance(syntax_element,ast.While):
            # Conditional or loop... Recursive call
            _capture_assignments_as_variables(syntax_element.body,instrs,prevglobals,newglobals)
            _capture_assignments_as_variables(syntax_element.orelse,instrs,prevglobals,newglobals)
            pass
        if hasattr(ast,"TryExcept") and isinstance(syntax_element,ast.TryExcept):
            # Try...Except block... Recursive call (python2)
            _capture_assignments_as_variables(syntax_element.body,instrs,prevglobals,newglobals)
            _capture_assignments_as_variables(syntax_element.handlers,instrs,globals)
            _capture_assignments_as_variables(syntax_element.orelse,instrs,prevglobals,newglobals)
            pass
        if hasattr(ast,"TryFinally") and isinstance(syntax_element,ast.TryFinally):
            # Try...Except block... Recursive call (python2)
            _capture_assignments_as_variables(syntax_element.body,instrs,prevglobals,newglobals)
            _capture_assignments_as_variables(syntax_element.finalbody,instrs,prevglobals,newglobals)
            pass
        if hasattr(ast,"try") and isinstance(syntax_element,ast.Try):
            # python3
            _capture_assignments_as_variables(syntax_element.body,instrs,prevglobals,newglobals)
            _capture_assignments_as_variables(syntax_element.handlers,instrs,globals)
            _capture_assignments_as_variables(syntax_element.orelse,instrs,prevglobals,newglobals)
            _capture_assignments_as_variables(syntax_element.finalbody,instrs,prevglobals,newglobals)
            pass

        if isinstance(syntax_element,ast.With):
            # With block... Recursive call
            _capture_assignments_as_variables(syntax_element.body,instrs,prevglobals,newglobals)
            pass

        if isinstance(syntax_element,ast.Import) or isinstance(syntax_element,ast.ImportFrom):
            # import or from xxx import statement
            for imported in syntax_element.names:
                # should be an ast.alias object
                asname=imported.asname
                if asname is None:
                    asname=imported.name
                    pass
                pass
                
                # No overwrite warning in this case because 
                # importing the same thing in both Abaqus
                # and our script is both common 
                # and useful
                #
                # Just don't assign in that case, so we
                # get the native copy where possible
                if asname not in prevglobals:
                    newglobals[asname]=instrs.preexisting_variable(asname)
                    pass
                pass

            pass
            

        pass
    pass
