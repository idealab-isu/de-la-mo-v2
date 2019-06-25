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


import subprocess
import collections
import ast
import sys
import os
import os.path
import copy
import glob
import re
from collections import OrderedDict

import redbaron
from redbaron import RedBaron

from . import simple


try:
    basestring
    pass
except NameError:
    # python3
    basestring = str
    pass

# These represent the required and keyword parameters to
# bond_layers() from api.py
bond_layers_params = [ "DM","layer1","layer2" ]
if redbaron is not None:
    bond_layers_default_params = OrderedDict([
        ("defaultBC", redbaron.RedBaron("\"TIE\"")[0]),
        ("delamBC", redbaron.RedBaron("\"CONTACT\"")[0]),
        ("delamRingBC", redbaron.RedBaron("\"NONE\"")[0]),
        ("CohesiveInteraction", redbaron.NameNode.from_fst({"type":"name","value":"None"})),
        ("ContactInteraction", redbaron.NameNode.from_fst({"type":"name","value":"False"})),
        ("delaminationlist", redbaron.NameNode.from_fst({"type":"name","value":"None"})),
        ("master_layer", redbaron.NameNode.from_fst({"type":"name","value":"None"})),
        ("cohesive_layer", redbaron.NameNode.from_fst({"type":"name","value":"None"})),
        ("delamo_sourceline", redbaron.NameNode.from_fst({"type":"name","value":"None"})),
        ("delamo_phase", redbaron.NameNode.from_fst({"type":"name","value":"None"})),
        ("delamo_basename", redbaron.NameNode.from_fst({"type":"name","value":"None"})),
    ])
else:
    output_filenames_default_params = None
    pass

# NOTE: Ideally in the RedBaron trees we would
# track the provenance of each line so when we
# reconstruct a source file and parse it into
# python we can repoint the line # references at
# the original source. 

def hasnamedparameter(callnode,name):
    commaproxylist=callnode.value
    for callargumentnode in commaproxylist:
        if not isinstance(callargumentnode,redbaron.CallArgumentNode):
            continue
        if callargumentnode.target is None:
            # unnamed parameter
            continue
        if callargumentnode.target.value==name:
            return True
        pass
    return False



def getparameter(callnode,name,index):
    commaproxylist=callnode.value
    listindex=0
    for callargumentnode in commaproxylist:
        assert(isinstance(callargumentnode,redbaron.CallArgumentNode))
        
        if callargumentnode.target is None:
            # unnamed parameter
            if (listindex==index):
                return callargumentnode.value
            pass
        elif callargumentnode.target.value==name:
            return callargumentnode.value
        index+=1
        pass
    return False

def troubleshoot_bond_layers(rb_tree):
    # bond_layers() calls not allowed inside functions or comprehensions or while loops...
    # This function checks for that

    
    # atomtrailers is a sequence of foo.bar.fubar(arg1,arg2)
    # where foo, bar and fubar would be a NameNode and (arg1,arg2) is a CallNode
    bondlayers_calls=rb_tree.find_all("atomtrailers",
                                      lambda node: (node[0].value=="bond_layers" and
                                                    isinstance(node[1],redbaron.CallNode)))
    
    # bondlayers_calls[:][0] is AtomTrailersNode; [1] is CallNode; [1][0] is first parameter, etc. 
    for bondlayers_call in bondlayers_calls:
        thisnode=bondlayers_call
        while thisnode.parent is not None:
            thisnode=thisnode.parent
            if isinstance(thisnode,redbaron.DefNode):
                raise ValueError("bond_layers() call inside function, line %d" % (bondlayers_call.absolute_bounding_box.top_left.line))
            if isinstance(thisnode,redbaron.ComprehensionLoopNode):
                raise ValueError("bond_layers() call inside comprehension, line %d" % (bondlayers_call.absolute_bounding_box.top_left.line))
            if isinstance(thisnode,redbaron.WhileNode):
                raise ValueError("bond_layers() call inside while loop, line %d" % (bondlayers_call.absolute_bounding_box.top_left.line))
            
            pass
        pass
    
    
    pass



def unwrap_loops(rb_tree):

    # Find all loops that contain bond_layers calls
    while True:        
        loop=rb_tree.find("for",
                          lambda node: (node.find("atomtrailers",
                                                  lambda atnode: (atnode[0].value=="bond_layers" and
                                                                  isinstance(atnode[1],redbaron.CallNode))) is not None))
        
        if loop is None:
            # No more loops
            break

        #rangeparams=None
        tupleentries=None
        # Check that loop target is a tuple or range(...)
        if (isinstance(loop.iterator,redbaron.NameNode) and
            isinstance(loop.target.value[0],redbaron.NameNode) and
            loop.target.value[0].value == "range" and
            isinstance(loop.target.value[1],redbaron.CallNode) and
            isinstance(loop.target.value[1].value,redbaron.CommaProxyList) and
            len(loop.target.value[1].value) <= 3 and
            all([ isinstance(rangeparam.value,redbaron.IntNode) for rangeparam in loop.target.value[1].value ])):
            # we have a range(...)
            # Extract list of range parameters. 
            rangeparams=[ int(rangeparam.value.value) for rangeparam in loop.target.value[1].value ]
            # Expand range into a tuple
            tupleentries = [ redbaron.IntNode({"section": "number","type": "int","value": ("%d" % (rangevalue)) }) for rangevalue in range(*rangeparams) ]
            pass
        elif (isinstance(loop.iterator,redbaron.NameNode) and
              (isinstance(loop.target,redbaron.TupleNode) or
               isinstance(loop.target,redbaron.ListNode))):
            # We are iterating over a tuple or list
            tupleentries=[ entry for entry in loop.target.value ]
            pass
        
        if tupleentries is None:
            raise ValueError("For loop containing bond_layers does not use simple list/tuple or integer range() as loop sequence, line %d" % (loop.absolute_bounding_box.top_left.line))

        iterator_var=loop.iterator.value  # variable used as iterator

        # check for continue or break or for..else statements
        # NOTE: This is more constraining than necessary
        # as a break or continue within a nested loop would be OK
        if loop.value.find("breaknode") is not None or loop.value.find("continuenode") is not None or getattr(loop,"else") is not None:
            raise ValueError("For loop containing bond_layers contains break or continue statement or for..else clause, line %d" % (loop.absolute_bounding_box.top_left.lie))
        
        
        
        
        # Unwrap loop
        template=loop.value.copy() #copy.deepcopy(loop.value)

        # de-indent template
        #deindent=len(template[0].indentation)-len(loop.indentation)
        # Indentation after copy() is incorrect (Redbaron bug)
        # so pull the indentation from the original NodeProxyList
        deindent=len(loop.value[0].indentation)-len(loop.indentation)
        template.decrease_indentation(deindent)
        
        loop_parent=loop.parent
        loop_parent_index=loop.index_on_parent

        # delete original
        del loop_parent[loop_parent_index]
        
        # iterate loop
        toinsert=[]
        
        for tupleentry in tupleentries:

            # iterator assignment
            toinsert.append(redbaron.AssignmentNode({
                "first_formatting": [],
                "operator": "",
                "second_formatting": "",
                "target": {"type":"name","value":iterator_var},
                "type": "assignment",
                "value": tupleentry.fst()}))
            
            # loop body
            for elem in template:
                toinsert.append(elem.copy())
                pass
            pass
        #for o in toinsert:
        #    loop_parent.insert(loop_parent_index,o);
        #    loop_parent_index+=1
        #    import pdb
        #    pdb.set_trace()
        loop_parent.insert_multiple(loop_parent_index,toinsert)
        #loop_parent.extend(toinsert)
        
        pass
    pass


def annotate_bond_layers_calls(basename,basename_withpath,output_directory,phase,rb_tree):
    # bond_layers() calls not allowed inside functions or comprehensions or while loops...
    # This function checks for that

        
    
    # atomtrailers is a sequence of foo.bar.fubar(arg1,arg2)
    # where foo, bar and fubar would be a NameNode and (arg1,arg2) is a CallNode
    bondlayers_calls=rb_tree.find_all("atomtrailers",lambda node: node[0].value=="bond_layers" and isinstance(node[1],redbaron.CallNode))

    for bondlayers_call in bondlayers_calls:
        callnode=bondlayers_call[1]

        insertpos=last_kwarg_index(callnode)

        linearg = callnode.value.find("callargument",lambda node: node.target is not None and node.target.value=="delamo_sourceline")
        if linearg is not None:
            # already have a kwarg for delamo_sourceline
            linearg.value.value="%d" % (callnode.absolute_bounding_box.top_left.line)
            pass
        else:
            # add a new kwarg for delamo_sourceline
            callnode.value.insert(insertpos,
                                  redbaron.CallArgumentNode({"type": "call_argument",
                                                             "first_formatting": [],
                                                             "second_formatting":[],
                                                             "target": {
                                                                 "type": "name",
                                                                 "value":"delamo_sourceline"},
                                                             "value": {
                                                                 "first_formatting": [],
                                                                 "second_formatting": [],
                                                                 "type": "int",
                                                                 "value": "%d" % (callnode.absolute_bounding_box.top_left.line)}}))
            
            # add a new kwarg for delamo_basename
            delamo_basename_arg=redbaron.CallArgumentNode({"type": "call_argument",
                                                           "first_formatting": [],
                                                           "second_formatting":[],
                                                           "target": {
                                                               "type": "name",
                                                               "value":"delamo_basename"},
                                                           "value": {
                                                               "first_formatting": [],
                                                               "second_formatting": [],
                                                               "type": "string",
                                                               "value": "\"%s\"" % (basename.encode('unicode_escape'))}})
        
            #delamo_basename_arg.value=copy.deepcopy(basename_value)
            callnode.value.insert(insertpos,delamo_basename_arg)

            # add a new kwarg for delamo_phase
            delamo_phase_arg=redbaron.CallArgumentNode({"type": "call_argument",
                                                        "first_formatting": [],
                                                        "second_formatting":[],
                                                        "target": {
                                                            "type": "name",
                                                            "value":"delamo_phase"},
                                                        "value": {
                                                            "first_formatting": [],
                                                            "second_formatting": [],
                                                            "type": "string",
                                                            "value": "\"%s\"" % (phase.encode('unicode_escape')) }})
            #delamo_phase_arg.value=copy.deepcopy(phase_value)
            callnode.value.insert(insertpos,delamo_phase_arg)
            pass
        
        pass

    pass



def linenumber_and_index_from_filename(delamfilename):
    matchobj=re.match(r"""layerdelam_([0-9]+)(_[0-9]+)?""",delamfilename,flags=re.IGNORECASE)
    if matchobj is None:
        raise ValueError("Unable to extract delamination line number and index from filename %s" % (delamfilename))
    (linenumstring,indexstring)=matchobj.groups()

    linenum=int(linenumstring)
    index=None
    if indexstring is not None:
        index=int(indexstring[1:])
        pass
    return (linenum,index)

def last_kwarg_index(callnode):
    numcallparams = len(callnode)

    insertpos=numcallparams
    kwargs_node=callnode.find("dictargumentnode",recursive=False)
    if kwargs_node is not None:
        # Place insertpos immediately prior to kwargs_node
        insertpos=kwarg_nodes.index_on_parent
        pass
    return insertpos

def add_delams_to_bond_layers_calls(basename,basename_withpath,output_directory,phase,rb_tree):
    # bond_layers() calls not allowed inside functions or comprehensions or while loops...


    # identify delamination files
    delams = glob.glob(os.path.join(output_directory,"layerdelam_*.csv"))

    delams_pathsplit = [ os.path.split(delam) for delam in delams ]

    # each element in delams_withnumber_and_index is:
    #  (delamdir, delamfile, (linenumber, index_or_None) )
    delams_withnumber_and_index = [ (delamdir,delamfile,linenumber_and_index_from_filename(delamfile)) for (delamdir,delamfile) in delams_pathsplit ]

    # Convert to dictionary
    delams_by_linenumber = {}
    for (delamdir, delamfile, (linenumber, index_or_None) ) in delams_withnumber_and_index:
        if linenumber not in delams_by_linenumber:
            delams_by_linenumber[linenumber]=[]
            pass

        delams_by_linenumber[linenumber].append( (delamdir, delamfile, (linenumber, index_or_None) ))
        pass
    

    # As we iterate through the bond_layers calls,
    # we will remove the corresponding entries from the
    # delams_by_linenumber dictionary, so that any
    # remaining dictionary entries are excess entries. 
    
    # atomtrailers is a sequence of foo.bar.fubar(arg1,arg2)
    # where foo, bar and fubar would be a NameNode and (arg1,arg2) is a CallNode
    bondlayers_calls=rb_tree.find_all("atomtrailers",lambda node: node[0].value=="bond_layers" and isinstance(node[1],redbaron.CallNode))

    sourcelines=[]

    
    for bondlayers_call in bondlayers_calls:

        paramvals,paramindexes = extract_parameter_values(bondlayers_call,
                                                          bond_layers_params,
                                                          bond_layers_default_params)
        
        callnode=bondlayers_call[1]
        numcallparams = len(callnode)
        #insertpos = last_kwarg_index(callnode)

        #linearg = callnode.value.find("callargument",lambda node: node.target is not None and node.target.value=="delamo_sourceline")
        #if linearg is None:
        #    raise ValueError("add_delams_to_bond_layers_calls(): bond_layers() call without delamo_sourceline parameter")
        if paramindexes["delamo_sourceline"] < 0:
            raise ValueError("add_delams_to_bond_layers_calls(): bond_layers() call without delamo_sourceline parameter")
        linearg=paramvals["delamo_sourceline"]
        
        
        # Extract value of sourceline
        delamo_sourceline=int(linearg.value)
        sourcelines.append(delamo_sourceline)
        # have a kwarg for delamo_sourceline... remove it!
        # ... modifications moved to later, after we have finished processing,
        # so the indexes don't change
        sourcelinenode=callnode.value[paramindexes["delamo_sourceline"]]
        #callnode.value.remove(callnode.value[paramindexes["delamo_sourceline"]])

        #delaminationlistarg = callnode.value.find("callargument",lambda node: node.target is not None and node.target.value=="delaminationlist")
        #if delaminationlistarg is not None:
        #    raise ValueError("add_delams_to_bond_layers_calls(): Call to bond_layers() already has a delaminationlist specified")
        if paramindexes["delaminationlist"] >= 0:
            raise ValueError("add_delams_to_bond_layers_calls(): Call to bond_layers() (sourceline specified as %d) already has a delaminationlist specified" % (delamo_sourceline))
        
        
        # add a new kwarg for delaminationlist
        delaminationlist = redbaron.CallArgumentNode({"type": "call_argument",
                                                         "first_formatting": [],
                                                         "second_formatting":[],
                                                         "target": {
                                                             "type": "name",
                                                             "value":"delaminationlist"},
                                                         "value": {
                                                             "first_formatting": [],
                                                             "second_formatting": [],
                                                             "third_formatting": [],
                                                             "fourth_formatting": [],
                                                             "type": "list",
                                                             "value": []}})

        if delamo_sourceline in delams_by_linenumber:
            # Found a delamination .csv file 
            delam_tuples = delams_by_linenumber[delamo_sourceline]

            for (delamdir, delamfile, (linenumber, index_or_None) ) in delam_tuples:
                # add this to the delamination list

                delaminationlist.value.insert(len(delaminationlist.value),redbaron.StringNode({"type": "string", "value": "\""+os.path.join(delamdir,delamfile).encode('unicode_escape')+"\"", 'first_formatting': [], 'second_formatting': []}))
                pass
            
            # Remove from delams_by_linenumber to indicate we've handled these files
            del delams_by_linenumber[delamo_sourceline]
            pass

        # insertion of delaminationlist moved to end
        
        

        # !!!*** REPLACE kwarg for delamo_phase
        #phasearg = callnode.value.find("callargument",lambda node: node.target is not None and node.target.value=="delamo_phase")
        #if phasearg is None:
        #    raise ValueError("add_delams_to_bond_layers_calls(): bond_layers() call without delamo_phase parameter")
        if paramindexes["delamo_phase"] < 0:
            raise ValueError("add_delams_to_bond_layers_calls(): bond_layers() call without delamo_phase parameter")
        oldphasearg=callnode.value[paramindexes["delamo_phase"]]
        # removal of old phase argument moved below
        
        delamo_phase_arg=redbaron.CallArgumentNode({"type": "call_argument",
                                                        "first_formatting": [],
                                                    "second_formatting":[],
                                                    "target": {
                                                        "type": "name",
                                                        "value":"delamo_phase"},
                                                    "value": {
                                                        "first_formatting": [],
                                                        "second_formatting": [],
                                                        "type": "name",
                                                        "value": "\"%s\"" % (phase.encode('unicode_escape'))}})


        
        # have a kwarg for delamo_sourceline... remove it!
        # we do all mods together here
        # so the indexes don't change
        
        callnode.value.remove(sourcelinenode)
        callnode.value.insert(last_kwarg_index(callnode),delaminationlist)
        callnode.value.remove(oldphasearg)  # remove old phase argument
        callnode.value.insert(last_kwarg_index(callnode),delamo_phase_arg)

        
        pass

    pass




def _apply_delam_outline(delamo_sourceline,outline,rb_tree):
    ol_copy=outline.copy()

    # find bond_layers() with this delamo_sourceline
    bondlayers_call=rb_tree.find("atomtrailers",
                                 lambda node: (node[0].value=="bond_layers" and
                                               isinstance(node[1],redbaron.CallNode) and
                                               node[1].find("callargumentnode",
                                                            lambda argnode: (isinstance(argnode.target,redbaron.NameNode) and
                                                                             argnode.target.value=="delamo_sourceline"
                                                                             and isinstance(argnode.value,redbaron.IntNode)
                                                                             and int(argnode.value.value)==delamo_sourceline)) is not None))
    
    if bondlayers_call is None:
        raise ValueError("No bond_layers() call found with sourceline marked as %d" % (delamo_sourceline))

    delam_outlines_3d=bondlayers_call[1].find("callargumentnode",
                                              lambda argnode: (isinstance(argnode.target,redbaron.NameNode) and
                                                               argnode.target.value=="delam_outlines_3d"))
    if delam_outlines_3d is not None:
        if not isinstance(delam_outlines_3d.value,redbaron.ListNode) :
            raise ValueError("Existing delam_outlines_3d parameter to bond_layers call on line %d is not a list" % (delam_outlines_3d.absolute_bounding_box.top_left.line))
        delam_outlines_3d.value.append(ol_copy)
        pass
    else:
        callnode=bondlayers_call[1]
        insertpos=last_kwarg_index(callnode)
        do3d_node=redbaron.CallArgumentNode({"type": "call_argument",
                                             "first_formatting": [],
                                             "second_formatting":[],
                                             "target": {
                                                 "type": "name",
                                                 "value":"delam_outlines_3d"},
                                            "value": {
                                                "first_formatting": [],
                                                "second_formatting": [],
                                                "third_formatting": [],
                                                "fourth_formatting": [],
                                                "type": "list",
                                                "value": [],
                                            }})
        do3d_node.value.append(ol_copy) # Add in outlin
        callnode.value.insert(insertpos,do3d_node)
 
        
        pass
    
    pass

# NOTE **** apply_delam_outlines is only used by delamo_test_processor,
# which should probably be updated to use the more modern routines !!!***
def apply_delam_outlines(outlines,rb_tree):
    """ Given text of outline declarations 
    and a processed rb_tree, add the outlines
    to the bond_layers calls in the rb_tree"""
    
    outlines_tree=RedBaron(outlines)
    for outline_node in outlines_tree:
        if (isinstance(outline_node,redbaron.EndlNode) or
            isinstance(outline_node,redbaron.CommentNode)):
            # blank line or comment... ignore
            pass
        elif (isinstance(outline_node,redbaron.AtomtrailersNode) and 
              isinstance(outline_node.value[0],redbaron.NameNode) and
              outline_node.value[0].value=="delam_outline" and
              isinstance(outline_node.value[1],redbaron.CallNode)):
            # A delam_outline() call
            delamo_sourceline=None
            outline=None
            for argnode in outline_node.value[1]:
                if (isinstance(argnode,redbaron.CallArgumentNode) and
                    isinstance(argnode.target,redbaron.NameNode) and
                    argnode.target.value=="delamo_sourceline" and
                    isinstance(argnode.value,redbaron.IntNode)):
                    # Got delamo_sourceline argument
                    delamo_sourceline=int(argnode.value.value)
                    pass
                elif (isinstance(argnode,redbaron.CallArgumentNode) and
                    isinstance(argnode.target,redbaron.NameNode) and
                    argnode.target.value=="outline"):
                    outline = argnode.value
                    pass
                else:
                    raise ValueError("Unknown or invalid delam_outline() parameter %s at line %d" % (argnode.dumps(),argnode.absolute_bounding_box.top_left.line))
                pass

            if delamo_sourceline is None:
                raise ValueError("delamo_sourceline parameter not given in delam_outlines at line %d" % (outline_node.absolute_bounding_box.top_left.line))

            if outline is not None:
                _apply_delam_outline(delamo_sourceline,outline,rb_tree)
                
                pass
            
            pass
        else:
            raise ValueError("Unknown or invalid syntax element in delam_outlines at line %d" % (outline_node.absolute_bounding_box.top_left.line))
        pass
    pass

def extract_parameter_values(function_call, params, default_params):
    """ Extract parameter values from a function call (i.e. atomtrailers
    ending in parameters) in the RedBaron tree. Params is a list of 
    required parameters (strings). default_params is an OrderedDict represeting
    optional parameter names (strings) and RedBaron representations 
    of their default values. 

    Returns a dictionary by parameter name of RedBaron representations
    of their parameter values and a dictionary by parameter name of the index of the parameters in the function call, or -1 for parameters not present"""
    
    # reset parameter and index values 
    paramvals = { paramname: None for paramname in params }
    indexvals = { paramname: -1 for paramname in params }
    paramvals.update(default_params)
    for default_param_name in default_params:
        indexvals[default_param_name]=-1
        pass
    

    # Extract parameter values from function call
        
    
    paramsnode=function_call[-1] # end of the atomtrailers
    argnum=0
    # Extract argument values
    for argumentnode in paramsnode.value:
        assert(isinstance(argumentnode,redbaron.nodes.CallArgumentNode))
        if argumentnode.target is not None:
            if argumentnode.target.value not in paramvals:
                raise ValueError("Argument %s unknown at line %d" % (argumentnode.target.value, argumentnode.absolute_bounding_box.top_left.line))
            
            paramvals[argumentnode.target.value]=argumentnode.value.copy()
            indexvals[argumentnode.target.value]=argnum
            pass
        elif argnum < len(params):
            paramvals[params[argnum]]=argumentnode.value.copy()
            indexvals[params[argnum]]=argnum
            pass
        elif argnum < len(params) + len(default_params):
            paramname=list(default_params.items())[argnum-len(params)][0]
            paramvals[paramname]=argumentnode.value.copy()
            indexvals[paramname]=argnum
            pass
        else:
            raise ValueError("Too many parameters in call %s line %d" % (function_call.dumps(), argumentnode.absolute_bounding_box.top_left.line))
            
        argnum+=1
        pass
    
    # Make sure all required parameters are specified
    for paramname in paramvals:
        if paramvals[paramname] is None:
            raise ValueError("Parameter %s not specified in call %s line %d" % (paramname,function_call.dumps(),argumentnode.absolute_bounding_box.top_left.line))
        pass
    return (paramvals,indexvals)

def substitute_function(rb_tree,fun_def,atomtrailers_nodematch):
    """ Substitute content for a function, 
    replacing params referenced within its content. 

    Note that the substitution is not entirely equivalent to calling
    the function: There is no new scope with the substitution, so 
    variable assignments within the function affect the calling scope. 

    
    rb_tree is the tree in which to substitute 
    atomtrailers_nodematch is a function that checks whether a particular 
    atomtrailers node in rb_tree is the function we want to substitute
    params is a list of parameter names. If None, use the 
    function name in fun_def
"""
    # turn fun_def into a RedBaron tree if it isn't already
    if isinstance(fun_def,basestring):
        fun_def=RedBaron(fun_def)
        pass

    if not isinstance(fun_def,redbaron.nodes.DefNode):
        
        def_node_indexes=[ index for index in range(len(fun_def)) if isinstance(fun_def[index],redbaron.nodes.DefNode) ]

        if len(def_node_indexes) != 1:
            raise ValueError("Found %d function definitions in fun_def" % (len(def_node_indexes)))
        
        def_node_index=def_node_indexes[0]

        fun_def = fun_def[def_node_index]
        pass
    
    fun_name=fun_def.name

    # If no node-matcher is specified, just find calls to this function
    if atomtrailers_nodematch is None:
        atomtrailers_nodematch = lambda node: (isinstance(node[0],redbaron.NameNode) and
                                               node[0].value==fun_name and
                                               isinstance(node[1],redbaron.CallNode))
        pass


    # assemble params and
    # default_params
    started_default_params=False

    params= [ ]  # just list of param names
    default_params=collections.OrderedDict()
    
    for argument in fun_def.arguments:
        assert(isinstance(argument.target,redbaron.nodes.NameNode))
        if argument.value is None:
            # parameter with no default value
            if started_default_params:
                raise ValueError("Parameter %s with no default value follows parameter with default value" % (argument.target.value))
            params.append(argument.target.value)
            pass
        else:
            # parameter with a default value
            started_default_params=True
            default_params[argument.target.value]=argument.value.copy()
            pass
        pass

    
    # Find calls to specified function
    function_calls = rb_tree.find_all("atomtrailers",
                                      atomtrailers_nodematch)

    # Iterate over matching function calls
    
    for function_call in function_calls:


        paramvals,paramindexes = extract_parameter_values(function_call, params, default_params)
        
        # OK.. Got all parameters. Now substitute them
        rb_code_indentation=fun_def.value[0].indentation
        call_indentation=function_call.indentation

        indentamount=len(call_indentation)-len(rb_code_indentation)

        rb_code=fun_def.value.copy()
        if indentamount > 0:
            rb_code.increase_indentation(indentamount)
            pass
        elif indentamount < 0:
            rb_code.decrease_indentation(-indentamount)
            pass
        
        
        for paramname in paramvals: # Each parameter
            paraminsts=rb_code.find_all("namenode",value=paramname)
            for paraminst in paraminsts:  # Each instance of each parameter
                paramparent=paraminst.parent
                paramparentindex=paraminst.index_on_parent

                # replace instance with value copy
                if paramparentindex is not None:
                    # if it is an indexed member, unless
                    # it is an atomtrailersnode where we are
                    # not the first emember
                    if (paramparentindex == 0 or
                        not isinstance(paramparent,redbaron.nodes.AtomtrailersNode)):
                        paramparent[paramparentindex] = paramvals[paramname].copy()
                        pass
                    pass
                else :
                    # if it is some named attribute, need to figure out which attribute
                    attrnames=[ attrname for attrname in dir(paramparent) if not attrname.startswith("_") and attrname != "next_rendered" and attrname != "absolute_bounding_box" and attrname != "root" and getattr(paramparent,attrname) is paraminst ]
                    assert(len(attrnames)==1)
                    setattr(paramparent,attrnames[0],paramvals[paramname].copy())
                    pass
                
                pass
            pass

        fc_parent=function_call.parent
        fc_parentindex=function_call.index_on_parent

        # replace simple.laminate call with generated code
        # fc_parent[fc_parentindex]=rb_code
        del fc_parent[fc_parentindex]
        index=fc_parentindex

        #for o in rb_code:
        #    print("Node to insert: %s" % (repr(o)))
        #    pass
        
        
        fc_parent.insert_multiple(index,rb_code)
        #fc_parent.extend(rb_code)
        #for o in rb_code:
        #    fc_parent.insert(index,o);
        #    index+=1

        
        
        #for code_el in rb_code:
        #    fc_parent.insert(index,code_el.copy())
        #    index+=1
        #    pass
        
        pass
    pass
        



def unwrap_simple(rb_tree):

    if simple.__file__.endswith(".pyc"):
        simple_file=os.path.splitext(simple.__file__)[0]+".py"
        pass
    else:
        simple_file=simple.__file__
        pass
    
    simple_module_tree = tree_from_file(simple_file)
    
    #code=
    def_laminate_nodematch=lambda node: node.name=="laminate"
    
    def_laminate_nodes = simple_module_tree.find_all("defnode",def_laminate_nodematch)

    if len(def_laminate_nodes) != 1:
        raise RuntimeError("Did not find single function definition of \"laminate\" in %s" % (simple_file))

    def_laminate_node=def_laminate_nodes[0]
    
    delamo_simple_nodematch=lambda node: (isinstance(node[0],redbaron.NameNode) and
                                          node[0].value=="simple" and
                                          isinstance(node[1],redbaron.NameNode) and
                                          node[1].value=="laminate" and
                                          isinstance(node[2],redbaron.CallNode))
    
    
    substitute_function(rb_tree,def_laminate_node,delamo_simple_nodematch)
    pass

def tree_from_str(fstr):
    rb_tree=RedBaron(fstr)

    # Verify that parsing and rendering round-trip does not change input
    # ... if this fails it is a baron/RedBaron bug. 
    assert(fstr==rb_tree.dumps())
    return rb_tree

    
def tree_from_file(fname):
    with open(fname,"r") as fh:
        fstr=fh.read()
        fh.close()

        # First parse with Python ast because of better error messages.
        ast.parse(fstr,fname)
        
        rb_tree=tree_from_str(fstr)
        return rb_tree
    pass

def tree_to_file(rb_tree,fname):
    with open(fname,"w") as fh:
        fh.write(rb_tree.dumps())
        fh.close()
        pass
    pass



def original_to_predamage_file(fname,phase):
    
    with open(fname,"r") as fh:
        fstr=fh.read()
        fh.close()
                       
        basename_withpath=os.path.splitext(fname)[0]                       
        (output_directory,basename)=os.path.split(basename_withpath)
                       
        rb_tree=original_to_predamage_str(basename,basename_withpath,output_directory,phase,fstr)
        return rb_tree
    pass

def original_to_predamage_str(basename,basename_withpath,output_directory,phase,fstr):
    rb_tree=RedBaron(fstr)

    # Verify that parsing and rendering round-trip does not change input
    # ... if this fails it is a baron/RedBaron bug. 
    assert(fstr==rb_tree.dumps())

    return original_to_predamage_tree(basename,basename_withpath,output_directory,phase,rb_tree)

def original_to_predamage_tree(basename,basename_withpath,output_directory,phase,rb_tree):

    unwrap_simple(rb_tree)
    
    troubleshoot_bond_layers(rb_tree)    
    unwrap_loops(rb_tree)

    annotate_bond_layers_calls(basename,basename_withpath,output_directory,phase,rb_tree)

    return rb_tree


def predamage_preclean(basename,basename_withpath,output_directory,phase):
    # Remove layerboundary outputs from this phase, if any
    for old_layerboundary_stl in glob.glob(os.path.join(output_directory,"layerboundary_%s_*.stl" % (phase))):
        print("  rm \'%s\'" % (old_layerboundary_stl))
        os.remove(old_layerboundary_stl)
        pass

    # remove any recorded delaminations (should regenerate after PREDAMAGE
    # completed)

    # ***!!!! Should postpone this to post-cleaning and not do anything
    # unless the line numbers have changed. 
    for old_layerdelam_csv in glob.glob(os.path.join(output_directory,"layerdelam_*.csv")):
        print("  rm \'%s\'" % (old_layerdelam_csv))
        os.remove(old_layerdelam_csv)
        pass
    
    
    pass

def predamage_to_damage_tree(basename,basename_withpath,output_directory,phase,rb_tree):
    add_delams_to_bond_layers_calls(basename,basename_withpath,output_directory,phase,rb_tree)

    return rb_tree


def run_damage_script(basename,basename_withpath,output_directory,phase,rb_tree):
    # rb_tree is from PRIOR phase
    from . import process

    # Extract damage_script parameter as in process.get_process_phases()
    output_filenames_call = process.get_output_filenames_call(rb_tree)

    paramvals,paramindexes = extract_parameter_values(output_filenames_call,
                                                      process.output_filenames_params,
                                                      process.output_filenames_default_params)

    apply_damage_script_rb=paramvals["apply_damage_script"]
    if isinstance(apply_damage_script_rb,redbaron.NameNode) and apply_damage_script_rb.value=="None":
        raise ValueError("apply_damage_script parameter not specified in call to output_filenames()")
    elif isinstance(apply_damage_script_rb,redbaron.StringNode):
        apply_damage_script = apply_damage_script_rb.to_python() 
    else:
        raise ValueError("Base script should have call to output_filenames() with apply_damage_script set to a string")
        
    directory = os.path.split(basename_withpath)[0]
    if len(directory) > 0:
        print("  cd %s" % (directory))
        pass
    params = [sys.executable]
    params.append(apply_damage_script)

    print("  %s" % (" ".join(params)))

    if len(directory)==0:
        directory=None
        pass
    
    try:
        subprocess.check_call(params,cwd=directory)
        pass
    except subprocess.CalledProcessError:
        print("")
        print("Error executing apply damage script")
        pass

    pass

    


