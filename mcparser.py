"""McStas parser for blinggen"""
import re
import numpy as np 
from math import cos, sin, pi, sqrt, log, atan
from collections import namedtuple
import itertools
import pycparser


class Component:
    """Class to hold the information for a component"""
    def __init__(self, ctype, name, pos, rot, **kwargs):
        self.type = ctype
        self.name = name
        self.pos = pos
        self.rot = rot
        for key, value in kwargs.items():
            setattr(self, key, value)


    def unfoldcopy(self, copied_comp, n=1):
        """Make a copy of Component copied_comp in component comp, with the proper
        placement and rotation. Notice that both must be declared"""
        exclude= ["rot", "pos", "name"]
        for var in vars(copied_comp):
            if var not in exclude:
                vars(self)[var] = vars(copied_comp)[var]
        self.name = "COPY_"+copied_comp.name+str(n) 



def parsearguments(rawargs):
    """parse from string rawargs of a component in order to get the relevant argument values"""
    keywords = ("w1", "w2", "h1", "h2", "l", "mleft", "mright", "mtop", "mbottom", "m", 
                "xwidth", "yheight", "dimensionsAt", "linxw", "linyh", "loutxw", "loutyh")
    
    # Change all commas inside brackets to spaces. This helps numpy and avoids confusion
    # This is done with regular expressions, not sure abou robustness!"
    O_ONE = re.compile('(\[[!\]]*),([!\]]*\])')
    arguments = O_ONE.sub("\g<1> \g<2>", rawargs)

    args = arguments.split(",")
    inp = {}   # empty dictionary
    for arg in args:
        kw = arg.split("=")[0].strip()
        if kw in keywords:
            val = arg.split("=")[1].strip()
            inp[kw] = val
    return inp


def parseloc(parsedstr, varlib = {}):
    """ Auxiliar parsing from McStas."""
    if parsedstr is None:
        return(0, 0, 0, None)  #TODO: Maybe we don't actually want defaults here
    else:  pass  
    coords = parsedstr.split("(")[1].split(")")[0]
    coords = coords.split(",")
    xyz = [0, 0, 0]
    for i in [0,1,2]:
        try:
            xyz[i] = float(coords[i])
        except:
            xyz[i] = parseop(coords[i], varlib)
    (x, y, z) = xyz[0:3]
    # HANDLE THE RELATIVE TO part.
    if "RELATIVE" in parsedstr:
        ref = parsedstr.split("RELATIVE")[1].strip()
    elif "ABSOLUTE" in parsedstr:
        ref = None
    return (z, x, y, ref)  # Yes, this is on purpose to have the MCNP guide pointing to +X


def parseAT(McStas_AT, varlib={}):
    """ Parse the AT statement, and return a displacement array """
    coords = parseloc(McStas_AT, varlib)
    location = namedtuple("location", "x y z reference")
    result = location(coords[0], coords[1], coords[2], coords[3])
    return result

    
def parseROT(McStas_ROT, varlib={}):
    """ Parse the ROTATED statement, and return a rotation matrix"""
    axis = parseloc(McStas_ROT, varlib)
    rotation = namedtuple("rotation", "Ux Uy Uz reference")
    rotang = rotation(axis[0], axis[1], axis[2], axis[3])
    Ux = rotang.Ux*pi/180
    Uy = rotang.Uy*pi/180
    Uz = rotang.Uz*pi/180
    # Define the rotation matrixes and multiply
    rotx = np.array([[1, 0, 0], [0, cos(Ux), -sin(Ux)], [0, sin(Ux), cos(Ux)]])
    roty = np.array([[cos(Uy), 0, -sin(Uy)], [0, 1, 0], [sin(Uy), 0, cos(Uy)]])
    rotz = np.array([[cos(Uz), sin(Uz), 0], [-sin(Uz), cos(Uz), 0], [0, 0, 1]])
    rot = np.matmul(rotx, roty)
    rot = np.matmul(rot, rotz)
    return rot, rotang[3]


def getcomps(infile, varlib={}):
    """ Get the components of a guide from an *instr infile"""
    comps = []
    with open(infile, "r") as instrfile:
        rawdata = instrfile.readlines()
    # Remove comments and blank lines, we don't want to deal with them later
    data = []
    for i, line in enumerate(rawdata):
        if not line:
            continue
        newline = line.split("//")[0]  # Remove inline comment
        if newline.strip().strip("\t") == "":
            continue
        elif "/*" in newline:
            if newline.strip()[0:2] == "/*":
                continue
        data.append(newline)
    # Split the file into strings containing a component each
#    print(data)
    spltline = []
    for i, line in enumerate(data):
        if re.match("COMPONENT", line.strip()):
            spltline.append(i)
    print(spltline)
    rawcomps = [data[i:j] for i, j in zip(spltline, spltline[1:]+[len(data)])]
    # Curate this
    components = []
    for rawcomp in rawcomps:
        ROTline = ""  # Default if no line defined
        for line in rawcomp:
#            if line.strip().strip("\t") == "":
#                print("skipping")  #I don't think this is needed anymore
#                continue
            if re.match("COMPONENT", line.strip()):
                print(line)
                name = line.split("=")[0].split()[1]
                comptype = line.split("=")[1].split("(")[0]  
                if "(" in line:
                    rawargs = line.split("(")[1].split(")")[0]  # To enclose all
                else:
                    rawargs = ""  #TODO: Come up with something less horrible. 
            elif re.match("AT", line.strip()) is not None:  
                ATline = line.strip()
            elif re.match("ROTATED", line.strip()) is not None:
                ROTline = line.strip()
            else:
                rawargs = rawargs+line.split(")")[0]
        AT = parseAT(ATline, varlib)
        if not ROTline:  # From McStas 2.6.1 manual sect. 4.3.6
            ROT = [np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]), AT.reference]
        else:
            ROT = parseROT(ROTline, varlib)
        args = parsearguments(rawargs)
        for key in args:
            args[key] = parseop(args[key], varlib)
        comp = Component(comptype, name, AT, ROT, **args)
        components.append(comp)
    return components 


def doBinOp(Node, varlib=None):
    """Do a BinaryOp Node and return the result"""
    if isinstance(Node, (pycparser.c_ast.Constant, pycparser.c_ast.UnaryOp, pycparser.c_ast.ID)):
        return getconstval(Node, varlib)
    mut = lambda  a, b : a*b
    div = lambda  a, b : a/b
    add = lambda  a, b : a+b
    subs = lambda a, b : a-b
    operator = {"/": div, "*": mut, "+": add, "-": subs}
    result = operator[Node.op](doBinOp(Node.left, varlib), doBinOp(Node.right, varlib))
    return result

def getconstval(Node, varlib=None):
    """Get constant value from a Constant, ID or UnaryOpo (sub)node. 
    ID nodes need a varlib to refer"""
    if isinstance(Node, pycparser.c_ast.Constant):
        if Node.type in ("double", "float"):
            return float(Node.value)
        if Node.type in ("int"):
            return int(Node.value)
    elif isinstance(Node, pycparser.c_ast.UnaryOp):
        if Node.op == '-':
            return(-getconstval(Node.expr))
        else:
            print("Unimplemented operator")
            return None
    elif isinstance(Node, pycparser.c_ast.ID):
        try:
            return varlib[Node.name]
        except:
            raise ValueError("Variable {0} not in library {1}\n".format(Node.name, varlib))

def getfuncvars(Node, varlib=None):
    """Get variables declared or initialized in child node of function
    Only works with = assignments for now, probably very limited and
    will need expanding"""
    if Node[1].op != "=":
        print("Operation not implemented")
        return None
    name = Node[1].lvalue.name
    if isinstance(Node[1].rvalue, (pycparser.c_ast.Constant, pycparser.c_ast.UnaryOp,
                                   pycparser.c_ast.ID)):
        value = getconstval(Node[1].rvalue, varlib)
    elif isinstance(Node[1].rvalue, pycparser.c_ast.BinaryOp):
        value = doBinOp(Node[1].rvalue, varlib)
    return name, value
    

def getvariables(infile):
    """ 
    Get the variables from a curated instr file. This is an instr file
    where the DECLARE instruction, and the component definition
    are removed. a main() function must exist as well"""
    ast = pycparser.parse_file(infile, use_cpp=True)
    variables = []
    for ele in ast.ext:
        if isinstance(ele, pycparser.c_ast.Decl):
            variables.append(ele)
    varlib = {}
    for var in variables:  # First loop over declared
        if var.init is None:
            continue
        if isinstance(var.init, pycparser.c_ast.UnaryOp):  # negative constant?
            varlib[var.name] = getconstval(var.children()[1][1].expr)
        if isinstance(var.init, pycparser.c_ast.Constant):  # constant?
            varlib[var.name] = getconstval(var.children()[1][1])
        if isinstance(var.init, pycparser.c_ast.InitList):  # array
            varlib[var.name] = getarraynode(var)
        
    for var in ast.ext:  # Now loop over all ext
        if isinstance(var, pycparser.c_ast.FuncDef): # function
            for child in var.children()[1][1].children():
                varlib[getfuncvars(child, varlib)[0]] = getfuncvars(child, varlib)[1]
    return varlib


def getarraynodedim(Node):
    """ Get the array dimensions from a node. Called recursively"""
    
    if hasattr(Node.children()[0][1], 'dim'):
        a = int(Node.children()[0][1].dim.value)
        b = getarraynodedim(Node.children()[0][1])      
        b.append(a)
        return b
    else:
        return []
    

def getarraynode(Node):
    """ Get the array values from a node. Arrays are assumed to be numerical 
    and will fail otherwise"""
    dims = getarraynodedim(Node) 
    dims.reverse()  # TODO getarraydim() should probably just return it right
    A = np.zeros(dims)
    for idx in itertools.product(*[range(s) for s in dims]):
        Subnode = Node.children()[1][1]
        for i in idx:
            Subnode = Subnode.exprs[i]
        A[idx] = getconstval(Subnode)
    return A


def parsearray(arrstr, varlib):
    """ get the value of array element defined in arrstr string from varlib"""
    try:
        result = float(arrstr)
        return result
    except: 
        if re.match("\"", arrstr) is not None:
            return arrstr.strip("\"")
        arr = re.split('\[|\]| ', arrstr)  # Notice this also strips arrstr
        while '' in arr:
            arr.remove('')  # Null strings may be generated by re, don't want them

        if arr[0] not in varlib.keys():
            raise ValueError("Variable {0} unknown".format(arr[0]))
            return None
        elif len(arr) == 1:
            return varlib[arr[0]]  
        else: 
            idx = [int(i) for i in arr[1:]]
            result = varlib[arr[0]]
            for i in idx:
                result = result[i]
            return result


def parseop(opstr, varlib):
    """ Get the result of the operation defined by opstr, with variables
    from varlib. Return another string with highest priority operation done.
    If only one or no operation, return value.
    Only basics operators implemented """
    if not opstr:
        return 0  #Return 0 for void arguments for instance in case of A=-B
    #TODO: We probably should define a lambda function for the ops...
    hi_ops = "\*|\/"  # Hi priority ops
    low_ops = "\+|\-"  # Low priority ops
    
    if not re.findall(low_ops, opstr):
        if not re.findall(hi_ops, opstr):
            return parsearray(opstr, varlib)  # No actual operations
        else:
            op = re.findall(hi_ops, opstr)[-1]
            token = re.split(hi_ops, opstr)
            token1 = token[-1]
            token0 = opstr[:-len(token1)-1]  # TODO: Check this hack....
            if op == "*":
                return parseop(token0, varlib) * parseop(token1, varlib)
            if op == "/":
                return parseop(token0, varlib) / parseop(token1, varlib)
    else:
        op = re.findall(low_ops, opstr)[-1]
        token = re.split(low_ops, opstr)
        token1 = token[-1]
        token0 = opstr[:-len(token1)-1]  # TODO: Check this hack....
        if op == "+":
            return parseop(token0, varlib) + parseop(token1, varlib)
        if op == "-":
            return parseop(token0, varlib) - parseop(token1, varlib)


def composepos(pos, ref):
    """compose position pos with Component ref, to change its reference
    to that of ref"""
    if ref.rot[1] not in ["ABSOLUTE", None]:
        print("Can not position an instrument with not oriented reference")
        return pos
    location = namedtuple("location", "x y z reference")
    TR = np.linalg.inv(ref.rot[0])

    x = ref.pos.x + pos.x*TR[0][0]+pos.y*TR[0][1]+pos.z*TR[0][2]
    y = ref.pos.y + pos.x*TR[1][0]+pos.y*TR[1][1]+pos.z*TR[1][2]
    z = ref.pos.z + pos.x*TR[2][0]+pos.y*TR[2][1]+pos.z*TR[2][2]
    result = location(x, y, z, ref.pos.reference)
    return result 

def composerot(rot, ref):
    """compose rotation rot with Component ref"""
    TR = np.matmul(ref.rot[0], rot[0])
    result = [TR, ref.rot[1]]
    return result


def makeinstr(comps):
    """Make a coherent instrument from array of Components comp.
    This implies unfolding the copies, and setting al positions and
    rotation to absolute values"""
    for i, comp in enumerate(comps):
        if comp.rot[1] not in [None, "ABSOLUTE"]:
            if comp.rot[1] == "PREVIOUS":
                if i == 0:  #special case for undefined 1st component rotation
                    comp.rot = (comp.rot[0], None)
                else:
                    comp.rot = composerot(comp.rot, comps[i-1])
            else:
                ref = [r for r in comps if comp.rot[1] == r.name][0] 
                comp.rot = composerot(comp.rot, ref)
        if comp.pos.reference not in [None, "ABSOLUTE"]:
            if comp.pos.reference == "PREVIOUS":
                comp.pos = composepos(comp.pos, comps[i-1])
            else:
                ref = [r for r in comps if comp.pos.reference == r.name][0] 
                comp.pos = composepos(comp.pos, ref)


def getbinvariables(infile, nsegments):
    """ 
    Gets the variables from the compiled McStas file infile,
    assuming there are nsegments guide segments. Currently unused,
    maybe not needed at all.
     """
    import ctypes
    import itertools
    instr = ctypes.CDLL(infile)
    guideflat = ctypes.cast(instr.guide, ctypes.POINTER(ctypes.c_double*10*nsegments))
    mflat = ctypes.cast(instr.m, ctypes.POINTER(ctypes.c_double*4*nsegments))
    w1 = [guideflat.contents[i][0] for i in range(nsegments)]
    w2 = [guideflat.contents[i][2] for i in range(nsegments)]
    h1 = [guideflat.contents[i][1] for i in range(nsegments)]
    h2 = [guideflat.contents[i][3] for i in range(nsegments)]
    l = [guideflat.contents[i][4]*guideflat.contents[i][9] for i in range(nsegments)]
    r = guideflat.contents[:][8]
    mu = mflat.contents[2::4]
    md = mflat.contents[3::4]
    ml = mflat.contents[0::4]
    mr = mflat.contents[1::4]
    m = [[mflat.contents[i][2], mflat.contents[i][3], mflat.contents[i][0], 
         mflat.contents[i][1]] for i in range(nsegments)]
    h = [a for b in zip(h1, h2) for a in b]  
    w = [a for b in zip(w1, w2) for a in b]  

    return {'m': m, 'l': l, 'h': h, 'w': w, 'r': r}


