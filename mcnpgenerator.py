"""BeamLINe Geometry GENerator"""
import re
import linecache
import numpy as np 
# from pylab import *
from math import cos, sin, pi, sqrt, log, atan
from array import array
from collections import namedtuple
import itertools
import pycparser


def BoundBox(Up, Down, Left, Right):
    S=[]
    S.append("999990 PZ {0}\n".format(Up))
    S.append("999991 PZ -{0}\n".format(Down))
    S.append("999992 PY -{0}\n".format(Left))
    S.append("999993 PY {0}\n".format(Right))
    return S

def calc_ellipse(x1, x2, y, posy):
    """Calculate the center, minor axis, and major axis of the sphere, given by focus
    placement x1 and x2, ellipse width y, and where that ellipse y is measured, posy."""
    f = abs(x1-x2)  # focal distance
    a = sqrt(y**2 + (x1-posy)**2) + sqrt(y**2 + (x2-posy)**2)
    b = sqrt(a**2- f**2)
    center = (x1 +x2)/2
    major = a/2
    minor = b/2
    return{"center": center, "major": major, "minor": minor}

def AddSurf(Type, N, TR, params, Slist):
    """
    Adds a surface of type Type with number N and Transform card TR to Slist, as defined by 
    params. Valid types are:
    PX,PY,PZ: The obvious
    P: Params: Point and normal 
    C/X, C/Y, C/Z: Params are the Origin and radius (like MCNP card)
    SQ: A B C D E F G, X0, Y0, Z0. (Like MCNP Card) 
    """
    if Type in ["PX","PY","PZ"]:
        Slist.append("{0} {1} {2} {3}".format(N,TR,Type,params)) 
    elif Type in ["P"]:  
        if len(params)<6:
            print("not enough parameters for general plane calculation in this function")
            #TODO: Implement an alternative way to do it with just 4 numbers, i.e: Normal+displacement
            return N
        Mag=sqrt(sum(i**2 for i in params[3:6]))
        Normal=[i/Mag for i in params[3:6]]
        D=sum(i*j  for i,j  in zip(Normal,params[0:3]))
        Slist.append("{0} {1} P {2:8F} {3:8F} {4:8F} {5:6E}".format(N,TR,Normal[0],Normal[1],Normal[2],D))
    elif Type in ["C/X","C/Y","C/Z"]:
        if len(params)<3:
            print("not enough parameters for general plane calculation in this function")
            return N
        Slist.append("{0} {1} {2} {3:6E} {4:.5f} {5:.5f}".format(N,TR,Type,params[0],params[1],params[2]))
    elif Type in ["TX","TY","TZ"]:
        if len(params)<6:
            print("not enough parameters to describe an Axis-parallel Torus")
            return N
        Slist.append("{0} {1} {2} {3} {4} {5} {6} {7} {8}".format(N,TR,Type,params[0],params[1],params[2],
        params[3],params[4],params[5]))    
    elif Type in ["SQ"]:
        Slist.append("{0} {1} {2} {3} {4} {5} {6} {7}".format(N, TR, Type, 
                     params[0], params[1], params[2], params[3], params[4]))
        
        Slist.append("      {0} {1} {2} {3} {4}".format(params[5], params[6], 
                     params[7], params[8], params[9]))
    else:
        print("Unknown or unimplemented surface type")
        return N
    N+=1
    return N


def AddCell (N, Mat, Ro, InSurf, ExSurf, Clist, IMP=1):
    """
    Adds a cell number N of Material Mat with density Ro to Clist. Surfaces indicated by InSurf 
    and ExSurf array/list, in the form (InSurf[0] InSurf[1]...InSurf[n] (ExSurf[0]:ExSurf[1]..ExSurf[n])
    THE SURF ARRAY MUST HAVE THE +/- SIGN!
    """
    import textwrap
    if Ro==0:
        if Mat!=0:
            print("WARNING: Zero density given, changing material to void.")
        Mat=0

    C="{0} {1} ".format(N, Mat)
    if Mat!=0:
        C=C+"{0} ".format(Ro)
    C=C+"("
    
    for S in InSurf:
        C=C+"{0} ".format(S)
    if not not ExSurf:
        C=C+"("  
        for S in ExSurf:
            C=C+"{0}:".format(S)
        C=C.strip(":") # REMOVE THE FINAL COLON
        C=C+")"
    C=C+") IMP:N={0}".format(IMP)
    
    C=textwrap.fill(C,75,subsequent_indent="      ") # Yeah, could be 80 but I do not feel like taking a risk here
    Clist.append(C)
    N+=1
    return N

def CreateMirrors(Comp, N, substhick=1):
    """Create and return the surface cards for the mirrors of component Comp. 
       N is the number of component, used for surface and TR numbering"""
    # Coating thickness dictionary. 
    coat_thick=	{
        1.5: 7E-5,
        2: 1.3E-4,
        2.5: 2.1E-4,
        3: 3.4E-4,
        4: 6E-4,
        5: 9E-4,
        6: 1.6E-3
    }

    # Previous checks
    Surf=[]
    Dump=[]  #Trash array for surfaces that will NOT be written
    # Build the mirror thickness variable:
    if hasattr(Comp, "m"):  # Single m-value for all surfaces
        t = [coat_thick[Comp.m]]*4
    else:  # Will fail if mtop, mdown, etc not defined
        t = [coat_thick[vars(Comp)[place]] for place in ["mtop", "mbottom", "mleft", "mright"]]
    snum = 100*N
    l = 100*Comp.l  
    snum = AddSurf("PX", snum, N, 0, Surf)
    snum = snum+10-snum%10
    if Comp.type.strip() in ["Guide", "Guide_gravity", "Guide_gravity_polar"]:  #Straight guide
        h1 = 100*Comp.h1
        h2 = 100*Comp.h2
        w1 = 100*Comp.w1
        w2 = 100*Comp.w2
        if Comp.h1 == Comp.h2:
            snum = AddSurf("PZ", snum, N, h1/2, Surf)
            snum = AddSurf("PZ", snum, N, h1/2+t[0], Surf)
            snum = AddSurf("PZ", snum, N, h1/2+t[0]+substhick, Surf)
            snum = snum+10-snum%10
            snum = AddSurf("PZ", snum, N, -h1/2, Surf)
            snum = AddSurf("PZ", snum, N, -h1/2-t[1], Surf)
            snum = AddSurf("PZ", snum, N, -h1/2-t[1]-substhick, Surf)
            snum = snum+10-snum%10
        else: 
            alpha = atan((Comp.h2-Comp.h1)/(2*l))
            snum = AddSurf("P", snum, N, [0, 0, h1/2, -sin(alpha), 0 ,cos(alpha)], Surf)
            snum = AddSurf("P", snum, N, [0, 0, h1/2+t[0], -sin(alpha), 0, cos(alpha)], Surf)
            snum = AddSurf("P", snum, N, [0, 0, h1/2+t[0]+substhick, -sin(alpha), 0, 
                           cos(alpha)], Surf)
            snum = snum+10-snum%10
            snum = AddSurf("P", snum, N, [0, 0, -h1/2, sin(alpha), cos(alpha), 0], Surf)
            snum = AddSurf("P", snum, N, [0, 0, -h1/2+t[0], sin(alpha), cos(alpha), 0], Surf)
            snum = AddSurf("P", snum, N, [0, 0, -h1/2+t[0]+substhick, sin(alpha), cos(alpha),
                           0], Surf)
            snum = snum+10-snum%10
        if Comp.w1 == Comp.w2:
            snum = AddSurf("PY", snum, N, -w1/2, Surf)
            snum = AddSurf("PY", snum, N, -(w1/2+t[2]), Surf)
            snum = AddSurf("PY", snum, N, -(w1/2+t[2]+substhick), Surf)
            snum = snum+10-snum%10
            snum = AddSurf("PY", snum, N, w1/2, Surf)
            snum = AddSurf("PY", snum, N, w1/2+t[3], Surf)
            snum = AddSurf("PY", snum, N, w1/2+t[3]+substhick, Surf)
            snum = snum+10-snum%10
        else: 
            alpha = atan((Comp.h2-Comp.h1)/(2*l))
            snum = AddSurf("P", snum, N, [0, 0, -w1/2, sin(alpha), cos(alpha), 0], Surf)
            snum = AddSurf("P", snum, N, [0, 0, -(w1/2+t[0]), sin(alpha), cos(alpha), 0], Surf)
            snum = AddSurf("P", snum, N, [0, 0, -(w1/2+t[0]+substhick), sin(alpha),  
                           cos(alpha), 0], Surf)
            snum = snum+10-snum%10
            snum = AddSurf("P", snum, N, [0, 0, w1/2, -sin(alpha), cos(alpha), 0], Surf)
            snum = AddSurf("P", snum, N, [0, 0, w1/2+t[0], -sin(alpha), cos(alpha), 0], Surf)
            snum = AddSurf("P", snum, N, [0, 0, w1/2+t[0]+substhick, -sin(alpha), cos(alpha),
                           0], Surf)

    elif Comp.type.strip() in ["Elliptic_guide_gravity", "Elliptic_Guide"]:  #Elliptical guide
        y1 = -100*Comp.linxw
        y2 = 100*Comp.l + 100*Comp.loutxw
        y = 100*Comp.xwidth
        z1 = -100*Comp.linyh
        z2 = 100*Comp.l + 100*Comp.loutyh
        z = 100*Comp.yheight
        if Comp.dimensionsAt == "entrance":
            posy = y1
            posz = z1
        elif Comp.dimensionsAt == "exit":
            posy = y2
            posz = z2
        elif Comp.dimensionsAt == "mid":
            posy = (y1+y2)/2
            posz = (z1+z2)/2
        else:
            raise Exception("Unknown definition of dimensionsAt for {0}\n".format(Comp.name))
        yell = calc_ellipse(y1, y2, z, posz)
        zell = calc_ellipse(z1, z2, y, posy)
        snum = AddSurf("REC", snum, N, [[(z1+z2)/2, 0, 0], [0, 1, 0], [zell["major"], 0, 0],
                                        [0, 0, zell["minor"]]], Surf)
        snum = AddSurf("REC", snum, N, [[(z1+z2)/2, 0, 0], [0, 1, 0], [zell["major"], 0, 0],
                                        [0, 0, zell["minor"]+t[0]]], Surf)
        snum = AddSurf("REC", snum, N, [[(z1+z2)/2, 0, 0], [0, 1, 0], [zell["major"], 0, 0],
                                        [0, 0, zell["minor"]+t[0]+substhick]], Surf)
        snum = snum+10-snum%10
        snum+=10
        snum = AddSurf("REC", snum, N, [[(y1+y2)/2, 0, 0], [0, 0, 1], [yell["major"], 0, 0],
                                        [0, yell["minor"], 0]], Surf)
        snum = AddSurf("REC", snum, N, [[(y1+y2)/2, 0, 0], [0, 0, 1], [yell["major"], 0, 0], 
                                        [0, yell["minor"]+t[0], 0]], Surf)
        snum = AddSurf("REC", snum, N, [[(y1+y2)/2, 0, 0], [0, 0, 1], [yell["major"], 0, 0],
                                        [0, yell["minor"]+t[0]+substhick, 0]], Surf)
        snum = snum+10-snum%10
        snum+=10
    snum = AddSurf("C/X", snum, N, [0,0,10], Surf)   #Housing
    snum = AddSurf("C/X", snum, N, [0,0,10.5], Surf)   #Housing
    snum = snum +10 - snum%10  
    snum = AddSurf("PX", snum, N, l, Surf)
    return Surf


def CreateShieldSurf(Transforms,SThick,CThick,radius,FLen,Sections):
    """
    Creates the surface cards for Shielding, using the Transformations for the axis 
    and Steel and Concrete given by SThick and Cthick.Radius array  and a final segment with 
    length FLen are used as in the mirror surface generator. Only does the geometry for now."""
 
    # Previous checks
    AxisLen=len(Transforms)
    if (len(SThick)<AxisLen) or (len(SThick)<AxisLen) or (len(radius)<AxisLen):
        print("Not enough sizes given, check variables")
        return
    Surf=[]
    if len(radius)==AxisLen-1:
        radius.append(0)
    for i,j in enumerate(Transforms):
        Section=Sections[i]
        Snum=1000*(i+1)+100
        TRN=i+1
        #UP
        Snum=AddSurf("PZ",Snum,TRN,22.5,Surf)
        Snum=AddSurf("PZ",Snum,TRN,22.5+SThick[i],Surf)
        Snum=AddSurf("PZ",Snum,TRN,22.5+SThick[i]+CThick[i],Surf)
        #Down
        Snum=Snum+10-Snum%10 

        if (Section=="D"):
            Snum=AddSurf("PZ",Snum,TRN,-12.5,Surf) # Lower edge of steel
            Snum=AddSurf("PZ",Snum,TRN,-112.5,Surf) # Ground
            Snum=AddSurf("PZ",Snum,TRN,-162.5,Surf) # Ground
            Snum=AddSurf("PZ",Snum,TRN,-212.5,Surf) # Ground
        elif (Section=="E"):
            Snum=AddSurf("PZ",Snum,TRN,-5.5,Surf) 
            Snum=AddSurf("PZ",Snum,TRN,-12.5,Surf) 
            Snum=AddSurf("PZ",Snum,TRN,-26.5,Surf) 
            Snum=AddSurf("PZ",Snum,TRN,-112.5,Surf) # Ground
        if radius[i]==0:
            #Left
            Snum=Snum+10-Snum%10 
            Snum=AddSurf("PY",Snum,TRN,-22.5,Surf)
            Snum=AddSurf("PY",Snum,TRN,-22.5-SThick[i],Surf)
            if (Section=="D"):
                Snum=AddSurf("PY",Snum,TRN,-33.5,Surf)
                Snum=AddSurf("PY",Snum,TRN,-33.5-CThick[i],Surf)
            elif (Section=="E"):
                Snum=AddSurf("PY",Snum,TRN,-50,Surf)
                Snum=AddSurf("PY",Snum,TRN,-83.5,Surf)
                Snum=AddSurf("PY",Snum,TRN,-83.5-CThick[i],Surf)

            #Right
            Snum=Snum+10-Snum%10 
            Snum=AddSurf("PY",Snum,TRN,22.5,Surf)
            Snum=AddSurf("PY",Snum,TRN,22.5+SThick[i],Surf)
            if (Section=="D"):
                Snum=AddSurf("PY",Snum,TRN,33.5,Surf)
                Snum=AddSurf("PY",Snum,TRN,33.5+CThick[i],Surf)
            elif (Section=="E"):
                Snum=AddSurf("PY",Snum,TRN,50,Surf)
                Snum=AddSurf("PY",Snum,TRN,83.5,Surf)
                Snum=AddSurf("PY",Snum,TRN,83.5+CThick[i],Surf)
        else:
            R=radius[i]
            #Left
            Snum=Snum+10-Snum%10 
            Snum=AddSurf("C/Z",Snum,TRN,[0,R,R+22.5],Surf)
            Snum=AddSurf("C/Z",Snum,TRN,[0,R,R+22.5+SThick[i]],Surf)
            if (Section=="D"):
                Snum=AddSurf("C/Z",Snum,TRN,[0,R,R+33.5],Surf)
                Snum=AddSurf("C/Z",Snum,TRN,[0,R,R+33.5+CThick[i]],Surf)
            if (Section=="E"):
                Snum=AddSurf("C/Z",Snum,TRN,[0,R,R+50],Surf)
                Snum=AddSurf("C/Z",Snum,TRN,[0,R,R+83.5],Surf)
                Snum=AddSurf("C/Z",Snum,TRN,[0,R,R+83.5+CThick[i]],Surf)
            #Right
            Snum=Snum+10-Snum%10 
            Snum=AddSurf("C/Z",Snum,TRN,[0,R,R-22.5],Surf)
            Snum=AddSurf("C/Z",Snum,TRN,[0,R,R-22.5-SThick[i]],Surf)
            if (Section=="D"):
                Snum=AddSurf("C/Z",Snum,TRN,[0,R,R-33.5],Surf)
                Snum=AddSurf("C/Z",Snum,TRN,[0,R,R-33.5-CThick[i]],Surf)
            if (Section=="E"):
                Snum=AddSurf("C/Z",Snum,TRN,[0,R,R-50],Surf)
                Snum=AddSurf("C/Z",Snum,TRN,[0,R,R-83.5],Surf)
                Snum=AddSurf("C/Z",Snum,TRN,[0,R,R-83.5-CThick[i]],Surf)
    return Surf 

def createguidecells(Comp, N):
    """
    Creates the guide cells for Component Comp. N is the corresponding number, used
    for cell numbering"""
    Coat_Mat=1
    Subs_Mat=2 # Constants for material definition. Kinda hackish but I don't think these need change often?
    Hous_Mat=3
    S_Mat=3
    C_Mat=4
    W_Mat=5
    Coat_Ro=-6.45
    Subs_Ro=-2.3 # TODO: Yeah, I just made up the numbers. Fix. 
    Hous_Ro=-7.8
    S_Ro=-7.8
    C_Ro=-2.5
    W_Ro=-2.3
    BB=BoundBox(150,250,150,150)
#    I=[B.rsplit()[0] for B in BB] #Global bounding box
#    print(I)
#    I[0]="-"+I[0]
#    I[3]="-"+I[3]
#    Sec=Section[i]
#    CType=Type[i]
    SPrefix = 100*N
    CNum = 100*N
    CList = []
    if Comp.type.strip() in ["Guide_gravity", "Guide", "Guide_gravity_polar"]:
        I0 = [SPrefix, -(SPrefix+10), SPrefix+20, SPrefix+30, -(SPrefix+40), -(SPrefix+60)] 
        I1 = [SPrefix, -(SPrefix+11), SPrefix+21, SPrefix+31, -(SPrefix+41), -(SPrefix+60)] 
        I2 = [SPrefix, -(SPrefix+12), SPrefix+22, SPrefix+32, -(SPrefix+42), -(SPrefix+60)] 
        I3 = [SPrefix, -(SPrefix+50), -(SPrefix+60)] 
        I4 = [SPrefix, SPrefix+50, -(SPrefix+51), -(SPrefix+60)] 
        E0 = []
        E1 = [SPrefix+10, -(SPrefix+20), -(SPrefix+30), SPrefix+40] 
        E2 = [SPrefix+11, -(SPrefix+21), -(SPrefix+31), SPrefix+41] 
        E3 = [SPrefix+12, -(SPrefix+22), -(SPrefix+32), SPrefix+42] 
        E4 = []
    elif Comp.type.strip() in ["Elliptic_guide", "Elliptic_guide_gravity"]:
        I0 = [SPrefix, -(SPrefix+10), -(SPrefix+30), -(SPrefix+60)] 
        I1 = [SPrefix, -(SPrefix+11), -(SPrefix+31), -(SPrefix+60)] 
        I2 = [SPrefix, -(SPrefix+12), -(SPrefix+32), -(SPrefix+60)] 
        I3 = [SPrefix, -(SPrefix+50), -(SPrefix+60)] 
        I4 = [SPrefix, SPrefix+50, -(SPrefix+51), -(SPrefix+60)] 
        E0 = []
        E1 = [SPrefix+10, SPrefix+30] 
        E2 = [SPrefix+11, SPrefix+31] 
        E3 = [SPrefix+12, SPrefix+32] 
        E4 = []

    CNum = AddCell(CNum,0,0,I0,E0,CList)  # Void inside guide
    CNum = AddCell(CNum,Coat_Mat,Coat_Ro,I1,E1,CList)  # Coating
    CNum = AddCell(CNum,Subs_Mat,Subs_Ro,I2,E2,CList)  # Substrate
    CNum = AddCell(CNum,0,0,I3,E3,CList)  # Void outside
    CNum = AddCell(CNum,Hous_Mat,Hous_Ro,I4,E4,CList)  # Housing
    CNum= 100*N + 99
    Cnum= AddCell(CNum, 0, 0, [SPrefix, -(SPrefix+60), SPrefix+51], [], CList, IMP=0) # Outside
    return CList
       

def writeTR (Comp, N):
    """Write the TR card of component Comp."""
    TR = Comp.rot[0]
    TRcard = ["TR{0} {1} {2} {3}".format(N, 100*Comp.pos.x, 100*Comp.pos.y, 100*Comp.pos.z)]
    for j in range(0,3):
        TRcard.append("      {0} {1} {2}".format(TR[j][0], TR[j][1], TR[j][2]))
    return TRcard


def createrefcards(Comp, N, RFLAG=2):
    """Write the reflection cards for Comp, number N, with RFLAGS"""
    cards = ["C =============== Section {0}================".format(N)]
    celln = 100*N
    if hasattr(Comp, "m"):
        m = [Comp.m]*4
    else:
        m = [Comp.mtop, Comp.mbottom, Comp.mleft, Comp.mright]
    if Comp.type.strip() in ["Elliptic_guide_gravity", "Elliptic_guide"]:
        srange = range(1, 5, 2)
    elif Comp.type.strip() in ["Guide_gravity", "Guide", "Guide_gravity_polar"]:
        srange = range(1, 5)
    for i in srange:
        cardn = 100*N + i
        surfn = 100*N + 10*i
        mn = int(10*m[i-1])  # Notice that m should not have more than 1 decimal
        cards.append("REFLE{0} {1} {2} -{3}".format(cardn, surfn, mn, celln))
        cards.append("RFLAG{0} {1}".format(cardn, RFLAG))
    return cards
    
