""" MAGiC guide to test Blinggen"""

import mcparser
import mcnpgenerator 
import math
import geocell

headerline = "MAGiC beamline model created using blinggen.\n"
varlib = mcparser.getvariables("MAGiC_var.c")
comps = mcparser.getcomps("MAGiC_therm.instr", varlib)
mcparser.makeinstr(comps)

wanted_comps = ["Guide_gravity_polar", "Elliptic_guide_gravity", "Guide_gravity"]

surfs = []
cells = []
TRs = []
reffs = []
marray = []
writtencomp = []
lastx = [None for i in comps]
l = [200 for i in range(100)]
l[13] = 100  # Polarizer
for i, comp in enumerate(comps):
    if  comp.type.strip() in wanted_comps:
        writtencomp.append(i)
        compsurfs = mcnpgenerator.createguidesurfs(comp, i, length=l[i])
        nsegs = math.ceil(comp.l*100 / l[i]) 
        compcells = mcnpgenerator.createguidecells(comp, i, nsegs)
        compTR = mcnpgenerator.createTR(comp, i)
        compreffs = mcnpgenerator.createreflecards(comp, i, length=l[i])
        lastx[i] = 100*i + 59 + 2*nsegs  # Last X plane
        for param in ["m", "mtop", "mbottom", "mleft", "mright"]:
            if hasattr(comp, param):
                marray.append(vars(comp)[param])
# Safecheck againts null instruments
        for line in compsurfs:
            surfs.append(line+"\n")
        for line in compcells:
            cells.append(line+"\n")
        for line in compTR:
            TRs.append(line+"\n")
        for line in compreffs:
            reffs.append(line+"\n")
mlist = list(set(marray))  # prune duplicates
for i1, i2 in zip(writtencomp[:-1], writtencomp[1:]):
##  ADD FILLER CELLS
    filler = geocell.cell([lastx[i1], -(100*i1+50), -100*i2], 10*i1) 
    cells.append(geocell.mcnpcard(filler)+"\n")
    filler = geocell.cell([lastx[i1], 100*i1+50, -100*i2], 10*i1+1) 
    cells.append(geocell.mcnpcard(filler)+"\n")
    
with open("MAGiC_ent.i", "w") as entfile:
    entfile.write(headerline)
    entfile.writelines(cells)
    entfile.write("\n")
    entfile.writelines(surfs)
    entfile.write("\n")
    entfile.writelines(TRs)
    entfile.write("READ FILE=MATERIAL.i\n")
    for m in mlist:
        entfile.write(mcnpgenerator.createreffcard(m)+"\n")
    entfile.writelines(reffs)

