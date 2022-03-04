import sandy
import os, glob, sys, shutil

library = sys.argv[1].lower().replace(".", "").replace("/", "").replace("-", "_")
nuclide = int(sys.argv[2])
temperature = float(sys.argv[3])
fmt = sys.argv[4]
samples = int(sys.argv[5])

tape = sandy.get_endf6_file(library, 'xs', nuclide)
name = "file.endf6"
tape.to_file(name)
cline = f"{name} --samples {samples} --temperatures {temperature} --mf 33"
if fmt == "ACE":
    cline += " --acer"
out = sandy.sampling(cline.split())

pendf = glob.glob("*.pendf")[0]
basename = pendf.split(".")[0]
os.remove(pendf)
os.remove(name)
if fmt == "ACE":
    for f in glob.glob(f"{name}*"):
        os.remove(f)
else:
    for f in glob.glob(f"{name}*"):
        base, smp = f.split("-")
        shutil.move(f, "-".join([basename, smp]))