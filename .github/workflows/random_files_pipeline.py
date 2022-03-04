import sandy
import os, glob, sys, shutil

library = sys.argv[1].lower().replace(".", "").replace("/", "").replace("-", "_")
nuclide = int(sys.argv[2])
temperature = float(sys.argv[3])
fmt = sys.argv[4]
samples = int(sys.argv[5])

tape = sandy.get_endf6_file(library, 'xs', nuclide)
name = "file.endf6"
folder = "random_files"
tape.to_file(name)
cline = f"{name} --samples {samples} --temperatures {temperature} --mf 33 -D {folder}"
if fmt == "ACE":
    cline += " --acer"
out = sandy.sampling(cline.split())

pendf = glob.glob("*.pendf")[0]
basename = pendf.split(".")[0]
os.remove(pendf)
os.remove(name)
if fmt == "ACE":
    for f in glob.glob(f"**/{name}*", recursive=True):
        os.remove(f)
else:
    for f in glob.glob(f"**/{name}*", recursive=True):
        base, smp = f.split("-")
        path, base = os.path.split(base)
        new_name = "-".join([basename, smp])
        new_file = os.path.join(path, new_name)
        shutil.move(f, new_file)
