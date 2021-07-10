import argparse
import fileinput
import os
import subprocess
import sys

parser = argparse.ArgumentParser(description='Compile all GLSL shaders')
parser.add_argument('--glslang', type=str, help='path to glslangvalidator executable')
parser.add_argument('--g', action='store_true', help='compile with debug symbols')
args = parser.parse_args()

def findGlslang():
    def isExe(path):
        return os.path.isfile(path) and os.access(path, os.X_OK)

    if args.glslang != None and isExe(args.glslang):
        return args.glslang

    exe_name = "glslangvalidator"
    if os.name == "nt":
        exe_name += ".exe"

    for exe_dir in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(exe_dir, exe_name)
        if isExe(full_path):
            return full_path

    sys.exit("Could not find glslangValidator executable on PATH")

def findSpirvOpt():
    def isExe(path):
        return os.path.isfile(path) and os.access(path, os.X_OK)

    if args.glslang != None and isExe(args.glslang):
        return args.glslang

    exe_name = "spirv-opt"
    if os.name == "nt":
        exe_name += ".exe"

    for exe_dir in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(exe_dir, exe_name)
        if isExe(full_path):
            return full_path

    print("Could not find spirv-opt on PATH")
    return ""

glslang_path = findGlslang()
spirvopt_path = findSpirvOpt()
subdirs = ['double', 'float']
for dp in subdirs:
    os.chdir(dp)
    for root, dirs, files in os.walk(os.getcwd()):
        for file in files:
            if file.endswith(".vert") or file.endswith(".frag") or file.endswith(".comp") or \
            file.endswith(".geom") or file.endswith(".tesc") or file.endswith(".tese") or \
            file.endswith(".rgen") or file.endswith(".rchit") or file.endswith(".rmiss") or \
            file.endswith(".glsl"):

                file_short = file[:file.index(".glsl")]
                input_file = os.path.join(root, file)
                output_file = file_short + ".spv"

                add_params = ""
                if args.g:
                    add_params = "-g"

                if file.endswith(".rgen") or file.endswith(".rchit") or file.endswith(".rmiss"):
                    add_params = add_params + " --target-env vulkan1.2"

                res = subprocess.call("%s -V %s --target-env vulkan1.2 -o %s %s" % (glslang_path, input_file, output_file, add_params), shell=True)

                if res != 0:
                    sys.exit()
        
        if(spirvopt_path):
            print("Optimizing shaders...")
            for file in files:
                if(file.endswith(".spv")):
                    res_opt = subprocess.call("%s -O %s -o %s" % (spirvopt_path, file, file), shell=True)
                    if res_opt != 0:
                        sys.exit()
    os.chdir('..')