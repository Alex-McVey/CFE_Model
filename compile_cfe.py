import os
import glob
import shutil
import PyInstaller.__main__
import PyInstaller.config


print('********************************************************************')
print('               COMPILING CFE Tool for distribution')
print('********************************************************************\n\n')

# Clean up build and dist files if exist:
print('************************************ Cleaning up any previous builds\n\n')
if os.path.isdir("dist"):
    shutil.rmtree('dist')
if os.path.isdir("build"):
    shutil.rmtree('build')
specfiles = glob.glob("*.spec")
for f in specfiles:
    os.remove(f)

print('******************************** Calling PyInstaller to package .exe\n\n')
script = os.path.join(os.getcwd(),"cfe_gui.py")

# call pyinstaller and package:
# PyInstaller.__main__.run([
#     script,
#     '-nACES',
#     f'-i{os.path.join(os.getcwd(), "images", "aces1.ico")}',
#     '--onedir',
#     '--windowed'
# ])
PyInstaller.__main__.run([
    script,
    '-nCFE_Model',
    '--onedir',
    '--windowed'
])

print('\n\n************************* Copying over missing files and directories\n\n')

# --> images (copy)
print('copying images directory')
imgdr1 = os.path.join(os.getcwd(), 'images')
imgdr2 = os.path.join('.','dist','CFE_Model','images')
shutil.copytree(imgdr1,imgdr2)
print('copying data directory')
datadr1 = os.path.join(os.getcwd(), 'data')
datadr2 = os.path.join('.','dist','CFE_Model','_internal','data')
shutil.copytree(datadr1,datadr2)
print('copying missing imports from rasterio')
rastdr1 = "C://Users/AlexMcVey/anaconda3/envs/cfe_tool/Lib/site-packages/rasterio"
rastdr2 = os.path.join('.','dist','CFE_Model', '_internal', 'rasterio')
shutil.copy2(os.path.join(rastdr1, 'sample.py'),os.path.join(rastdr2, 'sample.py'))
shutil.copy2(os.path.join(rastdr1, 'vrt.py'),os.path.join(rastdr2, 'vrt.py'))
shutil.copy2(os.path.join(rastdr1, '_features.cp311-win_amd64.pyd'),os.path.join(rastdr2, '_features.cp311-win_amd64.pyd'))
shutil.copy2(os.path.join(rastdr1, '_warp.cp311-win_amd64.pyd'),os.path.join(rastdr2, '_warp.cp311-win_amd64.pyd'))

# --> copy ui file
shutil.copy2(os.path.join('.','cfe_model.ui'),os.path.join('.','dist','CFE_Model', '_internal','cfe_model.ui'))
shutil.copy2(os.path.join('.','model_assumptions.ui'),os.path.join('.','dist','CFE_Model', '_internal','model_assumptions.ui'))



