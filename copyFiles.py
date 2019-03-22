import os
import shutil

BACKS = 'src/assets/backs'

hash = os.listdir('dist/static/img')[0].split('.')[1]
for file in os.listdir(BACKS):
    newFileName = file.split('.')[0] + '.' + hash + '.png'
    shutil.copyfile(BACKS + '/' + file, 'dist/static/css/static/img/' + newFileName)
