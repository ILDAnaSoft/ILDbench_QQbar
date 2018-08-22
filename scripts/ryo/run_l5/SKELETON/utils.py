#!/usr/bin/python
import os,conf
def makeDirectory(dirname):
	if not os.path.exists(dirname):
		os.mkdir(dirname)

def appendFilteredList(filteredlist,filepath,pname):
        filename = filepath.split('/')[-1]
	suffix = filename.split('.')[-1]
	if (filename.count(pname) and suffix==conf.FILESUFFIX):
		filteredlist.append(filepath)

def makeFilteredList(filteredlist,listdir,pname):
	items = os.listdir(listdir)
	for item in items:
		newpath = listdir+"/"+item
		if (os.path.isdir(newpath)) :
			makeFilteredList(filteredlist,newpath,pname)
		if (os.path.isfile(newpath)) :
			appendFilteredList(filteredlist,newpath,pname)
