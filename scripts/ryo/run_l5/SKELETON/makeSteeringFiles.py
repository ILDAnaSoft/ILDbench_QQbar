#!/usr/bin/python
import os,string,shutil,sys,re,utils,conf


def writeXMLs(pname):
	# make a list that has only accepted file names.
	filteredlist = []
	utils.makeFilteredList(filteredlist,conf.LISTDIR,pname)
	nxmlfile = 0
	while (1):
		if (len(filteredlist)<=0):
			break
		else:
			nxmlfile += 1
			ofname = conf.XMLDIR + "/"+conf.OUTPUT_PREFIX+pname + "_" + str(nxmlfile) + ".xml"
			foutxml = open(ofname,'w')
			template = open(conf.TEMPLATEDIR+"/"+conf.TEMPLATEFILE,'r')
			for t in template:
				if re.search("__INPUTFILES__",t):
	
					for i in range(conf.nfilesInOneshot):
						if (len(filteredlist)==0):
							break
						if (conf.nfilesInOneshot>1 and i==0):
							fpath=conf.SPACE+filteredlist.pop()+"\n"
						elif (conf.nfilesInOneshot==1):
							fpath=conf.SPACE+filteredlist.pop()
						elif (i==conf.nfilesInOneshot-1) :
							fpath+=conf.SPACE+filteredlist.pop()
						else :
							if (len(filteredlist)==1):#last case 
								fpath+=conf.SPACE+filteredlist.pop()
							else:
								fpath+=conf.SPACE+filteredlist.pop()+"\n"
					foutxml.write(re.sub("__INPUTFILES__",fpath,t))
				elif re.search("__OUTPUTROOTDIR__/__OUTPUTROOT__",t):		
					tmp = re.sub("__OUTPUTROOTDIR__",conf.OUTDIR_ROOT,t)
					tmp2 = conf.SPACE + re.sub("__OUTPUTROOT__",conf.OUTPUT_PREFIX+pname + "_" + str(nxmlfile) + ".root",tmp)
					foutxml.write(tmp2)
				elif re.search("__GEARFILE__",t):		
					tmp = conf.SPACE + re.sub("__GEARFILE__",conf.GEARFILE,t)
					foutxml.write(tmp)
				elif re.search("__INPUTPREFIX__",t):		
					tmp = conf.SPACE + re.sub("__INPUTPREFIX__",conf.INPUT_PREFIX,t)
					foutxml.write(tmp)
				else:
					foutxml.write(t)


# main function

# prepare directories if necessary
utils.makeDirectory(conf.TEMPLATEDIR)
utils.makeDirectory(conf.XMLDIR)
utils.makeDirectory(conf.OUTDIR_ROOT)

# copy xml to this directory
shutil.copy(conf.TEMPLATEXMLORIGDIR+"/"+conf.TEMPLATEFILE,conf.TEMPLATEDIR+"/"+conf.TEMPLATEFILE)

for proc in conf.PROCESSES:
	writeXMLs(proc);

