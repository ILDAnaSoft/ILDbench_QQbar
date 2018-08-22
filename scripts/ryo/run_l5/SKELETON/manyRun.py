#!/usr/bin/python
import os,utils,conf

utils.makeDirectory(conf.LOGDIR)
xmls = os.listdir(conf.XMLDIR)

minIndex=0
maxIndex=len(xmls)
#minIndex=100
#maxIndex=200

for xml in xmls:
	suffix = '.' + xml.split('.')[-1] # suffix = .slcio
	basename = xml.split(suffix)[0]
	index = basename.split('_')[-1]
	if int(index) <= maxIndex and int(index) > minIndex :
		cmd = "bsub -q s " + "-o "+conf.LOGDIR+"/" + basename + ".log " + '\"Marlin ' + conf.XMLDIR + '/' + xml.strip() + '\"'
		os.system(cmd)
		#print cmd

