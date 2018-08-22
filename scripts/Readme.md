## scripts directory


Put your Marlin steering files here and explain to which processor they refer (if the file name is not obvious) !

exampleX_l/s5.xml : test steering files.  
lcfiplus : test weight files of LCFIPlus. This should be replaced to a central one at some point.  
ryo : test directory mainly for debugging purpose.  
      * init_ilcsoft_v02-00-01_default.sh for setting environment fro KEKCC.
      * run_l5 : 
	* SKELETON : python scripts to be copied into a working directory e.g. test. These scripts will help to run many samples.
	* XML_TMPLATES : steering files to be used in python scripts above.
	* Example usage :   
           * mkdir test
           * cd test
           * cp ../SKELTON/*.py .
	   * EDIT conf.py so that it can read correct input directory where sample files exist, and specify template steering file etc. 
	   * python makeSteeringFiles.py  (This creates steering files in a directory, generatedXMLs, which automatically created.)
	   * python manyRun.py (Submit jobs for the steering files in generatedXMLs.)
 	

