import unittest
import sys
import os
import logging
import filecmp
import shutil

import structureWrapperNew2
import structure
import structure_func
import ruffus
import itertools


from itertools import combinations


class filterTest(unittest.TestCase):

    def test_calc_run_1(self):
        
        #------------------------------------------------directory names-----------------------------------------------#
        filenames = ['mainparams','extracols','infile','label','locdata','missing','markernames','recessivealleles','mapdistances','numinds','numloci','onerowperind','phenotype','ploidy','popdata']
    

        #------------holds all paths for input directory--------------#
        mainpath = []

        #------------holds all paths for mainparams inside individual directory-------------------#
        filepath = []

        #-----------holds all paths for output directory-------------#
        outpath = []

        #-----------path to running directory---------------#
        outdir = "/home/staff/asukeshkall/Downloads/sk2-2018-test/"

        #-----------path to mainparas ref----------------# 
        outmain = "/home/staff/asukeshkall/Downloads/sk2-2018-test/mainparams"

        #-----------path to store the correct output-----------#
        respath= "/home/staff/asukeshkall/Downloads/sk2-2018-test/Testing/mainparams"

        for i in range(15):
            os.makedirs(os.path.join('Testing', filenames[i]))

        for i in filenames:
            path = "/home/staff/asukeshkall/Downloads/sk2-2018-test/files/"+i
            mainpath.append(path)
 
        for i in filenames:
            main = "/home/staff/asukeshkall/Downloads/sk2-2018-test/files/"+i+"/mainparams"
            filepath.append(main)

        for i in filenames:
            path = "/home/staff/asukeshkall/Downloads/sk2-2018-test/Testing/"+i
            outpath.append(path)

        #-----------function for structure call with incorrect input--------------#
        def structureCall(index):
            shutil.move(filepath[i],outdir)
            structure_func.run()
            shutil.move(outmain,mainpath[i])

            self.assertTrue(os.path.isfile('structureOutput.log'))
            shutil.move('structureOutput.log', outpath[i])

        #---------function for structure call with correct input-------------#       
        for index,path in enumerate(filepath):
            shutil.move(filepath[0],outdir)
            structure_func.run()
            shutil.move(outmain,mainpath[0])

            self.assertTrue(os.path.isfile('blandingsout_f'))
            self.assertTrue(os.path.isfile('structureOutput.log'))

        shutil.move('blandingsout_f', respath)
        shutil.move('structureOutput.log', respath)

  

        #--------call structure
        for i in range(1,len(filenames)):    
            structureCall(i)
 

if __name__ == "__main__":
    unittest.main()
