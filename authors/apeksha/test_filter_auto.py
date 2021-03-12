import unittest
import sys
import os
import logging
import filecmp
import shutil

import vcf_filter
import vcftools
import ruffus
import itertools
from itertools import combinations


class filterTest(unittest.TestCase):
 
#=================================function to make dir and sub -dir #================================#

    def test_filter_run_1(self):
        for i in range(1,61):
           os.makedirs(os.path.join('Run', 'test' + str(i))) 

#============================# Create path names to move the generated o/p files #===================#
        newpath = []
        filenamee = 'test'
        counter1 = 1
        counter2 = 2
        list1 = []
        list2 = range(61)    #random number, given len needed
        for x in list2:
            counter1 = str(counter1)
            full_name = (filenamee+counter1)
            list1.append(full_name)
            counter1 = counter2
            counter2+=1

        for x in list1:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            newpath.append(y)
        


#================================================# Generate all combinations  #=======================================================================#

 
        comb = []

        a=[['merged_chr1_1000.vcf.gz'],['--out-prefix','out'],['--out-prefix','res'],['--out-format','removed_sites'],['--out-format','kept_sites'],['--out-format','vcf'],['--out-format','bcf'],['--filter-min-alleles', '2'],['--filter-max-alleles', '4']]

        for i in range(len(a)):

            b = list(combinations(a, i+1))
            c = map(list, b)
            comb.append(c)

        whole_list = []

        for x in comb:
            for y in x:
                b = list(itertools.chain.from_iterable(y))
                whole_list.append(b)
                #print b
        m = whole_list

#=============================================# Filter level 1 function #========================================================================#

        filList1 = []


        def filFunc(name, lol, listname):
            for x in lol:
                for y in x:
                    if y== name:
                       listname.append(x) 

#============================= input =====================#

        filFunc('merged_chr1_1000.vcf.gz', m, filList1)
        fnameFilter = filList1
        for t in fnameFilter:
            print t
            #print "--------------"

#======================# Filter level 2 - to filter the redundant options #===================================#

        app = []
        for m in fnameFilter:  
            for x,y in itertools.combinations(m, 2):  
                if x == y:
                   app.append(m) 
        
        unmatched = [d for d in fnameFilter if d not in app]

        print unmatched
            



#============================================================================================================================#

        file_names = ['out.removed.sites', 'out.kept.sites','out.recode.vcf','out.recode.bcf','res.removed.sites','res.kept.sites','res.recode.vcf','res.recode.bcf']

        list_names = ['fnameFilter', 'outKeptSitesFilter', 'resRemSitesFilter', 'resKeptSitesFilter', 'outVcfFilter', 'outBcfFilter', 'resVcfFilter', 'resBcfFilter']
        
#============================# test function #=======================================================#

        def testFunc(lname,fn,np):
            vcf_filter.run(lname)
            self.assertTrue(os.path.isfile(file_names[fn]))
            self.assertTrue(os.path.isfile(file_names[fn]+'.log'))     
            shutil.move(file_names[fn], newpath[np])
            shutil.move(file_names[fn]+'.log', newpath[np])

#=================================================================================# test function calls #===================================================================================#

        # i/p 
        testFunc(unmatched[0], 0, 0) 
        
        # i/p ; --out- prefix - out
        testFunc(unmatched[1], 0, 1)

        # i/p ; --out- prefix - res
        testFunc(unmatched[2], 4, 2)
        
        # i/p ; --out- format - rem sites
        testFunc(unmatched[3], 0, 3)

        # i/p ; --out- format - kept sites
        testFunc(unmatched[4], 1, 4)

        # i/p ; --out- format - vcf
        testFunc(unmatched[5], 2, 5)

        # i/p ; --out- format - bcf
        testFunc(unmatched[6], 3, 6)

        # i/p ; min alleles
        testFunc(unmatched[7], 0, 7)

        # i/p ; max alleles
        testFunc(unmatched[8], 0, 8)

        # i/p ; --out- prefix - out; --out- format - rem sites
        testFunc(unmatched[9], 0, 9)

        # i/p ; --out- prefix - out; --out- format - kept sites
        testFunc(unmatched[10], 1, 10)

        # i/p ; --out- prefix - out; --out- format - vcf
        testFunc(unmatched[11], 2, 11)

        # i/p ; --out- prefix - out; --out- format - bcf
        testFunc(unmatched[12], 3, 12)

        # i/p ; --out- prefix - out; min alleles 
        testFunc(unmatched[13], 0, 13)

        # i/p ; --out- prefix - out; max alleles
        testFunc(unmatched[14], 0, 14)

        # i/p ; --out- prefix - res; --out- format - rem sites
        testFunc(unmatched[15], 4, 15)

        # i/p ; --out- prefix - res; --out- format - kept sites
        testFunc(unmatched[16], 5, 16)

        # i/p ; --out- prefix - res; --out- format - vcf
        testFunc(unmatched[17], 6, 17)

        # i/p ; --out- prefix - res; --out- format - bcf
        testFunc(unmatched[18], 7, 18)

        # i/p ; --out- prefix - res; min alleles 
        testFunc(unmatched[19], 4, 19)

        # i/p ; --out- prefix - res; max alleles
        testFunc(unmatched[20], 4, 20)

      # i/p ; -out-format- rem sites ; min alleles
        testFunc(unmatched[21], 0, 21)

      # i/p ; -out-format- rem sites ; max alleles
        testFunc(unmatched[22], 0, 22)

      # i/p ; -out-format- kept sites ; min alleles
        testFunc(unmatched[23], 1, 23)

      # i/p ; -out-format- kept sites ; max alleles
        testFunc(unmatched[24], 1, 24)

      # i/p ;-out-format- vcf sites ; min alleles
        testFunc(unmatched[25], 2, 25)

      # i/p ;-out-format- vcf sites ; max alleles
        testFunc(unmatched[26], 2, 26)

      # i/p ;-out-format- bcf sites ; min alleles
        testFunc(unmatched[27], 3, 27)

      # i/p ;-out-format- bcf sites ; max alleles
        testFunc(unmatched[28], 3, 28)

      # i/p ;min alleles ; max alleles
        testFunc(unmatched[29], 0, 29)

      # i/p ;out - prefix - out ; -out-format - rem sites ; min alleles 
        testFunc(unmatched[30], 0, 30)

      # i/p ;out - prefix - out ; -out-format - rem sites ; max alleles 
        testFunc(unmatched[31], 0, 31)

      # i/p ;out - prefix - out ; -out-format - kept sites ; min alleles 
        testFunc(unmatched[32], 1, 32)

      # i/p ;out - prefix - out ; -out-format - kept sites ; max alleles 
        testFunc(unmatched[33], 1, 33)

      # i/p ;out - prefix - out ; -out-format - vcf ; min alleles 
        testFunc(unmatched[34], 2, 34)

      # i/p ;out - prefix - out ; -out-format - vcf ; max alleles 
        testFunc(unmatched[35], 2, 35)

      # i/p ;out - prefix - out ; -out-format - bcf ; min alleles 
        testFunc(unmatched[36], 3, 36)

      # i/p ;out - prefix - out ; -out-format - bcf ; max alleles 
        testFunc(unmatched[37], 3, 37)

      # i/p ;out - prefix - out ; min alleles ; max alleles
        testFunc(unmatched[38], 0, 38)

      # i/p ;out - prefix - res ; -out-format - rem sites ; min alleles
        testFunc(unmatched[39], 4, 39)

      # i/p ;out - prefix - res ; -out-format - rem sites ; max alleles
        testFunc(unmatched[40], 4, 40)

      # i/p ;out - prefix - res ; -out-format - kept sites ; min alleles
        testFunc(unmatched[41], 5, 41)

      # i/p ;out - prefix - res ; -out-format - kept sites ; max alleles
        testFunc(unmatched[42], 5, 42)

      #i/p ; out - prefix - res ; -out-format - vcf ; min alleles
        testFunc(unmatched[43], 6, 43)

      # i/p ;out - prefix - res ; -out-format - vcf ; max alleles
        testFunc(unmatched[44], 6, 44)

      # i/p ;out - prefix - res ; -out-format - bcf ; min alleles
        testFunc(unmatched[45], 7, 45)

      # i/p ;out - prefix - res ; -out-format - bcf ; max alleles
        testFunc(unmatched[46], 7, 46)

      # i/p ;out - prefix - res ; min alleles ; max alleles
        testFunc(unmatched[47], 4, 47)

      # i/p ;-out-format - rem sites ; min alleles ; max alleles
        testFunc(unmatched[48], 0, 48)

      # i/p ;-out-format - kept sites ; min alleles ; max alleles
        testFunc(unmatched[49], 1, 49)

      # i/p ;-out-format - vcf ; min alleles ; max alleles
        testFunc(unmatched[50], 2, 50)

      # i/p ;-out-format - bcf ; min alleles ; max alleles
        testFunc(unmatched[51], 3, 51)

      # i/p ;out - prefix - out ; -out-format - rem sites ; min alleles ; max alleles
        testFunc(unmatched[52], 0, 52)

      # i/p ;out - prefix - out ; -out-format - kept sites ; min alleles ; max alleles
        testFunc(unmatched[53], 1, 53)

      # i/p ;out - prefix - out ; -out-format - vcf ; min alleles ; max alleles
        testFunc(unmatched[54], 2, 54)

      # i/p ;out - prefix - out ; -out-format - bcf ; min alleles ; max alleles
        testFunc(unmatched[55], 3, 55)

      # i/p ;out - prefix - res ; -out-format - rem sites ; min alleles ; max alleles
        testFunc(unmatched[56], 4, 56)

      # i/p ;out - prefix - res ; -out-format - kept sites ; min alleles ; max alleles
        testFunc(unmatched[57], 5, 57)

      # i/p ;out - prefix - res ; -out-format - vcf ; min alleles ; max alleles
        testFunc(unmatched[58], 6, 58)

      # i/p ;out - prefix - res ; -out-format - bcf ; min alleles ; max alleles
        testFunc(unmatched[59], 7, 59)

   

if __name__ == '__main__':
    unittest.main()
