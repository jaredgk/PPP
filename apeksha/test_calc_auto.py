import unittest
import sys
import os
import logging
import filecmp
import shutil

import vcf_calc
import vcf_filter
import vcftools
import ruffus
import itertools
from itertools import combinations


class calcTest(unittest.TestCase):
 
#=================================function to make dir and sub -dir #================================#

    def test_calc_run_1(self):
        for i in range(1,80):
           os.makedirs(os.path.join('Run', 'test' + str(i))) 

#============================# Create path names to move the generated o/p files #===================#
        defaultPath = []
        winWeirPath = []
        weirPath = []
        tajimaPath = []
        piPath = []
        freqPath = []
        hetPath = []
        outWinWeirPath = []
        outWeirPath = []
        outTajimaPath = []
        outPiPath = []
        outFreqPath = []
        outHetPath = []
        resWinWeirPath = []
        resWeirPath = []
        resTajimaPath = []
        resPiPath = []
        resFreqPath = []
        resHetPath = []
        filename = 'test'
        counter1 = 1
        counter2 = 2
        list1 = []
        list2 = range(80)    #random number, given len needed
        for x in list2:
            counter1 = str(counter1)
            full_name = (filename+counter1)
            list1.append(full_name)
            counter1 = counter2
            counter2+=1

        for x in list1[:4]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            defaultPath.append(y)

        for x in list1[4:8]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            winWeirPath.append(y)

        for x in list1[8:12]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            weirPath.append(y)

        for x in list1[12:16]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            tajimaPath.append(y)

        for x in list1[16:20]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            piPath.append(y)

        for x in list1[20:24]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            freqPath.append(y)

        for x in list1[24:28]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            hetPath.append(y)

        for x in list1[28:32]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            outWinWeirPath.append(y)

        for x in list1[32:36]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            outWeirPath.append(y)

        for x in list1[36:40]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            outTajimaPath.append(y)

        for x in list1[40:44]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            outPiPath.append(y)

        for x in list1[44:48]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            outFreqPath.append(y)

        for x in list1[48:52]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            outHetPath.append(y)

        for x in list1[52:56]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            resWinWeirPath.append(y)

        for x in list1[56:60]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            resWeirPath.append(y)

        for x in list1[60:64]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            resTajimaPath.append(y)

        for x in list1[64:68]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            resPiPath.append(y)

        for x in list1[68:72]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            resFreqPath.append(y)

        for x in list1[72:76]:
            y = "/home/staff/asukeshkall/Downloads/Test/Run/"+ x
            resHetPath.append(y)




#================================================# Generate all combinations  #=======================================================================#

 
        comb = []

        a=[['merged_chr1_1000.vcf.gz'],['--out-prefix','out'],['--out-prefix','res'],['--pop-file', 'Paniscus.txt'],['--pop-file', 'Troglodytes.txt'],['--calc-statistic', 'weir-fst'],['--calc-statistic', 'windowed-weir-fst'],['--calc-statistic', 'TajimaD'],['--calc-statistic', 'pi'],['--calc-statistic', 'freq'],['--calc-statistic', 'het'], ['--statistic-window-size', '2'],['--statistic-window-step', '2']]

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
        filList2 = []
        filList3 = []


        def filFunc(name, lol, listname):
            for x in lol:
                for y in x:
                    if y== name:
                       listname.append(x) 

#============================= input =====================#

        filFunc('merged_chr1_1000.vcf.gz', m, filList1)
        fnameFilter = filList1
        #for t in fnameFilter:
           # print t
            #print "--------------"

        filFunc('Paniscus.txt', fnameFilter, filList2)
        panFilter = filList2
        #for t in panFilter:
            #print t
            #print "--------------"

        filFunc('Troglodytes.txt', panFilter, filList3)
        trogFilter = filList3
        #for t in trogFilter:
            #print t
            #print "--------------"

#======================# Filter level 2 - to filter the redundant options #===================================#

        #app = []
        #for m in trogFilter:  
            #for x,y in itertools.combinations(m, 2):  
                #if x == y:
                   #app.append(m) 

        #for x in trogFilter:
            #for y in x:
            #if x.count(--d)>2:
        unmatched = []

        for e in trogFilter:
            if e.count('--calc-statistic')>1 or e.count('--out-prefix')>1:
               unmatched.append(e)
        #print unmatched
               

        unmatched1 = [d for d in trogFilter if d not in unmatched]   
        #for t in unmatched1:
           # print t 
        
        #for t in app:
            #print t
        
        #unmatched = [d for d in trogFilter if d not in app]

        #for t in unmatched:
            #print t
            



#============================================================================================================================#

        file_names = ['out.windowed.weir.fst', 'out.weir.fst','out.recode.vcf','out.recode.bcf','res.removed.sites','res.kept.sites','res.recode.vcf','res.recode.bcf']

        #list_names = ['fnameFilter', 'outKeptSitesFilter', 'resRemSitesFilter', 'resKeptSitesFilter', 'outVcfFilter', 'outBcfFilter', 'resVcfFilter', 'resBcfFilter']
        
#============================# test function #=======================================================#   
        listnew = []
        defList = []
        winWeirList = []
        weirList = []
        tajimaList = []
        piList = []
        freqList = []
        hetList = []
        outWinWeirList = []
        outWeirList = []
        outTajimaList = []
        outPiList = []
        outFreqList = []
        outHetList = []
        resWinWeirList = []
        resWeirList = []
        resTajimaList = []
        resPiList = []
        resFreqList = []
        resHetList = []
        arr=[]

        def defaultFunc(name,path,pathname):
            self.assertTrue(os.path.isfile(name+'.windowed.weir.fst'))
            self.assertTrue(os.path.isfile(name+'.windowed.weir.fst.log'))
            shutil.move(name+'.windowed.weir.fst', pathname[path])
            shutil.move(name+'.windowed.weir.fst.log', pathname[path])

        def winWeirFunc(name,path,pathname):
            self.assertTrue(os.path.isfile(name+'.windowed.weir.fst'))
            self.assertTrue(os.path.isfile(name+'.windowed.weir.fst.log'))
            shutil.move(name+'.windowed.weir.fst', pathname[path])
            shutil.move(name+'.windowed.weir.fst.log', pathname[path])

        def weirFunc(name,path,pathname):
            self.assertTrue(os.path.isfile(name+'.weir.fst'))
            self.assertTrue(os.path.isfile(name+'.weir.fst.log'))
            shutil.move(name+'.weir.fst', pathname[path])
            shutil.move(name+'.weir.fst.log', pathname[path])

        def tajimaFunc(name,path,pathname):
            self.assertTrue(os.path.isfile(name+'.Tajima.D'))
            self.assertTrue(os.path.isfile(name+'.Tajima.D.log'))
            shutil.move(name+'.Tajima.D', pathname[path])
            shutil.move(name+'.Tajima.D.log', pathname[path])

        def piFunc(name,path,pathname):
            self.assertTrue(os.path.isfile(name+'.windowed.pi'))
            self.assertTrue(os.path.isfile(name+'.windowed.pi.log'))
            shutil.move(name+'.windowed.pi', pathname[path])
            shutil.move(name+'.windowed.pi.log', pathname[path])

        def freqFunc(name,path,pathname):
            self.assertTrue(os.path.isfile(name+'.frq'))
            self.assertTrue(os.path.isfile(name+'.frq.log'))
            shutil.move(name+'.frq', pathname[path])
            shutil.move(name+'.frq.log', pathname[path])

        def hetFunc(name,path,pathname):
            self.assertTrue(os.path.isfile(name+'.het'))
            self.assertTrue(os.path.isfile(name+'.het.log'))
            shutil.move(name+'.het', pathname[path])
            shutil.move(name+'.het.log', pathname[path])


#================================================================================================#


        
        def defaultFuncCall(index,):
               vcf_calc.run(defList[index])
               defaultFunc('out',index,defaultPath)

        def winWeirFuncCall(index):
               vcf_calc.run(winWeirList[index])
               winWeirFunc('out',index,winWeirPath)

               vcf_calc.run(outWinWeirList[index])
               winWeirFunc('out',index,outWinWeirPath)

               vcf_calc.run(resWinWeirList[index])
               winWeirFunc('res',index,resWinWeirPath)

        def weirFuncCall(index):

               vcf_calc.run(weirList[index])
               weirFunc('out',index,weirPath)

               vcf_calc.run(outWeirList[index])
               weirFunc('out',index,outWeirPath)

               vcf_calc.run(resWeirList[index])
               weirFunc('res',index,resWeirPath)

        def tajimaFuncCall(index):

               vcf_calc.run(tajimaList[index])
               tajimaFunc('out',index,tajimaPath)

               vcf_calc.run(outTajimaList[index])
               tajimaFunc('out',index,outTajimaPath)

               vcf_calc.run(resTajimaList[index])
               tajimaFunc('res',index,resTajimaPath)

        def piFuncCall(index):

               vcf_calc.run(piList[index])
               piFunc('out',index,piPath)

               vcf_calc.run(outPiList[index])
               piFunc('out',index,outPiPath)

               vcf_calc.run(resPiList[index])
               piFunc('res',index,resPiPath)

        def freqFuncCall(index):

               vcf_calc.run(freqList[index])
               freqFunc('out',index,freqPath)

               vcf_calc.run(outFreqList[index])
               freqFunc('out',index,outFreqPath)

               vcf_calc.run(resFreqList[index])
               freqFunc('res',index,resFreqPath)

        def hetFuncCall(index):

               vcf_calc.run(hetList[index])
               hetFunc('out',index,hetPath)

               vcf_calc.run(outHetList[index])
               hetFunc('out',index,outHetPath)

               vcf_calc.run(resHetList[index])
               hetFunc('res',index,resHetPath)


#=========================================================================================#


        def listOutputFunc(list1): 
                           
            for z,m in enumerate(list1):
                if not '--out-prefix' in m:
                   listnew.append(m) 

                elif all(x in m for x in ['out','windowed-weir-fst']):
                     print '===================out-win-weir===================='
                     print m
                     outWinWeirList.append(m) 

                elif all(x in m for x in ['out','weir-fst']):
                     print '===================out-weir===================='
                     print m
                     outWeirList.append(m)

                elif all(x in m for x in ['out','TajimaD']):
                     print '===================out-Tajima===================='
                     print m
                     outTajimaList.append(m) 

                elif all(x in m for x in ['out','pi']):
                     print '===================out-Pi===================='
                     print m
                     outPiList.append(m) 

                elif all(x in m for x in ['out','freq']):
                     print '===================out-Freq===================='
                     print m
                     outFreqList.append(m) 

                elif all(x in m for x in ['out','het']):
                     print '===================out-het===================='
                     print m
                     outHetList.append(m) 

                elif all(x in m for x in ['res','windowed-weir-fst']):
                     print '===================res-win-weir===================='
                     print m
                     resWinWeirList.append(m) 

                elif all(x in m for x in ['res','weir-fst']):
                     print '===================res-weir===================='
                     print m
                     resWeirList.append(m) 

                elif all(x in m for x in ['res','TajimaD']):
                     print '===================res-Tajima===================='
                     print m
                     resTajimaList.append(m) 

                elif all(x in m for x in ['res','pi']):
                     print '===================res-Pi===================='
                     print m
                     resPiList.append(m)

                elif all(x in m for x in ['res','freq']):
                     print '===================res-Freq===================='
                     print m
                     resFreqList.append(m) 

                elif all(x in m for x in ['res','het']):
                     print '===================res-Freq===================='
                     print m
                     resHetList.append(m)  
 
            for m in listnew:
                if not '--calc-statistic' in m:
                   print '===================def===================='
                   print m
                   defList.append(m)

                elif all(x in m for x in ['windowed-weir-fst']):
                     print '===================win-weir===================='
                     print m
                     winWeirList.append(m)

                elif all(x in m for x in ['weir-fst']):
                     print '===================weir===================='
                     print m
                     weirList.append(m)

                elif all(x in m for x in ['TajimaD']):
                     print '===================tajima===================='
                     print m
                     tajimaList.append(m)

                elif all(x in m for x in ['pi']):
                     print '===================pi===================='
                     print m
                     piList.append(m)

                elif all(x in m for x in ['freq']):
                     print '===================freq===================='
                     print m
                     freqList.append(m)

                elif all(x in m for x in ['het']):
                     print '===================het===================='
                     print m
                     hetList.append(m)

            

            for i in range(4):
                print i
                defaultFuncCall(i)       
                winWeirFuncCall(i)     
                weirFuncCall(i)
                tajimaFuncCall(i)
                piFuncCall(i)
                freqFuncCall(i)
                hetFuncCall(i)
             
  
        listOutputFunc(unmatched1)




#====================================================================================================#


if __name__ == '__main__':
    unittest.main()



