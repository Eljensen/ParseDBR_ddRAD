'''
Created on Jan 11, 2018

@author: karl
'''
import unittest
import os
import time
from ParseFastQ import * 

class Test(unittest.TestCase):

    def testOpenFile(self):
        p = ParseFastQ("/home/karl/Downloads/P_Bear_ddRAD1_R1_subset.fastq.gz")
        temp = p.next() 
        expected ="@D00430:276:CBV22ANXX:2:2202:1493:1969 1:N:0:1" 
        self.assertEqual(temp['id'],expected, "ID's Match")
        pass
    
    def testRead(self):
        f = ParseFastQ("/home/karl/Downloads/P_Bear_ddRAD1_R1_subset.fastq.gz")
        r = Read(f.next())
        expected ="@D00430:276:CBV22ANXX:2:2202:1493:1969 1:N:0:1" 
        self.assertEqual(r.id,expected, "ID's Match")

    def testReadFindSeq(self):
        f = ParseFastQ("/home/karl/Downloads/P_Bear_ddRAD1_R1_subset.fastq.gz")
        r = Read(f.next())
        idx = r.find_seq("AGAA")
        self.assertEqual(idx,5,"Extracted Enzyme Recognition Sequence")

    def testReadFindSeqErr(self):
        f = ParseFastQ("/home/karl/Downloads/P_Bear_ddRAD1_R1_subset.fastq.gz")
        r = Read(f.next())
        self.assertEqual(r.find_seq("CGAA",max_err=1),5,"Check Sequence")
        self.assertEqual(r.find_seq("CCAA",max_err=2),4,"Check Sequence")
        self.assertEqual(r.find_seq("ACAA",max_err=1),5,"Check Sequence")

    def testRead1(self):
        f = ParseFastQ("/home/karl/Downloads/P_Bear_ddRAD1_R1_subset.fastq.gz")
        r1 = Read1(f.next())
        self.assertEqual(r1.bc,"AATC","Extracted Barcode")
        self.assertEqual(r1.ers,"AGAA","Extracted Enzyme Recognition Sequence")
        self.assertEqual(r1.insert,"TGAGCCGCAACTTCGGGATGAAAATGCTCACAATGACAAATCTGTCCACGGAGTGCTTAATCCAACTTACCAAGCTGGGTTACGACGCGACGCCGTTCAACCAGATATTGAAGCAGA","Extracted Insert")
        pass

    def testRead1Exception(self):
        f = ParseFastQ("/home/karl/Downloads/P_Bear_ddRAD1_R1_subset.fastq.gz")
        fq = f.next()
        fq['seq']= 'MBCDEFG' + fq['seq'][1:]
        try:
            r1 = Read1(fq)
        except ValueError:
            caught_error = 1
        self.assert_(caught_error, "Catch exception")
        pass

    def testSearchErsRead1(self):
        f = ParseFastQ("/home/karl/Downloads/P_Bear_ddRAD1_R1_subset.fastq.gz")
        fq = f.next()
        # Shift the ers by adding an A to the barcode
        fq['seq']= fq['seq'][0] + 'A' + fq['seq'][1:]
        r1 = Read1(fq)
        self.assertEqual(r1.bc,"AAATC","Extracted Barcode")
        self.assertEqual(r1.ers,"AGAA","Extracted Enzyme Recognition Sequence")
        self.assertEqual(r1.insert,"TGAGCCGCAACTTCGGGATGAAAATGCTCACAATGACAAATCTGTCCACGGAGTGCTTAATCCAACTTACCAAGCTGGGTTACGACGCGACGCCGTTCAACCAGATATTGAAGCAGA","Extracted Insert")
        pass

    def testSearchErsFailRead1(self):
        f = ParseFastQ("/home/karl/Downloads/P_Bear_ddRAD1_R1_subset.fastq.gz")
        fq = f.next()
        # Insert some garbage into the ESR
        fq['seq']= fq['seq'][:6] + 'C' + fq['seq'][7:]
        try:
            r1 = Read1(fq)
        except ValueError:
            caught_error = 1
        self.assert_(caught_error, "Catch exception")
        pass

    def testRead1String(self):
        f = ParseFastQ("/home/karl/Downloads/P_Bear_ddRAD1_R1_subset.fastq.gz")
        fast1 = f.next()
        r1 = Read1(fast1)
        
        fp = open("temp.fast",'w')
        fp.write(r1.to_str())
        fp.close()
        
        f1 = ParseFastQ("temp.fast")
        fast2 = f1.next()
        #print fast1
        #print fast2
        self.assertEqual(fast1, fast2, "Compare to string file")
        pass
    
    def testRead2(self):
        f = ParseFastQ("/home/karl/Downloads/P_Bear_ddRAD1_R2_subset.fastq.gz")
        f.next()
        r2 = Read2(f.next())
        #G
        # DBR: "AAAAACAA"
        # Anchor: GGACG
        # Index: ACTTGA
        # ERS: AATT
        # Insert ATACTGAGAGATACTTGAGGCAGCACCAGAAACGTCACATGCAATAACAGGATAGGAAAGCATGCGAGGTATTTGTCTAGTCAACAAGTATTTATTGAGCAC
        self.assertEqual(r2.dbr,"AAAAACAA")
        pass

    def testRead2Exception(self):
        f = ParseFastQ("/home/karl/Downloads/P_Bear_ddRAD1_R1_subset.fastq.gz")
        fq = f.next()
        try:
            r2 = Read2(fq)
        except ValueError:
            caught_error = 1
        self.assert_(caught_error, "Catch exception")
        pass
    
    def testCompareRead2(self):
        f = ParseFastQ("/home/karl/Downloads/P_Bear_ddRAD1_R2_subset.fastq.gz")
        f.next()
        r2 = Read2(f.next())
        r2_1 = r2
        self.assert_(r2.compare(r2_1))
        r2_2 = Read2(f.next())
        self.assertFalse(r2.compare(r2_2), "Incorrect compare result")
        pass
    
    def testRead2StripDBR(self):
        f = open('temp.fastq','w')
        f.write("@D00430:276:CBV22ANXX:2:2202:1400:1979 2:N:0:1\n"
                "GAAAAACAAGGACGACTTGAAATTATACTGAGAGATACTTGAGGCAGCACCAGAAACGTCACATGCAATAACAGGATAGGAAAGCATGCGAGGTATTTGTCTAGTCAACAAGTATTTATTGAGCAC\n"
                "+\n"
                "BBABBGGGGGGCFEEGGGGGGGGGGGGGGGFG>FFFGGGGGGGGG<GGEEBGGGGFGFFF>BGGGGGGGCCFGFE0FC@@EC>DGEGGGGBGGBFGGGGGGGGGFC0@GGBBGGGCGEGGGGGGGF\n")
        f.close()
        fq = ParseFastQ('temp.fastq')
        r2 = Read2(fq.next())
        f1 = open('temp1.fastq','w')
        f1.write(r2.trim_dbr_anchor())
        f1.close()
        
        fq1 = ParseFastQ('temp1.fastq')
        test = fq1.next()
        self.assertEqual(r2.id, test['id'], "ID Match")
        self.assertEqual("ACTTGAAATTATACTGAGAGATACTTGAGGCAGCACCAGAAACGTCACATGCAATAACAGGATAGGAAAGCATGCGAGGTATTTGTCTAGTCAACAAGTATTTATTGAGCAC", test['seq'], "Seq Match")
        self.assertEqual("EGGGGGGGGGGGGGGGFG>FFFGGGGGGGGG<GGEEBGGGGFGFFF>BGGGGGGGCCFGFE0FC@@EC>DGEGGGGBGGBFGGGGGGGGGFC0@GGBBGGGCGEGGGGGGGF", test['conf'], "Confidence Match")
        pass

    def testStripDuplicates(self):
        #Read1.ers = 'AGAA'
        fn_r1 = "/home/karl/Downloads/P_Bear_ddRAD1_R1_subset.fastq.gz"
        fq = ParseFastQ(fn_r1)
        r = Read(fq.next())
        a = StripIds()
        a.append_drop(r.id)
        a.strip(fn_r1,"out.fastq")
        
        fq1 = ParseFastQ("out.fastq")
        r1 = Read(fq1.next())
        self.assertNotEqual(r1.id, r.id, "Failed to remove first id")
        pass
    
    def trialRun(self):
        fn = "/home/karl/Downloads/P_Bear_ddRAD1_R2_subset.fastq.gz"
        a = StripIds()
        count1 = 0
        for fq in ParseFastQ(fn):
            count1 += 1
            try:
                r2 = Read2(fq)
            except ValueError:
                r = Read(fq)
                a.append_drop(r.id)
        
        a.strip(fn,"out.fastq")  

        count2 = 0
        for fq in ParseFastQ("out.fastq"):
            count2 += 1
        
        print "Count1: " + str(count1) + " Count2: " + str(count2) 
        self.assertNotEqual(count1, count2, "Didn't remove anything from files... seems wrong")
        
    def trialRunDuplicates(self):
        fn = "/home/karl/Downloads/P_Bear_ddRAD1_R2_subset.fastq.gz"
        fout = "out.fastq"
        r2f = Read2File(fn)
        r2f.decruff("R2_decruff.fastq")

        r2f_1 = Read2File("R2_decruff.fastq")
        r2f_1.remove_duplicates(fout)
        
    def testRemoveDuplicates(self):
        fn = "/home/karl/Downloads/P_Bear_ddRAD1_R2_subset.fastq.gz"
        fout = open("dups.fastq",'w')
        count = 0
        drop_cnt = 0
        a = StripIds()
        for fq in ParseFastQ(fn):
            try:
                r2 = Read2(fq)
                fout.write(r2.to_str())
                count+=1
                if(count % 12 == 0) : 
                    fout.write(r2.to_str());
                    a.append_drop(r2.id)
                if(count > 1000) : break;
            except ValueError:
                drop_cnt += 1
        r2f = Read2File("dups.fastq") 
        t1 = time.time()
        num_removed = r2f.remove_duplicates("out.fastq")
        t2 = time.time()
        print "Drop Count: " + str(drop_cnt)
        print "Num Removed " + str(num_removed)
        print "exec time: " + str(t2-t1)
        exp = len(a.drop_list)
        self.assertEqual(num_removed, exp, "Number removed and number dropped don't match")
        
if __name__ == "__main__":
    import sys;sys.argv = ['', 'Test.testOpenFile', 
                           'Test.testRead',
                           'Test.testReadFindSeq',
                           'Test.testReadFindSeqErr',
                           'Test.testRead1',
                           'Test.testRead1Exception',
                           'Test.testRead2',
                           'Test.testRead2Exception',
                           'Test.testCompareRead2',
                           'Test.testSearchErsRead1',
                           'Test.testSearchErsFailRead1',
                           'Test.testRead1String',
                           'Test.testRead2StripDBR',
                           'Test.trialRun',
     #                      'Test.trialRunDuplicates',
                           'Test.testRemoveDuplicates',
                           'Test.testStripDuplicates' # This test is last because it changes static vars
                           ]
    unittest.main()