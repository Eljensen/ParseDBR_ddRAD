'''
Created on Jan 11, 2018
@author: karl
'''
import gzip
import re
import time
import argparse
import sys
from gzip import GzipFile

class ParseFastQ(object):
    ''' Accepts filename and provides basic incremental read functionality '''
    def __init__(self, filename):
        ''' Constructor '''
        self.fp = []
        if filename.endswith('.gz') :
            self.fp  = gzip.open(filename, 'rb')
        else :
            self.fp  = open(filename, 'rb')

        self.num_read = 0
        
    def __iter__(self):
        return self
    
    def next(self):
        ''' Get the next sequence object '''
        retval = {}
        retval['id'] = self.get_line() 
        retval['seq'] = self.get_line()
        self.get_line()
        retval['conf'] = self.get_line()
        self.num_read += 1
        return retval
    
    def tell(self):
        return self.fp.tell()
    
    def seek(self,x):
        self.fp.seek(x)

    def get_line(self):
        line = self.fp.readline().rstrip('\n')
        if( line == ''): 
            raise StopIteration
        return line
    
class Read(object):
    def __init__(self,FastQ):
        self.seq = FastQ['seq']
        self.id = FastQ['id']
        self.conf = FastQ['conf']

    def find_seq(self,seq,start=0,end=-1,max_err=0):
        retval = len(self.seq) + 1
        if max_err > 0 :
            dots = '.'*max_err
            for i in range(0,len(seq)-max_err + 1):
                pat = seq[:i] + dots + seq[i+max_err:]
                p = re.compile(pat)
                idx = p.search(self.seq[start:end])
                if idx:
                    if(retval >= idx.span()[0] + start):
                        retval = idx.span()[0] + start
                else:
                    raise ValueError("Sequence Not Found")
        else :
            retval = self.seq.find(seq,start,end)
            if retval < 0:
                raise ValueError("Sequence Not Found")

        return retval

    def check_id(self,seq_id):
        if self.id[:-7] == seq_id[:-7] :
            return True
        else : 
            return False

    def to_str(self):
        return self.id + '\n' + self.seq + '\n+\n' + self.conf + '\n'
    
class Read1(Read):
    ''' Parse a FastQ Dictionary into R1 components '''
    maxbclen = 8    # Max Barcode length
    minbclen = 4    # Min Barcode length
    
    ers = 'AGAA' # Enzyme recognition sequence
    anchor = 'N'
    
    def __init__(self,FastQ):
        '''
        Constructor
        '''
        super(Read1,self).__init__(FastQ)

        # check if ERS is present in the first 10 bp
        #idx = self.find_ers()
        idx = 0
        self.insert=self.seq[idx+len(Read1.ers):]
        self.bc=self.seq[len(Read1.anchor):idx]
        #if self.seq[:len(Read1.anchor)] != Read1.anchor : 
        #    raise ValueError("Invalid Anchor");

    def find_ers(self):
        # TODO: Support searching list of ERS
        low = len(Read1.anchor)+Read1.minbclen
        high = len(Read1.anchor) + Read1.maxbclen + len(Read1.ers)
        try:
            ers_idx = self.find_seq(Read1.ers,low,high,max_err=1)
        except ValueError:
            raise ValueError("Enzyme Recognition Sequence Not Found")

        return ers_idx
    

class Read2(Read):
    ''' Parse a FastQ Dictionary into R2 components '''
    hdrlen = 1    # First anchor length
    dbrlen = 8    # Degenerate Base Region length
    anchor = 'GGACG' # TODO: this is a list of possible anchors
    index = 'ACTTGA' # TODO: this is also a list of size 4
                     # Come up with way of latching on first read
    ers = 'AATT'  # Enzyme recognition sequence
    # TODO: the initial letter is always known but will change depending

    # Offsets to speed of runtime
    os_dbr = hdrlen
    os_anchor = os_dbr + dbrlen
    os_index = os_anchor + len(anchor)
    os_ers = os_index + len(index)
    os_insert = os_ers+len(ers)
    
    def __init__(self,FastQ):
        ''' Constructor '''
        super(Read2,self).__init__(FastQ)
        self.dbr = self.seq[Read2.os_dbr:Read2.os_anchor]
        self.insert= self.seq[Read2.os_insert:]
        if self.seq[0] != 'G' : 
            raise ValueError("Invalid Read 2 Sequence Start Character");
        
        if self.seq[Read2.os_ers:Read2.os_ers+len(Read2.ers)] != Read2.ers :
            raise ValueError("Invalid Enzyme Recognition Sequence")
        
        # Check the index
        try:
            self.find_seq(self.index,self.os_anchor,self.os_index+len(self.index))
        except ValueError:
            raise ValueError("Index not found")

    def compare(self,read): 
        ''' Check for matching dbr and inserts '''
        if self.dbr == read.dbr :
            return self.insert[:10] == read.insert[:10]
        else :
            return False
    
    def find_anchor(self):
        return 9

    def trim_dbr_anchor(self):
        idx = self.find_anchor()
        temp_str = self.id + '\n'
        temp_str += self.seq[idx+len(Read2.anchor):] + '\n'
        temp_str += '+\n'
        temp_str += self.conf[idx+len(Read2.anchor):] + '\n'
        return temp_str

class StripIds(object):
    ''' Takes a list of ID's and removes them from a FastQ file''' 
    def __init__(self):
        self.drop_list = []; 
        self.fp = ''
    
    def strip(self,filein,fileout):
        fp = open(fileout,'w')
        for fq in ParseFastQ(filein):
            r = Read(fq)
            if not r.id[:-7] in self.drop_list:
                fp.write(r.to_str())
        fp.close()
        
    def append_drop(self,seq_id):
        if not seq_id[:-7] in self.drop_list:
            self.drop_list.append(seq_id[:-7])

class Read2File(object):
    ''' Read2 actions class.  Does things like:
        + remove duplicates
        + remove malformed R2 sequences'''
    def __init__(self,filename):
        self.fn = filename
        
    def decruff(self,fileout):
        a = StripIds()
        for fq in ParseFastQ(self.fn):
            try:
                Read2(fq)
            except ValueError:
                r = Read(fq)
                a.append_drop(r.id)
        
        a.strip(self.fn,fileout)  
        
    def remove_duplicates(self,fileout):
        a = StripIds()
        masterf = ParseFastQ(self.fn)
        fq_sub = ParseFastQ(self.fn)
         # run time with fq_sub being reopened each time 1184.783 which saved about 40 seconds...
        for fq in masterf:
            # Save our current point in the file
            cur_loc = masterf.tell()
            # print "cur_loc: " + str(cur_loc)
            r2 = Read2(fq)
            
            # Jump to the same location as we were just at
            fq_sub.seek(cur_loc)
            
            # Run through the rest of the new file to find duplicates
            for fq_1 in fq_sub:
                r2_1 = Read2(fq_1)
                if r2_1.compare(r2) :
                    a.append_drop(r2_1.id)
        
        print "Num Duplicates Found: " + str(len(a.drop_list))    
        print "Stripping file..."
        a.strip(self.fn,fileout)
        return len(a.drop_list)
    
def VerifyRead1(fin,fout,fdrop):
    fd = open(fdrop,'w') 
    fo = open(fout,'w') 
    count = 0
    for fq in ParseFastQ(fin):
        try:
            r = Read1(fq)
            fo.write(r.to_str())
        except ValueError:
            r = Read(fq);
            fd.write(fq['id']+'\n')
            count += 1
            
    print "Dropped Read 1 sequences: " + str(count)
    fd.close()
    fo.close()

def StripFile(fdrop,fin,fout,fin_read_type=1,nozip=False):
    tbr = {}
    count = 0
    if fin_read_type == 1:
        fout_read_type = 2;
    else:
        fout_read_type = 1;
        
    #TODO: Refactor this if statement
    if(nozip):
        with open(fdrop) as f:
            for seq_id in f:
                pat = str(fout_read_type) + ":N:0:"
                tmp = re.sub(pat,str(fin_read_type) + ":N:0:",seq_id.rstrip())
                tbr[tmp] = 1
    else:
        with gzip.open(fdrop) as f:
            for seq_id in f:
                pat = str(fout_read_type) + ":N:0:"
                tmp = re.sub(pat,str(fin_read_type) + ":N:0:",seq_id.rstrip())
                tbr[tmp] = 1
    
    if(nozip):
        r1_out = open(fout,'w')
    else:
        r1_out = GzipFile(fout,'w')
    
    for fq in ParseFastQ(fin):
        if fq['id'] not in tbr:
            r1_out.write(Read(fq).to_str())
        else:
            count += 1
    print "Removed: " + str(count)
    r1_out.close()

def FindDuplicates(fin,fout,fdrop,nozip=False):
    count = 0
    bad = 0
    dup = 0
    if(nozip):
        fout = open(fout,'w')
        fdup = open(fdrop, 'w')
    else:
        fout = GzipFile(fout,'w')
        fdup = GzipFile(fdrop, 'w')
    t1 = time.time()
    dbr_table = {}
    for fq in ParseFastQ(fin):
        count += 1
        try:
            r = Read2(fq)
            if (r.dbr,r.insert) not in dbr_table:
                dbr_table[r.dbr,r.insert] = r.id
                fout.write(r.trim_dbr_anchor())
            else :
                dup += 1
                fdup.write(fq['id'] + '\n')
        except ValueError:
            fdup.write(fq['id'] + '\n')
            bad += 1
    fout.close()
    fdup.close()
            
    t2 = time.time()
    print "Total parsed: " + str(count)
    print "Total dropped: " + str(bad) 
    print "Total duplicate: " + str(dup) 
    print "Time: " + str(t2-t1)
    print "Time per read: " + str((t2-t1)/count)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="Reads and R1 and R2 FastQ file and removes the paired reads from both files if the R2 has duplicates"
        )
    parser.add_argument('--read1','-r', help='Read 1 fastq file', required=True,)
    parser.add_argument('--read2','-R', help='Read 2 fastq file', required=True)
    parser.add_argument('--index','-i', help='Read 2 index sequence', required=True)
    parser.add_argument('--enzyme','-e', help='Read 2 enzyme sequence', required=True)
    parser.add_argument('--out1','-n', help='Read 1 fastq output file', required=True)
    parser.add_argument('--out2','-N', help='Read 1 fastq output file', required=True)
    parser.add_argument('--drop', help='Drop file name',default="drop_list",required=False)
    parser.add_argument('--nogzip','-Z', help='Do not zip output files',action='store_true')
    parser.add_argument('--hdrlen','-l', help='Read 2 first anchor length', required=False,default=1)
    args = parser.parse_args()
    
    if(args.nogzip == False and args.drop == "drop_list"):
        args.drop = "drop_list.gz"
        
    # If read1 sequences need to be checked use this:
    check_read1 = 0

    Read2.hdrlen = args.hdrlen
    Read2.index = args.index
    Read2.ers = args.enzyme

    if(check_read1 == 1):
        VerifyRead1(args.read1, 'r1_tmp', 'r1_drop')
        StripFile('r1_drop', args.read2, 'r2_tmp', 2)
        FindDuplicates('r2_tmp', args.out2, args.drop)
        StripFile(args.drop, 'r1_tmp', args.out1, 1)
    else:
        FindDuplicates(args.read2, args.out2, args.drop,args.nogzip)
        StripFile(args.drop, args.read1, args.out1, 1, args.nogzip)
        
    sys.stdout.flush()
    sys.stderr.flush()
