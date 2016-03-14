## siyguo@indiana.edu
## May,2013

from collections import defaultdict
import sys
import re
import numpy


START_SYMBOL = '*'

class HMM():

    def __init__(self, trainFileName, testFileName):
        self.ftrain = trainFileName
        self.ftest = testFileName

    def get_counts(self,smoothing = 'NONE'): # get counts for deriving parameters
        self.wordtag = defaultdict(int) # emission freqs
        self.unitag = defaultdict(int) # unigram freqs of tags
        self.bitag = defaultdict(int) # bigram freqs of tags
        self.tritag = defaultdict(int) # trigram freqs of tags
        
        tag_penult = START_SYMBOL
        tag_last = START_SYMBOL
        tag_current = START_SYMBOL
        
        for line in open(self.ftrain, 'r'):
            #print "printing a line"
            #print line

            # removing trailing spaces
            line = line.strip()
            
            # for each line 
            tag_penult = START_SYMBOL
            tag_last = START_SYMBOL
            tag_current = START_SYMBOL
            if(line):
                for l in line.split(' '):
                    #print "printing a word tag pair"
                    #print l
                    l = l.strip()

                    word,tag = '',''
                    try:
                        word, tag = l.split('/')
                    except ValueError: # case when there are more words not separated with space but joined with '/'
                        split = l.split('/')
                        word, tag = "/".join(numpy.array(split)[[0,-2]]), split[len(split)-1]

                    tag_penult = tag_last
                    tag_last = tag_current
                    tag_current = tag
                    
                    # update emission freqs
                    self.wordtag[(word,tag)] += 1
                    
                    # update unitag freqs
                    self.unitag[tag] += 1
                    
                    # update bitag freqs
                    self.bitag[(tag_last, tag_current)] += 1
                    
                    # update tritag freqs
                    self.tritag[(tag_penult, tag_last, tag_current)] += 1
                    
                    # update starting bigrams
                    if tag_last == START_SYMBOL and tag_penult == START_SYMBOL:
                        self.bitag[(START_SYMBOL,START_SYMBOL)] += 1
                    # adding two tags * * at start of each sentence    
                    self.unitag[START_SYMBOL] += 2
        self.alltag = sum(self.unitag.itervalues())             
                
        


    def get_e(self,word,tag):
        return float(self.wordtag[(word,tag)])/self.unitag[tag]

    def get_q(self,tag_penult, tag_last, tag_current):
        return float(self.tritag[(tag_penult, tag_last, tag_current)])/self.bitag[(tag_penult, tag_last)]

  
    def get_laplace_q(self,tag_penult, tag_last, tag_current):
        return float((self.tritag[(tag_penult, tag_last, tag_current)]) + 1)/(self.bitag[(tag_penult, tag_last)]+ len(self.tags)**2)
    
    def get_parameters(self, method='UNK',smoothing='NONE'): # derive parameters from counts
        self.words = set([key[0] for key in self.wordtag.keys()])
        if method == 'UNK':
            self.UNK()
        elif method == 'MORPHO':
            self.MORPHO()
        self.words = set([key[0] for key in self.wordtag.keys()])
        self.tags = set(self.unitag.keys()).union({'*'})
        print str(self.tags)
        self.E = defaultdict(int)
        self.Q = defaultdict(int)
        for (word,tag) in self.wordtag:
            self.E[(word,tag)] = self.get_e(word,tag)

            #print str(self.tri_alpha)                

        for (tag_penult, tag_last, tag_current) in self.tritag:
            if(smoothing == 'LAPLACE'):
                self.Q[(tag_penult, tag_last, tag_current)] = self.get_laplace_q(tag_penult, tag_last, tag_current)        
            else:
                self.Q[(tag_penult, tag_last, tag_current)] = self.get_q(tag_penult, tag_last, tag_current)

    def tagger(self, outFileName, method='UNK', smoothing='NONE'):
        # train
        sys.stderr.write("get_counts\n")
        self.get_counts()
        sys.stderr.write("get_parameters\n")
        self.get_parameters(method,smoothing)
        # begin tagging
        self.sent = []
        fout = open(outFileName, 'w')
        cnt = 0
        #try:

        for line in open(self.ftest, 'r'):            
            line = line.strip()
            if(line):
                cnt += 1;
                if(cnt==1001):
                    break
                for l in line.split(' '):
                    l = l.strip()                    
                    if l == '.':                        
                        if self.sent:
                            #sys.stderr.write("generate path for "+str(self.sent)+'\n')
                            if cnt%100 == 0: print cnt
                            path = self.viterbi(self.sent,cnt, method)
                            #print path
                            for i in range(len(self.sent)):
                                fout.write(self.sent[i]+'/'+path[i]+' ')
                            self.sent = []                        
                    else:                        
                        self.sent.append(l)
            fout.write('./.\n')
        #except KeyError:
         #   sys.stderr.write(str(KeyError) + " in line " + str(cnt) +"\n")
        fout.close()

    def viterbi(self,sent,cnt,method='UNK'):

        V = {}
        path = {}
        # init
        V[0,START_SYMBOL,START_SYMBOL] = 1
        path[START_SYMBOL,START_SYMBOL] = []
        # run
        #print "sentence is"
        #print sent

        #sys.stderr.write("entering k loop\n")
        for k in range(1,len(sent)+1):
            temp_path = {}
            word = self.get_word(sent,k-1)
            ## handling unknown words in test set using low freq words in training set
            if word not in self.words:
                #print word
                if method=='UNK':
                    word = "<RARE>"
                elif method == 'MORPHO':
                    word = self.subcategorize(word)
            #sys.stderr.write("entering u loop "+str(k)+"\n")
            for u in self.get_possible_tags(k-1):
                #sys.stderr.write("entering v loop "+str(u)+"\n")
                for v in self.get_possible_tags(k):
                    V[k,u,v],prev_w = max([(V[k-1,w,u] * self.Q[w,u,v] * self.E[word,v],w) for w in self.get_possible_tags(k-2)])
                    temp_path[u,v] = path[prev_w,u] + [v]
            path = temp_path
           

        # last step
        #print "before last step " + str(cnt)
        prob,umax,vmax = max([(V[len(sent),u,v] * self.Q[u,v,'.'],u,v) for u in self.get_possible_tags(len(sent)-1) for v in self.get_possible_tags(len(sent))])
        return path[umax,vmax]

    def get_possible_tags(self,k):
        if k == -1:
            return set([START_SYMBOL])
        if k == 0:
            return set([START_SYMBOL])
        else:
            return self.tags

    def get_word(self,sent,k):
        if k < 0:
            return START_SYMBOL
        else:
            return sent[k]

    def UNK(self):
        new = defaultdict(int)
        # change words with freq <6 into unknown words ""<RARE>""
        for (word,tag) in self.wordtag:
            new[(word,tag)] = self.wordtag[(word,tag)]
            if self.wordtag[(word,tag)] < 6:
                new[("<RARE>",tag)] += self.wordtag[(word,tag)]
        self.wordtag = new

    def subcategorize(self,word):
        if not re.search(r'\w',word):
            return '<PUNCS>'
        elif re.search(r'[A-Z]',word):
            return '<CAPITAL>'
        elif re.search(r'\d',word):
            return '<NUM>'
        elif re.search(r'(ion\b|ty\b|ics\b|ment\b|ence\b|ance\b|ness\b|ist\b|ism\b)',word):
            return '<NOUNLIKE>'
        elif re.search(r'(ate\b|fy\b|ize\b|\ben|\bem)',word):
            return '<VERBLIKE>'
        elif re.search(r'(\bun|\bin|ble\b|ry\b|ish\b|ious\b|ical\b|\bnon)',word):
            return '<JJLIKE>'
        else:
            return '<OTHER>'

    def MORPHO(self):
        new = defaultdict(int)
        for (word,tag) in self.wordtag:
            new[(word,tag)] = self.wordtag[(word,tag)]
            if self.wordtag[(word,tag)] < 6:
                new[(self.subcategorize(word),tag)] += self.wordtag[(word,tag)]
        self.wordtag = new

