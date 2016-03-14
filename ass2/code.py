from collections import defaultdict

unigram = defaultdict(int)
bigram = defaultdict(int)
trigram = defaultdict(int)
word_tags = defaultdict(int)

def bi_count(line_tags):
    for i in xrange(len(line_tags) -1):
        u = line_tags[i]
        v = line_tags[i+1]
        bigram[u,v] += 1


def tri_count(line_tags):
    for i in xrange(len(line_tags) - 2):
        u = line_tags[i]
        v = line_tags[i+1]
        w = line_tags[i+2]
        trigram[u,v,w] += 1

def uni_count(line_tags):
    for i in xrange(len(line_tags)):
        u = line_tags[i]
        unigram[u] += 1


with open('Brown_tagged_train.txt','r') as f:
    for line in f:
        line_tags = []
        line = line.split()
        for w in line:
            x = w.split('/')
            line_tags.append(x[1])
            word_tags[x[0],x[1]] += 1
        line_tags = ['*'] + ['*'] + line_tags 
        uni_count(line_tags)
        bi_count(line_tags)
        tri_count(line_tags)

def E(x,y):
    return float(word_tags[x,y]/unigram[y])

def Q(u,v,w):
    if bigram[u,v] == 0:
        return 0
    return float(trigram[u,v,w]/bigram[u,v])

START = defaultdict(int)
END = defaultdict(int)
START['*'] = 1
END['.'] = 1


def viterbi(line):
    l = len(line)
    u = START
    v = START
    w = START
    bp = defaultdict()
    pi = defaultdict(int)
    pi[0,'*','*'] = 1
    tags = []

    for i in xrange(1,l+1):
        if i >=1 :
            w = unigram
        if i >=2 :
            v = unigram
        if i >= 3:
            u = unigram
        if i == l+1:
            w = END

        for z,val1 in  w.items():
            for y,val2 in v.items():
                max_val = -1;
                for x,val3 in u.items():
                    temp = pi[i-1,y,x] * Q(x,y,z) * E(line[i-1],z)
                    if temp > max_val:
                        max_val = temp
                        bp[i,y,z]  =x

                pi[i,z,y] = max_val

    max_val  = -1
    for i in range(0,l):
        tags.append('pp')

    for z, val1 in w.items():
            for y, val2 in v.items():
                now = (pi[l-1,y,z]* Q(y,z,'.'))
                if now > max_val:
                       ans = z,y

    print ans
    tags[l-1] = ans[0]
    tags[l-2] = ans[1]
    i = l-1
    while(i>=2):
        v = tags[i]
        u = tags[i-1]
        tags[i-2] = bp[i,u,v]
        print u,v
        i = i-1
    return tags


print unigram

    
